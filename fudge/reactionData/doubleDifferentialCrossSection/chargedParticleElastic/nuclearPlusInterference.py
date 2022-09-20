# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
"""
nuclearPlusInterference separates charged-particle elastic scattering into two terms:
 - (cross section) - (pure Rutherford cross section), integrated from mu = -1 up to a cut-off angle mu_cutoff
 - the distribution (P(mu|E) - P(mu|E)_Rutherford) / cross section), defined for mu = -1 to mu_cutoff
"""

import math

from LUPY import ancestry as ancestryModule

from fudge.reactionData import crossSection as crossSectionModule
from fudge.productData.distributions import angular as angularModule

from . import misc as miscModule

class CrossSection(miscModule.ChargedParticleElasticTerm):

    moniker = 'crossSection'

    allowedDataForms = (crossSectionModule.XYs1d, crossSectionModule.Regions1d)

class Distribution( miscModule.ChargedParticleElasticTerm):

    moniker = 'distribution'

    allowedDataForms = (angularModule.XYs2d, angularModule.Regions2d)


class NuclearPlusInterference( ancestryModule.AncestryIO ):

    moniker = "nuclearPlusInterference"

    def __init__( self, muCutoff, crossSection = None, distribution = None ):

        ancestryModule.AncestryIO.__init__( self )
        self.__muCutoff = muCutoff
        self.crossSection = crossSection
        self.distribution = distribution

    @property
    def muCutoff(self): return self.__muCutoff

    @property
    def crossSection( self ):
        return self.__crossSection

    @crossSection.setter
    def crossSection(self, value):
        if not isinstance(value, CrossSection):
            raise TypeError( '%s cross section must be instance of %s' % (self.moniker, CrossSection.moniker))
        value.setAncestor( self )
        self.__crossSection = value

    @property
    def distribution( self ):
        return self.__distribution

    @distribution.setter
    def distribution(self, value):
        if not isinstance(value, Distribution):
            raise TypeError( '%s distribution must be instance of %s' % (self.moniker, Distribution.moniker))
        value.setAncestor( self )
        self.__distribution = value

    @property
    def domainMin(self): return self.crossSection.data.domainMin

    @property
    def domainMax(self): return self.crossSection.data.domainMax

    @property
    def domainUnit(self): return self.crossSection.data.domainUnit

    @property
    def identicalParticles( self ) :

        return( self.ancestor.identicalParticles )

    def check( self, info ):
        """
        Things to check: mu_cutoff should be between 0.99 and 1.0,
        domain of effective xsc / distribution should match,
        if identical particles distribution should only be for 0->muCutoff, otherwise -1->muCutoff
        distribution should be normalized
        If the cross section goes negative, also check that the pure Coulomb portion is large enough to 'win'
        :return:
        """
        return []
        raise NotImplementedError

    def convertUnits( self, unitMap ):

        self.crossSection.convertUnits( unitMap )
        self.distribution.convertUnits( unitMap )

    def dSigma_dMu(self, energy, accuracy=1e-3, muMax=None, probability=False):
        """
        Returns d(Sigma)/d(mu) at the specified incident energy if probability is **False** and P(mu) otherwise.

        :param energy:      Energy of the projectile.
        :param accuracy:    Currently not used. Only need to be compatible with other *dSigma_dMu* methods.
        :param muMax:       Slices the upper domain mu to this value.
        :param probability: If **True** P(mu) is returned instead of d(Sigma)/d(mu).

        :return:            d(Sigma)/d(mu) at *energy* of P(mu) if probability is **True**.
        """

        distribution = self.distribution.data.evaluate( energy )
        
        if( muMax is not None ) : distribution = distribution.domainSlice( domainMax = muMax )

        _dSigma_dMu = distribution
        if not probability: _dSigma_dMu *= self.crossSection.data.evaluate(energy)
        _dSigma_dMu.axes = self.defaultAxes( self.crossSection.data.domainUnit )

        data = _dSigma_dMu.copyDataToXYs()
        if self.identicalParticles:
            if data[1][0] > 2e-8:
                data[0][0] = 1e-8
            else:
                del data[0]
            data = [[-1.0, 0.0], [0.0, 0.0]] + data

        miscModule.fixMuRange(data)
        _dSigma_dMu.setData(data)

        return _dSigma_dMu

    def evaluate( self, E, mu, phi = 0.0 ) :
        """
        :param E: incident energy
        :param mu: scattering angle cosine
        :return: differential cross section at E,mu (integrated over phi) in b/sr
        """

        if self.identicalParticles:
            mu = abs(mu)        # distribution is symmetric, so only positive mu is stored
        if mu > self.muCutoff: return 0
        if E < self.crossSection.data.domainMin: return 0
        if E > self.crossSection.data.domainMax:
            raise ValueError( "Attempted to evaluate at %s, outside of domain limit %s %s" % ( E, self.crossSection.data.domainMax, self.crossSection.data.domainUnit ) )
        angular = self.distribution.data.evaluate( E )
        return abs( self.crossSection.data.evaluate( E ) ) * angular.evaluate( mu ) / ( 2 * math.pi )

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent+kwargs.get('incrementalIndent','  ')
        xml = [ '%s<%s muCutoff="%s">' % (indent, self.moniker, self.muCutoff) ]
        xml += self.crossSection.toXML_strList( indent2, **kwargs )
        xml += self.distribution.toXML_strList( indent2, **kwargs )
        xml[-1] += '</%s>' % self.moniker
        return xml

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )

        Xsc = Dist = None
        for child in element:
            if child.tag == CrossSection.moniker:
                Xsc = CrossSection.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            elif child.tag == Distribution.moniker:
                Dist = Distribution.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            else:
                raise TypeError("Unknown child element '%s' encountered in %s" % (child.tag, element.tag))

        NPCI = cls(float(element.get('muCutoff')), crossSection=Xsc, distribution=Dist)

        xPath.pop()

        return NPCI

    @staticmethod
    def defaultAxes( energyUnit, crossSectionUnit = 'b' ) :

        return( miscModule.defaultAxes( energyUnit, crossSectionUnit ) )
