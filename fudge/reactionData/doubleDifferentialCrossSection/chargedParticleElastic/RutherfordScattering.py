# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import math

from pqu import PQU as PQUModule

from xData import ancestry as ancestryModule

from fudge.core.math import fudgemath as fudgemathModule
from fudge.productData.distributions import angular as angularModule

from . import misc as miscModule

class RutherfordScattering( ancestryModule.ancestry ):
    """
    Stores the pure-Coulomb contribution to charged-particle elastic scattering
    """

    moniker = "RutherfordScattering"

    def __init__(self, domainMin=None, domainMax=None, domainUnit=None ):

        ancestryModule.ancestry.__init__(self)

        self.domainMin = domainMin
        self.domainMax = domainMax
        self.domainUnit = domainUnit

    @property
    def etaCoefficient( self ) :

        return( self.ancestor.etaCoefficient )

    @property
    def spin( self ) :

        return( self.ancestor.spin )

    @property
    def kCoefficient( self ) :

        return( self.ancestor.kCoefficient )

    @property
    def identicalParticles( self ) :

        return( self.ancestor.identicalParticles )

    def convertUnits( self, unitMap ):

        if self.domainUnit in unitMap:
            newUnit = unitMap[ self.domainUnit ]
            factor = PQUModule.PQU( 1., self.domainUnit ).getValueAs( newUnit )
            self.domainMin *= factor
            self.domainMax *= factor
            self.domainUnit = newUnit

    def dSigma_dMu( self, energy, accuracy = 1e-3, muMax = 0.999 ) :
        """
        Returns d(Sigma)/d(mu) at the specified incident energy.

        :param energy:      Energy of the projectile.
        :param accuracy:    The accuracy of the returned *dSigma_dMu*.
        :param muMax:       Slices the upper domain mu to this value.

        :return:            d(Sigma)/d(mu) at *energy*.
        """

        class tester :

            def __init__( self, evaluate, energy, relativeTolerance, absoluteTolerance ) :

                self.evaluate = evaluate
                self.energy = energy
                self.relativeTolerance = relativeTolerance
                self.absoluteTolerance = absoluteTolerance

            def evaluateAtX( self, mu ) :

                return( self.evaluate( self.energy, mu ) )

        muMin = -1.0
        if( self.identicalParticles ) : muMin = 0.0
        _dSigma_dMu = [ [ muMin, self.evaluate( energy, muMin ) ], [ muMax, self.evaluate( energy, muMax ) ] ]
        average = sum( y for x, y in _dSigma_dMu ) / len( _dSigma_dMu )
        _tester = tester( self.evaluate, energy, accuracy, accuracy * average )
        _dSigma_dMu = 2.0 * math.pi * angularModule.XYs1d( fudgemathModule.thickenXYList( _dSigma_dMu, _tester, biSectionMax = 16 ) )
        _dSigma_dMu.axes = self.defaultAxes( self.domainUnit )

        return( _dSigma_dMu )

    def evaluate( self, energy, mu, phi = 0.0 ) :
        """
        Returns the dSigma/dOmega(energy,mu,phi) for Rutherford scattering. Note, Rutherford scattering is independent of phi
        but phi is listed as an argument to be consistent with others dSigma/dOmega that may depend on phi.

        :param energy: Incident energy.
        :param mu: Scattering angle cosine.
        :param phi: The scattering azimuthal angle.

        :return: Differential cross section dSigma/dOmega at point energy,mu,phi (in b/sr).
        """

        eta = self.etaCoefficient / math.sqrt( energy )
        k = self.kCoefficient * math.sqrt( energy )
        if( self.identicalParticles ) :
            twiceSpin = 2.0 * float( self.spin )
            term1 = ( 2 * eta**2 ) / ( k**2 * ( 1.0 - mu**2 ) )
            term2 = ( 1 + mu**2 ) / ( 1 - mu**2 ) + (-1)**twiceSpin / ( twiceSpin + 1 ) * math.cos( eta * math.log( ( 1 + mu ) / ( 1 - mu ) ) )
            return term1 * term2
        else:
            return eta**2 / (k**2 * (1-mu)**2 )

    def toXMLList( self, indent='', **kwargs ):

        attrStr = ''
        for attribute in ('domainMin','domainMax','domainUnit'):
            if getattr(self,attribute) is not None: attrStr += ' %s="%s"' % (attribute, getattr(self,attribute))
        return ['%s<%s%s/>' % (indent, self.moniker, attrStr)]

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        domainMin, domainMax = element.get('domainMin'), element.get('domainMax')
        if domainMin is not None:
            domainMin, domainMax = map(float, (domainMin, domainMax) )
        return RutherfordScattering( domainMin, domainMax, element.get('domainUnit') )

    @staticmethod
    def defaultAxes( energyUnit, crossSectionUnit = 'b' ) :

        return( miscModule.defaultAxes( energyUnit, crossSectionUnit ) )
