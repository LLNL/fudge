# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

r"""
This module contains the classes for representing the nuclear + interference term of the elastic scattering of two nuclei.
The main class is the :py:class:`NuclearPlusInterference` class which stores the double differential cross section as
the product :math:`\sigma(E) \, P(\mu|E)`. The :math:`\sigma(E)` most be an instance of :py:class:`CrossSection` and
the :math:`P(\mu|E)` most be an instance of :py:class:`Distribution`.

This module contains the following classes:
        
    +---------------------------+-----------------------------------------------------------------------------------+
    | Class                     | Description                                                                       |
    +===========================+===================================================================================+
    | CrossSection              | This class represents the :math:`\sigma(E)` part of :math:`\sigma(E) \, P(\mu|E)`.|
    +---------------------------+-----------------------------------------------------------------------------------+
    | Distribution              | This class represents the :math:`P(\mu|E)` part of :math:`\sigma(E) \, P(\mu|E)`. |
    +---------------------------+-----------------------------------------------------------------------------------+
    | NuclearPlusInterference   | This is an exception class that is raised when a method of the                    |
    |                           | :py:class:`Form` cannot perform its calculation.                                  |
    +---------------------------+-----------------------------------------------------------------------------------+
"""

import math

from LUPY import ancestry as ancestryModule

from fudge.reactionData import crossSection as crossSectionModule
from fudge.productData.distributions import angular as angularModule

from . import misc as miscModule

class CrossSection(miscModule.ChargedParticleElasticTerm):
    """
    This class represents the :math:`\sigma(E)` part of :math:`\sigma(E) \, P(\mu|E)`
    """

    moniker = 'crossSection'

    allowedDataForms = (crossSectionModule.XYs1d, crossSectionModule.Regions1d)

class Distribution( miscModule.ChargedParticleElasticTerm):
    """
    This class represents the :math:`P(\mu|E)` part of :math:`\sigma(E) \, P(\mu|E)`.
    """

    moniker = 'distribution'

    allowedDataForms = (angularModule.XYs2d, angularModule.Regions2d)


class NuclearPlusInterference( ancestryModule.AncestryIO ):
    """
    This class represents the nuclear + interference term of the elastic scattering of two nuclei and
    stores the double differential cross section as the product :math:`\sigma(E) \, P(\mu|E)`.
    In :math:`P(\mu|E)` the :math:`\mu` ranges from muMin to muCutoff where muCutoff is a member of this
    class and muMin is -muCutoff for identical particles and -1 otherwise.

    The following table list the primary members of this class:

    +-------------------+-------------------------------------------------------------------------------+
    | Member            | Description                                                                   |
    +===================+===============================================================================+
    | muCutoff          | The maximum :math:`\mu` value the data represent.                             |
    +-------------------+-------------------------------------------------------------------------------+
    | crossSection      | The cross section part of the nuclear + interference term.                    |
    +-------------------+-------------------------------------------------------------------------------+
    | distribution      | The distribution part of the nuclear + interference term.                     |
    +-------------------+-------------------------------------------------------------------------------+
    """

    moniker = "nuclearPlusInterference"

    def __init__( self, muCutoff, crossSection = None, distribution = None ):
        """
        :param muCutoff:        The maximum :math:`\mu` value the data represent.
        :param crossSection:    The cross section part of the nuclear + interference term.
        :param distribution:    The distribution part of the nuclear + interference term.
        """

        ancestryModule.AncestryIO.__init__( self )
        self.__muCutoff = muCutoff
        self.crossSection = crossSection
        self.distribution = distribution

    @property
    def muCutoff(self):
        """This methods returns the value of *muCutoff*."""

        return self.__muCutoff

    @property
    def crossSection( self ):
        """This methods returns a reference to the cross section member."""

        return self.__crossSection

    @crossSection.setter
    def crossSection(self, value):
        """
        This methods set the cross section member to *value*.

        :param value:   The cross section.
        """

        if not isinstance(value, CrossSection):
            raise TypeError( '%s cross section must be instance of %s' % (self.moniker, CrossSection.moniker))
        value.setAncestor( self )
        self.__crossSection = value

    @property
    def distribution( self ):
        """This methods returns a reference to the distribution member."""

        return self.__distribution

    @distribution.setter
    def distribution(self, value):
        """
        This methods set the distribution member to *value*.

        :param value:   The distribution.
        """

        if not isinstance(value, Distribution):
            raise TypeError( '%s distribution must be instance of %s' % (self.moniker, Distribution.moniker))
        value.setAncestor( self )
        self.__distribution = value

    @property
    def domainMin(self):
        """Returns the minimum projectile energy for *self*.""" 

        return self.crossSection.data.domainMin

    @property
    def domainMax(self):
        """Returns the maximum projectile energy for *self*.""" 

        return self.crossSection.data.domainMax

    @property
    def domainUnit(self):
        """Returns the energy unit of the projectile."""

        return self.crossSection.data.domainUnit

    @property
    def identicalParticles( self ) :
        """This function returns True if the projectile and target are the same paritcle type and False otherwise."""

        return( self.ancestor.identicalParticles )

    def check( self, info ):
        """
        This method currently does nothing.
        """

        return []
        raise NotImplementedError

    def convertUnits( self, unitMap ):
        """
        Converts all data in *self* per *unitMap*.
        
        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        self.crossSection.convertUnits( unitMap )
        self.distribution.convertUnits( unitMap )

    def dSigma_dMu(self, energy, accuracy=1e-3, muMax=None, probability=False):
        """
        Returns :math:`d\sigma / d\mu` at the specified projdctile energy.
                
        :param energy:          Energy of the projectile.
        :param accuracy:        This argument is not used.
        :param muMax:           Slices the upper domain mu to this value.
        :param probability:     If True P(mu) is returned instead of d(Sigma)/d(mu).

        :return:                A :py:class:`angularModule.XYs1d` instance.
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
        r"""
        Returns the :math:`d\sigma / d \Omega(energy,\mu,\phi) for Rutherford scattering. Note, Rutherford
        scattering is independent of phi but phi is listed as an argument to be consistent with others
        :math:`d\sigma / d\Omega` that may depend on phi.

        :param energy:  Projectile energy.
        :param mu:      Scattering angle cosine.
        :param phi:     The scattering azimuthal angle.

        :returns:       A float.
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
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent+kwargs.get('incrementalIndent','  ')
        xml = [ '%s<%s muCutoff="%s">' % (indent, self.moniker, self.muCutoff) ]
        xml += self.crossSection.toXML_strList( indent2, **kwargs )
        xml += self.distribution.toXML_strList( indent2, **kwargs )
        xml[-1] += '</%s>' % self.moniker
        return xml

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        """
        Parse *element* into an instance of *cls*.
        
        :param cls:         Form class to return.
        :param element:     Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.
        
        :return: an instance of *cls* representing *element*.
        """

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
        """
        This function returns an :py:class:`axesModule.Axes` instance for a double difference cross section.

        :param energyUnit:          The unit for energy.
        :oaram crossSectionUnit:    The unit for cross section.
        """

        return( miscModule.defaultAxes( energyUnit, crossSectionUnit ) )
