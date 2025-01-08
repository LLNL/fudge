# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
The double differential cross section for elastic scattering of two charged particle can be written
as the sum of three terms: pure Coulomb scattering, nuclear and an interference term (the 
interference term has real and imaginary components). This module contains classes for storing 
the nuclear and interference term as the class :py:class:`NuclearAmplitudeExpansion`.

This module contains the following classes:
        
    +-----------------------------------+-----------------------------------------------------------------------+
    | Class                             | Description                                                           |
    +===================================+=======================================================================+
    | NuclearTerm                       | This class stores the nuclear term.                                   |
    +-----------------------------------+-----------------------------------------------------------------------+
    | RealInterferenceTerm              | This class stores the real part of the interference term.             |
    +-----------------------------------+-----------------------------------------------------------------------+
    | ImaginaryInterferenceTerm         | This class stores the imaginary part of the interference term.        |
    +-----------------------------------+-----------------------------------------------------------------------+
    | NuclearAmplitudeExpansion         | This class contains the nuclear, and real and imaginary               |
    |                                   | interference terms as its members.                                    |
    +-----------------------------------+-----------------------------------------------------------------------+
"""

import math
import numpy

from LUPY import ancestry as ancestryModule

from xData import axes as axesModule
from xData import series1d as series1dModule

from fudge.core.math import fudgemath as fudgemathModule
from fudge.productData.distributions import angular as angularModule

from . import misc as miscModule


class NuclearTerm( miscModule.ChargedParticleElasticTerm ):
    """This class store the nuclear term."""

    moniker = 'nuclearTerm'

    allowedDataForms = (angularModule.XYs2d, angularModule.Regions2d)

class RealInterferenceTerm( miscModule.ChargedParticleElasticTerm ):
    """This class stores the real part of the interference term."""

    moniker = 'realInterferenceTerm'

    allowedDataForms = (angularModule.XYs2d, angularModule.Regions2d)

class ImaginaryInterferenceTerm( miscModule.ChargedParticleElasticTerm ):
    """This class stores the imaginary part of the interference term."""

    moniker = 'imaginaryInterferenceTerm'

    allowedDataForms = (angularModule.XYs2d, angularModule.Regions2d)


class NuclearAmplitudeExpansion( ancestryModule.AncestryIO):
    """
    This class contains the nuclear, and real and imaginary interference terms as its members.

    The following table list the primary members of this class: 
    
    +-----------------------+-----------------------------------------------------------+
    | Member                | Description                                               |
    +=======================+===========================================================+
    | nuclearTerm           | The data for the nuclear term.                            |
    +-----------------------+-----------------------------------------------------------+
    | realInterference      | The data for the real interference term.                  |
    +-----------------------+-----------------------------------------------------------+
    | imaginaryInterference | The data for the imaginary interference term.             |
    +-----------------------+-----------------------------------------------------------+
    """

    moniker = 'nuclearAmplitudeExpansion'

    def __init__( self, nuclearTerm = None, realInterference = None, imaginaryInterference = None ):
        """
        :param nuclearTerm:             The data for the nuclear term.
        :param realInterference:        The data for the real interference term.
        :param imaginaryInterference:   The data for the imaginary interference term.
        """

        ancestryModule.AncestryIO.__init__( self )
        self.nuclearTerm = nuclearTerm
        self.realInterferenceTerm = realInterference
        self.imaginaryInterferenceTerm = imaginaryInterference

    @property
    def etaCoefficient( self ) :
        r"""This method returns the parameter :math:`\eta \, \sqrt{E} = Z_1 \, Z_2 \, sqrt{\alpha^2 \mu m_1 / 2}`."""

        return( self.ancestor.etaCoefficient )

    @property
    def nuclearTerm( self ):
        """This method returns the nuclear term."""

        return self.__nuclear

    @nuclearTerm.setter
    def nuclearTerm(self, value):
        """
        This method sets the nuclear term to *value*.

        :param value:       The nuclear term.
        """

        if not isinstance(value, NuclearTerm):
            raise TypeError( '%s nuclear term must be instance of %s' % (self.moniker, NuclearTerm.moniker))
        value.setAncestor( self )
        self.__nuclear = value

    @property
    def realInterferenceTerm( self ):
        """This method returns the real part of the interference term."""

        return self.__realInterference

    @realInterferenceTerm.setter
    def realInterferenceTerm(self, value):
        """
        This method set the real part of the interference term to *value*.

        :param value:       The real part of the interference term.
        """

        if not isinstance(value, RealInterferenceTerm):
            raise TypeError( '%s real interference term must be instance of %s' %
                             (self.moniker, RealInterferenceTerm.moniker))
        value.setAncestor( self )
        self.__realInterference = value

    @property
    def imaginaryInterferenceTerm( self ):
        """This method returns the imaginary part of the interference term."""

        return self.__imaginaryInterference

    @imaginaryInterferenceTerm.setter
    def imaginaryInterferenceTerm(self, value):
        """
        This method set the imaginary part of the interference term to *value*.

        :param value:       The imaginary part of the interference term.
        """

        if not isinstance(value, ImaginaryInterferenceTerm):
            raise TypeError( '%s imaginary interference term must be instance of %s' %
                             (self.moniker, ImaginaryInterferenceTerm.moniker))
        value.setAncestor( self )
        self.__imaginaryInterference = value

    @property
    def domainMin(self) :
        """Returns the minimum projectile energy for *self*."""

        return self.nuclearTerm.data.domainMin

    @property
    def domainMax(self) :
        """Returns the maximum projectile energy for *self*."""

        return self.nuclearTerm.data.domainMax

    @property
    def domainUnit(self) :
        """Returns the energy unit of the projectile."""

        return self.nuclearTerm.data.domainUnit

    @property
    def identicalParticles( self ) :
        """This function returns True if the projectile and target are the same paritcle type and False otherwise."""

        return( self.ancestor.identicalParticles )

    def check( self, info ):
        """
        Check that incident energy domain of all three terms are the same
        FIXME: I'm not sure whether each term is meant to be normalized  (I see L=0 terms != 0)
        There may be more checks we can perform depending on whether particles are identical
        :return:
        """

        return []
        raise NotImplementedError()

    def convertUnits( self, unitMap ):
        """
        Converts all data in *self* per *unitMap*.
        
        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        self.nuclearTerm.convertUnits( unitMap )
        self.realInterferenceTerm.convertUnits( unitMap )
        self.imaginaryInterferenceTerm.convertUnits( unitMap )

    def dSigma_dMu(self, energy, accuracy=1e-3, muMax=0.999, probability=False):
        r"""
        This function returns :math:`d\sigma / d\mu` at the specified incident energy.

        :param energy:      Energy of the projectile.
        :param accuracy:    The accuracy of the returned *dSigma_dMu*.
        :param muMax:       Slices the upper domain mu to this value.
        :param probability: If **True** P(mu) is returned instead of d(Sigma)/d(mu).

        :return:            An instance of :pt:class:`angularModule.XYs1d`.
        """

        class Tester :
            r"""
            This class is used to added point to :math:`d\sigma / d\mu` until the desired lin-lin interpolation is met. 
            This class is for internal use.
            """

            def __init__( self, evaluate, energy, relativeTolerance, absoluteTolerance ) :
                """
                :param evaluate:            The function used to evaluate *self*.
                :param energy:              The energy of the projectile.
                :param relativeTolerance:   The relative tolerance used to determine if enough points have been added.
                :param absoluteTolerance:   The asbolute tolerance used to determine if enough points have been added.
                """

                self.evaluate = evaluate
                self.energy = energy
                self.relativeTolerance = relativeTolerance
                self.absoluteTolerance = absoluteTolerance

            def evaluateAtX( self, mu ) :
                r"""
                This function returns *self* evaluated at a projectile energy and outgoing :math:`\mu`.

                :param mu:  The mu point to evaluate *self* at.
                """

                return( self.evaluate( self.energy, mu ) )

        muMin = -1.0
        if( self.identicalParticles ) : muMin = 0.0

        _dSigma_dMu = []
        numberOfPoints = 20
        for i1 in range( numberOfPoints ) :
            mu = muMin + i1 * ( muMax - muMin ) / numberOfPoints
            _dSigma_dMu.append( [ mu, self.evaluate( energy, mu ) ] )
        _dSigma_dMu.append( [ muMax, self.evaluate( energy, muMax ) ] )
        maximum = max( [ abs( y ) for x, y in _dSigma_dMu ] )

        _tester = Tester( self.evaluate, energy, accuracy, accuracy * maximum )
        _dSigma_dMu = fudgemathModule.thickenXYList(_dSigma_dMu, _tester, biSectionMax=12)
        if self.identicalParticles:
            if _dSigma_dMu[1][0] > 2e-8:
                _dSigma_dMu[0][0] = 1e-8
            else:
                del _dSigma_dMu[0]
            _dSigma_dMu = [[-1.0, 0.0], [0.0, 0.0]] + _dSigma_dMu
        muEnd = _dSigma_dMu[-1][0]
        if muEnd < 1.0:
            muEnd = muEnd * (1 + 1e-6)
            if muEnd >= 1.0:
                muEnd = 0.5 * (_dSigma_dMu[-1][0] + 1.0)
            _dSigma_dMu += [[muEnd * (1 + 1e-6), 0.0], [1.0, 0.0]]
        _dSigma_dMu = 2.0 * math.pi * angularModule.XYs1d(_dSigma_dMu, axes=axesModule.Axes(2))
        _dSigma_dMu.axes = self.defaultAxes( self.domainUnit )

        if probability: _dSigma_dMu = _dSigma_dMu.normalize()

        return _dSigma_dMu

    def evaluate( self, E, mu, phi = 0.0 ) :
        """
        This function returns the differential cross section at E, mu (integrated over phi) in b/sr.

        :param E: incident energy
        :param mu: scattering angle cosine

        :return: differential cross section at E,mu (integrated over phi) in b/sr
        """

        if E < self.nuclearTerm.data.domainMin: return 0
        if E > self.nuclearTerm.data.domainMax:
            raise ValueError( "Attempted to evaluate at %s, outside of domain limit %s %s" % ( E, self.nuclearTerm.data.domainMax, self.nuclearTerm.data.domainUnit ) )


        eta = self.etaCoefficient / numpy.sqrt( E )

        # get Legendre series at desired energy, interpolating coefficients if necessary
        coefficients = {}
        for term in ( 'nuclearTerm', 'realInterferenceTerm', 'imaginaryInterferenceTerm' ) :
            data = getattr(self, term).data
            position, function1, function2, frac, interpolation, interpolationQualifier = data.getBoundingSubFunctions( E )
            if position == '<':
                coefficients[term] = [0]
            elif position == '=':
                coefficients[term] = function1.coefficients
            elif position == '>':
                raise ValueError("Attempted to evaluate outside of domain")
            else :                                  # Ignore interpolation as it may be "charged-particle" which does not make sense.
                coef1 = function1.coefficients[:]
                coef2 = function2.coefficients[:]
                coef1 += max( 0, len( coef2 ) - len( coef1 ) ) * [ 0 ]
                coef2 += max( 0, len( coef1 ) - len( coef2 ) ) * [ 0 ]
                coefficients[term] = [ ( 1.0 - frac ) * coef1[idx] + frac * coef2[idx] for idx in range( len( coef1 ) ) ]

        if self.identicalParticles:
            NL = len(coefficients['nuclearTerm'])
            NL_half = len(coefficients['realInterferenceTerm'])

            P_lAtMu = [series1dModule.Legendre(l, mu) for l in range(NL)]

            term1 = (1+mu) * numpy.exp(1j * eta * numpy.log((1-mu)/2))
            term2 = (1-mu) * numpy.exp(1j * eta * numpy.log((1+mu)/2))
            interferenceCoefficients = [ term1 + (-1)**ell * term2 for ell in range(NL) ]
            complexAs = [ real + 1j * imag for (real,imag) in zip(
                coefficients['realInterferenceTerm'], coefficients['imaginaryInterferenceTerm']) ]
            interferenceSum = numpy.sum( [interferenceCoefficients[ell] * (2*ell+1)/2 * complexAs[ell] * P_lAtMu[ell]
                                          for ell in range(NL_half)])
            interference = -2*eta / (1 - mu**2) * numpy.real( interferenceSum )

            nuclear = numpy.sum( [(4*ell+1)/2 * coefficients['nuclearTerm'][ell] * P_lAtMu[ell]
                                  for ell in range(NL)] )

        else:   # target != projectile
            RI = series1dModule.LegendreSeries( coefficients['realInterferenceTerm'] ).evaluate( mu )
            II = series1dModule.LegendreSeries( coefficients['imaginaryInterferenceTerm'] ).evaluate( mu )

            interference = -2*eta / (1-mu) * numpy.real(numpy.exp(1j * eta * numpy.log((1-mu)/2)) * (RI + 1j*II))

            nuclear = series1dModule.LegendreSeries( coefficients['nuclearTerm'] ).evaluate( mu )

        return( nuclear + interference )

    def toXML_strList(self, indent='', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent+kwargs.get('incrementalIndent','  ')
        xml = [ '%s<%s>' % (indent, self.moniker) ]
        xml += self.nuclearTerm.toXML_strList( indent2, **kwargs )
        xml += self.realInterferenceTerm.toXML_strList( indent2, **kwargs )
        xml += self.imaginaryInterferenceTerm.toXML_strList( indent2, **kwargs )
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
        for child in element:
            if child.tag == NuclearTerm.moniker:
                nuclear = NuclearTerm.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            elif child.tag == RealInterferenceTerm.moniker:
                real = RealInterferenceTerm.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            elif child.tag == ImaginaryInterferenceTerm.moniker:
                imag = ImaginaryInterferenceTerm.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            else:
                raise TypeError("Unknown child element '%s' encountered in %s" % (child.tag, element.tag))

        CoulExp = cls(nuclearTerm=nuclear, realInterference=real, imaginaryInterference=imag)

        xPath.pop()

        return CoulExp

    @staticmethod
    def defaultAxes( energyUnit, crossSectionUnit = 'b' ) :
        """
        This function returns an :py:class:`miscModule.defaultAxes` instance for a double difference cross section.

        :param energyUnit:          The unit for energy.
        :param crossSectionUnit:    The unit for cross section.
        """

        return( miscModule.defaultAxes( energyUnit, crossSectionUnit ) )
