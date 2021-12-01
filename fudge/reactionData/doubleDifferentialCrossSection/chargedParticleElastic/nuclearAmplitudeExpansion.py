# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
"""
The nuclearAmplitudeExpansion represents charged-particle elastic scattering using four terms:
 - pure Coulomb scattering  (i.e. RutherfordScattering),
 - pure nuclear scattering
 - complex interference between Coulomb and nuclear terms, broken into real and imaginary components
"""

import math
import numpy

from xData import ancestry as ancestryModule
from xData.series1d import Legendre, LegendreSeries

from fudge.core.math import fudgemath as fudgemathModule
from fudge.productData.distributions import angular as angularModule

from . import misc as miscModule

__metaclass__ = type

class nuclearTerm( miscModule.chargedParticleElasticTerm ):

    moniker = 'nuclearTerm'

    allowedDataForms = (angularModule.XYs2d, angularModule.regions2d)

class realInterferenceTerm( miscModule.chargedParticleElasticTerm ):

    moniker = 'realInterferenceTerm'

    allowedDataForms = (angularModule.XYs2d, angularModule.regions2d)

class imaginaryInterferenceTerm( miscModule.chargedParticleElasticTerm ):

    moniker = 'imaginaryInterferenceTerm'

    allowedDataForms = (angularModule.XYs2d, angularModule.regions2d)


class nuclearAmplitudeExpansion( ancestryModule.ancestry):

    moniker = 'nuclearAmplitudeExpansion'

    def __init__( self, nuclearTerm = None, realInterference = None, imaginaryInterference = None ):

        ancestryModule.ancestry.__init__( self )
        self.nuclearTerm = nuclearTerm
        self.realInterferenceTerm = realInterference
        self.imaginaryInterferenceTerm = imaginaryInterference

    @property
    def etaCoefficient( self ) :

        return( self.ancestor.etaCoefficient )

    @property
    def nuclearTerm( self ):
        return self.__nuclear

    @nuclearTerm.setter
    def nuclearTerm(self, value):
        if not isinstance(value, nuclearTerm):
            raise TypeError( '%s nuclear term must be instance of %s' % (self.moniker, nuclearTerm.moniker))
        value.setAncestor( self )
        self.__nuclear = value

    @property
    def realInterferenceTerm( self ):
        return self.__realInterference

    @realInterferenceTerm.setter
    def realInterferenceTerm(self, value):
        if not isinstance(value, realInterferenceTerm):
            raise TypeError( '%s real interference term must be instance of %s' %
                             (self.moniker, realInterferenceTerm.moniker))
        value.setAncestor( self )
        self.__realInterference = value

    @property
    def imaginaryInterferenceTerm( self ):
        return self.__imaginaryInterference

    @imaginaryInterferenceTerm.setter
    def imaginaryInterferenceTerm(self, value):
        if not isinstance(value, imaginaryInterferenceTerm):
            raise TypeError( '%s imaginary interference term must be instance of %s' %
                             (self.moniker, imaginaryInterferenceTerm.moniker))
        value.setAncestor( self )
        self.__imaginaryInterference = value

    @property
    def domainMin(self) :

        return self.nuclearTerm.data.domainMin

    @property
    def domainMax(self) :

        return self.nuclearTerm.data.domainMax

    @property
    def domainUnit(self) :

        return self.nuclearTerm.data.domainUnit

    @property
    def identicalParticles( self ) :

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

        self.nuclearTerm.convertUnits( unitMap )
        self.realInterferenceTerm.convertUnits( unitMap )
        self.imaginaryInterferenceTerm.convertUnits( unitMap )

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

        _dSigma_dMu = []
        numberOfPoints = 20
        for i1 in range( numberOfPoints ) :
            mu = muMin + i1 * ( muMax - muMin ) / numberOfPoints
            _dSigma_dMu.append( [ mu, self.evaluate( energy, mu ) ] )
        _dSigma_dMu.append( [ muMax, self.evaluate( energy, muMax ) ] )
        maximum = max( [ abs( y ) for x, y in _dSigma_dMu ] )

        _tester = tester( self.evaluate, energy, accuracy, accuracy * maximum )
        _dSigma_dMu = 2.0 * math.pi * angularModule.XYs1d( fudgemathModule.thickenXYList( _dSigma_dMu, _tester, biSectionMax = 12 ) )
        _dSigma_dMu.axes = self.defaultAxes( self.domainUnit )

        return( _dSigma_dMu )

    def evaluate( self, E, mu, phi = 0.0 ) :
        """
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

            P_lAtMu = [Legendre(l, mu) for l in range(NL)]

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
            RI = LegendreSeries( coefficients['realInterferenceTerm'] ).evaluate( mu )
            II = LegendreSeries( coefficients['imaginaryInterferenceTerm'] ).evaluate( mu )

            interference = -2*eta / (1-mu) * numpy.real(numpy.exp(1j * eta * numpy.log((1-mu)/2)) * (RI + 1j*II))

            nuclear = LegendreSeries( coefficients['nuclearTerm'] ).evaluate( mu )

        return( nuclear + interference )

    def toXMLList(self, indent='', **kwargs):

        indent2 = indent+kwargs.get('incrementalIndent','  ')
        xml = [ '%s<%s>' % (indent, self.moniker) ]
        xml += self.nuclearTerm.toXMLList( indent2, **kwargs )
        xml += self.realInterferenceTerm.toXMLList( indent2, **kwargs )
        xml += self.imaginaryInterferenceTerm.toXMLList( indent2, **kwargs )
        xml[-1] += '</%s>' % self.moniker
        return xml

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( element.tag )
        for child in element:
            if child.tag == nuclearTerm.moniker:
                nuclear = nuclearTerm.parseXMLNode( child, xPath, linkData )
            elif child.tag == realInterferenceTerm.moniker:
                real = realInterferenceTerm.parseXMLNode( child, xPath, linkData )
            elif child.tag == imaginaryInterferenceTerm.moniker:
                imag = imaginaryInterferenceTerm.parseXMLNode( child, xPath, linkData )
            else:
                raise TypeError("Unknown child element '%s' encountered in %s" % (child.tag, element.tag))

        CoulExp = nuclearAmplitudeExpansion( nuclearTerm=nuclear, realInterference=real, imaginaryInterference=imag)
        xPath.pop()
        return CoulExp

    @staticmethod
    def defaultAxes( energyUnit, crossSectionUnit = 'b' ) :

        return( miscModule.defaultAxes( energyUnit, crossSectionUnit ) )
