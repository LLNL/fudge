# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Photon incoherent double difference cross section form and its supporting classes.
"""

import math

from pqu import PQU as PQUModule

from numericalFunctions import integration as nf_integrationModule

import xData.standards as standardsModule

from fudge.core.math import fudgemath as fudgemathModule
from ... import crossSection as crossSectionModule
from .. import base as baseModule

__metaclass__ = type


class XYs1d( baseModule.XYs1d ) :

    ENDFMT = 504

class regions1d( baseModule.regions1d ) :

    ENDFMT = 504

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs1d, ) )

class form( baseModule.form ) :

    moniker = 'incoherentPhotonScattering'
    subformAttributes = ( 'scatteringFunction', )

    def __init__( self, pid, label, productFrame, scatteringFunction ) :

        if( not( isinstance( scatteringFunction, ( XYs1d, regions1d ) ) ) ) :
            raise Exception( 'Instance is class "%s" and not scatteringFunction' % scatteringFunction.__class__ )

        baseModule.form.__init__( self, pid, label, productFrame, ( scatteringFunction, ) )

    def check(self, info):
        return []

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        class Parameters :

            pass

        def integrand( mu, parameters ) :

            _x = min( parameters.energy_hc * math.sqrt( 0.5 * ( 1.0 - mu ) ), parameters.scatteringFunctionDomainMax )
            kPrime_k = 1.0 / ( 1 + parameters.k * ( 1.0 - mu ) )
            integral = kPrime_k * kPrime_k * ( 1 + mu * mu + kPrime_k * parameters.k**2 * ( 1 - mu )**2 )
            integral *= parameters.scatteringFunction.evaluate( _x )
            if( parameters.energyWeight ) : integral *= kPrime_k
            if( parameters.calcualateMomentum ) : integral *= mu
            return( integral )

        def averageEnergy( parameters ) :

            parameters.energyWeight = False
            norm = nf_integrationModule.adaptiveQuadrature_GnG( 4, integrand, parameters, -1.0, 1.0, parameters.tolerance, 20 )[0]
            parameters.energyWeight = True
            energyPrime = nf_integrationModule.adaptiveQuadrature_GnG( 4, integrand, parameters, -1.0, 1.0, parameters.tolerance, 20 )[0]

            return( energyPrime / norm )

        fraction = 1.2
        energyMin, energyMax = kwargs['multiplicity'].domain( )
        N = int( math.log( energyMax / energyMin ) / math.log( fraction ) )
        energy = energyMin
        energies = []
        for i in range( N ) :
            energies.append( energy )
            energy *= fraction
        if( 1.2 * energy < energyMax ) : energies.append( energy )
        energies.append( energyMax )

        factor_electronMass = 1.0 / PQUModule.PQU( 1.0, 'me * c**2' ).getValueAs( kwargs['incidentEnergyUnit'] )
        factor_E2x = PQUModule.PQU( 1.0, '%s / hplanck / c' % kwargs['incidentEnergyUnit'] ).getValueAs( self.scatteringFunction.axes[1].unit )

        parameters = Parameters( )
        parameters.tolerance = 1e-2 * kwargs['energyAccuracy']
        parameters.scatteringFunction = self.scatteringFunction            
        parameters.scatteringFunctionDomainMax = self.scatteringFunction.domainMax

        aveEnergy = []
        aveMomenta = []
        for energy in energies :
            parameters.k = energy * factor_electronMass
            parameters.energy_hc = energy * factor_E2x

            parameters.calcualateMomentum = False
            aveEnergy.append( [ energy, energy * averageEnergy( parameters ) ] )

            parameters.calcualateMomentum = True
            aveMomenta.append( [ energy, energy * averageEnergy( parameters ) ] )

        return( [ aveEnergy ], [ aveMomenta ] )

    def processMC_cdf( self, style, tempInfo, indent ) :

        return( None )

    def processMultiGroup( self, style, tempInfo, indent ) :

        from fudge.processing.deterministic import transferMatrices as transferMatricesModule
        from fudge.processing import group as groupModule

        verbosity = tempInfo['verbosity']
        if( verbosity > 2 ) : print('%sGrouping %s' % (indent, self.moniker))

        TM_1, TM_E = transferMatricesModule.comptonScattering( style, tempInfo, self.productFrame, self.scatteringFunction,
                comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )
        return( groupModule.TMs2Form( style, tempInfo, TM_1, TM_E ) )

    def crossSection( self, domainMin, domainMax, domainUnit, tolerance = 1e-3, interpolation = standardsModule.interpolation.linlinToken ) :
        """
        Integrates the double differential cross section to get the cross section :math:`\\sigma(E)` where :math:`E` is the projectile's energy.
        """

        class Parameters :

            pass

        class tester :

            def __init__( self, parameters, relativeTolerance, absoluteTolerance ) :

                self.parameters = parameters
                self.relativeTolerance = relativeTolerance
                self.absoluteTolerance = absoluteTolerance

            def evaluateAtX( self, energy ) :

                return( integratedCrossSection( energy, self.parameters ) )

        def integrand( mu, parameters ) :

            oneMinusMu = 1.0 - mu
            _x = parameters.energy_hc * math.sqrt( 0.5 * oneMinusMu )
            kPrime_k = 1.0 / ( 1.0 + parameters.k * oneMinusMu )

            integral = kPrime_k * kPrime_k * ( 1.0 + mu * mu + kPrime_k * parameters.k**2 * oneMinusMu**2 )
            integral *= parameters.scatteringFunction.evaluate( _x )

            return( integral )

        def integratedCrossSection( energy, parameters ) :

            parameters.energy_hc = energy * parameters.factor_E2x
            parameters.k = energy * factor_electronMass
            return( parameters.constant * nf_integrationModule.adaptiveQuadrature_GnG( 4, integrand, parameters, -1.0, 1.0, parameters.tolerance, 20 )[0] )

        factor_electronMass = 1.0 / PQUModule.PQU( 1.0, 'me * c**2' ).getValueAs( domainUnit )
        factor_E2x = PQUModule.PQU( 1.0, '%s / hplanck / c' % domainUnit ).getValueAs( self.scatteringFunction.axes[1].unit )

        parameters = Parameters( )
        tolerance = min( 0.1, max( tolerance, 1e-6 ) )
        parameters.tolerance = 1e-2 * tolerance
        parameters.factor_E2x = factor_E2x
        parameters.scatteringFunction = self.scatteringFunction
        electronRadius = PQUModule.PQU( 1.0, 'e**2 / ( 4 * pi * eps0 * me * c**2 )' ).getValueAs( '1e-12 * cm' )
        parameters.constant = math.pi * electronRadius**2

        domainMax = min( domainMax, self.scatteringFunction.domainMax / factor_E2x )
        energies = [ domainMin, math.sqrt( domainMin * domainMax ), domainMax ]
        energy_crossSection = []
        for energy in energies :
            energy_crossSection.append( [ energy, integratedCrossSection( energy, parameters ) ] )
        _tester = tester( parameters, tolerance, tolerance * energy_crossSection[-1][1] )
        energy_crossSection = fudgemathModule.thickenXYList( energy_crossSection, _tester, biSectionMax = 20, interpolation = interpolation )

        return( crossSectionModule.XYs1d( data = energy_crossSection, axes = crossSectionModule.defaultAxes( domainUnit ),
                interpolation = interpolation ) )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )

        data = element[0]
        if( data.tag == XYs1d.moniker ) :
            data = XYs1d.parseXMLNode( data, xPath, linkData )
        elif( data.tag == regions1d.moniker ) :
            data = regions1d.parseXMLNode( data, xPath, linkData )
        else :
            raise TypeError( 'Invalid data "%s" for "%s"' % ( data.tag, form.moniker ) )

        _form = form( element.get( 'pid' ), element.get( 'label' ), element.get( 'productFrame' ), data )
        xPath.pop( )
        return( _form )
