# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import math
import numpy

from fudge.core.math import fudgemath
from fudge.productData.distributions import angular as angularModule
from fudge.productData.distributions import uncorrelated as uncorrelatedModule

from . import base as baseModule

def process( self, label, energyMin, energyMax, temperature, kwargs ) :

    epsilon = min( 1e-3, max( 1e-15, kwargs['epsilon'] ) )
    incidentEnergyUnit = kwargs['incidentEnergyUnit']

    atomsPerMolecule = baseModule.getNumberOfAtomsPerMolecule(self)

    sigma_b = self.boundAtomCrossSection.getValueAs( 'b' )
    if( temperature < self.DebyeWallerIntegral.function1d.domainMin ) : temperature = self.DebyeWallerIntegral.function1d.domainMin
    if( temperature > self.DebyeWallerIntegral.function1d.domainMax ) : temperature = self.DebyeWallerIntegral.function1d.domainMax
    debyeWallerIntegralPrime = self.DebyeWallerIntegral.function1d.evaluate( temperature )
    momentumConstant = math.sqrt( 2.0 * kwargs['productMass'] )

    parameters = { 'debyeWallerIntegralPrime' : debyeWallerIntegralPrime, 'sigma_b' : sigma_b, 'momentumConstant' : momentumConstant }

    # initial energy grid: 4 points per decade
    bottom = numpy.floor(numpy.log10(energyMin))
    top = numpy.floor(numpy.log10(energyMax))
    decades = int(top - bottom)
    initialGrid = [E for E in numpy.logspace(bottom, top, decades*4 + 1) if E >= energyMin]
    if energyMax not in initialGrid: initialGrid.append(energyMax)

    crossSection = [ [ energy, crossSectionAtEnergy( energy, parameters ) ] for energy in initialGrid ]
    crossSection = fudgemath.thickenXYList( crossSection, Evaluator( crossSectionAtEnergy, parameters ), biSectionMax = 12 )

    averageProductEnergy = [ [ energy, averageProductEnergyAtEnergy( energy, parameters ) ] for energy in [ initialGrid[0], initialGrid[-1] ] ]

    averageProductMomentum = [ [ energy, averageProductMomentumAtEnergy( energy, parameters ) ] for energy in initialGrid ]
    averageProductMomentum = fudgemath.thickenXYList( averageProductMomentum, Evaluator( averageProductMomentumAtEnergy, parameters ), biSectionMax = 12 )
 
    crossSection, averageProductEnergy, averageProductMomentum = baseModule.processedToForms( label, crossSection, averageProductEnergy, averageProductMomentum, kwargs )
    crossSection /= atomsPerMolecule

    axes = angularModule.defaultAxes( incidentEnergyUnit )
    angular2d = angularModule.XYs2d( axes = axes )

    angularData = []
    for energy, sigma in crossSection :
        parameters['energy'] = energy
        PmuForEachE = [ [ mu, incoherentElasticPOfMuGivenE( mu, parameters ) ] for mu in [ -1., 1. ] ]
        PmuForEachE = fudgemath.thickenXYList( PmuForEachE, Evaluator( incoherentElasticPOfMuGivenE, parameters ), biSectionMax = 12 )
        PmuForEachE = angularModule.XYs1d( data = PmuForEachE, outerDomainValue = energy, axes = axes )
        PmuForEachE = PmuForEachE.normalize( )
        angular2d.append( PmuForEachE )
    angular2d = uncorrelatedModule.AngularSubform( angular2d )

    energy2d = baseModule.energyDelta2d( crossSection[0][0], crossSection[-1][0], epsilon, incidentEnergyUnit )
    energy2d = uncorrelatedModule.EnergySubform( energy2d )

    distribution = uncorrelatedModule.Form( label, self.productFrame, angular2d, energy2d )

    return ( crossSection, averageProductEnergy, averageProductMomentum, distribution )
 
def crossSectionAtEnergy( energy, parameters ) :
    """
    Calculates the incoherent elastic thermal scattering law cross section as:

    sigma(E) = sigma_b / 2 * ( 1 - exp[ -4 E W' ] ) / ( 2 E W' )

    where E in the neutron energy and W' is the DebyeWaller integral divided by the atomic mass.
    """

    debyeWallerIntegralPrime = parameters['debyeWallerIntegralPrime']
    sigma_b = parameters['sigma_b']
    x = 4. * energy * debyeWallerIntegralPrime
    return( sigma_b * ( 1. - math.exp( -x ) ) / x )

def averageProductEnergyAtEnergy( energy, parameters ) :
    """
    Calculates the incoherent elastic thermal scattering law average product energy as E'(E) = E.
    """

    return( energy )

def averageProductMomentumAtEnergy( energy, parameters ) :
    """
    Calculates the incoherent elastic thermal scattering law average product momentum as:

    p(E) = sqrt( 2 * m_n * E ) * sigma_b / 2 * ( 2 E W' - 1  + ( 1 + 2 E W' ) * exp[ -4 E W' ] ) / ( 2 E W' )**2 /
                    ( sigma_b / ( 2 * ( 1 - exp[ -4 E W' ] ) / ( 2 E W' ) ) )
         = sqrt( 2 * m_n * E ) * ( 2 E W' - 1  + ( 1 + 2 E W' ) * exp[ -4 E W' ] ) / ( 2 E W' ) / ( 1 - exp[ -4 E W' ] )
    
    where E in the neutron energy and W' is the DebyeWaller integral divided by the atomic mass.
    """

    debyeWallerIntegralPrime = parameters['debyeWallerIntegralPrime']
    TwoEW = 2.0 * energy * debyeWallerIntegralPrime
    exp_4EW = math.exp( -2.0 * TwoEW )
    angularFactor = ( TwoEW - 1.0 + ( 1.0 + TwoEW ) * exp_4EW ) / ( TwoEW * ( 1.0 - exp_4EW ) )
    
    return( parameters['momentumConstant'] * math.sqrt( energy ) * angularFactor )

def incoherentElasticPOfMuGivenE( mu, parameters) :

    energy = parameters['energy']
    debyeWallerIntegralPrime = parameters['debyeWallerIntegralPrime']
    return( math.exp( -2. * energy * debyeWallerIntegralPrime * ( 1. - mu ) ) )

class Evaluator :

    def __init__( self, func, parameters, relativeTolerance = 1e-3, absoluteTolerance = 1e-10 ) :

        self.func = func
        self.parameters = parameters
        self.relativeTolerance = relativeTolerance
        self.absoluteTolerance = absoluteTolerance

    def evaluateAtX( self, energy ) :

        return( self.func( energy, self.parameters ) )
