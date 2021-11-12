# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the reaction class.
"""

from xData import standards as standardsModule

from fudge.core.utilities import fudgeExceptions

from fudge import reactionProducts as reactionProductsModule
from fudge.reactionData import availableEnergy as availableEnergyModule
from fudge.reactionData import availableMomentum as availableMomentumModule
from fudge.reactionData.doubleDifferentialCrossSection.chargedParticleElastic import CoulombPlusNuclearElastic as CoulombPlusNuclearElasticModule
from ..reactionData.doubleDifferentialCrossSection.thermalNeutronScatteringLaw import \
    base as thermalNeutronScatteringLawBaseModule

from fudge.productData.distributions import angular as angularModule

from . import base as baseModule

__metaclass__ = type

class reaction( baseModule.base_reaction ) :
    """This is the class for a normal gnds reaction."""

    moniker = 'reaction'
    ancestryMembers = baseModule.base_reaction.ancestryMembers + ( 'availableEnergy', 'availableMomentum' )

    def __init__( self, genre, ENDF_MT, fissionGenre = None, documentation = None, label = None ) :
        """
        Creates a new reaction object. Reaction is two-body or uncorrelated-body, depending on
        the outputChannel genre. This class is only meant to be used for 'distinct' reactions (distinct reactions
        do not overlap any other reactions; altogether they sum up to total).
        To store a sum over these distinct reactions, use the product class.
        """

        baseModule.base_reaction.__init__( self, genre, ENDF_MT, fissionGenre, documentation, label = label )

        self.availableEnergy = availableEnergyModule.component( )
        self.availableEnergy.setAncestor( self )

        self.availableMomentum = availableMomentumModule.component( )
        self.availableMomentum.setAncestor( self )

        self.__reactionProducts = None

    def __eq__(self, other):
        if (not baseModule.isGNDSReaction( other )): return False
        selfParent, otherParent = self.getReactionSuite( ), other.getReactionSuite( )
        return( ( selfParent.projectile == otherParent.projectile ) and ( selfParent.target == otherParent.target ) 
            and ( self.outputChannel == other.outputChannel ) )
    
    def __cmp__( self, other ) :
        """Test if self is <, == or > other."""

        if( not baseModule.isGNDSReaction( other ) ) : raise fudgeExceptions.FUDGE_Exception( "Other not an reaction object." )
        selfParent, otherParent = self.getReactionSuite( ), other.getReactionSuite( )
        if( selfParent.projectile < otherParent.projectile ) : return( -1 )
        if( selfParent.projectile > otherParent.projectile ) : return(  1 )
        if( selfParent.target < otherParent.target ) : return( -1 )
        if( selfParent.target > otherParent.target ) : return(  1 )
        if( self.outputChannel < other.outputChannel ) : return( -1 )
        if( self.outputChannel > other.outputChannel ) : return(  1 )
        return( 0 )

    def getThreshold( self, unit ) :

        Q = self.getQ( unit = unit, final = False )
        if( Q >= 0. ) : return( 0. )
        reactionSuite = self.getReactionSuite( )
        projectile = reactionSuite.PoPs[reactionSuite.projectile]
        projectileMass = projectile.mass[0].float( 'amu' )
        targetID = reactionSuite.target
        if( targetID in reactionSuite.PoPs.aliases ) : targetID = reactionSuite.PoPs[targetID].pid
        target = reactionSuite.PoPs[targetID]
        try:
            targetMass = target.mass[0].float( 'amu' )
        except IndexError:
            targetMass = target.getMass('amu')
        return( -Q * ( 1. + projectileMass / targetMass ) )

    def isBasicReaction( self ) :

        return( True )

    def isCompleteReaction( self ) :

        return( True )

    def processCoulombPlusNuclearMuCutoff( self, style, excludeRutherfordScattering = False ) :

        for doubleDifferentialCrossSection in self.doubleDifferentialCrossSection :
            if( isinstance( doubleDifferentialCrossSection, CoulombPlusNuclearElasticModule.form ) ) :
                crossSection, angular = doubleDifferentialCrossSection.processCoulombPlusNuclearMuCutoff( style )

                self.crossSection.add( crossSection )

                product = self.outputChannel[0]
                if( product.pid != doubleDifferentialCrossSection.pid ) : raise Exception( 'First product out not one described in doubleDifferentialCrossSection data.' )
                product.distribution.add( angularModule.twoBodyForm( style.label, standardsModule.frames.centerOfMassToken, angular ) )

                residual = self.outputChannel[1]
                residual.distribution.add( angularModule.twoBodyForm( style.label, standardsModule.frames.centerOfMassToken, angularModule.recoil( angular ) ) )

                if( excludeRutherfordScattering ) :
                    crossSection, angular = doubleDifferentialCrossSection.processCoulombPlusNuclearMuCutoff( style, excludeRutherfordScattering = True )
                    if( crossSection is not None ) : return( crossSection, angular )

                return( None )

    def processThermalNeutronScatteringLaw( self, style, kwargs ) :

        for doubleDifferentialCrossSection in self.doubleDifferentialCrossSection :
            if( isinstance( doubleDifferentialCrossSection, thermalNeutronScatteringLawBaseModule.form ) ) :

                product = self.outputChannel[0]
                if( product.pid != doubleDifferentialCrossSection.pid ) : raise Exception( 'First product out not one described in doubleDifferentialCrossSection data.' )

                crossSection, averageProductEnergy, averageProductMomentum, distribution = doubleDifferentialCrossSection.processThermalNeutronScatteringLaw( style, kwargs )

                self.crossSection.add( crossSection )

                product.distribution.add( distribution )
                product.energyDeposition.add( averageProductEnergy )
                product.momentumDeposition.add( averageProductMomentum )

                axes = availableEnergyModule.defaultAxes( self.domainUnit )
                data = [ [ crossSection.domainMin, crossSection.domainMin ], [ crossSection.domainMax, crossSection.domainMax ] ]
                self.availableEnergy.add( availableEnergyModule.XYs1d( data = data, axes = axes, label = style.label ) )

                axes = availableMomentumModule.defaultAxes( self.domainUnit, kwargs['momentumUnit'] )
                data = data = [ [ crossSection.domainMin, 0.0 ], [ crossSection.domainMax, 0.0 ] ]
                self.availableMomentum.add( availableMomentumModule.XYs1d( data = data, axes = axes, label = style.label ) )

                break

    def reactionProducts( self ) :

        if( self.__reactionProducts is None ) : self.__reactionProducts = self.outputChannel.reactionProducts( reactionProductsModule.ReactionProducts( ) )
        if self.fissionGenre is not None:
            self.__reactionProducts[self.fissionGenre] = reactionProductsModule.ReactionProduct( reactionProductsModule.Category.process, 0 )
        return( self.__reactionProducts )

    def thermalNeutronScatteringLawTemperatures( self, temperatures ) :

        for doubleDifferentialCrossSection in self.doubleDifferentialCrossSection :
            if( isinstance( doubleDifferentialCrossSection, thermalNeutronScatteringLawBaseModule.form ) ) :
                doubleDifferentialCrossSection.temperatures( temperatures )
