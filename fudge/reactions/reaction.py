# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the reaction class.
"""

from xData import enums as xDataEnumsModule

from PoPs import IDs as popsIDsModule
from PoPs import database as PoPsDatabaseModule
from PoPs import specialNuclearParticleID as specialNuclearParticleID_PoPsModule

from fudge.core.utilities import fudgeExceptions

from .. import enums as enumsModule
from .. import outputChannel as outputChannelModule
from .. import reactionProducts as reactionProductsModule
from ..reactionData import availableEnergy as availableEnergyModule
from ..reactionData import availableMomentum as availableMomentumModule
from ..reactionData.doubleDifferentialCrossSection.chargedParticleElastic import CoulombPlusNuclearElastic as CoulombPlusNuclearElasticModule
from ..reactionData.doubleDifferentialCrossSection.thermalNeutronScatteringLaw import \
    base as thermalNeutronScatteringLawBaseModule

from fudge.productData.distributions import angular as angularModule

from . import base as baseModule

class Reaction(baseModule.Base_reaction2):
    """This is the class for a normal gnds reaction."""

    moniker = 'reaction'
    ancestryMembers = baseModule.Base_reaction2.ancestryMembers + ( 'availableEnergy', 'availableMomentum' )

    def __init__(self, label, genre, ENDF_MT, fissionGenre=enumsModule.FissionGenre.none):
        """
        Creates a new reaction object. Reaction is two-body or uncorrelated-body, depending on
        the outputChannel genre. This class is only meant to be used for 'distinct' reactions (distinct reactions
        do not overlap any other reactions; altogether they sum up to total).
        To store a sum over these distinct reactions, use the product class.
        """

        baseModule.Base_reaction2.__init__(self, label, genre, ENDF_MT, fissionGenre=fissionGenre)

        self.availableEnergy = availableEnergyModule.Component( )
        self.availableEnergy.setAncestor( self )

        self.availableMomentum = availableMomentumModule.Component( )
        self.availableMomentum.setAncestor( self )

        self.__reactionProducts = None

    def __eq__(self, other):
        if (not baseModule.isGNDSReaction( other )): return False
        selfParent, otherParent = self.getReactionSuite( ), other.getReactionSuite( )
        return( ( selfParent.projectile == otherParent.projectile ) and ( selfParent.target == otherParent.target ) 
            and ( self.outputChannel == other.outputChannel ) )
    
    def getThreshold(self, unit):

        Q = self.getQ(unit = unit, final = False)
        if Q >= 0.: return 0.
        reactionSuite = self.getReactionSuite()
        projectile = reactionSuite.PoPs[reactionSuite.projectile]
        projectileMass = projectile.getMass('amu')
        targetID = reactionSuite.target
        if targetID in reactionSuite.PoPs.aliases:
            targetID = reactionSuite.PoPs[targetID].pid
        target = reactionSuite.PoPs[targetID]
        targetMass = target.getMass('amu')
        return -Q * (1. + projectileMass / targetMass)

    def isBasicReaction( self ) :

        return( True )

    def isCompleteReaction( self ) :

        return( True )

    def processCoulombPlusNuclearMuCutoff( self, style, excludeRutherfordScattering = False ) :

        for doubleDifferentialCrossSection in self.doubleDifferentialCrossSection :
            if( isinstance( doubleDifferentialCrossSection, CoulombPlusNuclearElasticModule.Form ) ) :
                crossSection, angular = doubleDifferentialCrossSection.processCoulombPlusNuclearMuCutoff( style )

                self.crossSection.add( crossSection )

                product = self.outputChannel[0]
                if( product.pid != doubleDifferentialCrossSection.pid ) : raise Exception( 'First product out not one described in doubleDifferentialCrossSection data.' )
                product.distribution.add( angularModule.TwoBody( style.label, xDataEnumsModule.Frame.centerOfMass, angular ) )

                residual = self.outputChannel[1]
                residual.distribution.add( angularModule.TwoBody( style.label, xDataEnumsModule.Frame.centerOfMass, angularModule.Recoil( angular ) ) )

                if( excludeRutherfordScattering ) :
                    crossSection, angular = doubleDifferentialCrossSection.processCoulombPlusNuclearMuCutoff( style, excludeRutherfordScattering = True )
                    if( crossSection is not None ) : return( crossSection, angular )

        return( None )

    def processThermalNeutronScatteringLaw( self, style, kwargs ) :

        for doubleDifferentialCrossSection in self.doubleDifferentialCrossSection :
            if( isinstance( doubleDifferentialCrossSection, thermalNeutronScatteringLawBaseModule.Form ) ) :

                product = self.outputChannel[0]
                if( product.pid != doubleDifferentialCrossSection.pid ) : raise Exception( 'First product out not one described in doubleDifferentialCrossSection data.' )

                crossSection, averageProductEnergy, averageProductMomentum, distribution = doubleDifferentialCrossSection.processThermalNeutronScatteringLaw( style, kwargs )

                self.crossSection.add( crossSection )

                product.distribution.add( distribution )
                product.averageProductEnergy.add( averageProductEnergy )
                product.averageProductMomentum.add( averageProductMomentum )

                axes = availableEnergyModule.defaultAxes( self.domainUnit )
                data = [ [ crossSection.domainMin, crossSection.domainMin ], [ crossSection.domainMax, crossSection.domainMax ] ]
                self.availableEnergy.add( availableEnergyModule.XYs1d( data = data, axes = axes, label = style.label ) )

                axes = availableMomentumModule.defaultAxes( self.domainUnit, kwargs['momentumUnit'] )
                data = data = [ [ crossSection.domainMin, 0.0 ], [ crossSection.domainMax, 0.0 ] ]
                self.availableMomentum.add( availableMomentumModule.XYs1d( data = data, axes = axes, label = style.label ) )

                break

    def reactionProducts( self ) :

        if( self.__reactionProducts is None ) : self.__reactionProducts = self.outputChannel.reactionProducts( reactionProductsModule.ReactionProducts( ) )
        if self.fissionGenre != enumsModule.FissionGenre.none:
            self.__reactionProducts[str(self.fissionGenre)+'Fission'] = reactionProductsModule.ReactionProduct( reactionProductsModule.Category.process, 0 )
        return( self.__reactionProducts )

    def specialReactionLabel(self, outputMode, mode=specialNuclearParticleID_PoPsModule.Mode.familiar):
        """
        Returns the special label for *self*. The special label is a unique label for the reaction that is meant to 
        be the same for the same reaction type.

        :param outputMode:      An value from the reactionProducts.Mode enum.
        :param mode:            A value from the special particle id Mode enum.

        :return:                The special label for *self*.
        """

        pops = self.findAttributeInAncestry(PoPsDatabaseModule.Database.moniker)
        if self.label == str(enumsModule.Genre.sumOfRemainingOutputChannels):
            return self.label
        return self.reactionProducts().specialReactionLabel(outputMode, pops, mode=mode)

    def thermalNeutronScatteringLawTemperatures( self, temperatures ) :

        for doubleDifferentialCrossSection in self.doubleDifferentialCrossSection :
            if( isinstance( doubleDifferentialCrossSection, thermalNeutronScatteringLawBaseModule.Form ) ) :
                doubleDifferentialCrossSection.temperatures( temperatures )
