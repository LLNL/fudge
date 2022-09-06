# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Base form for all thermal neutron scattering law forms.
"""

from fudge.reactionData import crossSection as crossSectionModule
from fudge.productData import averageProductEnergy as averageProductEnergyModule
from fudge.productData import averageProductMomentum as averageProductMomentumModule
from fudge.productData.distributions import energy as energyModule

from .. import base as baseModule

class Form( baseModule.Form ) :

    @property
    def domainMin( self ) :
        "Returns the minimum projectile energy for which data are valid."

        return( self.rootAncestor.styles[0].projectileEnergyDomain.min )

    @property
    def domainMax( self ) :
        "Returns the maximum projectile energy for which data are valid."

        return( self.rootAncestor.styles[0].projectileEnergyDomain.max )

    def isThermalNeutronScatteringLaw( self ) :

        return( True )

    def energyMaximum( self ) :

        return( 0.0 )

def processedToForms( label, crossSection, averageProductEnergy, averageProductMomentum, kwargs ) :

    axes = crossSectionModule.defaultAxes( kwargs['incidentEnergyUnit'] )
    crossSection = crossSectionModule.XYs1d( data = crossSection, axes = axes, label = label )
        
    axes = averageProductEnergyModule.defaultAxes( kwargs['incidentEnergyUnit'] )
    averageProductEnergy = averageProductEnergyModule.XYs1d( data = averageProductEnergy, axes = axes, label = label )
        
    axes = averageProductMomentumModule.defaultAxes( kwargs['incidentEnergyUnit'], kwargs['momentumUnit'] )
    averageProductMomentum = averageProductMomentumModule.XYs1d( data = averageProductMomentum, axes = axes, label = label )

    return( crossSection, averageProductEnergy, averageProductMomentum )

def energyDelta2d( energyMin, energyMax, epsilon, incidentEnergyUnit ) :

    axes = energyModule.defaultAxes( incidentEnergyUnit )
    energy2d = energyModule.XYs2d( axes = axes )

    energy = energyMin
    data = [ [ energy * ( 1.0 - epsilon ), 0.0 ], [ energy, 1.0 / epsilon ], [ energy * ( 1.0 + epsilon ), 0.0 ] ]
    energy2d.append( energyModule.XYs1d( data = data, outerDomainValue = energy, axes = axes ).normalize( ) )

    energy = energyMax
    data = [ [ energy * ( 1.0 - epsilon ), 0.0 ], [ energy, 1.0 / epsilon ], [ energy * ( 1.0 + epsilon ), 0.0 ] ]
    energy2d.append( energyModule.XYs1d( data = data, outerDomainValue = energy, axes = axes ).normalize( ) )

    return( energy2d )

def getNumberOfAtomsPerMolecule( form ):
    """
    Like ENDF-6,  GNDS only stores the number of atoms per molecule in the incoherent inelastic section,
    but that value is also necessary for processing elastic. This routine grabs that value.
    @param form: any thermalNeutronScatteringLaw form (must have reactionSuite as an ancestor)
    @return: number of principal scattering atoms per molecule
    """
    from fudge import reactionSuite as reactionSuiteModule
    from fudge.reactionData.doubleDifferentialCrossSection.thermalNeutronScatteringLaw import incoherentInelastic as incoherentInelasticModule
    RS = form.findClassInAncestry(reactionSuiteModule.ReactionSuite)
    IIform = None
    for reaction in RS.reactions:
        if reaction.doubleDifferentialCrossSection:
            doubleDifferentialForm = reaction.doubleDifferentialCrossSection.evaluated
            if isinstance(doubleDifferentialForm, incoherentInelasticModule.Form):
                IIform = doubleDifferentialForm
                break
    if not IIform:
        print("WARNING: could not find incoherent inelastic reaction, so assuming number of atoms/molecule = 1")
        return 1
    return IIform.scatteringAtoms[0].numberPerMolecule
