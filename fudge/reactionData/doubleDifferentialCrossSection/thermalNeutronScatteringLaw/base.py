# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the base form for all thermal neutron scattering law forms and some common functions.

This module contains the following classes: 
        
    +-----------------------------------+-----------------------------------------------------------------------+
    | Class                             | Description                                                           |
    +===================================+=======================================================================+
    | Form                              | This is the base form for all thermal neutron scattering law forms.   |
    +-----------------------------------+-----------------------------------------------------------------------+

This module contains the following functions: 
        
    +-----------------------------------+-----------------------------------------------------------------------+
    | Class                             | Description                                                           |
    +===================================+=======================================================================+
    | processedToForms                  | This function takes cross section, average product energy and average |
    |                                   | product momentum data and puts them into the proper classes and       |
    |                                   | returns the created instances.                                        |
    +-----------------------------------+-----------------------------------------------------------------------+
    | energyDelta2d                     | This function creates a pointwise energy representation of elastic    |
    |                                   | scattering where :math:`P(E'|E) = \delta(E' - E)` by creating         |
    |                                   | small triangles at the minimum and maximum projectile energies.       |
    +-----------------------------------+-----------------------------------------------------------------------+
    | getNumberOfAtomsPerMolecule       | This function returns the number of atoms per molecule in the         |
    |                                   | incoherent inelastic reaction.                                        |
    +-----------------------------------+-----------------------------------------------------------------------+
"""

from fudge.reactionData import crossSection as crossSectionModule
from fudge.productData import averageProductEnergy as averageProductEnergyModule
from fudge.productData import averageProductMomentum as averageProductMomentumModule
from fudge.productData.distributions import energy as energyModule

from .. import base as baseModule

class Form( baseModule.Form ) :
    """
    This is the base form for all thermal neutron scattering law forms.
    """

    @property
    def domainMin( self ) :
        """Returns the minimum projectile energy for which data are valid."""

        return( self.rootAncestor.styles[0].projectileEnergyDomain.min )

    @property
    def domainMax( self ) :
        """Returns the maximum projectile energy for which data are valid."""

        return( self.rootAncestor.styles[0].projectileEnergyDomain.max )

    def isThermalNeutronScatteringLaw( self ) :
        """
        This method always returns True as any class that inherits this class contains thermal neutron scattering law data.
        """

        return( True )

    def energyMaximum( self ) :
        """
        This method only makes sense for incoherent inelastic scattering and is overwritten by its class.

        :returns:           This method always returns 0.
        """

        return( 0.0 )

def processedToForms( label, crossSection, averageProductEnergy, averageProductMomentum, kwargs ) :
    """
    This function takes cross section, average product energy and average product momentum data and puts
    them into the proper classes and returns the created instances.  This function is for internal use.

    :param label:                   The label of the style for the data.
    :param crossSection:            The cross section data.
    :param averageProductEnergy:    The average product energy data.
    :param averageProductMomentum:  The average product momentum data.
    :param kwargs:                  A dictionary with units for energy and momentum.

    :returns:                       An instance of :py:class:`crossSectionModule.XYs1d`, :py:class:`averageProductEnergyModule.XYs1d` 
                                    and :py:class:`averageProductMomentumModule.XYs1d`.
    """

    axes = crossSectionModule.defaultAxes( kwargs['incidentEnergyUnit'] )
    crossSection = crossSectionModule.XYs1d( data = crossSection, axes = axes, label = label )
        
    axes = averageProductEnergyModule.defaultAxes( kwargs['incidentEnergyUnit'] )
    averageProductEnergy = averageProductEnergyModule.XYs1d( data = averageProductEnergy, axes = axes, label = label )
        
    axes = averageProductMomentumModule.defaultAxes( kwargs['incidentEnergyUnit'], kwargs['momentumUnit'] )
    averageProductMomentum = averageProductMomentumModule.XYs1d( data = averageProductMomentum, axes = axes, label = label )

    return( crossSection, averageProductEnergy, averageProductMomentum )

def energyDelta2d( energyMin, energyMax, epsilon, incidentEnergyUnit ) :
    """
    This function creates a pointwise energy representation of elastic scattering where :math:`P(E'|E) = \delta(E' - E)` by
    creating small triangles at the minimum and maximum projectile energies.

    :param energyMin:               The minimum projectile energy.
    :param energyMax:               The maximum projectile energy.
    :param epsilon:                 The half width of the :math:`P(E')` triangle.
    :param incidentEnergyUnit:      The unit of the incident energy.

    :returns:                       A :py:class:`energyModule.XYs2d` instance.
    """

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
    This function returns the number of atoms per molecule in the incoherent inelastic reaction,
    Like ENDF-6, GNDS currently only stores the number of atoms per molecule in the incoherent inelastic reaction,
    but that value is also necessary for processing coherent and incoherent elastic scattering.

    :param form:    Any one of the thermalNeutronScatteringLaw forms (must have reactionSuite as an ancestor).

    :returns:       The number of principal scattering atoms per molecule
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
