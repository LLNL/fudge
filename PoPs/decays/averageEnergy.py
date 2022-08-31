# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module defines classes for storing the average outgoing energy of various types of particles
emitted during a decay.

Average energies are not as informative as detailed spectra, but for decay modes with very small probabilities
they may be the only available information.
These classes follow the example of the ENDF decay sub-library, providing average energies for a few general classes
of decay products, such as 'electroMagneticRadiation', 'lightParticles', 'heavyParticles', etc.

These classes are thin wrappers around the xData.physicalQuantity.PhysicalQuantity class.
"""

import abc

from xData import physicalQuantity as physicalQuantityModule
from xData.uncertainty.physicalQuantity import uncertainty as uncertaintyModule

from .. import suite as suiteModule

#
# FIXME Need a physicalQuantity class with keyName.
#
class AverageEnergy( physicalQuantityModule.PhysicalQuantity ) :

    moniker = 'averageEnergy'
    keyName = 'label'

    def __init__( self, value, unit ) :
        """See documentation in xData.physicalQuantity.PhysicalQuantity.__init__"""

        physicalQuantityModule.PhysicalQuantity.__init__( self, value, unit, self._label )

    # overrides required since __init__ arguments differ
    def copy( self ) :

        instance = self.__class__( self.value, self.unit )
        
        if self.uncertainty is not None: instance.uncertainty = self.uncertainty.copy( )
        return( instance )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )

        instance = cls( float( element.get( 'value' ) ), element.get( 'unit' ) )
        for child in element :
            if( child.tag == uncertaintyModule.Uncertainty.moniker ) :
                instance.uncertainty = uncertaintyModule.Uncertainty.parseNodeUsingClass(child, xPath, linkData, **kwargs)

        xPath.pop( )
        return( instance )

class LightParticles( AverageEnergy ) :

    _label = 'lightParticles'

class ElectroMagneticRadiation( AverageEnergy ) :

    _label = 'electroMagneticRadiation'

class HeavyParticles( AverageEnergy ) :

    _label = 'heavyParticles'

class BetaMinus( AverageEnergy ) :

    _label = 'betaMinus'

class BetaPlus( AverageEnergy ) :

    _label = 'betaPlus'

class AugerElectron( AverageEnergy ) :

    _label = 'AugerElectron'

class ConversionElectron( AverageEnergy ) :

    _label = 'conversionElectron'

class Gamma( AverageEnergy ) :

    _label = 'gamma'

class XRay( AverageEnergy ) :

    _label = 'xRay'

class InternalBremsstrahlung( AverageEnergy ) :

    _label = 'internalBremsstrahlung'

class Annihilation( AverageEnergy ) :

    _label = 'annihilation'

class Alpha( AverageEnergy ) :

    _label = 'alpha'

class Recoil( AverageEnergy ) :

    _label = 'recoil'

class SpontaneousFission( AverageEnergy ) :

    _label = 'spontaneousFission'

class FissionNeutrons( AverageEnergy ) :

    _label = 'fissionNeutrons'

class Proton( AverageEnergy ) :

    _label = 'proton'

class Neutrino( AverageEnergy ) :

    _label = 'neutrino'

class AverageEnergies( suiteModule.Suite ) :

    moniker = 'averageEnergies'

    def __init__( self ) :

        suiteModule.Suite.__init__( self, ( AverageEnergy, ) )

    def parseNode(self, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )

        for child in element:
            label = child.get('label')
            for subclass in (
                LightParticles, ElectroMagneticRadiation, HeavyParticles, BetaMinus, BetaPlus, AugerElectron,
                ConversionElectron, Gamma, XRay, InternalBremsstrahlung, Annihilation, Alpha, Recoil,
                SpontaneousFission, FissionNeutrons, Proton, Neutrino
            ):
                if label == subclass._label:
                    self.add( subclass.parseNodeUsingClass(child, xPath, linkData, **kwargs))
        xPath.pop( )
