# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
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

These classes are thin wrappers around the xData.physicalQuantity.physicalQuantity class.
"""

import abc

from xData import physicalQuantity as physicalQuantityModule
from xData.uncertainty.physicalQuantity import uncertainty as uncertaintyModule

from .. import suite as suiteModule

#
# FIXME Need a physicalQuantity class with keyName.
#
class averageEnergy( physicalQuantityModule.physicalQuantity ) :

    moniker = 'averageEnergy'
    keyName = 'label'

    def __init__( self, value, unit ) :
        """See documentation in xData.physicalQuantity.physicalQuantity.__init__"""

        physicalQuantityModule.physicalQuantity.__init__( self, value, unit, self._label )

    # overrides required since __init__ arguments differ
    def copy( self ) :

        instance = self.__class__( self.value, self.unit )
        instance.uncertainty = self.uncertainty.copy( )
        return( instance )

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):

        xPath.append( element.tag )

        instance = cls( float( element.get( 'value' ) ), element.get( 'unit' ) )
        for child in element :
            if( child.tag == uncertaintyModule.uncertainty.moniker ) :
                instance.uncertainty = uncertaintyModule.uncertainty.parseXMLNodeAsClass( child, xPath, linkData )

        xPath.pop( )
        return( instance )

class lightParticles( averageEnergy ) :

    _label = 'lightParticles'

class electroMagneticRadiation( averageEnergy ) :

    _label = 'electroMagneticRadiation'

class heavyParticles( averageEnergy ) :

    _label = 'heavyParticles'

class betaMinus( averageEnergy ) :

    _label = 'betaMinus'

class betaPlus( averageEnergy ) :

    _label = 'betaPlus'

class AugerElectron( averageEnergy ) :

    _label = 'AugerElectron'

class conversionElectron( averageEnergy ) :

    _label = 'conversionElectron'

class gamma( averageEnergy ) :

    _label = 'gamma'

class xRay( averageEnergy ) :

    _label = 'xRay'

class internalBremsstrahlung( averageEnergy ) :

    _label = 'internalBremsstrahlung'

class annihilation( averageEnergy ) :

    _label = 'annihilation'

class alpha( averageEnergy ) :

    _label = 'alpha'

class recoil( averageEnergy ) :

    _label = 'recoil'

class spontaneousFission( averageEnergy ) :

    _label = 'spontaneousFission'

class fissionNeutrons( averageEnergy ) :

    _label = 'fissionNeutrons'

class proton( averageEnergy ) :

    _label = 'proton'

class neutrino( averageEnergy ) :

    _label = 'neutrino'

class averageEnergies( suiteModule.suite ) :

    moniker = 'averageEnergies'

    def __init__( self ) :

        suiteModule.suite.__init__( self, ( averageEnergy, ) )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( element.tag )

        for child in element:
            label = child.get('label')
            for subclass in (
                lightParticles, electroMagneticRadiation, heavyParticles, betaMinus, betaPlus, AugerElectron,
                conversionElectron, gamma, xRay, internalBremsstrahlung, annihilation, alpha, recoil,
                spontaneousFission, fissionNeutrons, proton, neutrino
            ):
                if label == subclass._label:
                    self.add( subclass.parseXMLNode( child, xPath, linkData ) )
        xPath.pop( )
