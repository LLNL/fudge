# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
# <<END-copyright>>

"""
This module contains the average energy classes.
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

        physicalQuantityModule.physicalQuantity.__init__( self, value, unit, self._label )

    # overrides required since __init__ arguments differ
    def copy( self ) :

        instance = self.__class__( self.value, self.unit )
        instance.uncertainty = self.uncertainty
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

class proton( suiteModule.suite ) :

    _label = 'proton'

class neutrino( suiteModule.suite ) :

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
