#! /usr/bin/env python

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

from PoPs.quantities import quantity as quantityModule
from PoPs.quantities import mass as massModule
from PoPs.quantities import spin as spinModule
from PoPs.quantities import parity as parityModule
from PoPs.quantities import charge as chargeModule
from PoPs.quantities import halflife as halflifeModule
from PoPs.quantities import nuclearEnergyLevel as nuclearEnergyLevelModule

from PoPs.families import nucleus as nucleusModule 
from PoPs.families import nuclearLevel as nuclearLevelModule 
from PoPs.groups import isotope as isotopeModule

A = '16'
isotopeID = isotopeModule.isotopeIDFromElementIDAndA( 'O', A )
isotope = isotopeModule.suite( isotopeID, A )

data = [ [ '0', 15.99491461956, 0,       None, None, None, None ],
         [ '1',           None, 6049400, None, None, None, None ],
         [ '2',           None, 6129893, None, None, None, None ],
         [ '3',           None, 6917100, None, None, None, None ] ]

for index, mass, energy, charge, halflife, spin, parity in data :
    name = nucleusModule.nucleusNameFromNucleusNameAndIndex( isotopeID, index )
    nucleus = nucleusModule.particle( name, index )
    energy = nuclearEnergyLevelModule.double( 'base', energy, quantityModule.stringToPhysicalUnit( 'eV' ) )
    nucleus.energy.add( energy )

    name = nucleusModule.levelNameFromIsotopeNameAndIndex( isotopeID, index )
    level = nuclearLevelModule.particle( name, nucleus )

    if( mass is not None ) :
        mass = massModule.double( 'base', mass, quantityModule.stringToPhysicalUnit( 'amu' ) )
        level.mass.add( mass )

    if( charge is not None ) :
        charge = chargeModule.integer( 'base', charge, quantityModule.stringToPhysicalUnit( 'e' ) )
        level.charge.add( charge )

    if( halflife is not None ) :
        halflife = halflifeModule.double( 'base', halflife, quantityModule.stringToPhysicalUnit( 's' ) )
        level.halflife.add( halflife )

    if( spin is not None ) :
        spin = spinModule.fraction( 'base', spinModule.fraction.toValueType( spin ), quantityModule.stringToPhysicalUnit( 'hbar' ) )
        level.spin.add( spin )

    if( parity is not None ) :
        parity = parityModule.integer( 'base', parity, quantityModule.stringToPhysicalUnit( '' ) )
        level.parity.add( parity )

    isotope.add( level )

xmli1 = isotope.toXML( )
print xmli1
isotope2 = isotopeModule.suite.parseXMLStringAsClass( xmli1 )
if( xmli1 != isotope2.toXML( ) ) : raise Exception( 'Fix me.' )

print isotope.mass.pqu( )
print isotope.mass.pqu( 'MeV/c**2' )
print isotope.mass.float( 'MeV/c**2' )

print
print isotope.energy.pqu( )
print isotope.energy.pqu( 'MeV' )
print isotope.energy.float( 'MeV' )
