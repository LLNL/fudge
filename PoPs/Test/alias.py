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

from PoPs import database as databaseModule
from PoPs import alias as aliasModule

fIn = open( 'Answers/database.py.out' )
lines = ''.join( fIn.readlines( ) )
fIn.close( )

database = databaseModule.database.parseXMLStringAsClass( lines )

def checkFor( id ) :

    print '    %s %s' % ( id, id in database )

def getID( id ) :

    try :
        item = database[id]
        try :
            pid = item.pid
        except :
            pid = None
        print '    %s: %s --> %s' % ( id, item.id, pid )
    except KeyError :
        item = None
        print '    id %s not found' % id
    return( item )

checkFor( 'e-' )
checkFor( 'e-_anti' )
checkFor( 'e+' )
checkFor( 'e+_anti' )
checkFor( 'tau-_anti' )

#
# photon and some aliases.
#
print '\n========== photon and some aliases =========='
checkFor( 'photon' )
checkFor( 'gamma' )
checkFor( 'x-ray' )
print 'Now add gamma and x-ray'
database.add( aliasModule.particle( 'gamma', 'photon' ) )
database.add( aliasModule.particle( 'x-ray', 'gamma' ) )
checkFor( 'photon' )
checkFor( 'gamma' )
checkFor( 'x-ray' )

print '\n========== getting items =========='
getID( 'photon' )
gamma = getID( 'gamma' )
print '        mass %s: spin %s: parity %s: charge %s: halflife %s' % \
        ( gamma.mass[0].pqu( ), gamma.spin[0].pqu( ), gamma.parity[0].pqu( ), gamma.charge[0].pqu( ), gamma.halflife[0].value )
getID( 'x-ray' )

print '\n ---- positron ----'
database.add( aliasModule.particle( 'positron', 'e-_anti' ) )
database.add( aliasModule.particle( 'b+', 'e-_anti' ) )
print database.aliases.toXML( )

item = getID( 'e-' )
print '        mass %s: spin %s: parity %s: charge %s: halflife %s: isAnti %s' % \
        ( item.mass[0].pqu( ), item.spin[0].pqu( ), item.parity[0].pqu( ), item.charge[0].pqu( ), item.halflife[0].value, item.isAnti )

item = getID( 'e-_anti' )
print '        mass %s: spin %s: parity %s: charge %s: halflife %s: isAnti %s' % \
        ( item.mass[0].pqu( ), item.spin[0].pqu( ), item.parity[0].pqu( ), item.charge[0].pqu( ), item.halflife[0].value, item.isAnti )

item = getID( 'positron' )
print '        mass %s: spin %s: parity %s: charge %s: halflife %s: isAnti %s' % \
        ( item.mass[0].pqu( ), item.spin[0].pqu( ), item.parity[0].pqu( ), item.charge[0].pqu( ), item.halflife[0].value, item.isAnti )

item = getID( 'b+' )
print '        mass %s: spin %s: parity %s: charge %s: halflife %s: isAnti %s' % \
        ( item.mass[0].pqu( ), item.spin[0].pqu( ), item.parity[0].pqu( ), item.charge[0].pqu( ), item.halflife[0].value, item.isAnti )

print '\n========== meta stables =========='
print 'has Am', 'Am' in database
print 'has Am242', 'Am242' in database
print 'has Am242_e0', 'Am242_e0' in database
print 'has Am242_e1', 'Am242_e1' in database
print 'has Am242_e2', 'Am242_e2' in database
print 'has Am242_m1', 'Am242_m1' in database
print 'hasAlias Am242_m1', database.hasAlias( 'Am242_m1' )
item = getID( 'Am242_m1' )
print '        spin %s: parity %s: halflife %s: isAnti %s: energy %s' % \
        ( item.spin[0].pqu( ), item.parity[0].pqu( ), item.halflife[0].pqu( ), item.isAnti, item.energy.pqu( ) )
