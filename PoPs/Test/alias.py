#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from PoPs import database as databaseModule
from PoPs import alias as aliasModule

fIn = open( 'Answers/database.py.out' )
lines = ''.join( fIn.readlines( ) )
fIn.close( )

database = databaseModule.Database.parseXMLString(lines)

def checkFor( id ) :

    print( '    %s %s' % ( id, id in database ) )

def getID( id ) :

    try :
        item = database[id]
        try :
            pid = item.pid
        except :
            pid = None
        print( '    %s: %s --> %s' % ( id, item.id, pid ) )
    except KeyError :
        item = None
        print( '    id %s not found' % id )
    return( item )

checkFor( 'e-' )
checkFor( 'e-_anti' )
checkFor( 'e+' )
checkFor( 'e+_anti' )
checkFor( 'tau-_anti' )

#
# photon and some aliases.
#
print( '\n========== photon and some aliases ==========' )
checkFor( 'photon' )
checkFor( 'gamma' )
checkFor( 'x-ray' )
print( 'Now add gamma and x-ray' )
database.add(aliasModule.Alias('gamma', 'photon'))
database.add(aliasModule.Alias('x-ray', 'gamma'))
checkFor( 'photon' )
checkFor( 'gamma' )
checkFor( 'x-ray' )

print( '\n========== getting items ==========' )
getID( 'photon' )
gamma = getID( 'gamma' )
print( '        mass %s: spin %s: parity %s: charge %s: halflife %s' %
        ( gamma.mass[0].pqu( ), gamma.spin[0].pqu( ), gamma.parity[0].pqu( ), gamma.charge[0].pqu( ), gamma.halflife[0].value ) )
getID( 'x-ray' )

print( '\n ---- positron ----' )
database.add(aliasModule.Alias('positron', 'e-_anti'))
database.add(aliasModule.Alias('b+', 'e-_anti'))
print( database.aliases.toXML( ) )

item = getID( 'e-' )
print( '        mass %s: spin %s: parity %s: charge %s: halflife %s: isAnti %s' %
        ( item.mass[0].pqu( ), item.spin[0].pqu( ), item.parity[0].pqu( ), item.charge[0].pqu( ), item.halflife[0].value, item.isAnti ) )

item = getID( 'e-_anti' )
print( '        mass %s: spin %s: parity %s: charge %s: halflife %s: isAnti %s' %
        ( item.mass[0].pqu( ), item.spin[0].pqu( ), item.parity[0].pqu( ), item.charge[0].pqu( ), item.halflife[0].value, item.isAnti ) )

item = getID( 'positron' )
print( '        mass %s: spin %s: parity %s: charge %s: halflife %s: isAnti %s' %
        ( item.mass[0].pqu( ), item.spin[0].pqu( ), item.parity[0].pqu( ), item.charge[0].pqu( ), item.halflife[0].value, item.isAnti ) )

item = getID( 'b+' )
print( '        mass %s: spin %s: parity %s: charge %s: halflife %s: isAnti %s' %
        ( item.mass[0].pqu( ), item.spin[0].pqu( ), item.parity[0].pqu( ), item.charge[0].pqu( ), item.halflife[0].value, item.isAnti ) )

print( '\n========== meta stables ==========' )
print( 'has Am', 'Am' in database )
print( 'has Am242', 'Am242' in database )
print( 'has Am242_e0', 'Am242_e0' in database )
print( 'has Am242_e1', 'Am242_e1' in database )
print( 'has Am242_e2', 'Am242_e2' in database )
print( 'has Am242_m1', 'Am242_m1' in database )
print( 'hasAlias Am242_m1', database.hasAlias( 'Am242_m1' ) )
item = getID( 'Am242_m1' )
print( '        spin %s: parity %s: halflife %s: isAnti %s: energy %s' % \
        ( item.spin[0].pqu( ), item.parity[0].pqu( ), item.halflife[0].pqu( ), item.isAnti, item.energy.pqu( ) ) )
