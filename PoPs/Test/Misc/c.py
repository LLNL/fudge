# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from xml.etree import cElementTree

from PoPs import database as databaseModule
from PoPs import alias as aliasModule
from PoPs.families import gaugeBoson as gaugeBosonModule
from PoPs.families import baryon as baryonModule
from PoPs.families import nucleus as nucleusModule
from PoPs.groups import chemicalElement as chemicalElementModule

pops = databaseModule.database( 'LLNL', '0.0.1' )

element = cElementTree.parse( 'pops.xml' )
element = element.getroot( )

def aliases( element ) :

    for child in element :
        _alias = aliasModule.particle( child.get( 'id' ), child.get( 'pid' ) )
        pops.add( _alias )

def gaugeBosons( element ) :

    for child in element :
        pops.add( gaugeBosonModule.particle.parseXMLNodeAsClass( child , [], [] ) )
    
def baryons( element ) :

    for child in element :
        pops.add( baryonModule.particle.parseXMLNodeAsClass( child , [], [] ) )
 
def chemicalElements( element ) :

    for child in element :
        pops.add( chemicalElementModule.suite.parseXMLNodeAsClass( child, [], [] ) )
 
for child in element :
    if( child.tag == 'aliases' ) :
        aliases( child )
    elif( child.tag == 'gaugeBosons' ) :
        gaugeBosons( child )
    elif( child.tag == 'baryons' ) :
        baryons( child )
    elif( child.tag == 'chemicalElements' ) :
        chemicalElements( child )

for chemicalElement in pops.chemicalElements :
    for isotope in chemicalElement :
        for level in isotope :
            nucleusName = level.id
            nucleusName = nucleusName[0].lower( ) + nucleusName[1:]
            nucleus = nucleusModule.particle( nucleusName, '0' )
            level.nucleus = nucleus

fOut = open( 'p.xml', 'w' )
fOut.write( pops.toXML( ) + '\n' )
fOut.close( )

pops2 = databaseModule.database.readFile( 'p.xml' )
fOut = open( 'p2.xml', 'w' )
fOut.write( pops2.toXML( ) + '\n' )
fOut.close( )
