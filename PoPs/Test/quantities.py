#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from pqu import PQU as PQUModule

from PoPs.quantities import quantity as quantityModule
from PoPs.quantities import mass as massModule
from PoPs.quantities import spin as spinModule
from PoPs.quantities import parity as parityModule
from PoPs.quantities import charge as chargeModule
from PoPs.quantities import halflife as halflifeModule

@classmethod
def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

    suite = cls()
    suite.parseNode(node, xPath, linkData, **kwargs)

    return suite
quantityModule.Suite.parseNodeUsingClass = parseNodeUsingClass

print( '============== MASS ==============' )
m1 = massModule.Double( 'atomic', 2.12, quantityModule.stringToPhysicalUnit( 'amu' ) )
xmlm1 = m1.toXML( )
print( xmlm1 )
m2 = massModule.Double.parseXMLString(xmlm1)
if( xmlm1 != m2.toXML( ) ) : raise Exception( 'Fix me.' )

print()
suite = massModule.Suite( )
suite.add( m1 )

m2 = massModule.Double( 'nuclear', 2, quantityModule.stringToPhysicalUnit( 'amu' ) )
suite.add( m2 )

print( suite.toXML( ) )
suite2 = massModule.Suite.parseXMLString(suite.toXML())
if( suite2.toXML( ) != suite.toXML( ) ) : raise Exception( 'Fix me' )

print( '\n============== CHARGE ==============' )
c1 = chargeModule.Integer( 'nucleus', 3, quantityModule.stringToPhysicalUnit( 'e' ) )
xmlc1 = c1.toXML( )
print( xmlc1 )
c2 = chargeModule.Integer.parseXMLString(xmlc1)
if( xmlc1 != c2.toXML( ) ) : raise Exception( 'Fix me.' )

print()
suite = chargeModule.Suite( )
suite.add( c1 )
print( suite.toXML( ) )

suite2 = suite.parseXMLString(suite.toXML())
if( suite2.toXML( ) != suite.toXML( ) ) : raise Exception( 'Fix me' )

print( '\n============== HALFLIFE ==============' )
h1 = halflifeModule.Double( 'nucleus', 3.14e6, quantityModule.stringToPhysicalUnit( 'd' ) )
xmlh1 = h1.toXML( )
print( xmlh1 )
h2 = halflifeModule.Double.parseXMLString(xmlh1)
if( xmlh1 != h2.toXML( ) ) : raise Exception( 'Fix me.' )

print()
suite = halflifeModule.Suite( )
suite.add( h1 )
print( suite.toXML( ) )

suite2 = suite.parseXMLString(suite.toXML())
if( suite2.toXML( ) != suite.toXML( ) ) : raise Exception( 'Fix me' )

print()
print( h1.pqu( ) )
print( suite[0].pqu( ) )
print( suite[0].pqu( 's' ) )
print( PQUModule.floatToShortestString( suite[0].pqu( ).getValueAs( 's' ), 12 ) )

print()
suite = halflifeModule.Suite( )
h2 = halflifeModule.String( 'nucleus', 'stable', quantityModule.stringToPhysicalUnit( 's' ) )
suite.add( h2 )
print( suite.toXML( ) )

suite2 = suite.parseXMLString(suite.toXML())
if( suite2.toXML( ) != suite.toXML( ) ) : raise Exception( 'Fix me' )

print( '\n============== SPIN ==============' )
fraction = quantityModule.Fraction.toValueType( "5/2" )
s1 = spinModule.Fraction( 'nucleus', fraction, quantityModule.stringToPhysicalUnit( 'hbar' ) )
xmls1 = s1.toXML( )
print( xmls1 )
s2 = s1.parseXMLString(xmls1)
if( xmls1 != s2.toXML( ) ) : raise Exception( 'Fix me.' )

print()
suite = spinModule.Suite( )
suite.add( s1 )
print( suite.toXML( ) )

suite2 = suite.parseXMLString(suite.toXML())
if( suite2.toXML( ) != suite.toXML( ) ) : raise Exception( 'Fix me' )

print( '\n============== PARITY ==============' )
p1 = parityModule.Integer( 'nucleus', -1, quantityModule.stringToPhysicalUnit( '' ) )
xmlp1 = p1.toXML( )
print( xmlp1 )
p2 = p1.parseXMLString(xmlp1)

if( xmlp1 != p2.toXML( ) ) : raise Exception( 'Fix me.' )

print()
suite = parityModule.Suite( )
suite.add( p1 )
print( suite.toXML( ) )

suite2 = suite.parseXMLString(suite.toXML())
if( suite2.toXML( ) != suite.toXML( ) ) : raise Exception( 'Fix me' )
