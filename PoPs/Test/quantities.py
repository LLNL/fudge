#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
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

print( '============== MASS ==============' )
m1 = massModule.double( 'atomic', 2.12, quantityModule.stringToPhysicalUnit( 'amu' ) )
xmlm1 = m1.toXML( )
print( xmlm1 )
m2 = massModule.double.parseXMLStringAsClass( xmlm1 )
if( xmlm1 != m2.toXML( ) ) : raise Exception( 'Fix me.' )

print()
suite = massModule.suite( )
suite.add( m1 )

m2 = massModule.double( 'nuclear', 2, quantityModule.stringToPhysicalUnit( 'amu' ) )
suite.add( m2 )

print( suite.toXML( ) )

suite2 = massModule.suite.parseXMLStringAsClass( suite.toXML( ) )
if( suite2.toXML( ) != suite.toXML( ) ) : raise Exception( 'Fix me' )

print( '\n============== CHARGE ==============' )
c1 = chargeModule.integer( 'nucleus', 3, quantityModule.stringToPhysicalUnit( 'e' ) )
xmlc1 = c1.toXML( )
print( xmlc1 )
c2 = chargeModule.integer.parseXMLStringAsClass( xmlc1 )
if( xmlc1 != c2.toXML( ) ) : raise Exception( 'Fix me.' )

print()
suite = chargeModule.suite( )
suite.add( c1 )
print( suite.toXML( ) )

suite2 = suite.parseXMLStringAsClass( suite.toXML( ) )
if( suite2.toXML( ) != suite.toXML( ) ) : raise Exception( 'Fix me' )

print( '\n============== HALFLIFE ==============' )
h1 = halflifeModule.double( 'nucleus', 3.14e6, quantityModule.stringToPhysicalUnit( 'd' ) )
xmlh1 = h1.toXML( )
print( xmlh1 )
h2 = halflifeModule.double.parseXMLStringAsClass( xmlh1 )
if( xmlh1 != h2.toXML( ) ) : raise Exception( 'Fix me.' )

print()
suite = halflifeModule.suite( )
suite.add( h1 )
print( suite.toXML( ) )

suite2 = suite.parseXMLStringAsClass( suite.toXML( ) )
if( suite2.toXML( ) != suite.toXML( ) ) : raise Exception( 'Fix me' )

print()
print( h1.pqu( ) )
print( suite[0].pqu( ) )
print( suite[0].pqu( 's' ) )
print( PQUModule.floatToShortestString( suite[0].pqu( ).getValueAs( 's' ), 12 ) )

print()
suite = halflifeModule.suite( )
h2 = halflifeModule.string( 'nucleus', 'stable', quantityModule.stringToPhysicalUnit( 's' ) )
suite.add( h2 )
print( suite.toXML( ) )

suite2 = suite.parseXMLStringAsClass( suite.toXML( ) )
if( suite2.toXML( ) != suite.toXML( ) ) : raise Exception( 'Fix me' )

print( '\n============== SPIN ==============' )
fraction = quantityModule.fraction.toValueType( "5/2" )
s1 = spinModule.fraction( 'nucleus', fraction, quantityModule.stringToPhysicalUnit( 'hbar' ) )
xmls1 = s1.toXML( )
print( xmls1 )
s2 = s1.parseXMLStringAsClass( xmls1 )
if( xmls1 != s2.toXML( ) ) : raise Exception( 'Fix me.' )

print()
suite = spinModule.suite( )
suite.add( s1 )
print( suite.toXML( ) )

suite2 = suite.parseXMLStringAsClass( suite.toXML( ) )
if( suite2.toXML( ) != suite.toXML( ) ) : raise Exception( 'Fix me' )

print( '\n============== PARITY ==============' )
p1 = parityModule.integer( 'nucleus', -1, quantityModule.stringToPhysicalUnit( '' ) )
xmlp1 = p1.toXML( )
print( xmlp1 )
p2 = p1.parseXMLStringAsClass( xmlp1 )

if( xmlp1 != p2.toXML( ) ) : raise Exception( 'Fix me.' )

print()
suite = parityModule.suite( )
suite.add( p1 )
print( suite.toXML( ) )

suite2 = suite.parseXMLStringAsClass( suite.toXML( ) )
if( suite2.toXML( ) != suite.toXML( ) ) : raise Exception( 'Fix me' )
