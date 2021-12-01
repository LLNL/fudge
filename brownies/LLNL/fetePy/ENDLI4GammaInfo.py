#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from argparse import ArgumentParser

import brownies.legacy.converting.endf_endl as endf_endlModule
import brownies.legacy.endl.endlmisc as endlmiscModule

energyEps = 1e-6

usage = \
"""
Replaces all gamma data (i.e., yo07* and *c55* data) in an ENDL neutron database with the gamma data produced by fetePy.py.
"""

parser = ArgumentParser( description = usage )
parser.add_argument( 'file', type = str,
        help = 'Path to ENDL yo = 7, I = 4 file.' )
parser.add_argument( '--energyEps', default = energyEps, action = 'store', type = float,
        help = 'Estimate of fractional delta function width (default = %s).' % energyEps )
parser.add_argument( '-v', '--verbose', default = 0, action = 'count',
        help = 'Enable verbose output.' )

args = parser.parse_args( )

energyEps = 2 * args.energyEps

CS2MT = endf_endlModule.ENDLCS_To_ENDFMT( 1 )

def process( reaction, C, S ) :

    lines = reaction.split( '\n' )
    lines.pop( 0 )
    line = lines.pop( 0 )
    X1 = line[20:32]
    X1 = endlmiscModule.headerString2FunkyDouble( X1, 0, __file__ )
    X1 = float( X1 )
    MT = CS2MT.getMTFromCS( C, S )
    print( 'X1 = %s: for MT = %s' % ( X1, MT ) )
    energyInData = {}
    for line in lines :
        eIn, Ep, l, P = map( float, line.split( ) )
        if( l > 0 ) : break
        if( eIn not in energyInData ) : energyInData[eIn] = []
        energyInData[eIn].append( [ Ep, P ] )

    for energy in sorted( energyInData.keys( ) ) :
        print( '    incident energy = %s' % energy )
        data = energyInData[energy]
        n1 = len( data )
        sum = 0
        for i1 in range( 1, n1 - 1 ) :
            if( ( ( data[i1][0] - data[i1-1][0] ) < energyEps * data[i1][0] ) and
                ( ( data[i1+1][0] - data[i1][0] ) < energyEps * data[i1][0] ) ) :
                if( ( data[i1-1][1] < 1e-3 * data[i1][1] ) and
                    ( data[i1+1][1] < 1e-3 * data[i1][1] ) ) :
                        probability = 0.5 * ( ( data[i1][0]   - data[i1-1][0] ) * ( data[i1][1] - data[i1-1][1] ) +
                                              ( data[i1+1][0] - data[i1][0]   ) * ( data[i1][1] - data[i1+1][1] ) )
                        sum += probability
                        print( '        discrete: %-12s %.8f' % ( data[i1][0], probability ) )
        print( '        continuum: =========== %.8e' % ( 1 - sum ) )

fIn = open( args.file )
lines = ''.join( fIn.readlines( ) )
fIn.close( )

reactions = lines.split( '\n                                                                       1\n' )[:-1]
header2 = reactions[0].split( '\n' )[1]
C = int( header2[:2] )
I = int( header2[2:5] )
if( I != 4 ) : raise Exception( 'I=%d not supported. Only run with I=4 files' % I )
S = int( header2[5:8] )
for reaction in reactions : process( reaction, C, S )
