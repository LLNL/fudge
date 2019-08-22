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

import os
from argparse import ArgumentParser

import fudge.legacy.converting.endf_endl as endf_endlModule
import fudge.legacy.endl.endlmisc as endlmiscModule

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
    print 'X1 = %s: for MT = %s' % ( X1, MT )
    energyInData = {}
    for line in lines :
        eIn, Ep, l, P = map( float, line.split( ) )
        if( l > 0 ) : break
        if( eIn not in energyInData ) : energyInData[eIn] = []
        energyInData[eIn].append( [ Ep, P ] )

    for energy in sorted( energyInData.keys( ) ) :
        print '    incident energy = %s' % energy
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
                        print '        discrete: %-12s %.8f' % ( data[i1][0], probability )
        print '        continuum: =========== %.8e' % ( 1 - sum )

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
