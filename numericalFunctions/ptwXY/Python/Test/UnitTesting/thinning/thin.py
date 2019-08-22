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

import sys
sys.path.insert( 0, '../../../../../lib' )

import os
import pointwiseXY_C

if( 'CHECKOPTIONS' in os.environ ) :
    options = os.environ['CHECKOPTIONS'].split( )
    if( '-e' in options ) : print __file__

CPATH = '../../../../Test/UnitTesting/thinning'

os.system( 'cd %s; thin -v > v' % CPATH )
f = open( os.path.join( CPATH, 'v' ) )
ls = f.readlines( )
f.close( )

def getData( ls, hasLabel ) :

    i = 0
    for l in ls :
        if( l.strip( ) != '' ) : break
        i = i + 1
    ls = ls[i:]
    if( len( ls ) == 0 ) : return( None, None, None, None )
    label, accuracy = None, None
    if( hasLabel ) :
        label, ls = ls[0].strip( ), ls[1:]
        accuracy, ls = ls[0].strip( ), ls[1:]
        if( accuracy[:13] != '# accuracy = ' ) : raise Exception( '%s: line does not contain accuracy info: "%s"' % ( __file__, accuracy[:-1] ) )
        accuracy = float( accuracy.split( '=' )[1] )
    length, ls = ls[0], ls[1:]
    if( '# length = ' != length[:11] ) : raise Exception( '%s: line does not contain length info: "%s"' % ( __file__, length.strip( ) ) )
    length = int( length.split( '=' )[1] )
    data = [ map( float, ls[i].split( )[:2] ) for i in xrange( length ) ]
    return( ls[length:], label, accuracy, pointwiseXY_C.pointwiseXY_C( data, initialSize = 100, overflowSize = 10 ) )

def compareValues( label, i, v1, v2 ) :

    sv1, sv2 = '%.7g' % v1, '%.7g' % v2
    if( sv1 != sv2 ) : raise Exception( '%s: values %s %s diff at %d for label = %s' % ( __file__, v1, v2, i, label ) )

def thin( label, accuracy, original, data ) :

    thin = original.thin( accuracy )
    if( len( data ) != len( thin ) ) : raise Exception( '%s: len( data ) = %d != len( thin ) = %d for label = "%s"' % \
        ( __file__, len( data ), len( thin ), label ) )
    for i, xy in enumerate( data ) :
        xc, yc = xy
        xp, yp = thin[i]
        compareValues( label, i, xc, xp )
        compareValues( label, i, yc, yp )

while( 1 ) :
    ls, label, accuracy, original = getData( ls, True )
    if( ls is None ) : break
    ls, label, dummy, data = getData( ls, False )
    if( ls is None ) : raise Exception( '%s: missing thinned data for label = "%s"' % ( __file__, label ) )
    thin( label, accuracy, original, data )
