# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from fudge import *

import sys
f = open( sys.argv[1] );
ls = f.readlines( )
f.close( )

def getData( ls, datas, label ) :

    for i, l in enumerate( ls ) :
        if( l[:10] == "# label = " ) :
            label = l[10:].strip( )
            if( label not in datas ) : datas[label] = []
        elif( l[:11] == "# length = " ) :
            if( label is None ) : raise Exception( 'label is None, i = %d' % i )
            length = int( l[11:] )
            ls = ls[i+1:]
            data = []
            for i in xrange( length ) :
                data.append( map( float, ls[i].split( ) ) )
            datas[label].append( endl2dmath( data ) )
            ls = ls[length:]
            return( ls, label )
    return( None, label )

datas, label = {}, None
while( ls ) : ls, label = getData( ls, datas, label )

for key in datas :
    data = datas[key]
    if( len( data ) != 2 ) : raise Exception( 'len( datas[%s] ) = %d != 2' % ( key, len( data ) ) )
    orig = data[0]
    orig.label = 'original'
    thinned = data[1]
    thinned.label = 'thinned'
    d = orig - thinned
    d.label = 'diff'
    r = d.safeDivide( orig )
    r.label = 'rel. diff.'
    multiPlot( [ orig, thinned, d, r ], title = key )
