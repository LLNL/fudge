# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# <<END-copyright>>

from fudge import *

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
