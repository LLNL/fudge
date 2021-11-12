# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
f = open( sys.argv[1] )
ls = f.readlines( )
f.close( )

x1 = None
for l in ls :
    x2, y2 = map( float, l.split( ) )
    y2 *= 3
    if( x1 is not None ) : print("\drawline(%s, %s)(%s, %s)" % (x1, y1, x2, y2))
    x1, y1 = x2, y2
