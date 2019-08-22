#! /bin/env python

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

"""
Some simple stuff that is not repeatable as pointers are printed.
"""
import sys
sys.path.append( '../../../lib' )

import pointwiseXY_C

def setShow( f, x, y ) :
    print
    print 'Setting x, y to', x, y
    f.setValue( x, y )
    f.showInteralStructure( printPointersAsNull = True )

f = pointwiseXY_C.pointwiseXY_C( initialSize = 6, overflowSize = 3 )
f.showInteralStructure( printPointersAsNull = True )

setShow( f, 1., -1 )
setShow( f, 0.1, -2 )
setShow( f, 0.091, -1.2 )
setShow( f, 0.91, -1.3 )
setShow( f, 91, -1.6 )
setShow( f, 1e-3, -1.05 )
setShow( f, 1e-9, -1.15 )
setShow( f, 1e-6, -1.25 )
