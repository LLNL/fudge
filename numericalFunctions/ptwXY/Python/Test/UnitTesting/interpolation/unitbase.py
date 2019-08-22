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

import sys
sys.path.insert( 0, '../../Utilities' )
sys.path.insert( 0, '../../../../../lib' )

import os
import pointwiseXY_C
import utilities
options = utilities.getOptions( __file__ )

CPATH = '../../../../Test/UnitTesting/interpolation'
accuracy = 1e-3
count = 0

def checkUnitbase2( ls, XY1, XY2 ) :

    global count
    ls, w = utilities.getDoubleValue( 'w', ls )
    ls, wMin = utilities.getDoubleValue( 'wMin', ls )
    ls, wMax = utilities.getDoubleValue( 'wMax', ls )
    ls, XY = utilities.getXYData( ls )
    XYC = pointwiseXY_C.unitbaseInterpolate( w, wMin, XY1, wMax, XY2 )
    utilities.compareXYs( count, 'unitbase', XY, XYC )
    count += 1
    
    return( ls, XYC )

def checkUnitbase( ls ) :

    ls, XY1 = utilities.getXYData( ls )
    ls, XY2 = utilities.getXYData( ls )

    ls, XYl = checkUnitbase2( ls, XY1, XY2 )
    ls, XYr = checkUnitbase2( ls, XY1, XY2 )
    ls, XYm1 = checkUnitbase2( ls, XYl, XYr )
    ls, XYm2 = checkUnitbase2( ls, XY1, XY2 )
    return( ls )

os.system( 'cd %s; make -s clean; ./unitbase -v > v' % CPATH )
f = open( os.path.join( CPATH, 'v' ) )
ls = f.readlines( )
f.close( )

while( len( ls ) ) :
    ls = checkUnitbase( ls )
