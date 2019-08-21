#! /bin/env python

# <<BEGIN-copyright>>
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
