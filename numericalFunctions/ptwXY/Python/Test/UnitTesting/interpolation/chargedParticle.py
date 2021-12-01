#! /bin/env python

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import sys
sys.path.insert( 0, '../../Utilities' )

from numericalFunctions import pointwiseXY_C

import utilities

options = utilities.getOptions( __file__ )

CPATH = '../../../../Test/UnitTesting/interpolation'
accuracy = None, None

def getXYData( fileName, getParameters = False ) :

    global accuracy

    f = open( os.path.join( CPATH, fileName ) )
    ls = f.readlines( )
    f.close( )
    if( getParameters ) :
        if( 'accuracy' not in ls[0] ) : raise Exception( 'Sparse file missing "accuracy" data.' )
        accuracy = float( ls.pop( 0 ).split( '=' )[1] )
    xys = utilities.getXYData( ls )[1]
    return( xys )

os.system( 'cd %s; make -s clean; ./chargedParticle' % CPATH )
sparse = getXYData( 'curve_sparse.dat', True )
sparse = pointwiseXY_C.pointwiseXY_C( [ xy for xy in sparse ], interpolation = 'charged-particle' )
