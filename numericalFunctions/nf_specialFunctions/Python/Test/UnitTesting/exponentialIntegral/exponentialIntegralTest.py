#
# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
#

import os
import sys
import glob

import numericalFunctions

# Deployment via Makefiles places shared libraries in .../fudge/numericalFunctions/lib
if len(glob.glob(os.path.join(numericalFunctions.__path__[0], 'lib', '*specialFunctions*'))) > 0:
    from numericalFunctions.lib import specialFunctions

# Deployment via `pip install` places shared libraries in .../site-packages/numericalFunctions
else:
    from numericalFunctions import specialFunctions

options = []
if( 'CHECKOPTIONS' in os.environ ) : options = os.environ['CHECKOPTIONS'].split( )
for argv in sys.argv[1:] : options += [ argv ]

if( '-e' in options ) : print( __file__ )

f = open( '../../../../Test/UnitTesting/exponentialIntegral/exponentialIntegralTest.dat' )
ls = f.readlines( )
f.close( )

errors = 0
for l in ls :
    if( ';' in l ) : break
    if( ',' not in l ) : continue
    n, x, f = l.split( '{' )[1].split( '}' )[0].split( ',' )
    n = int( n )
    x = float( x )
    f = float( f )
    En = specialFunctions.exponentialIntegral( n, x )
    r = 1
    if( f != 0 ) : r = En / f - 1
    if( abs( r ) > 3e-14 ) :
        print( n, x, f, En, r )
        errors += 1

if( errors > 0 ) : raise Exception( "FAILED: %s" % __file__ )
