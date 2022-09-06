#
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
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

x = 1e-4
while( x < 5 ) :
    erf  = specialFunctions.erf( x )
    erfc = specialFunctions.erf( x, True )
    print( x, erf, erfc, erf+ erfc )
    x *= 1.2
