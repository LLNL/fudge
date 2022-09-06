# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os, sys, glob

CHECKOPTIONS = ''
if( 'CHECKOPTIONS' in os.environ  ) : CHECKOPTIONS = ' ' + os.environ['CHECKOPTIONS']
file = sys.argv[1]
os.system( './' + file + CHECKOPTIONS )
dats = glob.glob( 'Results/%s/*' % file )
for dat in dats :
    base = os.path.basename( dat )
    status = os.system( 'diff %s %s > /dev/null 2>&1' % ( dat, base ) )
    if( status ) :
        print( 'FAILURE for %s' % file )
