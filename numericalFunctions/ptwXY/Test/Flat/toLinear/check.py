# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys, os

file = sys.argv[1]

status = os.system( 'diff v %s.dat > /dev/null' % file )
if( status != 0 ) : print( "ERROR for %s" % file )
