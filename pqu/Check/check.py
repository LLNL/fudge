# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
import glob, os, shutil, filecmp

PYTHON = sys.executable

files = sorted( glob.glob( 't*.py' ) )

if( os.path.exists( 'Out' ) ) : shutil.rmtree( 'Out' )
os.mkdir( 'Out' )

for file in files :
    base = file[:-3]
    status = os.system( '%s %s > Out/%s.out' % ( PYTHON, file, base ) )
    if( status ) : print( '=========== %s ===========' % file )

outs = sorted( glob.glob( 'Out/t*.out' ) )
for out in outs :
    file = os.path.basename( out )
    if( not( filecmp.cmp( os.path.join( 'Out.checked', file ), out ) ) ) : print( 'ERROR: %s' % out )
