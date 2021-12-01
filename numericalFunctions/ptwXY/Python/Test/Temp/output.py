# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import glob, os, sys

stdouts = glob.glob( '*.stdout' )
stdouts.sort( )

errors = 0

def getThinFile( file ) :

    f = open( file )
    ls = f.readlines( )
    f.close( )
    ns = []
    for l in ls :
        if( l[0] != '#' ) : ns.append( l )
    return( ns )

for stdout in stdouts :
    stderr = stdout.replace( 'stdout', 'stderr' )
    errStat = os.stat( stderr )
    if( errStat.st_size == 0 ) : 
        os.system( 'rm -f %s' % stderr )
    else :
        errors += 1
        sys.stderr.write( 'FAILURE: %s is not empty\n' % stderr )
    result = os.path.join( '..', 'Results', stdout )
    ls1 = getThinFile( stdout )
    ls2 = getThinFile( result )
    if( ls1 == ls2 ) : 
        os.system( 'rm -f %s' % stdout )
    else :
        errors += 1
        sys.stderr.write( 'FAILURE: %s differs from reference\n' % stdout )

if( os.system( 'diff t6.out ../Results/t6.out 2> /dev/null' ) == 0 ) :
    os.system( 'rm -f t6.out' )
else :
    errors += 1
    sys.stderr.write( 'FAILURE: t6.out differs from reference\n' )

if( os.system( 'diff t7.orig.out t7.copy.out 2> /dev/null' ) == 0 ) :
    os.system( 'rm -f t7.copy.out' )
else :
    errors += 1
    sys.stderr.write( 'FAILURE: t7.orig.out and t7.copy.out differs\n' )

if( os.system( 'diff t7.orig.out ../Results/t7.orig.out 2> /dev/null' ) == 0 ) :
    os.system( 'rm -f t7.orig.out' )
else :
    errors += 1
    sys.stderr.write( 'FAILURE: t7.orig.out differs from reference\n' )

files = glob.glob( '*' )
files.remove( 'Makefile' )
files.remove( __file__ )
if( len( files ) > 0 ) :
    errors += 1
    sys.stderr.write( 'FAILURE: extra files not checked.\n' )
    sys.stderr.write( '    %s\n' % repr( files ) )

sys.exit( errors )
