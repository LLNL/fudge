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
    sys.stderr.write( '    %s\n' % `files` )

sys.exit( errors )
