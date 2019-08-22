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

import os, subprocess, glob, sys

def _getSTDStreamOpen( stdStream ) :

    if( stdStream is None ) : return( None )
    if( stdStream == subprocess.PIPE ) : return( stdStream )
    if( type( stdStream ) == type( 1 ) ) : return( stdStream )
    if( type( stdStream ) == type( '' ) ) : return( open( stdStream, 'w' ) )
    if( type( stdStream ) == type( sys.stdout ) ) : return( stdStream )
    raise Exception( 'Unsupported stream = "%s"' % type( stdStream ) )

def _getSTDStreamClose( stdStream, stdStream2, processStd ) :

    if( stdStream == subprocess.PIPE ) : return( processStd.readlines( ) )
    if( type( stdStream ) == type( '' ) ) :
        stdStream2.close( )
        return( stdStream )
    return( None )

def executeCommand( args, raiseOnError = True, useExecutable = False, stdout = subprocess.PIPE, stderr = subprocess.PIPE ) :

    stdout2, stderr2, stdin, shell = _getSTDStreamOpen( stdout ), _getSTDStreamOpen( stderr ), subprocess.PIPE, False
    os.environ.update( { 'PYTHONPATH' : ':'.join(os.sys.path ) } )
    try :
        if( useExecutable ) :
            executable = args[0]
            if( os.path.exists( executable ) ) : executable = os.path.realpath( args[0] )
            process = subprocess.Popen( args, shell = shell, stdin = stdin, stdout = stdout2, stderr = stderr2, executable = executable )
        else :
            process = subprocess.Popen( args, shell = shell, stdin = stdin, stdout = stdout2, stderr = stderr2 )
    except :
        print args
        print 'Execution of "%s" FAILED' % args[0]
        raise
    process.wait( )
    stdout_results = _getSTDStreamClose( stdout, stdout2, process.stdout )
    stderr_results = _getSTDStreamClose( stderr, stderr2, process.stderr )
    if( raiseOnError ) :
        if( process.returncode != 0 ) :
            if( type( stderr_results ) == type( [] ) ) :
                sys.stderr.write( ''.join( stderr_results ) + '\n' )
            elif( type( stderr_results ) == type( '' ) ) :
                sys.stderr.write( 'Error directed to file "%s"' % stderr_results )
            raise Exception( 'Execution of "%s" FAILED with status = %s' % ( args[0], process.returncode ) )
    return( process.returncode, stdout_results, stderr_results )

def spawn( args ) :

    os.environ.update( { 'PYTHONPATH' : ':'.join(os.sys.path ) } )
    sp = subprocess.Popen( args )
    return( sp.pid )

def deleteFilesUsingGlob( patterns ) :

    if( type( patterns ) == type( '' ) ) : patterns = [ patterns ]
    for pattern in patterns :
        files = glob.glob( pattern )
        for file in files : 
            if( os.path.isdir( file ) ) :
                pass
            else :
                os.remove( file )
