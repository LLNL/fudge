# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
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

def executeCommand( args, raiseOnError = True, useExecutable = False, stdout = subprocess.PIPE, stderr = subprocess.PIPE, bytesToStr = True ) :

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
        print( args )
        print( 'Execution of "%s" FAILED' % args[0] )
        raise
    process.wait( )
    stdout_results = _getSTDStreamClose( stdout, stdout2, process.stdout )
    stderr_results = _getSTDStreamClose( stderr, stderr2, process.stderr )
    if( bytesToStr and ( sys.version_info.major > 2 ) ) :
        if( stdout_results is not None ) : stdout_results = map( bytes.decode, stdout_results )
        if( stderr_results is not None ) : stderr_results = map( bytes.decode, stderr_results )
    if( raiseOnError ) :
        if( process.returncode != 0 ) :
            if( type( stderr_results ) == type( [] ) ) :
                sys.stderr.write( ''.join( map( bytes.decode, stderr_results ) ) + '\n' )
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
