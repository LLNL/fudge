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
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
