#! /usr/bin/env python3
 
# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys, os, time, multiprocessing
from argparse import ArgumentParser, FileType

description = """Runs multiprocess threading on a general script. For example:

python multiprocessGeneral.py ./convertEndfToGnd.py endl2009.2/endf/*/*.endf
"""

parser = ArgumentParser( description = description )
parser.add_argument( "code",                                                        help = "code to use" )
parser.add_argument( "files", nargs='+',                                            help = "list of file arguments to pass to code" )
parser.add_argument( "-o", "--options", default='',                                 help = "other options to pass to code (use quoted string for multiple args)" )
parser.add_argument( "-n", "--numberOfProcesses",type = int, default = 8,           help = "number of simultaneous processes" )
parser.add_argument( "-l", "--logFile", default = sys.stderr, type = FileType('w'), help = "summary log file" )
parser.add_argument( "-c", "--changeDir", action = 'store_true',                    help = "cd to same directory as file before executing the code" )
parser.add_argument( "-p", "--pythonInterpreter", default=None,                     help = "Python interpreter to use with args.code" )

args = parser.parse_args( )

def f( filename ) :

    CD = ''
    if args.changeDir :
        CD = 'cd %s;' % os.path.dirname( os.path.abspath(filename) )
        filename = os.path.basename(filename)
    if args.pythonInterpreter is not None:
        status = os.system( '%s %s %s %s %s' % ( CD, args.pythonInterpreter, args.code, filename, args.options) )
    else:
        status = os.system( '%s %s %s %s' % ( CD, args.code, filename, args.options) )
    if( status != 0 ) : status = 1
    sys.exit( status )

index = 0
def addProcess( processes, file ) :

    process = multiprocessing.Process( target = f, args = ( file, ) )
    process.start( )
    processes.append( [ index, file, process ] )

n1 = len( args.files )
processes = []
for i1 in range( 0, min( n1, args.numberOfProcesses ) ) :
    addProcess( processes, args.files[index] )
    index += 1

ik = 0
statusMessage = {}
while( ik < n1 ) :
    nextProcesses = []
    for i2, file, process in processes :
        if( process.is_alive( ) ) :
            nextProcesses.append( [ i2, file, process ] )
        else :
            s = "%3d of %3d : %s" % ( i2 + 1, n1, file )
            if( process.exitcode != 0 ) : s += ' ********* FAILED'
            statusMessage[i2] = s
            if( index < n1 ) :
                addProcess( nextProcesses, args.files[index] )
                index += 1
    while( ik in statusMessage ) :       # Print results ordered.
        args.logFile.write( statusMessage[ik] + '\n' )
        del statusMessage[ik]
        ik += 1
    time.sleep( 0.1 )
    processes = nextProcesses
