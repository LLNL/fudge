#! /usr/bin/env python3
 
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
import os
import time
import json
import multiprocessing
import argparse

description = """Runs a code on a list of files. For each file of the list, the code is run 
on a separate thread. For example, the following will execute the code "./convertEndfToGnd.py"
for each file matching "endl2009.2/endf/*/*.endf" on a unique thread:

python multiprocessGeneral.py ./convertEndfToGnd.py endl2009.2/endf/*/*.endf
"""

parser = argparse.ArgumentParser( description = description, formatter_class = argparse.RawTextHelpFormatter )
parser.add_argument( 'code',                                                        help = 'Code to use.' )
parser.add_argument( 'files', nargs='+',                                            help = 'List of files to be processed by the code. One code/file per thread.' )
parser.add_argument( '-o', '--options', default='',                                 help = 'Options to pass to code (use quoted string for multiple args).' )
parser.add_argument( "-j", '--optionsJSON', default=None,                           help = "location of JSON file with individual file command line arguments")
parser.add_argument( '-n', '--numberOfProcesses',type = int, default = 8,           help = 'Number of simultaneous processes (i.e., threads).' )
parser.add_argument( '-l', '--logFile', default = sys.stderr, type = argparse.FileType('w'), help = 'Log file.' )
parser.add_argument( '-p', '--pythonInterpreter', default = '',                     help = 'Python interpreter to use with code if it is a Python script.' )
parser.add_argument( '-v', '--verbose', action = 'count', default = 0,              help = 'Prints out the cmd executed for each file. This is meant for debugging.' )

group = parser.add_mutually_exclusive_group( )
group.add_argument( '-c', '--changeDir', action = 'store_true',                     help = 'Cd to same directory as file before executing the code. This options and "--cd" are mutually exclusive.')
group.add_argument( '--cd', action = 'store', default = '',                         help = 'Cd to the directory specified by the next argument. If the directory does not exist, it is created.' )

args = parser.parse_args( )

if args.optionsJSON is not None:
    with open(args.optionsJSON) as fileObject:
        individualFileOptions = json.load(fileObject)


def f( filename ) :

    fileOptions = args.options if args.optionsJSON is None else '%s %s' % (individualFileOptions[filename], args.options)

    CD = ''
    if args.cd != '':
        CD = 'cd %s;' % args.cd
        filename = os.path.realpath(filename)
    elif args.changeDir :
        CD = 'cd %s;' % os.path.dirname( os.path.abspath(filename) )
        filename = os.path.basename(filename)

    code = os.path.realpath(args.code)
    if args.pythonInterpreter is not None:
        cmd = '%s %s %s %s %s' % ( CD, args.pythonInterpreter, code, filename, fileOptions)
    else:
        cmd = '%s %s %s %s' % ( CD, args.code, filename, fileOptions)

    if args.verbose: print(cmd)
    status = os.system(cmd)
    if( status != 0 ) : status = 1
    sys.exit( status )

if args.cd != '':
    if not os.path.exists(args.cd): os.makedirs(args.cd)

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

indexFormat = '%%%dd' % len('%s' % n1)
format = '%s of %s: %%s' % ( indexFormat, indexFormat)
ik = 0
statusMessage = {}
while( ik < n1 ) :
    nextProcesses = []
    for i2, file, process in processes :
        if( process.is_alive( ) ) :
            nextProcesses.append( [ i2, file, process ] )
        else :
            s = format % ( i2 + 1, n1, file )
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
