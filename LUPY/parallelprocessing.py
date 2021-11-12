# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
import os
import multiprocessing
import time

def executeOne( command ) :
    """
    This function is mainly for an example of what the target argument of multiprocessing.Process may look.
    The main thing is to call sys.exit with a status and to properly handle the status returned by os.system.
    """

    status = os.system( command )
    if( status != 0 ) : status = 1
    sys.exit( status )

def addProcess( numberLaunched, processes, items, callback ) :
    """Added a process."""

    process = multiprocessing.Process( target = callback, args = ( items[numberLaunched], ) )
    process.start( )
    processes.append( [ numberLaunched, items[numberLaunched], process ] )
    return( numberLaunched + 1 )

def execute( numberOfParallelProcesses, items, callback ) :
    """
    This function loops over each item in *items* and processes it by calling user's callback function. At most *numberOfParallelProcesses* of
    processes can run in parallel. The user's callback must create a command to execute and call function **executeOne**.
    """

    numberLaunched = 0
    numberToDo = len( items )
    processes = []
    for i1 in range( 0, min( numberToDo, numberOfParallelProcesses ) ) : numberLaunched = addProcess( numberLaunched, processes, items, callback )

    numberDone = 0
    done = {}
    while( numberDone < numberToDo ) :
        nextProcesses = []
        for i1, item, process in processes :
            if( process.is_alive( ) ) :
                nextProcesses.append( [ i1, item, process ] )
            else :
                message = "%3d of %3d : %s" % ( i1 + 1, numberToDo, item )
                if( process.exitcode != 0 ) : message += ' ********* FAILED'
                done[i1] = message
                if( numberLaunched < numberToDo ) : numberLaunched = addProcess( numberLaunched, nextProcesses, items, callback )
        while( numberDone in done ) :       # Print results ordered.
            print(done[numberDone])
            del done[numberDone]
            numberDone += 1
        time.sleep( 0.1 )
        processes = nextProcesses
