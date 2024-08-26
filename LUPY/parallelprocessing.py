# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""         
This module contains functions for running tasks in parallel mainly through :py:func:`execute`
and :py:func:`executeMultiThreaded`.
                
    This module contains the following functions:

    +-----------------------+-----------------------------------------------------------------------+
    | Function              | Description                                                           |
    +=======================+=======================================================================+
    | executeOne            | This function is mainly for an example.                               | 
    +-----------------------+-----------------------------------------------------------------------+
    | addProcess            | This functions is internal use by :py:func:`execute`.                 |
    +-----------------------+-----------------------------------------------------------------------+
    | execute               | This function is used to run many processes in parallel.              |
    +-----------------------+-----------------------------------------------------------------------+
    | executeMultiThreaded  | This function is used to run many tasks in parallel as threads.       |
    +-----------------------+-----------------------------------------------------------------------+
"""     

import sys
import os
import multiprocessing
import time

def executeOne( command ) :
    """
    This function is mainly for an example of what the target argument of :py:class:`multiprocessing.Process` may look.
    The main thing is to call sys.exit with a status and to properly handle the status returned by os.system.

    :param command:     System comannd to execute.
    """

    status = os.system( command )
    if( status != 0 ) : status = 1
    sys.exit( status )

def addProcess( numberLaunched, processes, items, callback ) :
    """
    This function calls :py:func:`multiprocessing.Proces` with target of *callback* and args of items[numberLaunched] to create 
    a new process, and then calls start on that process. The created process is appended to *processes*.
    This function is for internal use by :py:class:`execute`.

    :param numberLaunched:  Index into *items* which contains the arguments for the callback function.
    :param processes:       The list of processes which the created process will be added to.
    :param items:           The list of args for all processes to create.
    :param callback:        That function called by :py:class:`multiprocessing.Process`.

    :returns:               *numberLaunched* incremented by one.
    """

    process = multiprocessing.Process( target = callback, args = ( items[numberLaunched], ) )
    process.start( )
    processes.append( [ numberLaunched, items[numberLaunched], process ] )
    return( numberLaunched + 1 )

def execute(numberOfParallelProcesses, items, callback, percentNotifications=100):
    """
    This function loops over each item in *items*, and creates and runs each with a :py:func:`multiprocessing.Proces` process.
    Also see function :py:func:`addProcess`.  At most *numberOfParallelProcesses* of
    processes can run in parallel. The user's callback must create a command to execute and a backcall function similar
    to the example :py:func:`executeOne`.

    percentNotifications is in percent of number to do and for every percent of number done a message is sent to sys.syserr.

    :param numberOfParallelProcesses:   The maximum number of process to have running at one time.
    :param items:                       A python list where each element contains the arguments for one process.
    :param callback:                    The target passed to :py:class:`multiprocessing.Process`.
    :param percentNotifications:        This parameter determines how often a progress message is printed.
    """

    numberOfFailures = 0
    numberLaunched = 0
    numberToDo = len( items )
    processes = []
    for i1 in range( 0, min( numberToDo, numberOfParallelProcesses ) ) : numberLaunched = addProcess( numberLaunched, processes, items, callback )

    syserrNotificationIncrement = numberToDo + 1
    if 0.99 < percentNotifications < 50.1:
        syserrNotificationIncrement = int(percentNotifications * numberToDo / 100.)
    syserrNotificationCounter = syserrNotificationIncrement

    numberDone = 0
    done = {}
    while( numberDone < numberToDo ) :
        nextProcesses = []
        for i1, item, process in processes :
            if( process.is_alive( ) ) :
                nextProcesses.append( [ i1, item, process ] )
            else :
                message = "%3d of %3d : %s" % ( i1 + 1, numberToDo, item )
                if( process.exitcode != 0 ) :
                    numberOfFailures += 1
                    message += ' ********* FAILED'
                done[i1] = message
                if( numberLaunched < numberToDo ) : numberLaunched = addProcess( numberLaunched, nextProcesses, items, callback )
        while( numberDone in done ) :       # Print results ordered.
            print(done[numberDone])
            del done[numberDone]
            numberDone += 1
            if numberDone == syserrNotificationCounter:
                syserrNotificationCounter += syserrNotificationIncrement
                sys.stderr.write('    %5d of %5d done\n' % ( numberDone, numberToDo ))
        time.sleep( 0.1 )
        processes = nextProcesses

    return numberOfFailures

def executeMultiThreaded(taskList, runner, nThreads=multiprocessing.cpu_count(), sleepTime=0.2, verbose=False):
    """
    The function runs the items in *taskList* as separate threads. Applies *runner* to each task in the taskList and
    returns list of completed runners once all tasks are complete.  Note: order of returned list may not match order of taskList.

    Due to the Python GIL, this should be most efficient for runners that spend most of their time
    in compiled code e.g., C extensions rather than in interpreted Python code.

    :param taskList:    List of tasks to run.
    :param runner:      Class (inheriting from threading.Thread) with a run() method that completes each task.
    :param nThreads:    Maximum number of concurrent threads.
    :param sleepTime:   How long (in seconds) to wait between polling for finished threads.
    :param verbose:     If True, print updates as each task starts.

    :returns:           List of finished tasks.
    """
#   TODO: add index to each thread in case we need to sort results? Could use thread.name.

    nTasks = len(taskList)
    if verbose:
        print(f"Executing {nTasks} tasks with multithreading")

    index = 0
    active = []
    maxThreads = min(nThreads, nTasks)
    while len(active) < maxThreads:
        if verbose:
            print(f"  Start task {index} of {nTasks}")
        active.append(runner(taskList.pop()))
        active[-1].start()
        index += 1

    finished = []
    while True:
        time.sleep(sleepTime)
        for t in active[::-1]:
            if not t.is_alive():
                t.join()
                finished.append(t)
                active.remove(t)
        while len(active) < maxThreads and len(taskList) != 0:
            if verbose:
                print(f"  Start task {index} of {nTasks}")
            active.append(runner(taskList.pop()))
            active[-1].start()
            index += 1

        if len(active)==0 and len(taskList)==0:
            break

    if verbose:
        print("All tasks completed")
    return finished
