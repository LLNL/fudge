# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

'''
This module contains the :class:`MultiExecuteFunction`. This class supports running many Python functions in parallel using the Python
multiprocessing module.
'''

# Comments:
#   -) Should files in outputDir be removed in run before processes are started.

import sys
import pathlib
import multiprocessing
import time
import traceback

def runProcess(function, outputFile, **kwargs):
    '''
    This function sets up **stdout** and calls the user function as:

        function(stdout=stdout, **kwargs)

    This function also tries to handle a **raise** executed by *function* in a good way.
    If a **raise** is executed by *function* then, the exceptions traceback information is
    written to *outputFile*.

    When calling **MultiExecuteFunction.run**, if *outputDir* is None, then *outputFile* is **sys.stdout**.

    :param function:        User function to run with keywork arguments *kwargs*.
    :param outputFile:      File name where *function* can write to using something like print('Blah', file=stdout).
    :param kwargs:          Keyword arguments to pass to *function*.
    '''

    stdout = sys.stdout
    if outputFile is not None:
        try:
            stdout = open(outputFile, 'w')
        except:
            traceback.print_exc()
            sys.exit(1)

    try:
        function(stdout=stdout, **kwargs)
        stdout.close()
    except:
        traceback.print_exc(file=stdout)
        stdout.close()
        sys.exit(1)

class MultiExecuteFunction:
    '''
    This class supports running many Python functions in parallel using the Python
    multiprocessing module.

    The following shows an example usage:

        processes = MultiExecuteFunction()
        for index in range(10):
            processes.add(myFunction, 'index = %s' % index, index=index)

        processes.run(numberOfConcurrentProcesses=4, outputDir='work')
    '''

    def __init__(self):

        self.functionsToRun = []

    def add(self, function, message, **kwargs):
        '''
        Adds a function to the list of functions to run in parallel.

        :param function:        Function to call.
        :param message:         Message to print when *function* is complete.
        :param kwargs:          Keyword arguments to pass to *function*.
        '''

        self.functionsToRun.append([function, message, kwargs])

    def run(self, outputDir=None, numberOfConcurrentProcesses=multiprocessing.cpu_count()):
        '''
        The methods runs all functions added to *self* in parallel using *numberOfConcurrentProcesses* concurrent processes.
        '''

        if outputDir is not None:
            outputDir = pathlib.Path(outputDir)
            if not outputDir.exists():
                outputDir.mkdir(parents=True)
            if not outputDir.is_dir():
                raise Exception('outputDir "%s" is not a directory.' % outputDir)

        numberOfProcessesToRun = len(self.functionsToRun)
        processesToRun = []
        for index, functionToRun in enumerate(self.functionsToRun):
            outputFile = outputDir
            if outputFile is not None:
                outputFile = outputDir / ('%.6d.out' % index)
            processesToRun.append(multiprocessing.Process(target=runProcess, args=(functionToRun[0], outputFile), kwargs=functionToRun[2]))
        processesActive = {}
        processesToPrint = {}
        processesDoneCounter = 0

        printIndex = 0
        while(processesDoneCounter < numberOfProcessesToRun):
            while(len(processesActive) < numberOfConcurrentProcesses and len(processesToRun) > 0):
                index = numberOfProcessesToRun - len(processesToRun)
                processesActive[index] = processesToRun.pop(0)
                processesActive[index].start()

            activeProcessIndices = sorted(processesActive)
            for activeProcessIndex in activeProcessIndices:
                process = processesActive[activeProcessIndex]
                if not process.is_alive():
                    processesActive.pop(activeProcessIndex)
                    processesToPrint[activeProcessIndex] = process

            processesToPrintIndices = list(sorted(processesToPrint))
            for processesToPrintIndex in processesToPrintIndices:
                process = processesToPrint[processesToPrintIndex]
                if processesToPrintIndex != printIndex:
                    break
                statusFlag = '********* FAILED' if process.exitcode else ''
                print('%5d of %5d completed: %s:%s' % (processesToPrintIndex + 1, numberOfProcessesToRun, self.functionsToRun[processesToPrintIndex][1], 
                        statusFlag), file=sys.stderr)
                processesToPrint.pop(processesToPrintIndex)
                processesDoneCounter += 1
                printIndex += 1

            time.sleep(0.01)

if __name__=='__main__':

    def myFunction(**kwargs):

        stdout = kwargs['stdout']
        if kwargs['index'] == 3:
            time.sleep(5)
        if kwargs['index'] == 5:
            time.sleep(2)
        if kwargs['index'] == 8:
            time.sleep(8)
        print('Hi', kwargs, file=stdout)
        if kwargs['index'] == 6:
            raise Exception('Raise test.')

    processes = MultiExecuteFunction()
    for index in range(10):
        processes.add(myFunction, 'index = %s' % index, index=index)

    processes.run(numberOfConcurrentProcesses=4, outputDir='work')
