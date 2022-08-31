# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import sys
import subprocess
import glob

def _getSTDStreamOpen(stdStream):
    """
    If stdStream is a string, open file of that name and return file handle.
    Otherwise check that it's an allowed type (int, buffer or None) and return as-is.
    """

    if stdStream is None:
        return None
    if stdStream == subprocess.PIPE:
        return stdStream
    if isinstance(stdStream, int):
        return stdStream
    if isinstance(stdStream, str):
        return open(stdStream, 'w')
    if isinstance(stdStream, type(sys.stdout)):
        return stdStream
    raise Exception('Unsupported stream = "%s"' % type(stdStream))

def _getSTDStreamClose(stdStream, stdStream2, processStd):
    """
    Close stream. If output is a file, write results and close the file. If output is a pipe, get result and decode to list of strings.
    """

    if stdStream == subprocess.PIPE:
        results = processStd.readlines()
        return map(bytes.decode, results)
    if isinstance(stdStream, str):
        stdStream2.close()
        return stdStream
    return None

def executeCommand(args, raiseOnError=True, useExecutable=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
    """
    Executes a command using subprocess. Supports redirecting stdout / stderr to files (or to each other).
    If stdout / stderr are strings, they will be treated as file names and output written to them.
    If they are pipes (default option), they will be returned as lists of strings along with the exit code.

    @param args: list of arguments to be executed, e.g. ['echo', 'Hello World!']
    @param raiseOnError: if True, raise Exception if the subprocess returns a non-zero exit code.
    @param useExecutable: if True, replace the executable with absolute path before running.
    @param stdout: where to send standard out. Options include file name, PIPE, or None (like sending to /dev/null).
    @param stderr: where to send standard error. Same options as stdout.
    @return: (return code, stdout, stderr)
    """

    stdout2, stderr2, stdin, shell = _getSTDStreamOpen(stdout), _getSTDStreamOpen(stderr), subprocess.PIPE, False
    os.environ.update({'PYTHONPATH': ':'.join(sys.path)})
    try:
        if useExecutable:
            executable = args[0]
            if os.path.exists(executable):
                executable = os.path.realpath(args[0])
            process = subprocess.Popen(args, shell=shell, stdin=stdin, stdout=stdout2, stderr=stderr2,
                                       executable=executable)
        else:
            process = subprocess.Popen(args, shell=shell, stdin=stdin, stdout=stdout2, stderr=stderr2)
    except Exception:
        print(args)
        print('Execution of "%s" FAILED' % args[0])
        raise

    process.wait()
    stdout_results = _getSTDStreamClose(stdout, stdout2, process.stdout)
    stderr_results = _getSTDStreamClose(stderr, stderr2, process.stderr)
    if raiseOnError:
        if process.returncode != 0:
            if isinstance(stderr_results, list):
                sys.stderr.write(''.join(stderr_results) + '\n')
            elif isinstance(stderr_results, str):
                sys.stderr.write('Error directed to file "%s"' % stderr_results)
            raise Exception('Execution of "%s" FAILED with status = %s' % (args[0], process.returncode))
    return process.returncode, stdout_results, stderr_results

def spawn(args):
    """
    Shortcut for launching a process with PYTHONPATH set. Returns process id.
    """

    os.environ.update({'PYTHONPATH': ':'.join(sys.path)})
    sp = subprocess.Popen(args)
    return sp.pid

def deleteFilesUsingGlob(patterns):
    """
    Deletes file(s) matching patterns. Skips deleting directories.
    """

    if isinstance(patterns, str):
        patterns = [patterns]
    for pattern in patterns:
        files = glob.glob(pattern)
        for file in files:
            if os.path.isdir(file):
                pass
            else:
                os.remove(file)
