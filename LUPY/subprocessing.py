# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains several functions that are useful for creating sub-processes using :py:class:`subprocess.Popen`.
"""

import os
import sys
import subprocess
import glob

def _getSTDStreamOpen(stdStream):
    """
    If *stdStream* is a string, a file is create for writing with that path *stdStream* and its file handle
    is returned.  Otherwise check that it's an allowed type (int, buffer or None) and return as-is. For internal use only.

    :param stdStream:   Any file handle type supported by :py:class:`subprocess.Popen` or a Python str.
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
    If stream is a file, write results to the file and close the file. 
    If stream is a pipe, get result and decode them to a list of strings.
    For internal use only.

    :param stdStream:       Instance passed to :py:func:`_getSTDStreamOpen`.
    :param stdStream2:      Instance returned by :py:func:`_getSTDStreamOpen` that is associated with *stdStream*.
    :param processStd:      A stdout or stderr instance associated with a :py:class:`subprocess.Popen` instance.
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
    This function executes a command using :py:class:`subprocess.Popen`. It supports redirecting stdout / stderr 
    to files (or to each other).  If stdout / stderr are strings, they will be treated as file paths. These
    paths are opened and output written to them.
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
    Shortcut for launching a process with PYTHONPATH set from the current sys.path. Calls :py:class:`subprocess.Popen` and 
    returns the id of the created process.

    :param args:        Argument passed to :py:class:`subprocess.Popen`.
    """

    os.environ.update({'PYTHONPATH': ':'.join(sys.path)})
    sp = subprocess.Popen(args)
    return sp.pid

def deleteFilesUsingGlob(patterns):
    """
    Deletes file(s) matching patterns. Skips deleting directories.
    The function is currently not used anywhere in FUDGE so probably should be deleted.

    :param patterns:        Any objected that can be passed to :py:func:`glob.glob`.
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
