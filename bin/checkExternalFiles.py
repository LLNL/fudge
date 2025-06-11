#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse
import pathlib

from LUPY import argumentsForScripts as argumentsForScriptsModule
from fudge import GNDS_file as GNDS_fileModule
from fudge import suites as suitesModule
from fudge.covariances import covarianceSuite as covarianceSuiteModule

summaryDocString__FUDGE = """Checks the externalFiles for consistency in one or more protares file."""

description = """
This scripts checks the externalFiles for consistency in one or more protares file. The path argument can refer to a map file,
or one or more protare files.
"""

parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('path', type=pathlib.Path, nargs='+',               help='The path to a GNDS map or protare file, or a list of protare files.')
parser.add_argument('-v', '--verbose', action='count', default=0,       help='The more the more verbosity.')

args = parser.parse_args( )

protarePath = None

def printMessage(message):
    """
    This function prints the protare's path if it has not be printed and then prints the message indented.

    :param message:     The message to print.
    """

    global protarePath

    if protarePath is not None:
        print(protarePath)
    protarePath = None

    if message is not None:
        print('    %s' % message)

def checkCovarianceExternals(covariancePath, protareRealPath):
    """
    This function checks the external files in a covarianceSuite.

    :param covariancePath:      The path to the covariance file to check.
    :param protareRealPath:     The real path to the protare file which contains an external reference to the covariance file.
    """

    covariance = GNDS_fileModule.preview(covariancePath, haltParsingMoniker=suitesModule.ExternalFiles.moniker)
    matchFound = False
    for externalFile in covariance.externalFiles:
        realpath = pathlib.Path(externalFile.realpath())
        try:
            if realpath.samefile(protareRealPath):
                matchFound = True
                break
        except FileNotFoundError:
            printMessage('Path for extenalFile in corariance does not exists: %s' % realpath)

    if not matchFound:
        printMessage('No external file in corariance pointing to a protare file: %s' % covariancePath)

def checkProtare(path):
    """
    This function checks all externalFile's in the file *path* and the externalFiles in their paths.

    :param path:    The path to the protare to check.
    """

    global protarePath

    protarePath = path
    if args.verbose > 0:
        printMessage(None)

    if pathlib.Path(path).exists():
        protare = GNDS_fileModule.preview(path, haltParsingMoniker=suitesModule.ExternalFiles.moniker)

        for externalFile in protare.externalFiles:
            realpath = externalFile.realpath()
            if pathlib.Path(realpath).exists():
                name, _ = GNDS_fileModule.type(realpath)
                if name == GNDS_fileModule.HDF5_values:
                    pass
                elif name == covarianceSuiteModule.CovarianceSuite.moniker:
                    checkCovarianceExternals(realpath, path)
                else:
                    printMessage('Unsupported external file type: %s' % name)
            else:
                printMessage('External file not found: %s' % realpath)
    else:
        printMessage('Protare not found.')

map1 = argumentsForScriptsModule.mapFromMapOrProtarePath(args.path)
for mapProtare in map1.iterate():
    checkProtare(mapProtare.fileName)
