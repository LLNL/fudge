# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import re
import os
import sys
import subprocess


versionOutputFile = os.path.join('fudge', 'fudgeVersion.py')


def runShellCommand(_shellCommand, returnSTDOUT=False):
    process = subprocess.Popen(_shellCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = [output.decode('utf-8') for output in process.communicate()]

    if process.returncode and stderr:
        print(stderr)

    if returnSTDOUT:
        return stdout


def defaultVersionNumber():
    if not os.path.isfile(versionOutputFile):
        output = ['FUDGE_MAJORVERSION = 6', 'FUDGE_MINORVERSION = 4', 'FUDGE_RELEASECANDIDATE = \'\'', 'FUDGE_POSTRELEASE = \'\'', 'FUDGE_REPOIDENTIFIER = \'\'']

        with open(versionOutputFile, 'w') as fileObject:
            fileObject.write('\n'.join(output))


def getVersionNumber():
    if os.path.isdir('.git'):
        remoteURL = runShellCommand(['git', 'config', '--get', 'remote.origin.url'], True).split('\n')[0]
        if 'zgitlab.llnl.gov' in remoteURL:
            gitDescribe = runShellCommand(['git', 'describe', '--abbrev=40', '--match', 'fudge[0-9].[0-9]*'], True)
            
            # see https://www.python.org/dev/peps/pep-0440/
            regexGitDescribe = re.compile('^fudge(?P<major>\d+).(?P<minor>\d+)(?:-(?P<releaseCondidate>rc\d+)(-(?P<postRelease>\d+)-(?P<repoIdentifier>.+))?)?$')
            regexMatch = regexGitDescribe.match(gitDescribe)
            if regexMatch:
                regexMatch = regexMatch.groupdict()

                if not regexMatch["repoIdentifier"]:
                    gitDescribe = runShellCommand(['git', 'rev-parse', 'HEAD'], True)
                    repoIdentifier = f'g{gitDescribe[:8]}'
                else:
                    repoIdentifier = regexMatch["repoIdentifier"]

                output = [f'FUDGE_MAJORVERSION = {regexMatch["major"]}', 
                          f'FUDGE_MINORVERSION = {regexMatch["minor"]}',
                          f'FUDGE_RELEASECANDIDATE = \'{regexMatch["releaseCondidate"] if regexMatch.get("releaseCondidate") else ""}\'',
                          f'FUDGE_POSTRELEASE = \'{regexMatch["postRelease"] if regexMatch.get("postRelease") else ""}\'',
                          f'FUDGE_REPOIDENTIFIER = \'{repoIdentifier}\'']

                with open(versionOutputFile, 'w') as fileObject:
                    fileObject.write('\n'.join(output))

            else:
                defaultVersionNumber()

        else:
            defaultVersionNumber()

    else:
        defaultVersionNumber()

    assert os.path.isfile(versionOutputFile)

    sys.path.append('fudge')
    import fudgeVersion as versionModule

    versionString = f'{versionModule.FUDGE_MAJORVERSION}.{versionModule.FUDGE_MINORVERSION}'
    if versionModule.FUDGE_RELEASECANDIDATE != '':
        versionString += f'{versionModule.FUDGE_RELEASECANDIDATE}post{versionModule.FUDGE_POSTRELEASE}'

    return versionString


if __name__ == '__main__':
    getVersionNumber()
