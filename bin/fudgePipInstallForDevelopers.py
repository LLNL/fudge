

import os
import argparse
import subprocess
import collections


fudgeDependencies = collections.OrderedDict([
  ('crossSectionAdjustForHeatedTarget', {
      'hasSubmodules': False, 'tag': None,
      'url': r'git+ssh://git@czgitlab.llnl.gov:7999/nuclear/fudge/crosssectionadjustforheatedtarget.git'}),
  ('Merced', {'hasSubmodules': False, 'tag': None,
              'url': r'git+ssh://git@czgitlab.llnl.gov:7999/nuclear/fudge/merced.git'}),
  ('pqu', {'hasSubmodules': False, 'tag': None, 'url': r'git+ssh://git@czgitlab.llnl.gov:7999/nuclear/common/pqu.git'}),
  ('brownies', {'hasSubmodules': False, 'tag': None,
                'url': r'git+ssh://git@czgitlab.llnl.gov:7999/nuclear/fudge/brownies.git'}),
  ('numericalFunctions', {'hasSubmodules': True, 'tag': 'fudge4.3-rc1',
                          'url': r'git+ssh://git@czgitlab.llnl.gov:7999/nuclear/common/numericalFunctions.git'}),
  ('xData', {'hasSubmodules': False, 'tag': None,
             'url': r'git+ssh://git@czgitlab.llnl.gov:7999/nuclear/common/xData.git'}),
  ('PoPs', {'hasSubmodules': False, 'tag': None, 'url': r'git+ssh://git@czgitlab.llnl.gov:7999/nuclear/pops/PoPs.git'})
])


def runShellCommand(_shellCommand, returnSTDOUT=False):
    process = subprocess.Popen(_shellCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = [output.decode('utf-8') for output in process.communicate()]

    if process.returncode and stderr:
        print(stderr)

    if returnSTDOUT:
        return stdout


def cloneRepository(url, _localPath, hasSubmodules=False):
    _shellCommand = ['git', 'clone', '-q']

    if hasSubmodules:
        _shellCommand.append('--recurse-submodules')

    _shellCommand += [url, _localPath]
    runShellCommand(_shellCommand)


def installFudgeModule(_localPath, _pythonInterpreter, repoTag, hasSubmodules):
    if repoTag is not None:
        workingFolder = os.getcwd()
        os.chdir(_localPath)
        _shellCommand = ['git', 'checkout', repoTag]
        runShellCommand(_shellCommand)

        if hasSubmodules:
            runShellCommand(['git', 'submodule', 'init'])
            runShellCommand(['git', 'submodule', 'update'])

        os.chdir(workingFolder)

    _shellCommand = [_pythonInterpreter, '-m', 'pip', 'install', '--upgrade', '--force-reinstall', '--no-deps', '-e',
                     _localPath]

    print('Executing: %s' % ' '.join(_shellCommand))
    runShellCommand(_shellCommand)


def createVirtualEnvironment(venvPath, _pythonInterpreter):
    # delete virtual environment if it already exists
    if os.path.isdir(venvPath):
        deleteVirtualEnvironment(venvPath)

    # create virtual environment
    _shellCommand = [_pythonInterpreter, '-m', 'venv', venvPath]
    runShellCommand(_shellCommand)

    # upgrade to the latest version of pip
    venvPythonInterprator = os.path.join(venvPath, 'bin/python')
    _shellCommand = [venvPythonInterprator, '-m', 'pip', 'install', '--upgrade', 'pip']
    runShellCommand(_shellCommand)

    # install numpy and matplotlib
    for pythonPackage in ['numpy', 'matplotlib']:
        _shellCommand = [venvPythonInterprator, '-m', 'pip', 'install', pythonPackage]
        runShellCommand(_shellCommand)


def deleteVirtualEnvironment(venvPath):
    _shellCommand = ['rm', '-rf', venvPath]
    runShellCommand(_shellCommand)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Pip installation script for FUDGE developers')
    parser.add_argument('topLevelFudge', type=str, help='Root folder for FUDGE installation')
    parser.add_argument('--fudgeRepoURL', type=str, help='URL to the FUDGE repo', 
                        default=r'git+ssh://git@czgitlab.llnl.gov:7999/nuclear/fudge/fudge.git')
    parser.add_argument('--fudgeRepoTag', metavar='fudgeRepoTag', type=str, default=None,
                        help='FUDGE Repo tag to use with `pip install`')
    parser.add_argument('--virtualEnvironmentPath', metavar='virtualEnvironmentPath', type=str, default=None,
                        help='Option to create a new python virtual environment')
    parser.add_argument('--pythonInterpreter', metavar='pythonInterpreter', type=str, default='python3',
                        help='Python interpreter name')
    parser.add_argument('--printArguments', default=False, action='store_true',
                        help='Option to print the input arguments')
    parser.add_argument('--cloneOnly', default=False, action='store_true', help='Option to clone without installing')

    args = parser.parse_args()
    if args.printArguments:
        for arg in vars(args):
            print(arg, '=', getattr(args, arg))

    pipInstallAfterCloning = not args.cloneOnly

    # create and activate new python virtual environment if necessary
    if args.virtualEnvironmentPath is not None:
        createVirtualEnvironment(args.virtualEnvironmentPath, args.pythonInterpreter)
        pythonInterpreter = os.path.join(args.virtualEnvironmentPath, 'bin/python')

    else:
        # default pythonInterpreter
        pythonInterpreter = args.pythonInterpreter

    # clone main fudge repository
    topLevelFudge = args.topLevelFudge
    cloneRepository(args.fudgeRepoURL, topLevelFudge)

    # clone and pip install fudge dependencies
    for packageName in fudgeDependencies.keys():
        localPath = os.path.join(topLevelFudge, packageName)
        cloneRepository(fudgeDependencies[packageName]['url'], localPath,
                        fudgeDependencies[packageName]['hasSubmodules'])

        if pipInstallAfterCloning:
            installFudgeModule(localPath, pythonInterpreter, fudgeDependencies[packageName]['tag'],
                               fudgeDependencies[packageName]['hasSubmodules'])

    # create symbolic link for statusMessageReporting
    shellCommand = ['ln', '-s', os.path.join(topLevelFudge, 'numericalFunctions/statusMessageReporting'),
                    os.path.join(topLevelFudge, 'statusMessageReporting')]
    runShellCommand(shellCommand)

    # finally install FUDGE
    if pipInstallAfterCloning:    
        installFudgeModule(topLevelFudge, pythonInterpreter, args.fudgeRepoTag, None)
