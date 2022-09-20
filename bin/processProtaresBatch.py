#! /usr/bin/env python3
 
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os, argparse, sys, datetime, math, copy, glob, json, re

from LUPY import subprocessing

from fudge import enums as enumsModule
from fudge import GNDS_file as GNDS_fileModule

""" 
Create batch scripts for all gnds processing tasks, mcf, ndf, ndf-heated, tdf and create one installable tarball of the resulting library 
    sample test run to process multigroup at 2 temperatures: 
        batchProcessProtares.py -library endl2009.2-rc2  -t 1.e-3 -t 1.e-4 -d MG
    sample small set run:  (using --query to check your list of files to process)
        batchProcessProtares.py --library endl2009.3-rc3 -t 1e-3 -d 'MG' --filelist './gnds/endl/za005*.xml' --tag 'testMG' --query -vv
"""
 
### environment specs
machineDefault = os.environ['LCSCHEDCLUSTER']
bankDefault = 'wbronze'
if os.environ['PWD'].find('lscratch') > 0 : 
    diskDefault = os.environ['PWD'].split(os.path.sep)[2]
else :
    diskDefault = ''
timeDefault = '24'
fudgeDefault = '/usr/apps/fudge/current/'
timestamp = str( datetime.datetime.now() ).replace(' ','_')
jobsDefault = 10

### processing specs
argumentsDefault = '-t 2.586e-8 -mc -mg -up --tag proc'
outputDefault = 'gnds'
filelistDefault = 'gnds/direct/*.xml'

parser = argparse.ArgumentParser( )
parser.add_argument( "library", type=str,                                                               help = "where to find the gnds files, (full path)" )
parser.add_argument( "--filelist", type=str, action='append', default = None,                           help = "choose a subset of files for testing, stated as a shell expandable regex. Default = %s"%filelistDefault  )
parser.add_argument( '-a', '--arguments', default = argumentsDefault,                                   help = 'Arguments to pass through to processProtare.py. Default is %s' % argumentsDefault )
parser.add_argument( '-f', '--argumentFiles', default = None, action='append',                          help = '''
    Alernative option to pass argument input files (instead of explicit arguments) to processProtare.py.
    Input is expected to be in the form "-f projectile=file_path" for each projecile,
      e.g. "-f n=/usr/gapps/data/nuclear/common/processProtare/processProtare.n.input -f TNSL=/usr/gapps/data/nuclear/common/processProtare/processProtare.TNSL.full.input"''')
parser.add_argument( '-w', '--writeArgsFile', action='append', default=None,                            help = '''
    Option to write processProtare.py arguments to a file for record keeping and rerunning processPortare.py
    Input is expected to be in the form "-w projectile=file_path" for each projectile,
      e.g. "-w n=processProtare.n.input -w TNSL=processProtare.TNSL.full.input"''')

parser.add_argument( "-n", "--nodes", type=int, default=1,                                              help = "Allow n nodes at once" )
parser.add_argument( "-nt", "--tasks", type=int, default=8,                                             help = "Allow nt tasks per node" )
parser.add_argument( "-c", "--cores", type=int, default=None,                                           help = "Allow c cores per task for OMP threading. Default is number of cores per node on the chosen machine " )
parser.add_argument( "-p", "--partition", type=str, default = machineDefault,                           help = "LC machine to use, default = %s"%machineDefault )
parser.add_argument( "-b", "--bank", type=str, default = bankDefault,                                   help = "bank on LC machine, default = %s"%bankDefault )
parser.add_argument( "-k", "--disk", type=str, default = diskDefault,                                   help = "lustre disk on LC machine, default = %s"%diskDefault )
parser.add_argument( "-hr", "--hours", type=float, default = timeDefault,                               help = "maximum job length, default = %s"%timeDefault )
parser.add_argument( "-q", "--queue", type=str, default = 'pbatch',                                     help = "which queue, pdebug or pbatch?" )
parser.add_argument( "-j", '--jobs', type=int, default = jobsDefault,                                   help = "maximum number of simultaneous jobs submitted to the queue. Default = %d"%jobsDefault)
                                                                  
parser.add_argument( "-v", "--verbose", action = "count", default = 0,                                  help = "enable verbose output" )
parser.add_argument( "--explicitPython3", action="store_true", default = False,                         help = "explicitely call python3 instead of just python" )
parser.add_argument( "--rerun", action="store_true", default = False,                                   help = "force a rerun of all jobs" )
parser.add_argument( "--frontend", action="store_true",                                                 help = "run everything on the frontend commandline. For testing or quick small jobs. " )
parser.add_argument( "--dryrun", action="store_true", default = False,                                  help = "build scripts but dont submit. for testing. " )
parser.add_argument( "--query", action="store_true",                                                    help = "dont build scripts. just list files not done. for testing. " )
parser.add_argument( "--fudgePath", dest="fudgePath", type=str,  default=fudgeDefault,                  help = "what fudge do you want to use? default = %s"%fudgeDefault )
parser.add_argument("--venvActivateScript", default=None,                                               help = "Path to virtual environment activate script")
args = parser.parse_args( )

cwd = os.getcwd()

pythonCommand = 'python3' if args.explicitPython3 else sys.executable
recordString = ' '.join(['### %s' % pythonCommand]+sys.argv)

try:
    tagsIndex = args.arguments.split().index('--tag')+1
    args.newDir = args.arguments.split()[tagsIndex].strip('."')
except:
    args.newDir = 'proc'
args.tag = args.newDir

outputArgs = {}

def main():
      
    if args.frontend: args.partition = os.environ['HOSTNAME'].strip('0123456789')
    if args.query and args.verbose < 1 : args.verbose = 1

    ### check for valid fudge
    try:
        sys.path.insert(0,args.fudgePath)
        import fudge
        args.fudgePath = os.path.dirname(os.path.dirname(fudge.__file__))
        if( args.verbose > 0 ) : print('path to fudge : ', fudge.__file__, '\n', args.fudgePath)
    except:
        raise("don't recognize your fudge path! Is it up-to-date? and rebuilt?")

    ### setup any working directories for stashing batch and log files
    os.chdir(args.library)
    if not os.path.exists('batch') : os.mkdir('batch')
    if not os.path.exists('logs') : os.mkdir('logs')
    
    ### find all files specified by filelist argument    
    unFinishedCount, newZAargs = 0, []
    if args.filelist == None: args.filelist = [filelistDefault] # [os.path.join(args.library,filelistDefault)] 
    for zaregex in args.filelist: 
        if( args.verbose > 1 ) : print(zaregex)
        newZAargs.extend(glob.glob(str(zaregex)))
    args.filelist = newZAargs
    if( args.verbose > 0 ) : print('number of files %d ' % len(args.filelist))

    evaluatedPath = os.path.dirname(os.path.dirname(args.filelist[0]))
    args.newDir = os.path.join(evaluatedPath,args.newDir)
    if not os.path.exists(args.newDir) : os.makedirs(args.newDir)
    if( args.verbose > 0 ) : print('processed path : %s' % args.newDir)

    ### exclude files as specified by options
    for zaregex in copy.copy(args.filelist): 
        fpath, fname = os.path.split(zaregex)
        if fpath.find('wkdr_') >=0 :            ### dont include the working dir copies
            args.filelist.remove(zaregex)
            continue

        fullfile = os.path.join(args.newDir, fname.replace('.xml','.%s.xml'%args.tag) )
        if os.path.exists(fullfile) and os.stat(fullfile).st_size != 0 and not args.rerun:   ### rerun will force jobs to be redone
            args.filelist.remove(zaregex)
            if( args.verbose > 3 ) : print('skipping this one %s' % fullfile)
        else:
            if( args.verbose > 2 ) : print('doing this one %s' % fullfile)
            unFinishedCount += 1

    if( args.verbose > 0 ) : print('jobs left : %d'%unFinishedCount)
    if args.query : return 
    args.filelist.sort()

    ### remove heating and upscatter from non-neutron arguments
    if args.argumentFiles is None:
        useExplicitArgs = True
        argumentsNonNeutron = args.arguments.split()
        while '-t' in  argumentsNonNeutron :
            aInd = argumentsNonNeutron.index('-t')
            del argumentsNonNeutron[aInd:aInd+2]
        if '-up' in argumentsNonNeutron: argumentsNonNeutron.remove('-up')
        argumentsNonNeutron = ' '.join(argumentsNonNeutron)

    else:
        useExplicitArgs = False
        argFileDict = {}
        argFileRegex = re.compile(r'^([^=]+)=([^=]+)$')
        for singleArgumentFile in args.argumentFiles:
            assert argFileRegex.match(singleArgumentFile), \
            'Arguments for the processProtare.py inputfile arguments are expected to be of the form %s!!!' % argFileRegex.pattern
            projectile, inputFile = argFileRegex.findall(singleArgumentFile)[0]
            argFileDict[projectile] = inputFile

    if args.writeArgsFile:
        tempOutputArgs = {}
        argFileRegex = re.compile(r'^([^=]+)=([^=]+)$')
        for singleArgumentFile in args.writeArgsFile:
            assert argFileRegex.match(singleArgumentFile), \
                'Arguments for the "--writeArgsFile" option are expected to be of the form "-w %s"!!!' % argFileRegex.pattern
            projectile, argumentFile = argFileRegex.findall(singleArgumentFile)[0]
            tempOutputArgs[projectile] = [argumentFile]

    
    ### split the fileNames List into two different lists for neutrons and non neutrons
    fileListNeutrons = []
    fileListOther = []
    fileListTNSL = []
    fileListLLNL_TNSL = []
    individualFileArguments = {}
    for filename in args.filelist:
        try : 
            passedFileName, d = GNDS_fileModule.type(filename)
        except : 
            print(filename)
            raise
        if 'projectile' in d :
            if d['projectile'] == 'n': 
                # Look at interaction type
                if d['interaction'] == enumsModule.Interaction.TNSL:
                    fileListTNSL.append(filename)
                    if args.writeArgsFile and len(tempOutputArgs['TNSL']) == 1:
                        tempOutputArgs['TNSL'].append(filename)
                elif d['interaction'] == enumsModule.Interaction.LLNL_TNSL:
                    fileListLLNL_TNSL.append(filename)
                    if args.writeArgsFile and len(tempOutputArgs['LLNL_TNSL']) == 1:
                        tempOutputArgs['LLNL_TNSL'].append(filename)
                else:
                    fileListNeutrons.append(filename)
                    if args.writeArgsFile and len(tempOutputArgs['n']) == 1:
                        tempOutputArgs['n'].append(filename)
            else:
                fileListOther.append(filename)
                if args.argumentFiles is not None:
                    individualFileArguments[filename] = '@%s' % argFileDict[d['projectile']]
                    if args.writeArgsFile and len(tempOutputArgs[d['projectile']]) == 1:
                        tempOutputArgs[d['projectile']].append(filename)
                        individualFileArguments[filename] += ' --printArgsFile %s' % os.path.join(args.library, tempOutputArgs[d['projectile']][0])
                elif args.writeArgsFile and len(tempOutputArgs[d['projectile']]) == 1:
                    tempOutputArgs['n'].append(filename)
    
    if args.writeArgsFile:
        for values in tempOutputArgs.values():
            if len(values) == 2:
                outputArgs[values[1]] = values[0]
            
    os.chdir(cwd)

    if args.cores: 
        args.procstr = '--NumProcesses %s'%args.cores
    else: 
        args.procstr = ''
    args.cores,args.numhours,args.numnodes = getQueueLims(
            args.partition,args.nodes,args.hours,procs=args.cores)      ### look up job limits on requested machine
    
    batsubname = '%s/batch/batchsubmit_%s.sh'%(args.library,os.path.split(args.library)[1])
    batsublines = []
    batsublines.append('#!/bin/env bash \n' )
    batsublines.append('### %s \n'%recordString)

    def neutronOrTNSLArgs(_projectileAlias, _filelist):
        if args.writeArgsFile:
            if tempOutputArgs[_projectileAlias][1] in _filelist:
                _useIndividualFileArguments = True
                for _filename in _filelist:
                    individualFileArguments[_filename] = args.arguments if useExplicitArgs else "@%s" % argFileDict[_projectileAlias]
                    if _filename == tempOutputArgs[_projectileAlias][1]:
                        individualFileArguments[_filename] += ' --printArgsFile %s' % os.path.join(args.library, tempOutputArgs[_projectileAlias][0])

                _arguments = None

            else:
                _useIndividualFileArguments = False
                _arguments = args.arguments if useExplicitArgs else "@%s" % argFileDict[_projectileAlias]

        else:
            _useIndividualFileArguments = False
            _arguments = args.arguments if useExplicitArgs else "@%s" % argFileDict[_projectileAlias]

        return _useIndividualFileArguments, _arguments


    #print('non-neutron list: ',len(fileListOther),fileListOther)
    for listNumber,filelist in enumerate([fileListNeutrons, fileListTNSL, fileListLLNL_TNSL, fileListOther]) :
        if len(filelist) == 0:
            continue
        
        if listNumber==0:
            batchNamePrefix = 'neutron'
            useIndividualFileArguments, arguments = neutronOrTNSLArgs('n', filelist)

        elif listNumber==1:
            batchNamePrefix = 'TNSL'
            useIndividualFileArguments, arguments = neutronOrTNSLArgs('TNSL', filelist)

        elif listNumber==2:
            batchNamePrefix = 'LLNL_TNSL'
            useIndividualFileArguments, arguments = neutronOrTNSLArgs('LLNL_TNSL', filelist)

        else:
            batchNamePrefix = 'other'
            if useExplicitArgs:
                useIndividualFileArguments = False
                arguments = argumentsNonNeutron
            else:
                useIndividualFileArguments = True
                arguments = None

        batchListDict = {}
        batchNames = []
        numFiles = len(filelist)
        if numFiles < 1 : 
            print('no files for %s'%batchNamePrefix)
            continue
        numJobs = args.jobs
        numTasks = args.tasks
        numLinks = (numFiles / float(numJobs*numTasks) )
        if numLinks < 1 :
            if args.verbose > 3: print('too few files for this many jobs, adjusting the striping:')
            if args.verbose > 3: print('Files', 'Jobs', 'Tasks', 'Links')
            if args.verbose > 3: print(numFiles, numJobs, numTasks, numLinks)
        while numLinks < 1 :
            if numTasks > 1 : 
                numTasks -= 1           
            if numJobs > 1 : 
                numJobs -= 1
            numLinks = int(numFiles / float(numJobs*numTasks) )
            if args.verbose > 3: print(numFiles, numJobs, numTasks, numLinks)
            
        ### walk the list of processing files and stripe them every job*tasks
        filenumber = -1 
        for filename in filelist: 
            #print(filename, os.path.exists(filename))
            #if os.path.exists(filename) and not args.rerun : continue
            
            filenumber += 1
            jobname = '%s%03d'%(batchNamePrefix,filenumber % (numJobs*numLinks))
            if jobname not in batchListDict : batchListDict[jobname] = []
            batchListDict[jobname].append(filename)
            
        ### write the individual batch files and store the names
        for batchNum,batchName in enumerate(batchListDict.keys()) :
            if useIndividualFileArguments:
                batchFolder = os.path.split(batsubname)[0]
                if not os.path.isdir(batchFolder):
                    os.mkdir(batchFolder)

                argumentsJSON = ('%s/%s.json' % (batchFolder, batchName), individualFileArguments)
                # fileObject = open(argumentsJSON, 'w')
                # json.dump(individualFileArguments, fileObject)
                # fileObject.close()

            else:
                argumentsJSON = None

            batchNames.append( getSeparateBatch(batchListDict[batchName],arguments,batchName,argumentsJSON) )
            #print(batchNames[-1])
            
        ### now stripe the chain submissions    
        for batchNum,batchname in enumerate(batchNames) :
            if args.frontend:
                batsublines.append('source %s \n' % (batchname))
            else:
                if batchNum//numJobs == 0:
                    batsublines.append('msub %s \n' % (batchname))
                else:
                    with open(batchNames[batchNum-numJobs], 'a') as fout:
                        fout.write( 'msub %s \n' % (batchname) )    

    with open(batsubname,'w') as batsub :
        batsub.writelines(batsublines)
    
    if not args.dryrun:  
        os.system('chmod u+x %s'%batsubname )
        os.system('. %s'%batsubname )
    else:
        print('submit to batch with ==> source %s'%batsubname)

def getSeparateBatch(filenames,arguments,batchName,argumentsJSON=None):
    thisnodes,thispart,thishours,thisprocs = args.numnodes,args.partition,args.numhours,args.cores
    workstrj = '%s_%s_%s_n%d'%(os.path.basename(args.library),os.path.basename(batchName),thispart,thisnodes)
    if( args.verbose > 0 ) : print('writing batch file %s' % workstrj)

    jout = open('%s/batch/batch_%s.sh' % (os.path.basename(args.library),workstrj), 'w') 
    thisname = jout.name
    print(thisname)
    jout.write( '#!/bin/bash \n')
    jout.write( '##Moab options \n')
    jout.write( '#MSUB -V \n')
    jout.write( '#MSUB -N %s \n'%workstrj )
    jout.write( '#MSUB -l partition=%s \n'%thispart )
    jout.write( '#MSUB -l walltime=%s \n'%thishours )
    jout.write( '#MSUB -l nodes=%d \n'%thisnodes )
    if len(args.disk)>0 : jout.write( '#MSUB -l gres=%s \n'%args.disk )
    if thispart in ('borax','boraxo','rztrona'): jout.write( '#MSUB -l ttc=%s \n'%thisprocs )
    jout.write( '#MSUB -A %s \n'%args.bank )
    jout.write( '#MSUB -q %s \n'%args.queue )
    jout.write( '#MSUB -j oe \n')
    jout.write( '#MSUB -o %s/logs/log_job_%s.out \n'%(args.library, batchName+'_'+timestamp) )
    if args.venvActivateScript is not None : jout.write('\nsource %s\n\n' % args.venvActivateScript)
    if not args.frontend : jout.write( 'echo "job_id = $SLURM_JOBID" \n')
    systemBinPath = os.path.join(sys.prefix, 'bin')
    fudgeBinPath = systemBinPath if os.path.exists(os.path.join(systemBinPath, 'processProtare.py')) else os.path.join(args.fudgePath, 'bin')
    assert os.path.exists(os.path.join(fudgeBinPath, 'processProtare.py')), 'The processProtare.py script was not found'
    jout.write( 'OMP_NUM_THREADS=%d \n'%thisprocs)
    jout.write( 'export PYTHONPATH=%s/:%s/:$PYTHONPATH\n' % (args.fudgePath, fudgeBinPath) )
    jout.write( 'export PATH=%s:$PATH\n'%fudgeBinPath )
    jout.write( 'CWD=$PWD \n' )
    if not args.frontend : jout.write( 'echo $PWD \n' )
    jout.write( 'echo "starting job : " ; date \n' )
    filelist = []
    jsonArgumentDictionary = {}
    for filename in filenames:
        fileBase = os.path.basename(filename)
        filePath = '%s/%s/working_%s/wkdr_%s'%(args.library,os.path.dirname(filename),args.tag,fileBase)
        jout.write( 'mkdir -p %s \n'%(filePath) )
        jout.write( 'cp %s/%s  %s/ \n'%(args.library,filename,filePath) )
        filelist.append( '%s/%s'%(filePath,fileBase) )

        if argumentsJSON is not None:
            jsonArgumentDictionary[filelist[-1]] = argumentsJSON[1][filename]

    # FIXME hard-coding time limits for now, should compute from walltime.  Also, timeout may not be installed on all systems:
    if argumentsJSON is None:
        jout.write( 'timeout -k 86100s 85800s %s %s %s %s -p %s -c --options "%s  > runLog_%s.txt 2>&1"   \n' %
                    (pythonCommand, os.path.join(fudgeBinPath, 'multiprocessGeneral.py'),
                     os.path.join(fudgeBinPath, 'processProtare.py'),
                     ' '.join(filelist), pythonCommand, ''.join(arguments), args.tag) )
    else:
        with open(argumentsJSON[0], 'w') as fileObject:
            json.dump(jsonArgumentDictionary, fileObject)

        jout.write( 'timeout -k 86100s 85800s %s %s %s %s -p %s -c --optionsJSON %s --options "  > runLog_%s.txt 2>&1"   \n' %
                    (pythonCommand, os.path.join(fudgeBinPath, 'multiprocessGeneral.py'),
                     os.path.join(fudgeBinPath, 'processProtare.py'),
                     ' '.join(filelist), pythonCommand, argumentsJSON[0], args.tag) )

    for filename in filenames:
        fileBase = os.path.basename(filename)
        filePath = '%s/%s/working_%s/wkdr_%s'%(args.library,os.path.dirname(filename),args.tag,fileBase)
        filetype = fileBase.split('.')[-1]
        processedFileName = fileBase[:-len(filetype)] + '%s.%s' % (args.tag, filetype)
        jout.write('ln %s/%s %s/%s/  \n' % (filePath, processedFileName, args.library, args.newDir))
    
    jout.close()
    return thisname
    
def getQueueLims(partition,nodes,hours,procs=None):
    def readQueueLimits(_machineDetails):
        _maximumTime = _machineDetails['partitions']['pbatch']['maxtime']
        if re.match('^(\d+-)*\d+:\d+:\d+$', _maximumTime):
            _maximumTime = datetime.datetime.strptime(_maximumTime, '%d-%H:%M:%S') if re.match('^\d+-\d+:\d+:\d+$', _maximumTime) else datetime.datetime.strptime(_maximumTime, '%H:%M:%S')
            _maximumHours = _maximumTime.day*24 + _maximumTime.hour + _maximumTime.minute/60
        else:
            _maximumHours = _maximumTime
    
        _numberProcesses = _machineDetails['partitions']['pbatch']['cpuspernode']
        _maximumNodes = _machineDetails['partitions']['pbatch']['maxnodes']
        if _maximumNodes == 'UNLIMITED':
            _maximumNodes = _machineDetails['partitions']['pbatch']['totalnodes']
    
        return _maximumHours, int(_numberProcesses), int(_maximumNodes)
    
    queueLims = {
            'borax':    (200., 36, 1),
            'boraxo':   (200., 36, 1),
            'cab':      (16., 16, 1),
            'syrah':    (16., 16, 128),
            'quartz':   (24., 36, 256),
            'ruby':     (24., 40, 520),
            'rzmerl':   (24., 16, 32),
            'rztopaz':  (24., 36, 64),
            'rztrona':  (168., 36, 1),
            'agate':    (200., 36, 1),
            'jade':     (24., 36, 256),
            'jadeita':  (24., 36, 256),
            'magma':    (24., 96, 256),
            'mica':     (8., 36, 290),
            }

    if queueLims.get(partition):
        numhours, numprocs, numnodes = queueLims[ partition ]
    else:
        _stdout = subprocessing.executeCommand(['/usr/global/tools/lorenz/bin/lora', '-k', '/clusters/batchdetails'])[1]
        _batchDetails = json.loads(''.join(list(_stdout)))

        assert _batchDetails['output'].get(partition), 'No information available for machine/partition %s' % partition       
        numhours, numprocs, numnodes = readQueueLimits(_batchDetails['output'][partition])

        if numhours == 'UNLIMITED':
            numhours = hours
        else:
            numhours = float(numhours)

    if nodes>numnodes:
        if( args.verbose > 1 ) : print('warning: too many nodes requested, using %d'%numnodes)
        nodes = numnodes
    if hours>numhours:
        if( args.verbose > 1 ) : print('warning: too many hours requested. using %d'%numhours)
        hours = numhours
    
    fracPart, intPart = math.modf(hours)
    hourStr = '%d:%02i:00'%(intPart,math.floor(fracPart*60))
    
    if procs: numprocs = procs  ### user overrides number of machine procs
    
    return numprocs,hourStr,nodes


if __name__=='__main__':
    main()

