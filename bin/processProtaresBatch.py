#!/bin/env python
 
# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
# <<END-copyright>>

import os, argparse, sys, datetime, math, copy, re, glob

""" 
Create batch scripts for all gnds processing tasks, mcf, ndf, ndf-heated, tdf and create one installable tarball of the resulting library 
    sample test run to process multigroup at 2 temperatures: 
        batchProcessProtares.py -library endl2009.2-rc2  -t 1.e-3 -t 1.e-4 -d MG
    sample small set run:  (using --query to check your list of files to process)
        batchProcessProtares.py --library endl2009.3-rc3 -t 1e-3 -d 'MG' --filelist './gnds/endl/za005*.xml' --tag 'testMG' --query -vv
"""

### environment specs
machineDefault = 'borax'
bankDefault = 'wbronze'
diskDefault = 'lscratche'
timeDefault = '24'
fudgeDefault = '/usr/apps/fudge/current/'
timestamp = str( datetime.datetime.now() ).replace(' ','_')
jobsDefault = 10

### processing specs

tempsDefault = 2.586e-8
energyUnitDefault = 'MeV'
tempUnitDefault = 'MeV/k'       ### k = 8.617e-5 eV/K
lMaxDefault = 9
fidDefault = 'LLNL_fid_1'
bdflsDefault = '/usr/gapps/data/nuclear/bdfls.archive/bdfls.Audi_etal.2003.12.22'
reconChoiceDefault = 'crossSection'
reconAccuracyDefault = 1e-6
tagDefault = 'proc'
outputDefault = 'gnds'
filelistDefault = 'gnds/direct/*.xml'

### default style prefixes
MGprefixDefault = 'MultiGroup'
MCprefixDefault = 'MonteCarlo'
UPprefixDefault = 'UpScatter'
HTprefixDefault = 'heated'
AEPprefixDefault = 'aep'
EVprefixDefault = 'eval'
RCprefixDefault = 'recon'

### Default group structures
gids = {'n'     : 'LLNL_gid_7',
        'H1'    : 'LLNL_gid_71', 
        'H2'    : 'LLNL_gid_71', 
        'H3'    : 'LLNL_gid_71', 
        'He3'   : 'LLNL_gid_71', 
        'He4'   : 'LLNL_gid_71', 
        'photon' : 'LLNL_gid_70' }
gidDefaults = ""
for gid in gids : gidDefaults += '%s="%s", ' % ( gid, gids[gid] )
gidDefaults = gidDefaults[:-2]   ### remove trailing comma

parser = argparse.ArgumentParser( )
parser.add_argument( "library", type=str,                                                             help = "where to find the gnds files, (full path)" )
parser.add_argument( "--filelist", type=str, action='append', default = None,                           help = "choose a subset of files for testing, stated as a shell expandable regex. Default = %s"%filelistDefault  )
parser.add_argument( "--tag", type=str, default=tagDefault,                                             help = "Any special tag you want to give the processed files. Default is %s " % tagDefault )
parser.add_argument( "--newDir", type=str, default=None,                                          help = "Any special tag you want to give the processed files. Default is %s " % tagDefault )
parser.add_argument( "-t", "--temperatures", type=float, action='append', default=None,                 help = "temperatures for heating, use successive -t arguments. Default is %s " % tempsDefault )
parser.add_argument( "-d", "--dataTypes", action='append', choices = ['all','MG','MC','UP'], default=[], help = "which data types to compute" )
parser.add_argument( '-o', '--output', default = outputDefault,                                         help = 'directory to write output file. If "SRC", writes to source directory. Default is %s' % outputDefault )
parser.add_argument( "--writeConvertedUnits", action = 'store_true', default = False,                   help = 'write data to file after units converted' )

parser.add_argument( "-n", "--nodes", type=int, default=1,                                              help = "Allow n nodes at once" )
parser.add_argument( "-nt", "--tasks", type=int, default=8,                                             help = "Allow nt tasks per node" )
parser.add_argument( "-c", "--cores", type=int, default=None,                                           help = "Allow c cores per task for OMP threading. Default is number of cores per node on the chosen machine " )
parser.add_argument( "-p", "--partition", type=str, default = machineDefault,                           help = "LC machine to use, default = %s"%machineDefault )
parser.add_argument( "-b", "--bank", type=str, default = bankDefault,                                   help = "bank on LC machine, default = %s"%bankDefault )
parser.add_argument( "-k", "--disk", type=str, default = diskDefault,                                   help = "lustre disk on LC machine, default = %s"%diskDefault )
parser.add_argument( "-hr", "--hours", type=float, default = timeDefault,                               help = "maximum job length, default = %s"%timeDefault )
parser.add_argument( "-q", "--queue", type=str, default = 'pbatch',                                     help = "which queue, pdebug or pbatch?" )
parser.add_argument( "-j", '--jobs', type=int, default = jobsDefault,                                   help = "maximum number of simultaneous jobs submitted to the queue. Default = %d"%jobsDefault)

parser.add_argument( "--energyUnit", type = str, default = energyUnitDefault,                           help = 'energy unit to convert to. Default is %s ' % energyUnitDefault)
parser.add_argument( "--temperatureUnit", type = str, default = tempUnitDefault,                        help = 'temperature unit to convert to, default is %s ' % tempUnitDefault )
parser.add_argument( "--bdfls", type = str, default = bdflsDefault,                                     help = "bdfls file to use. Default is %s " % bdflsDefault )
parser.add_argument( "--fluxID", type = str, default = fidDefault,                                      help = 'Flux ID from bdfls. Default is %s ' % fidDefault )
parser.add_argument( "-g", "--gid",  action='append', default = None,                                   help = "particle grouping schemes :  <particle>=LLNL_gid_<integer>. Default is %s." % gidDefaults )
parser.add_argument( "--legendreMax", type = int, default = lMaxDefault,                                help = "Maximum Legendre order for Sn prcessed data. Default is %s." % lMaxDefault )
parser.add_argument( "--reconAccuracy", type = float, default = reconAccuracyDefault,                   help = "Accuracy for reconstructing resonances. Default is %.1e" % reconAccuracyDefault )
parser.add_argument( "--reconstruct", type = str, choices = ['all','crossSection','angular'], default=reconChoiceDefault, 
                                                                                                        help = "What kind of reconstructing should be done ('all','crossSection','angular')" )

parser.add_argument( "-v", "--verbose", action = "count", default = 0,                                  help = "enable verbose output" )
parser.add_argument( "--frontend", action="store_true",                                                 help = "run everything on the frontend commandline" )
parser.add_argument( "--rerun", action="store_true", default = False,                                   help = "force a rerun of all jobs" )
parser.add_argument( "--dryrun", action="store_true", default = False,                                  help = "for testing. build scripts but dont submit. " )
parser.add_argument( "--query", action="store_true",                                                    help = "for testing. dont build scripts. just list files not done " )
parser.add_argument( "--elementals", action="store_true", default = False,                              help = "ignore elementals. Default is true " )
parser.add_argument( "--fudgePath", dest="fudgePath", type=str,  default=fudgeDefault,                  help = "what fudge do you want to use? default = %s"%fudgeDefault )
args = parser.parse_args( )

cwd = os.getcwd()

if ( args.temperatures == None ) : args.temperatures = [ tempsDefault ]
if ( args.newDir == None ) : args.newDir = args.tag

if not ( args.gid == None ) :
    for gidVal in args.gid:
        key, value = gidVal.split('=')
        if key not in gids: raise UserWarning("particle grouping specified for unknown particle name : %s " % (key) )
        gids[key] = value
gidDefaults = ""
for gid in gids : gidDefaults += '%s="%s", ' % ( gid, gids[gid] )
gidDefaults = gidDefaults[:-2]   ### remove trailing comma
 
def main():
      
    if args.frontend: args.partition = os.environ['HOSTNAME'].strip('0123456789')
    if args.query and args.verbose < 1 : args.verbose = 1

    ### check for valid fudge
    try:
        sys.path.insert(0,args.fudgePath)
        import fudge
        args.fudgePath = os.path.dirname(fudge.__file__)
        if( args.verbose > 0 ) : print 'path to fudge : ', fudge.__file__, '\n', args.fudgePath
    except:
        raise("don't recognize your fudge path! Is it up-to-date? and rebuilt?")

    ### setup any working directories for stashing batch and log files
    if not os.path.exists('batch') : os.mkdir('batch')
    if not os.path.exists('logs') : os.mkdir('logs')
    os.chdir(args.library)
    if not os.path.exists('logs') : os.mkdir('logs')
    if not os.path.exists(args.newDir) : os.makedirs(args.newDir)
    
    ### find all files specified by filelist argument    
    unFinishedCount, newZAargs = 0, []
    if args.filelist == None: args.filelist = [os.path.join(args.library,filelistDefault)] 
    print args.filelist
    for zaregex in args.filelist: 
        if( args.verbose > 1 ) : print zaregex
        newZAargs.extend(glob.glob(zaregex))
    args.filelist = newZAargs
    if( args.verbose > 0 ) : print 'number of files %d '%len(args.filelist)
    
    ### exclude files as specified by options
    for zaregex in copy.copy(args.filelist): 
        if zaregex.find('gnds') < 0:  ### ignore if file is not in gnds directory
            if( args.verbose > 2 ) : print 'ignoring non gnds file: %s '%zaregex
            args.filelist.remove(zaregex) 
            continue
                
        fpath, fname = os.path.split(zaregex)
        
        if fpath.find('wkdr_') >=0 :   
            args.filelist.remove(zaregex)
            continue
        
            
        fullfile = os.path.join(args.newDir, fname.replace('.xml','_%s.xml'%args.tag) )
        if os.path.exists(fullfile) or os.path.exists( zaregex.replace('.xml','.%s.xml'%args.tag) )  and not args.rerun:   ### rerun will force jobs to be redone
            args.filelist.remove(zaregex)
            if( args.verbose > 2 ) : print 'skipping this one %s' % fullfile 
        else:
            if( args.verbose > 1 ) : print 'doing this one %s' % fullfile 
            unFinishedCount += 1
                
    if args.query : 
        if( args.verbose > 0 ) : print 'jobs left : %d'%unFinishedCount
        return 
    args.filelist.sort()        

    os.chdir(cwd)

    if args.cores: 
        args.procstr = '--NumProcesses %s'%args.cores
    else: 
        args.procstr = ''
    args.cores,args.numhours,args.numnodes = getQueueLims(args.partition,args.nodes,args.hours,procs=args.cores)      ### look up job limits on requested machine

    arguments = []
    arguments.extend(["-t %s"%t for t in args.temperatures] )
    arguments.append("--tag %s"%args.tag )
    arguments.append("--fluxID %s"%args.fluxID )
    arguments.extend(['-g %s="%s"'%(gid, gids[gid]) for gid in gids] )
    arguments.append("--energyUnit %s"%args.energyUnit )
    arguments.append("--temperatureUnit %s"%args.temperatureUnit )
    arguments.append("--bdfls %s"%args.bdfls )
    if( outputDefault != args.output ) : arguments.append("--output %s"%args.output )
    arguments.append("--legendreMax %d"%args.legendreMax )
    if 'MG' in args.dataTypes: arguments.append("-mg" )
    if 'UP' in args.dataTypes: arguments.append("-mg -up" )
    if 'MC' in args.dataTypes: arguments.append("-mc" )
    if 'all' in args.dataTypes: arguments.append("-mc -mg" ) ### and -up when the style is added
    if args.writeConvertedUnits : arguments.append("--writeConvertedUnits" )
    arguments.append("--reconstruct %s"%args.reconstruct )
    arguments.append("--reconAccuracy %s"%reconAccuracyDefault )
    arguments.append("--prefixMultiGroup %s"%MGprefixDefault )
    arguments.append("--prefixMonteCarlo %s"%MCprefixDefault )
    arguments.append("--prefixUpScatter %s"%UPprefixDefault )
    arguments.append("--prefixHeated %s"%HTprefixDefault )
    arguments.append("--prefixAep %s"%AEPprefixDefault )
    arguments.append("--prefixEval %s"%EVprefixDefault )
    arguments.append("--prefixRecon %s"%RCprefixDefault )
    
   
    batsub = open('batchsubmit.sh','w')
    batsub.write('#!/bin/env bash \n' )
    
    filecounter = -1
    localFileList = []
    #print args.filelist
    for filenumber,filename in enumerate(args.filelist):        
        if os.path.exists(filename) and not args.rerun : continue
        
        ### split up the full list into jobs of size args.tasks
        localFileList.append(filename)
        if( ( (filenumber+1) % args.tasks == 0 ) or ( filenumber == len(args.filelist)-1 ) ): 
            
            batchname = getSeparateBatch(localFileList,arguments)
            filecounter += 1
            localFileList = []  
             
            if args.frontend:
                batsub.write('source %s \n' % (batchname))
            else:
                if filecounter%math.ceil(len(args.filelist)/(float(args.jobs)*float(args.tasks)))==0: ### every nth job is the head of a submission chain 
                    batsub.write('msub %s \n' % (batchname))
                else:
                    with open(lastbatchname, 'a') as fout:
                        fout.write( 'msub %s \n' % (batchname) )    
                lastbatchname = batchname
    batsub.close()
 
    if not args.dryrun:  os.system('source ./batchsubmit.sh' )

def getSeparateBatch(filenames,arguments):
    thisnodes,thispart,thishours,thisprocs = args.numnodes,args.partition,args.numhours,args.cores
    filenameStr = args.tag+'_'+os.path.basename(filenames[0]).replace('.xml','')+'-'+os.path.basename(filenames[-1]).replace('.xml','')
    workstrj = '%s_%s_%s_n%d'%(os.path.basename(args.library),os.path.basename(filenameStr),thispart,thisnodes)
    if( args.verbose > 0 ) : print 'writing batch file %s' % workstrj

    jout = open('batch/batch_%s.sh' % (workstrj), 'w') 
    thisname = jout.name
    jout.write( '#!/bin/bash \n')
    jout.write( '##Moab options \n')
    jout.write( '#MSUB -V \n')
    jout.write( '#MSUB -N %s \n'%workstrj )
    jout.write( '#MSUB -l partition=%s \n'%thispart )
    jout.write( '#MSUB -l walltime=%s \n'%thishours )
    jout.write( '#MSUB -l nodes=%d \n'%thisnodes )
    jout.write( '#MSUB -l gres=%s \n'%args.disk )
    if thispart == 'borax': jout.write( '#MSUB -l ttc=%s \n'%thisprocs )
    jout.write( '#MSUB -A %s \n'%args.bank )
    jout.write( '#MSUB -q %s \n'%args.queue )
    jout.write( '#MSUB -j oe \n')
    jout.write( '#MSUB -o %s/logs/log_job_%s.out \n'%(args.library, filenameStr+'_'+timestamp) )
    if not args.frontend : jout.write( 'echo "job_id = $SLURM_JOBID" \n')
    jout.write( 'OMP_NUM_THREADS=%d \n'%thisprocs)
    jout.write( 'export PYTHONPATH=%s/:$PYTHONPATH\n'%args.fudgePath )
    jout.write( 'CWD=$PWD \n' )
    if not args.frontend : jout.write( 'echo $PWD \n' )

    filelist = []
    for filename in filenames:
        filePath = '%s/%s'%(args.library,os.path.dirname(filename))
        fileBase = os.path.basename(filename)
        fileTag = '%s_%s'%(args.tag,fileBase)
        jout.write( 'mkdir -p %s/working/wkdr_%s \n'%(filePath,fileTag) )
        jout.write( 'cp %s/%s  %s/working/wkdr_%s/ \n'%(filePath,fileBase,filePath,fileTag) )   
        filelist.append( '%s/working/wkdr_%s/%s'%(filePath,fileTag,fileBase) )
      
    jout.write( 'multiprocessGeneral.py  %s --code "processProtare.py" --files "%s" --options "%s  > runLog_%s.txt 2>&1"   \n'%(args.procstr, ' '.join(filelist), ' '.join(arguments), args.tag ) ) 

    for filename in filenames:
        filePath = '%s/%s'%(args.library,os.path.dirname(filename))
        fileBase = os.path.basename(filename)
        fileTag = '%s_%s'%(args.tag,fileBase)
        jout.write( 'cp %s/working/wkdr_%s/%s %s/%s/  \n'%(filePath,fileTag,fileBase.replace('.xml','%s.xml'%args.tag),args.library,args.newDir) )   
    
    jout.close()
    return thisname
    
def getQueueLims(partition,nodes,hours,procs=None):
    
    if partition=='borax':
        numhours,numprocs,numnodes=200., 36, 1
    elif partition=='cab':
        numhours,numprocs,numnodes=16., 16, 258
    elif partition=='syrah':
        numhours,numprocs,numnodes=16., 16, 128
    elif partition=='quartz':
        numhours,numprocs,numnodes=16., 36, 256
    elif partition=='rzmerl':
        numhours,numprocs,numnodes=24., 16, 32
    elif partition=='rztopaz':
        numhours,numprocs,numnodes=200., 16, 1
        
    if nodes>numnodes:
        if( args.verbose > 1 ) : print 'warning: too many nodes requested, using %d'%numnodes
        nodes = numnodes
    if hours>numhours:
        if( args.verbose > 1 ) : print 'warning: too many hours requested. using %d'%numhours
        hours = numhours
    
    import math      
    fracPart, intPart = math.modf(hours)
    hourStr = '%d:%02i:00'%(intPart,math.floor(fracPart*60))
    
    if procs: numprocs = procs  ### user overrides number of machine procs
    
    return numprocs,hourStr,nodes


if __name__=='__main__':
    main()

