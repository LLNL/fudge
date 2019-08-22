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

import sys, os, glob, time, multiprocessing
from argparse import ArgumentParser

description = """Runs multiprocess threading on a general script"""
### Example
### python multiprocessGeneral.py --code ./convertEndfToGnd.py --glob "endl2009.2/endf/*/*.endf" 

parser = ArgumentParser( description = description )
parser.add_argument( "-c", "--code",        type = str, default = "ls -d",  help = "code to use" )                                                             
parser.add_argument( "-g", "--glob",        type = str, default = None,     help = "pattern to glob for files. args.files overrides this option." )            
parser.add_argument( "-f", "--files",       type = str, default = None,     help = "list of file arguments to pass to code" )                                  
parser.add_argument( "-o", "--options",     type = str, default = " ",      help = "other options to pass to code" )                                           
parser.add_argument( "-n", "--NumProcesses",type = int, default = 8,        help = "number of simultaneous processes" )      
parser.add_argument( "--no_changeDir", action = 'store_true', default = False, help = "do not change directories before executing the code" )      
parser.add_argument( "-v", "--verbose", action = "count", default = 0,      help = "enable verbose output" )                 

args = parser.parse_args( )

files = None
if args.glob : files = sorted( glob.glob( args.glob ) )
if args.files : files = args.files.split()
if args.verbose > 1 : 
    print( '\ncode : %s' % args.code )
    print( '\noptions : %s' % args.options )
    print( '\nfile list : ',files,'\n' )

def f( filename ) :

    baseName = os.path.basename( filename )
    dirName = os.path.dirname( filename )
    if args.no_changeDir : status = os.system( '%s %s  %s  ' % ( args.code, filename, args.options) )
    else : status = os.system( 'cd %s; %s %s  %s  ' % ( dirName, args.code, baseName, args.options) )
    if( status != 0 ) : status = 1
    sys.exit( status )

k = 0
def addProcess( processes ) :

    global k
    process = multiprocessing.Process( target = f, args = ( files[k], ) )
    process.start( )
    processes.append( [ k, files[k], process ] )
    k += 1

n1, processes = len( files ), []
for i in range( 0, min( n1, args.NumProcesses ) ) : addProcess( processes )

ik, done = 0, {}
while( ik < n1 ) :
    nextProcesses = []
    for j, zas, process in processes :
        if( process.is_alive( ) ) :
            nextProcesses.append( [ j, zas, process ] )
        else :
            s = "%3d of %3d : %s" % ( j + 1, n1, zas )
            if( process.exitcode != 0 ) : s += ' ********* FAILED'
            done[j] = s
            if( k < n1 ) : addProcess( nextProcesses )
    while( ik in done ) :       # Print results ordered.
        if args.verbose > 0: print( done[ik] )
        del done[ik]
        ik += 1
    time.sleep( 0.1 )
    processes = nextProcesses
