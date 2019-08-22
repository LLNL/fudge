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

import sys
import os
import multiprocessing
import time

def executeOne( command ) :
    """
    This function is mainly for an example of what the target argument of multiprocessing.Process may look.
    The main thing is to call sys.exit with a status and to properly handle the status returned by os.system.
    """

    status = os.system( command )
    if( status != 0 ) : status = 1
    sys.exit( status )

def addProcess( numberLaunched, processes, items, callback ) :
    """Added a process."""

    process = multiprocessing.Process( target = callback, args = ( items[numberLaunched], ) )
    process.start( )
    processes.append( [ numberLaunched, items[numberLaunched], process ] )
    return( numberLaunched + 1 )

def execute( numberOfParallelProcesses, items, callback ) :
    """
    This function loops over each item in *items* and processes it by calling user's callback function. At most *numberOfParallelProcesses* of
    processes can run in parallel. The user's callback must create a command to execute and call function **executeOne**.
    """

    numberLaunched = 0
    numberToDo = len( items )
    processes = []
    for i1 in xrange( 0, min( numberToDo, numberOfParallelProcesses ) ) : numberLaunched = addProcess( numberLaunched, processes, items, callback )

    numberDone = 0
    done = {}
    while( numberDone < numberToDo ) :
        nextProcesses = []
        for i1, item, process in processes :
            if( process.is_alive( ) ) :
                nextProcesses.append( [ i1, item, process ] )
            else :
                message = "%3d of %3d : %s" % ( i1 + 1, numberToDo, item )
                if( process.exitcode != 0 ) : message += ' ********* FAILED'
                done[i1] = message
                if( numberLaunched < numberToDo ) : numberLaunched = addProcess( numberLaunched, nextProcesses, items, callback )
        while( numberDone in done ) :       # Print results ordered.
            print done[numberDone]
            del done[numberDone]
            numberDone += 1
        time.sleep( 0.1 )
        processes = nextProcesses
