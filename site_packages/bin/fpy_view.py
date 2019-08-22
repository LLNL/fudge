#! /usr/bin/env python

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

import sys, os, argparse
binDir = os.path.dirname( os.path.abspath( __file__ ) )
sys.path.insert( 0, os.path.dirname( binDir ) )

def process_args():
    parser = argparse.ArgumentParser(description = 'Translate an ENDF file to the new GNDS format')
    parser.set_defaults( verbose = True )
    parser.add_argument("-v", action="store_true", dest='verbose', help="enable verbose output")
    parser.add_argument("-q", action="store_false", dest='verbose', help="disable verbose output")
    parser.add_argument("-ifpy", action="store_true", dest='ifpy', help="show independent (prompt) fission yields")
    parser.add_argument("-cfpy", action="store_true", dest='cfpy', help="show cummulative (delayed) fission yields")
    parser.add_argument("--za", default=None, type=int, dest='za', help='the za to attempt to learn about (optional)' )
    parser.add_argument("--state", default=0, type=int, dest='state', help='the excitation of the za to look for (default=0 for g.s.)' )
    parser.add_argument('--fudgepath', default='', dest='fudgepath', type=str, help="set the path to fudge (use if not already in your PYTHONPATH)")
    parser.add_argument("inFile", type=str, help="the input endf file (should have '.endf' extension)")
    return parser.parse_args()

def grok_yields( d ):
    isoList = d.keys()
    isoList.sort()
    print isoList 

def print_sum_by_A( d, iE=0 ): 
    result = 200*[0]
    dresult = 200*[0]
    for key in d:
        A = int( str( key[0] )[-3:] )
        result[A] += d[key][iE][1]
        dresult[A] += d[key][iE][2]
    #return result, dresult
    for i in range( 200 )[75:175]:
        print str(i).rjust(4), int( result[i]*200 )*'*'

def print_ZA_by_E( d, ZA, state=0 ): 
    for E, Y, DY in d[ ( ZA, state ) ]: print "{:<11} {:>11} +/- {:<5}%".format( E/1e6, Y*100.0, (DY/Y)*100.0 )

if __name__ == "__main__":
    args = process_args()
    if args.fudgepath != '' and args.fudgepath not in sys.path: sys.path.append( args.fudgepath )
    from fudge.legacy.converting.endfFileToGNDSMisc import parseENDFByMT_MF
    from fudge.legacy.converting.endfFileToGNDS import endfFileToGNDS
    header, MAT, theData = parseENDFByMT_MF( args.inFile )

    x, c, i = endfFileToGNDS( args.inFile )
    
    # ifpy = ( 454, 8 ), independent (prompt) fission yields
    # cfpy = ( 459, 8 ), cummulative (delayed) fission yields
    
    if args.ifpy: theData = i.independentFissionYields
    if args.cfpy: theData = i.cumulativeFissionYields

    if args.za != None: print_ZA_by_E( theData, args.za, state=args.state )
    else: print_sum_by_A( theData, iE=0 )
#    else: grok_yields( theData )
