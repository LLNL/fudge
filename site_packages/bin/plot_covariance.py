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

try:                import argparse
except ImportError: from fudge.core.utilities import argparse

from site_packages.BNL.plot_evaluation import plotio as io

# Process command line options
parser = argparse.ArgumentParser(description='Plot the covariance data in an ENDF file')
parser.add_argument('inFile', type=str, help='The ENDF file you want to check.' )
parser.add_argument('MT', type=int, help='The MT you want to plot.' )
parser.add_argument('--crossMT', type=int, default=None, help='The MT of a cross correlation you want to plot' )
parser.add_argument('--MF', type=int, default=33, help='MF to plot (Default is 33 for cross sections)' )
parser.add_argument('--abs', dest='abs', default=None, action='store_true', help='Change covariance to absolute' )
parser.add_argument('--rel', dest='rel', default=None, action='store_true', help='Change covariance to relative' )
parser.add_argument('--corr', dest='corr', default=None, action='store_true', help='Change covariance to a correlation matrix' )
parser.add_argument("--skipBadData", action="store_true", default=False, help="skip bad data, rather than throw an exception, when reading an ENDF file")
parser.add_argument('-v', dest='verbose', default=False, action='store_true', help='Enable verbose output' )
parser.add_argument('-l', dest='list', default=False, action='store_true', help='List covariance data in file' )
args = parser.parse_args()


#if True:
try:
    results = io.readEvaluation( args.inFile, verbose = args.verbose, skipBadData = args.skipBadData )
    evaluation = results[0]
    covariances = results[1]
except Exception as err:
    print( 'WARNING: ENDF READ HALTED BECAUSE '+str(err) )
    exit()
    
if args.list:
    print "\n\nList of covariance data"
    print "======================="
    for c in covariances.sections:
        row = c.rowData.attributes['ENDF_MFMT']
        if c.columnData != None:    col = c.columnData.attributes['ENDF_MFMT']
        else:                       col = row
        print c, row, col
    exit()

#evaluation.reconstructResonances( styleName='reconstructed', accuracy=0.001 )
made_a_plot=False
for c in covariances.sections: 
    if hasattr(c,'rowData') and c.rowData.attributes['ENDF_MFMT'] == '%i,%i' % (args.MF, args.MT): 
        if args.crossMT==None or \
                ( c.columnData != None and c.columnData.attributes['ENDF_MFMT'] == '%i,%i' % (args.MF, args.crossMT) ) or \
                ( c.columnData == None and  (args.MF, args.MT) == (args.MF, args.crossMT) ):
            otherMT = None
            if args.crossMT != None: otherMT=args.crossMT
            elif c.columnData == None: otherMT=args.MT
            else: otherMT = int( c.columnData.attributes['ENDF_MFMT'].split(',')[-1] )
            c2 = c.toCovarianceMatrix()
            c2.setAncestor( c )
            if args.abs:
                c2.toAbsolute().plot(title='(%i,%i) x (%i,%i)'  % (args.MF, args.MT, args.MF, otherMT))
            elif args.rel:
                c2.toRelative().plot(title='(%i,%i) x (%i,%i)'  % (args.MF, args.MT, args.MF, otherMT))
            elif args.corr:
                c2.toAbsolute().toCorrelationMatrix().plot(title='(%i,%i) x (%i,%i)' % (args.MF, args.MT, args.MF, otherMT))
            else:
                c.plot(title='(%i,%i) x (%i,%i)'   % (args.MF, args.MT, args.MF, otherMT))
            made_a_plot=True

if not made_a_plot: 
    if args.crossMT != None: otherMT=args.crossMT
    elif c.columnData == None: otherMT=args.MT
    raise KeyError( "(MF0,MT0)x(MF1,MT1) = (%s,%s)x(%s,%s) not found" % (args.MF, args.MT, args.MF, otherMT) )
