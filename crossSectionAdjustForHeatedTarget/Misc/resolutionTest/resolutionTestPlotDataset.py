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

#
# This routine calls resolutionTestFile.py.
#
# USAGE:
#
#   resolutionTestPlotDataset.py ENDL_file temperature dataset_number_within_ENDL_file
#
# OPTIONS:
#
#   ENDL_file                               ENDL cross section file to heat.
#   temperature                             Temperature to heat to in MeV.
#   dataset_number_within_ENDL_file         some ENDL files have more than 1 cross section (e.g., S = 1), this number is the index 
#                                               (0 based) of the cross section dataset to plot.
#
import sys, os, glob
from fudge.core.math.xData import axes, XYs
from fudge.vis.gnuplot import fudgeMultiPlots

if( len( sys.argv ) != 4 ) : raise Exception( 'need ENDL_file, temperature and dataset_number_within_ENDL_file' )

endlFile = sys.argv[1]
endlFileName = os.path.split( endlFile )[1]
status = os.system( 'resolutionTestFile.py %s %s' % ( endlFile, sys.argv[2] ) )
if( status != 0 ) : raise

nDataset = int( sys.argv[3] )
axes_ = axes.defaultAxes( labelsUnits = { 0 : ( 'energy_in', 'eV' ), 1 : ( 'crossSection' , 'b' ) } )

def getDataset( n, file ) :

    f = open( file )
    ls = f.readlines( )
    f.close( )
    for i in xrange( n + 1 ) :
        if( len( ls ) < 5 ) : raise Exception( 'file %s only has %d datasets' % ( file, i ) )
        ls = ls[2:]
        j = 0
        for l in ls :
            if( l == '                                                                       1\n' ) : break
            j += 1
        if( i == n ) :
            data = []
            for k in xrange( j ) : data.append( map( float, ls[k].split( ) ) )
            return( data )
        else :
            ls = ls[j+1:]

endlData = XYs.XYs( axes_, getDataset( nDataset, endlFile ), 1e-3 )
endlData.label = endlFileName
lowData = XYs.XYs( axes_, getDataset( nDataset, glob.glob( '%s_T*eV_low' % endlFileName )[0] ), 1e-3 )
lowData.label = 'low res.'
highData = XYs.XYs( axes_, getDataset( nDataset, glob.glob( '%s_T*eV_high' % endlFileName )[0] ), 1e-3 )
highData.label = 'high res.'

fudgeMultiPlots.multiPlot( [ endlData, lowData, highData ] )
