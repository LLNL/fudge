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
# This routine heats the data in an ENDL file to inputted temperature for a resolution of 2e-3 and 1e-4,
# then diffs the results to see if the first is accurate to 2e-3.
#
# USAGE:
#
# resolutionTestFile.py file temperature
#
# OPTIONS:
#
#   file          ENDL cross section file to heat.
#   temperature   temperatre to heat file in units of file's energy data.
#
import sys, os, glob
from fudge.core.math.xData import axes, XYs

if( len( sys.argv ) != 3 ) : raise Exception( 'need file and temperature' )
file = sys.argv[1]
T = sys.argv[2]

if( os.path.exists( 'nuc_xsec_adjust_for_heated_target_endl' ) ) : 
    heater = 'nuc_xsec_adjust_for_heated_target_endl'
elif( os.path.exists( '../nuc_xsec_adjust_for_heated_target_endl' ) ) : 
    heater = '../nuc_xsec_adjust_for_heated_target_endl'
else :
    raise Exception( 'cannot file executable nuc_xsec_adjust_for_heated_target_endl' )

os.system( 'rm -f %s_T*eV' % file )
eps1 = 2e-3
eps2 = 1e-4
axes_ = axes.defaultAxes( labelsUnits = { 0 : ( 'energy_in', 'eV' ), 1 : ( 'crossSection' , 'b' ) } )

def heat( file, eps, T, heater, suffix ) :

    infoName = 'resolutionTestFile.info%s' % suffix
    os.system( '%s %s %s %e > %s 2>&1' % ( heater, T, file, eps, infoName ) )
    f = open( infoName )
    ls = f.readlines( )
    f.close( )
    n = 0
    for l in ls :
        if( 'nuc_xsec_adjust_for_heated_target:' in l ) : n += 1
    if( n > 0 ) : print '   %d resolution or possible infinite loop issues for %s on %s' % ( n, suffix, file )
    fn = glob.glob( '%s_T*eV' % os.path.split( file )[1] )[0]
    f = open( fn )
    ls = f.readlines( )
    f.close( )
    datas = []
    n = len( ls )
    i = 2
    while( i < n ) :
        data = []
        while( ls[i] != '                                                                       1\n' ) :
            data.append( map( float, ls[i].split( ) ) )
            i += 1
        datas.append( XYs.XYs( axes_, data, eps, safeDivide = True ) )
        i += 3
    os.system( 'mv %s %s%s' % ( fn, fn, suffix ) )
    return( datas )

datas1 = heat( file, eps1, T, heater, '_low' )
datas2 = heat( file, eps2, T, heater, '_high' )
epsT = eps1 + eps2
errCount = 0
for i in xrange( len( datas1 ) ) :
    data1 = datas1[i]
    data2 = datas2[i]
    d = data1 - data2
    d.trim( )
    xMin = max( data1.xMin( ), data2.xMin( ) )
    xMax = min( data1.xMax( ), data2.xMax( ) )
    from fudge.core.utilities import brb
    d = d.domainSlice( xMin = xMin, xMax = xMax )
    r = d / data2
    yMax1 = data1.yMax( ) * 1e-10
    yMax2 = yMax1 * 1e-2
    yMax3 = yMax2 * 1e-3
    yMax4 = yMax3 * 1e-3
    yErrMax = 0.
    iii = 0
    for x, y in r :
        eps = epsT
        y1 = data1.getValue( x )
        if( y1 < yMax1 ) :
            eps = epsT * 2.
            if( y1 < yMax4 ) : continue
            if( y1 < yMax2 ) : eps = epsT * 10.
            if( y1 < yMax3 ) : eps = epsT * 100.
        if( abs( y ) > eps ) :
            if( abs( y ) > abs( yErrMax ) ) : xsecErr, xErr, yErrMax = y1, x, y
            if( errCount < 5 ) : print '   r = %e for xsec = %e and %e at E = %e for dataset %d' % ( y, y1, data2.getValue( x ), x, i + 1 )
            errCount += 1
        iii += 1
if( errCount > 0 ) :
    print '   errCount =', errCount, 'for file', file
    print '       xsecErr, xErr, yErrMax =', xsecErr, xErr, yErrMax
