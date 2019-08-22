#! /usr/bin/env python

# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
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
    d = d.xSlice( xMin = xMin, xMax = xMax )
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
