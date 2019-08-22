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

from fudge import *
sys.path.append( './' )
readOnly( )
import crossSectionAdjustForHeatedTarget
import time

endl_bdfls = bdfls.getDefaultBdfls( )

e = endlProject( )
projectileMass = endl_bdfls.mass( e.yi )

T = 1e-1

z = e.readZA( 94239 )
targetMass = endl_bdfls.mass( z.ZA )
massRatio = targetMass / projectileMass
z.read( )
a = z.findData( C = 15, I = 0 )
# f = open( "o", "w" )
# f.write( a.toString( ) )
# f.close( )
fi = 5e-3
Ts = [ 1e-8, 1e-6, 1e-4, 1e-2 ]
for T in Ts :
	print "Do T =", T
	t1 = time.clock( )
	hj = endl2dmath( crossSectionAdjustForHeatedTarget.crossSectionAdjustForHeatedTarget( massRatio, T, 1e-10, a.data, heatAllPoints = 0, interpolationAccuracy = fi))
	t2 = time.clock( )
	ha = endl2dmath( crossSectionAdjustForHeatedTarget.crossSectionAdjustForHeatedTarget( massRatio, T, 1e-10, a.data, heatAllPoints = 1, interpolationAccuracy = fi))
	t3 = time.clock( )
	print t1, t2, t3, t2 - t1, t3 - t2

#	q = qmultiPlot( ( a, hj, ha ), legends = ( "unheated", "judicial", "all" ), xylog = 3 )
#	q2 = qmultiPlot( ( a, hj, ha ), legends = ( "unheated", "judicial", "all" ), xylog = 3, xMin = .01 )

#	s = "T=%.1e" % T
#	f = open( "a" + s, "w" )
#	f.write( hj.toString( ) )
#	f.close( )

#	f = open( "b" + s, "w" )
#	f.write( ha.toString( ) )
#	f.close( )

	d = hj - ha
#	d.plot( xylog = 1 )
	d = d / hj
	d.plot( xylog = 1 )
	print len( a ), len( hj ), len( ha )

	Pause( )
