# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from brownies.legacy.endl import *
sys.path.append( './' )
readOnly( )
import crossSectionAdjustForHeatedTarget
import time

if hasattr(crossSectionAdjustForHeatedTarget, 'heat'):
	from crossSectionAdjustForHeatedTarget import heat as crossSectionAdjustForHeatedTarget

# make new instance from reading a file
endl_bdfls = bdfls.getDefaultBdfls( )

# initialize new endl project (probably from a file)
e = endlProject( )
projectileMass = endl_bdfls.mass( e.yi )

# ? Tmperature in MeV/k
T = 1e-1

z = e.readZA( 94239 ) # point to ENDL data folder
targetMass = endl_bdfls.mass( z.ZA )
massRatio = targetMass / projectileMass
z.read( ) # All data in folder
a = z.findData( C = 15, I = 0 )
# f = open( "o", "w" )
# f.write( a.toString( ) )
# f.close( )
fi = 5e-3
Ts = [ 1e-8, 1e-6, 1e-4, 1e-2 ]
for T in Ts :
	print( "Do T =", T )
	t1 = time.clock( )
	hj = endl2dmath( crossSectionAdjustForHeatedTarget.crossSectionAdjustForHeatedTarget( massRatio, T, 1e-10, a.data, heatAllPoints = 0, interpolationAccuracy = fi ) )
	t2 = time.clock( )
	ha = endl2dmath( crossSectionAdjustForHeatedTarget.crossSectionAdjustForHeatedTarget( massRatio, T, 1e-10, a.data, heatAllPoints = 1, interpolationAccuracy = fi ) )
	t3 = time.clock( )
	print( t1, t2, t3, t2 - t1, t3 - t2 )

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
	print( len( a ), len( hj ), len( ha ) )

	Pause( )
