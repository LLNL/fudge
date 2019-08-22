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
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
