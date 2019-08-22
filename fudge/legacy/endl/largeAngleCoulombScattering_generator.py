#! /bin/env python

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

import sys, math

yi = int( sys.argv[1] )
ZA = int( sys.argv[2] )
workDir = sys.argv[3]

from fudge.legacy.endl import bdfls
endl_bdfls =  bdfls.getDefaultBdfls ( template = '/usr/gapps/data/nuclear/bdfls.archive/bdfls.Audi_etal.2003.12.22')
from fudge.legacy import endl

def xSec( Es, z1, m1, z2, m2 ) :

    if( z2 == 0 ) : z2 = 1e-10
    m = ( m1 + m2 ) / m2
    C = 2.07323e-2 * math.pi / 2 * z1 * z1 * z2 * z2 * m * m * 16.16666667
    xsecs = []
    for E in Es : xsecs.append( [ E, C / E / E ] )
    return( xsecs )

number, EMin, EMax = 200, 1e-4, 30
frac = math.pow( EMax / EMin, 1. / number )
energies = [ EMin * math.pow( frac, i1 ) for i1 in range( number ) ]
energies.append( EMax )

angularStr = """-1.00000e+00  1.53490e-02
     -8.21610e-01  1.85030e-02
     -6.65520e-01  2.21340e-02
     -5.20570e-01  2.65540e-02
     -3.97930e-01  3.14180e-02
     -2.86440e-01  3.71000e-02
     -1.86090e-01  4.36430e-02
     -1.88510e-02  5.91470e-02
      1.26090e-01  8.03940e-02
      2.37590e-01  1.05630e-01
      3.37930e-01  1.40070e-01
      4.15980e-01  1.80010e-01
      4.94020e-01  2.39820e-01
      5.49770e-01  3.02890e-01
      5.94370e-01  3.73150e-01
      6.38970e-01  4.71030e-01
      6.72410e-01  5.72140e-01
      7.05860e-01  7.09630e-01
      7.39310e-01  9.03460e-01
      7.61610e-01  1.08040e+00
      7.83910e-01  1.31480e+00
      8.06210e-01  1.63480e+00
      8.28510e-01  2.08760e+00
      8.39660e-01  2.38810e+00
      8.61950e-01  3.22180e+00
      8.73100e-01  3.81280e+00
      8.84250e-01  4.58280e+00
      8.95400e-01  5.61180e+00
      9.06550e-01  7.03060e+00
      9.17700e-01  9.06470e+00
      9.28850e-01  1.21280e+01
      9.39900e-01  1.69980e+01
      9.40000e-01  0.00000e+00
      1.00000e+00  0.00000e+00"""

project = endl.endlProject( projectile = yi, workDir = workDir )
target = project.addZA( ZA )

halflife = endl_bdfls.halflife( ZA )
if( halflife is None ) :            # Happens when bdfls file does not have halflife for ZA.
    ZA_halflife = target.ZA
    if( target.ZA % 1000 == 0 ) : ZA_halflife = 0
    halflife = { 0 : 1e50, 99120 : 1e50, 99125 : 1e50, 24047 : 0.5, 28067 : 27., 30073 : 23.5 }[ZA_halflife]

xSecData = xSec( energies, project.Z, project.mass, target.Z, endl_bdfls.mass( target.ZA ) )
fileI0 = target.addFile( 0, 8, 0, 0 )
fileI0.addData( xSecData, halflife = halflife )

angularData = angularStr.split( '\n' )
angularData = [ map( float, line.split( ) ) for line in angularData ]
angularData = [ [ EMin, angularData ], [ EMax, angularData ] ]
fileI1 = target.addFile( yi, 8, 1, 0 )
fileI1.addData( angularData, halflife = halflife )

project.save( )
