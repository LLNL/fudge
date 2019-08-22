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

from fudge.core.utilities import brb

import site_packages.legacy.toENDF6.gndToENDF6 as gndToENDF6Module
import site_packages.legacy.toENDF6.endfFormats as endfFormatsModule
import fudge.gnd.productData.multiplicity as multiplicityModule
import fudge.gnd.productData.distributions.angular as angularModule

#
# constant
#
def toENDF6List( self, targetInfo ) :

    nPoints = 2
    interpolationFlatData = [ nPoints, 2 ]
    endfMult = endfFormatsModule.endfNdDataList( [ targetInfo['EMin'], self.value, targetInfo['EMax'], self.value ] )
    return( interpolationFlatData, nPoints, endfMult )

multiplicityModule.constant.toENDF6List = toENDF6List

#
# pointwise
#
def toENDF6List( self, targetInfo ) :

    nPoints = len( self )
    interpolationFlatData = [ nPoints, gndToENDF6Module.gndToENDFInterpolationFlag( self.interpolation ) ]
    endfMult = endfFormatsModule.endfNdDataList( self.copyDataToXYs( xUnitTo = 'eV' ) )
    return( interpolationFlatData, nPoints, endfMult )

multiplicityModule.pointwise.toENDF6List = toENDF6List

#
# piecewise
#
def toENDF6List( self, targetInfo ) :

    interpolationFlatData, multiplicityFlatData = [], []
    counter = 0
    lastX, lastY = None, None
    for region in self :
        ENDFInterpolation = gndToENDF6Module.gndToENDFInterpolationFlag( region.interpolation )
        data = region.copyDataToXYs( xUnitTo = 'eV' )
        if( lastX is not None ) :
            if( lastY == data[0][1] ) : data = data[1:]
        counter += len( data )
        interpolationFlatData.append( counter )
        interpolationFlatData.append( ENDFInterpolation )
        for xy in data : multiplicityFlatData += xy
        lastX, lastY = data[-1]
    return( interpolationFlatData, counter, endfFormatsModule.endfDataList( multiplicityFlatData ) )

multiplicityModule.piecewise.toENDF6List = toENDF6List

#
# polynomial
#
def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

    endfMFList[1][MT]  = [ endfFormatsModule.endfHeadLine( targetInfo['ZA'], targetInfo['mass'], 0, 1, 0, 0 ) ]
    endfMFList[1][MT] += [ endfFormatsModule.endfHeadLine( 0, 0, 0, 0, len( self ), 0 ) ]
    endfMFList[1][MT] += endfFormatsModule.endfDataList( self.coefficients ) + [ endfFormatsModule.endfSENDLineNumber( ) ]

multiplicityModule.polynomial.toENDF6 = toENDF6

def fissionNeutronsToENDF6( MT, multiplicity, endfMFList, flags, targetInfo ) :

    if( isinstance( multiplicity, multiplicityModule.polynomial ) ) :
        multiplicity.toENDF6( MT, endfMFList, flags, targetInfo )
    else :
        interpolation, n1, data = multiplicity.toENDF6List( targetInfo )
        if( MT not in endfMFList[1] ) :
            endfMFList[1][MT] = [ endfFormatsModule.endfHeadLine( targetInfo['ZA'], targetInfo['mass'], 0, 2, 0, 0 ) ]
        endfMFList[1][MT] += [ endfFormatsModule.endfHeadLine( 0, 0, 0, 0, len( interpolation ) / 2, n1 ) ]
        endfMFList[1][MT] += endfFormatsModule.endfInterpolationList( interpolation )
        endfMFList[1][MT] += data + [ endfFormatsModule.endfSENDLineNumber( ) ]
