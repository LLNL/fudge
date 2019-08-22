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

from fudge.gnd.productData import multiplicity as multiplicityModule
from fudge.gnd.productData.distributions import angular as angularModule

from .. import gndToENDF6 as gndToENDF6Module
from .. import endfFormats as endfFormatsModule

#
# constant1d
#
def toENDF6List( self, targetInfo ) :

    nPoints = 2
    interpolationFlatData = [ nPoints, 2 ]
    endfMult = endfFormatsModule.endfNdDataList( [ targetInfo['EMin'], self.constant, targetInfo['EMax'], self.constant ] )
    return( interpolationFlatData, nPoints, endfMult )

multiplicityModule.constant1d.toENDF6List = toENDF6List

#
# XYs1d 
#
def toENDF6List( self, targetInfo ) :

    nPoints = len( self )
    interpolationFlatData = [ nPoints, gndToENDF6Module.gndToENDFInterpolationFlag( self.interpolation ) ]
    endfMult = endfFormatsModule.endfNdDataList( self.copyDataToXYs( ) )
    return( interpolationFlatData, nPoints, endfMult )

multiplicityModule.XYs1d.toENDF6List = toENDF6List

#
# regions1d
#
def toENDF6List( self, targetInfo ) :

    interpolationFlatData, multiplicityFlatData = [], []
    counter = 0
    lastX, lastY = None, None
    for region in self :
        ENDFInterpolation = gndToENDF6Module.gndToENDFInterpolationFlag( region.interpolation )
        data = region.copyDataToXYs( )
        if( lastX is not None ) :
            if( lastY == data[0][1] ) : data = data[1:]
        counter += len( data )
        interpolationFlatData.append( counter )
        interpolationFlatData.append( ENDFInterpolation )
        for xy in data : multiplicityFlatData += xy
        lastX, lastY = data[-1]
    return( interpolationFlatData, counter, endfFormatsModule.endfDataList( multiplicityFlatData ) )

multiplicityModule.regions1d.toENDF6List = toENDF6List

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
