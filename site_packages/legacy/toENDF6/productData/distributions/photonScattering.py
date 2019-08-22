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

import site_packages.legacy.toENDF6.endfFormats as endfFormatsModule
import site_packages.legacy.toENDF6.gndToENDF6 as gndToENDF6Module
import fudge.gnd.productData.distributions.photonScattering as photonScatteringModule

#
# scatteringFactor.pointwise
#
def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

    Z = targetInfo['ZA'] / 1000
    endfInterpolation = gndToENDF6Module.gndToENDFInterpolationFlag( self.interpolation )
    data = []
    for xy in self.copyDataToXYs( xUnitTo = 'eV', yUnitTo = '' ) : data += xy
    endfMFList[27][self.ENDFMT] = [ endfFormatsModule.endfHeadLine( targetInfo['ZA'], targetInfo['mass'],  0, 0, 0, 0 ),
                                    endfFormatsModule.endfContLine( 0, Z,  0, 0, 1, len( data ) / 2 ) ] + \
                                    endfFormatsModule.endfInterpolationList( ( len( data ) / 2, endfInterpolation ) ) + \
                                    endfFormatsModule.endfDataList( data ) + [ endfFormatsModule.endfSENDLineNumber( ) ]

photonScatteringModule.scatteringFactor.pointwise.toENDF6 = toENDF6
#
# scatteringFactor.piecewise
#
def toENDF6( self, MT, endfMFList, flags, targetInfo, energyUnit = 'eV' ) :

    Z = targetInfo['ZA'] / 1000

    endfInterpolation, data = [], []
    counter, lastX, lastY = 0, None, None
    for region in self :
        ENDFInterpolation = gndToENDF6Module.gndToENDFInterpolationFlag( region.interpolation )
        regionData = region.copyDataToXYs( xUnitTo = energyUnit, yUnitTo = '' )
        if( lastX is not None ) :
            if( lastY == regionData[0][1] ) : del regionData[0]
        counter += len( regionData )
        endfInterpolation.append( counter )
        endfInterpolation.append( ENDFInterpolation )
        for xy in regionData : data += xy
        lastX, lastY = regionData[-1]

    endfMFList[27][self.ENDFMT] = [ endfFormatsModule.endfHeadLine( targetInfo['ZA'], targetInfo['mass'],  0, 0, 0, 0 ),
                                    endfFormatsModule.endfContLine( 0, Z,  0, 0, len( endfInterpolation ) / 2, len( data ) / 2 ) ] + \
                                    endfFormatsModule.endfInterpolationList( endfInterpolation ) + \
                                    endfFormatsModule.endfDataList( data ) + [ endfFormatsModule.endfSENDLineNumber( ) ]

photonScatteringModule.scatteringFactor.piecewise.toENDF6 = toENDF6

#
# coherent.form
#
def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

    if( self.formFactor is not None ) : self.formFactor.toENDF6( MT, endfMFList, flags, targetInfo, energyUnit = '1/Ang' )
    if( self.anomalousScatteringFactor_realPart is not None ) : self.anomalousScatteringFactor_realPart.toENDF6( MT, endfMFList, flags, targetInfo )
    if( self.anomalousScatteringFactor_imaginaryPart is not None ) : self.anomalousScatteringFactor_imaginaryPart.toENDF6( MT, endfMFList, flags, targetInfo )

photonScatteringModule.coherent.form.toENDF6 = toENDF6

#
# incoherent.form
#
def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

    self.scatteringFunction.toENDF6( MT, endfMFList, flags, targetInfo, energyUnit = '1/Ang' )

photonScatteringModule.incoherent.form.toENDF6 = toENDF6
