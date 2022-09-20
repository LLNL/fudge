# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from fudge.productData import multiplicity as multiplicityModule

from .. import gndsToENDF6 as gndsToENDF6Module
from .. import endfFormats as endfFormatsModule

#
# Constant1d
#
def toENDF6List( self, targetInfo ) :

    nPoints = 2
    interpolationFlatData = [ nPoints, 2 ]
    endfMult = endfFormatsModule.endfNdDataList( [ targetInfo['EMin'], self.value, targetInfo['EMax'], self.value ] )
    return( interpolationFlatData, nPoints, endfMult )

multiplicityModule.Constant1d.toENDF6List = toENDF6List

#
# XYs1d 
#
def toENDF6List( self, targetInfo ) :

    nPoints = len( self )
    interpolationFlatData = [ nPoints, gndsToENDF6Module.gndsToENDFInterpolationFlag( self.interpolation ) ]
    endfMult = endfFormatsModule.endfNdDataList( self.copyDataToXYs( ) )
    return( interpolationFlatData, nPoints, endfMult )

multiplicityModule.XYs1d.toENDF6List = toENDF6List

#
# Regions1d
#
def toENDF6List( self, targetInfo ) :

    interpolationFlatData, multiplicityFlatData = [], []
    counter = 0
    lastX, lastY = None, None
    for region in self :
        ENDFInterpolation = gndsToENDF6Module.gndsToENDFInterpolationFlag( region.interpolation )
        data = region.copyDataToXYs( )
        if( lastX is not None ) :
            if( lastY == data[0][1] ) : data = data[1:]
        counter += len( data )
        interpolationFlatData.append( counter )
        interpolationFlatData.append( ENDFInterpolation )
        for xy in data : multiplicityFlatData += xy
        lastX, lastY = data[-1]
    return( interpolationFlatData, counter, endfFormatsModule.endfDataList( multiplicityFlatData ) )

multiplicityModule.Regions1d.toENDF6List = toENDF6List

#
# Polynomial1d
#
def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

    endfMFList[1][MT]  = [ endfFormatsModule.endfHeadLine( targetInfo['ZA'], targetInfo['mass'], 0, 1, 0, 0 ) ]
    endfMFList[1][MT] += [ endfFormatsModule.endfHeadLine( 0, 0, 0, 0, len( self ), 0 ) ]
    endfMFList[1][MT] += endfFormatsModule.endfDataList( self.coefficients ) + [ endfFormatsModule.endfSENDLineNumber( ) ]

multiplicityModule.Polynomial1d.toENDF6 = toENDF6

def fissionNeutronsToENDF6( MT, multiplicity, endfMFList, flags, targetInfo ) :

    if( isinstance( multiplicity, multiplicityModule.Polynomial1d ) ) :
        multiplicity.toENDF6( MT, endfMFList, flags, targetInfo )
    else :
        interpolation, n1, data = multiplicity.toENDF6List( targetInfo )
        if( MT not in endfMFList[1] ) :
            endfMFList[1][MT] = [ endfFormatsModule.endfHeadLine( targetInfo['ZA'], targetInfo['mass'], 0, 2, 0, 0 ) ]
        endfMFList[1][MT] += [ endfFormatsModule.endfHeadLine( 0, 0, 0, 0, len( interpolation ) / 2, n1 ) ]
        endfMFList[1][MT] += endfFormatsModule.endfInterpolationList( interpolation )
        endfMFList[1][MT] += data + [ endfFormatsModule.endfSENDLineNumber( ) ]
