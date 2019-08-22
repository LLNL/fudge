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
    
import site_packages.legacy.toENDF6.gndToENDF6 as gndToENDF6Module
import site_packages.legacy.toENDF6.endfFormats as endfFormatsModule
import fudge.gnd.reactionData.crossSection as crossSectionModule

def toENDF6( self, MT, endfMFList, targetInfo, level, LR ) :
    """
    Convert self into ENDF format

    :param int MT: The ENDF reaction designator, MT
    :param endfMFList:
    :param targetInfo:
    :param level:
    :param LR:
    """

    ZA, mass, QI, QM = targetInfo['ZA'], targetInfo['mass'], targetInfo['Q'], targetInfo['QM']
    if( 'EFL' in targetInfo ) :
        QM = QI
        QI = targetInfo['EFL']
    else :
        if( QM is None ) :
            if( MT in ( 2, 5 ) ) :
                QM = QI
            elif( MT == 4 ) :                               # Q should be 0 except for excited-state targets:
                QM = 0
                if( hasattr( targetInfo['reactionSuite'].target, 'getLevelIndex' ) ) :
                    if( targetInfo['reactionSuite'].target.getLevelIndex() > 0 ) : QM = QI
            else :
                QM = QI + level
    interpolationFlatData, crossSectionFlatData = self[targetInfo['style']].toENDF6Data( MT, endfMFList, targetInfo, level )
    MF = targetInfo['crossSectionMF']
    endfMFList[MF][MT] = [ endfFormatsModule.endfHeadLine( ZA, mass, 0, 0, 0, 0 ) ]
    endfMFList[MF][MT].append( endfFormatsModule.endfContLine( QM, QI, 0, LR, len( interpolationFlatData ) / 2, len( crossSectionFlatData ) / 2 ) )
    endfMFList[MF][MT] += endfFormatsModule.endfInterpolationList( interpolationFlatData )
    endfMFList[MF][MT] += endfFormatsModule.endfDataList( crossSectionFlatData )
    endfMFList[MF][MT].append( endfFormatsModule.endfSENDLineNumber( ) )

crossSectionModule.component.toENDF6 = toENDF6

def toENDF6Data( self, MT, endfMFList, targetInfo, level ) :

    endfInterpolation = gndToENDF6Module.gndToENDFInterpolationFlag( self.interpolation )
    crossSectionFlatData = []
    for xy in self.copyDataToXYs( xUnitTo = 'eV', yUnitTo = 'b' ) : crossSectionFlatData += xy
    return( [ len( crossSectionFlatData ) / 2, endfInterpolation ], crossSectionFlatData )

crossSectionModule.pointwise.toENDF6Data = toENDF6Data

def toENDF6Data( self, MT, endfMFList, targetInfo, level ) :

    interpolationFlatData, crossSectionFlatData = [], []
    counter = 0
    lastX, lastY = None, None
    for region in self :
        ENDFInterpolation = gndToENDF6Module.gndToENDFInterpolationFlag( region.interpolation )
        data = region.copyDataToXYs( xUnitTo = 'eV', yUnitTo = 'b' )
        if( lastX is not None ) :
            if( lastY == data[0][1] ) :
                data = data[1:]
            elif( ( lastY == 0 ) and region.interpolation[4:] == 'log' ) :
                interpolationFlatData[-2] += 1
        counter += len( data )
        interpolationFlatData.append( counter )
        interpolationFlatData.append( ENDFInterpolation )
        for xy in data : crossSectionFlatData += xy
        lastX, lastY = data[-1]
    return( interpolationFlatData, crossSectionFlatData )

crossSectionModule.piecewise.toENDF6Data = toENDF6Data

def toENDF6Data( self, MT, endfMFList, targetInfo, level ) :

    return self.tabulatedData.toENDF6Data( MT, endfMFList, targetInfo, level )

crossSectionModule.resonancesWithBackground.toENDF6Data = toENDF6Data

def toENDF6Data( self, MT, endfMFList, targetInfo, level ) :

    endfInterpolation = gndToENDF6Module.axesToEndfInterpolationFlag( self.weights.axes )
    crossSectionFlatData = []
    for xy in self.weights.copyDataToXYs() : crossSectionFlatData += xy
    return( [ len( crossSectionFlatData ) / 2, endfInterpolation ], crossSectionFlatData )

crossSectionModule.weightedPointwise.toENDF6Data = toENDF6Data
