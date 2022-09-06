# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import numpy

from xData import xDataArray as arrayModule

from fudge.covariances import enums as covarianceEnumsModule
from fudge.covariances import covarianceMatrix as covarianceMatrixModule
from fudge.core.math import linearAlgebra as linearAlgebraModule

from .. import endfFormats as endfFormatsModule

def toENDF6(self, flags, targetInfo, inCovarianceGroup=False):

    endf = []
    conversionFlags = targetInfo['ENDFconversionFlags'].get(self,"")
    rowdat, coldat = targetInfo['dataPointer']
    MF,MT1 = list( map( int, rowdat.ENDF_MFMT.split(',') ) )
    if not inCovarianceGroup:
        # print header for this subsection (contains one NL sub-subsection)
        MAT1 = targetInfo['MAT1']
        XMF1,XLFS1,NC,NI = 0,0,0,1
        if coldat:
            MF1, MT1 = list( map( int, coldat.ENDF_MFMT.split(',') ) )
        if MF in (31,33):
            endf.append( endfFormatsModule.endfHeadLine(XMF1,XLFS1,MAT1,MT1,NC,NI) )
    # header for matrix:
    rows,cols = self.matrix.array.shape
    if isinstance( self.matrix.array, arrayModule.Diagonal ):
        LS = 0; LB = 1; NP = len(self.matrix.axes[2].values); NT = 2*NP
        if self.type == covarianceEnumsModule.Type.absolute:
            LB = 0
        if 'LB' in conversionFlags:
            LB = int( conversionFlags.split('=')[1] )
        matrixData = list( zip( self.matrix.axes[2].values, list(self.matrix.array.values) + [0] ) )
        matrixData = [val for sublist in matrixData for val in sublist] # flatten
    elif isinstance( self.matrix.array, arrayModule.Full ):
        LB = 5
        if 'LB' in conversionFlags:
            LB = int( conversionFlags.split('=')[1].split(',')[0] )
        if LB == 2:
            LS = 0; NP = len(self.matrix.axes[2].values); NT = 2*NP
            fullMatrix = self.matrix.array.constructArray()
            vals = numpy.sqrt( numpy.diagonal( fullMatrix ) )
            firstNonZero = numpy.nonzero(vals)[0][0]
            vals = numpy.copysign(vals, fullMatrix[firstNonZero])
            if 'firstNegative' in conversionFlags:
                negativeIndex = int(conversionFlags.split('firstNegativeIndex=')[-1])
                if vals[negativeIndex] > 0:
                    vals *= -1
            matrixData = list( zip( self.matrix.axes[2].values, list(vals) + [0] ) )
            matrixData = [val for sublist in matrixData for val in sublist] # flatten
        elif self.matrix.array.symmetry in (arrayModule.Symmetry.lower, arrayModule.Symmetry.upper):
            LS = 1; NT = (rows+1) + rows*(rows+1)/2; NP = rows+1
            arrayData = list( self.matrix.array.values )
            if self.matrix.array.symmetry == arrayModule.Symmetry.lower:
                arrayData = linearAlgebraModule.switchSymmetry( arrayData, upperToLower=False )
            matrixData = list(self.matrix.axes[2].values) + arrayData
        elif self.matrix.axes[1].isLink():
            LS = 0; NT = (rows+1) + rows*cols; NP = rows+1
            matrixData = list(self.matrix.axes[2].values) + list(self.matrix.array.values)
        else:
            LS = 0; LB = 6; NT = (rows+1) + (cols+1) + rows*cols; NP = rows+1
            matrixData = list(self.matrix.axes[2].values) + list(self.matrix.axes[1].values) + list(
                    self.matrix.array.values)
    else:
        raise NotImplemented
    if MF==35:  # header for fission spectra is different:
        slice_ = rowdat.slices[0]
        E1,E2 = (slice_.domainMin, slice_.domainMax)
        if LS: LB = 7
        else:
            raise Exception ("Unknown spectrum (MF35) covariance format")
        endf.append( endfFormatsModule.endfHeadLine( E1,E2,LS,LB,NT,NP ) )
    else:
        endf.append( endfFormatsModule.endfHeadLine( 0,0,LS,LB,NT,NP ) )
    endf += endfFormatsModule.endfDataList( matrixData )
    return endf

covarianceMatrixModule.CovarianceMatrix.toENDF6 = toENDF6
