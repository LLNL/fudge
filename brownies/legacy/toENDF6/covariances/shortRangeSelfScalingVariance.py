# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from xData import xDataArray as arrayModule

from fudge.covariances import shortRangeSelfScalingVariance as shortRangeSelfScalingVarianceModule

from .. import endfFormats as endfFormatsModule

def toENDF6(self, flags, targetInfo, inCovarianceGroup=False):

    endf = []
    rowdat, coldat = targetInfo['dataPointer']
    MF, MT1 = list( map( int, rowdat.ENDF_MFMT.split(',') ) )
    if not inCovarianceGroup:
        # print header for this subsection (contains one NL sub-subsection)
        MAT1 = targetInfo['MAT1']
        XMF1,XLFS1,NC,NI = 0,0,0,1
        if coldat:
            MF1, MT1 = list( map( int, coldat.ENDF_MFMT.split(',') ) )
        if MF in (31,33):
            endf.append( endfFormatsModule.endfHeadLine(XMF1,XLFS1,MAT1,MT1,NC,NI) )
    # header for matrix:
    if isinstance( self.matrix.array, arrayModule.diagonal ):
        LS = 0; LB = 8; NP = len(self.matrix.axes[2].values); NT = 2*NP
        if self.dependenceOnProcessedGroupWidth==shortRangeSelfScalingVarianceModule.directToken: LB = 9
        matrixData = list( zip( self.matrix.axes[2].values, list(self.matrix.array.values) + [ 0 ] ) )
        matrixData = [val for sublist in matrixData for val in sublist] # flatten
    else:
        raise NotImplementedError("shortRangeSelfScalingVariance with off-diagonal terms")
    endf.append( endfFormatsModule.endfHeadLine( 0,0,LS,LB,NT,NP ) )
    endf += endfFormatsModule.endfDataList( matrixData )
    return endf

shortRangeSelfScalingVarianceModule.shortRangeSelfScalingVariance.toENDF6 = toENDF6
