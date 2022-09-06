# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from fudge.covariances import summed as summedModule

from .. import endfFormats as endfFormatsModule

def toENDF6(self, flags, targetInfo, inCovarianceGroup=False):
    endf = []
    if not inCovarianceGroup:
        # print header for this subsection (contains one NL sub-subsection)
        XMF1,XLFS1,MAT1,NC,NI = 0,0,0,1,0
        rowdat, coldat = targetInfo['dataPointer']
        MT1 = int( rowdat.ENDF_MFMT.split(',')[1] )
        endf.append( endfFormatsModule.endfHeadLine(XMF1,XLFS1,MAT1,MT1,NC,NI) )
    # header:
    LTY=0
    endf.append( endfFormatsModule.endfHeadLine(0,0,0,LTY,0,0) )
    NCI = len(self)
    endf.append( endfFormatsModule.endfHeadLine(self.domainMin, self.domainMax, 0, 0, 2*NCI, NCI) )
    mtList = [ int( a.ENDF_MFMT.split(',')[1] ) for a in self ]
    coefficients = [a.coefficient for a in self]
    endf += endfFormatsModule.endfDataList( [i for j in zip(coefficients,mtList)
        for i in j] )
    return endf

summedModule.SummedCovariance.toENDF6 = toENDF6
