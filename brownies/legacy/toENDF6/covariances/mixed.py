# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from fudge.covariances import mixed as mixedModule

from .. import endfFormats as endfFormatsModule

def toENDF6(self, flags, targetInfo):
    endf = []
    XMF1,XLFS1,MAT1 = 0,0,0
    NI = len([cov for cov in self.components if hasattr(cov,'matrix')])
    NC = len(self.components) - NI
    rowdat, coldat = targetInfo['dataPointer']
    MF, MT1 = list( map(int, rowdat.ENDF_MFMT.split(',') ) )
    if coldat:
        MF1, MT1 = list( map(int, coldat.ENDF_MFMT.split(',') ) )
    if MF in (31,33):
        endf.append( endfFormatsModule.endfHeadLine(XMF1,XLFS1,MAT1,MT1,NC,NI) )
    for cov in self.components:
        endf += cov.toENDF6(flags, targetInfo, inCovarianceGroup=True)
    return endf

mixedModule.MixedForm.toENDF6 = toENDF6
