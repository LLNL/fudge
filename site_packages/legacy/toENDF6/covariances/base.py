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

from fudge.gnd.covariances.base import covarianceMatrix
from fudge.core.math import linearAlgebra

from .. import endfFormats
from xData import array as arrayModule
from xData import axes as axesModule

def toENDF6(self, flags, targetInfo, inCovarianceGroup=False):

    endf = []
    rowdat, coldat = targetInfo['dataPointer']
    MF,MT1 = map(int, rowdat['ENDF_MFMT'].split(','))
    if not inCovarianceGroup:
        # print header for this subsection (contains one NL sub-subsection)
        MAT1 = targetInfo['MAT1']
        XMF1,XLFS1,NC,NI = 0,0,0,1
        if coldat:
            MF1, MT1 = map(int, coldat['ENDF_MFMT'].split(','))
        if MF in (31,33):
            endf.append( endfFormats.endfHeadLine(XMF1,XLFS1,MAT1,MT1,NC,NI) )
    # header for matrix:
    rows,cols = self.matrix.array.shape
    if isinstance( self.matrix.array, arrayModule.diagonal ):
        LS = 0; LB = 1; NP = len(self.matrix.axes[2].grid); NT = 2*NP
        if self.type=='absolute': LB = 0
        if self.ENDFconversionFlag:
            LB = int( self.ENDFconversionFlag.split('=')[1] )
        matrixData = zip( self.matrix.axes[2].grid, list(self.matrix.array.values) + [0] )
        matrixData = [val for sublist in matrixData for val in sublist] # flatten
    elif isinstance( self.matrix.array, arrayModule.full ):
        if self.matrix.array.symmetry in (arrayModule.symmetryLowerToken, arrayModule.symmetryUpperToken):
            LS = 1; LB = 5; NT = (rows+1) + rows*(rows+1)/2; NP = rows+1
            arrayData = list( self.matrix.array.values )
            if self.matrix.array.symmetry == arrayModule.symmetryLowerToken:
                arrayData = linearAlgebra.switchSymmetry( arrayData, upperToLower=False )
            matrixData = list(self.matrix.axes[2].grid) + arrayData
        elif self.matrix.axes[1].style == axesModule.linkGridToken:
            LS = 0; LB = 5; NT = (rows+1) + rows*cols; NP = rows+1
            matrixData = list(self.matrix.axes[2].grid) + list(self.matrix.array.values)
        else:
            LS = 0; LB = 6; NT = (rows+1) + (cols+1) + rows*cols; NP = rows+1
            matrixData = list(self.matrix.axes[2].grid) + list(self.matrix.axes[1].grid) + list(
                    self.matrix.array.values)
    else:
        raise NotImplemented
    if MF==35:  # header for fission spectra is different:
        E1,E2 = [a.getValueAs('eV') for a in (rowdat.attributes['incidentEnergyLowerBound'],
                                              rowdat.attributes['incidentEnergyUpperBound'])]
        if LS: LB = 7
        else:
            raise Exception ("Unknown spectrum (MF35) covariance format")
        endf.append( endfFormats.endfHeadLine( E1,E2,LS,LB,NT,NP ) )
    else:
        endf.append( endfFormats.endfHeadLine( 0,0,LS,LB,NT,NP ) )
    endf += endfFormats.endfDataList( matrixData )
    return endf

covarianceMatrix.toENDF6 = toENDF6
