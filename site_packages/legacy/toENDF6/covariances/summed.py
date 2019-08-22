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

from fudge.gnd.covariances.summed import summedCovariance
from .. import endfFormats

def toENDF6(self, flags, targetInfo, inCovarianceGroup=False):
    endf = []
    if not inCovarianceGroup:
        # print header for this subsection (contains one NL sub-subsection)
        XMF1,XLFS1,MAT1,NC,NI = 0,0,0,1,0
        rowdat, coldat = targetInfo['dataPointer']
        MT1 = map(int, rowdat['ENDF_MFMT'].split(',')) [1]
        endf.append( endfFormats.endfHeadLine(XMF1,XLFS1,MAT1,MT1,NC,NI) )
    # header:
    LTY=0
    endf.append( endfFormats.endfHeadLine(0,0,0,LTY,0,0) )
    NCI = len(self.pointerList)
    endf.append( endfFormats.endfHeadLine(self.lowerBound.getValueAs('eV'),
        self.upperBound.getValueAs('eV'),0,0, 2*NCI,NCI) )
    mtList = [ map(int, a['ENDF_MFMT'].split(','))[1] for a in self.pointerList]
    coefficients = [a['coefficient'] for a in self.pointerList]
    endf += endfFormats.endfDataList( [i for j in zip(coefficients,mtList)
        for i in j] )
    return endf

summedCovariance.toENDF6 = toENDF6
