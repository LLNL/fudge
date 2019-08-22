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

from xData import standards as standardsModule

from PoPs.groups import misc as chemicalElementMiscPoPsModule
from PoPs.families import nuclide as nuclideModule

from fudge.legacy.converting import endf_endl as endf_endlModule

from fudge.gnds.covariances import covarianceSuite as covarianceSuiteModule
from fudge.gnds.covariances import mixed as covarianceMixedModule

from .. import endfFormats as endfFormatsModule
from .modelParameters import averageParametersToENDF6

def toENDF6(self, endfMFList, flags, targetInfo, verbosityIndent=''):
    """Convert to ENDF format."""

    reactionSuite = targetInfo['reactionSuite']

    ZAM, AWT = targetInfo['ZA'], targetInfo['mass']
    NIS, ABN = 1, 1.0
    ZAI = ZAM  # assuming one isotope/file
    MTL = 0 # mtl=1 sections are handled in lumpedCovariance

    if self.parameterCovariances:
        self.parameterCovariances.toENDF6(endfMFList, flags, targetInfo, verbosityIndent)

    sections = self.covarianceSections
    # sort covariances by MF/MT:
    mfmts = []
    for section_ in sections:
        mfmts.append(map(int,section_.rowData['ENDF_MFMT'].split(',')))
    if len(mfmts)==0: return

    mfs, mts = zip(*mfmts)
    zipList = zip(mfs,mts,sections)
    idx = 0
    while idx<len(zipList):
        mf,mt,covar = zipList[idx]
        thisMFMT = [a[2] for a in zipList if a[:2]==(mf,mt)]
        idx += len(thisMFMT)

        if mf in (31,33):
            endf = [endfFormatsModule.endfHeadLine( ZAM, AWT, 0, MTL, 0, len(thisMFMT) )]
        elif mf==34:
            LTT = 1
            NMT1 = 1    # cross-reaction terms not yet supported
            endf = [endfFormatsModule.endfHeadLine( ZAM, AWT, 0, LTT, 0, NMT1 )]
            MAT1 = 0
            MT1 = mt
            L1s, L2s = [],[]
            for section_ in thisMFMT:
                L1s.append( section_.rowData['L'] )
                if section_.columnData:
                    L2s.append( section_.columnData['L'] )
                else:
                    L2s.append( L1s[-1] )
            NL = len(set(L1s))
            NL1 = len(set(L2s))
            endf += [ endfFormatsModule.endfHeadLine( 0.0, 0.0, MAT1, MT1, NL, NL1 ) ]
        elif mf==35:
            endf = [endfFormatsModule.endfHeadLine( ZAM, AWT, 0, MTL, len(thisMFMT), 0 )]
        elif mf==40:
            endf = [endfFormatsModule.endfHeadLine( ZAM, AWT, 0, 0, len(thisMFMT), 0 )]
        for section_ in thisMFMT:
            MAT1 = 0
            form = section_[ targetInfo['style'] ]
            conversionFlags = targetInfo['ENDFconversionFlags'].get(form, "")
            if (section_.columnData is not None and
                    section_.columnData.root is not None and
                    section_.columnData.root[1:] in self.externalFiles):
                otherTarget = section_.columnData.root[1:]
                if otherTarget != 'reactions':
                    ZA, MAT1 = endf_endlModule.ZAAndMATFromParticleName( otherTarget )
            if mf==34:
                L1 = int(section_.rowData['L'])
                if section_.columnData:
                    L2 = int(section_.columnData['L'])
                else:
                    L2 = L1
                NI = 1
                if isinstance(form, covarianceMixedModule.mixedForm):
                    NI = len(form)
                    frame = form[0].productFrame
                else:
                    frame = form.productFrame
                LCT = {standardsModule.frames.labToken: 1,
                       standardsModule.frames.centerOfMassToken: 2}[frame]
                if 'LCT=0' in conversionFlags:
                    LCT = 0
                endf += [ endfFormatsModule.endfHeadLine( 0.0, 0.0, L1, L2, LCT, NI ) ]
            if mf==40:
                rowData = section_.rowData
                if( isinstance( rowData, str ) ) :
                    raise Exception( "Don't string me along!" ) # FIXME
                else:
                    quant = rowData.link
                    if mt == 18:
                        QM, QI = quant.findAttributeInAncestry('getQ')('eV'), 0
                        LFS, level, IZAP = 0,0,-1
                    else:
                        product = quant.findAttributeInAncestry('outputChannel')[0]
                        QI = quant.findAttributeInAncestry('getQ')('eV')
                        LFS, level = 0, 0.
                        particle = reactionSuite.PoPs[product.id]
                        if( isinstance( particle, nuclideModule.particle ) ) :
                            LFS = particle.index
                            level = particle.energy[0].float( 'eV' )
                        QM = QI + level
                        IZAP = chemicalElementMiscPoPsModule.ZA( reactionSuite.PoPs[product.id] )
                    NL = 1
                    endf += [endfFormatsModule.endfHeadLine( QM, QI, IZAP, LFS, 0, NL ) ]
                    XMF1, XLFS1, NC, NI = 10,LFS, 0,1
                    endf += [endfFormatsModule.endfHeadLine( XMF1,XLFS1,MAT1,mt,NC,NI )]
            targetInfo['MAT1'] = MAT1
            targetInfo['dataPointer'] = [section_.rowData,section_.columnData]
            endf += form.toENDF6(flags, targetInfo)
            targetInfo.pop( 'dataPointer' )
            targetInfo.pop( 'MAT1' )
        endf.append( endfFormatsModule.endfSENDLineNumber() )
        if mt not in endfMFList[mf]:
            endfMFList[mf][mt] = []
        endfMFList[mf][mt] += endf
    # also add ENDF-style pointers for lumped covariance data:
    for reactionSum in targetInfo['reactionSuite'].sums.crossSections:
        MT1 = reactionSum.ENDF_MT
        if MT1 not in range(851,872): continue
        for summand in reactionSum.summands:
            mt = summand.link.findAttributeInAncestry('ENDF_MT')
            endfMFList[33][mt] = [endfFormatsModule.endfHeadLine(ZAM,AWT,0,MT1,0,0),
                    endfFormatsModule.endfSENDLineNumber() ]
    return

covarianceSuiteModule.covarianceSuite.toENDF6 = toENDF6
