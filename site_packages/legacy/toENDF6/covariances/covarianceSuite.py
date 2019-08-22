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

from PoPs import misc as miscPoPsModule
from PoPs.families import nuclearLevel as nuclearLevelModule

from fudge.legacy.converting import endf_endl as endf_endlModule
from fudge.gnd.covariances import covarianceSuite as covarianceSuiteModule
from fudge.gnd.covariances import section as sectionModule
from fudge.gnd.covariances import distributions as distributionsModule

from .. import endfFormats as endfFormatsModule

def toENDF6(self, endfMFList, flags, targetInfo, verbosityIndent=''):
    """Convert to ENDF format."""

    reactionSuite = targetInfo['reactionSuite']

    ZAM, AWT = targetInfo['ZA'], targetInfo['mass']
    NIS, ABN = 1, 1.0;
    ZAI = ZAM  # assuming one isotope/file
    MTL = 0 # mtl=1 sections are handled in lumpedCovariance

    for section_ in self.modelParameterCovariances:
        section_.toENDF6(endfMFList, flags, targetInfo, verbosityIndent)

    sections = self.sections[:]
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
            if isinstance( thisMFMT[0][ targetInfo['style'] ], distributionsModule.LegendreOrderCovarianceForm ): LTT = 1
            NMT1 = len(thisMFMT)
            endf = [endfFormatsModule.endfHeadLine( ZAM, AWT, 0, LTT, 0, NMT1 )]
        elif mf==35:
            endf = [endfFormatsModule.endfHeadLine( ZAM, AWT, 0, MTL, len(thisMFMT), 0 )]
        elif mf==40:
            endf = [endfFormatsModule.endfHeadLine( ZAM, AWT, 0, 0, len(thisMFMT), 0 )]
        for section_ in thisMFMT:
            MAT1 = 0
            form = section_[ targetInfo['style'] ]
            if section_.columnData and isinstance(section_.columnData.link, sectionModule.externalReaction):
                otherTarget = section_.columnData.link.target
                ZA, MAT1 = endf_endlModule.ZAAndMATFromParticleName( otherTarget )
            if mf==34:
                L1s = [subsec.L1 for subsec in form]
                L2s = [subsec.L2 for subsec in form]
                NL = len( set(L1s) );  NL1 = len( set(L2s) )
                if section_.columnData: raise NotImplemented # cross-reaction or cross-material
                MT1 = mt
                endf += [ endfFormatsModule.endfHeadLine( 0.0, 0.0, MAT1, MT1, NL, NL1 ) ]
            if mf==40:
                rowData = section_.rowData
                if( isinstance( rowData, str ) ) :
                    raise Exception( "Don't string me along!" ) # FIXME
                else:
                    quant = rowData.link
                    product = quant.findAttributeInAncestry('outputChannel')[0]
                    QI = quant.findAttributeInAncestry('getQ')('eV')
                    LFS, level = 0, 0.
                    particle = reactionSuite.PoPs[product.id]
                    if( isinstance( particle, nuclearLevelModule.particle ) ) :
                        LFS = particle.intIndex
                        level = particle.energy[0].float( 'eV' )
                    QM = QI + level
                    IZAP = miscPoPsModule.ZA( reactionSuite.PoPs[product.id] )
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
        if MT1 not in xrange(851,872): continue
        for summand in reactionSum.summands:
            mt = summand.link.findAttributeInAncestry('ENDF_MT')
            endfMFList[33][mt] = [endfFormatsModule.endfHeadLine(ZAM,AWT,0,MT1,0,0),
                    endfFormatsModule.endfSENDLineNumber() ]
    return

covarianceSuiteModule.covarianceSuite.toENDF6 = toENDF6
