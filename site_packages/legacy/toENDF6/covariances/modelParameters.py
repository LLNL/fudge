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

import numpy

from pqu import PQU as PQUModule

from fudge.gnds.covariances import covarianceSuite as covarianceSuiteModule, modelParameters as modelParametersModule
from fudge.gnds import resonances as resonancesModule

from .. import endfFormats as endfFormatsModule
from site_packages.legacy.toENDF6 import resonances as resonancesRewriteModule

#
# helper methods:
#
def writeLCOMP2( matrix, NDIGIT, NNN ):
    import numpy
    nints = 56 // (NDIGIT + 1)  # how many numbers fit on each line?
    if NDIGIT == 3: nints = 13  # special case
    rsd = numpy.sqrt(matrix.diagonal())
    rsd[rsd == 0] = 1
    corr_mat = matrix / numpy.outer(rsd, rsd)
    corr_mat = numpy.rint(corr_mat * 10 ** NDIGIT)  # rint: round to nearest int
    # write lower-diagonal as sparse matrix using INTG format:
    endfCorrMat = []
    for i in range(len(corr_mat)):
        vals = corr_mat[i, :i]
        j = 0
        while j < i:
            if vals[j] != 0:
                endfCorrMat.append(endfFormatsModule.writeEndfINTG(
                    i + 1, j + 1, list(vals[j:j + nints]), NDIGIT))
                j += nints
            else:
                j += 1
    NM = len(endfCorrMat)
    endf = [endfFormatsModule.endfHeadLine(0, 0, NDIGIT, NNN, NM, 0)]
    endf += endfCorrMat
    return endf

def toENDF6(self, endfMFList, flags, targetInfo, verbosityIndent=''):
    """
    Write all unresolved covariance data to ENDF6
    """

    resolved, unresolved = [],[]
    for section_ in self:
        if isinstance(section_, modelParametersModule.averageParameterCovariance):
            unresolved.append(section_)
        else:
            resolved.append(section_)

    NER = len(resolved)
    if unresolved: NER += 1

    res = targetInfo['reactionSuite'].resonances
    if resolved:
        if isinstance( res.resolved.evaluated, resonancesModule.BreitWigner ):
            LFW = res.resolved.evaluated.resonanceParameters.table.getColumn('fissionWidth') is not None
        elif isinstance( res.resolved.evaluated, resonancesModule.RMatrix ):
            LFW = any( [RR.reactionLink.link.outputChannel.isFission() for RR in res.resolved.evaluated.resonanceReactions])
    else:
        LFW = any( [RR.reactionLink.link.outputChannel.isFission() for RR in res.unresolved.evaluated.resonanceReactions])

    # MF32 header information:
    ZAM, AWT = targetInfo['ZA'], targetInfo['mass']
    NIS, ABN, ZAI = 1, 1.0, ZAM  # assuming only one isotope per file
    endf = [endfFormatsModule.endfHeadLine(ZAM, AWT, 0, 0, NIS, 0)]
    endf.append(endfFormatsModule.endfHeadLine(ZAI, ABN, 0, LFW, NER, 0))

    endfMFList[32][151] = endf

    for section_ in resolved:
        section_[targetInfo['style']].toENDF6(endfMFList, flags, targetInfo, verbosityIndent)
    if unresolved:      # these need to all be done together
        averageParametersToENDF6(unresolved, endfMFList, flags, targetInfo, verbosityIndent)

    endfMFList[32][151].append( endfFormatsModule.endfSENDLineNumber() )

covarianceSuiteModule.parameterCovariances.toENDF6 = toENDF6


def toENDF6(self, endfMFList, flags, targetInfo, verbosityIndent=''):
    """
    Translate resolved resonance covariance back to ENDF-6
    """
    import numpy

    def swaprows( matrix, i1, i2, nrows ):
        # may need to rearrange parameters: ENDF often sorts first by L rather than by energy
        rows = matrix[i1:i1+nrows].copy()
        matrix[i1:i1+nrows] = matrix[i2:i2+nrows]; matrix[i2:i2+nrows] = rows
        cols = matrix[:,i1:i1+nrows].copy()
        matrix[:,i1:i1+nrows] = matrix[:,i2:i2+nrows]; matrix[:,i2:i2+nrows] = cols

    # need the resonance parameters as well as covariance matrix:
    res = targetInfo['reactionSuite'].resonances
    conversionFlags = targetInfo['ENDFconversionFlags'].get(self,"")

    LRF = 7
    if isinstance( res.resolved.evaluated, resonancesModule.BreitWigner ):
        LRF = {
            resonancesModule.BreitWigner.singleLevel: 1,
            resonancesModule.BreitWigner.multiLevel: 2
        }[ res.resolved.evaluated.approximation ]
    elif 'LRF3' in conversionFlags:
        LRF = 3

    endf = []
    EL, EH = res.resolved.domainMin, res.resolved.domainMax
    AWT, LRU, NRO = targetInfo['mass'], 1, 0
    NAPS = not res.resolved.evaluated.calculateChannelRadius
    endf.append(endfFormatsModule.endfHeadLine(EL, EH, LRU, LRF, NRO, NAPS))

    LCOMP = 1
    if 'LCOMP=0' in conversionFlags: LCOMP = 0
    elif 'LCOMP=2' in conversionFlags: LCOMP = 2

    if LRF in (1,2):
        RPs = res.resolved.evaluated.resonanceParameters.table
        NRes = len(RPs)

        SPI = targetInfo['spin']
        AP = res.resolved.evaluated.scatteringRadius.getValueAs('10*fm')

        sortByL = "sortByL" in conversionFlags
        Ls = RPs.getColumn('L')
        NLS = len(set(Ls))
        LAD = 0
        if LCOMP in (1,2) or not sortByL: NLS = 0
        ISR = any( [isinstance(parameter.link, resonancesModule.scatteringRadius) for parameter in self.parameters] )
        endf.append( endfFormatsModule.endfHeadLine( SPI,AP,LAD,LCOMP,NLS,ISR ) )
        MLS = 0
        if ISR:
            MLS = 1 # currently don't handle energy-dependent DAP
            DAP = PQUModule.PQU( self.matrix.constructArray()[0,0],
                self.parameters[0].link.form.axes[0].unit ).getValueAs('10*fm')
            endf.append( endfFormatsModule.endfContLine( 0,DAP,0,0,0,0 ) )

        # MF32 repeats the resonance parameter information.
        # Extract that info from reactionSuite.resonances:
        table = [RPs.getColumn('L'), RPs.getColumn('energy',unit='eV'), RPs.getColumn('J'),
                RPs.getColumn('totalWidth',unit='eV') or [0]*NRes,
                RPs.getColumn('neutronWidth',unit='eV'), RPs.getColumn('captureWidth',unit='eV'),
                RPs.getColumn('fissionWidth') or [0]*NRes]
        CS = RPs.getColumn('channelSpin')
        if CS is not None:  # ENDF hack: J<0 -> use lower available channel spin
            CS = [2*(cs-SPI) for cs in CS]
            Js = [v[0]*v[1] for v in zip(table[2],CS)]
            table[2] = Js
        table = zip(*table)
        matrix = self.matrix.constructArray()[MLS:,MLS:]
        # toss out extra rows/columns of zeros (for column 'L')
        MPAR2 = len(matrix) / len(table)
        index = []
        for idx in range(len(table)):
            index += [idx*MPAR2+1, idx*MPAR2+2, idx*MPAR2+3]    # indices to remove
        keep = numpy.array( sorted( set(range(len(matrix))).difference(index)) )
        matrix = matrix[keep.reshape(-1,1),keep]
        MPAR = len(matrix) / len(table)

        if sortByL:
            # reorder resonances, sorting first by L and second by energy:
            table.sort()

            elist1 = [(lis[1],lis[4],lis[5]) for lis in table]
            elist2 = zip( RPs.getColumn('energy',unit='eV'),
                    RPs.getColumn('neutronWidth',unit='eV'),
                    RPs.getColumn('captureWidth',unit='eV') )

            for i in range(len(elist1)):
                i2 = elist2.index( elist1[i] )
                if i2!=i:
                    swaprows( matrix, MPAR*i, MPAR*elist2.index( elist1[i] ), MPAR )
                    val = elist2[i]
                    elist2[i] = elist2[i2]; elist2[i2] = val

        if LCOMP==0:
            tableIndex = 0
            for L in set( Ls ):
                NRS = Ls.count(L)
                endf.append( endfFormatsModule.endfHeadLine( AWT, 0, L, 0, 18*NRS, NRS ) )
                for i in range(tableIndex, len(table)):
                    if table[i][0]!=L: break
                    endf.append( endfFormatsModule.endfDataLine( table[i][1:7] ) )
                    block = matrix[MPAR*i:MPAR*(i+1), MPAR*i:MPAR*(i+1)]
                    lis = [block[0,0], block[1,1], block[2,1], block[2,2]]
                    if MPAR==4:
                        lis += [block[3,1],block[3,2],block[3,3],0,0,0,0,0]
                    else:
                        lis += [0,0,0,0,0,0,0,0]
                    endf += endfFormatsModule.endfDataList( lis )
                tableIndex += NRS


        if LCOMP==1:
            NSRS, NLRS = 1,0    # short-range correlations only
            endf.append( endfFormatsModule.endfHeadLine( AWT, 0, 0, 0, NSRS, NLRS ) )
            NRB = NRes
            NVS = (NRB*MPAR)*(NRB*MPAR+1)/2 # length of the upper diagonal matrix
            endf.append( endfFormatsModule.endfHeadLine( 0,0, MPAR, 0, NVS+6*NRB, NRB ) )

            for res in table:
                endf.append( endfFormatsModule.endfDataLine( res[1:7] ) )

            dataList = []
            for i in range(len(matrix)): dataList.extend( list( matrix[i][i:] ) )
            endf += endfFormatsModule.endfDataList( dataList )

        elif LCOMP==2:
            QX, LRX = 0, 0  # haven't encountered any competitive widths yet
            dat = matrix.diagonal()
            NRes = 0
            matrixDat = []
            omitResonance = []
            for i in range(len(table)):
                params = table[i][1:7]
                uncerts = [dat[MPAR*i],0,0,dat[MPAR*i+1],dat[MPAR*i+2],0]
                if MPAR==4: uncerts[-1] = dat[MPAR*i+3]
                if not any(uncerts):
                    omitResonance.extend([True]*MPAR)  # resonance is not included in MF=32
                    continue
                omitResonance.extend([False]*MPAR)
                matrixDat += endfFormatsModule.endfDataList( params )
                matrixDat += endfFormatsModule.endfDataList( uncerts )
                NRes += 1

            endf.append(endfFormatsModule.endfHeadLine(AWT, QX, 0, LRX, 12 * NRes, NRes))
            endf += matrixDat

            # correlation matrix:
            if any(omitResonance):  # see Pt192
                omitResonance = numpy.array(omitResonance)
                matrix = matrix[~omitResonance][:,~omitResonance]
            NDIGIT = [a for a in conversionFlags.split(',') if a.startswith('NDIGIT')]
            NDIGIT = int( NDIGIT[0][-1] )
            endf += writeLCOMP2( matrix, NDIGIT, NRes*MPAR )

    elif LRF==3:
        conversionDetails = targetInfo['LRF3conversion']    # useful info saved when writing MF=2 back to ENDF6

        table = conversionDetails['table']
        sortedTable = conversionDetails['sortedTable']
        NRes = len(table['energies'])

        SPI = targetInfo['spin']
        AP = conversionDetails['AP']

        sortByL = "sortByL" in conversionFlags
        Ls = table['Ls']
        NLS = len(set(Ls))
        LAD = 0
        if LCOMP==2 and res.resolved.evaluated.supportsAngularReconstruction: LAD=1
        if LCOMP==2 or not sortByL: NLS = 0
        ISR = any( [isinstance(parameter.link, resonancesModule.scatteringRadius) for parameter in self.parameters] )
        endf.append( endfFormatsModule.endfHeadLine( SPI,AP,LAD,LCOMP,NLS,ISR ) )
        MLS = 0
        matrix = self.matrix.constructArray()
        if ISR:
            MLS = 1 # currently don't handle energy-dependent DAP
            DAP = PQUModule.PQU( matrix[0,0], self.parameters[0].link.form.axes[0].unit ).getValueAs('10*fm')
            endf.append( endfFormatsModule.endfHeadLine( 0,0,0,0,MLS,1 ) )
            endf.append( endfFormatsModule.endfDataLine( [DAP] ) )

        matrix = matrix[MLS:,MLS:]
        MPAR = len(matrix) / NRes

        if not sortByL:
            sortedTable.sort(key=lambda foo: foo[1])    # sort by resonance energy. Otherwise sorted first by L, then E

        elist1 = [(lis[1],lis[3],lis[4]) for lis in sortedTable]
        elist2 = zip( table['energies'], table['elastic'], table['capture'] )

        for i in range(len(elist1)):
            i2 = elist2.index( elist1[i] )
            if i2!=i:
                swaprows( matrix, MPAR*i, MPAR*elist2.index( elist1[i] ), MPAR )
                val = elist2[i]
                elist2[i] = elist2[i2]; elist2[i2] = val

        for i1 in range(len(elist1)):   # switch order of elastic and capture widths
            swaprows(matrix, MPAR * i1 + 1, MPAR * i1 + 2, 1)

        omitRow = numpy.sum(matrix, axis=1) == 0

        # check for empty rows/columns:
        if any(omitRow):
            omitResonance = []
            for idx in range(len(sortedTable)):
                if numpy.all(omitRow[MPAR*idx:MPAR*(idx+1)]):
                    omitResonance.append(idx)
                else:
                    omitRow[MPAR*idx:MPAR*(idx+1)] = 0  # only omit of all parameters are zero for this resonance
            for idx in omitResonance[::-1]:
                sortedTable.pop(idx)
            NRes = len(sortedTable)
            matrix = matrix[~omitRow][:, ~omitRow]

        if LCOMP==0:
            tableIndex = 0
            for L in set( Ls ):
                NRS = Ls.count(L)
                endf.append( endfFormatsModule.endfHeadLine( AWT, 0, L, 0, 18*NRS, NRS ) )
                for i in range(tableIndex, len(sortedTable)):
                    if sortedTable[i][0]!=L: break
                    endf.append( endfFormatsModule.endfDataLine( sortedTable[i][1:7] ) )
                    block = matrix[MPAR*i:MPAR*(i+1), MPAR*i:MPAR*(i+1)]
                    lis = [block[0,0], block[1,1], block[2,1], block[2,2]]
                    if MPAR==4:
                        lis += [block[3,1],block[3,2],block[3,3],0,0,0,0,0]
                    else:
                        lis += [0,0,0,0,0,0,0,0]
                    endf += endfFormatsModule.endfDataList( lis )
                tableIndex += NRS

        elif LCOMP==1:
            NSRS, NLRS = 1,0    # short-range correlations only
            endf.append( endfFormatsModule.endfHeadLine( AWT, 0, 0, 0, NSRS, NLRS ) )
            NRB = NRes
            NVS = (NRB*MPAR)*(NRB*MPAR+1)/2 # length of the upper diagonal matrix
            endf.append( endfFormatsModule.endfHeadLine( 0,0, MPAR, 0, NVS+6*NRB, NRB ) )

            for res in sortedTable:
                endf.append( endfFormatsModule.endfDataLine( res[1:] ) )

            dataList = []
            for i in range(len(matrix)): dataList.extend( list( matrix[i][i:] ) )
            endf += endfFormatsModule.endfDataList( dataList )

        elif LCOMP==2:
            QX, LRX = 0, 0  # haven't encountered any competitive widths yet
            dat = matrix.diagonal()
            NRes = 0
            matrixDat = []
            for i in range(len(sortedTable)):
                params = sortedTable[i][1:]
                uncerts = [dat[MPAR*i],0,dat[MPAR*i+1],dat[MPAR*i+2],0,0]
                if MPAR==5: uncerts[-2:] = [dat[MPAR*i+3], dat[MPAR*i+4]]
                if not any(uncerts): continue  # Some resonances may not be included in MF=32
                matrixDat += endfFormatsModule.endfDataList( params )
                matrixDat += endfFormatsModule.endfDataList( uncerts )
                NRes += 1

            endf.append(endfFormatsModule.endfHeadLine(AWT, QX, 0, LRX, 12 * NRes, NRes))
            endf += matrixDat

            # correlation matrix:
            NDIGIT = [a for a in conversionFlags.split(',') if a.startswith('NDIGIT')]
            NDIGIT = int( NDIGIT[0][-1] )
            endf += writeLCOMP2( matrix, NDIGIT, NRes*MPAR )

    else:   # LRF = 7
        import numpy
        matrix = self.matrix.constructArray()
        RML = res.resolved.evaluated

        IFG = int( RML.reducedWidthAmplitudes )
        NJS = len( self.parameters )
        ISR = 0 # FIXME: scattering radius uncertainty not yet handled

        if LCOMP==1:
            AWRI, NSRS, NLRS = 0, 1, 0    # FIXME: hard-coded
            NJSX = len( RML.spinGroups )
            NPARB = 0

            endf.append(endfFormatsModule.endfContLine(0, 0, 0, LCOMP, 0, ISR))
            endf.append( endfFormatsModule.endfContLine(AWRI, 0, 0, 0, NSRS, NLRS) )
            endf.append( endfFormatsModule.endfContLine(0,0,NJSX,0,0,0) )
            for spingrp in RML.spinGroups:
                NCH = len( spingrp.channels )
                NRB = len( spingrp.resonanceParameters.table )
                NX = (NCH//6 + 1)*NRB
                endf.append( endfFormatsModule.endfContLine(0,0,NCH,NRB,6*NX,NX) )
                for res in spingrp.resonanceParameters.table:
                    for jidx in range(NCH // 6 + 1):
                        endfLine = res[jidx * 6:jidx * 6 + 6]
                        while len(endfLine) < 6: endfLine.append(0)
                        endf.append(endfFormatsModule.endfDataLine(endfLine))
                if NRB == 0:
                    endf.append(endfFormatsModule.endfDataLine([0, 0, 0, 0, 0, 0]))
                NPARB += NRB * (NCH+1)

            # matrix header
            N = (NPARB * (NPARB+1))/2
            endf.append( endfFormatsModule.endfContLine(0,0,0,0,N,NPARB))
            dataList = []
            for i in range(len(matrix)): dataList.extend( list( matrix[i][i:] ) )   # upper-diagonal
            endf += endfFormatsModule.endfDataList( dataList )

        elif LCOMP==2:
            uncertainties = list(numpy.sqrt(matrix.diagonal()))
            endf.append(endfFormatsModule.endfContLine(0, 0, IFG, LCOMP, NJS, ISR))
            endf.extend( resonancesRewriteModule.writeRMatrixParticlePairs(RML, targetInfo) )

            uidx = 0
            NNN = 0
            for spingrp in RML.spinGroups:
                NRSA = len(spingrp.resonanceParameters.table)
                if NRSA==0: continue
                NCH, spinGroupHeader = resonancesRewriteModule.writeRMatrixSpinGroupHeader(RML, spingrp)
                NNN += NRSA * (NCH + 1) # +1 for resonance energy
                endf.extend( spinGroupHeader )

                # write resonance parameters (redundant with MF=2), followed by uncertainties
                NX = (NCH//6 + 1)*NRSA
                endf.append( endfFormatsModule.endfHeadLine( 0,0,0,NRSA,12*NX,NX ) )
                for res in spingrp.resonanceParameters.table:
                    for jidx in range(NCH//6+1):
                        endfLine = res[jidx*6:jidx*6+6]
                        while len(endfLine)<6: endfLine.append(0)
                        endf.append( endfFormatsModule.endfDataLine( endfLine ) )

                    uncertaintiesThisResonance = uncertainties[uidx:uidx + len(res)]
                    for jidx in range(NCH//6+1):
                        endfLine = uncertaintiesThisResonance[jidx*6:max(jidx*6+6, NCH)]
                        while len(endfLine)<6: endfLine.append(0)
                        endf.append( endfFormatsModule.endfDataLine( endfLine ) )
                    uidx += len(res)

                if NRSA==0:
                    endf.append( endfFormatsModule.endfDataLine( [0,0,0,0,0,0] ) )

            # correlation matrix:
            NDIGIT = [a for a in conversionFlags.split(',') if a.startswith('NDIGIT')]
            NDIGIT = int(NDIGIT[0][-1])
            endf += writeLCOMP2(matrix, NDIGIT, NNN)
        else:
            raise NotImplementedError("MF32 LRF7 with LCOMP=%d" % LCOMP)

    endfMFList[32][151] += endf

modelParametersModule.parameterCovarianceMatrix.toENDF6 = toENDF6


def averageParametersToENDF6( averageParameterSections, endfMFList, flags, targetInfo, verbosityIndent):
    """
    Unresolved resonance parameters need to be converted from multiple sections back into a single
    covariance matrix for ENDF-6
    """

    endf = []
    res = targetInfo['reactionSuite'].resonances
    URR = res.unresolved
    EL,EH = URR.domainMin, URR.domainMax
    LRU,LRF,NRO,NAPS = 2,1,0,0
    endf.append(endfFormatsModule.endfContLine(EL,EH,LRU,LRF,NRO,NAPS))

    SPI = targetInfo['spin']
    if res.unresolved.evaluated.scatteringRadius is not None:
        AP = res.unresolved.evaluated.scatteringRadius.getValueAs('10*fm')
    else:
        AP = res.scatteringRadius.getValueAs('10*fm')

    NLS = len(URR.evaluated.Ls)
    AWRI = targetInfo['mass']
    endf.append(endfFormatsModule.endfHeadLine(SPI, AP, 0, 0, NLS, 0))

    MPARs, diagonal, diagonalPointers = [], [], []
    for lsection in URR.evaluated.Ls:
        NJS = len( lsection.Js )
        endf.append(endfFormatsModule.endfHeadLine(AWRI,0,lsection.L,0,6*NJS,NJS))

        for jsection in lsection.Js:

            params = {'D':0, 'AJ':jsection.J, 'GNO':0, 'GG':0, 'GF':0, 'GX':0}
            MPARs.append( len(jsection.widths) )
            for width in jsection.widths:
                covarianceMatrix = width.data.uncertainty.data.link

                savedParams = targetInfo['ENDFconversionFlags'].get(covarianceMatrix,"").split(',')
                for p in savedParams:
                    key, value = p.split('=')
                    params[key] = float(value)

                if covarianceMatrix.matrix.array.shape != (1,1):
                    raise NotImplementedError("ENDF-6 format does not support energy-dependent unresolved covariances")
                diagonal.append( covarianceMatrix.matrix.array.values[0] )
                diagonalPointers.append( width )

            endf.append(endfFormatsModule.endfDataLine([params[key] for key in ('D','AJ','GNO','GG','GF','GX')]))

    assert len(set(MPARs)) == 1 # same number of widths for each L/J section
    MPAR = MPARs[0]
    NPAR = len(diagonal)

    endf.append(endfFormatsModule.endfContLine(0,0,MPAR,0,(NPAR*(NPAR+1))/2,NPAR))
    matrix = numpy.identity(NPAR) * diagonal

    crossTerms = [section for section in averageParameterSections if section.crossTerm]
    for crossTerm in crossTerms:
        ridx = diagonalPointers.index( crossTerm.rowData.link )
        cidx = diagonalPointers.index( crossTerm.columnData.link )
        matrix[ridx,cidx] = matrix[cidx,ridx] = crossTerm[targetInfo['style']].matrix.array.values[0]

    datalist = []
    for idx,row in enumerate(matrix):
        datalist += row[idx:].tolist()
    endf += endfFormatsModule.endfDataList(datalist)

    endfMFList[32][151] += endf
