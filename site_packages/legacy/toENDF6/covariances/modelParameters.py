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

from pqu import PQU
from fudge.gnd.covariances.modelParameters import resonanceParameterCovariance, inputParameter
from .. import endfFormats

def toENDF6(self, endfMFList, flags, targetInfo, verbosityIndent=''):
    """ go back to ENDF format """
    def swaprows( matrix, i1, i2, nrows ):
        # may need to rearrange parameters: ENDF often sorts first by L rather than by energy
        rows = matrix[i1:i1+nrows].copy()
        matrix[i1:i1+nrows] = matrix[i2:i2+nrows]; matrix[i2:i2+nrows] = rows
        cols = matrix[:,i1:i1+nrows].copy()
        matrix[:,i1:i1+nrows] = matrix[:,i2:i2+nrows]; matrix[:,i2:i2+nrows] = cols

    # need the resonance parameters as well as covariance matrix:
    res = targetInfo['reactionSuite'].resonances
    RPs = res.resolved.evaluated.resonanceParameters.table
    NRes = self.inputParameters[-1].nResonances

    # MF32 header information:
    ZAM, AWT = targetInfo['ZA'], targetInfo['mass']
    NIS, ABN, ZAI = 1, 1.0, ZAM  # assuming only one isotope per file
    endf = [endfFormats.endfHeadLine( ZAM, AWT, 0, 0, NIS, 0 )]
    LFW = RPs.getColumn('fissionWidthA') is not None; NER=1
    endf.append( endfFormats.endfHeadLine( ZAI,ABN,0,LFW,NER,0 ) )
    EL,EH = res.resolved.lowerBound.getValueAs('eV'), res.resolved.upperBound.getValueAs('eV')
    LRU,NRO =1,0
    LRF = {'SingleLevel_BreitWigner':1, 'MultiLevel_BreitWigner':2, 'Reich_Moore':3}[
            res.resolved.evaluated.moniker ]
    NAPS = not res.resolved.evaluated.calculateChannelRadius
    endf.append( endfFormats.endfHeadLine( EL,EH,LRU,LRF,NRO,NAPS ) )
    SPI = targetInfo['spin']
    AP = res.resolved.evaluated.scatteringRadius.getValueAs('10*fm')
    LCOMP=1
    if 'LCOMP=0' in self.attributes.get('endfConversionFlags',''): LCOMP=0
    elif 'LCOMP=2' in self.attributes.get('endfConversionFlags',''): LCOMP=2

    sortByL = ("sortByL" in self.attributes.get('endfConversionFlags',''))
    Ls = RPs.getColumn('L')
    NLS = len(set(Ls))
    if LCOMP==2 or not sortByL: NLS = 0
    ISR = int( isinstance(self.inputParameters[0], inputParameter) and
            ('scatteringRadius' in self.inputParameters[0].name) )
    endf.append( endfFormats.endfHeadLine( SPI,AP,0,LCOMP,NLS,ISR ) )
    MLS = 0
    if ISR:
        MLS = 1 # currently don't handle energy-dependent DAP
        DAP = PQU.PQU( self.matrix.data[0][0], self.inputParameters[0].unit ).getValueAs('10*fm')
        if LRF in (1,2):
            endf.append( endfFormats.endfDataLine( [0,DAP] ) )
        elif LRF==3:
            endf.append( endfFormats.endfHeadLine( 0,0,0,0,MLS,1 ) )
            endf.append( endfFormats.endfDataLine( [DAP] ) )
        else:
            raise Exception("ISR>0 not yet supported for LRF=%i!" % LRF)

    # MF32 repeats the resonance parameter information.
    # Extract that info from reactionSuite.resonances:
    table = [RPs.getColumn('L'), RPs.getColumn('energy',units='eV'), RPs.getColumn('J'),
            RPs.getColumn('totalWidth',units='eV') or [0]*NRes,
            RPs.getColumn('neutronWidth',units='eV'), RPs.getColumn('captureWidth',units='eV'),
            RPs.getColumn('fissionWidthA') or [0]*NRes,
            RPs.getColumn('fissionWidthB') or [0]*NRes]
    CS = RPs.getColumn('channelSpin')
    if CS is not None:  # ENDF hack: J<0 -> use lower available channel spin
        CS = [2*(cs-SPI) for cs in CS]
        Js = [v[0]*v[1] for v in zip(table[2],CS)]
        table[2] = Js
    table = zip(*table)
    matrix = self.matrix.data[MLS:,MLS:].copy()
    MPAR = len(matrix) / len(table)

    if sortByL:
        # reorder resonances, sorting first by L and second by energy:
        table.sort()

        elist1 = [(lis[1],lis[4],lis[5]) for lis in table]
        elist2 = zip( RPs.getColumn('energy',units='eV'),
                RPs.getColumn('neutronWidth',units='eV'),
                RPs.getColumn('captureWidth',units='eV') )

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
            endf.append( endfFormats.endfHeadLine( AWT, 0, L, 0, 18*NRS, NRS ) )
            for i in range(tableIndex, len(table)):
                if table[i][0]!=L: break
                endf.append( endfFormats.endfDataLine( table[i][1:7] ) )
                block = matrix[MPAR*i:MPAR*(i+1), MPAR*i:MPAR*(i+1)]
                lis = [block[0,0], block[1,1], block[2,1], block[2,2]]
                if MPAR==4:
                    lis += [block[3,1],block[3,2],block[3,3],0,0,0,0,0]
                else:
                    lis += [0,0,0,0,0,0,0,0]
                endf += endfFormats.endfDataList( lis )
            tableIndex += NRS


    if LCOMP==1:
        NSRS, NLRS = 1,0    # short-range correlations only
        endf.append( endfFormats.endfHeadLine( AWT, 0, 0, 0, NSRS, NLRS ) )
        MPAR = len( self.inputParameters[0].parametersPerResonance.split(',') )
        NRB = NRes
        NVS = (NRB*MPAR)*(NRB*MPAR+1)/2 # length of the upper diagonal matrix
        endf.append( endfFormats.endfHeadLine( 0,0, MPAR, 0, NVS+6*NRB, NRB ) )

        for res in table:
            if LRF in (1,2):
                endf.append( endfFormats.endfDataLine( res[1:7] ) )
            elif LRF==3:
                endf.append( endfFormats.endfDataLine( res[1:3] + res[4:8] ) )

        dataList = []
        for i in range(len(matrix)): dataList.extend( list( matrix[i][i:] ) )
        endf += endfFormats.endfDataList( dataList )

    elif LCOMP==2:
        import numpy
        QX, LRX = 0, 0  # haven't encountered any competitive widths yet
        endf.append( endfFormats.endfHeadLine( AWT,QX,0,LRX, 12*NRes, NRes ) )
        dat = matrix.diagonal()
        for i in range(len(table)):
            if LRF in (1,2):
                params = table[i][1:7]
                uncerts = [dat[MPAR*i],0,0,dat[MPAR*i+1],dat[MPAR*i+2],0]
                if MPAR==4: uncerts[-1] = dat[MPAR*i+3]
            elif LRF==3:
                params = table[i][1:3] + table[i][4:8]
                uncerts = [dat[MPAR*i],0,dat[MPAR*i+1],dat[MPAR*i+2],0,0]
                if MPAR==5: uncerts[-2:] = [dat[MPAR*i+3], dat[MPAR*i+4]]
            endf += endfFormats.endfDataList( params )
            endf += endfFormats.endfDataList( uncerts )

        # correlation matrix:
        NDIGIT = [a for a in self.attributes['endfConversionFlags'].split(',') if a.startswith('NDIGIT')]
        NDIGIT = int( NDIGIT[0][-1] )
        nints = 56 // (NDIGIT+1)    # how many numbers fit on each line?
        if NDIGIT==3: nints = 13    # special case
        rsd = numpy.sqrt( matrix.diagonal() )
        rsd[ rsd==0 ] = 1
        corr_mat = matrix / numpy.outer( rsd,rsd )
        corr_mat = numpy.rint( corr_mat * 10**NDIGIT )  # rint: round to nearest int
        # write lower-diagonal as sparse matrix using INTG format:
        endfCorrMat = []
        for i in range(len(corr_mat)):
            vals = corr_mat[i,:i]
            j = 0
            while j < i:
                if vals[j]!=0:
                    endfCorrMat.append( endfFormats.writeEndfINTG(
                        i+1, j+1, list(vals[j:j+nints]), NDIGIT ) )
                    j += nints
                else: j+=1
        NNN = NRes * MPAR
        NM = len(endfCorrMat)
        endf.append( endfFormats.endfHeadLine( 0,0, NDIGIT, NNN, NM, 0 ) )
        endf += endfCorrMat
    endf.append( endfFormats.endfSENDLineNumber() )
    endfMFList[32][151] = endf

resonanceParameterCovariance.toENDF6 = toENDF6
