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

import math

from pqu import PQU as PQUModule

from PoPs import misc as miscPoPsModule

import fudge.gnd.resonances as resonancesModule

from . import endfFormats as endfFormatsModule
from . import gndToENDF6 as gndToENDF6Module

#
# resonances
#
def toENDF6( self, endfMFList, flags, targetInfo, verbosityIndent = '' ) :

    ZAM, AWT = targetInfo['ZA'], targetInfo['mass']
    NIS, ABN = 1, 1.0; ZAI=ZAM  # assuming only one isotope per file

    # get target spin from the particle list:
    target = targetInfo['reactionSuite'].PoPs[targetInfo['reactionSuite'].target]
    if hasattr(target, 'nucleus'):
        targetInfo['spin'] = target.nucleus.spin[0].value
    else:
        targetInfo['spin'] = target.spin[0].value

    endf = [endfFormatsModule.endfHeadLine( ZAM, AWT, 0, 0, NIS, 0)]
    resolvedCount, unresolvedCount = 0, 0
    # resolved may have multiple energy regions:
    if self.resolved is not None: resolvedCount = max(1,len(self.resolved.regions))
    if self.unresolved is not None: unresolvedCount = 1

    # resonances may only contain a scattering radius:
    if not (resolvedCount + unresolvedCount) and self.scatteringRadius:
        scatRadius = self.scatteringRadius.form
        lowerBound, upperBound = scatRadius.domainMin, scatRadius.domainMax
        endf.append( endfFormatsModule.endfHeadLine( ZAM, ABN, 0,0,1,0 ) )
        endf.append( endfFormatsModule.endfHeadLine(lowerBound, upperBound, 0,0,0,0 ) )
        AP = PQUModule.PQU( self.scatteringRadius.form.constant, self.scatteringRadius.form.rangeUnit ).getValueAs('10*fm')
        endf.append( endfFormatsModule.endfHeadLine(targetInfo['spin'], AP, 0,0,0,0 ) )
        endf.append( endfFormatsModule.endfSENDLineNumber() )
        endfMFList[2][151] = endf
        return

    # For now I'm storing the LRF/LFW flags in xml, since they are tricky to compute
    # LFW is a pain: only applies to unresolved, but must be written at the start of MF2
    LRFurr, LFW = 0,0
    if unresolvedCount != 0:
        LRF_LFW = self.unresolved.evaluated.ENDFconversionFlag
        LRFurr, LFW = map(int, LRF_LFW.split('=')[-1].split(',') )
    NER = resolvedCount + unresolvedCount
    endf.append( endfFormatsModule.endfHeadLine( ZAI, ABN, 0, LFW, NER, 0 ) )
    for idx in range(resolvedCount):
        if resolvedCount==1: region = self.resolved
        else: region = self.resolved.regions[idx]
        LRU = 1 #resolved
        LRF = { resonancesModule.SLBW.moniker:1,
                resonancesModule.MLBW.moniker:2,
                resonancesModule.RMatrix.moniker:7
              }[region.evaluated.moniker]
        EL, EH = region.domainMin, region.domainMax
        if LRF==7:
            NRO = 0
            if region.evaluated.ENDFconversionFlag and 'LRF3' in region.evaluated.ENDFconversionFlag:
                LRF = 3
        else: NRO = region.evaluated.scatteringRadius.isEnergyDependent()
        NAPS = not( region.evaluated.calculateChannelRadius )
        endf.append(endfFormatsModule.endfHeadLine( EL,EH,LRU,LRF,NRO,NAPS ) )
        endf += region.evaluated.toENDF6( flags, targetInfo, verbosityIndent )
    if unresolvedCount != 0:
        LRU = 2 #unresolved
        region = self.unresolved
        EL, EH = region.domainMin, region.domainMax
        NRO, NAPS = 0,0
        if region.evaluated.scatteringRadius.isEnergyDependent(): NRO = 1
        endf.append(endfFormatsModule.endfHeadLine( EL,EH,LRU,LRFurr,NRO,NAPS ) )
        # pass LFW/LRF so we don't have to calculate twice:
        targetInfo['unresolved_LFW'] = LFW
        targetInfo['unresolved_LRF'] = LRFurr
        targetInfo['LSSF'] = not self.unresolved.reconstructCrossSection
        targetInfo['regionEnergyBounds'] = (region.domainMin, region.domainMax)
        endf += region.evaluated.toENDF6( flags, targetInfo, verbosityIndent )
    endf.append( endfFormatsModule.endfSENDLineNumber() )
    endfMFList[2][151] = endf

resonancesModule.resonances.toENDF6 = toENDF6

#
# resonanceFormalismBaseClass
#
def toENDF6( self, flags, targetInfo, verbosityIndent='' ):

    endf = []
    AP = getattr( self, 'scatteringRadius' )
    if AP.isEnergyDependent():
        scatRadius = AP.form
        NR, NP = 1, len(scatRadius)
        endf.append( endfFormatsModule.endfHeadLine( 0,0,0,0, NR,NP ) )
        endf += endfFormatsModule.endfInterpolationList( (NP,
            gndToENDF6Module.gndToENDFInterpolationFlag( scatRadius.interpolation ) ) )
        endf += endfFormatsModule.endfNdDataList( scatRadius.convertAxisToUnit( 0, '10*fm' ) )
        AP = 0
    else :
        AP = self.scatteringRadius.getValueAs( '10*fm' )
    L_list = self.resonanceParameters.table.getColumn('L')
    NLS = len( set(L_list) )
    LAD = getattr(self, 'computeAngularDistribution') or 0
    NLSC = getattr(self, 'LvaluesNeededForConvergence') or 0
    endf += [endfFormatsModule.endfHeadLine( targetInfo['spin'], AP, LAD, 0, NLS, NLSC )]

    table = [ self.resonanceParameters.table.getColumn('energy',unit='eV'),
            self.resonanceParameters.table.getColumn('J') ]
    NE = len(table[0])
    for attr in ('totalWidth','neutronWidth','captureWidth','fissionWidth'):
        column = self.resonanceParameters.table.getColumn( attr, unit='eV' )
        if not column: column = [0]*NE
        table.append( column )
    CS = self.resonanceParameters.table.getColumn('channelSpin')
    if CS is not None:  # ENDF hack: J<0 -> use lower available channel spin
        targetSpin = targetInfo['spin']
        CS = [2*(cs-targetSpin) for cs in CS]
        Js = [v[0]*v[1] for v in zip(table[1],CS)]
        table[1] = Js
    table = zip(*table)

    for L in set(L_list):
        APL = 0
        resonances = [table[i] for i in range(NE) if L_list[i] == L]
        NRS = len(resonances)
        endf.append( endfFormatsModule.endfHeadLine( targetInfo['mass'], APL, L, 0, 6*NRS, NRS ) )
        for row in resonances:
            endf.append( endfFormatsModule.endfDataLine( row ) )
    return endf

resonancesModule.resonanceFormalismBaseClass.toENDF6 = toENDF6

#
# helper functions for RMatrix
#
def getENDFtuple( spin, parity ):
    # ENDF combines spin & parity UNLESS spin==0. Then it wants (0,+/-1)
    if parity is None: raise ValueError("parity is None")
    if spin:
        return (spin * parity, 0)
    else:
        return (spin, parity)

def writeRMatrixParticlePairs( RMatrixLimited, targetInfo ):
    """
    Gets used for writing both MF=2 and MF=32
    :param RMatrixLimited: 
    :return: 
    """
    endf = []
    NPP = len(RMatrixLimited.resonanceReactions)
    endf.append( endfFormatsModule.endfHeadLine(0,0,NPP,0,12*NPP,2*NPP) )

    def MZIP( particle, ignoreMissingJpi=False ):  # helper: extract mass, z, spin and parity from particle list

        PoPs = targetInfo['reactionSuite'].PoPs
        try:
            nMass = PoPs['n'].getMass( 'amu' )
        except:
            nMass = 1.00866491578 # From ENDF102 manual, Appendix H.4

        mass = particle.getMass( 'amu' ) / nMass

        if hasattr(particle, 'nucleus'): particle = particle.nucleus
        Z = miscPoPsModule.ZAInfo( particle )[0]
        try:
            parity = particle.parity[0].value
            spin = particle.spin.float('hbar')
            I, P = getENDFtuple( spin, parity )
        except IndexError, err:
            if ignoreMissingJpi:
                I,P = (0,0) # ENDF omits spin/parity for compound nucleus from capture
            else:
                raise ValueError( "When formatting particle %s, encountered %s" % ( particle.id, err.message ) )
        return mass, Z, I, P

    for idx,pp in enumerate(RMatrixLimited.resonanceReactions):
        reaction = pp.reactionLink.link
# BRBBRB
        pA,pB = pp.label.split(' + ')
        # get the PoPs instances for pA and pB:
        pA,pB = targetInfo['reactionSuite'].getParticle(pA), targetInfo['reactionSuite'].getParticle(pB)
        MT = reaction.ENDF_MT
        MA, ZA, IA, PA = MZIP( pA )
        MB, ZB, IB, PB = MZIP( pB, ignoreMissingJpi=(MT == 102))
        PNT = RMatrixLimited.calculatePenetrability
        SHF = pp.computeShiftFactor
        if MT in (19,102): PNT = 0  # special case
        if pp.Q is not None: Q = pp.Q.getValueAs('eV')
        else:
            Q = reaction.getQ('eV')
            # getQ doesn't account for residual left in excited state:
            for particle in reaction.outputChannel:
                if( hasattr( particle, 'getLevelAsFloat' ) ) : Q -= particle.getLevelAsFloat('eV')
        if MT==102: Q = 0
        endf.append( endfFormatsModule.endfDataLine( [MA,MB,ZA,ZB,IA,IB] ) )
        endf.append( endfFormatsModule.endfDataLine( [Q,PNT,SHF,MT,PA,PB] ) )
        pp.index = idx+1    # 1-based index in ENDF

    return endf

def writeRMatrixSpinGroupHeader( RMatrixLimited, spingrp ):
    endf = []
    AJ, PJ = getENDFtuple(float(spingrp.spin), int(spingrp.parity))
    KBK = spingrp.background
    KPS = spingrp.applyPhaseShift
    NCH = len(spingrp.resonanceParameters.table.columns) - 1  # skip the 'energy' column
    try:
        endf.append(endfFormatsModule.endfHeadLine(AJ, PJ, KBK, KPS, 6 * NCH, NCH))
    except TypeError as err:
        raise TypeError("Got '%s' when formatting '%s'" % (err.message, str((AJ, PJ, KBK, KPS, 6 * NCH, NCH))))
    for chan in spingrp.channels:
        rreac = RMatrixLimited.resonanceReactions[chan.resonanceReaction]
        PPI = rreac.index
        L = chan.L
        if type(L) not in [int, float]: L = L.value
        SCH = chan.channelSpin
        if chan.boundaryConditionOverride is not None:
            BND = chan.boundaryConditionOverride
        elif RMatrixLimited.boundaryCondition == 'S':
            BND = 0
        elif RMatrixLimited.boundaryCondition == '-L':
            BND = -L

        APT = chan.scatteringRadius or rreac.scatteringRadius or PQUModule.PQU(0, 'fm')
        APE = chan.hardSphereRadius or rreac.hardSphereRadius or APT
        APT = APT.getValueAs('10*fm')
        APE = APE.getValueAs('10*fm')
        if rreac.reactionLink.link.ENDF_MT == 102:
            APT, APE = 0, 0
        endf.append(endfFormatsModule.endfDataLine([PPI, L, SCH, BND, APE, APT]))
    return NCH, endf

#
# RMatrix
#
def toENDF6( self, flags, targetInfo, verbosityIndent = '' ):

    if self.ENDFconversionFlag and "LRF3" in self.ENDFconversionFlag:
        return writeAsLRF3( self, flags, targetInfo, verbosityIndent = verbosityIndent )

    KRM = {'SLBW':1, 'MLBW':2, 'Reich_Moore':3, 'Full R-Matrix':4} [self.approximation]
    try:
        endf = [endfFormatsModule.endfHeadLine(0,0,self.reducedWidthAmplitudes,KRM,
            len(self.spinGroups),self.relativisticKinematics)]
    except TypeError as err:
        raise TypeError("Got '%s' when formatting '%s'"%(err.message,str((0,0,self.reducedWidthAmplitudes,KRM,
            len(self.spinGroups),self.relativisticKinematics))))

    endf.extend( writeRMatrixParticlePairs(self, targetInfo) )

    for spingrp in self.spinGroups:
        NCH, spinGroupHeader = writeRMatrixSpinGroupHeader(self, spingrp)
        endf.extend( spinGroupHeader )

        # resonances:
        NRS = len(spingrp.resonanceParameters.table)
        NX = (NCH//6 + 1)*NRS
        if NRS==0: NX=1 # special case
        endf.append( endfFormatsModule.endfHeadLine( 0,0,0,NRS,6*NX,NX ) )
        for res in spingrp.resonanceParameters.table:
            for jidx in range(NCH//6+1):
                endfLine = res[jidx*6:jidx*6+6]
                while len(endfLine)<6: endfLine.append(0)
                endf.append( endfFormatsModule.endfDataLine( endfLine ) )
        if NRS==0:
            endf.append( endfFormatsModule.endfDataLine( [0,0,0,0,0,0] ) )
    return endf

resonancesModule.RMatrix.toENDF6 = toENDF6

#
# LRF=3 is converted and stored in RMatrix. Logic below translates back to LRF=3:
#
def writeAsLRF3( RMatrix, flags, targetInfo, verbosityIndent = '' ):
    endf = []
    elastic, = [reac for reac in RMatrix.resonanceReactions if
               reac.reactionLink.link == targetInfo['reactionSuite'].getReaction('elastic')]
    AP = elastic.scatteringRadius
    if AP.isEnergyDependent():
        scatRadius = AP.form
        NR, NP = 1, len(scatRadius)
        endf.append( endfFormatsModule.endfHeadLine( 0,0,0,0, NR,NP ) )
        endf += endfFormatsModule.endfInterpolationList( (NP,
            gndToENDF6Module.gndToENDFInterpolationFlag( scatRadius.interpolation ) ) )
        endf += endfFormatsModule.endfNdDataList( scatRadius.convertAxisToUnit( 0, '10*fm' ) )
        AP = 0
    else :
        AP = AP.getValueAs( '10*fm' )

    table = {'Ls':[], 'energies':[], 'Js':[], 'elastic':[], 'capture':[], 'fission width_1':[], 'fission width_2':[]}
    APLs = {}
    for sg in RMatrix.spinGroups:
        elasticChan, = [chan for chan in sg.channels if chan.resonanceReaction == elastic.label]
        if elasticChan.scatteringRadius:
            APL = elasticChan.scatteringRadius.getValueAs( '10*fm' )
            if APL != AP: APLs[elasticChan.L] = APL
        J = float(sg.spin)
        if J!=0 and elasticChan.channelSpin < targetInfo['spin']: J *= -1

        energies = sg.resonanceParameters.table.getColumn('energy',unit='eV')
        NE = len(energies)
        table['energies'] += energies

        table['Ls'] += [elasticChan.L] * NE
        table['Js'] += [J] * NE
        table['elastic'] += sg.resonanceParameters.table.getColumn('elastic width',unit='eV')
        table['capture'] += sg.resonanceParameters.table.getColumn('capture width',unit='eV')
        for column in ('fission width_1','fission width_2'):
            vals = sg.resonanceParameters.table.getColumn(column,unit='eV')
            if not vals: vals = [0]*NE
            table[column] += vals

    L_list = sorted(set(table['Ls']))
    NLS = len(L_list)
    if APLs and 0 not in APLs:
        APLs[0] = AP
    LAD = int(RMatrix.reconstructAngular)
    NLSC = 0
    if 'LvaluesNeededForConvergence' in RMatrix.ENDFconversionFlag:
        NLSC = int( RMatrix.ENDFconversionFlag.split('LvaluesNeededForConvergence=')[1].split(',')[0] )
    endf += [endfFormatsModule.endfHeadLine(targetInfo['spin'], AP, LAD, 0, NLS, NLSC)]

    sortedTable = sorted(zip(*(
        table['Ls'],table['energies'],table['Js'],table['elastic'],table['capture'],
        table['fission width_1'],table['fission width_2']
    )))

    defaultAP = 0
    if 'explicitAPL' in RMatrix.ENDFconversionFlag: defaultAP = AP
    for L in L_list:
        APL = APLs.get(L,defaultAP)
        resonances = [sortedTable[i][1:] for i in range(len(sortedTable)) if sortedTable[i][0] == L]
        NRS = len(resonances)
        endf.append( endfFormatsModule.endfHeadLine( targetInfo['mass'], APL, L, 0, 6*NRS, NRS ) )
        for row in resonances:
            endf.append( endfFormatsModule.endfDataLine( row ) )

    targetInfo['LRF3conversion'] = {    # save useful information for covariance re-translation
        'table': table,
        'sortedTable': sortedTable,
        'AP': AP,
        'APLs': APLs
    }
    return endf


#
# unresolvedTabulatedWidths
#
def toENDF6( self, flags, targetInfo, verbosityIndent = ''):

    endf = []
    AP = self.scatteringRadius
    if AP.isEnergyDependent():
        scatRadius = AP.form
        NR, NP = 1, len(scatRadius)
        endf.append( endfFormatsModule.endfHeadLine( 0,0,0,0, NR,NP ) )
        endf += endfFormatsModule.endfInterpolationList( (NP,
            gndToENDF6Module.gndToENDFInterpolationFlag( scatRadius.interpolation ) ) )
        endf += endfFormatsModule.endfNdDataList( scatRadius.convertAxisToUnit( 0, '10*fm' ) )
        AP = 0
    else :
        AP = AP.getValueAs( '10*fm' )

    NLS = len(self.L_values)
    LFW = targetInfo['unresolved_LFW']; LRF = targetInfo['unresolved_LRF']

    def v(val): # get value back from PhysicalQuantityWithUncertainty
        if type(val)==type(None): return
        return val.getValueAs('eV')

    if LFW==0 and LRF==1:   # 'Case A' from ENDF 2010 manual pg 70
        endf.append( endfFormatsModule.endfHeadLine( targetInfo['spin'], AP,
            targetInfo['LSSF'],0,NLS,0) )
        for Lval in self.L_values:
            NJS = len(Lval.J_values)
            endf.append(endfFormatsModule.endfHeadLine( targetInfo['mass'], 0, Lval.L, 0, 6*NJS, NJS ))
            for Jval in Lval.J_values:
                # here we only have one width per J:
                ave = Jval.constantWidths
                endf.append( endfFormatsModule.endfDataLine([v(ave.levelSpacing),Jval.J,
                    Jval.neutronDOF,v(ave.neutronWidth),v(ave.captureWidth),0]) )

    elif LFW==1 and LRF==1: # 'Case B'
        energies = self.L_values[0].J_values[0].energyDependentWidths.getColumn('energy',unit='eV')
        NE = len(energies)
        endf.append( endfFormatsModule.endfHeadLine( targetInfo['spin'], AP,
            targetInfo['LSSF'], 0, NE, NLS ) )
        nlines = int(math.ceil(NE/6.0))
        for line in range(nlines):
            endf.append( endfFormatsModule.endfDataLine( energies[line*6:line*6+6] ) )
        for Lval in self.L_values:
            NJS = len(Lval.J_values)
            endf.append( endfFormatsModule.endfHeadLine( targetInfo['mass'], 0, Lval.L, 0, NJS, 0 ) )
            for Jval in Lval.J_values:
                cw = Jval.constantWidths
                endf.append( endfFormatsModule.endfHeadLine(0,0,Lval.L,Jval.fissionDOF,NE+6,0) )
                endf.append( endfFormatsModule.endfDataLine([v(cw.levelSpacing),Jval.J,Jval.neutronDOF,
                    v(cw.neutronWidth),v(cw.captureWidth),0]) )
                fissWidths = Jval.energyDependentWidths.getColumn('fissionWidthA',unit='eV')
                for line in range(nlines):
                    endf.append( endfFormatsModule.endfDataLine( fissWidths[line*6:line*6+6] ) )

    elif LRF==2:            # 'Case C', most common in ENDF
        endf.append( endfFormatsModule.endfHeadLine( targetInfo['spin'], AP,
            targetInfo['LSSF'],0,NLS,0) )
        INT = gndToENDF6Module.gndToENDFInterpolationFlag( self.interpolation )
        for Lval in self.L_values:
            NJS = len(Lval.J_values)
            endf.append( endfFormatsModule.endfHeadLine( targetInfo['mass'], 0, Lval.L, 0, NJS, 0 ))
            for Jval in Lval.J_values:
                NE = len( Jval.energyDependentWidths )
                useConstant = not NE
                if useConstant: NE = 2
                endf.append( endfFormatsModule.endfHeadLine( Jval.J, 0, INT, 0, 6*NE+6, NE ) )
                endf.append( endfFormatsModule.endfDataLine([0,0,Jval.competitiveDOF,
                    Jval.neutronDOF,Jval.gammaDOF,Jval.fissionDOF]) )
                cws = Jval.constantWidths
                if useConstant:
                    # widths are all stored in 'constantWidths' instead. get energies from parent class
                    NE = 2; useConstant=True
                    eLow,eHigh = targetInfo['regionEnergyBounds']
                    for e in (eLow,eHigh):
                        endf.append(endfFormatsModule.endfDataLine([v(e),v(cws.levelSpacing),
                            v(cws.competitiveWidth),v(cws.neutronWidth),v(cws.captureWidth),
                            v(cws.fissionWidthA)]) )

                else:
                    table = [ Jval.energyDependentWidths.getColumn('energy',unit='eV') ]
                    for attr in ('levelSpacing','competitiveWidth','neutronWidth','captureWidth',
                            'fissionWidthA'):
                        # find each attribute, in energy-dependent or constant width section
                        column = ( Jval.energyDependentWidths.getColumn( attr, unit='eV' ) or
                                [v(getattr( cws, attr ))]*NE )
                        if not any(column): column = [0]*NE
                        table.append( column )
                    for row in zip(*table):
                        endf.append( endfFormatsModule.endfDataLine( row ) )
    return endf

resonancesModule.unresolvedTabulatedWidths.toENDF6 = toENDF6
