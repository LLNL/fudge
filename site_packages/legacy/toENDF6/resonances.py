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

import math
import endfFormats as endfFormatsModule
import gndToENDF6
from pqu import PQU as PQUModule
import fudge.gnd.resonances as resonancesModule


#
# resonances
#
def toENDF6( self, endfMFList, flags, targetInfo, verbosityIndent = '' ) :

    ZAM, AWT = targetInfo['ZA'], targetInfo['mass']
    NIS, ABN = 1, 1.0; ZAI=ZAM  # assuming only one isotope per file
        
    # get target spin from the particle list:
    target = targetInfo['reactionSuite'].getParticle( targetInfo['reactionSuite'].target.name )
    targetInfo['spin'] = target.getSpin().value
        
    endf = [endfFormatsModule.endfHeadLine( ZAM, AWT, 0, 0, NIS, 0)]
    resolvedCount, unresolvedCount = 0, 0
    # resolved may have multiple energy regions:
    if self.resolved is not None: resolvedCount = max(1,len(self.resolved.regions))
    if self.unresolved is not None: unresolvedCount = 1
        
    # resonances may only contain a scattering radius:
    if not (resolvedCount + unresolvedCount) and self.scatteringRadius:
        scatRadius = self.scatteringRadius.form
        lowerBound, upperBound = scatRadius.bounds
        endf.append( endfFormatsModule.endfHeadLine( ZAM, ABN, 0,0,1,0 ) )
        endf.append( endfFormatsModule.endfHeadLine(lowerBound.getValueAs('eV'),
            upperBound.getValueAs('eV'), 0,0,0,0 ) )
        AP = scatRadius.value.getValueAs('10*fm')
        endf.append( endfFormatsModule.endfHeadLine(targetInfo['spin'], AP, 0,0,0,0 ) )
        endf.append( endfFormatsModule.endfSENDLineNumber() )
        endfMFList[2][151] = endf
        return

    # For now I'm storing the LRF/LFW flags in xml, since they are tricky to compute
    # LFW is a pain: only applies to unresolved, but must be written at the start of MF2
    LRFurr, LFW = 0,0
    if unresolvedCount != 0:
        LRF_LFW = self.unresolved.nativeData.ENDFconversionFlag
        LRFurr, LFW = map(int, LRF_LFW.split('=')[-1].split(',') )
    NER = resolvedCount + unresolvedCount
    endf.append( endfFormatsModule.endfHeadLine( ZAI, ABN, 0, LFW, NER, 0 ) )
    for idx in range(resolvedCount):
        if resolvedCount==1: region = self.resolved
        else: region = self.resolved.regions[idx]
        LRU = 1 #resolved
        LRF = { resonancesModule.SLBW.moniker:1,
                resonancesModule.MLBW.moniker:2,
                resonancesModule.RM.moniker:3,
                resonancesModule.RMatrix.moniker:7
              }[region.nativeData.moniker]
        EL, EH = region.lowerBound.getValueAs('eV'), region.upperBound.getValueAs('eV')
        if LRF==7: NRO = 0
        else: NRO = region.nativeData.scatteringRadius.isEnergyDependent()
        NAPS = not( region.nativeData.calculateChannelRadius )
        endf.append(endfFormatsModule.endfHeadLine( EL,EH,LRU,LRF,NRO,NAPS ) )
        endf += region.nativeData.toENDF6( flags, targetInfo, verbosityIndent )
    if unresolvedCount != 0:
        LRU = 2 #unresolved
        region = self.unresolved
        EL, EH = region.lowerBound.getValueAs('eV'), region.upperBound.getValueAs('eV')
        NRO, NAPS = 0,0
        endf.append(endfFormatsModule.endfHeadLine( EL,EH,LRU,LRFurr,NRO,NAPS ) )
        # pass LFW/LRF so we don't have to calculate twice:
        targetInfo['unresolved_LFW'] = LFW
        targetInfo['unresolved_LRF'] = LRFurr
        targetInfo['regionEnergyBounds'] = (region.lowerBound, region.upperBound)
        endf += region.nativeData.toENDF6( flags, targetInfo, verbosityIndent )
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
            gndToENDF6.gndToENDFInterpolationFlag( scatRadius.interpolation ) ) )
        endf += endfFormatsModule.endfNdDataList( scatRadius.convertAxisToUnit( 0, '10*fm' ) )
        AP = 0
    else :
        AP = self.scatteringRadius.getValueAs( '10*fm' )
    L_list = self.resonanceParameters.table.getColumn('L')
    NLS = len( set(L_list) )
    LAD = getattr(self, 'computeAngularDistribution') or 0
    NLSC = getattr(self, 'LvaluesNeededForConvergence') or 0
    endf += [endfFormatsModule.endfHeadLine( targetInfo['spin'], AP, LAD, 0, NLS, NLSC )]
        
    table = [ self.resonanceParameters.table.getColumn('energy',units='eV'),
            self.resonanceParameters.table.getColumn('J') ]
    NE = len(table[0])
    if isinstance( self, (resonancesModule.SLBW, resonancesModule.MLBW) ):
        attrList = ('totalWidth','neutronWidth','captureWidth','fissionWidthA')
    elif isinstance( self, resonancesModule.RM ):
        attrList = ('neutronWidth','captureWidth','fissionWidthA','fissionWidthB')
    for attr in attrList:
        column = self.resonanceParameters.table.getColumn( attr, units='eV' )
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
        if self.scatteringRadius.isLdependent() and L in self.scatteringRadius.form.lvals:
            APL = self.scatteringRadius.getValueAs('10*fm', L=L)
        resonances = [table[i] for i in range(NE) if L_list[i] == L]
        NRS = len(resonances)
        endf.append( endfFormatsModule.endfHeadLine( targetInfo['mass'], APL, L, 0, 6*NRS, NRS ) )
        for row in resonances:
            endf.append( endfFormatsModule.endfDataLine( row ) )
    return endf

resonancesModule.resonanceFormalismBaseClass.toENDF6 = toENDF6

#
# RMatrix
#
def toENDF6( self, flags, targetInfo, verbosityIndent = '' ):

    KRM = {'SLBW':1, 'MLBW':2, 'Reich_Moore':3, 'Full R-Matrix':4} [self.approximation]
    endf = [endfFormatsModule.endfHeadLine(0,0,self.reducedWidthAmplitudes,KRM,
        len(self.spinGroups),self.relativisticKinematics)]
        
    # first describe all the particle pairs (two-body output channels)
    NPP = len(self.channels)
    endf.append( endfFormatsModule.endfHeadLine(0,0,NPP,0,12*NPP,2*NPP) )

    def getENDFtuple( spin, parity ):
        # ENDF combines spin & parity UNLESS spin==0. Then it wants (0,+/-1)
        if spin.value: return (spin.value * parity.value, 0)
        else: return (spin.value, parity.value)
        
    def MZIP( p ):  # helper: extract mass, z, spin and parity from particle list
        mass = p.getMass( 'amu' ) / targetInfo['reactionSuite'].getParticle( 'n' ).getMass( 'amu' )
        Z = p.getZ_A_SuffixAndZA()[0]
        I,P = getENDFtuple( p.getSpin(), p.getParity() )
        return mass, Z, I, P
        
    for pp in self.channels:
        pA,pB = pp.name.split(' + ')
        # get the xParticle instances for pA and pB:
        pA,pB = targetInfo['reactionSuite'].getParticle(pA), targetInfo['reactionSuite'].getParticle(pB)
        MA, ZA, IA, PA = MZIP( pA )
        MB, ZB, IB, PB = MZIP( pB )
        MT = pp.ENDF_MT
        PNT = pp.calculatePenetrability
        if PNT is None: PNT = self.calculatePenetrability
        SHF = pp.calculateShift
        if SHF is None: SHF = self.calculateShift
        if MT in (19,102): PNT = 0  # special case
        try: Q = pp.Q.inUnitsOf('eV').value
        except:
            Q = targetInfo['reactionSuite'].getReaction( pp.channel ).getQ('eV')
            # getQ doesn't account for residual left in excited state:
            for particle in targetInfo['reactionSuite'].getReaction( pp.channel ).outputChannel.particles:
                if( hasattr( particle, 'getLevelAsFloat' ) ) : Q -= particle.getLevelAsFloat('eV')
        if MT==102: Q = 0
        endf.append( endfFormatsModule.endfDataLine( [MA,MB,ZA,ZB,IA,IB] ) )
        endf.append( endfFormatsModule.endfDataLine( [Q,PNT,SHF,MT,PA,PB] ) )
        
    for spingrp in self.spinGroups:
        AJ,PJ = getENDFtuple( spingrp.spin, spingrp.parity )
        KBK = spingrp.background
        KPS = spingrp.applyPhaseShift
        NCH = len(spingrp.resonanceParameters.table.columns)-1    # skip the 'energy' column
        endf.append( endfFormatsModule.endfHeadLine( AJ,PJ,KBK,KPS,6*NCH,NCH ) )
        for chan in spingrp.resonanceParameters.table.columns:
            # which open channel does it correspond to?
            if chan.name == 'energy': continue
            name = chan.name.split(' width')[0]
            PPI, openChannel = [ch for ch in enumerate(self.channels) if ch[-1].name == name][0]
            PPI += 1  # 1-based index in ENDF
            attr = chan.attributes
            L = attr['L']
            SCH = attr['channelSpin']
            # some data may have been moved up to the channel list:
            BND = attr.get('boundaryCondition') or self.boundaryCondition

            channelOverride = resonancesModule.channelOverride( openChannel.label )
            if openChannel.label in spingrp.overrides:
                channelOverride = spingrp.overrides[ openChannel.label ]

            APT = channelOverride.scatteringRadius or openChannel.scatteringRadius or PQUModule.PQU(0,'fm')
            APE = channelOverride.effectiveRadius or openChannel.effectiveRadius or APT
            APT = APT.getValueAs('10*fm')
            APE = APE.getValueAs('10*fm')
            if openChannel.ENDF_MT==102:
                APT, APE = 0,0
            endf.append( endfFormatsModule.endfDataLine( [PPI,L,SCH,BND,APE,APT] ) )
        # resonances:
        NRS = len(spingrp.resonanceParameters.table)
        NX = (NCH//6 + 1)*NRS
        if NRS==0: NX=1 # special case
        endf.append( endfFormatsModule.endfHeadLine( 0,0,0,NRS,6*NX,NX ) )
        for res in spingrp.resonanceParameters.table:
            #resproperties = [res.energy.inUnitsOf('eV').value] + [
            #        w.inUnitsOf('eV').value for w in res.widths]
            for jidx in range(NCH//6+1):
                endfLine = res[jidx*6:jidx*6+6]
                while len(endfLine)<6: endfLine.append(0)
                endf.append( endfFormatsModule.endfDataLine( endfLine ) )
        if NRS==0:
            endf.append( endfFormatsModule.endfDataLine( [0,0,0,0,0,0] ) )
    return endf

resonancesModule.RMatrix.toENDF6 = toENDF6

#
# unresolvedTabulatedWidths
#
def toENDF6( self, flags, targetInfo, verbosityIndent = ''):

    AP = self.scatteringRadius.form.value.inUnitsOf('10*fm').value
    NLS = len(self.L_values)
    LFW = targetInfo['unresolved_LFW']; LRF = targetInfo['unresolved_LRF']
        
    def v(val): # get value back from PhysicalQuantityWithUncertainty
        if type(val)==type(None): return
        return val.getValueAs('eV')
        
    if LFW==0 and LRF==1:   # 'Case A' from ENDF 2010 manual pg 70
        endf = [endfFormatsModule.endfHeadLine( targetInfo['spin'], AP, 
            self.forSelfShieldingOnly,0,NLS,0) ]
        for Lval in self.L_values:
            NJS = len(Lval.J_values)
            endf.append(endfFormatsModule.endfHeadLine( targetInfo['mass'], 0, Lval.L, 0, 6*NJS, NJS ))
            for Jval in Lval.J_values:
                # here we only have one width per J:
                ave = Jval.constantWidths
                endf.append( endfFormatsModule.endfDataLine([v(ave.levelSpacing),Jval.J.value,
                    Jval.neutronDOF,v(ave.neutronWidth),v(ave.captureWidth),0]) )
        
    elif LFW==1 and LRF==1: # 'Case B'
        energies = self.L_values[0].J_values[0].energyDependentWidths.getColumn('energy',units='eV')
        NE = len(energies)
        endf = [endfFormatsModule.endfHeadLine( targetInfo['spin'], AP,
            self.forSelfShieldingOnly, 0, NE, NLS )]
        nlines = int(math.ceil(NE/6.0))
        for line in range(nlines):
            endf.append( endfFormatsModule.endfDataLine( energies[line*6:line*6+6] ) )
        for Lval in self.L_values:
            NJS = len(Lval.J_values)
            endf.append( endfFormatsModule.endfHeadLine( targetInfo['mass'], 0, Lval.L, 0, NJS, 0 ) )
            for Jval in Lval.J_values:
                cw = Jval.constantWidths
                endf.append( endfFormatsModule.endfHeadLine(0,0,Lval.L,Jval.fissionDOF,NE+6,0) )
                endf.append( endfFormatsModule.endfDataLine([v(cw.levelSpacing),Jval.J.value,Jval.neutronDOF,
                    v(cw.neutronWidth),v(cw.captureWidth),0]) )
                fissWidths = Jval.energyDependentWidths.getColumn('fissionWidthA',units='eV')
                for line in range(nlines):
                    endf.append( endfFormatsModule.endfDataLine( fissWidths[line*6:line*6+6] ) )
        
    elif LRF==2:            # 'Case C', most common in ENDF
        endf = [endfFormatsModule.endfHeadLine( targetInfo['spin'], AP, 
            self.forSelfShieldingOnly,0,NLS,0) ]
        INT = gndToENDF6.gndToENDFInterpolationFlag( self.interpolation )
        for Lval in self.L_values:
            NJS = len(Lval.J_values)
            endf.append( endfFormatsModule.endfHeadLine( targetInfo['mass'], 0, Lval.L, 0, NJS, 0 ))
            for Jval in Lval.J_values:
                NE = len( Jval.energyDependentWidths )
                useConstant = not NE
                if useConstant: NE = 2
                endf.append( endfFormatsModule.endfHeadLine( Jval.J.value, 0, INT, 0, 6*NE+6, NE ) )
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
                    table = [ Jval.energyDependentWidths.getColumn('energy',units='eV') ]
                    for attr in ('levelSpacing','competitiveWidth','neutronWidth','captureWidth',
                            'fissionWidthA'):
                        # find each attribute, in energy-dependent or constant width section
                        column = ( Jval.energyDependentWidths.getColumn( attr, units='eV' ) or
                                [v(getattr( cws, attr ))]*NE )
                        if not any(column): column = [0]*NE
                        table.append( column )
                    for row in zip(*table):
                        endf.append( endfFormatsModule.endfDataLine( row ) )
    return endf

resonancesModule.unresolvedTabulatedWidths.toENDF6 = toENDF6
