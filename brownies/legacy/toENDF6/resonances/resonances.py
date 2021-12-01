# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from pqu import PQU as PQUModule

from fudge.resonances import resonances as resonancesModule, common as commonResonancesModule, \
    resolved as resolvedModule

from brownies.legacy.toENDF6 import endfFormats as endfFormatsModule

def toENDF6( self, endfMFList, flags, targetInfo, verbosityIndent = '' ) :
    """
    Write gnds.resonances.resonances to ENDF6

    :param self:
    :param endfMFList:
    :param flags:
    :param targetInfo:
    :param verbosityIndent:
    :return:
    """

    ZAM, AWT = targetInfo['ZA'], targetInfo['mass']
    NIS, ABN = 1, 1.0; ZAI=ZAM  # assuming only one isotope per file

    # get target spin from the particle list:
    reactionSuite = targetInfo['reactionSuite']
    targetID = reactionSuite.target
    if( targetID in reactionSuite.PoPs.aliases ) : targetID = reactionSuite.PoPs[targetID].pid
    targetInfo['targetID'] = targetID
    target = reactionSuite.PoPs[targetID]
    if hasattr(target, 'nucleus'):
        targetInfo['spin'] = target.nucleus.spin[0].value
    else:
        targetInfo['spin'] = target.spin[0].value

    endf = [endfFormatsModule.endfHeadLine( ZAM, AWT, 0, 0, NIS, 0)]
    resolvedCount, unresolvedCount = 0, 0
    # resolved may have multiple energy regions:
    if self.resolved is not None:
        resolvedCount = 1
        if isinstance(self.resolved.evaluated, commonResonancesModule.energyIntervals):
            resolvedCount = len(self.resolved.evaluated)
    if self.unresolved is not None: unresolvedCount = 1

    # resonances may only contain a scattering radius:
    if not (resolvedCount + unresolvedCount) and self.scatteringRadius:
        scatRadius = self.scatteringRadius.form
        lowerBound, upperBound = scatRadius.domainMin, scatRadius.domainMax
        endf.append( endfFormatsModule.endfHeadLine( ZAM, ABN, 0,0,1,0 ) )
        endf.append( endfFormatsModule.endfHeadLine(lowerBound, upperBound, 0,0,0,0 ) )
        AP = PQUModule.PQU( self.scatteringRadius.form.value, self.scatteringRadius.form.rangeUnit ).getValueAs('10*fm')
        endf.append( endfFormatsModule.endfHeadLine(targetInfo['spin'], AP, 0,0,0,0 ) )
        endf.append( endfFormatsModule.endfSENDLineNumber() )
        endfMFList[2][151] = endf
        return

    # For now I'm storing the LRF/LFW flags in xml, since they are tricky to compute
    # LFW is a pain: only applies to unresolved, but must be written at the start of MF2
    LRFurr, LFW = 0,0
    if unresolvedCount != 0:
        LRF_LFW = targetInfo['ENDFconversionFlags'].get(self.unresolved.evaluated)
        if LRF_LFW is not None:
            LRFurr, LFW = LRF_LFW.split(',')
            LRFurr = int(LRFurr.replace('LRF',''))
            LFW = int(LFW.replace('LFW',''))
    NER = resolvedCount + unresolvedCount
    endf.append( endfFormatsModule.endfHeadLine( ZAI, ABN, 0, LFW, NER, 0 ) )
    for idx in range(resolvedCount):
        if resolvedCount==1: region = self.resolved
        else: region = self.resolved.evaluated[idx]
        LRU = 1 #resolved
        if isinstance(region.evaluated, resolvedModule.BreitWigner):
            LRF = {
                resolvedModule.BreitWigner.singleLevel: 1,
                resolvedModule.BreitWigner.multiLevel: 2
            }[ region.evaluated.approximation ]
        elif isinstance(region.evaluated, resolvedModule.RMatrix):
            LRF = 7
        EL, EH = region.domainMin, region.domainMax
        if LRF==7:
            NRO = 0
            if 'LRF3' in targetInfo['ENDFconversionFlags'].get(region.evaluated,''):
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
        targetInfo['LSSF'] = self.unresolved.evaluated.useForSelfShieldingOnly
        targetInfo['regionEnergyBounds'] = (region.domainMin, region.domainMax)
        endf += region.evaluated.toENDF6( flags, targetInfo, verbosityIndent )
    endf.append( endfFormatsModule.endfSENDLineNumber() )
    endfMFList[2][151] = endf

resonancesModule.resonances.toENDF6 = toENDF6
