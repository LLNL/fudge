# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from brownies.legacy.toENDF6 import endfFormats as endfFormatsModule
from PoPs.atomic import atomic as atomicModule

from brownies.legacy.converting.ENDFToGNDS.ENDF_ITYPE_3_6_Misc import AtomicConfigurations_MT

def toENDF6( self, MT, endfMFList, flags, info, verbosityIndent = '' ) :

    NSS = len(self.configurations)
    endfMFList[28][MT] = [endfFormatsModule.endfHeadLine( info['ZA'], info['AWR'], 0, 0, NSS, 0 )]

    for configuration in self.configurations:
        SUBI = AtomicConfigurations_MT[configuration.subshell] - 533
        bindingEnergy = configuration.bindingEnergy.float('eV')

        data = [endfFormatsModule.endfContLine(bindingEnergy, configuration.electronNumber, 0, 0, 0, 0)]
        for decayMode in configuration.decayData.decayModes:
            finalStates = decayMode.decayPath[0].products[-1].pid.split('{')[-1].split('}')[0].split(',')  # FIXME make standard qualifiers parser
            finalSUBIs = []
            bindingEnergyDiff = bindingEnergy
            for finalState in finalStates:
                finalSUBIs.append(AtomicConfigurations_MT[finalState] - 533)
                bindingEnergyDiff -= self.configurations[finalState].bindingEnergy.float('eV')

            while len(finalSUBIs) < 2: finalSUBIs.append(0)
            data.append( endfFormatsModule.endfDataLine( finalSUBIs + [ bindingEnergyDiff, decayMode.probability.float( '' ), 0, 0 ] ) )

        endfMFList[28][MT] += [endfFormatsModule.endfContLine(SUBI, 0, 0, 0, 6 * len(data), 0)] + data

    endfMFList[28][MT].append( endfFormatsModule.endfSENDLineNumber() )

atomicModule.atomic.toENDF6 = toENDF6
