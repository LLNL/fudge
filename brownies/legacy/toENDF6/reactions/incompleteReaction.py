# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from fudge.reactions import incompleteReaction as incompleteReactionModule

from .. import endfFormats as endfFormatsModule
from .. import gndsToENDF6 as gndsToENDF6Module

from . import base

#
# incompleteReaction
#
def toENDF6( self, endfMFList, flags, targetInfo, verbosityIndent = '' ) :

    MT = self.ENDF_MT

    if MT == 18:  # sub-actinide fission, write back to MF 8/10
        if flags['verbosity'] >= 10:
            print('%sincompleteReaction: %s' % (verbosityIndent, self.outputChannel.toString(simpleString=True)))
        ZA, mass = targetInfo['ZA'], targetInfo['mass']
        Qval = self.outputChannel.Q.evaluated.value
        ZAP, LMF = -1, 10
        endfMFList[8][MT] = [endfFormatsModule.endfHeadLine(ZA, mass, 0, 0, 1, 1),
                             endfFormatsModule.endfHeadLine(ZAP, Qval, LMF, 0, 0, 0)]
        endfMFList[8][MT].append(endfFormatsModule.endfSENDLineNumber())

        # MF10 portion:
        interpolationFlatData, flatData = gndsToENDF6Module.getForm( targetInfo['style'], self.crossSection ).toENDF6Data( MT, endfMFList, targetInfo, 0 )
        endfMFList[LMF][MT] = [endfFormatsModule.endfHeadLine(ZA, mass, 0, 0, 1, 0),
            endfFormatsModule.endfHeadLine(Qval, 0, ZAP, 0, len(interpolationFlatData) / 2, len(flatData) / 2) ]
        endfMFList[LMF][MT] += endfFormatsModule.endfInterpolationList(interpolationFlatData)
        endfMFList[LMF][MT] += endfFormatsModule.endfDataList(flatData)

        endfMFList[LMF][MT].append(endfFormatsModule.endfSENDLineNumber())

    elif MT in list(range(201,208)) + [525]:
        base.toENDF6( self, endfMFList, flags, targetInfo, verbosityIndent )

    else:
        if MT == 102:
            return

        raise NotImplementedError("incompleteReaction only supports fission right now")

incompleteReactionModule.IncompleteReaction.toENDF6 = toENDF6
