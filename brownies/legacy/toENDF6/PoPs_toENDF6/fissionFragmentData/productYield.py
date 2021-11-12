# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from brownies.legacy.toENDF6 import endfFormats as endfFormatsModule

from PoPs.fissionFragmentData import productYield as productYieldModule

from brownies.legacy.converting import massTracker as massTrackerModule
from brownies.legacy.endl import misc as miscENDLModule


def toENDF6( self, endfMFList, flags, info, verbosityIndent = '' ) :

    header = [endfFormatsModule.endfHeadLine( info['ZA'], info['AWR'], 1, 0, 0, 0 )]
    productInfo = list(map(miscENDLModule.getZ_A_suffix_andZAFromName, self.nuclides.data))
    productZAs = [v[-1] for v in productInfo]
    def computeFPS(val):
        if not val: return 0
        return int(val.replace('m',''))
    productFPs = [computeFPS(v[-2]) for v in productInfo]

    # spontaneous yields are broken into two sections: prompt yields go in MT=454, delayed go in MT=459
    for duration in self.durations:
        yieldList = list(duration.yields.values)
        uncertaintyList = list(duration.yields.uncertainty.form().matrix().constructArray().diagonal())

        datalist = header + [endfFormatsModule.endfHeadLine(0,0,0,0, 4*len(productZAs), len(productZAs))]

        reorderedList = [y for x in zip(productZAs, productFPs, yieldList, uncertaintyList) for y in x]
        datalist += endfFormatsModule.endfDataList(reorderedList)
        datalist.append(endfFormatsModule.endfSENDLineNumber())
        
        if duration.time[0].value == 0:
            MT = 454
        else:
            MT = 459
        endfMFList[8][MT] = datalist

productYieldModule.productYield.toENDF6 = toENDF6
