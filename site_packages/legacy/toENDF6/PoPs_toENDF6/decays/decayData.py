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

from math import sqrt

from site_packages.legacy.toENDF6 import endfFormats as endfFormatsModule
from fudge.legacy.converting.ENDFToGNDS import ENDF_ITYPE_4

from PoPs.decays import decayData as decayDataModule
from PoPs.decays import spectrum as spectrumModule

RTYPdict = {}
for key, value in ENDF_ITYPE_4.decayType.items():
    RTYPdict[value] = key

STYPdict = {}
for key, value in ENDF_ITYPE_4.STYPProduct.items():
    STYPdict[value] = key

def toENDF6( self, MT, endfMFList, flags, info, verbosityIndent = '' ) :

    for decayMode in self.decayModes:
        if ',' in decayMode.mode:
            RTYPlist = [str(RTYPdict[key]) for key in decayMode.mode.split(',')]
            RTYPlist.insert(1,'.')
            RTYP = float( ''.join(RTYPlist) )
        else:
            RTYP = RTYPdict[ decayMode.mode ]

        RFS = 0     # isomer decay product?
        for decay in decayMode.decayPath:
            for product in decay.products:
                if product.pid in info['PoPs'].aliases: RFS = 1

        Q = decayMode.Q.float('eV')
        Q_uncert = 0
        if decayMode.Q[0].uncertainty is not None:
            Q_uncert = decayMode.Q[0].uncertainty.form.value.float('eV')

        BR = decayMode.probability.float('')
        BR_uncert = 0
        if decayMode.probability[0].uncertainty is not None:
            BR_uncert = decayMode.probability[0].uncertainty.form.value.float('')

        endfMFList[8][457].append(endfFormatsModule.endfDataLine([RTYP,RFS,Q,Q_uncert,BR,BR_uncert]))

    # 2nd loop to get all the outgoing spectra:
    for decayMode in self.decayModes:
        for spectrum in decayMode.spectra:
            discretes, continuums = [],[]
            for emission in spectrum:
                if isinstance(emission, spectrumModule.discrete): discretes.append(emission)
                elif isinstance(emission, spectrumModule.continuum): continuums.append(emission)
                else:
                    raise NotImplementedError("Unknown emission mode %s" % emission.moniker)

            STYP = STYPdict[ spectrum.label ]
            if len(continuums) > 0 and len(discretes) > 0: LCON = 2
            elif len(continuums) > 0: LCON = 1
            else: LCON = 0
            LCOV = 0    # FIXME hard-coded for now, Cf258 SF probably has non-zero
            NER = len(spectrum)

            endfMFList[8][457].append(endfFormatsModule.endfContLine(0,STYP,LCON,LCOV,6,NER))
            spectrumData = []
            average_energy, d_average_energy = 0,0

            for discrete in discretes:
                ER = discrete.energy.value
                dER = 0
                if discrete.energy.uncertainty is not None:
                    dER = discrete.energy.uncertainty.form.value.float('eV')

                TYPE = 0
                if discrete.type is not None:
                    TYPE = spectrumModule.transitionType.types.index( discrete.type ) + 1

                RI = discrete.intensity.value
                dRI = 0
                if discrete.intensity.uncertainty is not None:
                    dRI = discrete.intensity.uncertainty.form.value.float('')

                average_energy += ER * RI
                d_average_energy += (ER*RI)**2 * ((dER/ER)**2 + (dRI/RI)**2)

                RIS, dRIS = 0,0
                if STYP == 2:
                    if discrete.positronEmissionIntensity is not None:
                        RIS = discrete.positronEmissionIntensity.value
                        if discrete.positronEmissionIntensity.uncertainty is not None:
                            dRIS = discrete.positronEmissionIntensity.uncertainty.form.value.float('')
                else:
                    if discrete.internalPairFormationCoefficient is not None:
                        RIS = discrete.internalPairFormationCoefficient.value
                        if discrete.internalPairFormationCoefficient.uncertainty is not None:
                            dRIS = discrete.internalPairFormationCoefficient.uncertainty.form.value.float('')

                dataList = [RTYP, TYPE, RI, dRI, RIS, dRIS]

                if len(discrete.internalConversionCoefficients) > 0:
                    for shell in discrete.internalConversionCoefficients:
                        coeff = shell.value
                        dcoeff = 0
                        if shell.uncertainty is not None:
                            dcoeff = shell.uncertainty.form.value.float('')
                        dataList += [coeff, dcoeff]

                NT = len(dataList)
                spectrumData.append( endfFormatsModule.endfContLine(ER,dER,0,0,NT,0) )
                spectrumData += endfFormatsModule.endfDataList(dataList)

            for continuum in continuums:
                pass

            endfMFList[8][457].append( endfFormatsModule.endfDataLine([1,0,average_energy,sqrt(d_average_energy),0,0]) )
            endfMFList[8][457] += spectrumData

    endfMFList[8][457].append( endfFormatsModule.endfSENDLineNumber() )


decayDataModule.decayData.toENDF6 = toENDF6
