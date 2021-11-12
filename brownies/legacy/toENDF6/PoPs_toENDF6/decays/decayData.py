# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from math import sqrt

from brownies.legacy.toENDF6 import endfFormats as endfFormatsModule
from brownies.legacy.converting.ENDFToGNDS import ENDF_ITYPE_4

from PoPs.decays import decayData as decayDataModule
from PoPs.decays import spectrum as spectrumModule

RTYPdict = {}
for key, value in ENDF_ITYPE_4.decayType.items():
    RTYPdict[value] = key

STYPdict = {}
for key, value in ENDF_ITYPE_4.STYPProduct.items():
    STYPdict[value] = key

def toENDF6( self, MT, endfMFList, flags, info, verbosityIndent = '' ) :

    RTYPs = []
    for decayMode in self.decayModes:
        if ',' in decayMode.mode:
            RTYPlist = [str(RTYPdict[key]) for key in decayMode.mode.split(',')]
            RTYPlist.insert(1,'.')
            RTYP = float( ''.join(RTYPlist) )
        else:
            RTYP = RTYPdict[ decayMode.mode ]
        RTYPs.append( RTYP )

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
    for RTYP, decayMode in zip(RTYPs, self.decayModes):
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
