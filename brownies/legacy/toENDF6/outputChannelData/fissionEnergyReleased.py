# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from fudge.outputChannelData.fissionFragmentData import fissionEnergyRelease as fissionEnergyReleaseModule

from .. import endfFormats as endfFormatsModule
from .. import gndsToENDF6 as gndsToENDF6Module

#
# component
#
def toENDF6( self, endfMFList, flags, targetInfo ) :

    order = 0

    energyReleaseTerms = (
        'promptProductKE',      'promptNeutronKE',      'delayedNeutronKE', 'promptGammaEnergy',
        'delayedGammaEnergy',   'delayedBetaEnergy',    'neutrinoEnergy',   'nonNeutrinoEnergy',
        'totalEnergy')

    for term in energyReleaseTerms:
        form = getattr( self, term ).data
        if isinstance( form, fissionEnergyReleaseModule.Polynomial1d ):
            order = len( form.coefficients )-1

    n, data = order+1, {}
    tabulatedTerms = []
    for i in range( n ) :
        data[i] = []
        for term in energyReleaseTerms:
            form = getattr( self, term ).data
            if isinstance( form, fissionEnergyReleaseModule.Polynomial1d ):
                data[i] += [ form.coefficients[i], form.uncertainty.data.coefficients[i] ]
            elif isinstance( form, fissionEnergyReleaseModule.XYs1d ):
                data[i] += [ form.evaluate(1e-5), 0 ]
                tabulatedTerms.append( term )
            else:
                raise NotImplementedError("Unknown fission energy release form %s" % type(form))
    datalist = []
    for i in range( n ) : datalist += data[i]

    LFC = len(tabulatedTerms) > 0
    endfMFList[1][458] = [ endfFormatsModule.endfContLine( targetInfo['ZA'], targetInfo['mass'], 0, LFC, 0, len(tabulatedTerms) ) ]
    endfMFList[1][458].append( endfFormatsModule.endfContLine( 0, 0, 0, order, 18 * n, 9 * n ) )
    endfMFList[1][458] += endfFormatsModule.endfDataList( datalist )

    for term in tabulatedTerms:
        LDRV = 2 if term == 'promptProductKE' else 1    # FIXME, hard-coded
        IFC = energyReleaseTerms.index(term) + 1
        form = getattr( self, term ).data

        endfInterpolation = gndsToENDF6Module.gndsToENDFInterpolationFlag(form.interpolation)
        flatData = []
        for xy in form.copyDataToXYs(): flatData += xy
        interpolations = [len(flatData) / 2, endfInterpolation]
        endfMFList[1][458].append(endfFormatsModule.endfContLine(0, 0, LDRV, IFC, len(interpolations) / 2,
                                                                                len(flatData) / 2))
        endfMFList[1][458] += endfFormatsModule.endfInterpolationList(interpolations)
        endfMFList[1][458] += endfFormatsModule.endfDataList(flatData)

    endfMFList[1][458].append( endfFormatsModule.endfSENDLineNumber( ) )


fissionEnergyReleaseModule.FissionEnergyRelease.toENDF6 = toENDF6
