# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module converts a GNDS reactionSuite into an ACE file.
"""

from xData import date as dateModule

from fudge.processing.montecarlo import fudge2dEqualProbableBinning as fudge2dEqualProbableBinningModule

from fudge.reactionData.doubleDifferentialCrossSection.thermalNeutronScatteringLaw import coherentElastic as coherentElasticModule
from fudge.reactionData.doubleDifferentialCrossSection.thermalNeutronScatteringLaw import incoherentElastic as incoherentElasticModule
from fudge.reactionData.doubleDifferentialCrossSection.thermalNeutronScatteringLaw import incoherentInelastic as incoherentInelasticModule

from fudge.productData.distributions import uncorrelated as uncorrelatedModule
from fudge.productData.distributions import angular as angularModule
from fudge.productData.distributions import energy as energyModule
from fudge.productData.distributions import energyAngularMC as energyAngularMCModule

from . import GNDS_toACE_misc as GNDS_toACE_miscModule

floatFormatOthers = '%20.12E'
floatFormatEnergyGrid = '%21.13E'
floatFormatEnergyGrid = floatFormatOthers

standardsData = {   'H1'    : [  125, 0.999167,     [  1001 ] ], 
                    'HC'    : [  125, 0.999167,     [  1001, 6012 ] ],
                    'H2'    : [  128, 1.996800,     [  1002 ] ],
                    'Be'    : [  425, 8.934780,     [  4009 ] ], 
                    'C'     : [  625, 11.893650,    [  6012, 6000 ] ], 
                    'N'     : [  725, 13.882780,    [  7014 ] ], 
                    'O'     : [  825, 15.857510,    [  8016 ] ], 
                    'SiO2'  : [ 825,  15.857510,    [  8016, 8017, 14028 ] ], 
                    'Al'    : [ 1325, 26.749750,    [ 13027 ] ], 
                    'Si'    : [ 1425, 27.737000,    [ 14028, 14029, 14030 ] ], 
                    'Fe'    : [ 2631, 55.454430,    [ 26056, 26054, 26057 ] ], 
                    'Y'     : [ 3925, 88.142100,    [ 39089 ] ], 
                    'Zr'    : [ 4025, 89.132400,    [ 40090, 40091, 40092 ] ],
                    'U'     : [ 9237, 236.005800,   [ 92238 ] ] }

targetMapping = {   'Be-metal'              : [ 'be-met',   'Be' ],
                    'BeinBeO'               : [ 'be-beo',   'Be' ],
                    'tnsl-Al27'             : [ 'al-27',    'Al' ],
                    'benzene'               : [ 'benz',     'HC' ],
                    'CinSiC'                : [ 'c-sic',    'C' ],
                     'DinD2O'               : [ 'd-d2o',    'H2' ],
                    'tnsl-Fe56'             : [ 'fe-56',    'Fe' ],
                    'crystalline-graphite'  : [ 'grph',     'C' ],
                    'reactor-graphite-10P'  : [ 'grph10',   'C' ],
                    'reactor-graphite-30P'  : [ 'grph30',   'C' ],
                    'HinH2O'                : [ 'h-h2o',    'H1' ],
                    'HinYH2'                : [ 'h-yh2',    'H1' ],
                    'HinZrH'                : [ 'h-zrh',    'H1' ],
                    'NinUN'                 : [ 'n-un',     'N' ],
                    'OinBeO'                : [ 'o-beo',    'O' ],
                    'OinD2O'                : [ 'o-d2o',    'O' ],
                    'OinUO2'                : [ 'o-uo2',    'O' ],
                    'SiinSiC'               : [ 'si-sic',   'Si' ],
                    'UinUN'                 : [ 'u-un',     'U' ],
                    'UinUO2'                : [ 'u-uo2',    'U' ],
                    'YinYH2'                : [ 'y-yh2',    'Y' ],
                    'ZrinZrH'               : [ 'zr-zrh',   'Zr' ],
                    'HinC5O2H8'             : [ 'h-luci',   'H1' ],
                    'HinCH2'                : [ 'h-poly',   'H1' ],
                    'SiO2-alpha'            : [ 'sio2',     'SiO2' ] }

def getEqualProbableBins(xs_pdf_cdf, numberOfBins, check=False):
    """Returns the equal probable bins of xs_pdf_cdf as a list."""

    pdf = xs_pdf_cdf.to_pdf_and_cdf()[0]
    epbs = pdf.equalProbableBins(numberOfBins)
    if check:
        numberOfErrors = pdf.checkEqualProbableBinsResult(epbs)
        if numberOfErrors != 0:
            pdf.checkEqualProbableBinsResult(epbs, printResults=True)
            raise Exception( 'checkEqualProbableBinsResult failed for outerDomainValue = %s' % xs_pdf_cdf.outerDomainValue)
    return epbs[1:-1]

def toACE(self, args, styleLabel, cdf_style, fileName, evaluationId, addAnnotation):

    NILm1 = args.NILm1
    NCL = args.NCL
    IFENG = 2
    checkEqualProbableBinning = False

    targetMap = targetMapping[self.target]
    standardData = standardsData[targetMap[1]]

    strRecords = []
    date = dateModule.Date.parse(cdf_style.date).date
    date = '%s/%s/%s' % (date.month, date.day, date.year % 100)
    
    line = '%6s.%2.2dt%12.6f%12.5e %10s' % ( targetMap[0], evaluationId, standardData[1], cdf_style.temperature.getValueAs('MeV/k'), date)
    strRecords.append(line)

    line = '%-70s' % ('%s at %12.2fK' % (self.target, cdf_style.temperature.getValueAs('K'))) + '%10s' % ('mat%4d' % standardData[0])
    strRecords.append(line)

    IZ = 16 * [ 0 ]
    for index, ZA in enumerate(standardData[2]): IZ[index] = ZA
    for index1 in range(4):
        line = ''
        for index2 in range(4): line += '%7d%11.0f' % (IZ[index1*2+index2], 0.0)
        strRecords.append(line)
            

    NXS = 16 * [ 0 ]        # NXS[1]: NXSS, length of the XSS block.
                            # NXS[2]: IDPNI, inelastic scattering mode. Currently, only allowed value is 3.
                            # NXS[3]: NIL, inelastic dimensioning parameter.
                            # NXS[4]: NIEB, number of inelastic exiting energies.
                            # NXS[5]: IDPNC, elastic scattering mode (coherent=4, incoherent!=4).
                            # NXS[6]: NCL, elastic dimensioning parameter.
                            # NXS[7]: IFENG, secondary energy mode (discrete=0, skewed=1, continuous=2).

    JXS = 32 * [ 0 ]        # JXS[1]: ITIE, location of inelastic energy table
                            # JXS[2]: ITIX, location of inelastic cross sections
                            # JXS[3]: ITXE, location of inelastic energy/ angle distributions
                            # JXS[4]: ITCE, location of elastic energy table
                            # JXS[5]: ITCX, location of elastic cross sections
                            # JXS[6]: ITCA, location of elastic angular distributions

    IDPNC = 0
    ITCE_X = None
    ITCA_X = []
    ITIE_X = None
    ITIEA_X = []

    for reaction in self.reactions:
        doubleDifferential = reaction.doubleDifferentialCrossSection[0]
        if isinstance(doubleDifferential, coherentElasticModule.form):          # IDPNC = 4.
            IDPNC = 4
            gridded2d = doubleDifferential.S_table.gridded2d
            temperatures = gridded2d.axes[2]
            temperature = cdf_style.temperature.getValueAs(temperatures.unit)
            deltas = [abs(value - temperature) for value in temperatures.values]
            index = deltas.index(min(deltas))
            energies = gridded2d.axes[1].values
            array = [float(value) for value in gridded2d.array.constructArray()[index]]
            ITCE_X = [energies, array]
            NXS[6-1] = -1

        elif isinstance(doubleDifferential, incoherentElasticModule.form):      # IDPNC != 4.
            crossSection = reaction.crossSection[cdf_style.label]
            ITCE_X = crossSection.copyDataToXsAndYs( )
            NXS[6-1] = NCL
            if NCL == -1:
                IDPNC = 4
                for index, energy in enumerate(ITCE_X[0]): ITCE_X[1][index] *= energy
            else:
                NXS[6-1] = NCL - 1
                IDPNC = 3
                energies = ITCE_X[0]

                distribution = reaction.outputChannel.products[0].distribution[cdf_style.label]
                if not isinstance( distribution, uncorrelatedModule.form ):
                    raise Exception('Incoherent elastic distribution is "%s" but needs to be "%s".' % (distribution.moniker, uncorrelatedModule.form.moniker))

                angular2d = distribution.angularSubform.data
                if not isinstance(angular2d, angularModule.XYs2d):
                    raise Exception('Incoherent elastic angular distribution is "%s" but needs to be "%s".' % (angular2d.moniker, angularModule.XYs2d.moniker))

                for index, angular1d in enumerate(angular2d):
                    if not isinstance(angular1d, angularModule.xs_pdf_cdf1d):
                        raise Exception('Incoherent elastic 1d angular distribution is "%s" but needs to be "%s".' % (angular1d.moniker, angularModule.xs_pdf_cdf1d.moniker))
                    if ( abs(energies[index] - angular1d.outerDomainValue) > 1e-7 * energies[index] ):
                        raise Exception('Incoherent elastic cross section and angular energies differ at %s and %s.' % (energies[index], angular1d.outerDomainValue))
                    ITCA_X += getEqualProbableBins(angular1d, 2 * NCL, checkEqualProbableBinning)[::2]

        elif isinstance(doubleDifferential, incoherentInelasticModule.form):
            crossSection = reaction.crossSection[cdf_style.label]
            ITIE_X = crossSection.copyDataToXsAndYs( )

            distribution = reaction.outputChannel.products[0].distribution[cdf_style.label]
            if not isinstance( distribution, energyAngularMCModule.form ):
                raise Exception('Incoherent inelastic distribution is "%s" but needs to be "%s".' % (distribution.moniker, energyAngularMCModule.form.moniker))
            
            energy2d = distribution.energy.data
            if not isinstance(energy2d, energyModule.XYs2d):
                raise Exception('Inelastic energy distribution is "%s" but needs to be "%s".' % (energy2d.moniker, energyModule.XYs2d.moniker))

            angular3d = distribution.energyAngular.data
            if not isinstance(angular3d, energyAngularMCModule.XYs3d):
                raise Exception('Inelastic angular distribution is "%s" but needs to be "%s".' % (angular3d.moniker, energyAngularMCModule.XYs2d.moniker))

            if len(energy2d) != len(angular3d): raise Exception('Inelastic energy and angular have different number of data: %s vs %s' % (len(energy2d), len(angular3d)))
            NXS[2-1] = 3                        # Equally-likely cosines; currently the only scattering mode allowed for inelastic angular distributions.
            NXS[3-1] = NILm1 + 1
            NXS[4-1] = 0
            NXS[7-1] = IFENG

            for index1, energy1d in enumerate(energy2d):
                if not isinstance(energy1d, energyModule.xs_pdf_cdf1d):
                    raise Exception('Incoherent inelastic outgoing energy1d is %s but needs to be "%s".' % (energy1d.moniker, energyModule.xs_pdf_cdf1d.moniker))

                angular2d = angular3d[index1]
                if not isinstance(angular2d, energyAngularMCModule.XYs2d):
                    raise Exception('Incoherent inelastic outgoing angular2d is %s but needs to be "%s".' % (angular2d.moniker, energyAngularMCModule.energyAngularMC.XYs2d.moniker))

                if(len(energy1d) != len(angular2d)):
                    raise Exception('Incoherent inelastic energy1d and angular2d lenghts differ: %s vs %s' % (len(energy1d), len(angular2d)))

                if ( abs(energy1d.outerDomainValue - angular2d.outerDomainValue) > 1e-7 * energy1d.outerDomainValue ):
                    raise Exception('Incoherent inelastic energy1d and angular2d energies differ at %s and %s.' % (energy1d.outerDomainValue, angular2d.outerDomainValue))

                xs  = energy1d.xs.values
                pdf = energy1d.pdf.values
                cdf = energy1d.cdf.values
                subITIEA_X = []
                for index2, angular1d in enumerate(angular2d):
                    if not isinstance(angular1d, energyAngularMCModule.xs_pdf_cdf1d):
                        raise Exception('Incoherent inelastic outgoing angular1d is %s but needs to be "%s".' % (angular1d.moniker, energyAngularMCModule.xs_pdf_cdf1d.moniker))

                    subITIEA_X += [xs[index2], pdf[index2], cdf[index2]]
                    subITIEA_X += getEqualProbableBins(angular1d, 2 * NILm1, checkEqualProbableBinning)[::2]

                ITIEA_X.append([len(energy1d), subITIEA_X])
        else:
            raise Exception( 'Unsupport ENDF_MT = %d' % reaction.ENDF_MT )

    NXS[5-1] = IDPNC

    XSS = []
    annotates = []

    if ITIE_X is not None:
        JXS[1-1] = len(XSS) + 1
        energyGrid, crossSection = ITIE_X
        GNDS_toACE_miscModule.updateXSSInfo('Number of inelastic energies', annotates, XSS, [ len( energyGrid ) ])
        GNDS_toACE_miscModule.updateXSSInfo('Inelastic energy', annotates, XSS, energyGrid)
        JXS[2-1] = len(XSS) + 1
        GNDS_toACE_miscModule.updateXSSInfo('Inelastic cross section', annotates, XSS, crossSection)

        numberOfIncidentEnergies = len(energyGrid)
        offset = len(XSS)
        JXS[3-1] = offset + 1
        incidentEnergyLocator = numberOfIncidentEnergies * [ 0 ]
        GNDS_toACE_miscModule.updateXSSInfo('Inelastic distribution incident energy locators', annotates, XSS, incidentEnergyLocator)
        outgoingEnergyLengths = numberOfIncidentEnergies * [ 0 ]
        GNDS_toACE_miscModule.updateXSSInfo('Inelastic distribution outgoing energy lenghts', annotates, XSS, outgoingEnergyLengths)

        for index, number_subITIEA_X in enumerate(ITIEA_X):
            incidentEnergyLocator[index] = len(XSS)
            outgoingEnergyLengths[index] = number_subITIEA_X[0]
            subITIEA_X = number_subITIEA_X[1]
            GNDS_toACE_miscModule.updateXSSInfo('Inelastic distribution data for incident energy index %s' % index, annotates, XSS, subITIEA_X)

        XSS[offset:offset+numberOfIncidentEnergies] = incidentEnergyLocator
        offset += numberOfIncidentEnergies
        XSS[offset:offset+numberOfIncidentEnergies] = outgoingEnergyLengths

    if ITCE_X is not None:
        JXS[4-1] = len(XSS) + 1
        energyGrid, crossSection = ITCE_X
        GNDS_toACE_miscModule.updateXSSInfo('Number of elastic energy', annotates, XSS, [ len( energyGrid ) ])
        GNDS_toACE_miscModule.updateXSSInfo('Elastic energy', annotates, XSS, energyGrid)
        JXS[5-1] = len(XSS) + 1
        GNDS_toACE_miscModule.updateXSSInfo('Elastic cross section', annotates, XSS, crossSection)

        if IDPNC != 4:
            JXS[6-1] = len(XSS) + 1
            GNDS_toACE_miscModule.updateXSSInfo('Elastic angular equal probably bins', annotates, XSS, ITCA_X)

    NXS[1-1] = len(XSS)

    strRecords += GNDS_toACE_miscModule.intArrayToRecords(NXS)
    strRecords += GNDS_toACE_miscModule.intArrayToRecords(JXS)
    strRecords += GNDS_toACE_miscModule.XSSToStrings(annotates, XSS, addAnnotation, len(energyGrid))

    strRecords.append('')
    fOut = open(fileName, 'w')
    fOut.write('\n'.join(strRecords))
    fOut.close()