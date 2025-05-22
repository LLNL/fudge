# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import textwrap
import time

from PoPs import IDs as IDsPoPsModule
from PoPs.quantities import halflife as halflifePoPsModule
from PoPs.families import gaugeBoson as gaugeBosonPoPsModule
from PoPs.families import lepton as leptonPoPsModule
from PoPs.families import baryon as baryonPoPsModule
from PoPs.families import nuclide as nuclidePoPsModule
from PoPs.families import nucleus as nucleusPoPsModule
from PoPs.chemicalElements import misc as chemicalElementMiscPoPsModule

from brownies.legacy.converting import endf_endl as endf_endlModule, massTracker as massTrackerModule

from fudge import documentation as documentationModule
from fudge import reactionSuite as reactionSuiteModule
from fudge import styles as stylesModule

from . import endfFormats as endfFormatsModule
from . import gndsToENDF6 as gndsToENDF6Module
from . import ENDFconversionFlags as ENDFconversionFlagsModule
from .productData import multiplicity as multiplicityModule


def toENDF6(self, styleLabel, flags, verbosityIndent='', covarianceSuite=None, useRedsFloatFormat=False,
            lineNumbers=True, **kwargs):
    if self.domainUnit != 'eV':
        self = self.parseXMLString(self.toXML(), lazyParsing=True)
        self.convertUnits({self.domainUnit: 'eV'})

    _useRedsFloatFormat = endfFormatsModule.useRedsFloatFormat
    endfFormatsModule.useRedsFloatFormat = useRedsFloatFormat

    style = self.styles[styleLabel]
    if not isinstance(style, (stylesModule.Evaluated, stylesModule.CrossSectionReconstructed, stylesModule.Heated)):
        raise TypeError('Invalid style via label "%s".' % styleLabel)
    evaluatedStyle = style
    while True:
        if evaluatedStyle is None: raise Exception('No evaluated style found via label "%s".' % styleLabel)
        if isinstance(evaluatedStyle, stylesModule.Evaluated): break
        evaluatedStyle = evaluatedStyle.derivedFromStyle

    if flags == {}: flags['verbosity'] = 0
    if flags['verbosity'] >= 10:
        print('%s%s' % (verbosityIndent, self.inputParticlesToReactionString(suffix=" -->")))
    verbosityIndent2 = verbosityIndent + ' ' * (len(self.inputParticlesToReactionString(suffix=" -->")) + 1)

    projectileZA = chemicalElementMiscPoPsModule.ZA(self.PoPs[self.projectile])
    IPART = projectileZA
    if self.projectile == 'e-': IPART = 11

    targetInfo = {}
    targetID = self.target
    if targetID in self.PoPs.aliases: targetID = self.PoPs[targetID].pid

    targetInfo['ENDFconversionFlags'] = {}
    if 'LLNL' in self.applicationData:
        LLNL_institution = self.applicationData['LLNL']
        ENDF_conversionFlags = [data for data in LLNL_institution.data if
                                data.moniker == ENDFconversionFlagsModule.ENDFconversionFlags.moniker]
        if len(ENDF_conversionFlags) > 0:
            XML_string = '\n'.join(ENDF_conversionFlags[0].toXML_strList())
            conversionFlags = ENDFconversionFlagsModule.ENDFconversionFlags.parseXMLString(XML_string)
            for link in conversionFlags.flags:
                if covarianceSuite and '/covarianceSuite/' in link.path:
                    linkTarget = link.follow(covarianceSuite)
                else:
                    linkTarget = link.follow(self)
                targetInfo['ENDFconversionFlags'][linkTarget] = link.flags

    if self.isThermalNeutronScatteringLaw():
        EMAX = self.styles.getEvaluatedStyle().projectileEnergyDomain.max
        MAT, ZA = targetInfo['ENDFconversionFlags'][self].split(',')
        MAT = int(MAT.replace('MAT=', ''))
        targetZA = int(ZA.replace('ZA=', ''))
        target = None
        targetInfo['ZA'] = MAT + 100
    else:
        EMAX = max([reaction.crossSection.domainMax for reaction in self.reactions])
        if 'MAT=' in targetInfo['ENDFconversionFlags'].get(self, ""):
            # use MAT number from conversion flags instead of recomputing
            MAT = int(targetInfo['ENDFconversionFlags'][self].replace('MAT=', ''))
            targetZA = chemicalElementMiscPoPsModule.ZAInfo_fromString(targetID)[2]
        else:
            targetZA, MAT = endf_endlModule.ZAAndMATFromParticleName(targetID)
        targetInfo['ZA'] = targetZA
        try:
            target = self.PoPs[targetID]
        except:
            target = self.PoPs.chemicalElements.getSymbol(targetID)

    targetZ, targetA = divmod(targetZA, 1000)

    targetInfo['massTracker'] = massTrackerModule.MassTracker()

    try:
        neutronMass = self.PoPs[IDsPoPsModule.neutron].getMass('amu')
    except:
        neutronMass = 1.00866491578  # From ENDF102 manual, Appendix H.4
    targetInfo['massTracker'].addMassAMU(1, neutronMass)

    for particle in self.PoPs:
        if particle.id == IDsPoPsModule.neutron:  # Already added above.
            continue
        if isinstance(particle, (gaugeBosonPoPsModule.Particle, leptonPoPsModule.Particle, baryonPoPsModule.Particle,
                                 nuclidePoPsModule.Particle)):
            ZA = chemicalElementMiscPoPsModule.ZA(particle)
            if isinstance(particle, nuclidePoPsModule.Particle):
                if particle.index != 0: continue
            if len(particle.mass) > 0: targetInfo['massTracker'].addMassAMU(ZA, particle.getMass('amu'))

    for chemicalElement in self.PoPs.chemicalElements:
        ZA = 1000 * chemicalElement.Z
        try:
            targetInfo['massTracker'].getMassAMU(ZA)
        except:  # If not present, i.e., a raise, add.
            targetInfo['massTracker'].addMassAMU(ZA, targetInfo['massTracker'].getElementalMassAMU(ZA))

    targetInfo['style'] = styleLabel
    targetInfo['reactionSuite'] = self

    levelIndex = 0
    levelEnergy_eV = 0
    STA = 0

    if isinstance(target, nuclidePoPsModule.Particle):  # isomer target
        levelIndex = target.index
        levelEnergy_eV = target.energy[0].float('eV')
        if len(target.nucleus.halflife) > 0:
            if target.nucleus.halflife[0].value == halflifePoPsModule.UNSTABLE: STA = 1
    if levelIndex > 0: STA = 1

    targetInfo['mass'] = targetInfo['massTracker'].getMassAWR(targetZA, levelEnergyInEv=levelEnergy_eV)
    if self.isThermalNeutronScatteringLaw():
        targetInfo['mass'] = self.PoPs[self.target].mass.float('amu') / targetInfo['massTracker'].neutronMass

    targetInfo['LIS'] = levelIndex
    targetInfo['metastables'] = {}
    targetInfo['LISO'] = 0
    if levelIndex > 0:
        targetInfo['LISO'] = self.PoPs[self.target].metaStableIndex
    # BRBBRB
    for alias in self.PoPs.aliases:
        if hasattr(alias, 'metaStableIndex'):
            targetInfo['metastables'][alias.pid] = alias
    MAT += targetInfo['LISO']
    if self.MAT is not None: MAT = self.MAT

    ITYPE = 0
    for reaction in self.reactions:
        if 500 <= reaction.ENDF_MT < 573: ITYPE = 3
    targetInfo['crossSectionMF'] = {0: 3, 3: 23}[ITYPE]

    if ITYPE == 3:
        targetInfo['EFL'] = 0

    targetInfo['delayedRates'] = []
    targetInfo['MTs'], targetInfo['MF8'], targetInfo['LRs'] = {}, {}, {}
    endfMFList = {1: {451: []}, 2: {}, 3: {}, 4: {}, 5: {}, 6: {}, 8: {}, 9: {}, 10: {}, 12: {}, 13: {},
                  14: {}, 15: {}, 23: {}, 26: {}, 27: {}, 31: {}, 32: {}, 33: {}, 34: {}, 35: {}, 40: {}}
    if self.resonances is not None:  # Add resonances, independent of reaction channels
        self.resonances.toENDF6(endfMFList, flags, targetInfo, verbosityIndent=verbosityIndent2)

    targetInfo['production_gammas'] = {}

    for multiplicitySum in self.sums.multiplicitySums:
        if multiplicitySum.ENDF_MT == 452:
            targetInfo['totalNubar'] = multiplicitySum.multiplicity.evaluated
        elif multiplicitySum.ENDF_MT == 455:
            targetInfo['totalDelayedNubar'] = multiplicitySum.multiplicity.evaluated

    for reaction in self:
        reaction.toENDF6(endfMFList, flags, targetInfo, verbosityIndent=verbosityIndent2)
    gndsToENDF6Module.upDateENDFMF8Data(endfMFList, targetInfo)
    for MT, production_gammas in targetInfo['production_gammas'].items():
        MF, production_gammas = production_gammas[0], production_gammas[1:]
        for productionReaction in production_gammas:
            gammas = [gamma for gamma in productionReaction.outputChannel]
            targetInfo['crossSection'] = gndsToENDF6Module.getForm(targetInfo['style'], productionReaction.crossSection)
            gndsToENDF6Module.gammasToENDF6_MF12_13(MT, MF, endfMFList, flags, targetInfo, gammas)

    for particle in self.PoPs:
        if isinstance(particle, (nucleusPoPsModule.Particle, nuclidePoPsModule.Particle)):
            if len(particle.decayData.decayModes) > 0:
                doGammaDecay = True
                for baseMT in [50, 600, 650, 700, 750, 800, 900,
                               1000]:  # 1000 causes raise in endf_endlModule.ENDF_MTZAEquation.
                    if baseMT == 1000:
                        doGammaDecay = False
                        break
                    residualZA = endf_endlModule.ENDF_MTZAEquation(projectileZA, targetZA, baseMT)[0][-1]
                    if chemicalElementMiscPoPsModule.ZA(particle) == residualZA: break
                if doGammaDecay:
                    addDecayGamma(self, particle, baseMT, endfMFList, flags, targetInfo)
                else:
                    print('WARNING: Ignoring photon branching data for particle "%s".' % particle.id)

    if 'totalNubar' in targetInfo:
        multiplicityModule.fissionNeutronsToENDF6(452, targetInfo['totalNubar'], endfMFList, flags, targetInfo)
    if 'promptNubar' in targetInfo:
        multiplicityModule.fissionNeutronsToENDF6(456, targetInfo['promptNubar'], endfMFList, flags, targetInfo)
    if 'totalDelayedNubar' in targetInfo:
        endfMFList[1][455] = [endfFormatsModule.endfHeadLine(targetZA, targetInfo['mass'], 0, 2, 0,
                                                             0)]  # Currently, only LDG = 0, LNU = 2 is supported.
        endfMFList[1][455] += [endfFormatsModule.endfHeadLine(0, 0, 0, 0, len(targetInfo['delayedRates']), 0)]
        endfMFList[1][455] += endfFormatsModule.endfDataList(targetInfo['delayedRates'])

        multiplicityModule.fissionNeutronsToENDF6(455, targetInfo['totalDelayedNubar'], endfMFList, flags, targetInfo)

        if 455 in endfMFList[5]:
            MF5MT455s = endfMFList[5][455]

            MF5MT455List = [endfFormatsModule.endfHeadLine(targetZA, targetInfo['mass'], 0, 0, len(MF5MT455s), 0)]
            for MF5MT455 in MF5MT455s: MF5MT455List += MF5MT455
            if len(MF5MT455s) == 0:
                del endfMFList[5][455]
            else:
                endfMFList[5][455] = MF5MT455List + [endfFormatsModule.endfSENDLineNumber()]

    if covarianceSuite:
        if covarianceSuite.domainUnit != 'eV':
            covarianceSuite = covarianceSuite.parseXMLString(covarianceSuite.toXML(), lazyParsing=True)
            covarianceSuite.convertUnits({covarianceSuite.domainUnit: 'eV'})
        if styleLabel not in covarianceSuite.styles: targetInfo['style'] = covarianceSuite.styles[
            0].label  # FIXME, this is a kludge.
        covarianceSuite.toENDF6(endfMFList, flags, targetInfo)
        targetInfo['style'] = styleLabel

    endfDoc = ''
    otherDocs = ''
    for style1 in self.styles:
        if isinstance(style1, stylesModule.Evaluated):
            if endfDoc == '':
                endfDoc = style1.documentation.endfCompatible.body
            elif style1.documentation.endfCompatible.body != '':
                otherDocs = '\n\n' + style1.documentation.endfCompatible.body + otherDocs
    if otherDocs != '':
        lines = endfDoc.split('\n')
        endfDoc = '\n'.join(lines[:5]) + otherDocs + '\n'.join(lines[5:])

    if endfDoc == '':
        evaluatedDate = evaluatedStyle.date.split('-')
        month = {1: 'JAN', 2: 'FEB', 3: 'MAR', 4: 'APR', 5: 'MAY', 6: 'JUN', 7: 'JUL', 8: 'AUG', 9: 'SEP',
                 10: 'OCT', 11: 'NOV', 12: 'DEC'}[int(evaluatedDate[1])]
        year = evaluatedDate[0][2:]
        evaluation = self.evaluation[:15]
        docHeader2 = [' %2d-%-2s-%3d LLNL       EVAL-%3s%2s Unknown' % (
        targetZ, chemicalElementMiscPoPsModule.symbolFromZ[targetZ], targetA, month, year),
                      '                      DIST-%s                       %8s   ' % (
                      time.strftime("%b%y").upper(), time.strftime("%Y%m%d")),
                      '----%15s   MATERIAL %4d' % (evaluation, MAT),
                      '-----INCIDENT %s DATA' %
                      {1: 'NEUTRON', 1001: 'PROTON', 1002: 'DEUTERON', 1003: 'TRITON', 2003: 'HELION', 2004: 'ALPHA'}[
                          projectileZA],
                      '------ENDF-6 FORMAT']
        endfDoc = ['Translated via GNDS to ENDF6 by FUDGE.',
                   '' ' ************************ C O N T E N T S ***********************']
    else:
        lines = endfDoc.split('\n')
        endfDoc = []
        for line in lines:
            if len(line) < 67:
                endfDoc.append(line)
            else:
                endfDoc += textwrap.wrap(line, 66, replace_whitespace=False, drop_whitespace=False)
        docHeader2 = []

        # update the documentation, including metadata on first 4 lines:
    if len([reac for reac in self.reactions if reac.isFission()]) > 0:
        LFI = True
    else:
        LFI = False
    LRP = -1
    if evaluatedStyle is not style:
        LRP = 2
    else:
        if self.resonances is not None:
            if self.resonances.resolved is None and self.resonances.unresolved is None:
                LRP = 0  # scattering radius only
            elif self.resonances.reconstructCrossSection:
                LRP = 1
            elif self.resonances.unresolved is not None:
                LRP = 1  # Self-shielding only.
                if self.resonances.resolved is not None:
                    LRP = 2
            else:
                LRP = 2  # Resonances provided for info only

    temperature = style.temperature.getValueAs('K')
    library = evaluatedStyle.library
    version = evaluatedStyle.version
    if library == 'ENDL':  # Additional ENDF meta-data. If the library is unknown, use NLIB = -1
        NVER, LREL, NMOD = 1, 1, 1
        NLIB = -1
    elif library == 'ecpl':  # For ECPL translations
        NVER, LREL, NMOD = 8, 1, 1
        NLIB = +1  # ENDF/A for ENDL/B-VIII.1
        temperature = 0.
    else:
        try:
            NVER, LREL, NMOD = map(int, version.split('.'))  # Version stored as '7.2.1'
        except Exception:
            print('Could not extract NVER/LREL/NMOD from version "%s". Defaulting to "0"' % version)
            NVER, LREL, NMOD = 0, 0, 0

        NLIB = -1
        for key, value in endf_endlModule.NLIBs.items():
            if value == library:
                NLIB = key

    NFOR = 6  # ENDF-6 format
    NSUB = 10 * IPART + ITYPE
    if self.isThermalNeutronScatteringLaw():
        NSUB = 12
    LDRV = 0
    if not isinstance(style, stylesModule.Evaluated):
        LDRV = 1

    if kwargs.get('NLIB', -1) != -1: NLIB = kwargs.get('NLIB')

    projectileMass = targetInfo['massTracker'].getMassAWR(projectileZA, asTarget=False)
    if self.projectile == IDsPoPsModule.photon:
        projectileMass = 0  # FIXME fix massTracker confusion in photo-atomic libraries between photon and electron
    docHeader = [endfFormatsModule.endfHeadLine(targetInfo['ZA'], targetInfo['mass'], LRP, LFI, NLIB, NMOD),
                 endfFormatsModule.endfHeadLine(levelEnergy_eV, STA, levelIndex, targetInfo['LISO'], 0, NFOR),
                 endfFormatsModule.endfHeadLine(projectileMass, EMAX, LREL, 0, NSUB, NVER),
                 endfFormatsModule.endfHeadLine(temperature, 0, LDRV, 0, len(docHeader2 + endfDoc), -1)]
    new_doc = documentationModule.Documentation('endf', '\n'.join(docHeader + docHeader2 + endfDoc))
    endfMFList[1][451] += endfFormatsModule.toEndfStringList(new_doc)

    endfFormatsModule.useRedsFloatFormat = _useRedsFloatFormat

    return endfFormatsModule.endfMFListToFinalFile(endfMFList, MAT, lineNumbers=lineNumbers)


reactionSuiteModule.ReactionSuite.toENDF6 = toENDF6


def addDecayGamma(reactionSuite, particle, baseMT, endfMFList, flags, targetInfo):
    MF = 12
    LP = 0
    MT = baseMT + particle.index
    gammaData = []
    levelEnergy_eV = particle.energy[0].float('eV')
    for decayMode in particle.decayData.decayModes:
        IDs = [product.pid for decay in decayMode.decayPath for product in decay.products]
        IDs.remove(IDsPoPsModule.photon)
        if len(IDs) != 1: raise Exception('Do not know how to handle this.')
        probability = decayMode.probability[0].value
        finalEnergy_eV = reactionSuite.PoPs[IDs[0]].energy[0].float('eV')
        _data = [finalEnergy_eV, probability]
        if decayMode.photonEmissionProbabilities:
            _data.append(decayMode.photonEmissionProbabilities[0].value)
        gammaData.append(_data)

    gammaData.sort(reverse=True)
    nGammas = len(gammaData)
    LGp = len(gammaData[0])
    endfMFList[MF][MT] = [
        endfFormatsModule.endfHeadLine(targetInfo['ZA'], targetInfo['mass'], 2, LGp - 1, MT - baseMT, 0),
        endfFormatsModule.endfHeadLine(levelEnergy_eV, 0., LP, 0, LGp * nGammas, nGammas)]

    endfMFList[MF][MT] += endfFormatsModule.endfNdDataList(gammaData)
    endfMFList[MF][MT].append(endfFormatsModule.endfSENDLineNumber())

    # Currently, assume all distributions are isotropic
    endfMFList[14][MT] = [endfFormatsModule.endfHeadLine(targetInfo['ZA'], targetInfo['mass'], 1, 0, nGammas, 0)]
    endfMFList[14][MT].append(endfFormatsModule.endfSENDLineNumber())
