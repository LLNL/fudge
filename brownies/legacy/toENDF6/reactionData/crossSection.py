# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import math

from xData import enums as xDataEnumsModule
from PoPs import IDs as IDsPoPsModule
from PoPs.families import nuclide as nuclideModule

from fudge import enums as enumsModule
from fudge import outputChannel as outputChannelModule
import fudge.reactionData.crossSection as crossSectionModule
import fudge.reactionData.doubleDifferentialCrossSection.chargedParticleElastic as CPElasticModule

from .. import gndsToENDF6 as gndsToENDF6Module
from .. import endfFormats as endfFormatsModule

def toENDF6(self, MT, endfMFList, targetInfo, level, LR):
    """
    Convert self into ENDF format

    :param int MT: The ENDF reaction designator, MT
    :param endfMFList:
    :param targetInfo:
    :param level:
    :param LR:
    """

    reactionSuite = targetInfo['reactionSuite']
    ZA, mass, QI, QM = targetInfo['ZA'], targetInfo['mass'], targetInfo['Q'], targetInfo['QM']
    if targetInfo['crossSectionMF'] == 23 and QI != 0:
        QI = -QI    # MF23 stores binding energy instead of Q-value
    if 'EFL' in targetInfo:
        QM = QI
        QI = targetInfo['EFL']
    else:
        if QM is None:
            QM = 0
            if MT in (2, 5):
                QM = QI
            elif MT == 4 and reactionSuite.projectile == IDsPoPsModule.neutron:
                # Q should be 0 except for excited-state targets:
                targetID = reactionSuite.target
                if targetID in reactionSuite.PoPs.aliases:
                    targetID = reactionSuite.PoPs[targetID].pid
                target = reactionSuite.PoPs[targetID]
                if isinstance(target, nuclideModule.Particle):
                    if target.index != 0:
                        QM = QI
            else:
                if targetInfo['reaction'] is not None:
                    if targetInfo['reaction'].outputChannel.genre == enumsModule.Genre.twoBody:
                        QM = QI + level
                    else:
                        QM = targetInfo['reaction'].thresholdQAs('eV', final=True)
                else:
                    QM = QI

    if targetInfo['reaction'] is not None:
        fissionEnergyReleases = targetInfo['reaction'].outputChannel.fissionFragmentData.fissionEnergyReleases
        if len(fissionEnergyReleases) > 0:
            QM = fissionEnergyReleases[0].nonNeutrinoEnergy.data.coefficients[0]
            QI = QM
    interpolationFlatData, crossSectionFlatData = gndsToENDF6Module.getForm(targetInfo['style'], self).toENDF6Data(MT, endfMFList, targetInfo, level)
    if interpolationFlatData is not None:
        MF = targetInfo['crossSectionMF']
        endfMFList[MF][MT] = [endfFormatsModule.endfHeadLine(ZA, mass, 0, 0, 0, 0)]
        endfMFList[MF][MT].append(endfFormatsModule.endfContLine(QM, QI, 0, LR, len(interpolationFlatData) / 2, len(crossSectionFlatData) / 2))
        endfMFList[MF][MT] += endfFormatsModule.endfInterpolationList(interpolationFlatData)
        endfMFList[MF][MT] += endfFormatsModule.endfDataList(crossSectionFlatData, crossSectionData=True)
        endfMFList[MF][MT].append(endfFormatsModule.endfSENDLineNumber())

crossSectionModule.Component.toENDF6 = toENDF6

def toENDF6Data( self, MT, endfMFList, targetInfo, level ) :

    endfInterpolation = gndsToENDF6Module.gndsToENDFInterpolationFlag( self.interpolation )
    crossSectionFlatData = []
    for xy in self.copyDataToXYs( ) : crossSectionFlatData += xy
    return( [ len( crossSectionFlatData ) / 2, endfInterpolation ], crossSectionFlatData )

crossSectionModule.XYs1d.toENDF6Data = toENDF6Data

def toENDF6Data( self, MT, endfMFList, targetInfo, level ) :

    interpolationFlatData, crossSectionFlatData = [], []
    counter = 0
    lastX, lastY = None, None
    for region in self :
        ENDFInterpolation = gndsToENDF6Module.gndsToENDFInterpolationFlag( region.interpolation )
        data = region.copyDataToXYs( )
        if( lastX is not None ) :
            isLogY = region.interpolation in [xDataEnumsModule.Interpolation.loglin, xDataEnumsModule.Interpolation.loglog]
            if lastY == data[0][1]:
                data = data[1:]
            elif lastY == 0 and isLogY:
                interpolationFlatData[-2] += 1
            elif ENDFInterpolation == lastENDFInterpolation:
                interpolationFlatData = interpolationFlatData[:-2]
        counter += len( data )
        interpolationFlatData.append( counter )
        interpolationFlatData.append( ENDFInterpolation )
        for xy in data : crossSectionFlatData += xy
        lastX, lastY = data[-1]
        lastENDFInterpolation = ENDFInterpolation
    return( interpolationFlatData, crossSectionFlatData )

crossSectionModule.Regions1d.toENDF6Data = toENDF6Data

def toENDF6Data( self, MT, endfMFList, targetInfo, level ) :

    regions = []
    axes = None
    mergeFlag = None
    for term in (self.background.resolvedRegion, self.background.unresolvedRegion, self.background.fastRegion):
        if term is None: continue

        if axes is None:
            axes = term.data.axes
        elif term.data.axes != axes:
            raise NotImplementedError("Inconsistent axes/units inside background cross section")

        if isinstance(term.data, crossSectionModule.XYs1d):
            regions.append(term.data.copy())
        elif isinstance(term.data, crossSectionModule.Regions1d):
            regions.append(term.data[0].copy())

        if mergeFlag is not None:
            reg2 = regions.pop()
            reg1 = regions[-1]
            if 'deletePoint' in mergeFlag:
                data = reg1.copyDataToXYs()[:-1] + reg2.copyDataToXYs()[1:]
            else:
                data = reg1.copyDataToXYs()[:-1] + reg2.copyDataToXYs()
            reg1.setData(data)

        if isinstance(term.data, crossSectionModule.Regions1d):
            for region in term.data[1:]:
                regions.append(region.copy())

        mergeFlag = targetInfo['ENDFconversionFlags'].get(term)

    regions_ = crossSectionModule.Regions1d(axes=axes)
    for region in regions:
        regions_.append( region )

    return regions_.toENDF6Data( MT, endfMFList, targetInfo, level )

crossSectionModule.ResonancesWithBackground.toENDF6Data = toENDF6Data

def toENDF6Data( self, MT, endfMFList, targetInfo, level ) :

    if isinstance( self.link, CPElasticModule.CoulombPlusNuclearElastic.Form ):
        if isinstance( self.link.data, CPElasticModule.nuclearAmplitudeExpansion.NuclearAmplitudeExpansion ):
            # actual cross section is in the Coulomb expansion, MF=3 stores a dummy value:
            return( [ 2, 2 ], [ self.link.domainMin,1, self.link.domainMax,1] )
        elif isinstance(self.link.data, CPElasticModule.nuclearPlusInterference.NuclearPlusInterference):
            # ENDF-6 uses different convention than GNDS or ENDL: factor of 2pi smaller
            crossSection = crossSectionModule.XYs1d( data = self.link.data.crossSection.data / (2 * math.pi) )
            return crossSection.toENDF6Data( MT, endfMFList, targetInfo, level )
    else :
        return( [ None, None ] )

crossSectionModule.Reference.toENDF6Data = toENDF6Data
