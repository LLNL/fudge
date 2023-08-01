# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import sys
import subprocess
import math
import copy

from pqu import PQU as PQUModule

from PoPs import IDs as IDsPoPsModule

from LUPY import locateFudgeBin as locateBinFolderModule

from LUPY import subprocessing as subprocessingModule
from LUPY import times as timesModule

from xData import enums as xDataEnumsModule
from xData import axes as axesModule
from xData import XYs1d as XYs1dModule
from xData import multiD_XYs as multiD_XYsModule
from xData import regions as regionsModule

from fudge.productData.distributions import angular as angularModule
from fudge.productData.distributions import energy as energyModule
from fudge.productData.distributions import energyAngular as energyAngularModule

from . import specialCases as specialCasesModule

linlin = xDataEnumsModule.Interpolation.linlin
linlog = xDataEnumsModule.Interpolation.linlog
loglin = xDataEnumsModule.Interpolation.loglin
loglog = xDataEnumsModule.Interpolation.loglog
flat = xDataEnumsModule.Interpolation.flat

GNDS2ProcessingInterpolationQualifiers = {
    xDataEnumsModule.InterpolationQualifier.none: 'direct', 
    xDataEnumsModule.InterpolationQualifier.unitBase: xDataEnumsModule.InterpolationQualifier.unitBase,
    xDataEnumsModule.InterpolationQualifier.unitBaseUnscaled: xDataEnumsModule.InterpolationQualifier.unitBaseUnscaled,
    xDataEnumsModule.InterpolationQualifier.correspondingPoints: 'cumulativepoints',
    xDataEnumsModule.InterpolationQualifier.correspondingPointsUnscaled: 'cumulativepoints-unscaled' }

versionStr = "xndfgenTransferMatrix: version 1.0"
doubleFmt = "%19.12e"
lowerEps = 1e-8
upperEps = 1e-8
startOfNewData = "\n# Start data\n"
startOfNewSubData = "# Start sub data"

transferMatrixExecute = locateBinFolderModule.locateMerced()

def twoBodyTransferMatrix( style, tempInfo, productFrame, crossSection, angularData, Q, weight = None, comment = None ) :
    """
    Generate input and call processing code to generate a transfer matrix for two-body angular distribution.
    If the distribution is actually made up of two different forms in different energy regions, this function
    calls itself in the two regions and sums the result.
    """

    if( isinstance( angularData, angularModule.Isotropic2d ) ) : angularData = angularData.toPointwise_withLinearXYs( )

    if( isinstance( angularData, angularModule.XYs2d ) ) :
        TM1, TME = twoBodyTransferMatrix2( style, tempInfo, crossSection, angularData, Q, productFrame, comment = comment )
    elif( isinstance( angularData, angularModule.Regions2d ) ) :
        if len(angularData) == 1: return twoBodyTransferMatrix(style, tempInfo, productFrame, crossSection, angularData[0], Q, weight, comment)
        TM1s, TMEs = [], []
        lowestBound, highestBound = angularData[0].domainMin, angularData[-1].domainMax
        weightAxes = axesModule.Axes(2)

        reactionSuite = tempInfo['reactionSuite']
        projectileName = reactionSuite.projectile
        maxGroupEnergy = style.transportables[projectileName].group.boundaries.values[-1]
        for iRegion, region in enumerate( angularData ) :
            if region.domainMin >= maxGroupEnergy:
                break
            if( iRegion == 0 ) :
                weightData = [ [ lowestBound, 1 ], [ region.domainMax, 0 ], [ highestBound, 0 ] ]
            elif( iRegion == len( angularData ) - 1 ) :
                weightData = [ [ lowestBound, 0 ], [ region.domainMin, 1 ], [ highestBound, 1 ] ]
            else :
                weightData = [ [ lowestBound, 0 ], [ region.domainMin, 1 ], [ region.domainMax, 0 ], [ highestBound, 0 ] ]
            _weight = XYs1dModule.XYs1d(data=weightData, axes=weightAxes, interpolation=flat)

            tempInfo['workFile'].append( 'r%s' % iRegion )
            try :
                TM1, TME = twoBodyTransferMatrix2( style, tempInfo, crossSection, region, Q,
                        productFrame, weight = _weight, comment = comment )
            except :
                del tempInfo['workFile'][-1]
                raise
            del tempInfo['workFile'][-1]
            TM1s.append( TM1 )
            TMEs.append( TME )
        TM1 = addTMs( TM1s )
        TME = addTMs( TMEs )
    else :
        raise Exception( 'Unsupported P(mu|E) = %s' % angularData.moniker )

    return( TM1, TME )
        
def twoBodyTransferMatrix2( style, tempInfo, crossSection, angularData, Q, productFrame, weight = None, comment = None ) :
    """
    Helper function for twoBodyTransferMatrix. This should be called separately for every interpolation 
    region within a Regions2d form.
    """

    reactionSuite = tempInfo['reactionSuite']

#
# This next section handles the case where Merced is having issues for photo-nuclear data at threshold and for small (<1e-9 MeV) outgoing energies.
#
    modifiedProductGroupIndex = 0
    modifiedProductGroup = None
    if reactionSuite.projectile == IDsPoPsModule.photon and tempInfo['productName'] == IDsPoPsModule.neutron:

        cutoffEnergy = PQUModule.PQU( 1e-9, 'MeV' ).getValueAs(tempInfo['incidentEnergyUnit'])
        productName = tempInfo['productName']
        productGroupBoundaries = style.transportables[productName].group
        for modifiedProductGroupIndex, boundary in enumerate(productGroupBoundaries.boundaries.values):
            if boundary > cutoffEnergy: break
        if modifiedProductGroupIndex > 0:
            modifiedProductGroup = productGroupBoundaries.copy( )
            modifiedProductGroupBoundaries = modifiedProductGroup.boundaries.values.values[:modifiedProductGroupIndex+3] # 3 should be 3.    # 2 because the first one returned by Merced extends to 0.0 product energy.
            modifiedProductGroup.boundaries.values.values = modifiedProductGroup.boundaries.values.values[modifiedProductGroupIndex:]

    logFile = tempInfo['logFile']
    workDir = tempInfo['workDir']

    if( isinstance( angularData, angularModule.Recoil ) ) : angularData = angularData.getNumericalDistribution( )

    s = versionStr + '\n'
    if( isinstance( angularData[0], angularModule.XYs1d ) ) :
        s += "Process: 'two body transfer matrix'\n"
    elif( isinstance( angularData[0], angularModule.Legendre ) ) :
        s += "Process: 'Legendre two body transfer matrix'\n"
    else :
        raise Exception( 'Unsupported P(mu) = %s' % angularData[0].moniker )

    s += "Reaction's Q value: %s\n" % PQUModule.floatToShortestString( ( Q ), 12 )

    s += commonDataToString( comment, style, tempInfo, crossSection, productFrame, photonFrame = xDataEnumsModule.Frame.centerOfMass,
            modifiedProductGroup=modifiedProductGroup)
    s += angularToString( angularData, crossSection, weight = weight, twoBody = True )
    TM1, TME = executeCommand( logFile, transferMatrixExecute, s, workDir, tempInfo['workFile'], tempInfo['restart'],
            productOffset = modifiedProductGroupIndex )

    if modifiedProductGroup is None: return TM1, TME

    projectileName = reactionSuite.projectile
    projectileGroupBoundaries = style.transportables[projectileName].group.boundaries.values
    for projectileIndex in range(len(TM1)):
        row = TM1[projectileIndex]
        data = []
        for productIndex in row: data += list(map(abs, row[productIndex]))
        total = sum(data)
        if total != 0.0:
            dE = projectileGroupBoundaries[projectileIndex+1] - projectileGroupBoundaries[projectileIndex]  # This ignores any flux dependency on energy.
            groupFlux = 1 / dE
            TM1atThreshold, TMEatThreshold = specialCasesModule.twoBodyPhotoNuclearAtThreshold(tempInfo['masses'], Q, 
                modifiedProductGroupBoundaries, crossSection, angularData, tempInfo['legendreMax'], groupFlux)
            norm = TM1atThreshold.pop(-1)[0]
            if norm != 0.0: norm = row[len(TM1atThreshold)][0] / norm
            TMEatThreshold.pop(-1)
            for index, TM1atThresholdRow in enumerate(TM1atThreshold):
                TM1atThresholdRow[0] *= norm
                row[index] = TM1atThresholdRow
                TME[projectileIndex][index] = TMEatThreshold[index]
            break

    return TM1, TME

def wholeAtomScattering( style, tempInfo, productFrame, formFactor, realAnomalousFactor = None, imaginaryAnomalousFactor = None, comment = None ) :

    logFile = tempInfo['logFile']
    workDir = tempInfo['workDir']
    incidentEnergyUnit = tempInfo['incidentEnergyUnit']
    waveLengthUnit = formFactor.axes[1].unit
    inverseWaveLengthToEnergyFactor = PQUModule.PQU( 1, 'hplanck * c * %s' % waveLengthUnit ).getValueAs( incidentEnergyUnit )
    electronRadius = PQUModule.PQU( 1, "e**2 / ( 4 * pi * eps0 * me * c**2 )" )
    ThomsonScatteringCrossSection = 8 * math.pi * electronRadius**2 / 3

    s1  = versionStr + '\n'
    s1 += "Process: 'coherent scattering'\n"
    s1 += 'inverseWaveLengthToEnergyFactor: %s\n' % PQUModule.floatToShortestString( inverseWaveLengthToEnergyFactor, 12 )
    s1 += 'ThompsonScattering: %s\n' % PQUModule.floatToShortestString( ThomsonScatteringCrossSection.getValueAs( 'b' ), 12 )

    s1 += commonDataToString( comment, style, tempInfo, None, productFrame )
    s1 += startOfNewData
    s1 += '\n'.join( twoDToString( "FormFactorData", formFactor ) )
# BRB: FIXME For realAnomalousFactor and imaginaryAnomalousFactor need to convert energy to proper unit.
    if( realAnomalousFactor is not None ) : 
        s1 += startOfNewData
        s1 += '\n'.join( twoDToString( "realAnomalousFactor", realAnomalousFactor ) )
    if( imaginaryAnomalousFactor is not None ) :
        s1 += startOfNewData
        s1 += '\n'.join( twoDToString( "imaginaryAnomalousFactor", imaginaryAnomalousFactor ) )
    return( executeCommand( logFile, transferMatrixExecute, s1, workDir, tempInfo['workFile'], tempInfo['restart'] ) )

def comptonScattering( style, tempInfo, productFrame, scatteringFactor, comment = None ) :

    logFile = tempInfo['logFile']
    workDir = tempInfo['workDir']
    incidentEnergyUnit = tempInfo['incidentEnergyUnit']
    waveLengthUnit = scatteringFactor.axes[1].unit
    inverseWaveLengthToEnergyFactor = PQUModule.PQU( 1, 'hplanck * c * %s' % waveLengthUnit ).getValueAs( incidentEnergyUnit )
    electronRadius = PQUModule.PQU( 1, "e**2 / ( 4 * pi * eps0 * me * c**2 )" )
    ThomsonScatteringCrossSection = 8 * math.pi * electronRadius**2 / 3

    s1  = versionStr + '\n'
    s1 += "Process: 'Compton scattering'\n"
    s1 += 'inverseWaveLengthToEnergyFactor: %s\n' % PQUModule.floatToShortestString( inverseWaveLengthToEnergyFactor, 12 )
    s1 += 'ThompsonScattering: %s\n' % PQUModule.floatToShortestString( ThomsonScatteringCrossSection.getValueAs( 'b' ), 12 )
    s1 += 'Electron mass: %s\n' % PQUModule.floatToShortestString( PQUModule.PQU( 1, "me * c**2" ).getValueAs( incidentEnergyUnit ), 14 )

    s1 += commonDataToString( comment, style, tempInfo, None, productFrame )
    s1 += '\n'.join( twoDToString( "ScatteringFactorData", scatteringFactor ) )
    return( executeCommand( logFile, transferMatrixExecute, s1, workDir, tempInfo['workFile'], tempInfo['restart'] ) )

def ENDFEMuEpP_TransferMatrix( style, tempInfo, productFrame, crossSection, angularEnergyData, multiplicity, comment = None ) :
    """This is ENDF MF = 6, LAW = 7 type data."""

    logFile = tempInfo['logFile']
    workDir = tempInfo['workDir']

    s  = versionStr + '\n'
    s += "Process: ENDF Double differential EMuEpP data\n"

    s += commonDataToString( comment, style, tempInfo, crossSection, productFrame, multiplicity = multiplicity )
    s += EMuEpPDataToString( angularEnergyData )
    return( executeCommand( logFile, transferMatrixExecute, s, workDir, tempInfo['workFile'], tempInfo['restart'] ) )

def ENDLEMuEpP_TransferMatrix( style, tempInfo, crossSection, productFrame, angularData, EMuEpPData, multiplicity, comment = None ) :
    """This is LLNL I = 1, 3 type data."""

    logFile = tempInfo['logFile']
    workDir = tempInfo['workDir']

    s  = versionStr + '\n'
    s += "Process: 'Double differential EMuEpP data transfer matrix'\n"

    s += commonDataToString( comment, style, tempInfo, crossSection, productFrame, multiplicity = multiplicity )
    s += angularToString( angularData, crossSection )
    s += EMuEpPDataToString( EMuEpPData )
    return( executeCommand( logFile, transferMatrixExecute, s, workDir, tempInfo['workFile'], tempInfo['restart'] ) )

def ELEpP_TransferMatrix( style, tempInfo, crossSection, productFrame, LEEpPData, multiplicity, comment = None ) :
    """
    This is for ENDL I = 4 data with l > 0. This form is deprecated.
    """

    logFile = tempInfo['logFile']
    workDir = tempInfo['workDir']

    s  = versionStr + '\n'
    s += "Process: 'Legendre EEpP data transfer matrix'\n"

    s += commonDataToString( comment, style, tempInfo, crossSection, productFrame, multiplicity = multiplicity )
    s += LEEpPDataToString( LEEpPData )
    return( executeCommand( logFile, transferMatrixExecute, s, workDir, tempInfo['workFile'], tempInfo['restart'] ) )

def EEpMuP_TransferMatrix( style, tempInfo, productFrame, crossSection, energyAngularData, multiplicity, comment = None ) :
    """This is ENDF MT=6 LAW=1, LANG=1,11-15 type data."""

    if( isinstance( energyAngularData[0][0], energyAngularModule.Legendre ) ) :
        return( Legendre_TransferMatrix( style, tempInfo, productFrame, crossSection, 
                energyAngularData, multiplicity, comment = comment ) )
    elif( isinstance( energyAngularData[0][0], energyAngularModule.XYs1d ) ) :
        return( doubleDifferential_EEpMuP(  style, tempInfo, productFrame, crossSection,
                energyAngularData, multiplicity, comment = comment ) )
    else :
        raise Exception( "Unsupported P(mu) for P(E',mu|E) = %s" % energyAngularData[0][0].moniker )

def Legendre_TransferMatrix( style, tempInfo, productFrame, crossSection, LegendreData, multiplicity, comment = None ) :

    logFile = tempInfo['logFile']
    workDir = tempInfo['workDir']

    s  = versionStr + '\n'
    s += "Process: Legendre energy-angle data\n"

    s += commonDataToString( comment, style, tempInfo, crossSection, productFrame, multiplicity = multiplicity )
    s += LegendreDataToString( LegendreData )

    return( executeCommand( logFile, transferMatrixExecute, s, workDir, tempInfo['workFile'], tempInfo['restart'] ) )

def doubleDifferential_EEpMuP( style, tempInfo, productFrame, crossSection, EEpMuPData, multiplicity, comment = None ) :
    
    logFile = tempInfo['logFile']
    workDir = tempInfo['workDir']
    
    s  = versionStr + '\n'
    s += "Process: pointwise energy-angle data\n"
    
    s += commonDataToString( comment, style, tempInfo, crossSection, productFrame, multiplicity = multiplicity )
    s += EEpMuPDataToString( EEpMuPData )
    
    return( executeCommand( logFile, transferMatrixExecute, s, workDir, tempInfo['workFile'], tempInfo['restart'] ) )

def uncorrelated_EMuP_EEpP_TransferMatrix( style, tempInfo, crossSection, productFrame, angularData, energyData, 
        multiplicity, comment = None, weight = None ) :

    logFile = tempInfo['logFile']
    workDir = tempInfo['workDir']

    legendreMax = None
    if( angularData.isIsotropic() ) :
        if productFrame == xDataEnumsModule.Frame.lab: legendreMax = 0

    sSpecific = ''
    if( isinstance( energyData, energyModule.FunctionalBase ) ) :
        sProcess, sSpecific, sData = energyFunctionToString( energyData, weight = weight )
        weight = None
    elif( isinstance( energyData, energyModule.WeightedFunctionals ) ) :
        TM1s, TMEs = [], []
        for weight in energyData :
            TM1, TME = uncorrelated_EMuP_EEpP_TransferMatrix( style, tempInfo, crossSection, productFrame, angularData, 
                    weight.functional, multiplicity, comment = comment, weight = weight.weight )
            TM1s.append( TM1 )
            TMEs.append( TME )
        TM1 = addTMs( TM1s )
        TME = addTMs( TMEs )
        return( TM1, TME )
    elif( isinstance( energyData, energyModule.Regions2d ) ) :
        TM1s, TMEs = [], []
        axes = axesModule.Axes(2)
        for i1, region in enumerate( energyData ) :
            weight = XYs1dModule.XYs1d(data=[[region.domainMin, 1.], [region.domainMax, 1.]], axes=axes, interpolation=flat)
            tempInfo['workFile'].append( 'r%s' % i1 )
            try :
                TM1, TME = uncorrelated_EMuP_EEpP_TransferMatrix( style, tempInfo, crossSection, productFrame, angularData,
                        region, multiplicity, comment = comment, weight = weight )
            except :
                del tempInfo['workFile'][-1]
                raise
            del tempInfo['workFile'][-1]
            TM1s.append( TM1 )
            TMEs.append( TME )
        TM1 = addTMs( TM1s )
        TME = addTMs( TMEs )
        return( TM1, TME )
    elif( isinstance( energyData, energyModule.PrimaryGamma ) ) :
        return( primaryGammaAngularData( style, tempInfo, crossSection, energyData, angularData, multiplicity = multiplicity, comment = comment ) )
    else :
        if( angularData.isIsotropic() ) :
            sProcess = "Process: isotropic table\n"
            sData = EEpPDataToString( energyData )
        else:
            sProcess = "Process: 'Uncorrelated energy-angle data transfer matrix'\n"
            sData = angularToString( angularData, crossSection )
            sData += EEpPDataToString( energyData )

    sCommon = commonDataToString( comment, style, tempInfo, crossSection, productFrame, multiplicity = multiplicity, weight = weight, legendreMax = legendreMax )
    s = versionStr + '\n' + sProcess + sSpecific + sCommon + sData

    return( executeCommand( logFile, transferMatrixExecute, s, workDir, tempInfo['workFile'], tempInfo['restart'] ) )

def discreteGammaAngularData(style, tempInfo, gammaEnergy, crossSection, angularData, multiplicity, comment = None):
    """Currently, only isotropic (i.e., l = 0) data are returned. That is, lMax and angularData are ignored. This routine is also used
    for pair-production which pass angularData as None."""

    from fudge.processing import miscellaneous as miscellaneousModule

    reactionSuite = tempInfo['reactionSuite']
    projectileName = reactionSuite.projectile
    projectileGroupBoundaries = style.transportables[projectileName].group.boundaries.values
    productName = tempInfo['productName']
    productGroupBoundaries = style.transportables[productName].group.boundaries.values

    nProj = len(projectileGroupBoundaries) - 1
    nProd = len(productGroupBoundaries) - 1
    TM_1, TM_E = {}, {}
    for i1 in range(nProj):
        TM_1_i1, TM_E_i1 = {}, {}
        for i2 in range(nProd):
            TM_1_i1[i2] = [0.]
            TM_E_i1[i2] = [0.]
        TM_1[i1] = TM_1_i1
        TM_E[i1] = TM_E_i1

    indexEp = -1
    for Ep in productGroupBoundaries :
        if gammaEnergy <= Ep:
            break
        indexEp += 1
    indexEp = min(max(indexEp, 0), nProd - 1)

    if multiplicity.domainMin < productGroupBoundaries[-1]:
        xsecTimesMult = miscellaneousModule.groupTwoFunctionsAndFlux(style, tempInfo, crossSection, multiplicity, norm = tempInfo['groupedFlux'])
        for indexE in range(nProj):
            x = xsecTimesMult[indexE]
            TM_1[indexE][indexEp] = [x]
            TM_E[indexE][indexEp] = [x * gammaEnergy]

    return TM_1, TM_E

def primaryGammaAngularData( style, tempInfo, crossSection, energyData, angularData, multiplicity = 1, comment = None ) :
    """
    Currently, only isotropic (i.e., l = 0) data are returned. That is, lMax and angularData are ignored. massRatio is the
    target mass divide by the sum of the projectile and target masses.

    This function perform the integration 

        TM = \int_g dE \int_h dE' S(E) M(E) f(E) P(E -> E') / \int_g dE f(E)

    where \int_g is the integral of E from E_i to E_{i+1}, \int_g is the integral of E' from E'_j to E'_{j+1}, S(E) is the cross section,
    M(E) is the products multiplicity, f(E) is the flux weighting and P(E -> E') is the probability for a projectile of energy E producing
    a primary gamma of energy E' for binding energy bindingEnergy. This function assumes that M(E) is a constant. For primary gamma's captured
    into binding energy bindingEnergy P(E -> E') = deltaFunction( E' - ( bindingEnergy + massRatio E ) ) where massRatio = mt / ( mp + mt ),
    mp is the projectile's mass and mt is the target's mass. Note, this formula is from the ENDF manual which ignores the recoil of the residual
    nucleus.
    """

    reactionSuite = tempInfo['reactionSuite']
    projectileName = reactionSuite.projectile
    projectileGroupBoundaries = style.transportables[projectileName].group.boundaries.values
    productName = tempInfo['productName']
    productGroupBoundaries = style.transportables[productName].group.boundaries.values
    flux0 = style.flux.getFluxAtLegendreOrder( 0 )
    groupedFlux = tempInfo['groupedFlux']

    bindingEnergy = energyData.value * energyData.axes[1].unitConversionFactor( tempInfo['incidentEnergyUnit'] )
    massRatio = energyData.massRatio

    nProj = len( projectileGroupBoundaries ) - 1
    nProd = len( productGroupBoundaries ) - 1
    TM_1, TM_E = {}, {}
    for i1 in range( nProj ) :
        TM_1[i1] = {}
        TM_E[i1] = {}
        for i2 in range( nProd ) :
            TM_1[i1][i2] = [ 0. ]
            TM_E[i1][i2] = [ 0. ]
    Eg2 = bindingEnergy + massRatio * projectileGroupBoundaries[0]
    for indexEo, Eo in enumerate( productGroupBoundaries ) :
        if( Eg2 <= Eo ) : break
    indexEo = min( max( indexEo - 1, 0 ), nProd - 1 )

    EMin, EMax = crossSection.domainMin, crossSection.domainMax
    axes = axesModule.Axes(2, labelsUnits = {0: ('energy_out', tempInfo['incidentEnergyUnit']), 1: (crossSection.axes[1].label, crossSection.axes[1].unit)})
    Egp = XYs1dModule.XYs1d( data = [ [ EMin, bindingEnergy + massRatio * EMin ], [ EMax, bindingEnergy + massRatio * EMax ] ], axes = axes )

    for indexEi in range( nProj ) :
        Ei2 = projectileGroupBoundaries[indexEi + 1]
        Eg2 = bindingEnergy + massRatio * Ei2
        EiMin = projectileGroupBoundaries[indexEi]
        while( True ) :
            incrementIndexEo, EiMax = 0, Ei2
            if( indexEo < ( nProd - 1 ) ) :
                if( Eg2 > productGroupBoundaries[indexEo + 1] ) :
                    incrementIndexEo = 1
                    EiMax = ( productGroupBoundaries[indexEo + 1] - bindingEnergy ) / massRatio
            TM_1[indexEi][indexEo][0] = float( crossSection.integrateTwoFunctions( flux0, 
                    domainMin = EiMin, domainMax = EiMax ) / groupedFlux[indexEi] )
            TM_E[indexEi][indexEo][0] = float( crossSection.integrateThreeFunctions( flux0, Egp, 
                    domainMin = EiMin, domainMax = EiMax ) / groupedFlux[indexEi] )
            if( incrementIndexEo == 0 ) : break
            EiMin = EiMax
            if( indexEo < ( nProd - 1 ) ) : indexEo += 1
    return( TM_1, TM_E )

def NBodyPhaseSpace( style, tempInfo, crossSection, numberOfProducts, mTotal, Q, multiplicity = 1, comment = None ) :

    logFile = tempInfo['logFile']
    workDir = tempInfo['workDir']

    s  = versionStr + '\n'
    s += "Process: phase space spectrum\n"
    s += "Number of particles: %s\n" % numberOfProducts
    s += "Total mass: %s\n" % PQUModule.floatToShortestString( mTotal, 14 )
    s += "Q value: %s\n" % PQUModule.floatToShortestString( Q, 12 )
    productFrame = xDataEnumsModule.Frame.centerOfMass
    s += commonDataToString( comment, style, tempInfo, crossSection, productFrame, multiplicity = multiplicity )
    return( executeCommand( logFile, transferMatrixExecute, s, workDir, tempInfo['workFile'], tempInfo['restart'] ) )

def KalbachMann_TransferMatrix( style, tempInfo, crossSection, particlesData, KalbachMannData, multiplicity = 1, comment = None ) :

    logFile = tempInfo['logFile'], 
    workDir = tempInfo['workDir']

    energy_in_unit = 'MeV'
    s  = versionStr + '\n'
    s += "Process: Kalbach spectrum\n"
    s += "Projectile's ZA: %s\n" % particlesData['projectile']['ZA']
    s += "Target's ZA: %s\n" % particlesData['target']['ZA']
    s += "Product's ZA: %s\n" % particlesData['product']['ZA']
    s += "Compound's mass: %s\n" % PQUModule.floatToShortestString( particlesData['compound']['mass'], 14 )
    productFrame = xDataEnumsModule.Frame.centerOfMass
    s += commonDataToString( comment, style, tempInfo, crossSection, productFrame, multiplicity = multiplicity, energy_in_unit = energy_in_unit )
    s += KalbachMannDataToString( KalbachMannData, energy_in_unit )
    return( executeCommand( logFile, transferMatrixExecute, s, workDir, tempInfo['workFile'], tempInfo['restart'] ) )

def executeCommand( logFile, executable, cmd, workDir, workFile, restart, productOffset = 0 ) :
    # FIXME remove unused argument logFile?
    """
    Runs executable (presumably Merced) on a single input, parses resulting output file and returns two transfer matrices:
     - 1st has constant weight (conserves number of particles)
     - 2nd is weighted by E' (conserves energy)

    If executable fails, create a '.err' file in the working directory (for easier searching)

    @param logFile: currently unused! Remove argument?
    @param executable: compiled executable to run, generally FUDGE/bin/merced
    @param cmd: full text of executable input file
    @param workDir: directory where executable will be run. processProtare creates new workDir for each projectile/target/temperature
    @param workFile: list of form [reaction identifier, product identifier], e.g. ['r0001','n'] or ['r0045','photon__b']. Used to generate unique file names.
    @param restart: boolean. If True, check whether executable was previously run for this input and if so skip re-running (useful for resuming jobs that timed out).
    """

    def checkNegative_l0( TM_EEpL, weight, f ) :

        negative_l0Counter, largestNegative = 0, 0.
        for ig in sorted( TM_EEpL.keys( ) ) :
            TM_EpL = TM_EEpL[ig]
            for ih in sorted( TM_EpL.keys( ) ) :
                l0Cell = TM_EpL[ih][0]
                if( l0Cell < 0 ) :
                    negative_l0Counter += 1
                    largestNegative = min( largestNegative, l0Cell )
                    fOut.write( '    negative l=0 for weight %s at row %3d column %3d: %s\n' % ( weight, ig, ih, l0Cell ) )
        return( negative_l0Counter )

    if executable is None : raise Exception( 'Merced executable not found.' )

    if not os.path.exists( workDir ): os.makedirs( workDir )
    workFile = '_'.join( workFile ) + '_in'
    fullFileName = os.path.join( workDir, workFile )
    outputFile = fullFileName + ".out"
    infoFile = fullFileName + ".info"
    fOut = open( infoFile, 'a' )

    if restart:
        if (os.path.exists(fullFileName) and os.path.exists(outputFile) and os.stat(outputFile).st_size != 0
                and not os.path.exists(fullFileName + ".err")):
            with open(fullFileName) as fin:
                previousInput = fin.read()
            if previousInput == cmd:
                # This input has already been successfully run (e.g. by previous processProtare run that timed out).
                # Ensure variables are initialized:
                fOut.write("  Reading previously computed results\n")
                t0 = timesModule.Times( )
                status = 0
            else:
                restart = False
        else:
            restart = False

    if not restart:
        with open( fullFileName, 'w' ) as dataFile:
            dataFile.write( cmd )
        t0 = timesModule.Times( )
        try :
            status, stdout, stderr = subprocessingModule.executeCommand( [ executable, '-output', outputFile+".tmp", fullFileName ],
                stdout = infoFile, stderr = subprocess.STDOUT )
            os.rename(outputFile+".tmp", outputFile)
        except :
            fErr = open( fullFileName + ".err", "w" )
            fErr.close( )
            raise

    # executable is finished, parse results
    fOut.write( str( t0 ) + "\n" )
    if( status != 0 ) :
        print("status = ", status)
        raise Exception( 'Transfer failed for %s %s' % ( executable, fullFileName ) )
    TM1, TME, N = parseOutputFile( fullFileName + '.out', productOffset )
    negative_l0Counter = 0
    if( TM1 is not None ) : negative_l0Counter = checkNegative_l0( TM1, "0", fOut )
    if( TME is not None ) : negative_l0Counter += checkNegative_l0( TME, "E", fOut )
    if( negative_l0Counter > 0 ) : print('    WARNING: %d negative l=0 elements found in transfer matrix' % negative_l0Counter)
    fOut.close( )

    if( workDir is None ) :
        os.remove( fullFileName + '.out' )
        os.remove( fullFileName + '.info' )
        if( workFile == [] ) :
            dataFile.delete( )
        else :
            dataFile.close( )
            os.remove( fullFileName )
    return( TM1, TME )

def parseOutputFile( file, productOffset, firstLine = None ) :

    def checkKeyValue( lineNumber, keyValue, _key ) :

        try :
            key, value = keyValue.split( '=' )
        except :
            raise Exception( 'bad key/value = "%s" at line %d while looking for key = "%d"' % ( keyValue, lineNumber, _key ) )
        key = key.strip( )
        if( key != _key ) : raise KeyError( 'bad key = "%s" at line %d' % ( key, lineNumber ) )
        try :
            return( int( value ) )
        except :
            raise TypeError( 'bad value type %s at line %d' % ( type(value), lineNumber ) )

    def parseCrossSection( paras, lineNumber, ls, file ) :

        s1 = ls[lineNumber].split( ':' )                        # This line should start as 'Cross section: n = '
        n1 = int( s1[1].split( '=' )[1] )
        s1 = ls[lineNumber].split( ':' )                        # This line should start as 'Interpolation: '
        interpolation = s1[1].strip( )
        lineNumber += n1 + 1
        return( lineNumber )

    def parseTransferMatrix( paras, lineNumber, ls, file, productOffset ) :

        s1 = ls[lineNumber].split( ":" )                      # This line should start as 'Integrals, weight = ...'
        EBins = checkKeyValue( lineNumber, s1[1], 'numEinBins' )
        if( EBins > paras['numEinBins'] ) : 
            raise Exception( 'Bad transfer numEinBins = %d != %d in file %s\n%s' % ( EBins, paras['numEinBins'], file, ls[lineNumber] ) )
        lMax = paras['outputLegendreOrder']
        lMaxPlus1 = lMax + 1
        TM = {}
        for i1 in range( paras['numEinBins'] ) : 
            TMEp = {}
            for i2 in range( paras['numEoutBins'] + productOffset ) : TMEp[i2] = lMaxPlus1 * [ 0. ]
            TM[i1] = TMEp
        for i1 in range( EBins ) :
            lineNumber += 1
            s1 = ls[lineNumber].split( ":" )                  # This should start as 'EinBin = ...'
            EBin = checkKeyValue( lineNumber, s1[0], 'EinBin' )
            if( not( 0 <= EBin < EBins ) ) : raise Exception( 'Bad transfer EinBin at lineNumber = %d in file %s\n%s' % ( lineNumber, file, ls[lineNumber] ) )
            EpBins = checkKeyValue( lineNumber, s1[1], 'numEoutBins' )
            if( EpBins > paras['numEoutBins'] ) : raise Exception( 'Bad transfer numEoutBins at lineNumber = %d in file %s\n%s' % ( lineNumber, file, ls[lineNumber] ) )
            for i2 in range( EpBins ) :
                lineNumber += 1
                rowLs = list( map( float, ls[lineNumber].split( ) ) )
                if( len( rowLs ) > lMaxPlus1 ) : raise Exception( 'Bad transfer lRows at lineNumber = %d in file %s\n%s' % ( lineNumber, file, ls[lineNumber] ) )
                for i3, Cl in enumerate( rowLs ) : TM[i1][i2+productOffset][i3] = Cl
        return( TM, lineNumber )

    crossSection = None
    TM1 = None
    TME = None
    try :
        f = open( file )
    except :
        raise Exception( 'Could not open transfer file = %s' % file )
    ls = f.readlines( )
    f.close( )
    rightFirstLine = "merced: version 1\n"
    if isinstance(firstLine, str) : rightFirstLine = firstLine
    if( ls[0] != rightFirstLine ) : raise Exception( 'Bad transfer file version in file %s\n%s' % ( file, ls[0] ) )
    lineNumber = 1
    n = len( ls )
    paras = { 'outputLegendreOrder' : None, 'numEinBins' : None, 'numEoutBins' : None }
    while( lineNumber < n ) :
        s = ls[lineNumber].split( ":" )
        key = s[0]
        if( key in [ 'outputLegendreOrder', 'numEinBins', 'numEoutBins' ] ) :
            if( paras[key] is not None ) : raise Exception( "Bad transfer file, multiple key = %s lines in file %s" % ( key, file ) )
            try :
                paras[key] = int( s[1] )
            except :
                raise Exception( "Bad transfer file, could not convert data for key = %s to integer\n%s in file %s" % ( key, ls[lineNumber], file ) )
        elif( key == "Comment" ) : 
            pass
        elif( key == "maxIncidentEnergyGroup" ) : 
            paras[key] = int( s[1] )
        elif( key == "Cross section" ) : 
            if( crossSection is not None ) : raise Exception( "Bad transfer file, multiple cross section datasets in file %s" % file )
            lineNumber = parseCrossSection( paras, lineNumber, ls, file )
        elif( key == "Integrals, weight = 1" ) :
            if( TM1 is not None  ) : raise Exception( "Bad transfer file, multiple weight = 1 datasets in file %s" % file )
            TM1, lineNumber = parseTransferMatrix( paras, lineNumber, ls, file, productOffset )
        elif( key == "Integrals, weight = E'" ) :
            if( TME is not None  ) : raise Exception( "Bad transfer file, multiple weight = E' datasets in file %s" % file )
            TME, lineNumber = parseTransferMatrix( paras, lineNumber, ls, file, productOffset )
        else :
            raise Exception( 'Bad transfer file key in file %s\n%s' % ( file, ls[lineNumber] ) )
        lineNumber += 1
    return( TM1, TME, paras )

def commonDataToString( comment, style, tempInfo, crossSection, productFrame, multiplicity = None, energy_in_unit = None, weight = None,
        photonFrame = xDataEnumsModule.Frame.lab, legendreMax = None, modifiedProductGroup = None ) :

    reactionSuite = tempInfo['reactionSuite']
    projectileName = reactionSuite.projectile
    projectileGroupBoundaries = style.transportables[projectileName].group
    productName = tempInfo['productName']
    productGroupBoundaries = style.transportables[productName].group
    if( modifiedProductGroup is not None ): productGroupBoundaries = modifiedProductGroup

    if( legendreMax is None ) : legendreMax = tempInfo['legendreMax']
    s  = "outputLegendreOrder: %d\n" % legendreMax
    if( comment is not None ) : s += "Comment: %s\n" % comment

    masses = tempInfo['masses']
    for particle in sorted( masses.keys( ) ) :
        mass = masses[particle]
        if( mass is not None ) : s += "%s's mass: %s\n" % ( particle, PQUModule.floatToShortestString( ( masses[particle] ), 14 ) )
    s += "Projectile Frame: %s\n" % reactionSuite.projectileFrame
    s += GBToString( 'Projectile', projectileGroupBoundaries, energy_in_unit )
    maximumEnergy_in = projectileGroupBoundaries.boundaries.values[-1]
    s += GBToString( 'Product', productGroupBoundaries, energy_in_unit )

    s += fluxToString( style.flux, energy_in_unit )
    if( crossSection is not None ) : s += crossSectionToString( style, crossSection, energy_in_unit )
    if( multiplicity is not None ) : s += multiplicityToString( style, multiplicity, energy_in_unit = energy_in_unit )
    if( productName == IDsPoPsModule.photon ) : productFrame = photonFrame
    s += "\nProduct Frame: %s\n" % productFrame
    if( weight is not None ) : s += startOfNewData + '\n'.join( twoDToString( "weight", weight, addExtraBlankLine = False ) )
    return( s )

def GBToString( name, gb, energy_in_unit ) :

    boundaries = gb.boundaries
    factor = 1
    if( energy_in_unit is not None ) : factor = boundaries.unitConversionFactor( energy_in_unit )
    s  = startOfNewData
    s += "%s's group boundaries: n = %d\n" % ( name, len( boundaries.values ) )
    s += oneDToString( boundaries.values, factor )
    return( s )

def fluxToString( flux, energy_in_unit ) :

    flux0 = flux.getFluxAtLegendreOrder( 0 )
    if( energy_in_unit is not None ) : flux0 = flux0.convertAxisToUnit( 1, energy_in_unit )
    interpolationStr, s = twoDToString( "", flux0, addHeader = False )
    s  = [ "Fluxes: n = %d" % len( flux0 ) ]
    s.append( interpolationStr )
    for Ein, Cls in flux0 :
        s.append( 'Ein: ' + ( doubleFmt % Ein ) + ': n = 1' )
        s.append( '  ' + ( doubleFmt % Cls ) )
    return( startOfNewData + '\n'.join( s ) )

def crossSectionToString( style, crossSection, energy_in_unit = None ) :

    if( energy_in_unit is not None ) : crossSection = crossSection.copy( ).convertAxisToUnit( 1, energy_in_unit )
    return( startOfNewData + '\n'.join( twoDToString( "Cross section", crossSection ) ) )

def multiplicityToString( style, _multiplicity, energy_in_unit = None ) :

    form = style.findFormMatchingDerivedStyle( _multiplicity )
    multiplicity = form.toPointwise_withLinearXYs( accuracy = 1e-5, upperEps = 1e-8 )
    if( energy_in_unit is not None ) : multiplicity = multiplicity.copy( ).convertAxisToUnit( 1, energy_in_unit )
    return( startOfNewData + '\n'.join( twoDToString( "Multiplicity", multiplicity ) ) )

def angularToString( angularData, crossSection, weight = None, twoBody = False, changeInterpolationQualifierWarning = False ) :

    weightString = ''
    if( weight is not None ) : weightString = startOfNewData + '\n'.join( twoDToString( "weight", weight, addExtraBlankLine = False ) ) + '\n'

    if( isinstance( angularData, angularModule.Isotropic2d ) ) :
        s = "Angular data: n = 2\n"
        s += 'Incident energy interpolation: %s %s\n' % ( linlin, GNDS2ProcessingInterpolationQualifiers[xDataEnumsModule.InterpolationQualifier.none] )
        s += 'Outgoing cosine interpolation: %s\n' % linlin
        for x in [ crossSection.domainMin, crossSection.domainMax ] : s += " %s : n = 2\n   -1 0.5\n   1 0.5\n" % ( doubleFmt % x )

    elif( isinstance( angularData, angularModule.XYs2d ) ) :
        if( isinstance( angularData[0], angularModule.XYs1d ) ) :
            s = "Angular data: n = %s\n" % len( angularData )

            interpolationQualifier = angularData.interpolationQualifier
            if( angularData.interpolationQualifier == xDataEnumsModule.InterpolationQualifier.none and not twoBody ) :    # Don't use direct interpolation!
                interpolationQualifier = xDataEnumsModule.InterpolationQualifier.unitBase

            if( angularData.interpolation not in [ linlin, flat ] ) :
                _angularData = angularData.copy( )
                _angularData.interpolation = linlin
                if( changeInterpolationQualifierWarning ) :
                    print('WARNING: interpolation %s is ignored, treated as %s' % ( angularData.interpolation, _angularData.interpolation ))

            if( twoBody ) :
                s += 'Incident energy interpolation: %s %s\n' % \
                        ( angularData.interpolation, GNDS2ProcessingInterpolationQualifiers[interpolationQualifier] )
            else :
                s += 'Incident energy interpolation: %s %s\n' % ( angularData.interpolation,
                        GNDS2ProcessingInterpolationQualifiers[interpolationQualifier] )

            interpolations = set( POfMu.interpolation for POfMu in angularData )
            linearizeXYs = False
            if( len( interpolations ) > 1 ) :
                linearizeXYs = True
            elif( ( linlin not in interpolations ) or ( flat not in interpolations ) ) :
                linearizeXYs = True
            interpolation = angularData[0].interpolation
            if( linearizeXYs ) : interpolation = linlin
            s += 'Outgoing cosine interpolation: %s\n' % interpolation

            s += threeDListToString( angularData, linearizeXYs = linearizeXYs )

        elif( isinstance( angularData[0], angularModule.Legendre ) ) :
            interpolation = angularData.interpolation
            if( angularData.interpolationQualifier != xDataEnumsModule.InterpolationQualifier.none ) :
                raise Exception( 'interpolation qualifier = "%s" is not supported.' % angularData.interpolationQualifier )
            if( interpolation not in [ linlin, flat ] ) :
                _angularData = angularData.copy( )
                _angularData.interpolation = linlin
                if( changeInterpolationQualifierWarning ) :
                    print('WARNING: ignoring interpolation "%s" and using "%s" instead' %
                          ( angularData.interpolation, _angularData.interpolation ) )
                angularData = _angularData

            s = [ "Legendre coefficients: n = %s" % len( angularData ) ]

            if( twoBody ) :
                interpolationStr = 'Interpolation: %s' % angularData.interpolation
            else :
                interpolationStr = 'Interpolation: %s' % angularData.interpolation
            s.append( interpolationStr )

            for energy_in in angularData :
                s.append( ' Ein: %s:  n = %d' % ( PQUModule.floatToShortestString( energy_in.outerDomainValue, 12 ), len( energy_in ) ) )
                for coefficient in energy_in : s.append( "%s" % PQUModule.floatToShortestString( coefficient, 12 ) )
            s = '\n'.join( s ) + '\n'
        else :
            raise Exception( "angular data = %s not supported" % angularData[0].moniker )

    else :
        raise Exception( "angular data = %s not supported" % angularData.moniker )

    return( startOfNewData + weightString + s )

def energyFunctionToString( energyData, weight = None ) :

    def getParameter( data, label ) :

        _linlin = data.toPointwise_withLinearXYs( lowerEps = lowerEps, upperEps = upperEps )
        sData = [ startOfNewData, "%s: n = %d" % ( label, len( _linlin ) ) ]
        sData.append( 'Interpolation: %s' % linlin )
        for x, y in _linlin : sData.append( '%s %s' % ( PQUModule.floatToShortestString( x, 12 ), PQUModule.floatToShortestString( y, 12 ) ) )
        return( sData )

    sData = []
    if( isinstance( energyData, ( energyModule.SimpleMaxwellianFission, energyModule.Evaporation, energyModule.Watt ) ) ) :
        sData.append( 'U: %s' % PQUModule.floatToShortestString( energyData.U.getValueAs( energyData.parameter1.data.axes[-1].unit ), 12 ) )

    parameter1 = energyData.parameter1.data.toPointwise_withLinearXYs( accuracy = 1e-5, lowerEps = lowerEps, upperEps = upperEps )

    if( energyData.parameter2 is not None ) :
        parameter2 = energyData.parameter2.data.toPointwise_withLinearXYs( accuracy = 1e-5, lowerEps = lowerEps, upperEps = upperEps )

    if( weight is None ) :
        axes = axesModule.Axes(2)
        weight = XYs1dModule.XYs1d(data=[[parameter1.domainMin, 1.], [parameter1.domainMax, 1.]], axes=axes, interpolation=flat)
    sData += twoDToString( "weight", weight, addExtraBlankLine = False )

    if( isinstance( energyData, energyModule.GeneralEvaporation ) ) :
        sProcess = 'Process: general evaporation\n'
        sSpecific = ''
        sData += getParameter( parameter1, 'Theta' )
        sData += getParameter( parameter2, 'g' )
    elif( isinstance( energyData, energyModule.SimpleMaxwellianFission ) ) :
        sProcess = 'Process: Maxwell spectrum\n'
        sSpecific = ''
        sData += getParameter( parameter1, 'Theta' )
    elif( isinstance( energyData, energyModule.Evaporation ) ) :
        sProcess = 'Process: Evaporation spectrum\n'
        sSpecific = ''
        sData += getParameter( parameter1, 'Theta' )
    elif( isinstance( energyData, energyModule.Watt ) ) :
        sProcess = 'Process: Watt spectrum\n'
        sSpecific = 'Conserve: number\n'
        sData += getParameter( parameter1, 'a' )
        sData += getParameter( parameter2, 'b' )
    elif( isinstance( energyData, energyModule.MadlandNix ) ) :
        sProcess = 'Process: Madland-Nix spectrum\n'
        sSpecific = 'Conserve: number\n'
        sSpecific += 'EFL: %s\n' % PQUModule.floatToShortestString( energyData.EFL.getValueAs( energyData.parameter1.data.axes[-1].unit ), 12 )
        sSpecific += 'EFH: %s\n' % PQUModule.floatToShortestString( energyData.EFH.getValueAs( energyData.parameter1.data.axes[-1].unit ), 12 )
        sData += getParameter( parameter1, 'TM' )
    else :
        raise Exception( 'Unsupport energy functional = %s' % energyData.moniker )

    sData.append( '' )

    return( sProcess, sSpecific, startOfNewData + '\n'.join( sData ) )

def EEpPDataToString( EEpPData, changeInterpolationQualifierWarning = False ) :

    s  = startOfNewData
    s += "EEpPData: n = %d\n" % len( EEpPData )

    qualifier = EEpPData.interpolationQualifier
    if qualifier == xDataEnumsModule.InterpolationQualifier.none:
        if EEpPData.interpolation != flat:
            if( changeInterpolationQualifierWarning ) : print('    WARNING: changing interpolation qualifier from None to unitbase')
            qualifier = xDataEnumsModule.InterpolationQualifier.unitBase
    s += 'Incident energy interpolation: %s %s\n' % ( EEpPData.interpolation, GNDS2ProcessingInterpolationQualifiers[qualifier] )

    linearizeXYs = False
    for POfEp in EEpPData :
        if( not( isinstance( POfEp, energyModule.XYs1d ) ) ) : linearizeXYs = True
    if( not( linearizeXYs ) ) :
        interpolations = set( POfEp.interpolation for POfEp in EEpPData )
        if( len( interpolations ) > 1 ) :
            linearizeXYs = True
        else :
            interpolation = interpolations.pop( )
            if( interpolation not in ( linlin, flat ) ) :
                linearizeXYs = True
        interpolation = EEpPData[0].interpolation
    if( linearizeXYs ) : interpolation = linlin
    s += 'Outgoing energy interpolation: %s\n' % interpolation

    s += threeDToString( EEpPData, linearizeXYs = linearizeXYs )
    return( s )

def LEEpPDataToString( LEEpPData ) :

    s  = startOfNewData
    s += "LEEpPData: n = %d\n" % len( LEEpPData )
    s += 'Incident energy interpolation: %s %s\n' % (linlin, GNDS2ProcessingInterpolationQualifiers[xDataEnumsModule.InterpolationQualifier.unitBase])
    s += 'Outgoing energy interpolation: %s\n' % linlin
    for l_EEpP in LEEpPData :
        s += startOfNewSubData + '\n'
        s += "  l = %d: n = %d\n" % ( int( l_EEpP.outerDomainValue ), len( l_EEpP ) )
        s += threeDToString( l_EEpP )
    return( s )

def LegendreDataToString( LegendreData, changeInterpolationQualifierWarning = False ) :

    s = [ 'Legendre data by incident energy:  n = %d' % len( LegendreData ) ]


    qualifier = LegendreData.interpolationQualifier
    if qualifier == xDataEnumsModule.InterpolationQualifier.none:
        qualifier = xDataEnumsModule.InterpolationQualifier.unitBase
    s.append( 'Incident energy interpolation: %s %s' % ( LegendreData.interpolation, GNDS2ProcessingInterpolationQualifiers[qualifier] ) )

    interpolation = LegendreData[0].interpolation
    if( interpolation in [ linlin, flat ] ) :
        s.append( 'Outgoing energy interpolation: %s' % interpolation )
    else :
        if( changeInterpolationQualifierWarning ) :
            print ('    WARNING: converting interpolation from "%s" to "%s".' % (interpolation, linlin))
        s.append( 'Outgoing energy interpolation: %s' % linlin )

    for EEpCl in LegendreData :
        s.append( " Ein: %s:  n = %d" % ( PQUModule.floatToShortestString( EEpCl.outerDomainValue, 12 ), len( EEpCl ) ) )
        for EpCl in EEpCl :
            s.append( startOfNewSubData )
            s.append( "  Eout: %s:  n = %d" % ( PQUModule.floatToShortestString( EpCl.outerDomainValue, 12 ), len( EpCl ) ) )
            for Cl in EpCl : s.append( "   %s" % PQUModule.floatToShortestString( Cl, 12 ) )
    s.append( '' )
    return( startOfNewData + '\n'.join( s ) )

def EMuEpPDataToString( EMuEpPData ) :

    s  = startOfNewData
    s += "EMuEpPData: n = %d\n" % len( EMuEpPData )

    qualifier = EMuEpPData.interpolationQualifier
    if qualifier not in [xDataEnumsModule.InterpolationQualifier.none, xDataEnumsModule.InterpolationQualifier.unitBase]:
        raise Exception( 'interpolation qualifier = %s is not supported' % qualifier )
    if( EMuEpPData.interpolation != linlin ) :
        raise Exception( 'interpolation = "%s" is not supported' % EMuEpPData.interpolation )
    s += 'Incident energy interpolation: %s %s\n' % (linlin, GNDS2ProcessingInterpolationQualifiers[xDataEnumsModule.InterpolationQualifier.unitBase])

    s += 'Outgoing cosine interpolation: %s %s\n' % (linlin, GNDS2ProcessingInterpolationQualifiers[xDataEnumsModule.InterpolationQualifier.unitBase])

    s += 'Outgoing energy interpolation: %s\n' % linlin

    for muEpP in EMuEpPData :
        s += startOfNewSubData + '\n'
        E = muEpP.outerDomainValue
        s += "  E = %s: n = %d\n" % ( doubleFmt % E, len( muEpP ) )
        s += threeDToString( muEpP )
    return( s )

def EEpMuPDataToString( EEpMuPData ) :

    s  = startOfNewData
    s += "EEpMuPData: n = %d\n" % len( EEpMuPData )

    qualifier = EEpMuPData.interpolationQualifier
    if qualifier in [xDataEnumsModule.InterpolationQualifier.none, xDataEnumsModule.InterpolationQualifier.correspondingPoints]:
        print('WARNING: converting interpolation qualifier from "%s" to "%s".' % (qualifier, xDataEnumsModule.InterpolationQualifier.unitBase))
        qualifier = xDataEnumsModule.InterpolationQualifier.unitBase
    s += 'Incident energy interpolation: %s %s\n' % ( EEpMuPData.interpolation, qualifier )

    outgoingInterpolations = set([xys2d.interpolation for xys2d in EEpMuPData])
    if len(outgoingInterpolations) > 1: raise Exception( 'Multiple outgoing energy interpolation(s) %s not supported' % outgoingInterpolations)
    s += 'Outgoing energy interpolation: %s direct\n' % EEpMuPData[0].interpolation

    s += 'Outgoing cosine interpolation: %s\n' % EEpMuPData[0][0].interpolation

    for EpMuP in EEpMuPData :
        s += startOfNewSubData + '\n'
        s += "  E = %s: n = %d\n" % ( doubleFmt % EpMuP.outerDomainValue, len( EpMuP ) )
        s += threeDToString( EpMuP )
    return( s )

def KalbachMannDataToString( KalbachMannData, energy_in_unit ) :

    nForm = 4
    if( KalbachMannData.aSubform.data is None ) : nForm = 3
    if( nForm == 4 ) : raise Exception( 'This is currently not implemented' )

    fSubform = KalbachMannData.fSubform.data
    rSubform = KalbachMannData.rSubform.data

    factor = fSubform.domainUnitConversionFactor( energy_in_unit ) # Energies passed to get_transfer must be in energy_in_unit, traditionally MeV.
    _fSubform = multiD_XYsModule.XYs2d( )
    for xys in fSubform :
        _xys = XYs1dModule.XYs1d( data = xys.nf_pointwiseXY.scaleOffsetXAndY( xScale = factor, yScale = 1 / factor ), outerDomainValue = factor * xys.outerDomainValue )
        _fSubform.append( _xys )

    _rSubform = multiD_XYsModule.XYs2d( )
    for xys in rSubform :
        _xys = XYs1dModule.XYs1d( data = xys.nf_pointwiseXY.scaleOffsetXAndY( xScale = factor ), outerDomainValue = factor * xys.outerDomainValue )
        _rSubform.append( _xys )

    s  = startOfNewData
    s += "Kalbach probabilities: n = %s\n" % len( _fSubform )
    s += "Incident energy interpolation: %s %s\n" % (linlin, GNDS2ProcessingInterpolationQualifiers[xDataEnumsModule.InterpolationQualifier.unitBase])
    s += "Outgoing energy interpolation: %s\n" % fSubform[0].interpolation
    s += threeDListToString( _fSubform )

    s += "Kalbach r parameter: n = %s\n" % len( rSubform )
    s += "Incident energy interpolation: %s %s\n" % (rSubform.interpolation, GNDS2ProcessingInterpolationQualifiers[xDataEnumsModule.InterpolationQualifier.unitBaseUnscaled])
    s += "Outgoing energy interpolation: %s\n" % rSubform[0].interpolation
    s += threeDListToString( _rSubform )

    return( s )

def oneDToString( data, factor = 1 ) :

    i = 1
    a = [ "" ]
    for d in data :
        a.append( doubleFmt % ( factor * d ) )
        if( ( i % 10 ) == 0 ) : a.append( "\n" )
        i += 1
    if( ( i % 10 ) != 1 ) : a.append( "\n" )
    return( " ".join( a ) )

def twoDToString( label, data, addHeader = True, addExtraBlankLine = True ) :

    if( isinstance( data, XYs1dModule.XYs1d ) ) :
        if( data.interpolation not in [ linlin, flat ] ) :
            data = data.toPointwise_withLinearXYs( lowerEps = lowerEps, upperEps = upperEps )
        interpolationStr = 'Interpolation: %s' % data.interpolation
    elif( isinstance( data, regionsModule.Regions1d ) ) :
        data = data.toPointwise_withLinearXYs( lowerEps = lowerEps, upperEps = upperEps )
        interpolationStr = 'Interpolation: %s' % data.interpolation
    else :
        print ('ERROR from %s' % __file__)
        print(type(data))
        raise TypeError( 'Unsupported data type: see prior output' )
    fmt = " %s %s" % ( doubleFmt, doubleFmt )
    a = [ fmt % ( x, y ) for x, y in data ]
    if( not( addHeader ) ) : return( interpolationStr, a )
    a.insert( 0, interpolationStr )
    a.insert( 0, "%s: n = %d" % ( label, len( data ) ) )
    if( addExtraBlankLine ) : a.append( '' )
    return( a )

def threeDToString( data, linearizeXYs = False ) :

    a = []
    wFmt = " %s : n = %%d" % ( doubleFmt )
    xyFmt = "   %s %s" % ( doubleFmt, doubleFmt )
    for w_xy in data :
        a.append( startOfNewSubData )
        if( type( w_xy ) == type( [] ) ) :       # ???? This is a kludge until Legendre EEpP data is represented like angular EMuP data.
            raise Exception( 'need to update for new gnds forms' )
            x, yz = xyz
        else :
            w = w_xy.outerDomainValue
            xy = w_xy
            if( linearizeXYs ) : xy = w_xy.toPointwise_withLinearXYs( accuracy = 1e-3, upperEps = 1e-8 )
        a.append( wFmt % ( w, len( xy ) ) )
        for x, y in xy : a.append( xyFmt % ( x, y ) )
    a.append( "" )
    return( "\n".join( a ) )

def threeDListToString( data, linearizeXYs = False ) :

    a = []
    wFmt = " %s : n = %%d" % ( doubleFmt )
    xyFmt = "   %s %s" % ( doubleFmt, doubleFmt )
    for w_xy in data :
        a.append( startOfNewSubData )
        w = w_xy.outerDomainValue
        if( linearizeXYs ) : w_xy = w_xy.toPointwise_withLinearXYs( accuracy = 1e-3, upperEps = 1e-8 )
        a.append( wFmt % ( w, len( w_xy ) ) )
        for x, y in w_xy : a.append( xyFmt % ( x, y ) )
    a.append( "" )
    return( "\n".join( a ) )

def addTMs( TMs ) :

    TM = TMs.pop( 0 )
    for i1 in sorted( TM.keys( ) ) :
        TM_energy_in = TM[i1]
        for i2 in sorted( TM_energy_in.keys( ) ) :
            TM_energy_in_out = TM_energy_in[i2]
            for TM2 in TMs :
                for i3 in range( len( TM_energy_in_out ) ) : TM_energy_in_out[i3] += TM2[i1][i2][i3]
                
    return( TM )
