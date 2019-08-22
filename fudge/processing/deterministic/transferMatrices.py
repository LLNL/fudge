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

import os, subprocess

from fudge.core.utilities import brb
from fudge.core.utilities import fudgeFileMisc, subprocessing, times

from xData import standards as standardsModule
from xData import axes as axesModule
from xData import XYs as XYsModule
from xData import multiD_XYs as multiD_XYsModule
from xData import regions as regionsModule

from fudge.gnd.reactionData import crossSection as crossSectionModule

from fudge.gnd.productData import multiplicity as multiplicityModule

from fudge.gnd.productData.distributions import angular as angularModule
from fudge.gnd.productData.distributions import energyAngular as energyAngularModule

linlin = standardsModule.interpolation.linlinToken
linlog = standardsModule.interpolation.linlogToken
loglin = standardsModule.interpolation.loglinToken
loglog = standardsModule.interpolation.loglogToken
flat = standardsModule.interpolation.flatToken

GND2ProcessingInterpolationQualifiers = { 
    '' : 'direct', 
    standardsModule.interpolation.unitBaseToken : standardsModule.interpolation.unitBaseToken,
    standardsModule.interpolation.unitBaseUnscaledToken : standardsModule.interpolation.unitBaseUnscaledToken,
    standardsModule.interpolation.correspondingPointsToken : 'cumulativepoints',
    standardsModule.interpolation.correspondingPointsUnscaledToken : 'cumulativepoints-unscaled' }

versionStr = "xndfgenTransferMatrix: version 1.0"
doubleFmt = "%19.12e"
lowerEps = 1e-8
upperEps = 1e-8
startOfNewData = "\n# Start data\n"
startOfNewSubData = "# Start sub data"

transferMatrixExecute = 'getTransferMatrix'
if( os.path.exists( transferMatrixExecute ) ) :
    srcPath = os.path.abspath( './' )
else :
    srcPath = os.path.split( os.path.abspath( __file__ ) )[0]
transferMatrixExecute = '%s/%s' % ( srcPath, transferMatrixExecute )

def twoBodyTransferMatrix( style, tempInfo, productFrame, crossSection, angularData, Q, weight = None, comment = None ) :
    """
    Generate input and call processing code to generate a transfer matrix for two-body angular distribution.
    If the distribution is actually made up of two different forms in different energy regions, this function
    calls itself in the two regions and sums the result.
    """

    if( isinstance( angularData, angularModule.isotropic ) ) : angularData = angularData.toPointwise_withLinearXYs( )

    if( isinstance( angularData, angularModule.XYs2d ) ) :
        TM1, TME = twoBodyTransferMatrix2( style, tempInfo, crossSection, angularData, Q, productFrame, comment = comment )
    elif( isinstance( angularData, angularModule.regions2d ) ) :
        TM1s, TMEs = [], []
        lowestBound, highestBound = angularData[0].domainMin( ), angularData[-1].domainMax( )
        weightAxes = axesModule.axes( )
        for iRegion, region in enumerate( angularData ) :
                if( iRegion == 0 ) :
                    weightData = [ [ lowestBound, 1 ], [ region.domainMax( ), 0 ], [ highestBound, 0 ] ]
                elif( iRegion == len( angularData ) - 1 ) :
                    weightData = [ [ lowestBound, 0 ], [ region.domainMin( ), 1 ], [ highestBound, 1 ] ]
                else :
                    weightData = [ [ lowestBound, 0 ], [ region.domainMin( ), 1 ], [ region.domainMax( ), 0 ], [ highestBound, 0 ] ]
                _weight = XYsModule.XYs1d( data = weightData, axes = weightAxes, accuracy = 1e-6 )

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
    region within a regions2d form.
    """

    from fudge.gnd.productData import distributions

    logFile = tempInfo['logFile']
    workDir = tempInfo['workDir']

    if( isinstance( angularData, distributions.angular.recoil ) ) : angularData = angularData.getNumericalDistribution( )

    s = versionStr + '\n'
    if( isinstance( angularData[0], angularModule.XYs1d ) ) :
        s += "Process: 'two body transfer matrix'\n"
    elif( isinstance( angularData[0], angularModule.Legendre ) ) :
        s += "Process: 'Legendre two body transfer matrix'\n"
        s += "Quadrature method: adaptive\n"
    else :
        raise Exception( 'Unsupported P(mu) = %s' % angularData[0].moniker )

    if( comment is not None ) : s += "Comment: %s\n" % comment
    s += "Reaction's Q value: %s\n" % floatToString( Q )

    s += commonDataToString( style, tempInfo, crossSection, productFrame )
    s += angularToString( angularData, crossSection, weight = weight, twoBody = True )
    return( executeCommand( logFile, transferMatrixExecute, s, workDir, tempInfo['workFile'] ) )

def wholeAtomScattering( processInfo, projectile, product, formFactor, comment = None ) :

    workDir = tempInfo['workDir']

    s  = versionStr + '\n'
    s += "Process: 'whole atom scattering'\n"

    if( comment is not None ) : s += "Comment: %s\n" % comment
    s += commonDataToString( style, tempInfo, crossSection, productFrame )
    s += startOfNewData
    s += '\n'.join( twoDToString( "FormFactorData", formFactor ) )
    return( executeCommand( logFile, transferMatrixExecute, s, workDir, tempInfo['workFile'] ) )

def comptonScattering( processInfo, projectile, product, formFactor, comment = None ) :

    logFile = tempInfo['logFile']
    workDir = tempInfo['workDir']

    s  = versionStr + '\n'
    s += "Process: 'Compton scattering'\n"

    if( comment is not None ) : s += "Comment: %s\n" % comment
    s += commonDataToString( style, tempInfo, crossSection, standardsModule.frames.labToken )
    s += '\n'.join( twoDToString( "ScatteringFactorData", formFactor ) )
    return( executeCommand( logFile, transferMatrixExecute, s, workDir, tempInfo['workFile'] ) )

def ENDFEMuEpP_TransferMatrix( style, tempInfo, productFrame, crossSection, angularEnergyData, multiplicity, comment = None ) :
    """This is ENDF MF = 6, LAW = 7 type data."""

    logFile = tempInfo['logFile']
    workDir = tempInfo['workDir']

    s  = versionStr + '\n'
    s += "Process: ENDF Double differential EMuEpP data\n"

    if( comment is not None ) : s += "Comment: %s\n" % comment
    s += "Quadrature method: adaptive\n"
    s += commonDataToString( style, tempInfo, crossSection, productFrame, multiplicity = multiplicity )
    s += EMuEpPDataToString( angularEnergyData )
    return( executeCommand( logFile, transferMatrixExecute, s, workDir, tempInfo['workFile'] ) )

def EMuEpP_TransferMatrix( processInfo, projectile, product, crossSection, angularData, EMuEpPData, multiplicity, comment = None ) :
    """This is LLNL I = 1, 3 type data."""

    logFile = tempInfo['logFile']
    workDir = tempInfo['workDir']

    s  = versionStr + '\n'
    s += "Process: 'Double differential EMuEpP data transfer matrix'\n"

    if( comment is not None ) : s += "Comment: %s\n" % comment
    s += commonDataToString( style, tempInfo, crossSection, productFrame, multiplicity = multiplicity )
    s += angularToString( angularData, crossSection )
    s += EMuEpPDataToString( EMuEpPData )
    return( executeCommand( logFile, transferMatrixExecute, s, workDir, tempInfo['workFile'] ) )

def ELEpP_TransferMatrix( processInfo, projectile, product, crossSection, LEEpPData, multiplicity, comment = None ) :
    """
    This is for ENDL I = 4 data with l > 0. This form is deprecated.
    """

    logFile = tempInfo['logFile']
    workDir = tempInfo['workDir']

    s  = versionStr + '\n'
    s += "Process: 'Legendre EEpP data transfer matrix'\n"
    if( comment is not None ) : s += "Comment: %s\n" % comment

    s += commonDataToString( style, tempInfo, crossSection, productFrame, multiplicity = multiplicity )
    s += LEEpPDataToString( LEEpPData )
    return( executeCommand( logFile, transferMatrixExecute, s, workDir, tempInfo['workFile'] ) )

def EEpMuP_TransferMatrix( style, tempInfo, productFrame, crossSection, energyAngularData, multiplicity, comment = None ) :
    """This is ENDF MT=6 LAW=1, LANG=1,11-15 type data."""

    if( isinstance( energyAngularData[0][0], energyAngularModule.Legendre ) ) :
        return( Legendre_TransferMatrix( style, tempInfo, productFrame, crossSection, 
                energyAngularData, multiplicity, comment = comment ) )
    else :
        raise Exception( "Unsupported P(mu) for P(E',mu|E) = %s" % energyAngularData[0][0].moniker )

def Legendre_TransferMatrix( style, tempInfo, productFrame, crossSection, LegendreData, multiplicity, comment = None ) :

    logFile = tempInfo['logFile']
    workDir = tempInfo['workDir']

    s  = versionStr + '\n'
    s += "Process: Legendre energy-angle data\n"
    if( comment is not None ) : s += "Comment: %s\n" % comment

    s += "Quadrature method: adaptive\n"
    s += commonDataToString( style, tempInfo, crossSection, productFrame, multiplicity = multiplicity )
    s += LegendreDataToString( LegendreData )

    return( executeCommand( logFile, transferMatrixExecute, s, workDir, tempInfo['workFile'] ) )

def uncorrelated_EMuP_EEpP_TransferMatrix( style, tempInfo, crossSection, productFrame, angularData, energyData, 
        multiplicity, comment = None, weight = None ) :

    from fudge.gnd.productData.distributions import energy as energyModule

    logFile = tempInfo['logFile']
    workDir = tempInfo['workDir']

    sSpecific = ''
    if( isinstance( energyData, energyModule.functionalBase ) ) :
        sProcess, sSpecific, sData = energyFunctionToString( energyData, weight = weight )
        weight = None
    elif( isinstance( energyData, energyModule.weightedFunctionals ) ) :
        TM1s, TMEs = [], []
        for weight in energyData :
            TM1, TME = uncorrelated_EMuP_EEpP_TransferMatrix( style, tempInfo, crossSection, productFrame, angularData, 
                    weight.functional, multiplicity, comment = comment, weight = weight )
            TM1s.append( TM1 )
            TMEs.append( TME )
        TM1 = addTMs( TM1s )
        TME = addTMs( TMEs )
        return( TM1, TME )
    elif( isinstance( energyData, energyModule.regions2d ) ) :
        TM1s, TMEs = [], []
        axes = axesModule.axes( )
        for i1, region in enumerate( energyData ) :
            weight = XYsModule.XYs1d( data = [ [ region.domainMin( ), 1. ], [ region.domainMax( ), 1. ] ], axes = axes, accuracy = 1e-6 )
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
    elif( isinstance( energyData, energyModule.primaryGamma ) ) :
        return( primaryGammaAngularData( style, tempInfo, crossSection, energyData, angularData, multiplicity = multiplicity, comment = comment ) )
    else :
        sProcess = "Process: 'Uncorrelated energy-angle data transfer matrix'\n"
        sData = angularToString( angularData, crossSection )
        sData += EEpPDataToString( energyData )
    if( comment is not None ) : sProcess += "Comment: %s\n" % comment
    sCommon = commonDataToString( style, tempInfo, crossSection, productFrame, multiplicity = multiplicity, weight = weight )
    s = versionStr + '\n' + sProcess + sSpecific + sCommon + sData
    return( executeCommand( logFile, transferMatrixExecute, s, workDir, tempInfo['workFile'] ) )

def discreteGammaAngularData( style, tempInfo, gammaEnergy, crossSection, angularData, multiplicity = 1, comment = None ) :
    """Currently, only isotropic (i.e., l = 0) data are returned. That is, lMax and angularData are ignored. This routine is also used
    for pair-production which pass angularData as None."""

    from fudge.processing import miscellaneous as miscellaneousModule

    reactionSuite = tempInfo['reactionSuite']
    projectileName = reactionSuite.projectile.name
    projectileGroupBoundaries = style.transportables[projectileName].group.boundaries.values
    productName = tempInfo['productName']
    productGroupBoundaries = style.transportables[productName].group.boundaries.values

    nProj = len( projectileGroupBoundaries ) - 1
    nProd = len( productGroupBoundaries ) - 1
    TM_1, TM_E = {}, {}
    for i1 in range( nProj ) :
        TM_1_i1, TM_E_i1 = {}, {}
        for i2 in range( nProd ) :
            TM_1_i1[i2] = [ 0. ]
            TM_E_i1[i2] = [ 0. ]
        TM_1[i1] = TM_1_i1
        TM_E[i1] = TM_E_i1

    indexEp = -1
    for Ep in productGroupBoundaries :
        if( gammaEnergy <= Ep ) : break
        indexEp += 1
    indexEp = min( max( indexEp, 0 ), nProd - 1 )
    form = style.findFormMatchingDerivedStyles( crossSection )
    f1 = form
    if( isinstance( form, crossSectionModule.reference ) ) : form = style.findFormMatchingDerivedStyles( form.link )
    if( isinstance( form, crossSectionModule.regions1d ) ) : form = form.toPointwise_withLinearXYs( 1e-8, 1e-8 )
    xsec = miscellaneousModule.groupOneFunctionAndFlux( style, tempInfo, form )
    for indexE in range( nProj ) :
        x = multiplicity * xsec[indexE]
        TM_1[indexE][indexEp] = [ x ]
        TM_E[indexE][indexEp] = [ x * gammaEnergy ]

    return( TM_1, TM_E )

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
    projectileName = reactionSuite.projectile.name
    projectileGroupBoundaries = style.transportables[projectileName].group.boundaries
    productName = tempInfo['productName']
    productGroupBoundaries = style.transportables[productName].group.boundaries
    flux = style.flux
    groupedFlux = tempInfo['groupedFlux']

    bindingEnergy = energyData.value.getValueAs( tempInfo['incidentEnergyUnit']  )
    massRatio = energyData.massRatio

    nProj = len( projectileGroupBoundaries.values ) - 1
    nProd = len( productGroupBoundaries.values ) - 1
    TM_1, TM_E = {}, {}
    for i1 in range( nProj ) :
        TM_1[i1] = {}
        TM_E[i1] = {}
        for i2 in range( nProd ) :
            TM_1[i1][i2] = [ 0. ]
            TM_E[i1][i2] = [ 0. ]
    Eg2 = bindingEnergy + massRatio * projectileGroupBoundaries.values[0]
    for indexEo, Eo in enumerate( productGroupBoundaries.values ) :
        if( Eg2 <= Eo ) : break
    indexEo = min( max( indexEo - 1, 0 ), nProd - 1 )

    crossSectionForm = style.findFormMatchingDerivedStyles( crossSection )
    crossSectionForm = crossSectionForm.toPointwise_withLinearXYs( 1e-8, 1e-8 )
    EMin, EMax = crossSectionForm.domainMin( ), crossSectionForm.domainMax( )
    axes = axesModule.axes( labelsUnits = { 0 : ( 'energy_out', tempInfo['incidentEnergyUnit'] ), 
                                            1 : ( crossSectionForm.axes[1].label, crossSectionForm.axes[1].unit ) } )
    Egp = XYsModule.XYs1d( data = [ [ EMin, bindingEnergy + massRatio * EMin ], [ EMax, bindingEnergy + massRatio * EMax ] ], 
            axes = axes, accuracy = 1e-12 )

    for indexEi in range( nProj ) :
        Ei2 = projectileGroupBoundaries.values[indexEi + 1]
        Eg2 = bindingEnergy + massRatio * Ei2
        EiMin = projectileGroupBoundaries.values[indexEi]
        while( True ) :
            incrementIndexEo, EiMax = 0, Ei2
            if( indexEo < ( nProd - 1 ) ) :
                if( Eg2 > productGroupBoundaries.values[indexEo + 1] ) :
                    incrementIndexEo = 1
                    EiMax = ( productGroupBoundaries.values[indexEo + 1] - bindingEnergy ) / massRatio
            TM_1[indexEi][indexEo][0] = float( crossSectionForm.integrateTwoFunctions( flux[0], 
                    domainMin = EiMin, domainMax = EiMax ) / groupedFlux[indexEi] )
            TM_E[indexEi][indexEo][0] = float( crossSectionForm.integrateThreeFunctions( flux[0], Egp, 
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
    if( comment is not None ) : s += "Comment: %s\n" % comment
    s += "Number of particles: %s\n" % numberOfProducts
    s += "Total mass: %s\n" % mTotal
    s += "Q value: %s\n" % Q
    s += 'Quadrature method: Gauss6\n'
    productFrame = standardsModule.frames.centerOfMassToken
    s += commonDataToString( style, tempInfo, crossSection, productFrame, multiplicity = multiplicity )
    return( executeCommand( logFile, transferMatrixExecute, s, workDir, tempInfo['workFile'] ) )

def KalbachMann_TransferMatrix( style, tempInfo, crossSection, particlesData, KalbachMannData, multiplicity = 1, comment = None ) :

    logFile = tempInfo['logFile'], 
    workDir = tempInfo['workDir']

    energy_in_unit = 'MeV'
    s  = versionStr + '\n'
    s += "Process: Kalbach spectrum\n"
    if( comment is not None ) : s += "Comment: %s\n" % comment
    s += "Projectile's ZA: %s\n" % particlesData['projectile']['ZA']
    s += "Target's ZA: %s\n" % particlesData['target']['ZA']
    s += "Product's ZA: %s\n" % particlesData['product']['ZA']
    s += "Compound's mass: %s\n" % particlesData['compound']['mass']
    s += "Quadrature method: adaptive\n"
    productFrame = standardsModule.frames.centerOfMassToken
    s += commonDataToString( style, tempInfo, crossSection, productFrame, multiplicity = multiplicity, energy_in_unit = energy_in_unit )
    s += KalbachMannDataToString( KalbachMannData, energy_in_unit )
    return( executeCommand( logFile, transferMatrixExecute, s, workDir, tempInfo['workFile'] ) )

def executeCommand( logFile, file, cmd, workDir, workFile ) :

    def checkNegative_l0( TM_EEpL, weight, f ) :

        negative_l0Counter, largestNegative = 0, 0.
        for ig in sorted( TM_EEpL.keys( ) ) :
            TM_EpL = TM_EEpL[ig]
            for ih in sorted( TM_EpL.keys( ) ) :
                l0Cell = TM_EpL[ih][0]
                if( l0Cell < 0 ) :
                    negative_l0Counter += 1
                    largestNegative = min( largestNegative, l0Cell )
                    fOut.write( 'negative l=0 for weight %s at row %3d column %3d: %s\n' % ( weight, ig, ih, l0Cell ) )
        return( negative_l0Counter )

    if( workFile == [] ) :
        dataFile = fudgeFileMisc.fudgeTempFile( dir = workDir, suffix = "_in" )
        fullFileName = dataFile.name
    else :
        workFile = '_'.join( workFile ) + '_in'
        if( workDir is not None ) :
            if( not( os.path.exists( workDir ) ) ) : os.makedirs( workDir )
            fullFileName = os.path.join( workDir, workFile )
        dataFile = open( fullFileName, 'w' )
    dataFile.write( cmd )
    dataFile.close( )
    infoFile = '%s.info' % fullFileName
    t0 = times.times( )
    try :
        status, stdout, stderr = subprocessing.executeCommand( [ file, '-output', fullFileName + '.out', fullFileName ], 
            stdout = infoFile, stderr = subprocess.STDOUT )
    except :
        fErr = open( fullFileName + ".err", "w" )
        fErr.close( )
        raise
    fOut = open( infoFile, 'a' )
    fOut.write( str( t0 ) + "\n" )
    if( status != 0 ) :
        print "status = ", status
        print s
        raise Exception( 'Transfer failed for %s %s' % ( file, fullFileName ) )
    TM1, TME = parseOutputFile( fullFileName + '.out' )
    negative_l0Counter = 0
    if( TM1 is not None ) : negative_l0Counter = checkNegative_l0( TM1, "0", fOut )
    if( TME is not None ) : negative_l0Counter += checkNegative_l0( TME, "E", fOut )
    if( negative_l0Counter > 0 ) : print 'WARNING: %d negative l=0 elements found in transfer matrix' % negative_l0Counter
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

def parseOutputFile( file ) :

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
            raise TypeError( 'bad value type at line %d' % ( value, lineNumber ) )

    def parseTransferMatrix( paras, lineNumber, ls, file ) :

        s1 = ls[lineNumber].split( ":" )                      # This line should start as 'Integrals, weight = ...'
        EBins = checkKeyValue( lineNumber, s1[1], 'numEinBins' )
        if( EBins > paras['numEinBins'] ) : 
            raise Exception( 'Bad transfer numEinBins = %d != %d in file %s\n%s' % ( EBins, paras['numEinBins'], file, ls[lineNumber] ) )
        lMax = paras['outputLegendreOrder']
        lMaxPlus1 = lMax + 1
        TM = {}
        for i1 in range( paras['numEinBins'] ) : 
            TMEp = {}
            for i2 in range( paras['numEoutBins'] ) : TMEp[i2] = lMaxPlus1 * [ 0. ]
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
                rowLs = map( float, ls[lineNumber].split( ) )
                if( len( rowLs ) > lMaxPlus1 ) : raise Exception( 'Bad transfer lRows at lineNumber = %d in file %s\n%s' % ( lineNumber, file, ls[lineNumber] ) )
                for i3, Cl in enumerate( rowLs ) : TM[i1][i2][i3] = Cl
        return( TM, lineNumber )

    TM1 = None
    TME = None
    try :
        f = open( file )
    except :
        raise Exception( 'Could not open transfer file = %s' % file )
    ls = f.readlines( )
    f.close( )
    if( ls[0] != "xndfgenGetTransfer: version 1\n" ) : raise Exception( 'Bad transfer file version in file %s\n%s' % ( file, ls[0] ) )
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
        elif( key == "Integrals, weight = 1" ) :
            if( TM1 is not None  ) : raise Exception( "Bad transfer file, multiple weight = 1 datasets in file %s" % file )
            TM1, lineNumber = parseTransferMatrix( paras, lineNumber, ls, file )
        elif( key == "Integrals, weight = E'" ) :
            if( TME is not None  ) : raise Exception( "Bad transfer file, multiple weight = E' datasets in file %s" % file )
            TME, lineNumber = parseTransferMatrix( paras, lineNumber, ls, file )
        else :
            raise Exception( 'Bad transfer file key in file %s\n%s' % ( file, ls[lineNumber] ) )
        lineNumber += 1
    return( TM1, TME )

def commonDataToString( style, tempInfo, crossSection, productFrame, multiplicity = None, energy_in_unit = None, weight = None ) :

    reactionSuite = tempInfo['reactionSuite']
    projectileName = reactionSuite.projectile.name
    projectileGroupBoundaries = style.transportables[projectileName].group
    productName = tempInfo['productName']
    productGroupBoundaries = style.transportables[productName].group
    flux = style.flux

    if( energy_in_unit is not None ) :
        flux = [ fluxOrder.convertAxisToUnit( 1, energy_in_unit ) for fluxOrder in flux ]

    s  = "outputLegendreOrder: %d\n" % style.lMax
    masses = tempInfo['masses']
    for particle, mass in masses.items( ) :
        if( mass is not None ) : s += "%s's mass: %s\n" % ( particle, floatToString( masses[particle] ) )
    s += "Projectile Frame: %s\n" % reactionSuite.projectileFrame
    s += GBToString( 'Projectile', projectileGroupBoundaries, energy_in_unit )
    maximumEnergy_in = projectileGroupBoundaries.boundaries.values[-1]
    s += GBToString( 'Product', productGroupBoundaries, energy_in_unit )
    s += fluxToString( flux )
    if( crossSection is not None ) : s += crossSectionToString( style, crossSection, energy_in_unit )
    if( multiplicity is not None ) : s += multiplicityToString( style, multiplicity, energy_in_unit = energy_in_unit )
    if( productName == 'gamma' ) : productFrame = standardsModule.frames.labToken     # BRB, hardwired??????
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

def fluxToString( flux ) :

    if( len( flux ) != 1 ) : raise Exception( 'FIXME' )
    interpolationStr, s = twoDToString( "", flux[0], addHeader = False )
    s  = [ "Fluxes: n = %d" % len( flux[0] ) ]
    s.append( interpolationStr )
    for Ein, Cls in flux[0] :
        s.append( 'Ein: ' + ( doubleFmt % Ein ) + ': n = 1' )
        s.append( '  ' + ( doubleFmt % Cls ) )
    return( startOfNewData + '\n'.join( s ) )

def crossSectionToString( style, crossSection, energy_in_unit = None ) :

    form = style.findFormMatchingDerivedStyles( crossSection )
    if( isinstance( form, crossSectionModule.reference ) ) : form = style.findFormMatchingDerivedStyles( form.link )
    if( isinstance( form, crossSectionModule.regions1d ) ) : form = form.toPointwise_withLinearXYs( 1e-8, 1e-8 )
    if( energy_in_unit is not None ) : form = form.convertAxisToUnit( 1, energy_in_unit )
    return( startOfNewData + '\n'.join( twoDToString( "Cross section", form ) ) )

def multiplicityToString( style, _multiplicity, energy_in_unit = None ) :

    form = style.findFormMatchingDerivedStyles( _multiplicity )
    multiplicity = form.toPointwise_withLinearXYs( 1e-8, 1e-8 )
    if( energy_in_unit is not None ) : multiplicity = multiplicity.convertAxisToUnit( 1, energy_in_unit )
    return( startOfNewData + '\n'.join( twoDToString( "Multiplicity", multiplicity ) ) )

def angularToString( angularData, crossSection, weight = None, twoBody = False ) :

    weightString = ''
    if( weight is not None ) : weightString = startOfNewData + '\n'.join( twoDToString( "weight", weight, addExtraBlankLine = False ) ) + '\n'

    if( isinstance( angularData, angularModule.isotropic ) ) :
        s = "Angular data: n = 2\n"
        s += 'Incident energy interpolation: %s %s\n' % ( linlin, GND2ProcessingInterpolationQualifiers[""] )
        s += 'Outgoing cosine interpolation: %s\n' % linlin
        for x in [ crossSection.domainMin( ), crossSection.domainMax( ) ] : s += " %s : n = 2\n   -1 0.5\n   1 0.5\n" % ( doubleFmt % x )

    elif( isinstance( angularData, angularModule.XYs2d ) ) :
        if( isinstance( angularData[0], angularModule.XYs1d ) ) :
            s = "Angular data: n = %s\n" % len( angularData )

            if( angularData.interpolationQualifier != '' ) :
                raise Exception( 'interpolation qualifier = %s is not supported' % angularData.interpolationQualifier )

            if( angularData.interpolation not in [ linlin, flat ] ) :
                _angularData = angularData.copy( )
                _angularData.interpolation = linlin
                print 'WARNING: interpolation %s is ignored, treated as %s' % ( angularData.interpolation, _angularData.interpolation )

            if( twoBody ) :
                s += 'Incident energy interpolation: %s %s\n' % \
                        ( angularData.interpolation, GND2ProcessingInterpolationQualifiers[angularData.interpolationQualifier] )
            else :
                s += 'Incident energy interpolation: %s %s\n' % ( angularData.interpolation,
                        GND2ProcessingInterpolationQualifiers[standardsModule.interpolation.unitBaseToken] ) # Ignoring data's interpolation qualifier.

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
            if( angularData.interpolationQualifier != '' ) :
                raise Exception( 'interpolation qualifier = %s is not supported' % angularData.interpolationQualifier )
            if( interpolation not in [ linlin, flat ] ) :
                _angularData = angularData.copy( )
                _angularData.interpolation = linlin
                print 'WARNING: ignoring interpolation "%s" and using "%s" instead' % ( angularData.interpolation, _angularData.interpolation )
                angularData = _angularData

            s = [ "Legendre coefficients: n = %s" % len( angularData ) ]

            if( twoBody ) :
                interpolationStr = 'Interpolation: %s' % angularData.interpolation
            else :
                interpolationStr = 'Interpolation: %s' % angularData.interpolation
            s.append( interpolationStr )

            for energy_in in angularData :
                s.append( ' Ein: %s:  n = %d' % ( energy_in.value, len( energy_in ) ) )
                for coefficient in energy_in : s.append( "%s" % coefficient )
            s = '\n'.join( s ) + '\n'
        else :
            raise Exception( "angular data = %s not supported" % angularData[0].moniker )

    else :
        brb.objectoutline( angularData )
        raise Exception( "angular data = %s not supported" % angularData.moniker )

    return( startOfNewData + weightString + s )

def energyFunctionToString( energyData, weight = None ) :

    def getParameter( data, label ) :

        _linlin = data.toPointwise_withLinearXYs( lowerEps = lowerEps, upperEps = upperEps )
        sData = [ startOfNewData, "%s: n = %d" % ( label, len( _linlin ) ) ]
        sData.append( 'Interpolation: %s' % linlin )
        for x, y in _linlin : sData.append( '%s %s' % ( x, y ) )
        return( sData )

    from fudge.gnd.productData.distributions import energy

    sData = []
    if( isinstance( energyData, ( energy.simpleMaxwellianFissionSpectrum, energy.evaporationSpectrum, energy.WattSpectrum ) ) ) :
        sData.append( 'U: %s' % energyData.U.getValueAs( energyData.parameter1.data.axes[-1].unit ) )

    if( isinstance( energyData.parameter1.data, regionsModule.regions1d ) ) :
        parameter1 = energyData.parameter1.data.toPointwise_withLinearXYs( 1e-5, lowerEps = lowerEps, upperEps = upperEps )
    else :
        parameter1 = energyData.parameter1.data.toPointwise_withLinearXYs( lowerEps = lowerEps, upperEps = upperEps )

    if( energyData.parameter2 is not None ) :
        if( isinstance( energyData.parameter2.data, regionsModule.regions1d ) ) :
            parameter2 = energyData.parameter2.data.toPointwise_withLinearXYs( 1e-5, lowerEps = lowerEps, upperEps = upperEps )
        else :
            parameter2 = energyData.parameter2.data.toPointwise_withLinearXYs( lowerEps = lowerEps, upperEps = upperEps )

    if( weight is None ) :
        axes = axesModule.axes( )
        weight = XYsModule.XYs1d( data = [ [ parameter1.domainMin( ), 1. ], [ parameter1.domainMax( ), 1. ] ], axes = axes, accuracy = 1e-6 )
    else :
        weight = weight.weight
    sData += twoDToString( "weight", weight, addExtraBlankLine = False )

    if( isinstance( energyData, energy.generalEvaporationSpectrum ) ) :
        sProcess = 'Process: general evaporation\n'
        sSpecific = ''
        sData += getParameter( parameter1, 'Theta' )
        sData += getParameter( parameter2, 'g' )
    elif( isinstance( energyData, energy.simpleMaxwellianFissionSpectrum ) ) :
        sProcess = 'Process: Maxwell spectrum\n'
        sSpecific = ''
        sData += getParameter( parameter1, 'Theta' )
    elif( isinstance( energyData, energy.evaporationSpectrum ) ) :
        sProcess = 'Process: Evaporation spectrum\n'
        sSpecific = 'Interpolate Eout integrals: false\n'
        sSpecific += 'Quadrature method: Gauss6\n'
        sData += getParameter( parameter1, 'Theta' )
    elif( isinstance( energyData, energy.WattSpectrum ) ) :
        sProcess = 'Process: Watt spectrum\n'
        sSpecific = 'Interpolate Eout integrals: false\n'
        sSpecific += 'Conserve: number\n'
        sData += getParameter( parameter1, 'a' )
        sData += getParameter( parameter2, 'b' )
    elif( isinstance( energyData, energy.MadlandNix ) ) :
        sProcess = 'Process: Madland-Nix spectrum\n'
        sSpecific  = 'Quadrature method: adaptive\n'
        sSpecific += 'Interpolate Eout integrals: false\n'
        sSpecific += 'Conserve: number\n'
        sSpecific += 'EFL: %s\n' % energyData.EFL.getValueAs( energyData.parameter1.data.axes[-1].unit )
        sSpecific += 'EFH: %s\n' % energyData.EFH.getValueAs( energyData.parameter1.data.axes[-1].unit )
        sData += getParameter( parameter1, 'TM' )
    else :
        raise Exception( 'Unsupport energy functional = %s' % energyData.moniker )
    sData.append( '' )
    return( sProcess, sSpecific, startOfNewData + '\n'.join( sData ) )

def EEpPDataToString( EEpPData ) :

    from fudge.gnd.productData.distributions import energy as energyModule

    s  = startOfNewData
    s += "EEpPData: n = %d\n" % len( EEpPData )

    qualifier = EEpPData.interpolationQualifier
    if( qualifier == standardsModule.interpolation.noneQualifierToken ) : qualifier = standardsModule.interpolation.unitBaseToken
    s += 'Incident energy interpolation: %s %s\n' % ( EEpPData.interpolation, GND2ProcessingInterpolationQualifiers[qualifier] )

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
    s += 'Incident energy interpolation: %s %s\n' % (
            linlin, GND2ProcessingInterpolationQualifiers[standardsModule.interpolation.unitBaseToken] )
    s += 'Outgoing energy interpolation: %s\n' % linlin
    for l_EEpP in LEEpPData :
        s += startOfNewSubData + '\n'
        EEpP = l_EEpP.EpP
        s += "  l = %d: n = %d\n" % ( l_EEpP.l, len( EEpP ) )
        s += threeDToString( EEpP )
    return( s )

def LegendreDataToString( LegendreData ) :

    s = [ 'Legendre data by incident energy:  n = %d' % len( LegendreData ) ]


    qualifier = LegendreData.interpolationQualifier
    if( qualifier == standardsModule.interpolation.noneQualifierToken ) : qualifier = standardsModule.interpolation.unitBaseToken
    s.append( 'Incident energy interpolation: %s %s' % ( LegendreData.interpolation, GND2ProcessingInterpolationQualifiers[qualifier] ) )

    interpolation = LegendreData[0].interpolation
    if( interpolation in [ linlin, flat ] ) :
        s.append( 'Outgoing energy interpolation: %s' % interpolation )
    else :
        print '    WARNING: converting interpolation from "%s" to "%s".' % ( interpolation, linlin )
        s.append( 'Outgoing energy interpolation: %s' % linlin )

    for EEpCl in LegendreData :
        s.append( " Ein: %s:  n = %d" % ( EEpCl.value, len( EEpCl ) ) )
        for EpCl in EEpCl :
            s.append( startOfNewSubData )
            s.append( "  Eout: %s:  n = %d" % ( EpCl.value, len( EpCl ) ) )
            for Cl in EpCl : s.append( "   %s" % Cl )
    s.append( '' )
    return( startOfNewData + '\n'.join( s ) )

def EMuEpPDataToString( EMuEpPData ) :

    s  = startOfNewData
    s += "EMuEpPData: n = %d\n" % len( EMuEpPData )

    qualifier = EMuEpPData.interpolationQualifier
    if( qualifier not in [ standardsModule.interpolation.noneQualifierToken, standardsModule.interpolation.unitBaseToken ] ) :
        raise Exception( 'interpolation qualifier = %s is not supported' % qualifier )
    if( EMuEpPData.interpolation != linlin ) :
        raise Exception( 'interpolation = "%s" is not supported' % EMuEpPData.interpolation )
    s += 'Incident energy interpolation: %s %s\n' % ( linlin,
            GND2ProcessingInterpolationQualifiers[standardsModule.interpolation.unitBaseToken] )

    s += 'Outgoing cosine interpolation: %s %s\n' % ( linlin, 
            GND2ProcessingInterpolationQualifiers[standardsModule.interpolation.unitBaseToken] )

    s += 'Outgoing energy interpolation: %s\n' % linlin

    for muEpP in EMuEpPData :
        s += startOfNewSubData + '\n'
        E = muEpP.value
        s += "  E = %s: n = %d\n" % ( doubleFmt % E, len( muEpP ) )
        s += threeDToString( muEpP )
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
        _xys = XYsModule.XYs1d( data = xys.scaleOffsetXAndY( xScale = factor, yScale = 1 / factor ), value = factor * xys.value )
        _fSubform.append( _xys )

    _rSubform = multiD_XYsModule.XYs2d( )
    for xys in rSubform :
        _xys = XYsModule.XYs1d( data = xys.scaleOffsetXAndY( xScale = factor ), value = factor * xys.value )
        _rSubform.append( _xys )

    s  = startOfNewData
    s += "Kalbach probabilities: n = %s\n" % len( _fSubform )
    s += "Incident energy interpolation: %s %s\n" % ( linlin,
            GND2ProcessingInterpolationQualifiers[standardsModule.interpolation.unitBaseToken] )
    s += "Outgoing energy interpolation: %s\n" % fSubform[0].interpolation
    s += threeDListToString( _fSubform )

    s += "Kalbach r parameter: n = %s\n" % len( rSubform )
    s += "Incident energy interpolation: %s %s\n" % ( rSubform.interpolation, 
            GND2ProcessingInterpolationQualifiers[standardsModule.interpolation.unitBaseUnscaledToken] )
    s += "Outgoing energy interpolation: %s\n" % rSubform[0].interpolation
    s += threeDListToString( _rSubform )

    return( s )

def floatToString( f ) :

    if( f == 0 ) :
        s = "0.0"
    elif( ( f > 1000. ) or ( f < .01 ) ) :
        s = "%19.12e" % f
    else :
        s = "%.15f" % f
        while( ( len( s ) > 2 ) and ( s[-1] == '0' ) ) :
            s = s[:-1]
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

    if( isinstance( data, XYsModule.XYs1d ) ) :
        if( data.interpolation not in [ linlin, flat ] ) :
            data = data.toPointwise_withLinearXYs( lowerEps = lowerEps, upperEps = upperEps )
        interpolationStr = 'Interpolation: %s' % data.interpolation
    else :
        print '======== HI =========', __file__
        brb.objectoutline( data )
        raise 'hell'
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
            raise Exception( 'need to update for new gnd forms' )
            x, yz = xyz
        else :
            w = w_xy.value
            xy = w_xy
            if( linearizeXYs ) : xy = w_xy.toPointwise_withLinearXYs( accuracy = 1e-3, lowerEps = 1e-8, upperEps = 1e-8 )
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
        w = w_xy.value
        if( linearizeXYs ) : w_xy = w_xy.toPointwise_withLinearXYs( accuracy = 1e-3, lowerEps = 1e-8, upperEps = 1e-8 )
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
