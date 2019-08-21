# <<BEGIN-copyright>>
# <<END-copyright>>

import sys
import math

from fudge import gnd
from fudge.core import fudgemisc
from fudge.core.utilities import brb, fudgeZA
from pqu import physicalQuantityWithUncertainty
from fudge.core.math import fudgemath, matrix as gndMatrix
from fudge.core.math.xData import axes, XYs, LegendreSeries, regions, W_XYs

from fudge.processing import processingInfo
from fudge.gnd.productData import distributions

import endfFormats
import endf_endl
import endfFileToGNDMisc
import endlToGND

ENDFInterpolationToGND = { 1 : axes.flatToken, 2 : axes.linearToken }
frames = { 1 : axes.labToken, 2 : axes.centerOfMassToken }
ENDF_Accuracy = endfFileToGNDMisc.ENDF_Accuracy
FUDGE_EPS = endfFileToGNDMisc.FUDGE_EPS
productNameToZA = { 'n' : 1, 'H1' : 1001, 'H2' : 1002, 'H3' : 1003, 'He3' : 2003, 'He4' : 2004, 'gamma' : 0 }
lightIsotopeNames = [ 'n', 'H1', 'H2', 'H3', 'He3', 'He4' ]

class logFiles :

    def __init__( self, toStdOut = False, toStdErr = False, logFile = None, defaultIsStderrWriting = True ) :

        def addIfDifferent( newOutput ) :

            if( newOutput.isatty( ) ) :
                for output in self.outputs :
                    if( output.isatty( ) ) : return( False )
            return( True )

        self.defaultIsStderrWriting = defaultIsStderrWriting
        self.outputs = []
        if( toStdOut ) : self.outputs.append( sys.stdout )
        if( toStdErr and addIfDifferent( sys.stderr ) ) : self.outputs.append( sys.stderr )
        if( logFile ) :
            if( type( logFile ) == type( '' ) ) :
                logFile = open( logFile, 'w' )
            elif( type( logFile ) != type( sys.stdout ) ) :
                raise Exception( 'logFile type = %s is no supported' % ( type( logFlie ) ) )
            if( addIfDifferent( logFile ) ) : self.outputs.append( logFile )

    def write( self, msg, stderrWriting = None ) :

        if( stderrWriting is None ) : stderrWriting = self.defaultIsStderrWriting
        for output in self.outputs :
            if( ( output == sys.stderr ) and ( not( stderrWriting ) ) ) : continue
            output.write( msg )

# iterator that keeps track of line number:
class myIter:
    def __init__(self, iterable):
        self.index = 0
        self.iterable = iter(iterable)
        self.length = len(iterable)
    def next(self):
        next = self.iterable.next()
        self.index += 1
        return next

# two useful shortcuts for reading ENDF data:
funkyF = endfFileToGNDMisc.sixFunkyFloatStringsToFloats
def funkyFI( a, logFile = sys.stderr ) :   # read ENDF line with 2 floats and 4 ints

    return( endfFileToGNDMisc.sixFunkyFloatStringsToIntsAndFloats( a, [ 2, 3, 4, 5 ], logFile = logFile ) )

# create some custom Exceptions:
class BadResonances( Exception ): pass
class BadCovariance( Exception ): pass
    
def calculateZA( ZACompound, ZAOther, minus = True ) :
    """This function handles the removal (or addition) of ZAOther to ZACompound include natural compound (but not a natural other)."""

    if( ( ZACompound % 1000 ) == 0 ) : ZAOther = 1000 * ( ZAOther // 1000 )
    if( minus ) : return( ZACompound - ZAOther )
    return( ZACompound + ZAOther )

def getTypeNameGamma( info, ZA, level = None, levelIndex = None ) :

    if( not( level is None ) ) :
        if( fudgemath.isNumber( level ) ) :
            if( level < 0 ) :
                if( ( level > -100 ) and ( levelIndex == 0 ) ) :
                    level = 0
                else :
                    message = 'Negative excitation level = %s for ZA = %s and levelIndex = %s is not allowed' % ( level, ZA, levelIndex )
                    sys.stderr.write( message+'\n' )
                    info.doRaise.append( message )
                    level = 0
        elif( type( level ) != type( '' ) ) :
            raise Exception( 'level (type ="%s") must be a number or a string' % ( brb.getType( level ), ZA ) )
    if( ZA == 0 ) : ZA = 7                             # Special case for gammas which are yo = 7 for ENDL.
    p =  endlToGND.getTypeName( info, ZA, level = level, levelIndex = levelIndex, levelUnit = 'eV' )
    return( p )

def getTypeNameENDF( info, ZA, undefinedLevelInfo ) :

    levelIndex, level = None, 0.
    if( not( undefinedLevelInfo is None ) ) :
        if( not( undefinedLevelInfo['ZA'] is None ) ) :
            if( ( undefinedLevelInfo['ZA'] == ZA ) and ( ZA not in [ 1001, 1002, 1003, 2003, 2004 ] ) ) :
                undefinedLevelInfo['count'] += 1
                if( undefinedLevelInfo['count'] > 1 ) : raise Exception( "undefinedLevelInfo['count'] > 1 for ZA = %s" % ZA )
                levelIndex, level = undefinedLevelInfo['levelIndex'], undefinedLevelInfo['level']
    return( getTypeNameGamma( info, ZA, level = level, levelIndex = levelIndex ) )

def returnConstantQ( Q_eV ) :

    return( gnd.channelData.Q.component( gnd.channelData.Q.constant( physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( Q_eV, 'eV' ) ) ) )

def countZAMasses( info, ZA, AWR ) :

    if( ZA != 0 ) : 
        if( ZA not in info.ZA_AWRMasses ) : info.ZA_AWRMasses[ZA] = {}
        ZAData = info.ZA_AWRMasses[ZA]
        if( AWR not in ZAData ) : ZAData[AWR] = 0
        ZAData[AWR] += 1

def smearXYs( data, eps = FUDGE_EPS ) :
    """data must be an python array of [ [ x0, y0 ], [ x1, y1 ], ...  [ xn, yn ] ]."""

    if( eps == 0. ) : raise Exception( 'Need non-zeor eps' )
    if( eps > 0 ) :
        xPrior = None
        for xy in data :
            x = xy[0]
            if( xPrior is not None ) :
                if( x < xPrior ) : raise Exception( 'x = %s < xPrior = %s' % ( x, x < xPrior ) )
                if( x == xPrior ) :
                    x *= ( 1 + FUDGE_EPS )
                    xy[0] = x
            xPrior = x
    else :
        raise Exception( 'eps = %s < 0 is not implemented' % eps )

def getCrossSectionLinearOrPointwise( axes_, data ) :

    cls = gnd.reactionData.crossSection.pointwise
    independent, dependent, qualifier = axes_[0].interpolation.getInterpolationTokens( )
    if( ( axes.linearToken, axes.linearToken, None ) == ( independent, dependent, qualifier ) ) :
        cls = gnd.reactionData.crossSection.linear
    return( cls( axes_, data = data, accuracy = ENDF_Accuracy ) )

def getMultiplicityPointwiseOrPieceWise( data, warningList, yUnit = 'multiplicity' ) :

    energyInterpolation, multiplicityInterpolation, dummy = data[0].axes[0].interpolation.getInterpolationTokens( )
    axes_ = gnd.productData.multiplicity.pointwise.defaultAxes( energyInterpolation = energyInterpolation, multiplicityName = yUnit, \
        multiplicityInterpolation = multiplicityInterpolation )
    if( len( data  ) == 1 ) :
        data = data[0]
        xPrior = -1 * data[0][0] - 1                        # Make sure its less than first point
        doWarning = True
        for x, y in data :
            if( xPrior == x ) :
                if( doWarning ) : warningList.append( '       WARNING: same multiplicity energies, second one being incremented' )
                doWarning = False
                data[i] = [ ( 1. + FUDGE_EPS ) * x, y ]
            xPrior = x
        multiplicity = gnd.productData.multiplicity.pointwise( axes_, data, ENDF_Accuracy )
    else :
        doKludge = False
        if( ( len( data ) == 2 ) and ( data[0][-1] == data[1][0] ) ) :
            if( data[0].axes[0].interpolation == data[1].axes[0].interpolation ) : doKludge = True
        if( not( doKludge ) ) : raise Exception( 'piecewise multiplicity is not supported' )
        xys = data[1].copyDataToXYs( )
        xys[0][0] *= ( 1 + FUDGE_EPS )
        xys = data[0].copyDataToXYs( ) + xys
        multiplicity = gnd.productData.multiplicity.pointwise( axes_, xys, ENDF_Accuracy )
    return( multiplicity )

def getTotalOrPromptFission( info, MT456Data, totalOrPrompt, warningList ) :

    ZA, AWR, dummy, LNU, dummy, dummy = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( MT456Data[0], logFile = info.logs )
    countZAMasses( info, ZA, AWR )
    ZA = int( ZA )
    LNU = int( LNU )
    info.logs.write( '     %s fission neutron data: LNU = %d\n' % ( totalOrPrompt, LNU ) )
    if( LNU == 1 ) :
        dataLine, polynomial = endfFileToGNDMisc.getList( 1, MT456Data, logFile = info.logs )
        axes_ = gnd.productData.multiplicity.polynomial.defaultAxes( )
        fissionMultiplicity = gnd.productData.multiplicity.polynomial( axes_, polynomial['data'] )
    else :
        axes_ = gnd.productData.multiplicity.pointwise.defaultAxes( )
        dataLine, TAB1, multiplicityRegions = endfFileToGNDMisc.getTAB1Regions( 1, MT456Data, axes_, logFile = info.logs )
        if( len( multiplicityRegions ) > 1 ) : raise Exception( "Currently only one interpolation flag is supported" )
        fissionMultiplicity = getMultiplicityPointwiseOrPieceWise( multiplicityRegions, warningList )
    return( fissionMultiplicity )

def getDelayedFission( info, MT455Data, warningList ) :

    info.logs.write( '     Delayed fission neutron data (MT=455)' )
    MT455DataMF1 = MT455Data[1]
    ZA, AWR, LDG, LNU, dummy, dummy = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( MT455DataMF1[0], logFile = info.logs )
    countZAMasses( info, ZA, AWR )
    LDG, LNU = int( LDG ), int( LNU )
    info.logs.write( ' LDG=%s LNU=%s' % ( LDG, LNU ) )
    if( LDG != 0 ) : raise Exception( "Only energy-independent delayed fission neutrons are supported" )
    if( LNU != 2 ) : raise Exception( "Only tables of delayed fission neutrons are supported" )

    dataLine, decayRateData = endfFileToGNDMisc.getList( 1, MT455DataMF1, logFile = info.logs )
    NNF = int( decayRateData['NPL'] )
    decayRates = decayRateData['data']

    axes_ = gnd.productData.multiplicity.pointwise.defaultAxes( )
    dataLine, TAB1, multiplicityRegions = endfFileToGNDMisc.getTAB1Regions( dataLine, MT455DataMF1, axes_, logFile = info.logs )
    if( len( multiplicityRegions ) > 1 ) : raise Exception( "Currently only one interpolation flag is supported" )
    energyInterpolation, multiplicityInterpolation, dummy = multiplicityRegions[0].axes[0].interpolation.getInterpolationTokens( )

    totolDelayedMultiplicity = multiplicityRegions[0]
    if( 5 in MT455Data ) :
        delayedNeutronEnergies, weights = readMF5( info, 455, MT455Data[5], warningList, delayNeutrons = True )
        weightsSum = XYs.XYs( weights[0].axes, data = [], accuracy = ENDF_Accuracy )        # Renormalize weights to sum to 1.
        for weight in weights : weightsSum = weightsSum + weight
        for weight in weights : 
            for i, xy in enumerate( weight ) : XYs.pointwiseXY.__setitem__( weight, i, xy[1] / weightsSum.getValue( xy[0] ) )
    else :
        delayedNeutronEnergies = len( decayRates ) * [ None ]
        totolDelayedMultiplicity_ = totolDelayedMultiplicity / len( decayRates )
        nuBar = getMultiplicityPointwiseOrPieceWise( [ totolDelayedMultiplicity_ ], warningList )
        if( len( decayRates ) > 1 ) : warningList.append( '       WARNING: More than one delayed fission neutron decay time but no MF = 5 data' )
    decayChannel = gnd.channels.NBodyDecayChannel( returnConstantQ( 0. ) )     # Dummy channel, so Q's value does not matter.
    for i, decayRate in enumerate( decayRates ) :
        energyComponent = delayedNeutronEnergies[i]
        if( energyComponent is not None ) :
            weights_energyInterpolation, weights_multiplicityInterpolation, dummy = weights[i].axes[0].interpolation.getInterpolationTokens( )
            if( weights_energyInterpolation != energyInterpolation ) : raise \
                Exception( 'Total and weight energy interpolations differ: %s vs %s' % ( weights_energyInterpolation, energyInterpolation ) )
            if( weights_multiplicityInterpolation != multiplicityInterpolation ) : raise \
                Exception( 'Total and weight multiplicity interpolations differ: %s vs %s' % ( weights_multiplicityInterpolation, multiplicityInterpolation ) )
            if( weights_energyInterpolation != axes.linearToken ) : 
                raise Exception( 'For energy only linear interpolations supported: %s' % weights_energyInterpolation )
            if( weights_multiplicityInterpolation != axes.linearToken ) : 
                raise Exception( 'For multiplicity only linear interpolations supported: %s' % weights_multiplicityInterpolation )
            totolDelayedM, weight = totolDelayedMultiplicity, weights[i]
            if( totolDelayedMultiplicity.domainMax( ) > weight.domainMax( ) ) : totolDelayedM = totolDelayedMultiplicity.xSlice( xMax = weight.domainMax( ) )
            if( weight.domainMax( ) > totolDelayedM.domainMax( ) ) : weight = weight.xSlice( xMax = totolDelayedM.domainMax( ) )
            if( weight.domainMin( ) < totolDelayedM.domainMin( ) ) : 
                if( len( weight ) == 2 ) :
                    if( weight[0][1] == weight[1][1] ) : weight = weight.xSlice( xMin = totolDelayedM.domainMin( ) )
            multiplicity = totolDelayedM * weight
            nuBar = getMultiplicityPointwiseOrPieceWise( [ multiplicity ], warningList )
        for xy in nuBar:
            m,e = math.frexp(xy[-1])
            nuBar.setValue( xy[0], round(xy[-1], 10-int(e/3.321928095)) )
            #nuBar.setValue( xy[0], round(xy[-1], 10+int(abs(math.log10(xy[-1])))) )
        particle = endlToGND.newGNDParticle( info, getTypeNameGamma( info, 1 ), multiplicity = nuBar )
        particle.addAttribute( 'emissionMode', 'delayed' )
        particle.addAttribute( 'decayRate', physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty(decayRate, "1/s") )
        if( not( energyComponent is None ) ) :
            form = distributions.angular.isotropic( frames[1] )     # MF = 5 data is always in lab frame.
            angularComponent = distributions.angular.component( form.moniker )
            angularComponent.addForm( form )
            component = distributions.uncorrelated.component( angularComponent, energyComponent )
            angularComponent.setParent( component ); energyComponent.setParent( component )
            component.parent = particle.distributions
            particle.addDistributionComponent( component )
            particle.distributions.setNativeData( component.moniker )
        decayChannel.addProduct( particle )
    info.logs.write( '\n' )
    return( decayChannel )

def getFissionEnergies( info, MF458Data ) :
    """
    For NPLY = 0 this data consists of pairs ( energy, standard deviation ).  The energies are:
    EFR: kinetic energy of all fission fragments
    ENP: energy of prompt fission neutrons
    END: energy of delayed fission neutrons
    EGP: energy of prompt fission gammas
    EGD: energy of delayed fission gammas
    EB: energy of delayed fission betas
    ENU: energy of fission neutrinos
    ER: EFR + ENP + END + EGP + EGD + EB, fission Q minus the neutrinos
    ET: ER + ENU, total fission Q
    For NPLY > 0 the fission energies are polynomials of degree NPLY in the incident energy,
    and the structure of the above table is repeated, once for each polynomial coefficient.
    """

    MF1Data = MF458Data[1]
    dataLine = 0
    ZA, AWR, dummy, dummy, dummy, dummy = endfFileToGNDMisc.sixFunkyFloatStringsToIntsAndFloats( MF1Data[dataLine], intIndices = [ 0 ], logFile = info.logs )
    countZAMasses( info, ZA, AWR )

    dataLine += 1
    dummy, dummy, dummy, NPLY, lengthData, numEnergies = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( MF1Data[ dataLine ], logFile = info.logs )
    NPLY = int( NPLY )               # order of the polynomial representation of energy produced
    lengthData = int( lengthData )   # how much data
    numEnergies = int( numEnergies ) # number of fission energies (9)

    dataLine += 1
    energies = endfFileToGNDMisc.nFunkyFloatStringsToFloats( numEnergies, dataLine, MF1Data, dimension = 2, logFile = info.logs )

    labels = gnd.channelData.fissionEnergyReleased.polynomial.labels
    fissionEnergyData = {}
    for label in labels : fissionEnergyData[label] = []
    for i in xrange( 0, numEnergies, 9 ) :
        for j, label in enumerate( labels ) : fissionEnergyData[label].append( energies[i + j] )

    fissionEnergies = gnd.channelData.fissionEnergyReleased.component( )
    fissionEnergies.addFormAsNativeData( gnd.channelData.fissionEnergyReleased.polynomial( NPLY, fissionEnergyData, 'eV', hasUncertainties = True ) )

    return fissionEnergies

def angularLegendrePiecewiseToPointwiseIfPossible( piecewiseForm ) :

    if( len( piecewiseForm ) != 1 ) : return( piecewiseForm )
    region = piecewiseForm[0]
    energyInterpolation, energyFunctionInterpolation, dummy = region.axes[0].interpolation.getInterpolationTokens( )
    frame = piecewiseForm.axes[1].frame
    axes_ = distributions.angular.LegendrePointwise.defaultAxes( frame, energyInterpolation, energyFunctionInterpolation )
    pointwiseForm = distributions.angular.LegendrePointwise( axes_ )

    for i, energy in enumerate( region ) :
        pointwiseForm[i] = LegendreSeries.XYs_LegendreSeries( axes_[1].getUnit( ), energy.coefficients, index = i, value = energy.value )
    return( pointwiseForm )

def angularLegendreToPointwiseOrPiecewiseLegendre( MT, frame, angularData, warningList, MF, msg ) :

    axes_ = distributions.angular.LegendrePiecewise.defaultAxes( frame )
    formLegendre = distributions.angular.LegendrePiecewise( axes_ )
    index, start, lastRegion = 0, 0, None
    lists = angularData['Lists']
    for end, interpolationFlag in angularData['interpolationInfo'] :
        if( interpolationFlag > 5 ) : raise Exception( 'Unsupported interpolation = %s for MF=%s, MT=%s' % ( interpolationFlag, MF, MT ) )
        interpolationE_in, interpolationCl, interpolationQualifier = endfFileToGNDMisc.ENDFInterpolationToGNDAxes3plusd( interpolationFlag )
        regionsAxes = axes.interpolationAxes( 0, axes.interpolationXY( interpolationE_in, interpolationCl ), None )
        region = LegendreSeries.W_XYs_LegendreSeries( regionsAxes, index = index, parent = formLegendre )
        indexEnergyIn, priorEnergy = 0, -1
        if( not( lastRegion is None ) ) :
            if( lists[start]['C2'] != lastRegion[0] ) :
                region[indexEnergyIn] = LegendreSeries.XYs_LegendreSeries( axes_[1].getUnit( ), lastRegion[1], index = indexEnergyIn, value = lastRegion[0] )
                indexEnergyIn += 1
        for i in xrange( start, end ) :
            energy, coefficients = lists[i]['C2'], [ 1.0 ] + lists[i]['data']
            if( energy == priorEnergy ) :                           # This fixes a problem with some data having two same energy values.
                energy += FUDGE_EPS * energy
                warningList.append( '       WARNING: same energies, second one being incremented for MT = %d, MF = %d, %s' % ( MT, MF, msg ) )
            priorEnergy = energy
            region[indexEnergyIn] = LegendreSeries.XYs_LegendreSeries( axes_[1].getUnit( ), coefficients, index = indexEnergyIn, value = energy )
            indexEnergyIn += 1
        formLegendre[index] = region
        lastRegion = [ energy, coefficients ]
        index += 1
        start = end
    return( angularLegendrePiecewiseToPointwiseIfPossible( formLegendre ) )

def convertNuclearPlusInterferenceDataToPiecewise( MT, frame, angularData, warningList, MF, msg, identicalParticles ) :
    """Return distributions.LegendreNuclearPlusCoulombInterference instance. This in turn contains
    Legendre sections for both the nuclear and interference contributions."""

    axes_ = distributions.angular.LegendrePiecewise.defaultAxes( frame )
    nuclear = distributions.angular.LegendrePiecewise( axes_ )
    interferenceReal = distributions.angular.LegendrePiecewise( axes_ )
    interferenceImaginary = distributions.angular.LegendrePiecewise( axes_ )
    index, start, lastRegion = 0, 0, None
    lists = angularData['Lists']
    for end, interpolationFlag in angularData['interpolationInfo'] :
        interpolationE_in, interpolationCl, interpolationQualifier = endfFileToGNDMisc.ENDFInterpolationToGNDAxes3plusd( interpolationFlag )
        if( interpolationFlag > 6 ) : raise Exception( 'Unsupported interpolation = %s for MF=%s, MT=%s' % ( interpolationFlag, MF, MT ) )
        interpolationE_in, interpolationCl, interpolationQualifier = endfFileToGNDMisc.ENDFInterpolationToGNDAxes3plusd( interpolationFlag )
        regionsAxes = axes.interpolationAxes( 0, axes.interpolationXY( interpolationE_in, interpolationCl ), None )
        region_Nuc = LegendreSeries.W_XYs_LegendreSeries( regionsAxes, index = index, parent = nuclear )
        region_IntReal = LegendreSeries.W_XYs_LegendreSeries( regionsAxes, index = index, parent = interferenceReal )
        region_IntImaginary = LegendreSeries.W_XYs_LegendreSeries( regionsAxes, index = index, parent = interferenceImaginary )
        indexEnergyIn, priorEnergy = 0, -1
        if( not( lastRegion is None ) ) :
            if( lists[start]['C2'] != lastRegion ) : # ensure no gaps between regions
                region_Nuc[indexEnergyIn] = LegendreSeries.XYs_LegendreSeries( axes_[1].getUnit( ), nuclear_term,      index = indexEnergyIn, value = lastRegion )
                region_IntReal[indexEnergyIn] = LegendreSeries.XYs_LegendreSeries( axes_[1].getUnit( ), interference_termReal, index = indexEnergyIn, value = lastRegion )
                region_IntImaginary[indexEnergyIn] = LegendreSeries.XYs_LegendreSeries( axes_[1].getUnit( ), interference_termImaginary, index = indexEnergyIn, value = lastRegion )
                indexEnergyIn += 1
        for i in xrange( start, end ) :
            list = lists[i]
            energy = list['C2']
            if( energy == priorEnergy ) :                           # This fixes a problem with some data having two same energy values.
                energy += FUDGE_EPS * energy
                warningList.append( '       WARNING: same energies, second one being incremented for MT = %d, MF = %d, %s' % ( MT, MF, msg ) )
            priorEnergy = energy
            if( identicalParticles ) :
                splitPoint = list['N2'] + 1
            else :
                splitPoint = list['N2'] * 2 + 1
            nuclear_term = list['data'][:splitPoint]
            interference_term, interference_termReal, interference_termImaginary = list['data'][splitPoint:], [], []
            for i in range( 0, len( interference_term ), 2 ) :
                interference_termReal.append( interference_term[i] )
                interference_termImaginary.append( interference_term[i+1] )
            region_Nuc[indexEnergyIn] = LegendreSeries.XYs_LegendreSeries( axes_[1].getUnit( ), nuclear_term,      index = indexEnergyIn, value = energy )
            region_IntReal[indexEnergyIn] = LegendreSeries.XYs_LegendreSeries( axes_[1].getUnit( ), interference_termReal, index = indexEnergyIn, value = energy )
            region_IntImaginary[indexEnergyIn] = LegendreSeries.XYs_LegendreSeries( axes_[1].getUnit( ), interference_termImaginary, index = indexEnergyIn, value = energy )
            indexEnergyIn += 1
        lastRegion = energy
        nuclear[index] = region_Nuc
        interferenceReal[index] = region_IntReal
        interferenceImaginary[index] = region_IntImaginary
        index += 1
        start = end
    nuclear = angularLegendrePiecewiseToPointwiseIfPossible( nuclear )
    interferenceReal = angularLegendrePiecewiseToPointwiseIfPossible( interferenceReal )
    interferenceImaginary = angularLegendrePiecewiseToPointwiseIfPossible( interferenceImaginary )
    # nuclear and interference are both angularLegendrePiecewise objects, but need to specify which is which:
    nuclear.moniker = distributions.base.LegendrePiecewise_NuclearFormToken
    interferenceReal.moniker = distributions.base.LegendrePiecewise_CoulombInterferenceRealFormToken
    interferenceImaginary.moniker = distributions.base.LegendrePiecewise_CoulombInterferenceImaginaryFormToken
    return( distributions.angular.nuclearPlusCoulombInterference( nuclear, interferenceReal, interferenceImaginary ) )

def convertAngularToPointwiseOrPiecewise( MT, LANG, frame, angularLOrT, warningList, allowDuplicateEnergiesInMF4Table = False ) :

    interpolationE_in = endfFileToGNDMisc.ENDFInterpolationToGND3plusd( angularLOrT['interpolationInfo'][0][1] )
    if( LANG > 0 ) :
        interpolationFlagMin = LANG - 10
        interpolationFlagMax = interpolationFlagMin
    else :
        interpolationFlagMin = interpolationFlagMax = angularLOrT['TAB1s'][0]['interpolationInfo'][0][1]
        for data in angularLOrT['TAB1s'] :
            if( data['NR'] != 1 ) : raise Exception( "Currently only one interpolation flag is supported for MT = %d, MF = 4, LTT = 2" % MT )
            interpolationFlagMin = min( interpolationFlagMin, data['interpolationInfo'][0][1] )
            interpolationFlagMax = min( interpolationFlagMax, data['interpolationInfo'][0][1] )

    if( interpolationFlagMin != interpolationFlagMax ) :
        warningList.append( '       WARNING: Some interpolation flags changing %s for MF = 4, LTT = 2' % interpolationFlagMin )
    interpolationx, interpolationy = endfFileToGNDMisc.ENDFInterpolationToGND2d( interpolationFlagMin )
    axes_ = axes.axes( dimension = 3 )
    axes_[0] = axes.axis( 'energy_in', 0, 'eV', frame = axes.labToken, interpolation = axes.interpolationXY( axes.byRegionToken, axes.byRegionToken ) )
    axes_[1] = axes.axis( 'mu', 1, '', frame = frame, interpolation = axes.interpolationXY( axes.byRegionToken, axes.byRegionToken ) )
    axes_[2] = axes.axis( 'P(energy_in|mu)', 2, '', frame = frame )
    piecewiseForm = distributions.angular.piecewise( axes_, IKnowWhatIAmDoing = True )
    index, start, lastRegion = 0, 0, None
    if( LANG > 0 ) :
        lOrTs = angularLOrT['Lists']
    else :
        lOrTs = angularLOrT['TAB1s']
    doInterpolationChangeWarning = True
    for i, endInterpolationFlag in enumerate( angularLOrT['interpolationInfo'] ) :
        end, interpolationFlag = endInterpolationFlag
        if( LANG > 0 ) :
            mu_pdf_Interpolation = LANG - 10
        else :
            mu_pdf_Interpolation = lOrTs[i]['interpolationInfo'][0][1]
            for data in lOrTs : mu_pdf_Interpolation = min( mu_pdf_Interpolation, data['interpolationInfo'][0][1] )
            for data in lOrTs :
                if( mu_pdf_Interpolation != data['interpolationInfo'][0][1] ) :
                    warningList.append( '       WARNING: Some interpolation changed for angular data from %s to %s.' % ( mu_pdf_Interpolation, 
                        data['interpolationInfo'][0][1] ) )
                    break
        mu_pdf_InterpolationString = endfFileToGNDMisc.interpolationString( mu_pdf_Interpolation )
        region = distributions.angular.piecewiseRegion( index, endfFileToGNDMisc.interpolationString( interpolationFlag ), mu_pdf_InterpolationString )
        if( not( lastRegion is None ) ) :
            if( lOrTs[start]['C2'] != lastRegion[0] ) : region.append( lastRegion[0], lastRegion[1] )
        priorEnergy = -1.
        for i in xrange( start, end ) :
            lOrT = lOrTs[i]
            if( LANG > 0 ) :
                if( lOrT['L1'] != LANG ) : raise Exception( "LANG = %d changed from %d: LAW = 2, MT = %d" % ( LANG, MT ) )
            if( LANG > 0 ) :
                list = lOrT['data']
                data = []
                for i in xrange( 0, len( list ), 2 ) : data.append( [ list[i], list[i+1] ] )
            else :
                if( ( mu_pdf_Interpolation != lOrT['interpolationInfo'][0][1] ) and doInterpolationChangeWarning ) : 
                    warningList.append( '       WARNING: Some interpolation flags changing to %s for MF = 4, LTT = 2 or 3' % mu_pdf_Interpolation )
                    doInterpolationChangeWarning = False
                data = lOrT['data']
            energy = lOrT['C2']
            if( energy == priorEnergy ) :
                if not allowDuplicateEnergiesInMF4Table:
                    warningList.append( '       WARNING: At Ein = %s in MF = 4, LTT = 2, MT = %s, found duplicate energy' % ( energy, MT ) )
                    continue
                piecewiseForm.append( region )
                index += 1
                region = distributions.angular.piecewiseRegion( index, endfFileToGNDMisc.interpolationString( interpolationFlag ), mu_pdf_InterpolationString )
            region.append( lOrT['C2'], data )       # Currently, only supports one interpolation mu, P region.
            priorEnergy = energy
        lastRegion = [ lOrT['C2'], data ]
        piecewiseForm.append( region )
        index += 1
        start = end
    if( len( piecewiseForm ) == 1 ) :
        energyInterpolation, energyFunctionInterpolation, interpolationQualifier = \
            endfFileToGNDMisc.ENDFInterpolationToGNDAxes3plusd( angularLOrT['interpolationInfo'][0][1] )
        piecewiseRegion = piecewiseForm[0]
        muInterpolation, probabilityInterpolation = piecewiseRegion.yzInterpolationString.split( ',' )
        axes_ = distributions.angular.pointwise.defaultAxes( energyInterpolation = energyInterpolation, 
            energyFunctionInterpolation = energyFunctionInterpolation, muProbabilityFrame = frame, energyInterpolationQualifier = interpolationQualifier,
            muInterpolation = muInterpolation, probabilityInterpolation = probabilityInterpolation )
        pointwiseForm = distributions.angular.pointwise( axes_ )
        axesMuP = axes.referenceAxes( pointwiseForm )
        xi, yi, = piecewiseRegion.yzInterpolationString.split( ',' )
        axesMuP[0].setInterpolation( axes.interpolationXY( xi, yi ) ) # This should be axes logic?????
        for i, energy in enumerate( piecewiseRegion.energies ) :
            pointwiseForm[i] = XYs.XYs( axesMuP, energy.xys, accuracy = ENDF_Accuracy, value = energy.value, index = i, parent = pointwiseForm )
        piecewiseForm = pointwiseForm
    else :
        print LANG, LTT, piecewiseForm.data
        raise Exception( 'piecewise len( %s ) > 1 not supported: MT=%i' % ( len( piecewiseForm ), MT ) )
    return( piecewiseForm )

def toPointwiseEnergy( data, wInterpolation, dependentInterpolation, interpolationQualifier, xInterpolation, yInterpolation, frame ) :

    axes_ = axes.axes( dimension = 3 )
    axes_[0] = axes.axis( 'energy_in', 0, 'eV', frame = axes.labToken, interpolation = axes.interpolationXY( wInterpolation, dependentInterpolation,
        interpolationQualifier ) )
    axes_[1] = axes.axis( 'energy_out', 1, 'eV', frame = frame, interpolation = axes.interpolationXY( xInterpolation, yInterpolation ) )
    axes_[2] = axes.axis( 'P(energy_in|energy_out)', 2, '1/eV', frame = frame )
    form = distributions.energy.pointwise( axes_ )
    axes_xy = axes.axes( )
    axes_xy[0] = axes_[1]
    axes_xy[1] = axes_[2]
    energyPrior = -1
    for energy, EoutP in data :
        if( EoutP[-1][0] == EoutP[-2][0] ) : EoutP[-2][0] *= ( 1. - FUDGE_EPS )
        if( energy == energyPrior ) : energy *= ( 1. + FUDGE_EPS )
        form.append( XYs.XYs( axes_xy, EoutP, accuracy = ENDF_Accuracy, value = energy ) )
        energyPrior = energy
    return( form )

def convertToPointwiseOrPiecewiseEnergy( MT, frame, data ) :

    if( data['NR'] > 1 ) : raise Exception( 'Currently only one interpolation flag is supported for MT=%d, MF=5, LF=1 data' % ( MT ) )
    interpolationE_in, interpolationF, interpolationQualifier = endfFileToGNDMisc.ENDFInterpolationToGNDAxes3plusd( data['interpolationInfo'][0][1] )
    doPointwise, interpolation = True, None
    for energy, region in data['TAB1s'] :
        doPointwise &= ( len( region ) == 1 )
        interpolation2 = region[0].axes[0].interpolation
        if( interpolation is None ) : interpolation = interpolation2
        doPointwise &= ( interpolation == interpolation2 )
    if( doPointwise ) :
        data_ = [ [ energy, xys[0] ] for energy, xys in data['TAB1s'] ]
        independent, dependent, dummy = interpolation.getInterpolationTokens( )
        form = toPointwiseEnergy( data_, interpolationE_in, interpolationF, interpolationQualifier, independent, dependent, frame )
    else :
        for NBT, interpolation in data['interpolationInfo'] :
            energyInterpolation, energyFunctionInterpolation, interpolationQualifier = endfFileToGNDMisc.ENDFInterpolationToGNDAxes3plusd( interpolation )
        axes_ = distributions.energy.semiPiecewise.defaultAxes( frame )
        form = distributions.energy.semiPiecewise( axes_ )
        regionsAxes = axes.referenceAxes( form )
        for energy, regions_ in data['TAB1s'] :
            regions__ = regions.regionsXYs( regionsAxes, parent = form, isPrimaryXData = False )
            for i, region in enumerate( regions_ ) : regions__[i] = region
            form.append( energy, regions__ )
    return( form )

def readMF2( info, MF2, warningList ) :
    """
    parse MF2 into resonances class (and sub-classes)
    cmattoon, 11/09/2010
    """
    PhysQuant = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty
    from fudge.core.math.table import table as gndTable, columnHeader as gndColumnId

    # store MT #s for all reactions that need to include resonance data:
    resonanceMTs = set()
    
    def readResonanceSection( LRU, LRF, NRO, NAPS ):
        """ helper function, read in resonance info for one energy range """
        
        if NRO!=0:  # energy-dependent scattering radius
            if NAPS==2: raise BadResonances("NAPS=2 option not yet supported!")
            line1 = mf2.next()
            dum, dum, dum, dum, NR, NP = funkyFI( line1, logFile = info.logs )
            nLines = NR//3 + bool(NR%3)  +  NP//3 + bool(NP%3)
            data = [line1] + [mf2.next() for i in range(nLines)]
            axes_ = axes.defaultAxes()
            dataLine, TAB1, regions = endfFileToGNDMisc.getTAB1Regions( 0, data, axes_, logFile = info.logs )
            if TAB1['NR']!=1 or TAB1['interpolationInfo'][0][1]!=2:
                raise BadResonances("scattering radius must be pointwise linear")
            data = regions[0]
            data.axes[0].label="energy_in"; data.axes[0].unit="eV"; data.axes[0].frame="lab"
            data.axes[1].label="radius"; data.axes[1].unit="10*fm"; data.axes[1].frame="lab"
            data = data.convertAxisToUnit(1,'fm')
            data.name = "scatteringRadius"
            data.isPrimaryXData = True; data.index = None
            scatRadius = gnd.resonances.scatteringRadius( data )
        
        if LRU==0:  # scattering radius only. Note AP given in 10*fm
            SPI, AP, dum, dum, NLS, dum = funkyFI( mf2.next(), logFile = info.logs )
            info.particleSpins['target'] = (gnd.xParticle.spin(SPI), 0) # no parity information
            scatRadius = gnd.resonances.scatteringRadius( PhysQuant(AP*10,"fm"), EL, EH)
            return scatRadius

        elif LRU==1 and LRF in (1,2):   #SLBW or MLBW form
            SPI, AP, dum, dum, NLS, dum = funkyFI( mf2.next(), logFile = info.logs )
            info.particleSpins['target'] = (gnd.xParticle.spin(SPI), 0)
            if NRO==0:
                scatRadius = gnd.resonances.scatteringRadius( PhysQuant(AP*10,"fm") )
            resList = []
            negativeJs = False
            for lidx in range(NLS):
                AWRI, QX, L, LRX, tmp, NRS = funkyFI( mf2.next(), logFile = info.logs )
                if tmp!=6*NRS:
                    raise BadResonances( "incorrect values in resonance section line %i" % mf2.index )
                for line in range(NRS):
                    e,j,gtot,gn,gg,gf = funkyF( mf2.next(), logFile = info.logs )
                    if j<0: negativeJs = True
                    resList.append( [e,L,j,gtot,gn,gg,gf] )
            if negativeJs: raise BadResonances("Encountered negative J-values for SLBW/MLBW")

            table = gndTable( columns = [
                gndColumnId( 0, name="energy", units="eV" ),
                gndColumnId( 1, name="L", units="" ),
                gndColumnId( 2, name="J", units="" ),
                gndColumnId( 3, name="totalWidth", units="eV" ),
                gndColumnId( 4, name="neutronWidth", units="eV" ),
                gndColumnId( 5, name="captureWidth", units="eV" ),
                gndColumnId( 6, name="fissionWidthA", units="eV" ), ],
                data = sorted(resList, key=lambda(res): res[0]) )   # sort by energy only
            for column in ("totalWidth","fissionWidthA"):
                if not any( table.getColumn(column) ):
                    table.removeColumn(column)

            if LRF==1:
                return gnd.resonances.SLBW( table, scatteringRadius=scatRadius,
                        calculateChannelRadius=not(NAPS) )
            else:
                return gnd.resonances.MLBW( table, scatteringRadius=scatRadius,
                        calculateChannelRadius=not(NAPS) )
        
        elif LRU==1 and LRF==3:     # Reich-Moore form
            SPI, AP, LAD, dum, NLS, NLSC = funkyFI( mf2.next(), logFile = info.logs )
            info.particleSpins['target'] = (gnd.xParticle.spin(SPI), 0) # store spin in GND particle list
            if NRO==0:
                scatRadius = gnd.resonances.scatteringRadius( PhysQuant(AP*10,"fm") )
            resList = []
            LdependentAP = {}
            negativeJs = False  # in Reich-Moore, j<0 indicates channel spin
            for lidx in range(NLS):
                AWRI, APL, L, dum, tmp, NRS = funkyFI( mf2.next(), logFile = info.logs )
                if tmp!=6*NRS:
                    raise BadResonances("incorrect values in resonance section line %i" % mf2.index)
                if APL: # and APL != AP:
                    LdependentAP[L] = PhysQuant(APL*10,"fm")
                for line in range(NRS):
                    e,j,gn,gg,gf1,gf2 = funkyF( mf2.next(), logFile = info.logs )
                    resList.append( [e,L,j,gn,gg,gf1,gf2] )
                    if j<0: negativeJs = True

            columns = [
                gndColumnId( 0, name="energy", units="eV" ),
                gndColumnId( 1, name="L", units="" ),
                gndColumnId( 2, name="J", units="" ),
                gndColumnId( 3, name="neutronWidth", units="eV" ),
                gndColumnId( 4, name="captureWidth", units="eV" ),
                gndColumnId( 5, name="fissionWidthA", units="eV" ),
                gndColumnId( 6, name="fissionWidthB", units="eV" ), ]

            if negativeJs:
                # for incident neutrons, channel spin s = vector sum(target spin, 1/2)
                # in Reich-Moore, the channel spin is encoded in the sign of J
                e,L,j,gn,gg,gf1,gf2 = zip(*resList)
                channelSpins = [SPI + math.copysign(0.5,jval) for jval in j]
                jlist = [abs(jval) for jval in j]
                resList = [list(row) for row in zip( e,L,jlist,channelSpins,gn,gg,gf1,gf2 )]
                columns.insert(3, gndColumnId( 3, name="channelSpin", units="" ) )
                for i in range(len(columns)): columns[i].index = i

            table = gndTable( columns = columns,
                    data = sorted(resList, key=lambda(res): res[0]) ) # sort by resonance energy
            for column in ('fissionWidthA','fissionWidthB'):
                if not any( table.getColumn(column) ):
                    table.removeColumn(column)

            return gnd.resonances.RM( table, scatteringRadius=scatRadius,
                    calculateChannelRadius=not(NAPS), computeAngularDistribution=LAD,
                    LvaluesNeededForConvergence=NLSC, LdependentScatteringRadii=LdependentAP )

        elif LRU==1 and LRF==4:     # Adler-Adler, not currently supported
            raise BadResonances( "Adler-Adler resonance formalism not yet supported!" )

        elif LRU==1 and LRF==7:     # R-Matrix Limited
            dum,dum,IFG,KRM,NJS,KRL = funkyFI( mf2.next(), logFile = info.logs )
            if KRM==3:
                approximation = 'Reich_Moore'
            else:
                raise BadResonances( "R-Matrix with KRM=%i not yet implemented!\n" % KRM )
            
            dum,dum,NPP,dum,tmp1,tmp2 = funkyFI( mf2.next(), logFile = info.logs )
            if tmp1!=12*NPP or tmp2!=2*NPP:
                raise BadResonances( "incorrect LRF7 header!" )
            
            # some helper functions:
            def getOutgoingParticles( MT, targZA, projZA ):
                reacStr = endf_endl.ENDF_MTZAEquation(projZA,targZA, MT)[1]
                outgoing = reacStr.split('->')[1].strip()
                pA, pB = outgoing.split()[::2]
                return pA, pB
            
            def translateENDFJpi(I,P):
                # endf uses weird convention for Jpi. Translate to simpler version:
                spin = abs(I)
                if I: parity = abs(I)/I
                else: parity = P or 1   # if (I,P)=(0,0) treat as '0+'
                return gnd.xParticle.spin( spin ), gnd.xParticle.parity( parity )
                
            # ENDF R-Matrix starts by listing outgoing particle pairs
            # these are referred back to later on. Store in 'IPPlist':
            openChannels = []
            for i in range(NPP):
                MA, MB, ZA, ZB, IA, IB = funkyF( mf2.next(), logFile = info.logs )
                Q, PNT, SHF, MT, PA, PB = funkyF( mf2.next(), logFile = info.logs )
                MT = int(MT)
                resonanceMTs.add(MT)
                
                Qvalue = PhysQuant(Q,'eV')
                calculatePenetrability = False
                if PNT==1 or (PNT==0 and MT in (19,102)):
                    calculatePenetrability = True
                
                # identify the channel using ZA and MT:
                projectileZA = info.reactionSuite.projectile.getZ_A_SuffixAndZA()[-1]
                pA,pB = getOutgoingParticles(MT,int(ZAM),projectileZA)
                # get target spin. In future, this should already be present in particle list
                info.particleSpins[pA] = translateENDFJpi(IA,PA)
                info.particleSpins[pB] = translateENDFJpi(IB,PB)
                
                channel = "%s + %s" % (pA,pB)
                openChannels.append( gnd.resonances.openChannel( i,channel=channel,ENDF_MT=MT,Qvalue=Qvalue,
                    calculatePenetrability=calculatePenetrability,calculateShift=bool(SHF)) )
                
            
            # next we have NJS spin groups, each containing channels and resonances
            spinGroups = []
            for spinGroupIndex in range(NJS):
                # read scattering radius, binding, etc:
                AJ, PJ, KBK, KPS, tmp, NCH = funkyFI( mf2.next(), logFile = info.logs )
                if tmp!=6*NCH:
                    raise BadResonances("incorrect LRF7 header, line %i" % mf2.index)
                idx = 0
                channelData = [ gndColumnId(idx, name="energy", units="eV") ]
                channelNames = []
                for i in range(NCH):
                    idx += 1
                    IPP, L, SCH, BND, APE, APT = funkyF( mf2.next(), logFile = info.logs )
                    thisChannel = openChannels[int(IPP)-1]
                    channelName = "%s width" % thisChannel.channel
                    jdx = 2
                    while True:
                        if channelName not in channelNames:
                            channelNames.append( channelName ); break
                        channelName = '%s width_%i' % (thisChannel.channel, jdx)
                        jdx += 1
                    scatteringRadius = gnd.resonances.scatteringRadius( PhysQuant(APT*10,"fm") )
                    effectiveRadius = gnd.resonances.scatteringRadius( PhysQuant(APE*10,"fm") )
                    channelData.append( gndColumnId(idx, name=channelName, units="eV",
                        L=int(L), channelSpin=SCH,scatteringRadius=scatteringRadius,
                                effectiveRadius=(effectiveRadius if not APT==APE else None),
                                boundaryCondition=BND) )
                
                # resonances for this J:
                dum,dum,dum,NRS,tmp,NX = funkyFI( mf2.next(), logFile = info.logs )
                if tmp!=6*NX:
                    raise BadResonances("incorrect LRF7 header, line %i" % mf2.index)
                if NRS==0: mf2.next()   # skip empty line
                resonances = []
                for i in range(NRS):
                    nlines = int(math.ceil( NCH/6.0 ))
                    vals = []
                    for j in range(nlines):
                        vals += funkyF( mf2.next(), logFile = info.logs )
                    resonances.append( vals[:NCH+1] )
                
                table = gndTable( columns=channelData, data=resonances )
                # done with this spin group:
                J, pi = translateENDFJpi(AJ,PJ)
                spinGroups.append( gnd.resonances.spinGroup(spinGroupIndex, J, pi, KBK, KPS, table) )
            
            # fix up the data a bit: in ENDF, scatteringRadius etc are defined multiple times
            # but are usually all identical. In GND, define only once for the whole resonance region
            # only include in specific channels if we need to override the default value
            kwargs = {}
            chanList = [chan for sg in spinGroups for chan in sg.resonanceParameters.columns]
            for option in ('scatteringRadius','boundaryCondition'):
                vals = [chan.attributes.get(option, None) for chan in chanList]
                valSet = [v for v in set(vals) if v is not None]
                if not valSet: continue
                counts = [vals.count(v) for v in valSet]
                mostCommon = valSet[ counts.index( max(counts) ) ]
                kwargs[option] = mostCommon
                for chan in chanList:
                    opt = chan.attributes.get(option,None)
                    if opt == mostCommon or not opt: chan.attributes[ option ] = None
            # also fix-up the penetrability and phase shift: define only once at the top if possible:
            for option in ('calculatePenetrability','calculateShift'):
                vals = [getattr(chan, option, None) for chan in openChannels]
                valSet = [v for v in set(vals)]
                mostCommon = valSet[ counts.index( max(counts) ) ]
                kwargs[option] = mostCommon
                for chan in openChannels:
                    opt = getattr(chan,option,None)
                    if opt == mostCommon: setattr(chan, option, None)
            
            # end of spin groups. write RMatrix class:
            #spinGroups.sort()   # sort by Jpi. Disable for easier comparison with ENDF-6
            return gnd.resonances.RMatrix( openChannels, spinGroups, approximation=approximation,
                    relativisticKinematics=bool(KRL), reducedWidthAmplitudes=bool(IFG),
                    calculateChannelRadius=not(NAPS), **kwargs )
        
        elif LRU==2: # unresolved
            L_list = []
            LRF_LFW = "LRF,LFW=%i,%i" % (LRF, LFW)
            
            if LFW==0 and LRF==1:   # 'Case A', see ENDF 2010 manual page 70
                SPI,AP,LSSF,dum,NLS,dum = funkyFI( mf2.next(), logFile = info.logs )
                info.particleSpins['target'] = (gnd.xParticle.spin(SPI), 0)
                scatRadius = gnd.resonances.scatteringRadius( PhysQuant(AP*10,"fm") )
                for lidx in range(NLS):
                    AWRI, dum, L, dum, tmp, NJS = funkyFI( mf2.next(), logFile = info.logs )
                    if tmp!=6*NJS:
                        raise BadResonances("bad unresolved flag, line %i" % mf2.index)
                    J_list = []
                    for jidx in range(NJS):
                        D,AJ,AMUN,GNO,GG,dum = funkyF( mf2.next(), logFile = info.logs )
                        if AMUN.is_integer(): AMUN = int(AMUN)
                        D,GNO,GG,GF = [PhysQuant(a,'eV') for a in (D,GNO,GG,0)]
                        J_list.append( gnd.resonances.URR_Jsection( gnd.xParticle.spin(AJ), eDepWidths=[],
                            constantWidths={'levelSpacing':D,'neutronWidth':GNO,'captureWidth':GG, 'fissionWidthA':GF},
                            neutronDOF=AMUN, gammaDOF=0, competitiveDOF=0, fissionDOF=0 ) )
                    L_list.append( gnd.resonances.URR_Lsection( L, J_list ) )
                return gnd.resonances.unresolvedTabulatedWidths( L_list,
                        interpolation = '%s,%s' % (axes.linearToken, axes.linearToken),
                        scatteringRadius=scatRadius, forSelfShieldingOnly=bool(LSSF), 
                        ENDFconversionFlag=LRF_LFW)
            
            elif LFW==1 and LRF==1: # 'Case B'
                SPI,AP,LSSF,dum,NE,NLS = funkyFI( mf2.next(), logFile = info.logs )
                info.particleSpins['target'] = (gnd.xParticle.spin(SPI), 0)
                scatRadius = gnd.resonances.scatteringRadius( PhysQuant(AP*10,"fm") )
                nlines = int(math.ceil(NE/6.0))
                energyList = []
                for i in range(nlines):
                    #energyList += [PhysQuant(a,'eV') for a in funkyF(mf2.next(), logFile = info.logs )]
                    energyList += funkyF(mf2.next(), logFile = info.logs)
                for lidx in range(NLS):
                    AWRI,dum,L,dum,NJS,dum = funkyFI( mf2.next(), logFile = info.logs )
                    J_list = []
                    for jidx in range(NJS):
                        eDepWidths = []
                        dum,dum,L,MUF,tmp,dum = funkyFI( mf2.next(), logFile = info.logs )
                        if tmp!=NE+6:
                            raise BadResonances("Bad unresolved flag, line %i" % mf2.index)
                        D,AJ,AMUN,GNO,GG,dum = funkyF( mf2.next(), logFile = info.logs )
                        if AMUN.is_integer(): AMUN = int(AMUN)
                        D,GNO,GG,GX = [PhysQuant(a,'eV') for a in (D,GNO,GG,0)]
                        widthList = []
                        for i in range(nlines):
                            #widthList += [PhysQuant(a,'eV') for a in funkyF( mf2.next(), logFile = info.logs )]
                            widthList += funkyF( mf2.next(), logFile = info.logs )
                        table = gndTable( [
                            gndColumnId( 0, name="energy", units="eV" ),
                            gndColumnId( 1, name="fissionWidthA", units="eV" ) ] )
                        for e in range(NE):
                            table.addRow( [energyList[e], widthList[e]] )
                            #eDepWidths.append({'energy':energyList[e], 'fissionWidthA':widthList[e]})
                        J_list.append( gnd.resonances.URR_Jsection( gnd.xParticle.spin(AJ), table,
                            constantWidths = { 'levelSpacing' : D, 'neutronWidth' : GNO, 'captureWidth' : GG, 'competitiveWidth' : GX },
                            neutronDOF = AMUN, fissionDOF = MUF, competitiveDOF = 0 ) )
                    L_list.append( gnd.resonances.URR_Lsection( L, J_list ) )
                return gnd.resonances.unresolvedTabulatedWidths( L_list,
                        interpolation = '%s,%s' % (axes.linearToken, axes.linearToken),
                        scatteringRadius=scatRadius, forSelfShieldingOnly=bool(LSSF),
                        ENDFconversionFlag=LRF_LFW)
                
            elif LRF==2:            # 'Case C', most common in ENDF-VII.1
                SPI,AP,LSSF,dum,NLS,dum = funkyFI( mf2.next(), logFile = info.logs )
                info.particleSpins['target'] = (gnd.xParticle.spin(SPI), 0)
                scatRadius = gnd.resonances.scatteringRadius( PhysQuant(AP*10,"fm") )
                interpolations = []
                for Lidx in range(NLS):
                    J_list = []
                    AWRI,dum,L,dum,NJS,dum = funkyFI( mf2.next(), logFile = info.logs )
                    for jidx in range(NJS):
                        resList = []
                        AJ,dum,INT,dum,tmp,NE = funkyFI( mf2.next(), logFile = info.logs )
                        interpolations.append(INT)
                        if tmp!=6*NE+6:
                            raise BadResonances("bad unresolved flag, line %i" % mf2.index)
                        dof = {}
                        dum,dum,dof['AMUX'],dof['AMUN'],dof['AMUG'],dof['AMUF'] = funkyF( mf2.next(), logFile = info.logs )
                        for key in ('AMUX','AMUN','AMUG','AMUF'):
                            if type(dof[key]) is float and dof[key].is_integer(): dof[key] = int(dof[key])
                        for i in range(NE):
                            resList.append( funkyF( mf2.next(), logFile = info.logs ) )
                        table = gndTable( columns= [
                            gndColumnId( 0, name="energy", units="eV" ),
                            gndColumnId( 1, name="levelSpacing", units="eV" ),
                            gndColumnId( 2, name="competitiveWidth", units="eV" ),
                            gndColumnId( 3, name="neutronWidth", units="eV" ),
                            gndColumnId( 4, name="captureWidth", units="eV" ),
                            gndColumnId( 5, name="fissionWidthA", units="eV" ), ],
                            data = resList )
                        J_list.append( gnd.resonances.URR_Jsection( gnd.xParticle.spin(AJ), table,
                                constantWidths={}, competitiveDOF=dof['AMUX'], neutronDOF=dof['AMUN'],
                                gammaDOF=dof['AMUG'], fissionDOF=dof['AMUF'] ) )
                        # sometimes energy-dependent flag is used, but widths are constant in energy:
                        J_list[-1].eliminateRedundantInfo()
                    L_list.append( gnd.resonances.URR_Lsection( L, J_list ) )
                if len(set(interpolations))>1:
                    warningList.append( '       WARNING: inconsistent interpolations in unresolved region will be ignored!' )
                return gnd.resonances.unresolvedTabulatedWidths( L_list,
                        interpolation=endfFileToGNDMisc.interpolationString(interpolations[0]),
                        scatteringRadius=scatRadius, forSelfShieldingOnly=bool(LSSF),
                        ENDFconversionFlag=LRF_LFW)
        
        else:
            info.logs.write( "Unexpected LRU=%i, LRF=%i encountered\n" % ( LRU, LRF ) )
    
    # end of helper functions.
    # now read MF2 data:
    mf2 = myIter(MF2) # mf2.next() to get each line
    
    scatteringRadius = None
    resolvedList = []
    unresolvedList = []
    
    # store particle spins, will be stored in the particle list:
    info.particleSpins = {}
    
    # read MF2 header:
    ZAM, AWR, dum, dum, NIS, dum = funkyFI( mf2.next(), logFile = info.logs )
    countZAMasses( info, ZAM, AWR )
    if NIS!=1: info.logs.write( "careful, more than one isotope in MF2!" )
    ZAI, ABN, dum, LFW, NER, dum = funkyFI( mf2.next(), logFile = info.logs )
    
    for erange in range(NER):
        # each energy range
        EL, EH, LRU, LRF, NRO, NAPS = funkyFI( mf2.next(), logFile = info.logs )
        EL, EH = PhysQuant(EL,"eV"), PhysQuant(EH,"eV")
        nativeData = readResonanceSection( LRU, LRF, NRO, NAPS )
        if LRU==0:
            scatteringRadius = nativeData
        else:
            resonanceMTs.update( [1,2,102] )
            if ZAI//1000>=90: resonanceMTs.update([18,19])
        if LRU==1:
            resolvedList.append( (nativeData,EL,EH) )
        elif LRU==2:
            unresolvedList.append( (nativeData,EL,EH) )
    if not resolvedList: resolved = None
    elif len(resolvedList)==1:
        resolved = gnd.resonances.resolved( *resolvedList[0] )
    else:
        warningList.append( "       WARNING: multiple resolved energy intervals are deprecated!" )
        resolved = gnd.resonances.resolved( multipleRegions=True )
        idx = 0
        for nativeData, EL, EH in resolvedList:
            resolved.regions.append( gnd.resonances.energyInterval(idx,nativeData,EL,EH) )
            idx += 1
    if not unresolvedList: unresolved = None
    elif len(unresolvedList)==1:
        unresolved = gnd.resonances.unresolved( *unresolvedList[0] )
    else:
        raise BadResonances, "multiple unresolved regions not supported"
    
#    if mf2.index != mf2.length: raise BadResonances("Not all resonance data converted!")
    if mf2.index != mf2.length: warningList.append("WARNING: Not all resonance data converted!")
    resonances = gnd.resonances.resonances( scatteringRadius, resolved, unresolved )
    return resonances, sorted(resonanceMTs)


def readMF3( info, MT, MF3Data, warningList ) :

    axes_ = gnd.reactionData.crossSection.pointwise.defaultAxes( )
    dataLine, TAB1, crossSectionRegions = endfFileToGNDMisc.getTAB1Regions( 1, MF3Data, axes_, allowInterpolation6 = True, logFile = info.logs )
    QM, QI, LR = TAB1['C1'], TAB1['C2'], int( TAB1['L2'] )
    breakupProducts = None
    if(   LR == 0 ) : 
        pass
    elif( LR in [ 22, 23, 24, 25, 28, 29, 30, 32, 33, 34, 35, 36 ] ) : 
        info.logs.write( ' : MF=3, LR=%s' % LR )
        breakupProducts, productCounts = {}, endf_endl.endfMTtoC_ProductLists[LR].productCounts
        for product in productCounts :
            if( productCounts[product] != 0 ) : breakupProducts[product] = productCounts[product]
        breakupProducts[info.projectile] -= 1
        if( breakupProducts[info.projectile] == 0 ) : del breakupProducts[info.projectile]
    elif( LR == 31 ) :
        warningList.append( '       WARNING: Invalid LR = %s for MT = %s is being ignored' % ( LR, MT ) )
    elif( LR in [ 1, 39, 40 ] ) : 
        if( LR == 40 ) :
            warningList.append( '       WARNING: LR = %s for MT = %s is being ignored' % ( LR, MT ) )
        else :
            raise Exception( "Breakup LR = %s is not supported: MT = %s" % ( LR, MT ) )
    else :
        raise Exception( "Invalide breakup flag LR %s: MT = %d" % ( LR, MT ) )

    crossSection = gnd.reactionData.crossSection.component( )
    if( len( crossSectionRegions ) == 1 ) :         # Store as pointwise:
        ptwsXSec = getCrossSectionLinearOrPointwise( crossSectionRegions[0].axes, crossSectionRegions[0] )
        crossSection.addFormAsNativeData( ptwsXSec )
    else:
        axes_ = gnd.reactionData.crossSection.piecewise.defaultAxes( )
        piecewise = gnd.reactionData.crossSection.piecewise( axes_ )
        for region in crossSectionRegions : 
            if( len( region ) > 1 ) : piecewise.append( region )
        crossSection.addFormAsNativeData( piecewise )

    if( MT in info.sumCrossSections ) :
        info.sumCrossSections[MT]['total'] = crossSectionRegions
    else :
        for totalMT in info.sumCrossSections :
            if( MT in info.sumCrossSections[totalMT]['MTs'] ) :
                info.sumCrossSections[totalMT]['partialPresent'] = True
                break

    return( QM, QI, crossSection, breakupProducts )

def readMF4( info, product, MT, MF4Data, componentClass, warningList, addAngularComponent = True ) :

    ZA, AWR, LVT, LTT, dummy, dummy = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( MF4Data[0], logFile = info.logs )
    countZAMasses( info, ZA, AWR )
    ZA = int( ZA )
    LVT = int( LVT )                # 1: transformation matrix given. Must be 0 for endf/b6 format but not older formats.
    LTT = int( LTT )                # 0: isotropic, 1: Legendre, 2: table, 3: Legendre for low E and table for high E.

    dummy, AWR, LI, LCT, NK, NM = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( MF4Data[1], logFile = info.logs )
    LI = int( LI )                  # if 1, gammas isotropic
    LCT = int( LCT )                # 1 for lab frame, 2 for center of mass
    NK = int( NK )                  # number of entries in transformation matrix
    NM = int( NM )                  # maximum Legendre order
    if( ( LCT != 2 ) and ( componentClass == distributions.angular.twoBodyComponent ) ): 
        raise Exception( "Discrete two-body must be in the center-of-mass frame: LCT = %d MT = %d." % ( LCT, MT ) )

    firstDataLine = 2
    if( LVT != 0 ) : 
        warningList.append( '       WARNING: MF = 4, MT = 2 contains obsolete matrix used to transform Legendre coefficients between frames.' )
        firstDataLine += ( NK + 5 ) / 6

    info.logs.write( ' : MF=4, LTT = %s' % LTT )
    frame = frames[LCT]
    if( LTT == 0 ) :                # Purely isotropic angular distributions
        form = distributions.angular.isotropic( frame )
        component = componentClass( form.moniker )
    elif( LTT == 1 ) :              # Legendre polynomial coefficient
        nextDataLine, angularData = endfFileToGNDMisc.getTAB2_Lists( firstDataLine, MF4Data, logFile = info.logs )
        form = angularLegendreToPointwiseOrPiecewiseLegendre( MT, frame, angularData, warningList, 4, 'LTT = 1' )
        component = componentClass( form.moniker )
    elif( LTT == 2 ) :              # Tabulated probability distribution
        nextDataLine, angularTable = endfFileToGNDMisc.getTAB2_TAB1s( firstDataLine, MF4Data, logFile = info.logs )
        form = convertAngularToPointwiseOrPiecewise( MT, 0, frame, angularTable, warningList )
        component = componentClass( form.moniker )
    elif( LTT == 3 ) :              # Mixed Legendre and Tabulated probability distribution
        nextDataLine, angularData = endfFileToGNDMisc.getTAB2_Lists( firstDataLine, MF4Data, logFile = info.logs )
        formLegendre = angularLegendreToPointwiseOrPiecewiseLegendre( MT, frame, angularData, warningList, 4, 'LTT = 3' )

        nextDataLine, angularTable = endfFileToGNDMisc.getTAB2_TAB1s( nextDataLine, MF4Data, logFile = info.logs )
        formPointwise = convertAngularToPointwiseOrPiecewise( MT, 0, frame, angularTable, warningList )

        form = distributions.angular.mixedRanges( formLegendre, formPointwise )
        component = componentClass( form.moniker )

    component.addForm( form )
    form.setParent( component )
    if( addAngularComponent ) : product.addDistributionComponent( component )
    return( component )

def readMF5( info, MT, MF5Data, warningList, delayNeutrons = False, product = None ) :

    ZA, AWR, dummy, dummy, NK, dummy = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( MF5Data[0], logFile = info.logs )
    countZAMasses( info, ZA, AWR )
    NK = int( NK )                 # number of partial energy distributions
    dataLine = 1
    energyComponents, energyForms, weights = [], [], []
    for k in xrange( NK ) :
        dataLine, productData = endfFileToGNDMisc.getTAB1( dataLine, MF5Data, oneNR = False, logFile = info.logs )
        if( productData['NR'] > 1 ) :
            oldInterpolations = productData['interpolationInfo']
            if( oldInterpolations == [ [ 3, 1 ], [ 4, 2 ], [ 5, 1 ] ] ) :       # This is a kludge for about 5 data sets, but does a good job.
                productData['NR'] = 1
                productData['data'].insert( 1, [ ( 1. - FUDGE_EPS ) * productData['data'][1][0], productData['data'][0][1] ] )
                productData['interpolationInfo'] = [ [ len( productData['data'] ), 2 ] ]
            else :
                raise Exception( "Currently only one interpolation flag is supported" )
        LF = int( productData['L2'] )           # breakup flag
        xPrior, addWarnging = None, True
        for xy in productData['data'] :
            if( xPrior is not None ) :
                if( xy[0] < xPrior ) : raise Exception( 'xy[0] = %s < xPrior = %s for MT=%d, MF=5' % ( xy[0], xPrior, MT ) )
                if( xy[0] == xPrior ) : 
                    xy[0] *= ( 1 + FUDGE_EPS )
                    if( addWarnging ) : warningList.append( '       WARNING: weights have same energies, second one being incremented for MT=%d, MF=5' % MT )
                    addWarnging = False
            xPrior = xy[0]
        weight = None
        if( ( NK > 1 ) and not( delayNeutrons ) ) :
            if( productData['NR'] != 1 ) : raise Exception( "Currently only one interpolation flag is supported for weight for MF=6, MT=%s" % MT )
            interpolationx, interpolationy = endfFileToGNDMisc.ENDFInterpolationToGND2d( productData['interpolationInfo'][0][1] )
            axes_ = axes.axes( )
            axes_[0] = axes.axis( 'energy_in', 0, 'eV', frame = axes.labToken, interpolation = axes.interpolationXY( interpolationx, interpolationy ) )
            axes_[1] = axes.axis( 'weight',    1, '',   frame = axes.labToken )
            weight = XYs.XYs( axes_,  productData['data'], accuracy = ENDF_Accuracy )
        energyInterpolation, multiplicityInterpolation = endfFileToGNDMisc.ENDFInterpolationToGND2d( productData['interpolationInfo'][0][1] )
        axes_ = gnd.productData.multiplicity.pointwise.defaultAxes( energyInterpolation = energyInterpolation,
            multiplicityName = '', multiplicityInterpolation = multiplicityInterpolation )
        weights.append( XYs.XYs( axes_, data = productData['data'], accuracy = ENDF_Accuracy ) )   # weights is only used for delayed nu_bar data

        frame = frames[1]                       # Always lab frame.
        info.logs.write( ' : MF=5, LF=%s' % LF )
        # 'upper energy limit':
        U = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( productData['C1'], 'eV' )
        if( LF == 1 ) :
            axes_ = axes.axes( )
            axes_[0] = axes.axis( 'energy_out',              0, 'eV',   frame = frame, interpolation = axes.interpolationXY( axes.linearToken, axes.linearToken ) )
            axes_[1] = axes.axis( 'P(energy_out|energy_in)', 1, '1/eV', frame = frame )
            dataLine, EEpETable = endfFileToGNDMisc.getTAB2_TAB1s( dataLine, MF5Data, asRegions = True, axes = axes_, logFile = info.logs )
            form = convertToPointwiseOrPiecewiseEnergy( MT, frame, EEpETable )
            component = distributions.energy.component( form.moniker )
        elif( LF == 5 ) :
            dataLine, thetas = endfFileToGNDMisc.getTAB1( dataLine, MF5Data, logFile = info.logs )
            thetas = endfFileToGNDMisc.toEnergyFunctionalData( 5, 'theta', 'eV', thetas )
            dataLine, gs = endfFileToGNDMisc.getTAB1( dataLine, MF5Data, logFile = info.logs )
            gs = endfFileToGNDMisc.toEnergyFunctionalData( 5, 'g', '1/eV', gs )
            form = distributions.energy.generalEvaporationSpectrum( U, thetas, gs )
            component = distributions.energy.component( form.moniker )
        elif( LF == 7 ) :
            dataLine, thetas = endfFileToGNDMisc.getTAB1( dataLine, MF5Data, logFile = info.logs )
            thetas = endfFileToGNDMisc.toEnergyFunctionalData( 7, 'theta', 'eV', thetas )
            form = distributions.energy.simpleMaxwellianFissionSpectrum( U, thetas )
            component = distributions.energy.component( form.moniker )
        elif( LF == 9 ) :
            dataLine, thetas = endfFileToGNDMisc.getTAB1( dataLine, MF5Data, oneNR = False, logFile = info.logs )
            if( thetas['NR'] > 1 ) :
                oldInterpolations = thetas['interpolationInfo']
                iMin = oldInterpolations[0][1]
                for n, i in oldInterpolations :
                    if( i < iMin ) : iMin = i
                thetas['NR'] = 1
                thetas['interpolationInfo'] = [ [ len( thetas['data'] ), iMin ] ]
                warningList.append( '       WARNING: multiple interpolation = %s reduce to %s for particle %d of %d for MF = 5, LF = %d, MT = %d' % \
                    ( oldInterpolations, thetas['interpolationInfo'], k + 1, NK, LF, MT ) )
            thetas = endfFileToGNDMisc.toEnergyFunctionalData( 9, 'theta', 'eV', thetas )
            form = distributions.energy.evaporationSpectrum( U, thetas )
            component = distributions.energy.component( form.moniker )
        elif( LF == 11 ) :
            dataLine, a = endfFileToGNDMisc.getTAB1( dataLine, MF5Data, logFile = info.logs )
            a = endfFileToGNDMisc.toEnergyFunctionalData( 11, 'a', 'eV', a )
            dataLine, b = endfFileToGNDMisc.getTAB1( dataLine, MF5Data, logFile = info.logs )
            b = endfFileToGNDMisc.toEnergyFunctionalData( 11, 'b', '1/eV', b )
            form = distributions.energy.WattSpectrum( U, a, b )
            component = distributions.energy.component( form.moniker )
        elif( LF == 12 ) :
            dataLine, Ts = endfFileToGNDMisc.getTAB1( dataLine, MF5Data, logFile = info.logs )
            EFL, EFH = [physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty(tmp,'eV') for tmp in (Ts['C1'], Ts['C2'])]
            Ts = endfFileToGNDMisc.toEnergyFunctionalData( 11, 'T_M', 'eV', Ts )
            form = distributions.energy.MadlandNix( EFL, EFH, Ts )
            component = distributions.energy.component( form.moniker )
        else :
            raise Exception( "Unsupported LF = %d" % LF )
        if( not( delayNeutrons ) and ( NK > 1 ) ) : form.weight = weight
        component.addForm( form )

        energyForms.append( form )
        energyComponents.append( component )

    if( not( delayNeutrons ) ) :
        if( NK > 1 ) :
            info.logs.write( 'using energy.weightedFunctionals form' )
            form = distributions.energy.weightedFunctionals( )
            for functional in energyForms :
                weight = distributions.energy.weighted( XYs.XYs( functional.weight.axes, functional.weight,
                    ENDF_Accuracy, isPrimaryXData = True ), functional )
                form.addWeight( weight )
                del functional.weight
            energyComponents = distributions.energy.component( form.moniker )
            energyComponents.addForm( form )
        else :
            energyComponents = energyComponents[0]
    return( energyComponents, weights )

def readMF6( MT, info, MF6Data, productList, warningList, undefinedLevelInfo, isTwoBody ) :

    ZA, AWR, dummy, LCT, NK, dummy = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( MF6Data[0], logFile = info.logs )
    countZAMasses( info, ZA, AWR )
    LCT = int( LCT )
    LCTLight, LCTWeight = LCT, LCT
    if( LCT == 3 ) : LCTLight, LCTWeight = 2, 1
    NK = int( NK )                  # number of outgoing particle data sets

    dataLine, discreteGammas, discretePrimaryGammas = 1, {}, []
    info.logs.write( ' : MF=6' )
    for outGoing in xrange( NK ) :
        ifLegendreConvertedToEnergy = False
        ZAP = int( endfFileToGNDMisc.sixFunkyFloatStringsToFloats( MF6Data[ dataLine ], logFile = info.logs )[0] )
        axes_ = gnd.productData.multiplicity.pointwise.defaultAxes( )
        dataLine, productData, multiplicityRegions = endfFileToGNDMisc.getTAB1Regions( dataLine, MF6Data, axes_, oneNR = ( ZAP != 0 ), logFile = info.logs )
        # ZAP is ZA for outgoing particle; AWP is its atomic mass, LIP: 0 for residual in ground state, 1 for first excited state, etc
        ZAP, AWP, LIP, LAW, NP = int( productData['C1'] ), productData['C2'], productData['L1'], productData['L2'], productData['NR']
        countZAMasses( info, ZAP, AWP )
        if( ZAP not in info.ZAMasses ) : 
            info.ZAMasses[ZAP] = AWP * info.ZAMasses[1]
        elif( info.ZAMasses[ZAP] is None ) :
            info.ZAMasses[ZAP] = -AWP * info.ZAMasses[1]
        ZAP = int( ZAP )            # ZA for outgoing particle; AWP is its atomic mass
        LCTp = LCTLight
        if( ZAP % 1000 > 4 ) : LCTp = LCTWeight
        LAW = int( LAW )
        frame = frames[LCTp]

        info.logs.write( ' : ZAP=%s, LAW=%s' % ( ZAP, LAW ) )
        isTwoBodyGamma = False
        genre = distributions.base.NBodyGenre
        if( LAW == 0 ) :
            if( MT == 5 ) :                 # Need to check if this should be done for all reactions???????//
                genre = distributions.base.unknownComponentToken
                component = distributions.base.unknownComponent( )
                form = None
            else :
                genre = distributions.base.noneComponentToken
        elif( LAW == 1 ) :
            dummy, dummy, LANG, LEP, NR, NE = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( MF6Data[ dataLine ], logFile = info.logs )
            LANG = int( LANG )          # identifies the type of data
            info.logs.write( ', LANG=%s' % LANG )
            LEP = int( LEP )            # interpolation type for outgoing energy
            if( LEP not in [ 1, 2 ] ) : 
                raise Exception( "Currently only LEP = 1 or 2 (flat or linear) interpolation is allowed for MF = 6 data; LEP = %d." % LEP )
            interpolationE_out = axes.linearToken                   # Only works because of LEP check a few lines above.
            interpolationC_l = ENDFInterpolationToGND[LEP]  # Only works because of LEP check a few lines above.
            if( LANG == 1 ) :
                NR = int( NR )              # number of interpolation regions for incident energy
                if( NR != 1 ) : raise Exception( "Currently only one interpolation flag is supported for MF = 6, LAW = 1, LANG = 2; MT = %s" % MT )
                NE = int( NE )              # number of incident energies
        
                dataLine += 1
                EinInterpolationTypes = endfFileToGNDMisc.nStringsToInts( NR, dataLine, MF6Data, dimension = 2 )
                interpolationE_in, interpolationF, interpolationQualifier = endfFileToGNDMisc.ENDFInterpolationToGNDAxes3plusd( EinInterpolationTypes[0][1] )
                dataLine += 1 + (NR-1)/3    # the next data is energy-angle distributions
                axes_ = distributions.Legendre.pointwise.defaultAxes( interpolationE_in, interpolationF, energyInterpolationQualifier = interpolationQualifier,
                    frame = frame, energy_outInterpolation = interpolationE_out, C_lInterpolation = interpolationC_l )
                form = distributions.Legendre.pointwise( axes_ )
                axes__ = axes.referenceAxes( form )
                massRatio = AWR / ( 1. + AWR )
                maxLegendre = 0
                for EinCount in xrange( NE ) :
                    dummy, Ein, ND, NA, NW, NEP = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( MF6Data[ dataLine ], logFile = info.logs )
                    ND = int( ND )          # number of discrete gammas (nonzero only for gammas)
                    NA = int( NA )          # number of angular parameters (i.e., lMax).
                    NW = int( NW )          # number of data values for this incident energy
                    NEP = int( NEP )        # number of outgoing energy values
                    maxLegendre = max( maxLegendre, NA )
                    dataLine += 1
                    if( ND != 0 ) :
                        discreteGammasAtE, EoutData = endfFileToGNDMisc.readDiscreteAndLegendre( ND, NEP-ND, dataLine, MF6Data, dimension = 2+NA, logFile = info.logs )
                        discreteEnergies = [discrete[0] for discrete in discreteGammasAtE]
                        if len(discreteEnergies) != len(set(discreteEnergies)):
                            # Have two or more gammas with identical energy (Fm255 MT91 for example).
                            # For now, just add a small epsilon
                            tmp = []
                            for idx in range(len(discreteEnergies)):
                                while discreteGammasAtE[idx][0] in tmp:
                                    discreteGammasAtE[idx][0] -= FUDGE_EPS
                                tmp.append(discreteGammasAtE[idx][0])
                        discretePrimaryGammasAtE = []
                        for Eg, P in discreteGammasAtE :
                            if( Eg < 0 ) :
                                if( Ein < 1e-2 ) :              # Generally, Eg ~ 1e6 and first Ein ~ 1e-5, so ignore Ein correction here.
                                    Epg = -Eg
                                else :
                                    Epg = -Eg - massRatio * Ein
                                discretePrimaryGammasAtE.append( [ Epg, [ Ein, P ] ] )
                            else :
                                if( Eg not in discreteGammas ) : discreteGammas[Eg] = []
                                discreteGammas[Eg].append( [ Ein, P ] )
                        if( len( discretePrimaryGammasAtE ) > 0 ) :
                            if( len( discretePrimaryGammas ) == 0 ) :
                                discretePrimaryGammas = discretePrimaryGammasAtE
                            else :      # Now we need to match level energies (aka binding energies) from prior with Ein's with current Ein.
                                        # This is needed since for different Ein's, the calculation of Epg will vary slightly (hopefully less than 1e-4).
                                if( len( discretePrimaryGammas ) != len( discretePrimaryGammasAtE ) ) :
                                    raise Exception( 'number of primary gammas at Ein = %s is different then first incident energy' % Ein )
                                for index, Eg in enumerate( discretePrimaryGammas ) :
                                    if( abs( Eg[0] - discretePrimaryGammasAtE[index][0] ) > 1e-4 * Eg[0] ) : raise \
                                        Exception( 'primary energy of %s is not close to primary energy %s' % ( Eg[0], discretePrimaryGammasAtE[index][0] ) )
                                    Eg.append( discretePrimaryGammasAtE[index][1] )
                    else :
                        EoutData = endfFileToGNDMisc.nFunkyFloatStringsToFloats( NEP, dataLine, MF6Data, dimension = 2 + NA, logFile = info.logs )
                    if( len( EoutData ) > 0 ) :
                        # test for some common TENDL problems:
                        # 1: all outgoing energies==0 for continuum gammas.
                        # This usually only affects one incident energy, so just replace that energy with empty outgoing distribution:
                        energy_out_list = [e[0] for e in EoutData]
                        if sum( energy_out_list ) == 0:
                            warningList.append("Format error! At Ein=%s eV, continuum gamma outgoing energies are all 0.0 (MT=%i, ZAP=%i)!"
                                    % (Ein,MT,ZAP))
                            info.doRaise.append( warningList[-1] )
                            EoutData = [[0.0,0.0],[1.0,0.0]]
                        # 2: trouble with duplicate outgoing energies:
                        elif (max( [energy_out_list.count(a) for a in energy_out_list] ) > 2
                                or energy_out_list.count(0.0) > 1):
                            warningList.append("Too many duplicate outgoing energies for Ein=%s eV in W_XYs_LegendreSeries (MT=%i, ZAP=%i)!"
                                    % (Ein,MT,ZAP))
                            info.doRaise.append( warningList[-1] )
                            tmp = []
                            i = 0
                            if energy_out_list.count(0.0)>1:
                                finalZeroIndex = energy_out_list.count(0.0)-1
                                energy_out_list = energy_out_list[ finalZeroIndex: ]
                                EoutData = EoutData[ finalZeroIndex: ]
                            while i<len(energy_out_list):
                                eout = energy_out_list[i]
                                eoutCount = energy_out_list.count(eout)
                                if eoutCount > 2:
                                    tmp.extend( [EoutData[i], EoutData[i+eoutCount-1]] )
                                    i += eoutCount
                                else:
                                    tmp.append( EoutData[i] )
                                    i += 1
                            EoutData = tmp
                        # end of TENDL-specific tests
                        w_xys_LegendreSeries = LegendreSeries.W_XYs_LegendreSeries( axes__, index = EinCount, value = Ein )
                        nm1 = len( EoutData ) - 1
                        for i, EpCs in enumerate( EoutData ) :
                            e_out, data = EpCs[0], EpCs[1:]
                            if( i != nm1 ) :
                                if( e_out == EoutData[i+1][0] ) : e_out *= ( 1 - FUDGE_EPS )
                            w_xys_LegendreSeries[i] = LegendreSeries.XYs_LegendreSeries( axes_[2].getUnit( ), data, index = i, value = e_out )
                        form.append( w_xys_LegendreSeries )
                    dataLine += 1 + ( NW -  1 ) / 6    # Offset for the next incident energy
                component = distributions.Legendre.component( form.moniker )

                if( ( maxLegendre == 0 ) and len( form ) ) :        # Only have l=0 for each outgoing energy. Convert this to
                                                                    # uncorrelated with P(E_out|E_in) and isotropic angular distribution.
                    ifLegendreConvertedToEnergy = True
                    componentAngular = distributions.angular.component( distributions.base.isotropicFormToken )
                    componentAngular.addForm( distributions.angular.isotropic( frame ) )
                    componentEnergy = distributions.energy.component( distributions.base.pointwiseFormToken )
                    pointwiseData = []
                    for w_xys_LegendreSeries in form :
                        energy_in = w_xys_LegendreSeries.value
                        pointwiseData.append( [ energy_in, [ [ lco.value, lco.coefficients[0] ] for lco in w_xys_LegendreSeries ] ] )
                        smearXYs( pointwiseData[-1][1] )
                    form2 = toPointwiseEnergy( pointwiseData, interpolationE_in, interpolationF, interpolationQualifier, interpolationE_out, interpolationC_l, frame )
                    componentEnergy.addForm( form2 )
                    component = distributions.uncorrelated.component( componentAngular, componentEnergy )

            elif( LANG == 2 ) :             # Kalbach-Mann data
                if( LCTp != 2 ) : raise Exception( 'LCT = %s != 2 as required for Kalbach-Mann data for MF=6, MT=%s' % ( LCTp, MT ) )
                dataLine, KalbachMannData = endfFileToGNDMisc.getTAB2_Lists( dataLine, MF6Data, logFile = info.logs )
                if( KalbachMannData['NR'] != 1 ) :
                    raise Exception( "Currently only one interpolation flag is supported for MF = 6, LAW = 1, LANG = 2; MT = %s" % MT )
                interpolationE_in, interpolationF, interpolationQualifier = \
                    endfFileToGNDMisc.ENDFInterpolationToGNDAxes3plusd( KalbachMannData['interpolationInfo'][0][1] )
                axes_ = axes.axes( dimension = 3 )
                axes_[0] = axes.axis( 'energy_in',  0, 'eV', frame = axes.labToken, interpolation = axes.interpolationXY( interpolationE_in, interpolationF, interpolationQualifier ) )
                axes_[1] = axes.axis( 'energy_out', 1, 'eV', frame = frames[2], interpolation = axes.interpolationXY( interpolationE_out, interpolationC_l ) )
                axes_[2] = axes.axis( 'f',          2, '1/eV', frame = frames[2] )
                KMForm = distributions.energyAngular.KalbachMann_Form_fr_Token
                NA = int( KalbachMannData['Lists'][0]['L2'] )
                if( NA == 2 ) : KMForm = distributions.energyAngular.KalbachMann_Form_fra_Token
                form = distributions.energyAngular.KalbachMann( KMForm, axes_ )
                for data in KalbachMannData['Lists'] :
                    # test for another common TENDL problem: repeated incident energy in Kalbach-Mann:
                    try:
                        form.append( distributions.energyAngular.KalbachMannCoefficients( 0, data['C2'], data['data'] ) )
                    except Exception:
                        warningList.append( "Format error! Duplicate KalbachMann incident energy for MT%i, ZAP%i" % (MT,ZAP) )
                        info.doRaise.append( warningList[-1] )
                component = distributions.energyAngular.component( form.moniker )
            else :
                raise Exception( "Unsupported LANG = %d for continuum energy-angle distribution MF = 6: ZAP = %d, LAW = %d: MT = %d" % ( LANG, ZAP, LAW, MT ) )
        elif( LAW == 2 ) :
            if( LCT != 2 ): raise Exception( "Discrete two-body must be in the center-of-mass frame: LCT = %d MT = %d." ( LCT, MT ) )
            isTwoBody = True
            if( ZAP == 0 ) : isTwoBodyGamma = True
            dataLine, angularData = endfFileToGNDMisc.getTAB2_Lists( dataLine, MF6Data, logFile = info.logs )
            LANG = int( angularData['Lists'][0]['L1'] )
            info.logs.write( ', LANG=%s' % LANG )
            if( angularData['NR'] != 1 ) :
                raise Exception( "Currently only one interpolation flag is supported for MF = 6, LAW = 2, LANG = %s; MT = %s" % ( LANG, MT ) )
            interpolationE_in = endfFileToGNDMisc.ENDFInterpolationToGND3plusd( angularData['interpolationInfo'][0][1] )
            if( ( ZAP == 0 ) and ( AWP != 0 ) ) :
                if( LANG != 0 ) : raise Exception( "primary gamma data for LAW=2 with LANG=%s != 0" % ( LANG ) )
                form = []           # This only works as later its length is check and if 0, set to none.
                discretePrimaryGammas = [ AWP ]
                for angularDatum in angularData['Lists'] : discretePrimaryGammas.append( [ angularDatum['C2'], [ 1.0 ] + angularDatum['data'] ] )
                discretePrimaryGammas = [ discretePrimaryGammas ]
            elif( LANG == 0 ) :
                form = angularLegendreToPointwiseOrPiecewiseLegendre( MT, frame, angularData, warningList, 6, 'LAW = 2, LANG = 0' )
                component = distributions.angular.twoBodyComponent( form.moniker )
            elif( LANG in [ 12, 14 ] ) :
                form = convertAngularToPointwiseOrPiecewise( MT, LANG, frame, angularData, warningList )
                component = distributions.angular.twoBodyComponent( form.moniker )
            else :
                raise Exception( "LANG = %d for LAW = %d not supported: MT = %d" % ( LANG, LAW, MT ) )
            genre = distributions.base.angularTwoBodyGenre
        elif( LAW == 3 ) :
            form = distributions.angular.isotropic( frame )
            component = distributions.angular.twoBodyComponent( form.moniker )
            genre = distributions.base.angularTwoBodyGenre
        elif( LAW == 4 ) :
            form = distributions.angular.recoil( product )
            component = distributions.angular.twoBodyComponent( form.moniker )
            genre = distributions.base.angularTwoBodyGenre
        elif( LAW == 5 ) :  # charged-particle elastic scattering
            if( LCT != 2 ): raise Exception( "Charged-particle elastic must be in the center-of-mass frame: LCT = %d MT = %d." ( LCT, MT ) )
            dataLine, angularData = endfFileToGNDMisc.getTAB2_Lists( dataLine, MF6Data, logFile = info.logs )
            SPI = angularData['C1']
            LIDP = angularData['L1']
            LTP = int( angularData['Lists'][0]['L1'] )
            info.logs.write( ', LTP=%s' % LTP )
            interpolationE_in = endfFileToGNDMisc.ENDFInterpolationToGND3plusd( angularData['interpolationInfo'][0][1] )
            # LTP flag changes interpretation of the data:
            if( LTP == 1 ) :
                form = convertNuclearPlusInterferenceDataToPiecewise( MT, frame, angularData, warningList, 6, 'LAW = 5, LTP = %i'%LTP, LIDP )
            elif( LTP == 2 ) :
                raise Exception( "MF=6 LAW=5 LTP=2 not yet implemented (MT%i)!" % MT )
            elif( LTP in ( 12, 14, 15 ) ) :
                form = convertAngularToPointwiseOrPiecewise( MT, LTP, frame, angularData, warningList )
                if( not( isinstance( form, distributions.angular.pointwise ) ) ) : raise Exception( 'Have we hit this yet. If not, we have not tested it: LTP = %d' % LTP )
            else:
                raise Exception( "unknown LTP encountered for MF=6, LAW=5, MT=%s" % MT )
            component = distributions.angular.CoulombElasticComponent( form.moniker, identicalParticles=bool(LIDP), spin=gnd.xParticle.spin(SPI) )
            genre = distributions.base.angularTwoBodyGenre
        elif( LAW == 6 ) :
            APSX, dummy, dummy, dummy, dummy, NPSX = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( MF6Data[ dataLine ], logFile = info.logs )
            dataLine += 1
            APSX *= info.masses.getMassFromZA( 1 )
            componentAngular = distributions.angular.component( distributions.base.isotropicFormToken )
            componentAngular.addForm( distributions.angular.isotropic( frames[2] ) )        # Some data has the wrong frame, should always be com.
            componentEnergy = distributions.energy.component( distributions.base.NBodyPhaseSpaceFormToken )
            form = distributions.energy.NBodyPhaseSpace( int( NPSX ), physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( APSX, 'amu' ) )
            componentEnergy.addForm( form )
            component = distributions.uncorrelated.component( componentAngular, componentEnergy )

        elif( LAW == 7 ) :
            dataLine, NEData = endfFileToGNDMisc.getTAB2Header( dataLine, MF6Data, logFile = info.logs )
            NR = int( NEData['NR'] )                    # number of interpolation regions for this incident energy
            if( NR != 1 ) : raise Exception( "Currently only one interpolation flag is supported for MF = 6, LAW = 7; MT = %s" % MT )
            for iE in xrange( int( NEData['NZ'] ) ) :   # Loop over incident energies
                dataLine, muEpETable = endfFileToGNDMisc.getTAB2_TAB1s( dataLine, MF6Data, logFile = info.logs )
                if( iE == 0 ) :
                    energy_inInterpolation, energy_inFunctionInterpolation, energy_inInterpolationQualifier = \
                        endfFileToGNDMisc.ENDFInterpolationToGNDAxes3plusd( NEData['interpolationInfo'][0][1] )
                    interpolationXYFlag = muEpETable['interpolationInfo'][0][1]
                    muInterpolation, dummy1, dummy2 = endfFileToGNDMisc.ENDFInterpolationToGNDAxes3plusd( interpolationXYFlag )
                    if( dummy1 != axes.linearToken ) : raise Exception( 'mu functional interpolation not linear but %s' % dummy1 )
                    if( dummy2 is not None ) : raise Exception( 'mu interpolation qualifier not None but %s' % dummy2 )
                    energy_outInterpolation, probabilityInterpolation, dummy2 = \
                        endfFileToGNDMisc.ENDFInterpolationToGNDAxes3plusd( muEpETable['TAB1s'][0]['interpolationInfo'][0][1] )
                    if( dummy2 is not None ) : raise Exception( 'energy_out interpolation qualifier not None but %s' % dummy2 )
                    axes_ = distributions.angularEnergy.pointwise.defaultAxes( energyInterpolation = energy_inInterpolation, 
                        energyFunctionInterpolation = energy_inFunctionInterpolation, energyInterpolationQualifier = energy_inInterpolationQualifier,
                        muInterpolation = muInterpolation, energy_outInterpolation = energy_outInterpolation, 
                        probabilityInterpolation = probabilityInterpolation, othersFrame = frame )
                    form = distributions.angularEnergy.pointwise( axes_ )
                    axesW_XY = axes.referenceAxes( form, dimension = 3 )
                    axesXY = axes.referenceAxes( form )
                if( muEpETable['NR'] != 1 ) : raise Exception( "Currently only one interpolation flag is supported for MF = 6, LAW = 7; MT = %s" % MT )
                w_xys = W_XYs.W_XYs( axesW_XY, index = iE, value = muEpETable['C2'] )
                for iMu, TAB1 in enumerate( muEpETable['TAB1s'] ) :
                    if( TAB1['NR'] != 1 ) : raise Exception( "Currently only one interpolation flag is supported for MF = 6, LAW = 7; MT = %s" % MT )
                    if( TAB1['interpolationInfo'][0][1] != interpolationXYFlag ) : 
                        raise Exception( "interpolationXYFlag = %s changing to %s is not allowed" % ( interpolationXYFlag, TAB1['interpolationInfo'][0][1] ) )
                    w_xys.append( XYs.XYs( axesXY, TAB1['data'], accuracy = ENDF_Accuracy, value = TAB1['C2'], index = iMu, parent = w_xys ) )
                form.append( w_xys )
            component = distributions.angularEnergy.component( form.moniker )
        else :
            raise Exception( "LAW = %d not implemented: MT = %d." %  ( LAW, MT ) )

        multiplicityData = productData['data']
        yMin, intergerMultiplicity = multiplicityData[0][1], False
        if( ( ZAP != 0 ) and ( int( yMin ) == yMin ) ) :
            intergerMultiplicity = True
            for x, y in multiplicityData : intergerMultiplicity &= ( yMin == y )
        if( intergerMultiplicity ) :
                multiplicity = int( yMin )
        else :
            multiplicity = getMultiplicityPointwiseOrPieceWise( multiplicityRegions, warningList )
        if( isinstance( multiplicity, gnd.productData.multiplicity.pointwise ) and ( multiplicity.yMin( ) < 0 ) ) :
            warningList.append( "       WARNING: Negative multiplicity encountered for MF6, MT%i %s" % 
                ( MT, getTypeNameENDF( info, ZAP, undefinedLevelInfo ) ) )

        if( ZAP == 0 ) :

            def addGammaProduct( component, form, firstGamma_Multiplicity, attributes = {} ) :

                firstGamma, multiplicity = firstGamma_Multiplicity
                if( firstGamma is not None ) : multiplicity = gnd.productData.multiplicity.reference( firstGamma )
                product = endlToGND.newGNDParticle( info, getTypeNameENDF( info, ZAP, None ), multiplicity = multiplicity, attributes = attributes )
                if( firstGamma is None ) : firstGamma_Multiplicity[0] = product
                component.addForm( form )
                product.distributions.setNativeData( component.moniker )
                product.addDistributionComponent( component )
                productList.append( product )
                return( product )       # May be required for LAW = 4 logic.

            def addGammaWeightedMultiplicity( component, form, firstGamma_Multiplicity, attributes = {} ) :
                """Have Legendre section with only L=0 for discrete gammas. 
                Convert to a weighted multiplicity (with link to continuum gamma mult) along with isotropic angular dist."""

                firstGamma, multiplicity = firstGamma_Multiplicity
                if( firstGamma is not None ) :
                    if( isinstance( form, distributions.angular.LegendrePointwise ) ) :
                        weights = [ [ lco.value, lco.coefficients[0] ] for lco in form ]
                    else :
                        weights = [ [ lco.value, lco.coefficients[0] ] for lco in form.data ]
                    axes_ = gnd.productData.multiplicity.pointwise.defaultAxes( multiplicityName = 'weight' )
                    axes_[1].frame = frame
                    weights = XYs.XYs( axes_, weights, accuracy = ENDF_Accuracy, isPrimaryXData = True )
                    multiplicity = gnd.productData.multiplicity.weightedReference( firstGamma, weights )
                product = endlToGND.newGNDParticle( info, getTypeNameENDF( info, ZAP, None ), multiplicity = multiplicity, attributes = attributes )
                if( firstGamma is None ) : firstGamma_Multiplicity[0] = product
                component.addForm( distributions.angular.isotropic( frame ) )
                component.nativeData = distributions.base.isotropicFormToken
                product.distributions.setNativeData( component.moniker )
                product.addDistributionComponent( component )
                productList.append( product )
                return( product )       # May be required for LAW = 4 logic.

            def addGamma( distinctType, Eg, ELegenres, firstGamma_Multiplicity, axes ) :

                form = distributions.angular.LegendrePointwise( axes )
                for i, ELegenre in enumerate( ELegenres ) :
                    form[i] = LegendreSeries.XYs_LegendreSeries( axes[1].getUnit( ), ELegenre[1], index = i, value = ELegenre[0] )
                attributes = { distinctType : physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( Eg, 'eV' ) }
                if( isTwoBodyGamma ) :
                    component = distributions.angular.twoBodyComponent( form.moniker )
                else :
                    component = distributions.angular.component( form.moniker )

                maxLegendre = max( [ len( e[1] ) for e in ELegenres ] )
                if( maxLegendre == 1 ) and ( firstGamma_Multiplicity[0] is not None  ) : # Only have L_0 for each incident energy, convert to weighted multiplicity
                    return( addGammaWeightedMultiplicity( component, form, firstGamma_Multiplicity, attributes ) )

                return( addGammaProduct( component, form, firstGamma_Multiplicity, attributes ) )

            interpolationQualifier = None
            if( form != [] ) :
                dgInterpolationE_in, dgInterpolationCl, dummy = form.axes[0].interpolation.getInterpolationTokens( )
            else :  
                dgInterpolationE_in, dgInterpolationCl = axes.linearToken, axes.linearToken
            axes_ = distributions.angular.LegendrePointwise.defaultAxes( frame, dgInterpolationE_in, dgInterpolationCl, 
                interpolationQualifier = interpolationQualifier )
            firstGamma_Multiplicity = [ None, multiplicity ]
            if( len( form ) != 0 ) : product = addGammaProduct( component, form, firstGamma_Multiplicity )
            for Eg in sorted( discreteGammas ) : product = addGamma( 'discrete', Eg, discreteGammas[Eg], firstGamma_Multiplicity, axes_ )
            for EgEinPs in discretePrimaryGammas : product = addGamma( 'primary', EgEinPs[0], EgEinPs[1:], firstGamma_Multiplicity, axes_ )
        else :
            thisParticle = getTypeNameENDF( info, ZAP, undefinedLevelInfo )
            if( LIP > 0 ) :
                try:
                    thisParticle = thisParticle.levels[ int(LIP) ]
                except:
                    # this excited level is missing from MF8. Make a level with unknown energy:
                    excitedStateName = thisParticle.name + "_e%i" % int(LIP)
                    unknownLevel = gnd.miscellaneous.undefinedLevel( physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty(0,'eV') )
                    thisParticle.levels[ int(LIP) ] = gnd.xParticle.nuclearLevel( excitedStateName, unknownLevel,
                            int(LIP), groundState=thisParticle )
                    thisParticle = thisParticle.levels[ int(LIP) ]
                    warningList.append( "       WARNING: unknown excited state %s encountered in MF6" % thisParticle )
            product = endlToGND.newGNDParticle( info, thisParticle, multiplicity = multiplicity )
            if( genre != distributions.base.noneComponentToken ) :
                if( genre != distributions.base.unknownComponentToken ) : component.addForm( form )
                product.distributions.setNativeData( component.moniker )
                product.addDistributionComponent( component )
            productList.append( product )
            if( ifLegendreConvertedToEnergy ) : product.addAttribute( 'ENDFconversionFlag', 'MF6' )

        if( ( len( discreteGammas ) > 0 ) and ( productData['interpolationInfo'][0][1] != 2 ) ) :
            info.logs.write( 'interpolation %s is not linear for gamma multiplicity' % productData['interpolationInfo'][0][1] )

    for product in productList :
        if( product.getName( ) in [ 'n', 'gamma' ] ) : product.addAttribute( 'ENDFconversionFlag', 'MF6' )
    return( isTwoBody )

def readMF8( info, MT, MTData, warningList ) :

    if MT in [ 454, 459 ]: # this is fission yield data
        MF8Data = MTData[8]
        ZA, AWR, LE1 = endfFileToGNDMisc.sixFunkyFloatStringsToIntsAndFloats( MF8Data[0], intIndices = [ 0, 2 ], logFile = info.logs )[0:3]
        iLine = 1
        allData = {}
        for iLE in range( LE1 ):
            iLine, thisData = endfFileToGNDMisc.getList( iLine, MF8Data, logFile = info.logs )
            NFP = thisData[ 'N2' ]
            E = thisData[ 'C1' ]
            for iFP in range( NFP ):
                ZAFP = int( thisData[ 'data' ][ iFP*4 ] )
                FPS =  int( thisData[ 'data' ][ iFP*4 + 1 ] )
                Y =    thisData[ 'data' ][ iFP*4 + 2 ]
                DY =   thisData[ 'data' ][ iFP*4 + 3 ]
                key = ( ZAFP, FPS )
                if key not in allData: allData[key] = []
                allData[key].append( ( E, Y, DY ) )
        return allData

    else: # this is regular decay data
        firstLMF, radioactiveDatas = None, []
        if( 8 in MTData ) :
            dataLine, MF8Data = 1, MTData[8]
            ZA, AWR, LIS, LISO, NS, NO = endfFileToGNDMisc.sixFunkyFloatStringsToIntsAndFloats( MF8Data[0], intIndices = [ 0, 2, 3, 4, 5 ], logFile = info.logs )
            countZAMasses( info, ZA, AWR )
            MF9Data = readMF9or10( info, MT, MTData, 9, LIS )
            MF10Data = readMF9or10( info, MT, MTData, 10, LIS )
            residualXSecWeight = []
            for i in xrange( NS ) :
                ZAP, ELFS, LMF, LFS, ND6, dummy = endfFileToGNDMisc.sixFunkyFloatStringsToIntsAndFloats( MF8Data[dataLine], intIndices = [ 0, 2, 3, 4 ], logFile = info.logs )
                if( firstLMF is None ) : firstLMF = LMF
                if( LMF != firstLMF ) : raise Exception( 'LMF changing from %s to %s is not allowed' % ( firstLMF, LMF ) )
                dataLine += 1
                level = ELFS
                try :
                    residual = getTypeNameGamma( info, ZAP, level = level, levelIndex = LFS )
                except :
                    info.logs.write( '\nMT = %s\n' % MT )
                    raise
                crossSection, weight = None, None
                QM, QI = None, None
                if( LMF == 3 ) :
                    pass
                elif( LMF == 6 ) :
                    if ND6 != 0:
                        raise Exception( 'LMF = 6 not supported for MF = 8 data, MT = %s' % MT )
                elif( LMF in [ 9, 10 ] ) :
                    if( LMF == 9 and not MF9Data or LMF == 10 and not MF10Data) :
                        LMF1 = (10 if LMF==9 else 9)
                        warningList.append( 'LMF = %i, but no MF=%i found for MT%i. Trying LMF=%i instead' % ( LMF, LMF, MT, LMF1 ) )
                        info.doRaise.append( warningList[-1] )
                        LMF = LMF1
                    if( LMF == 9 ):
                        TAB1, weight = MF9Data[i]
                    else :
                        TAB1, crossSection = MF10Data[i]
                    QM, QI, IZAP, LFS = TAB1['C1'], TAB1['C2'], int( TAB1['L1'] ), int( TAB1['L2'] )
                    ELFS9or10 = QM - QI
                    if( abs( ELFS - ELFS9or10 ) > 2e-4 * abs( ELFS ) ) : 
                        warningList.append( "MF8 residual level energy = %s for level %s of ZA = %d not close to MF%s's value = %s for MT = %s"
                            % ( ELFS, LIS, ZAP, LMF, ELFS9or10, MT ) )
                        info.doRaise.append( warningList[-1] )
                        ELFS9or10 = ELFS    # assume MF=8 has the correct value
                radioactiveDatas.append( [ ZAP, residual, crossSection, weight, LFS, ELFS, QM, QI ] )

                if( NO == 0 ) : dataLine += ( ND6 + 5 ) / 6
        else :
            radioactiveDatas.append( [ None, None, None, None, 0, 0, 0, 0 ] )
        return( firstLMF, radioactiveDatas )

def readMF9or10( info, MT, MTData, MF, targetLIS ) :

    if( MF not in MTData.keys( ) ) : return( None )
    dataLine, MFData, MF9or10 = 1, MTData[MF], []
    ZA, AWR, LIS, dummy, NS, dummy = endfFileToGNDMisc.sixFunkyFloatStringsToIntsAndFloats( MFData[0], intIndices = [ 0, 2, 4 ], logFile = info.logs )
    countZAMasses( info, ZA, AWR )
    if( LIS != targetLIS ) : raise Exception( "residual's LIS = %s does not match target's LIS = %s" % ( LIS, targetLIS ) )
    if( MF == 10 ) :
        axes_ = gnd.reactionData.crossSection.pointwise.defaultAxes( crossSectionUnit = 'b' )
    else:
        axes_ = axes.defaultAxes(labelsUnits={0:('energy_in','eV'), 1:('weight','')})
    for i in xrange( NS ) :
        QM, QI, IZAP, LFS, NR, NP = endfFileToGNDMisc.sixFunkyFloatStringsToIntsAndFloats( MFData[0], intIndices = [ 0, 2, 4 ], logFile = info.logs )
        dataLine, TAB1, regions = endfFileToGNDMisc.getTAB1Regions( dataLine, MFData, axes_, oneNR = ( MF == 9 ), logFile = info.logs )
        if( MF == 10 ) :
            if( len( regions ) == 1 ) :         # Store as pointwise:
                XSec = getCrossSectionLinearOrPointwise( regions[0].axes, regions[0] )
            else:
                axes_ = gnd.reactionData.crossSection.piecewise.defaultAxes( crossSectionUnit = 'b' )
                XSec = gnd.reactionData.crossSection.piecewise( axes_ )
                for i, region in enumerate( regions ) : XSec[i] = region
        else :
            XSec = regions
        MF9or10.append( [ TAB1, XSec ] )
    return( MF9or10  )

def readMF12_13( info, MT, MTData, productList, warningList, gammaBRTolerance = 1e-6 ) :

    def addMF12_13GammaToList( gList, EGk, ESk, LP, LF, region ) :
        """EGk is the gamma's energy and ESk is the gamma's origination level energy."""

        if( EGk in gList ) : raise Exception( 'Gammas with the same energy (%s) are not supported: MT=%s' % ( EGk, MT ) )
        gList.append( { 'EGk' : EGk, 'ESk' : ESk, 'LP' : LP, 'LF' : LF, 'yield' : region, 'angular' : None, 'energy' : None } )
        #if( ( EGk != 0. ) and ( len( region ) > 1 ) ) :     # There are many continuum gammas with multiple yield regions
        #    raise Exception( 'Gamma yield with multiple regions not implemented: MT=%s EGk=%s' % ( MT, EGk ) )

    def checkAngularData( gamma ) :

        angularForm = gamma['angular']
        if( angularForm is None ) :
            angularForm = distributions.angular.isotropic( frames[1] )             # Why lab frame????????
        angularComponent = distributions.angular.component( angularForm.moniker )
        angularComponent.addForm( angularForm )
        gamma['angular'] = angularComponent

    def addGammaProduct( info, MF, gamma, attributes, component, productList, warningList, ESk ) :

        yields_ = gamma['yield']
        if( MF == 12 ) :
            if( len( yields_ ) == 1 ) :
                multiplicity = getMultiplicityPointwiseOrPieceWise( yields_, warningList, yUnit = yields_[0].axes[1].getLabel( ) )
            else :
                multiplicity = gnd.productData.multiplicity.piecewise( yields_[0].axes )
                for i, yield_ in enumerate( yields_ ) : multiplicity[i] = yield_
        else :
            if( len( yields_ ) == 1 ) :
                mult = yields_[0].copy()
                xsc = info.crossSection.toPointwiseLinear()
                for index, xy in enumerate(mult):
                    xsc_i = xsc.getValue( xy[0] )
                    if xsc_i != 0: mult[ index ] = [ xy[0], xy[1] / xsc_i ]
                multiplicity = getMultiplicityPointwiseOrPieceWise( [mult], warningList,
                        yUnit = yields_[0].axes[1].getLabel( ) )
            else:
                raise Exception( 'Multiple MF=13 regions for MT=%s not yet supported' % MT )
        if( gamma['ESk'] != 0 ) : attributes['originationLevel'] = "%s" % physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( ESk, " eV" )
        product = endlToGND.newGNDParticle( info, getTypeNameENDF( info, 0, None ), multiplicity = multiplicity, attributes = attributes )
        product.distributions.setNativeData( component.moniker )
        product.addDistributionComponent( component )
        if MF==13: product.attributes['ENDFconversionFlag'] = 'MF13'
        productList.append( product )

    if( 12 in MTData ) :
        if( 13 in MTData ) : raise Exception( 'MF = 12 and 13 present for MT=%s, this is not suppored' % MT )
        MF = 12
    elif( 13 in MTData ) :
        MF = 13
    elif( ( 14 in MTData ) or ( 15 in MTData ) ) :
        raise Exception( 'MF 14 and/or 15 data and no MF 12 or 13 data: MT=%s MFs=%s' % ( MT, MTData.keys( ) ) )
    else :
        return

    MF12_13Data = MTData[MF]
    ZA, AWR, LO, LG, NK, dummy = endfFileToGNDMisc.sixFunkyFloatStringsToIntsAndFloats( MF12_13Data[0], intIndices = [ 0, 2, 3, 4 ], logFile = info.logs )
    countZAMasses( info, ZA, AWR )
    info.logs.write( ' : MF=%s LO=%s : ZAP=0 ' % ( MF, LO ) )
    dataLine, continuousGamma, discreteGammas, primaryGammas, branchingGammas = 1, [], [], [], []
    if( ( ( MF == 12 ) and ( LO == 1 ) ) or ( ( MF == 13 ) and ( LO == 0 ) ) ) :
        axes_ = gnd.productData.multiplicity.pointwise.defaultAxes( multiplicityName = 'yield' )
        if( NK > 1 ) : dataLine, TAB1, regions = endfFileToGNDMisc.getTAB1Regions( dataLine, MF12_13Data, axes_, logFile = info.logs ) # Total is not stored.
        keys = []
        for i in xrange( NK ) :
            dataLine, TAB1, regions = endfFileToGNDMisc.getTAB1Regions( dataLine, MF12_13Data, axes_, logFile = info.logs )
            EGk, ESk, LP, LF = TAB1['C1'], TAB1['C2'], int( TAB1['L1'] ), int( TAB1['L2'] )
            if( EGk == 0. ) :
                if( LP not in [ 0, 1 ] ) : raise Exception( 'LP = %s != 0 for continuous gamma for MT = %s' % ( LP, MT ) )
                if( len( continuousGamma ) == 0 ) :
                    addMF12_13GammaToList( continuousGamma, EGk, ESk, LP, LF, regions )
                else :
                    raise Exception( 'continuous gamma information for MF=%s, MT = %s already exist' % ( MF, MT ) )
            elif( LP == 2 ) :
                addMF12_13GammaToList( primaryGammas, EGk, ESk, LP, LF, regions )
            else :
                addMF12_13GammaToList( discreteGammas, EGk, ESk, LP, LF, regions )
    elif( ( MF == 12 ) and ( LO == 2 ) ) :
        dataLine, LO2 = endfFileToGNDMisc.getList( dataLine, MF12_13Data, logFile = info.logs )
        LP = int(LO2['L1'])
        if LP==2 and MT not in (91,649,699,749,799,849):
            warningList.append("       WARNING: Incorrect 'primary gamma' flag for MF12 MT%i" % MT)
        NT, LGp = LO2['N2'], LG + 1
        NK = NT
        for i in xrange( NT ) :
            parentEnergy, finalEnergy = LO2['C1'], LO2['data'][i*LGp]
            if parentEnergy==0:
                raise Exception("Gamma decay from ground state in MF12 MT%i" % MT)
            if abs((parentEnergy-finalEnergy)/parentEnergy)<0.001:
                raise Exception("Zero-energy gamma from %f eV to %f eV in MF12 MT%i" % (parentEnergy,finalEnergy,MT))
            branchingGammas.append( {'ES' : LO2['C1'], 'EGk' : 0, 'ESk' : LO2['data'][i*LGp], 'angular' : None, 'branching' : LO2['data'][i*LGp+1:i*LGp+LGp]} )
        gammaBRList = [ g['branching'][0] for g in branchingGammas ]
        sumGammaBRList = sum( gammaBRList )
        if abs( sumGammaBRList - 1.0 ) > gammaBRTolerance: 
            warningList.append( "       WARNING: sum of gamma BR's for MT="+str(MT)+" MF=12 is " + str(sumGammaBRList)+' != 1.0' )
            # leave re-normalization for checker/fixer codes
            #for i in xrange( len( branchingGammas ) ): branchingGammas[ i ][ 'branching' ][ 0 ] /= sumGammaBRList
        info.MF12_LO2[MT] = branchingGammas
    else :
        raise Exception( 'LO=%s is not valid for MF=%s, MT=%s' % ( LO, MF, MT ) )

    readMF14( info, MT, MTData, MF, NK, warningList, discreteGammas, primaryGammas, continuousGamma, branchingGammas )
    readMF15( info, MT, MTData, continuousGamma, warningList )

    for gamma in branchingGammas :
        if( not( gamma['angular'] is None ) ) :
            info.logs.write( 'NON-isotropic gamma' )
            break
    if( 14 in MTData ) : info.logs.write( ': MF=14 ' )
    if( 15 in MTData ) : info.logs.write( ': MF=15 ' )

    products = []
    if( len( continuousGamma ) ) :
        gamma = continuousGamma[0]
        checkAngularData( gamma )
        component = distributions.uncorrelated.component( gamma['angular'], gamma['energy'] )
        addGammaProduct( info, MF, gamma, {}, component, productList, warningList, gamma['ESk'] )
    for gamma in discreteGammas :
        checkAngularData( gamma )
        attributes = { 'discrete' : physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( gamma['EGk'], " eV" ) }
        addGammaProduct( info, MF, gamma, attributes, gamma['angular'], productList, warningList, gamma['ESk'] )
    for gamma in primaryGammas :
        checkAngularData( gamma )
        attributes = { 'primary' : physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( gamma['EGk'], " eV" ) }
        addGammaProduct( info, MF, gamma, attributes, gamma['angular'], productList, warningList, gamma['ESk'] )
    for gamma in branchingGammas : checkAngularData( gamma )

def readMF14( info, MT, MTData, MF, NK, warningList, discreteGammas, primaryGammas, continuousGamma, branchingGammas ) :

    allGammas = discreteGammas + primaryGammas + continuousGamma + branchingGammas

    def addAngularData( EGk, ESk, form ) :

        possibleGammas = []
        for gamma in allGammas :
            if( ( abs( EGk - gamma['EGk'] ) < 1e-4 ) and ( abs( ESk - gamma['ESk'] ) < 1e-4 ) ) : possibleGammas.append( gamma )
        if( len( possibleGammas ) == 0 ) :
            sys.stderr.write( "MF = 14 gamma data: EGk = %s, ESk = %s\n" % ( EGk, ESk ) )
            sys.stderr.write( 'MF = %s gamma list is,\n' % MF )
            for gamma in allGammas : sys.stderr.write( "    EGk = %s, ESk = %s\n" % ( gamma['EGk'], gamma['ESk'] ) )
            raise Exception( 'No matching gamma from MF = %s found for MF = 14 gamma: MT=%s: see above' % ( MF, MT ) )
        if( len( possibleGammas ) != 1 ) : raise Exception( 'Multiple possible gammas found for MF=14 data: EGk=%s, ESk=%s, MT=%s' % ( EGk, ESk, MT ) )
        gamma = possibleGammas[0]
        if( gamma['angular'] is None ) :
            gamma['angular'] = form
        else :
            raise Exception( 'Gamma already has MF=14 angular data: EGk=%s, ESk=%s, MT=%s' % ( EGk, ESk, MT ) )

    if( not( 14 in MTData ) ) : return
    MF14Data = MTData[14]
    ZA, AWR, LI, LTT, NK14, NI = endfFileToGNDMisc.sixFunkyFloatStringsToIntsAndFloats( MF14Data[0], intIndices = [ 0, 2, 3, 4, 5 ], logFile = info.logs )
    countZAMasses( info, ZA, AWR )
    if( NK14 != NK ) :
        if( not( ( NK14 == 1 ) and ( LTT == 0 ) and info.ignoreBadNK14 ) ) :
            warningList.append( '       WARNING: MF14 NK = %s != MF12/13 NK = %s for MT = %s' % ( NK14, NK, MT ) )
            info.doRaise.append( warningList[-1] )
    dataLine, frame = NI + 1, axes.labToken
    if( LTT == 0 ) :                                     # All distributions are isotropic
        pass
    elif( LTT == 1 ) :
        for i in xrange( NK14 - NI ) :
            dataLine, angularData = endfFileToGNDMisc.getTAB2_Lists( dataLine, MF14Data, logFile = info.logs )
            if( angularData['NR'] != 1 ) : raise Exception( 'Currently only one interpolation flag is supported: NR=%s, MT=%s' % ( angularData['NR'], MT ) )
            EGk, ESk = angularData['C1'], angularData['C2']

            dgInterpolationE_in, dgInterpolationCl, interpolationQualifier = endfFileToGNDMisc.ENDFInterpolationToGNDAxes3plusd( angularData['interpolationInfo'][0][1] )
            axes_ = distributions.angular.LegendrePointwise.defaultAxes( frame, dgInterpolationE_in, dgInterpolationCl )
            lists = angularData['Lists']
            form = distributions.angular.LegendrePointwise( axes_ )
            for i, ELegenre in enumerate( lists ) :
                form[i] = LegendreSeries.XYs_LegendreSeries( axes_[1].getUnit( ), [ 1.0 ] + ELegenre['data'], index = i, value = ELegenre['C2'] )
            addAngularData( EGk, ESk, form )
    else :
        raise Exception( 'MF=14, LI=%s not implemented' % LI )

def readMF15( info, MT, MTData, continuousGamma, warningList ) :

    if( not( 15 in MTData ) ) :
        if( len( continuousGamma ) ) : raise Exception( 'Continous gamma with no MF=15 data: MT=%s' % MT )
        return
    if( len( continuousGamma ) == 0 ) : raise Exception( 'MF=15 data and no continous gamma MF=%s data: MT=%s' % ( MT, 13 ) )
    MF15Data = MTData[15]
    ZA, AWR, dummy, dummy, NC, dummy = endfFileToGNDMisc.sixFunkyFloatStringsToIntsAndFloats( MF15Data[0], intIndices = [ 0, 4 ], logFile = info.logs )
    countZAMasses( info, ZA, AWR )
    if( NC != 1 ) : raise Exception( 'NC = %s > 1 is not supported: MT=%s' % ( NC, MT ) )
    dataLine, weights = endfFileToGNDMisc.getTAB1( 1, MF15Data, logFile = info.logs )    # Why are we not checking weights here???????
    LF = int( weights['L2'] )
    if( LF != 1 ) : raise Exception( 'For MF=15 data, only LF=1 is supported, not LF=%s : MT=%s' % ( LF, MT ) )
    axes_ = axes.axes( )
    axes_[0] = axes.axis( 'energy_out',              0, 'eV', frame = frames[1], interpolation = axes.interpolationXY( axes.linearToken, axes.linearToken ) )
    axes_[1] = axes.axis( 'P(energy_out|energy_in)', 1, '1/eV', frame = frames[1] )
    dataLine, EEpETable = endfFileToGNDMisc.getTAB2_TAB1s( dataLine, MF15Data, axes = axes_, logFile = info.logs )
    pointwiseData, interpolation = [], EEpETable['TAB1s'][0]['interpolationInfo'][0][1]
    interpolationE_out, interpolationP = endfFileToGNDMisc.ENDFInterpolationToGND2d( interpolation )
    interpolationE_in, interpolationF, interpolationQualifier = endfFileToGNDMisc.ENDFInterpolationToGNDAxes3plusd( EEpETable['interpolationInfo'][0][1] )
    if( ( info.projectile == 'n' ) and ( info.target.getName( ) == 'Ni59' ) and ( MT == 102 ) ) :
        specialNi59 = ( EEpETable['TAB1s'][0]['interpolationInfo'][0][1] == 1 )
        for tab1 in EEpETable['TAB1s'] :
            if( len( tab1['interpolationInfo'] ) != 1 ) :
                specialNi59 = False
            elif( tab1['interpolationInfo'][0][1] not in [ 1, 2 ] ) :
                specialNi59 = False
        if( specialNi59 ) :
            warningList.append( '       WARNING: Changing linear interpolation to flat for special Ni59 gamma energy data' )
            for tab1 in EEpETable['TAB1s'] : tab1['interpolationInfo'][0][1] = 1
    for tab1 in EEpETable['TAB1s'] :
        if( len( tab1['interpolationInfo'] ) != 1 ) :
            raise Exception( 'Need to implement piecewise energy here. Interpolations differ, %s vs %s, MF=15, MT=%s' % \
                ( interpolation, tab1['interpolationInfo'], MT ) )
        if( interpolation != tab1['interpolationInfo'][0][1] ) :
            raise Exception( 'Need to implement piecewise energy here. Interpolations differ, %s vs %s, MF=15, MT=%s' % \
                ( interpolation, tab1['interpolationInfo'], MT ) )
        pointwiseData.append( [ tab1['C2'], tab1['data'] ] )
    form = toPointwiseEnergy( pointwiseData, interpolationE_in, interpolationF, interpolationQualifier, interpolationE_out, interpolationP, frames[1] )
    component = distributions.energy.component( form.moniker )
    component.addForm( form )
    continuousGamma[0]['energy'] = component

# for covariances we need a unique id for each section, and also need link info.
# this is messy: lots of special cases
def genID( cov_info, MT, MF, MT2=None, MF2=None, MAT2=None, QI=None, QI2=None ):
    MTdict = cov_info['MTdict']
    def getReaction( MT, MF, QI=None ):
        # find the <reaction> corresponding to this covariance info, if any
        if MT in (452,455,456):
            MTold = MT
            MT = 18 # nubar actually refers to fission
        if MT in MTdict:
            chThisMT = MTdict[MT]
            if len(chThisMT)==1:
                return chThisMT[0]
            elif MF==40:   # MF40 section may contain covariances for multiple excited states of residual
                thisChannel = [ch for ch in chThisMT if ch.getQ('eV')==QI and isinstance(ch, gnd.production.production)]
                if len(thisChannel)==1: return thisChannel[0]
                elif MT in cov_info:
                    # Q-values disagree between MF10/MF40. Assume 1st covariance is for lowest-energy residual, etc
                    index = cov_info[MT]
                    cov_info[MT] += 1
                else:
                    index = 0; cov_info[MT] = 1
                return thisChannel[ index ]
            elif MF==33:
                # residual must be in ground state unless MT in 51-90, etc
                thisChannel = [ch for ch in chThisMT if (isinstance(ch, gnd.reaction.reaction) and
                        sum( [p.getLevelAsFloat('eV') for p in ch.outputChannel] )==0) or
                        isinstance(ch, gnd.summedReaction.summedReaction)]
                if len(thisChannel) != 1: raise BadCovariance("MF33 MT%i covariance doesn't correspond to any channel!"
                        % MT )
                return thisChannel[0]
            else: raise BadCovariance("Can't determine which reaction this covariance (MF%i MT%i) corresponds to!"
                    % (MF,MT))
        return

    rowReaction = getReaction( MT, MF, QI )

    def makeID(MT, reaction):
        if reaction is not None:
            #Id = str( reaction.outputChannel )
            Id = str( reaction )
        elif MT in (452,455,456):
            Id = {452:'total', 455:'delayed', 456:'prompt'}[MT]
        elif MT in (1,4,103,104,105,106,107):
            Id = {1:'total', 4:'sum(z,n)', 103:'sum(z,p)', 104:'sum(z,d)',
                    105:'sum(z,t)', 106:'sum(z,He3)', 107:'sum(z,a)'}[MT]
        elif 850 < MT < 871:
            Id = "lump%i" % (MT-851)
        else: Id = "MF%i_MT%i" % (MF,MT)
        return Id

    rowId = makeID( MT, rowReaction )

    if MAT2:
        # cross-material covariance
        colReaction = None  # must create an 'externalReaction' and link to it below
        versus = ' vs. '
        colId = "MAT%i_MF%i_MT%i" % (MAT2,MF2,MT2)
    elif (MT2 and MT2!=MT):
        MF2 = MF2 or MF
        colReaction = getReaction( MT2, MF2, QI2 )
        versus = ' vs. '
        colId = makeID( MT2, colReaction )
    else:
        colId = versus = ''
    qualifier = {31:' [nubar]', 33:'', 34:' [angular distribution]', 35:' [spectrum]', 40:''}[MF]
    ID = rowId + versus + colId + qualifier

    """ also, create links from covariances to data. If we're trying to point to nonexistant data
    (ie for lumped channels), make a 'fake channel' """
    filename = ''   # should hold file for covariances
    def makeLink( MT, reaction, Id, MAT2 ):
        link_ = gnd.link.link(ENDF_MFMT="%i,%i" % (MF,MT))
        if reaction is not None:
            if MF in (33,40): link_.link = reaction.crossSection
            elif MF==31: link_.path = reaction.outputChannel[0].toXLink() + '/multiplicity'
            elif MF==34: link_.link = reaction.outputChannel[0].distributions['angular']
            elif MF==35:
                product = reaction.outputChannel[0]
                if 'energy' in product.distributions.components:
                    link_.link = product.distributions['energy']
                elif 'uncorrelated' in product.distributions.components:
                    link_.link = product.distributions['uncorrelated'].energyComponent
        elif MAT2:
            # cross-material covariance: get other isotope from the MAT number
            otherTarget = endf_endl.getParticleNameFromMAT( MAT2 )
            quant = gnd.covariances.externalReaction(id=Id, target=otherTarget, ENDF_MFMT=(MF,MT))
            cov_info['externalReactions'][(MAT2,MT,MF)] = quant
            link_.link = quant
        else:
            if 850 < MT < 871: pass
            elif MT==4: cov_info['MTL_2'][(MT,MF)] = zip(range(50,92),[33]*41)
            elif MT==103: cov_info['MTL_2'][(MT,MF)] = zip(range(600,650),[33]*49)
            elif MT==104: cov_info['MTL_2'][(MT,MF)] = zip(range(651,700),[33]*49)
            elif MT==105: cov_info['MTL_2'][(MT,MF)] = zip(range(701,750),[33]*49)
            elif MT==106: cov_info['MTL_2'][(MT,MF)] = zip(range(751,800),[33]*49)
            elif MT==107: cov_info['MTL_2'][(MT,MF)] = zip(range(801,850),[33]*49)
            elif MT==1: cov_info['MTL_2'][(MT,MF)] = zip(range(2,1000),[33]*998)
            if (MT,MF) not in cov_info['lumpedChannels']:
                quant = gnd.covariances.reactionSum(id=Id, ENDF_MFMT=(MF,MT))
                cov_info['lumpedChannels'][(MT,MF)] = quant
            quant = cov_info['lumpedChannels'][(MT,MF)]

            link_.link = quant

        return link_

    rowData = makeLink( MT, rowReaction, rowId, MAT2 )
    rowData.label = 'rowData'
    colData = None
    if (MT2 and MT2!=MT) or MAT2:
        colData = makeLink( MT2, colReaction, colId, MAT2 )
        colData.label = 'columnData'
    return ID, [rowData, colData]

def readMatrix( info, LS,LB, NT,NP, dat ):
    """ matrix format is very similar for MF31, MF33 and MF35, so we can
    generalize parsing the matrix """
    
    nlines = int(math.ceil(NT/6.0))
    subsec = []
    for line in range(nlines): subsec += funkyF(dat.next(), logFile = info.logs)
    # LB flag tells how to interpret this data:
    if LB in (0,1,2,8): # diagonal
        subsec = subsec[:NT]
        energyBounds = [subsec[::2]]
        data = subsec[1::2]
        matrix = []
        for i in range(NP-1): # writing full matrix for now
            matrix.append( [0,] * (NP-1) )
            matrix[-1][i] = data[i]
    elif LB==5:
        if LS==1: #symmetric upper-diagonal
            energyBounds = [subsec[:NP],]
            data = subsec[NP:NT]
            matrix = []
            start, length = 0, NP-1
            for i in range(NP-1):
                matrix.append([0]*i + data[start:start+length])
                start = start+length; length = length-1
            for i in range(NP-1):   # copy to lower-diagonal portion
                for j in range(i,NP-1):
                    matrix[j][i] = matrix[i][j]
        else: # asymmetric, but same energy grids for both axes:
            energyBounds = [subsec[:NP]]
            matrix = subsec[NP:]
            NP = NP-1
            matrix = [matrix[NP*i:NP*i+NP] for i in range(NP)]
    elif LB==6: # asymmetric
        NER = NP
        NEC = (NT-1)//NER
        energyBounds = [subsec[:NER], subsec[NER:NER+NEC]]
        matrix = subsec[NER+NEC:NT]
        NEC-=1; NER-=1  # matrix dimensions
        matrix = [matrix[NEC*i:NEC*i+NEC] for i in range(NER)]
    else:
        return None,None

    axes_ = [0,0,0]
    axes_[0] = gnd.covariances.covarianceAxis( index=0, label="row_energy_bounds", unit="eV",
            interpolation=axes.interpolationXY('linear','flat'), data=energyBounds[0] )
    axes_[1] = gnd.covariances.covarianceAxis( index=1, label="column_energy_bounds", unit="eV",
            interpolation=axes.interpolationXY('linear','flat') )
    if len(energyBounds)==2: axes_[1].data = energyBounds[1]
    else: axes_[1].mirrorOtherAxis = True
    axes_[2] = gnd.covariances.covarianceAxis( index=2, label="matrix_elements", unit="" )

    return matrix, axes_

def readMF31_33( info, dat, mf, mt, cov_info, warningList ):
    """ nubar and cross section covariances have basically the same form
    in ENDF, so we can treat them the same way: """
    
    endfPrecision = 6
    
    dat = myIter(dat)
    # dat contains one MFMT section, but in gnd each cross-correlation gets its own section:
    sectionList, linkData = [], []
    
    ZA,AWR,dum,MTL,dum,NL = funkyFI(dat.next(), logFile = info.logs)
    if MTL!=0:
        # MTL!=0 implies section is a placeholder pointing to a lumped channel
        if (MTL,mf) in cov_info['MTL']:
            cov_info['MTL'][(MTL,mf)].append((mt,mf))
        else: cov_info['MTL'][(MTL,mf)] = [(mt,mf)]
    for subsection in range(NL):
        XMF1,XLFS1,MAT1,MT1,NC,NI = funkyFI(dat.next(), logFile = info.logs)
        if MAT1 == info.MAT:
            # Some files explicitly give MAT1==MAT for internal covariance.
            # for simplicity, just set MAT1=0 unless it is actually a different material
            MAT1 = 0
        XMF1,XLFS1 = int(XMF1),int(XLFS1)

        covarsThisSection = []

        if XMF1!=0 or XLFS1!=0:
            #info.logs.write( "non-zero XMF1/XLFS1!!!" ) # may not be properly dealt with
            raise BadCovariance( "non-zero XMF1/XLFS1 in covariances not currently handled!" )

        for NCdx in range(NC):
            dum,dum,dum,LTY,dum,dum = funkyFI(dat.next(), logFile = info.logs)
            if LTY==0:
                E1,E2,dum,dum,NCI2,NCI = funkyFI(dat.next(), logFile = info.logs)
                subsec = []
                nlines = int(math.ceil(NCI2/6.0))
                for line in range(nlines): subsec += funkyF(dat.next(), logFile = info.logs)
                #coefs = subsec[:NCI2][::2]
                pointerList = []
                for coef,mtnum in zip(subsec[:NCI2][::2],subsec[:NCI2][1::2]):
                    Id, pointers = genID( cov_info, int(mtnum), mf )
                    link = pointers[0]
                    link.attributes['coefficient'] = coef
                    link.label = "summand"
                    pointerList.append( link )
                covarsThisSection.append( gnd.covariances.summedCovariance(
                        lowerBound = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty(E1,'eV'),
                        upperBound = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty(E2,'eV'),
                        pointerList=pointerList) )
                cov_info['NC_data'].append( covarsThisSection[-1] )
            else:
                warningList.append( '       WARNING, non-zero LTY in MF33' )    

        for NIdx in range(NI):
            dum,dum,LS,LB,NT,NP = funkyFI(dat.next(), logFile = info.logs)
            matrix, axes_ = readMatrix( info, LS,LB,NT,NP, dat )
            if LB not in (0,1,2,5,6,8):
                warningList.append( '       WARNING: skipping LB%i section for MF%i MT%i' % ( LB, mf, mt ) )
                continue
            Type='relative'
            matrix = gndMatrix.matrix( matrix, precision=endfPrecision )
            if LS==1: matrix.form='symmetric'
            elif LB in (0,1,2,8):
                matrix.form='diagonal'
                if LB==0:
                    Type='absolute'
                    axes_[-1].unit='b**2'
            elif LB in (5,6): matrix.form='asymmetric'
            
            ENDFconversionFlag = None
            if LB in (2,3,4,8,): ENDFconversionFlag = "LB=%i" % LB

            covarsThisSection.append( gnd.covariances.covarianceMatrix(type=Type, axes=axes_, matrix=matrix,
                ENDFconversionFlag=ENDFconversionFlag) )
        
        # create unique id for each section:
        idNow, pointers = genID( cov_info, mt, mf, MT2=MT1, MF2=(XMF1 or mf), MAT2=MAT1 )
        rowdat, coldat = pointers
        section = gnd.covariances.section( id=idNow, rowData=rowdat, columnData=coldat )

        if len(covarsThisSection)>1:
            form = gnd.covariances.mixedForm( covarsThisSection )
            for idx in range(len(form)): form[idx].index = idx
        elif len(covarsThisSection)==1:
            form = covarsThisSection[0]
        else:
            #raise Exception("Encountered empty covariance section!!!")
            info.logs.write("Missing covariance data from section!")
            continue
        section.nativeData = form.moniker
        section.addForm( form )

        sectionList.append( section )
        linkData.append( (mt,mf,MT1,XMF1, idNow) )
        # end loop over NL subsections

    if dat.index != dat.length: raise BadCovariance("Not all covariance data converted, MF%i MT%i" % (mf,mt))
    return sectionList, linkData

def readMF32( info, dat, mf, mt, cov_info, warningList ) :
    # MF=32 resonance parameter covariances. Must be synchronized with MF=2 resonance parameters

    try: import numpy
    except ImportError:
        warningList.append("       WARNING! Skipping MF32 since NumPy is unavailable")
        return [],[]
    endfPrecision = 6
    resonances = cov_info['resonances']

    def swaprows( matrix, i1, i2, nrows ):
        # matrix rows may be out-of-order and need resorting
        rows = matrix[i1:i1+nrows].copy()
        matrix[i1:i1+nrows] = matrix[i2:i2+nrows]; matrix[i2:i2+nrows] = rows
        cols = matrix[:,i1:i1+nrows].copy()
        matrix[:,i1:i1+nrows] = matrix[:,i2:i2+nrows]; matrix[:,i2:i2+nrows] = cols

    dat = myIter(dat)
    ZA,AWR,dum,dum,NIS,dum = funkyFI(dat.next(), logFile = info.logs)
    if (NIS!=1): raise BadCovariance( "Can't handle multi-isotope file 32!" )
    ZAI,ABN,dum,LFW,NER,dum = funkyFI(dat.next(), logFile = info.logs)

    sections = []
    for subsection in range(NER):
        EL,EH,LRU,LRF,NRO,NAPS = funkyFI(dat.next(), logFile = info.logs)
        if (NRO!=0): raise BadCovariance( "Can't handle non-zero NRO in MF32!" )
        # format is determined mainly by combination of LCOMP and LRU/LRF
        if LRU==1:  # resolved resonance covariance section
            attributes = {'endfConversionFlags':''}
            if LRF==7:
                warningList.append("       WARNING: skipping LRF=7 resonance covariances!")
                continue
            if LRF in (1,2,3):  # Breit-Wigner and simplified Reich-Moore formats are similar
                mf2_elist = zip( resonances.resolved.nativeData.resonanceParameters.getColumn('energy'),
                        resonances.resolved.nativeData.resonanceParameters.getColumn('neutronWidth'),
                        resonances.resolved.nativeData.resonanceParameters.getColumn('captureWidth') )
                SPI,AP,dum,LCOMP,NLS,ISR = funkyFI(dat.next(), logFile = info.logs)
                DAP = []
                if ISR>0:   # scattering radius uncertainty
                    MLS = 1
                    if LRF==3:
                        dum,dum,dum,dum,MLS,one = funkyFI(dat.next(), logFile = info.logs)
                    for idx in range( int( math.ceil(MLS/6.0) ) ):
                        DAP.extend( funkyF(dat.next(), logFile = info.logs) )
                    DAP = [10*dap for dap in DAP[:MLS]] # convert from 10*fm to fm
                if LCOMP==0:
                    # internal correlations given for each resonance, no cross-resonance terms
                    attributes['endfConversionFlags'] = 'LCOMP=0'
                    AWRI, dum, L, dum, tmp, NRS = funkyFI(dat.next(), logFile = info.logs)
                    mf32_resonances, mf32_covars = [],[]
                    for i in range(NRS):
                        mf32_resonances.append( funkyF(dat.next(), logFile = info.logs) )
                        mf32_covars.append( funkyF(dat.next(), logFile = info.logs) + funkyF(dat.next(), logFile = info.logs) )

                    dEsq, dNsq, dNdG, dGsq, dNdF, dGdF, dFsq, dJdN, dJdG, dJdF, dJsq, dum = zip(*mf32_covars)
                    MPAR = 3
                    if any(dFsq): MPAR = 4
                    if any(dJsq): raise BadCovariance("Encountered uncertainty on J in MF32!")
                    matrix = numpy.zeros((MPAR*NRS,MPAR*NRS))
                    for i in range(NRS):
                        matrix[i*MPAR,i*MPAR] = dEsq[i]
                        matrix[i*MPAR+1,i*MPAR+1] = dNsq[i]
                        matrix[i*MPAR+2,i*MPAR+1:i*MPAR+3] = [dNdG[i], dGsq[i]]
                        if MPAR==4:
                            matrix[i*MPAR+3,i*MPAR+1:i*MPAR+4] = [dNdF[i], dGdF[i], dFsq[i]]
                    # symmetrize:
                    for i in range( MPAR*NRS ):
                        matrix[i,i:] = matrix[i:,i]

                    mf32_elist = [(lis[0],lis[3],lis[4]) for lis in mf32_resonances]
                    Type="absolute"
                    GNDmatrix = gndMatrix.matrix( matrix, form="sparse_symmetric" )

                if LCOMP==1:
                    AWRI,dum,dum,dum,NSRS,NLRS = funkyFI(dat.next(), logFile = info.logs)
                    dum,dum,MPAR,dum,tmp,NRB = funkyFI(dat.next(), logFile = info.logs)
                    dim = NRB * MPAR  # num. of resonances * num. parameters per resonance
                    matrixSize = dim * (dim+1) // 2
                    if matrixSize + 6*NRB != tmp:
                        raise BadCovariance("Incorrect dimension for the matrix!")

                    # resonanances are listed again (redundant!):
                    mf32_resonances = [ funkyF(dat.next(), logFile = info.logs) for i in range(NRB) ]
                    if LRF in (1,2):
                        mf32_elist = [(lis[0],lis[3],lis[4]) for lis in mf32_resonances]
                    else:
                        mf32_elist = [(lis[0],lis[2],lis[3]) for lis in mf32_resonances]
                    # then the matrix:
                    data = []
                    nLines, rem = divmod( matrixSize, 6 )
                    for i in range(nLines): data += funkyF(dat.next(), logFile = info.logs)
                    if rem: data += funkyF(dat.next(), logFile = info.logs)[:rem]
                    matrix = numpy.zeros((dim,dim))
                    start, length = 0, dim
                    for i in range(dim):
                        # data stores upper-diagonal matrix. Symmetrize:
                        matrix[i,i:] = matrix[i:,i] = data[start:start+length]
                        start = start+length; length = length-1
                    Type="absolute"
                    GNDmatrix = gndMatrix.matrix( matrix, form="symmetric", precision=endfPrecision )

                elif LCOMP==2:
                    attributes['endfConversionFlags'] = 'LCOMP=2'
                    AWRI, QX, dum, LRX, tmp, NRSA = funkyFI(dat.next(), logFile = info.logs)
                    # resonance parameters + uncertainties:
                    mf32_resonances = [funkyF(dat.next(), logFile = info.logs) for i in range(NRSA*2)]
                    if LRF in (1,2):
                        mf32_elist = [(lis[0],lis[3],lis[4]) for lis in mf32_resonances[::2]]
                    else:
                        mf32_elist = [(lis[0],lis[2],lis[3]) for lis in mf32_resonances[::2]]
                    # for LCOMP==2, off-diagonal terms are given as correlation matrix:
                    dum,dum,NDIGIT,NNN,NM,dum = funkyFI(dat.next(), logFile = info.logs)
                    MPAR = NNN/NRSA
                    matrix = numpy.eye(NNN)
                    diagonal = []
                    for i in range(NRSA):
                        if LRF in (1,2):
                            dE,dum,dum,dGammaN,dGammaG,dGammaF = mf32_resonances[2*i+1]
                            diagonal.extend( [dE,dGammaN,dGammaG] )
                            if MPAR==4: diagonal.append( dGammaF )
                        elif LRF==3:
                            dE,dum,dGammaN,dGammaG,dGammaF1,dGammaF2 = mf32_resonances[2*i+1]
                            diagonal.extend( [dE,dGammaN,dGammaG] )
                            if MPAR==5: diagonal.extend( [dGammaF1,dGammaF2] )
                    if len(diagonal)!=NNN: raise BadCovariance( "Incorrect dimensions for LCOMP=2 matrix!" )

                    # off-diagonal parts of matrix are stored as sparse correlation matrix:
                    attributes['endfConversionFlags'] += ',NDIGIT=%i' % NDIGIT
                    matrix = numpy.eye(NNN, dtype="float") * 10**NDIGIT
                    for i in range(NM):
                        row,col,vals = endfFileToGNDMisc.readEndfINTG( dat.next(), NDIGIT )
                        vals = vals[:row-col]
                        # go to 0-based index:
                        row -= 1
                        col -= 1
                        if row>=NNN or col>=row:
                            raise BadCovariance("Matrix indices out of range for MF32 LCOMP=2 matrix")
                        matrix[row,col:col+len(vals)] = vals
                    for i in range(NNN):    # symmetrize
                        matrix[i,i:] = matrix[i:,i]
                    matrix /= float(10**NDIGIT)

                    # convert correlation -> relative covariance matrix
                    rsd = numpy.sqrt( numpy.array( diagonal ) )
                    matrix = matrix * numpy.outer(rsd,rsd)
                    Type="absolute"
                    GNDmatrix = gndMatrix.matrix( matrix, form="sparse_symmetric" )

                # mf32 may cover only the low-energy portion of the resonance region:
                if len(mf2_elist) != len(mf32_elist):
                    mf2_elist = mf2_elist[:len(mf32_elist)]

                # check whether resonances are sorted by L rather than by energy. If so, rearrange:
                if mf32_elist != mf2_elist or LCOMP==0:
                    if attributes['endfConversionFlags'] == '': attributes['endfConversionFlags'] = 'sortByL'
                    else: attributes['endfConversionFlags'] += ',sortByL'
                if sorted(mf32_elist, key=lambda res: res[0])!=mf2_elist:
                    raise BadCovariance("MF32 resonance parameters don't match MF2 parameters!")

                for i in range(len(mf2_elist)):
                    i2 = mf32_elist.index( mf2_elist[i] )
                    if i2 != i:
                        swaprows( GNDmatrix.data, MPAR*i, MPAR*mf32_elist.index( mf2_elist[i] ), MPAR )
                        # also need to swap values in elist2:
                        val = mf32_elist[i]
                        mf32_elist[i] = mf32_elist[i2]; mf32_elist[i2] = val

                if DAP: # scattering radius uncertainty was specified. Expand matrix to include it:
                    if len(DAP) > 1:
                        raise BadCovariances("Energy-dependent scat. radius uncertainty not yet handled!")
                    dim = len(matrix) + len(DAP)
                    new_matrix = numpy.zeros( (dim,dim) )
                    for i in range(len(DAP)): new_matrix[i,i] = DAP[i]
                    new_matrix[ len(DAP):, len(DAP): ] = matrix
                    GNDmatrix = gndMatrix.matrix( new_matrix, form=GNDmatrix.form,
                            precision=GNDmatrix.precision )

                # switch to diagonal matrix if possible (much more compact):
                if numpy.all( matrix==( numpy.identity(len(matrix)) * matrix.diagonal() ) ):
                    GNDmatrix.form = "diagonal"

                # store into GND:
                resData = resonances.resolved.nativeData.resonanceParameters
                nResonances = len( mf32_elist )
                parametersPerResonance = "energy,neutronWidth,captureWidth"
                if MPAR>=4: parametersPerResonance = parametersPerResonance + ',fissionWidthA'
                if MPAR==5: parametersPerResonance = parametersPerResonance + ',fissionWidthB'
                inputParameters = [ gnd.covariances.loopOverResonanceParameters( nResonances,
                        parametersPerResonance, resData ) ]
                if DAP:
                    inputParameters.insert( 0, gnd.covariances.inputParameter("scatteringRadius",
                        resonances.resolved.nativeData.scatteringRadius, unit="fm") )
                sections.append( gnd.covariances.resonanceParameterCovariance( None, inputParameters,
                    GNDmatrix, type=Type, **attributes ) )
        else:
            # unresolved resonance parameters
            SPI,AP,dum,dum,NLS,dum = funkyFI( dat.next(), logFile = info.logs )
            for lval in range(NLS):
                AWRI,dum,L,dum,tmp,NJS = funkyFI( dat.next(), logFile = info.logs )
                if tmp!=6*NJS: raise BadCovariance( "Incorrect header in MF32 unresolved section!" )
                for jval in range(NJS):
                    D,AJ,GNO,GG,GF,GX = funkyFI( dat.next(), logFile = info.logs )

            # matrix:
            dum,dum,MPAR,dum,tmp,NPAR = funkyFI( dat.next(), logFile = info.logs )
            if tmp != (NPAR*(NPAR+1))/2: raise BadCovariance( "Incorrect header in MF32 unresolved section!" )
            nlines = int(math.ceil(tmp/6.0))
            data = []
            for line in range(nlines): data += funkyF(dat.next(), logFile = info.logs)
            matrix = numpy.zeros((NPAR,NPAR))
            start, length = 0, NPAR
            for i in range(NPAR):
                matrix[i,i:] = matrix[i:,i] = data[start:start+length]
                start = start+length; length = length-1
            if numpy.all( matrix==0 ):
                warningList.append("       WARNING: ignoring empty unresolved covariance matrix!")
            else:
                warningList.append("       WARNING: non-empty unresolved covariance!")

    return sections, []

def readMF34( info, dat, mf, mt, cov_info, warningList ):
    """ angular distribution covariances: """
    endfPrecision = 6
    # dat contains one MFMT section
    dat = myIter(dat)
    sectionList, linkData = [], []
    
    ZA,AWR,dum,LTT,dum,NMT = funkyFI(dat.next(), logFile = info.logs)
    
    form = gnd.covariances.LegendreOrderCovarianceForm()
    groupIndex = 0
    
    for subsection in range(NMT):
        dum,dum,MAT1,MT1,NL,NL1 = funkyFI(dat.next(), logFile = info.logs)
        if MT1 == mt:
            NSS = NL*(NL+1)/2
        else:
            NSS=NL*NL1
        if MAT1 != 0: raise BadCovariance( "Cross material angular distribution covariance is not allowed in ENDF format.  Found MAT1 =",str(MAT1) )
        if MT1 != mt: raise NotImplementedError( "Cross reaction covariances in angular distribution covariance data not supported" )

        for iNSS in range(NSS):
            dum, dum, L, L1, LCT, NI = funkyFI(dat.next(), logFile = info.logs)
            frame = ["frameOfMF4","lab","centerOfMass"][LCT]
            if frame == 'frameOfMF4': 
                warningList.append( "       WARNING: did not look up what frame MF4 data given in, so frame of covariance unknown" )

            Lsection = gnd.covariances.LegendreLValue(L,L1,frame)
            covarsThisL = []

            for NIdx in range(NI):
                dum,dum,LS,LB,NT,NE = funkyFI(dat.next(), logFile = info.logs)
                if not (LS in [ 1, 0 ] and LB<7):
                    #info.logs.write( "Unexpected LS%i LB%i in MF35" % LB )
                    raise BadCovariance( "Unexpected LS%i LB%i in MF34" % LB )
                # pretend it's MF33 file:
                matrix, axes_ = readMatrix( info, LS,LB,NT,NE, dat )
                matrix = gndMatrix.matrix( matrix, precision=endfPrecision )
                if LS==1: matrix.form='symmetric'
                elif LB in (0,1,2,8):
                    matrix.form='diagonal'
                    if LB==0: raise Exception("LB=0 in Legendre covariances")
                else: matrix.form='asymmetric'
                covarsThisL.append( gnd.covariances.covarianceMatrix(type='relative', axes=axes_, matrix=matrix ) )

            if len(covarsThisL)>1:
                sectionForm = gnd.covariances.mixedForm( covarsThisL )
                for idx in range(len(sectionForm)): sectionForm[idx].index = idx
            elif len(covarsThisL)==1:
                sectionForm = covarsThisL[0]

            Lsection.nativeData = sectionForm.moniker
            Lsection.forms[ Lsection.nativeData ] = sectionForm

            # end loop over subsubsections
            form.lvalues.append( Lsection )
       
        if dat.index != dat.length: raise BadCovariance("Not all covariance data converted, MF%i MT%i" % (mf,mt))
        
        # add unique id to the section:
        idNow, pointers = genID( cov_info, mt, mf )
        rowdat, coldat = pointers
        section = gnd.covariances.section( id=idNow, rowData=rowdat, columnData=coldat)
        section.nativeData = form.moniker
        section.addForm( form )
    
        sectionList.append( section )
        linkData.append( (mt,mf,mt,mf, idNow) )
    return sectionList, linkData


def readMF35( info, dat, mf, mt, cov_info, warningList ):
    """ spectra covariances are fairly simple: """
    endfPrecision = 6
    # dat contains one MFMT section
    dat = myIter(dat)
    sectionList, linkData = [], []
    
    ZA,AWR,dum,dum,NK,dum = funkyFI(dat.next(), logFile = info.logs)
    
    form = gnd.covariances.energyIntervalForm()
    groupIndex = 0
    
    for subsection in range(NK):
        E1,E2,LS,LB,NT,NE = funkyFI(dat.next(), logFile = info.logs)
        if not (LS==1 and LB==7):
            #info.logs.write( "Unexpected LS%i LB%i in MF35" % LB )
            raise BadCovariance( "Unexpected LS%i LB%i in MF35" % LB )
        # pretend it's MF33 file:
        LS=1; LB=5
        E1,E2 = [physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty(a,'eV') for a in (E1,E2)]
        
        matrix, axes_ = readMatrix( info, LS, LB, NT, NE, dat )
        matrix = gndMatrix.matrix( matrix, precision=endfPrecision )
        if LS==1: matrix.form='symmetric'
        else: matrix.form='asymmetric'
        covar = gnd.covariances.covarianceMatrix(type='relative', axes=axes_, matrix=matrix, energyBounds=(E1,E2) )
        
        covar.index = groupIndex; groupIndex+=1
        form.addComponent( covar )
        # end loop over NK subsections
    
    if dat.index != dat.length: raise BadCovariance("Not all covariance data converted, MF%i MT%i" % (mf,mt))
    
    # add unique id to the section:
    idNow, pointers = genID( cov_info, mt, mf )
    rowdat, coldat = pointers
    section = gnd.covariances.section( id=idNow, rowData=rowdat, columnData=coldat)
    section.nativeData = form.moniker
    section.addForm( form )

    sectionList.append( section )
    linkData.append( (mt,mf,mt,mf, idNow) )
    return sectionList, linkData

def readMF40(info,dat,mf,mt,cov_info,warningList):
    """ production of radioactive isotopes. Also very similar to MF33 """

    endfPrecision = 6

    dat = myIter(dat)
    sectionList, linkData = [],[]
    ZA,AWR,LIS,dum,NS,dum = funkyFI(dat.next(), logFile = info.logs)
    # each subsection represents different excited state of residual
    for subsection in range(NS):
        try:
            QM,QI,dum,LFS,dum,NL = funkyFI(dat.next(), logFile = info.logs)
        except StopIteration:
            warningList.append('       WARNING: MF40 MT%i lists incorrect number of subsections!' % mt )
            info.doRaise.append( warningList[-1] )
            break
        for subsubsection in range(NL):
            # each subsubsection is a single matrix
            XMF1,XLFS1,MAT1,MT1,NC,NI = funkyFI(dat.next(), logFile = info.logs)
            XMF1,XLFS1 = int(XMF1),int(XLFS1) # XLFS1: level index

            covarsThisSection = []
            if XMF1 not in (0,10):
                raise BadCovariance( "non-zero XMF1/XLFS1 in covariances not currently handled!" )
            if MAT1!=0:
                warningList.append( "       WARNING: cross-material covariance with MAT=%i" % MAT1 )

            for NCdx in range(NC):
                dum,dum,dum,LTY,dum,dum = funkyFI(dat.next(), logFile = info.logs)
                if LTY==0:
                    E1,E2,dum,dum,NCI2,NCI = funkyFI(dat.next(), logFile = info.logs)
                    subsec = []
                    nlines = int(math.ceil(NCI2/6.0))
                    for line in range(nlines): subsec += funkyF(dat.next(), logFile = info.logs)
                    #coefs = subsec[:NCI2][::2]
                    pointerList = [
                            gnd.link.link( 'summand', genID( cov_info, int(mtnum), mf )[0],
                            attributes={'ENDF_MFMT':"%i,%i"%(mf,mtnum), 'coefficient':coef})
                            for coef,mtnum in zip(subsec[:NCI2][::2],subsec[:NCI2][1::2])
                            ]
                    covarsThisSection.append( gnd.covariances.summedCovariance(
                            lowerBound = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty(E1,'eV'),
                            upperBound = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty(E2,'eV'),
                            pointerList=pointerList) )
                else:
                    warningList.append( '       WARNING, non-zero LTY in MF40' )

            for NIdx in range(NI):
                dum,dum,LS,LB,NT,NP = funkyFI(dat.next(), logFile = info.logs)
                matrix, axes_ = readMatrix( info, LS,LB,NT,NP, dat )
                if LB not in (0,1,5,6):
                    warningList.append( '       WARNING: skipping LB%i section for MF%i MT%i' % ( LB, mf, mt ) )
                    continue
                Type='relative'
                matrix = gndMatrix.matrix( matrix, precision=endfPrecision )
                if LS==1: matrix.form='symmetric'
                elif LB in (0,1):
                    matrix.form='diagonal'
                    if LB==0:
                        Type='absolute'
                        axes_[-1].unit='b**2'
                elif LB in (5,6): matrix.form='asymmetric'

                covarsThisSection.append( gnd.covariances.covarianceMatrix(type=Type, axes=axes_, matrix=matrix) )

            if MAT1!=0:
               continue

            # create unique id for each section:
            idNow, pointers = genID( cov_info, mt, mf, MT2=MT1, MF2=(XMF1 or mf), MAT2=MAT1, QI=QI )
            rowdat, coldat = pointers
            section = gnd.covariances.section( id=idNow, rowData=rowdat, columnData=coldat )

            if len(covarsThisSection)>1:
                form = gnd.covariances.mixedForm( covarsThisSection )
                for idx in range(len(form)): form[idx].index = idx
            elif len(covarsThisSection)==1:
                form = covarsThisSection[0]
            else:
                #raise Exception("Encountered empty covariance section!!!")
                info.logs.write("Missing covariance data from section!")
                continue
            section.nativeData = form.moniker
            section.addForm( form )

            sectionList.append( section )
            linkData.append( (mt,mf,MT1,XMF1, idNow) )
            # end loop over NL subsections

    if dat.index != dat.length:
        warningList.append( "Not all covariance data converted for MF%i MT%i" % (mf,mt) )
        info.doRaise.append( warningList[-1] )
    return sectionList, linkData

def fillRemainingProductsResidualForBreakup( info, decayChannel, lightIsotopeNames, breakupProducts, residualZA ) :

    residualZA2 = residualZA
    for productName in lightIsotopeNames :
        if( productName in breakupProducts ) :
            multiplicity = breakupProducts[productName]
            decayChannel.addProduct( endlToGND.newGNDParticle( info, getTypeNameGamma( info, productNameToZA[productName] ), multiplicity = multiplicity ) )
            if( ( residualZA % 1000 ) > 0 ) :
                residualZA2 -= multiplicity * productNameToZA[productName]
            else :
                residualZA2 -= multiplicity * ( 1000 * ( productNameToZA[productName] / 1000 ) )
    if( residualZA2 != 0 ) : decayChannel.addProduct( endlToGND.newGNDParticle( info, getTypeNameGamma( info, residualZA2 ) ) )

def parseReaction( info, target, projectileZA, targetZA, MT, MTData, parseCrossSectionOnly = False ) :
    "Currently, the following logic only supports neutron as projectile."

    warningList, productList = [], []
    info.newIndices( )
    MFKeys = MTData.keys( )
    info.logs.write( '    %3d %s' % ( MT, sorted( MFKeys ) ) )

    for MF in [ 8, 9, 10, 31, 32, 33, 34, 35, 40, 45 ] :
        if( MF in MFKeys ) : MFKeys.remove( MF )

    QM, QI, crossSection, breakupProducts = readMF3( info, MT, MTData[3], warningList )
    MFKeys.remove( 3 )
    if( parseCrossSectionOnly ) : return( QM, QI, crossSection, gnd.channels.NBodyOutputChannel( returnConstantQ( QM ) ), warningList, MFKeys )

    try :
        fissionGenre = { 18 : gnd.channels.fissionGenreTotal, 19 : gnd.channels.fissionGenreFirstChance, 
            20 : gnd.channels.fissionGenreSecondChance, 21 : gnd.channels.fissionGenreThirdChance, 38 : gnd.channels.fissionGenreFourthChance }[MT]
    except :
        fissionGenre = None

    neutronMFs = []
    for MF in [ 4, 5, 6 ] :
        if( MF in MFKeys ) : neutronMFs.append( MF )

    endfMTProductList = endf_endl.endfMTtoC_ProductLists[MT]
    compoundZA = calculateZA( targetZA, projectileZA, minus = False )
    lightIsotopeZAs = sorted( [ productNameToZA[product] for product in lightIsotopeNames ] )
    lightIsotopeZAsMultiplicity = {}
    for product in lightIsotopeNames : lightIsotopeZAsMultiplicity[productNameToZA[product]] = endfMTProductList.productCounts[product]

    if( ( 4 in neutronMFs ) or ( ( MT == 18 ) and ( neutronMFs == [ 5 ] ) ) ) :          # MT == 18 and neutronMFs == [ 5 ] is a special case for bad data.
        ZAP = 1
        if( MT not in [ 2, 18, 19, 20, 21, 38 ] )  :                # Not elastic or fission.
            for product in lightIsotopeNames :
                if( endfMTProductList.productCounts[product] > 0 ) : break
            ZAP = productNameToZA[product]

    levelIndex, decayChannel, twoBodyResidualZA = None, None, None
    if(    50 <= MT < 91 ) :
        levelIndex = MT - 50
    elif( 600 <= MT < 649 ) :
        levelIndex = MT - 600
    elif( 650 <= MT < 699 ) :
        levelIndex = MT - 650
    elif( 700 <= MT < 749 ) :
        levelIndex = MT - 700
    elif( 750 <= MT < 799 ) :
        levelIndex = MT - 750
    elif( 800 <= MT < 849 ) :
        levelIndex = MT - 800
    elif( 875 <= MT < 891 ) :
        levelIndex = MT - 875
    if( MT in [ 91, 649, 699, 749, 799, 849, 891 ] ) : levelIndex = 'c'
    if( MT in [ 103, 104, 105, 106, 107 ] ) : levelIndex = 's'
    level = QM - QI                                                 # If level > 0., residual is in an excited state.
    if( breakupProducts is not None ) :
        if( 50 <= MT < 91 ) :
            level = -QI
        elif( MT == 91 ) :
            pass
        else :
            print breakupProducts
            raise Exception( 'breakup for MT %s is not supported' % MT )
    isUndefinedTwoBody = ( MT == 91 ) or ( 102 < MT <= 107 ) or ( MT in [ 649, 699, 749, 799, 849 ] )
    isTwoBody = ( MT == 2 ) or ( 50 <= MT < 91 ) or ( ( 600 <= MT < 849 ) and not( isUndefinedTwoBody ) )

    if( isTwoBody or isUndefinedTwoBody ) :
        if( MT == 2 ) :
            ZAP = projectileZA
        else :
            for productName in endfMTProductList.productCounts :
                if( endfMTProductList.productCounts[productName] != 0 ) : break
            if( productName == 'gamma' ) :
                ZAP = 0
            else :
                ZAP = productNameToZA[productName]
        twoBodyResidualZA = calculateZA( compoundZA, ZAP )
    undefinedLevelInfo = { 'ZA' : twoBodyResidualZA, 'level' : level, 'levelIndex' : levelIndex, 'count' : 0 }
    if( neutronMFs == [ 4 ] ) :                     # This is a two-body reaction with only angular data.
        if( not( isTwoBody ) ) : raise Exception( 'With only MF = 4 data present, reaction is assumed to be two-body and it is not for MT = %s' % MT )
        product = endlToGND.newGNDParticle( info, getTypeNameGamma( info, ZAP ) )
        component = readMF4( info, product, MT, MTData[4], distributions.angular.twoBodyComponent, warningList )
        MFKeys.remove( 4 )
        product.distributions.setNativeData( component.moniker )
        component.setParent( product.distributions )
        productList.append( product )
    elif( ( neutronMFs == [ 4, 5 ] ) or ( ( neutronMFs == [ 5 ] ) and ZAP == 1 ) ) : 
        # cmattoon, Mar2011: don't check ZAP if MT=5. Currently this combination, MT=5, MF=4/5 appears only for incident gammas
        if( MT!=5 and ZAP != 1 ) : raise Exception( 'ZAP = %d != 1 for MFs = [ 4, 5 ] for MT = %d' % ( ZAP, MT ) )
        multiplicity = 1
        if( MT not in [ 2, 18, 19, 20, 21, 38 ] )  :                # Not elastic or fission.
            for product in lightIsotopeNames :
                if( endfMTProductList.productCounts[product] > 0 ) : break
            ZAP = productNameToZA[product]
            multiplicity = endfMTProductList.productCounts[product]
        else :
            if( MT != 2 ) : multiplicity = -1
        product = endlToGND.newGNDParticle( info, getTypeNameENDF( info, ZAP, undefinedLevelInfo ), multiplicity = multiplicity )

        if( neutronMFs == [ 5 ] ) :
            form = distributions.angular.isotropic( frames[1] )             # MF = 5 data is always in lab frame.
            angularComponent = distributions.angular.component( form.moniker )
            angularComponent.addForm( form )
        else :
            angularComponent = readMF4( info, product, MT, MTData[4], distributions.angular.component, warningList, addAngularComponent = False )
            MFKeys.remove( 4 )

        energyComponent, weights = readMF5( info, MT, MTData[5], warningList, product = product )
        MFKeys.remove( 5 )

        component = distributions.uncorrelated.component( angularComponent, energyComponent )
        angularComponent.setParent( component ); energyComponent.setParent( component )
        component.parent = product.distributions
        product.addDistributionComponent( component )
        product.distributions.setNativeData( component.moniker )
        productList.append( product )
    elif( neutronMFs == [ 6 ] ) :
        isTwoBody = readMF6( MT, info, MTData[6], productList, warningList, undefinedLevelInfo, isTwoBody )
        MFKeys.remove( 6 )
    elif( neutronMFs == [] ) : 
        if( isTwoBody and False ) :                 # ????????? Why False
            raise Exception( 'How did we get here.' )
            product = endlToGND.newGNDParticle( info, getTypeNameGamma( info, ZAP ) )
            residualZA = calculateZA( compoundZA, ZAP )
            if( not( levelIndex is None ) ) :
                if( ( levelIndex <= info.targetLevel ) and ( info.target.getZ_A_SuffixAndZA( )[3] == residualZA ) ) : levelIndex -= 1
            if( QI != QM ) :     # Residual is in an excited state.
                decayChannel = gnd.channels.NBodyDecayChannel( )
                decayChannel.addProduct( endlToGND.newGNDParticle( info, getTypeNameGamma( info, residualZA ) ) )
            residual = endlToGND.newGNDParticle( info, getTypeNameGamma( info, residualZA, level = level ), decayChannel = decayChannel )
            productList.append( product )
            productList.append( residual )
    else :
        pass

    info.crossSection = crossSection
    readMF12_13( info, MT, MTData, productList, warningList )
    del info.crossSection

    for MF in [ 12, 13, 14, 15 ] :
        if( MF in MFKeys ) : MFKeys.remove( MF )
    specialBe9n2nExcited, specialBe9n2nExcitedLevel = False, 0
    if( 875 <= MT < 892 ) :                                         # Special case for (z,2n[?])
        if( targetZA == 4009 ) :                                    # Special case for Be9(n,2n[?]He4)He4
            specialBe9n2nExcited = True
            specialBe9n2nExcitedLevel = ( QM - QI ) / 1e6

    if( MT == 5 ) :
        if( QM != QI ) : info.logs.write( '    --QM %s != QI = %s\n' % ( QM, QI ) )
        outputChannel = gnd.channels.sumOfRemainingOutputChannels( returnConstantQ( QM ) )
    elif( ( MT == 102 ) and not( isTwoBody ) ) :
        residualIndex, gammaMissing = -1, False
        for index, product in enumerate( productList ) :
            if( product.getName( ) != 'gamma' ) : residualIndex = index
            gammaMissing = ( product.getName( ) == 'gamma' ) or gammaMissing
        if( residualIndex == -1 ) : 
            productList.insert( 0, endlToGND.newGNDParticle( info, getTypeNameENDF( info, calculateZA( compoundZA, 0 ), undefinedLevelInfo ) ) )
        if( residualIndex > 0 ) : productList.insert( 0, productList.pop( residualIndex ) )
        if( not( gammaMissing ) ) : productList.append( endlToGND.newGNDParticle( info, getTypeNameENDF( info, 0, undefinedLevelInfo ) ) )
        outputChannel = gnd.channels.NBodyOutputChannel( returnConstantQ( QM ) )        # Q????? What about QI?
    elif( isTwoBody ) :
        if( ( QI == 0 ) and ( QM != 0 ) ) : raise Exception("QI = 0, QM = %f for MT=%i" % (QM,MT))
        outputChannel = gnd.channels.twoBodyOutputChannel( returnConstantQ( QI ) )
        if( len( productList ) == 0 ) :
            for ZA in lightIsotopeZAs :
                if( lightIsotopeZAsMultiplicity[ZA] != 0 ) :
                    productList.append( endlToGND.newGNDParticle( info, getTypeNameENDF( info, ZA, undefinedLevelInfo ) ) )
                    break
            if( len( productList ) == 0 ) :
                if( MT != 2 ) : raise Exception( "product data for reaction MT = %s needs to be implemented" % MT )
                productList.append( endlToGND.newGNDParticle( info, getTypeNameENDF( info, projectileZA, undefinedLevelInfo ) ) )
        decayProductList = productList[1:]
        productList = productList[:1]                               # Assume first product is "b" in "a + A -> b + B" where B is the larger product.
        ZA = productList[0].particle.getZ_A_SuffixAndZA( )[3]
        residualZA = calculateZA( compoundZA, ZA )
        levelIndex = undefinedLevelInfo['levelIndex']
        if( not( levelIndex is None ) ) :
            if( ( levelIndex <= info.targetLevel ) and ( info.target.getZ_A_SuffixAndZA( )[3] == residualZA ) ) : levelIndex -= 1
        undefinedLevelInfo['levelIndex'] = levelIndex
        for index, product in enumerate( decayProductList ) :
            ZA = product.particle.getZ_A_SuffixAndZA( )[3]
            if( residualZA == ZA ) :
                productList.append( decayProductList.pop( index ) )
                break
        if( len( productList ) < 2 ) : 
            if( MT == 2 ) :
                productList.append( endlToGND.newGNDParticle( info, target ) )
            else :
                if( ZA == undefinedLevelInfo['ZA'] ) : undefinedLevelInfo['ZA'] = None
                productList.append( endlToGND.newGNDParticle( info, getTypeNameENDF( info, residualZA, undefinedLevelInfo ) ) )
        decayZAs, decayGammaList, decayNonGammaList = 0, [], []
        for decayProduct in decayProductList : 
            decayZAs += decayProduct.particle.getZ_A_SuffixAndZA( )[-1]
            if( decayProduct.getName( ) == 'gamma' ) :
                decayGammaList.append( decayProduct )
            else :
                decayNonGammaList.append( decayProduct )
        if( decayZAs == 0 ) :
            if( len( decayGammaList ) != 0 ) :
                if( len( decayNonGammaList ) == 0 ) : 
                    decayNonGammaList.append( endlToGND.newGNDParticle( info, getTypeNameENDF( info, residualZA, None ) ) )
            elif( len( decayNonGammaList ) != 0 ) : 
                if( len( decayGammaList ) == 0 ) : decayGammaList.append( endlToGND.newGNDParticle( info, getTypeNameENDF( info, 0, None ) ) )
            decayProductList = decayNonGammaList + decayGammaList
        else :
            raise Exception( "decayZAs = %d != 0" % decayZAs )

        if( breakupProducts is not None ) :
            if( decayChannel is not None ) : raise Exception( 'breakupProducts and decayChannel both not None' )
            decayChannel = gnd.channels.NBodyDecayChannel( returnConstantQ( QM ) )
            fillRemainingProductsResidualForBreakup( info, decayChannel, lightIsotopeNames, breakupProducts, 
                productList[1].particle.getZ_A_SuffixAndZA( )[3] )
            productList[1].addDecayChannel( decayChannel )
        elif( len( decayProductList ) > 0 ) :                         # At this point, both two bodies are in productList and second one is redisual.
            if( QI >= QM ) : raise Exception( "Negative decay Q-value for MT%i" % MT )
            decayChannel = gnd.channels.NBodyDecayChannel( returnConstantQ( QM - QI ) )     # Q????? Not right.
            for decayProduct in decayProductList : decayChannel.addProduct( decayProduct )
            productList[1].addDecayChannel( decayChannel )

    elif( endfMTProductList.isFission ) :
        if hasattr( info, 'fissionEnergies' ): 
            ER = info.fissionEnergies.forms[info.fissionEnergies.nativeData].data[ 'nonNeutrinoEnergy' ] # gets whole polynomial + uncertianty
            useThisQM = ER[ 0 ][ 0 ] # grab constant term and take the value, not the uncertainty
            if abs( QM - useThisQM )/QM > 1e-7: warningList.append( "WARNING: Fission QM inconsistent with energy release data for MT = " + str( MT ) )
        else: useThisQM = QM
        outputChannel = gnd.channels.fissionChannel( returnConstantQ( useThisQM ), fissionGenre = fissionGenre )
        if( MT == 18 ) :
            if( hasattr( info, 'fissionEnergies' ) ) : outputChannel.addFissionEnergyReleased( info.fissionEnergies )
            if( len( productList ) > 0 ) : outputChannel.addProduct( productList.pop( 0 ) )
            if( len( outputChannel ) == 0 ) :
                multiplicity = gnd.productData.multiplicity.unknown( )
                product = endlToGND.newGNDParticle( info, getTypeNameGamma( info, 1 ), multiplicity = multiplicity )
                outputChannel.addProduct( product )
            else :
                for product in outputChannel :
                    if( product.getName( ) == 'n' ) : break
                info.firstFissionNeutron = product
                if( 'prompt' in info.totalOrPromptFissionNeutrons ) :
                    product.setMultiplicity( info.totalOrPromptFissionNeutrons['prompt'] )
                    product.multiplicity.removeForm('constant')
                    product.addAttribute( 'emissionMode', 'prompt' )
                    if( 'total' in info.totalOrPromptFissionNeutrons ) :
                        warningList.append( '       WARNING: have prompt fission nu_bar so not including total' )
                elif( 'total' in info.totalOrPromptFissionNeutrons ) :
                    product.setMultiplicity( info.totalOrPromptFissionNeutrons['total'] )
                    product.multiplicity.removeForm('constant')
                    product.addAttribute( 'emissionMode', 'total' )
                if( hasattr( info, 'delayedFissionDecayChannel' ) ) :
                    for delayedNeutron in info.delayedFissionDecayChannel : outputChannel.addProduct( delayedNeutron )
        else :
            if( neutronMFs == [] ) :
                if( hasattr( info, 'firstFissionNeutron' ) ) :
                    multiplicity = gnd.productData.multiplicity.reference( info.firstFissionNeutron )
                else :                                              # When singleMTOnly is fission MT != 18.
                    multiplicity = gnd.productData.multiplicity.unknown( )
                product = endlToGND.newGNDParticle( info, getTypeNameGamma( info, 1 ), multiplicity = multiplicity )
                if( hasattr( info, 'firstFissionNeutron' ) ) :
                    component = distributions.base.referenceComponent( info.firstFissionNeutron.distributions )
                    product.addDistributionComponent( component )
                    product.distributions.setNativeData( component.moniker )
                outputChannel.addProduct( product )

        # July 2011: some files have distributions for 1stChanceFission etc, but should still link to total nubar:
        for product in productList:
            if( ( product.multiplicity.nativeData == gnd.tokens.constantFormToken ) and ( product.multiplicity.getConstant() == -1 ) ) :
                if hasattr( info, 'firstFissionNeutron' ):
                    multiplicity = gnd.productData.multiplicity.reference( info.firstFissionNeutron )
                    product.multiplicity.addFormAsNativeData( multiplicity )

        while( len( productList ) > 0 ) : outputChannel.addProduct( productList.pop( 0 ) )
    else :
        outputChannel = gnd.channels.NBodyOutputChannel( returnConstantQ( QI ) )                    # Q?????
        if( MT not in [ 18, 19, 20, 21, 38 ] ) :
            residualZA, ZAsMultiplicities, productAsResidual, biggestProduct = compoundZA, {}, None, 0
            for index, product in enumerate( productList ) :
                if( product.getName( ) == 'gamma' ) : continue
                Z, A, suffix, ZA = product.particle.getZ_A_SuffixAndZA( )
                try :
                    m = product.multiplicity.getConstant( )
                except :
                    info.logs.write( '\n\nIncorrect multiplicity in ENDF file! MT = %s\n' % MT )
                    info.logs.write( 'Multiplicity should be constant but is (%s).\n' % product.multiplicity.nativeData )
                    info.logs.write( '%s\n' % product.multiplicity.getFormByToken( product.multiplicity.nativeData ) )
                    raise
                if( ZA in lightIsotopeZAs ) :
                    residualZA = calculateZA( residualZA, m * ZA, minus = True )
                        # If we have different distributions for both neutrons in (n,2n), n shows up twice in the productList.
                    if ZA in ZAsMultiplicities: ZAsMultiplicities[ZA] += m
                    else: ZAsMultiplicities[ZA] = m
                else :
                    if( not( productAsResidual is None ) ) : 
                        raise Exception( 'multiple residuals for MT = %, %s %s' % ( MT, productAsResidual.getToken( ), product.getToken( ) ) )
                    productAsResidual = product

            _residualZA = compoundZA
            for ZA in lightIsotopeZAsMultiplicity : _residualZA = calculateZA( _residualZA, lightIsotopeZAsMultiplicity[ZA] * ZA, minus = True )
            if( residualZA != 0 ) :
                for ZA in lightIsotopeZAs :
                    if( ZA not in ZAsMultiplicities ) : ZAsMultiplicities[ZA] = 0
                    if( ZAsMultiplicities[ZA] == lightIsotopeZAsMultiplicity[ZA] ) : continue       # All this ZA accounted for.
                    if( ZAsMultiplicities[ZA] > lightIsotopeZAsMultiplicity[ZA] ) :
                        raise Exception( 'negative multiplicity for ZA = %s for MT = %s' % ( ZA, MT ) )
                    multiplicity = lightIsotopeZAsMultiplicity[ZA] - ZAsMultiplicities[ZA]
                    productList.append( endlToGND.newGNDParticle( info, getTypeNameENDF( info, ZA, None ), multiplicity = multiplicity ) )
                    residualZA = calculateZA( residualZA, multiplicity * ZA, minus = True )
                if( productAsResidual is None ) :
                    if( residualZA > 0 ) : productList.append( endlToGND.newGNDParticle( info, getTypeNameENDF( info, residualZA, undefinedLevelInfo ) ) )

            if( MT in [ 103, 104, 105, 106, 107, 91, 649, 699, 749, 799, 849, 891 ] ) :
                gammaIndices = []
                for index, product in enumerate( productList ) :
                    if( product.getName( ) == 'gamma' ) : gammaIndices.append( index )
                if( len( gammaIndices ) > 0 ) :
                    decayChannel = gnd.channels.NBodyDecayChannel( returnConstantQ( -QI ) )          # Q?????
                    finalResidual = endlToGND.newGNDParticle( info, getTypeNameENDF( info, residualZA, None ) )
                    decayChannel.addProduct( finalResidual )
                    if( productList[-1].getAttribute( 'ENDFconversionFlag' ) ) : finalResidual.addAttribute( 'ENDFconversionFlag', productList[-1].getAttribute( 'ENDFconversionFlag' ) )
                    for index in gammaIndices : decayChannel.addProduct( productList[index] )
                    gammaIndices.reverse( )
                    for index in gammaIndices : del productList[index]
                    productList[-1].addDecayChannel( decayChannel )
                    if( productList[-1].getDistributionNativeData( ) != 'none' ) :  # One must be careful here as this assumes that distributions
                        u_distributions = productList[-1].distributions             # can be simply traded.
                        productList[-1].distributions = finalResidual.distributions
                        finalResidual.distributions = u_distributions
            if( breakupProducts is not None ) : 
                if( MT == 91 ) :
                    if( decayChannel is not None ) : raise Exception( 'breakupProducts and decayChannel both not None' )
                    decayChannel = gnd.channels.NBodyDecayChannel( returnConstantQ( QM ) )
                    fillRemainingProductsResidualForBreakup( info, decayChannel, lightIsotopeNames, breakupProducts, 
                        productList[1].particle.getZ_A_SuffixAndZA( )[3] )
                    productList[1].addDecayChannel( decayChannel )
                else :
                    raise Exception( 'breakup not supported for MT %d' % MT )

    for product in productList : outputChannel.addProduct( product )

    return( QM, QI, crossSection, outputChannel, warningList, MFKeys )

def parseCovariances( info, MTDatas, MTdict, singleMTOnly=None, resonances=None ):

    covarianceSuite = gnd.covariances.covarianceSuite()
    linkData = []    # which mf/mts need to be linked for each covariance?
    if( singleMTOnly ) : return( covarianceSuite, linkData )

    # make list of available covariance information:
    warningList = []
    cov_info = {'MTL':{}, 'MTL_2':{}, 'lumpedChannels':{}, 'externalReactions':{}, 'mfmts':[], 'MTdict':MTdict,
            'resonances':resonances, 'NC_data':[]}
    for mt in MTDatas.keys():
        if (singleMTOnly is not None) and (mt!=singleMTOnly): continue
        for mf in MTDatas[mt].keys():
            if mf >= 30:    # all covariance-type data
                cov_info['mfmts'].append((mf,mt))
    cov_info['mfmts'].sort()  # sorting first by MF, then by MT

    for mf,mt in cov_info['mfmts']:
        try:
            if mf in (31,33):
                covars, tmp = readMF31_33( info, MTDatas[mt][mf], mf, mt, cov_info, warningList )
            elif mf == 32:
                covars, tmp = readMF32( info, MTDatas[151][32], mf, mt, cov_info, warningList )
            elif mf == 34:
                covars, tmp = readMF34( info, MTDatas[mt][mf], mf, mt, cov_info, warningList )
            elif mf == 35:
                covars, tmp = readMF35( info, MTDatas[mt][mf], mf, mt, cov_info, warningList )
            elif mf == 40:
                covars, tmp = readMF40( info, MTDatas[mt][mf], mf, mt, cov_info, warningList )
            else:
                warningList.append( '       WARNING: MF%i not yet supported' % mf)
                continue
            for cov in covars:
                if mf == 32: covarianceSuite.addModelParameterCovariance(cov)
                else: covarianceSuite.addSection(cov)
            linkData += tmp
        except BadCovariance, e:
            warningList.append('       WARNING: MF%i MT%i covariance conversion failed with message "%s"' % (mf,mt,e) )
            info.doRaise.append( warningList[-1] )

    # fix links for summed matrices:
    for summedMatrix in cov_info['NC_data']:
        for pointer in summedMatrix.pointerList:
            pointed_to = [sec for sec in covarianceSuite.sections if sec.columnData is None
                    and sec.rowData['ENDF_MFMT'] == pointer['ENDF_MFMT']]
            if len(pointed_to) != 1:
                raise Exception("Can't resolve links for summed covariance matrix!")
            pointer.link = pointed_to[0]

    # fix lumped channel covariances (MT851-871) and summed channels (MT1,4,103-107)
    summedReactions = cov_info['MTL'].copy();  summedReactions.update( cov_info['MTL_2'] )
    for (mt,mf) in summedReactions:
        lumpedChannels = cov_info['lumpedChannels'][(mt,mf)]
        for (mt2,mf2) in summedReactions[(mt,mf)]:
            if mt not in range(851,872) and mt2 not in cov_info['MTdict']: continue
            label, pointers = genID(cov_info, mt2, mf2)
            lumpedChannels.reactions.append( pointers[0] )
    for lc in sorted(cov_info['lumpedChannels']):
        covarianceSuite.addReactionSum( cov_info['lumpedChannels'][lc] )
    for exReac in sorted(cov_info['externalReactions']):
        covarianceSuite.addExternalReaction( cov_info['externalReactions'][exReac] )

    # add labels
    for label,section in enumerate(covarianceSuite.sections + covarianceSuite.modelParameterCovariances):
        section.label = str(label)

    sys.stdout.flush( )
    for warning in warningList : info.logs.write( warning + '\n', stderrWriting = True )
    return covarianceSuite, linkData

def endfFileToGND( fileName, xenslIsotopes = None, useFilesQAlways = True, singleMTOnly = None, parseCrossSectionOnly = False, 
        toStdOut = True, toStdErr = True, logFile = None, skipBadData = False, deprecatedOptions = {} ) :

    logs = logFiles( toStdOut = toStdOut, toStdErr = toStdErr, logFile = logFile, defaultIsStderrWriting = False )
    header, MAT, MTDatas = endfFileToGNDMisc.parseENDFByMT_MF( fileName, logFile = logs )

    targetZA, targetMass, LRP, LFI, NLIB, NMOD = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( MTDatas[451][1][0], logFile = logs )
    targetZA = int( targetZA )      # Target's ZA
    LRP = int( LRP )            # Resonance parameter data info
    LFI = int( LFI )            # Is fission present
    NLIB = int( NLIB )          # What library (e.g., 0 = ENDF/B
    NMOD = int( NMOD )          # Version modification flag
    isNaturalTarget = ( targetZA % 1000 ) == 0

    targetExcitationEnergy, STA, LIS, LISO, dummy, NFOR = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( MTDatas[451][1][1], logFile = logs )
    STA = int( STA )            # Is nucleus unstable
    LIS = int( LIS )            # Excitation number
    LISO = int( LISO )          # Isomeric state number
    NFOR = int( NFOR )          # Must be 6 for ENDF/B6 format

    projectileMass, dummy, LREL, dummy, NSUB, NVER = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( MTDatas[451][1][2], logFile = logs )
    NSUB = int( NSUB )          # 10 * ZA + iType for projectile
    NVER = int( NVER )          # Evaluation version number
    LREL = int( LREL )          # Evaluation sub-version number
    projectileZA, ITYPE = NSUB / 10, NSUB % 10
    #if( ITYPE != 0 ) : raise Exception( "Currently, only ITYPE = 0 files are supported" )
    #if( projectileZA != 1 ) : raise Exception( "Currently, only neutron as projectile is supported" )

    targetTemperature, dummy, LDRZ, dummy, NWD, NXC = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( MTDatas[451][1][3], logFile = logs )
    LDRZ = int( LDRZ )          # Primary or special evaluation of this material
    NWD = int( NWD )            # 
    NXC = int( NXC )            #

    #evaluation = MTDatas[451][1][6][4:22].strip( )
    evaluation = {\
        0: "ENDF/B", \
        1: "ENDF/A", \
        2: "JEFF", \
        3: "EFF",\
        4: "ENDF/B (HE)",\
        5: "CENDL", \
        6: "JENDL", \
        21: "SG-23",\
        31: "INDL/V", \
        32: "INDL/A", \
        33: "FENDL",\
        34: "IRDF", \
        35: "BROND (IAEA version)",\
        36: "INGDB-90",\
        37: "FENDL/A",\
        41: "BROND",\
    }.get( NLIB, 'Unknown' )
    documentation = gnd.documentation.documentation( 'endfDoc', '\n'.join( MTDatas[451][1][4:4+NWD] ) )
    evaluatedStyle = gnd.miscellaneous.style( 'evaluated', attributes = { 'library' : evaluation, 'version' : "%i.%i.%i" % (NVER,LREL,NMOD) } )
    transportables = [ 'n', 'gamma' ]                          # ???? Check this with Gerry?
    info = endlToGND.infos( xenslIsotopes, transportables = transportables )
    info.ignoreBadNK14 = False
    deprecatedOptionList = [ 'ignoreBadNK14' ]
    for options in deprecatedOptions :
        if( options not in deprecatedOptionList ) : raise Exception( 'invalid deprecated option "%s"' % options )
        setattr( info, options, deprecatedOptions[options] )
    info.doRaise = []
    info.logs = logs
    info.projectile = { 0 : 'g', 1 :  'n', 1001 : 'H1', 1002 : 'H2', 1003 : 'H3', 2003 : 'He3', 2004 : 'He4' }[projectileZA]
    info.ZA_AWRMasses = {}
    info.ZA_AWRMasses[projectileZA] = { projectileMass : 1 }
    info.ZA_AWRMasses[targetZA] = { targetMass : 1 }
    info.ZAMasses[projectileZA] = info.masses.getMassFromZA( projectileZA )
    info.ZAMasses[targetZA] = targetMass * info.masses.getMassFromZA( 1 )
    info.ZAMasses[1] = info.masses.getMassFromZA( 1 ) # always need neutron mass in table
    
    info.sumCrossSections = { 4 : { 'total' : None, 'MTs' : range( 50, 91 ), 'partialPresent' : False }, 
        103 : { 'total' : None, 'MTs' : range( 600, 650 ), 'partialPresent' : False },
        104 : { 'total' : None, 'MTs' : range( 650, 700 ), 'partialPresent' : False },
        105 : { 'total' : None, 'MTs' : range( 700, 750 ), 'partialPresent' : False },
        106 : { 'total' : None, 'MTs' : range( 750, 800 ), 'partialPresent' : False },
        107 : { 'total' : None, 'MTs' : range( 800, 850 ), 'partialPresent' : False } }
    info.MF12_LO2 = {}

    try: Date = endfFileToGNDMisc.getENDFDate( MTDatas[451][1][4][22:33] )
    except Exception as e:
        info.doRaise.append( str(e) )
        import datetime
        Date = datetime.datetime.today().strftime("%Y-%m-%d")
    author = MTDatas[451][1][4][33:66]

    projectile = getTypeNameGamma( info, projectileZA )
    level = targetExcitationEnergy
    levelIndex = None
    if( LIS  != 0 ) : levelIndex = LIS
    target = getTypeNameGamma( info, targetZA, level = level, levelIndex = levelIndex )
    if STA!=0 and not isinstance(target,gnd.xParticle.nuclearLevel): target.attributes['unstable'] = True
    reactionSuite = gnd.reactionSuite.reactionSuite( projectile, target, particleList = info.particleList, style = evaluatedStyle, documentation = documentation )
    reactionSuite.setAttribute( 'temperature', physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( targetTemperature, 'K' ) )
    if( LISO != 0 ) :
        metaStableName = target.getName( ).split( '_' )[0]
        metaStableName += '_m%d' % LISO
        reactionSuite.addAlias( metaStableName, target.getName( ) )
    info.setReactionSuite( reactionSuite )
    info.target = reactionSuite.target
    info.targetLevel = LIS

    MTDatas[451][1] = MTDatas[451][1][:4+NWD]
    warningList = []

    info.totalOrPromptFissionNeutrons = {}
    if( 452 in MTDatas ) :
        info.totalOrPromptFissionNeutrons['total'] = getTotalOrPromptFission( info, MTDatas[452][1], 'total', warningList )
        #MTDatas.pop( 452 ) # don't remove these yet, still need the covariance info
    if( 455 in MTDatas ) :
        info.indices.delayedFissionNeutron = True
        info.delayedFissionDecayChannel = getDelayedFission( info, MTDatas[455], warningList )
        info.indices.delayedFissionNeutron = False
        #MTDatas.pop( 455 )
    if( 456 in MTDatas ) :
        info.totalOrPromptFissionNeutrons['prompt'] = getTotalOrPromptFission( info, MTDatas[456][1], 'prompt', warningList )
        #MTDatas.pop( 456 )
    if( 458 in MTDatas ) :
        info.fissionEnergies = getFissionEnergies( info, MTDatas[458] )
        #MTDatas.pop( 458 )
    if ( 454 in MTDatas ) : 
        info.independentFissionYields = readMF8( info, 454, MTDatas[454], warningList )
    if ( 459 in MTDatas ) : 
        info.cumulativeFissionYields = readMF8( info, 459, MTDatas[459], warningList )
    sys.stdout.flush( )
    for warning in warningList : logs.write( warning + '\n', stderrWriting = True )

    if( 3 in MTDatas ) :
        # MT3 is a redundant quantity ('nonelastic'). Remove any product information, and store only the cross section
        removed = False
        for mf in range(4,16):
            if mf in MTDatas[3]:
                MTDatas[3].pop(mf)
                removed = True
        if removed:
            logs.write( "       WARNING: distributions for MT=3 (nonelastic) are not supported and have been ignored\n", stderrWriting = True )

    MTList = endfFileToGNDMisc.niceSortOfMTs( MTDatas.keys( ), verbose = False, logFile = logs )
    haveTotalFission = (18 in MTList)
    fissionMTs = [mt for mt in MTList if mt in (19,20,21,38)]
    iChannel, totalReactions, MT5Reaction, summedReactions, fissionComponents, productions = 0, {}, None, [], [], []
    for MT in MTList :
        if( not( singleMTOnly is None ) ) and ( MT != singleMTOnly ) : continue
        MTData = MTDatas[MT]

        if 3 in MTData: # normal reaction, with cross section and distributions
            QM, QI, crossSection, outputChannel, warningList, MFKeys = parseReaction( info, target, projectileZA,
                    targetZA, MT, MTData, parseCrossSectionOnly = parseCrossSectionOnly )
            logs.write( '\n' )
            sys.stdout.flush( )
            for warning in warningList : logs.write( warning + '\n', stderrWriting = True )
            if( len( MFKeys ) ) : logs.write( '       WARNING: For reaction MT = %d, the following MFs were not converted: %s\n' % ( MT, MFKeys ) )
            if( outputChannel is None ) : break
            reactionOutputChannel = outputChannel

            reaction = gnd.reaction.reaction( reactionOutputChannel, "%s" % iChannel, ENDF_MT = MT, crossSection = crossSection )
            reaction.setAttribute( 'date', Date )
#            if( useFilesQAlways or ( MT in [ 18, 19, 20, 21, 38 ] ) or ( isNaturalTarget and ( QM != 0. ) ) ) : 
#                reaction.setAttribute( 'Q', physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( QM, 'eV' ) )  # Q????? The meaning of Q needs to be clarified.
            if( MT in info.sumCrossSections ) :
                if( MT not in totalReactions ) : totalReactions[MT] = []
                totalReactions[MT].append( reaction )
            elif( MT == 5 ) :
                MT5Reaction = reaction
            elif MT in (1,3):
                channelId = {1:'total', 3:'nonelastic'}
                reaction = gnd.summedReaction.summedReaction( channelId[MT], "%s" % iChannel, ENDF_MT = MT, crossSection = crossSection, Q = reactionOutputChannel.Q )
                reaction.setAttribute( 'date', Date )
                summedReactions.append( reaction )
            elif MT in fissionMTs and haveTotalFission: # this is 1st, 2nd, etc fission but total is also present
                fissionComponents.append( reaction )
            else :
                reactionSuite.addReaction( reaction, iChannel )
                iChannel += 1
        
        # get radioactive production data (if any) from MF 8-10. x-section form depends on value of LMF:
        # LMF=3 says just use MF3, LMF=9 says multiply MF3 by MF9 weights, LMF=10 explicitly gives x-section
        LMF, radioactiveDatas = readMF8( info, MT, MTData, warningList )
        if LMF:
            for radioactiveData in radioactiveDatas:
                if radioactiveData[-1] == None:
                    Q = None
                else:
                    Q = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( radioactiveData[-1], 'eV' )
                    Q = gnd.channelData.Q.component( gnd.channelData.Q.constant( Q ) )
                productionCrossSection = gnd.reactionData.crossSection.component()
                if LMF==3:
                    Q = outputChannel.Q
                    axes_ = axes.defaultAxes(labelsUnits={0:('energy_in','eV'), 1:('weight','')})
                    axes_[0].frame, axes_[1].frame = axes.labToken, axes.labToken
                    weights = XYs.XYs( axes_, (crossSection.getDomain(), (1,1)), ENDF_Accuracy,
                            moniker="weight", isPrimaryXData=True, dataForm="XsAndYs" )
                    productionCrossSection.addForm( gnd.reactionData.crossSection.weightedPointwise( crossSection, weights ) )
                elif LMF==6:
                    # product yield = MF6 multiplicity * MF3 xsc
                    Q = outputChannel.Q
                    prod = [p for p in outputChannel if
                            p.particle.getZ_A_SuffixAndZA()[-1] == radioactiveData[0] and
                            p.particle.getLevelIndex() == radioactiveData[4]]
                    if len( prod ) == 1:
                        axes_ = axes.defaultAxes(labelsUnits={0:('energy_in','eV'), 1:('weight','')})
                        axes_[0].frame, axes_[1].frame = axes.labToken, axes.labToken
                        weights = XYs.XYs( axes_, prod[0].multiplicity.forms[ prod[0].multiplicity.nativeData ].copy(),
                                ENDF_Accuracy, moniker="weight", isPrimaryXData=True, dataForm="XYs" )
                        productionCrossSection.addForm( gnd.reactionData.crossSection.weightedPointwise(
                            crossSection, weights ) )
                    else:
                        raise Exception( 'Cannot find unique product in MF6 corresponding to MT%i LMF6 data!' % MT )
                elif LMF==9:
                    radioactiveWeight = radioactiveData[3]
                    if len( radioactiveWeight ) == 1:
                        weights = XYs.XYs( radioactiveWeight[0].axes, radioactiveWeight[0], ENDF_Accuracy,
                                moniker="weight", isPrimaryXData=True )
                        weights.axes[0].frame, weights.axes[1].frame = axes.labToken, axes.labToken
                        productionCrossSection.addForm( gnd.reactionData.crossSection.weightedPointwise(
                            crossSection, weights ) )
                    else:
                        raise Exception( 'Piecewise weighted cross section (MF=9) not supported: MT=%s' % MT )
                elif LMF==10:
                    productionCrossSection.addForm( radioactiveData[2] )
                else:
                    raise Exception("Unknown LMF=%i encountered in MF=8 for MT=%i" % (LMF,MT))
                production = gnd.production.production( radioactiveData[1], label = -1, ENDF_MT = MT, crossSection = productionCrossSection, Q = Q )
                productions.append( production )

    # end loop over MT sections

    totalReactionMTs = endfFileToGNDMisc.niceSortOfMTs( totalReactions.keys( ), verbose = False, logFile = logs )
    for MT in totalReactionMTs :
        if( not( info.sumCrossSections[MT]['partialPresent'] ) ) :
            for totalReaction in totalReactions[MT] :
                totalReaction.setLabel( "%s" % iChannel )
                reactionSuite.addReaction( totalReaction, iChannel )
                iChannel += 1
        else:
            for totalReaction in totalReactions[MT] :
                channelId = {4:'(z,n)', 103:'(z,p)', 104:'(z,d)', 105: '(z,t)', 106:'(z,He3)', 107:'(z,alpha)'}
                summands = [gnd.link.link('summand', link=r) for r in reactionSuite.reactions
                        if int(r.attributes['ENDF_MT']) in info.sumCrossSections[MT]['MTs']]
                reaction = gnd.summedReaction.summedReaction( channelId[MT], "%s" % iChannel, ENDF_MT = MT,
                        crossSection = totalReaction.crossSection, Q = totalReaction.outputChannel.Q, summands = summands)
                reaction.setAttribute( 'date', Date )
                summedReactions.append( reaction )
    if( not( MT5Reaction is None ) ) :
        MT5Reaction.setLabel( "%s" % iChannel )
        reactionSuite.addReaction( MT5Reaction, iChannel )
        iChannel += 1
    for summedReac in summedReactions:  # [total], [nonelastic]
        summedReac.setLabel( "%s" % iChannel )
        reactionSuite.addSummedReaction( summedReac, iChannel )
        iChannel += 1
    for fissionComponent in fissionComponents:  # 1st-chance, 2nd-chance, etc
        fissionComponent.setLabel( "%s" % iChannel )
        fissionComponent.__class__ = gnd.fissionComponent.fissionComponent
        fissionComponent.moniker = gnd.reactions.base.fissionComponentToken
        reactionSuite.addFissionComponent( fissionComponent, iChannel )
        iChannel += 1
    for production in productions:
        production.setLabel( "%s" % iChannel )
        production.setAttribute( 'date', Date )
        reactionSuite.addProductionReaction( production, iChannel )
        iChannel += 1

    warningList = []
    # parse resonance section:
    try:
        if 151 not in MTDatas: mf2 = None   # no resonance data available
        else: mf2 = MTDatas.get(151).get(2)
        if( mf2 and ( singleMTOnly is None ) ) :
            resonances, resonanceMTs = readMF2(info, mf2, warningList)
            resonances.reconstructCrossSection = (LRP==1)   # LRP was read in from first line of ENDF file
            if resonances.unresolved and not resonances.resolved:
                if resonances.unresolved.tabulatedWidths.forSelfShieldingOnly:
                    resonances.reconstructCrossSection = False
            reactionSuite.addResonances( resonances )

            # add spins appearing in resonance region to the particle list
            for particle, spinParity in info.particleSpins.items():
                if particle=='target': particle = reactionSuite.target
                else: particle = reactionSuite.getParticle( particle )
                if isinstance(particle, gnd.xParticle.xParticle) and particle.name not in ('gamma','n','H1'):
                    # spin should be associated with ground state level:
                    particle.levels[0] = gnd.xParticle.nuclearLevel(name=particle.name+'_e0',
                        energy=physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty('0 eV'), label=0)
                    particle = particle.levels[0]
                particle.attributes['spin'] = spinParity[0]
                if spinParity[1]: particle.attributes['parity'] = spinParity[1]

            if resonances.reconstructCrossSection:
                # modify cross sections for relevant channels to indicate resonance contribution is needed:
                resonanceLink = gnd.link.link( label = "resonanceRegion", link = resonances, path = reactionSuite.resonances.toXLink( ) )

                for MT in resonanceMTs :
                    MTChannels = [ a for a in reactionSuite.reactions + reactionSuite.summedReactions
                            + reactionSuite.fissionComponents if( a.getENDL_CS_ENDF_MT()['MT'] == MT ) ]
                    if( len( MTChannels ) == 0 ) :
                        warningList.append( '       WARNING: unable to find channel corresponding to resonance data for MT%i' % MT )
                    elif( len( MTChannels ) == 1 ) :
                        crossSectionComponent = MTChannels[0].crossSection
                        backgroundForm = crossSectionComponent.forms[ crossSectionComponent.nativeData ]
                        crossSectionComponent.removeForm( backgroundForm, force = True )
                        crossSectionComponent.addFormAsNativeData( gnd.reactionData.crossSection.resonancesWithBackground( backgroundForm, resonanceLink ) )
                    else :
                        crossSectionComponent = MTChannels[0].crossSection
                        backgroundComponent = crossSectionComponent.forms[ crossSectionComponent.nativeData ].crossSection
                        backgroundForm = backgroundComponent.forms[ backgroundComponent.nativeData ]
                        backgroundComponent.removeForm( backgroundForm, force = True )
                        referredXSecForm = gnd.reactionData.crossSection.resonancesWithBackground( backgroundForm, resonanceLink )
                        backgroundComponent.addFormAsNativeData( referredXSecForm )
    except BadResonances as e:
        warningList.append( '       ERROR: unable to parse resonances! Error message: %s' % e )
        info.doRaise.append( warningList[-1] )

    try:
        """ parse covariances. This also requires setting up links from data to covariances, so we
        must ensure the links are synchronized """

        MTdict = {}
        info.MAT = MAT
        for ch in reactionSuite._getBaseAndDerivedReactions():
            MT = int(ch.attributes['ENDF_MT'])
            if MT in MTdict: MTdict[MT].append( ch )
            else: MTdict[MT] = [ch]
        covarianceSuite, linkData = parseCovariances( info, MTDatas, MTdict, singleMTOnly=singleMTOnly,
                resonances=getattr(reactionSuite,'resonances',None) )
        covarianceSuite.target = str(info.target)
        covarianceSuite.projectile = str(info.projectile)
        covarianceSuite.styles = {evaluatedStyle.name: evaluatedStyle}
        #covarianceSuite.removeExtraZeros() # disable for easier comparison to ENDF
    except Exception as e:
        warningList.append("       WARNING: Couldn't parse covariances! Error message: %s" % e )
        info.doRaise.append( warningList[-1] )
        covarianceSuite = gnd.covariances.covarianceSuite()        # create empty section
    
    for ZA in info.ZA_AWRMasses :
        mostCount = 0
        for AWR in info.ZA_AWRMasses[ZA] :
            if( info.ZA_AWRMasses[ZA][AWR] > mostCount ) :
                mostCount = info.ZA_AWRMasses[ZA][AWR]
                mostAWR = AWR
        info.ZA_AWRMasses[ZA] = mostAWR * info.masses.getMassFromZA( 1 )

    if level>0: # AWR is for isomer mass. Adjust info.ZAMasses to GS mass:
        info.ZA_AWRMasses[targetZA] -= physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty(level,'eV/c**2').getValueAs('amu')

    for name in reactionSuite.particles :                        # Fix up any particle whose mass is not defined (i.e., is None).
        if( name == 'gamma' ) : continue
        ZA = fudgeZA.gndNameToZ_A_Suffix( name )[3]
        particle = reactionSuite.particles[name]
        mass = particle.getMass( 'amu' )
        if( ZA in info.ZA_AWRMasses ) :
            mass = info.ZA_AWRMasses[ZA]
        elif( ( mass is None ) or ( ZA in info.ZAMasses ) ) :
            if( info.ZAMasses[ZA] is None ) :
                mass = info.masses.getMassFromZA( ZA )
            else :
                mass = abs( info.ZAMasses[ZA] )
        particle.setMass( physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( mass, 'amu' ) )
        # also add ground state level if there are discrete excited levels:
        if 1 in reactionSuite.particles[name].levels and 0 not in reactionSuite.particles[name].levels:
            reactionSuite.particles[name].levels[0] = gnd.xParticle.nuclearLevel(name=name+'_e0',
                    energy=physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty('0 eV'), label=0)

    MF12BaseMTsAndRange = [ [ 50, 92 ], [ 600, 650 ], [ 650, 700 ], [ 700, 750 ], [ 750, 800 ], [ 800, 850 ] ]

    if( singleMTOnly is None ) :
        branchingData = None
        #if( len( info.MF12_LO2 ) > 0 ) : reactionSuite.gammaBranching = {}
        for MTLO2, MF12_LO2 in sorted(info.MF12_LO2.items()) :  # The logic below assumes MTs are in ascending order per base MT.
            branchingBaseMT = None
            for MTBase, MTEnd in MF12BaseMTsAndRange :             # Determine base MT for this MTLO2
                if( MTBase < MTLO2 < MTEnd ) :
                    branchingBaseMT = MTBase
                    break
            if( not( branchingBaseMT is None ) ) :
                residualZA = endf_endl.ENDF_MTZAEquation( projectileZA, targetZA, branchingBaseMT )[0][-1]
                residual = getTypeNameENDF( info, residualZA, None )
                residualName = residual.getName( )
                level = MTLO2 - branchingBaseMT
                levelName, levelEnergy = '_e%d' % level, MF12_LO2[0]['ES']
                fullName = residualName + levelName
                levelEnergy_eV = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( levelEnergy, 'eV' )
                # compare this value to level energy from the particle list (from MF3 Q-value):
                particleLevelEnergy_eV = reactionSuite.getParticle(fullName).energy
                if levelEnergy_eV.value != particleLevelEnergy_eV.value:
                    if abs( levelEnergy_eV - particleLevelEnergy_eV ) / particleLevelEnergy_eV < 1e-4:
                        # MF12/MF14 value wins over MF3 value for small disagreements:
                        reactionSuite.getParticle(fullName).energy = levelEnergy_eV
                    else:
                        raise Exception ("MF12 parent level energy (%s eV) doesn't match known level" % particleLevelEnergy_eV )
                for MF12 in MF12_LO2 :
                    finalLevelEnergy = MF12['ESk']
                    if( finalLevelEnergy > 0. ) :
                        # find particle in the particleList with energy == finalLevelEnergy
                        finalParticles = [lev for lev in reactionSuite.getParticle( residualName ).levels.values()
                                if lev.energy.getValueAs('eV')==finalLevelEnergy]
                        if len(finalParticles)==1: finalParticle = finalParticles[0]
                        else:   # no exact match, look for levels within .01% of the exact value:
                            idx = 0
                            while True:
                                idx += 1
                                finalParticleName = residualName+'_e%i'%idx
                                if not reactionSuite.hasParticle(finalParticleName):
                                    raise Exception ("MF12 final level energy (%s eV) doesn't match known level" % finalLevelEnergy )
                                thisLevelEnergy = reactionSuite.getParticle(finalParticleName).energy.getValueAs('eV')
                                if abs( thisLevelEnergy - finalLevelEnergy ) < 1e-4 * finalLevelEnergy:
                                    finalParticle = reactionSuite.getParticle(finalParticleName)
                                    break   # found it
                    else: finalParticle = reactionSuite.getParticle(residualName+'_e0')
                    gammaTransition = 1.
                    if( len( MF12['branching'] ) > 2 ) : gammaTransition = MF12['branching'][1]
                    gamma = gnd.xParticle.nuclearLevelGamma( finalParticle, MF12['angular'], MF12['branching'][0], 1-gammaTransition )
                    reactionSuite.getParticle(fullName).addGamma( gamma )
            else :
                raise Exception( "Could not determine base MT for MF=12's MT=%s" % MTLO2 )
    sys.stdout.flush( )
    for warning in warningList : logs.write( warning + '\n', stderrWriting = True )

    if( len( info.doRaise ) > 0 and not skipBadData ) :
        info.logs.write( '\nRaising due to following errors:\n' )
        for err in info.doRaise : info.logs.write( err + '\n' )
        raise Exception( 'len( info.doRaise ) > 0' )

    return( reactionSuite, covarianceSuite )

if( __name__ == '__main__' ) :

    flags = processingInfo.tempInfo( )
    flags['verbosity'] = 31
    x, c = endfFileToGND( sys.argv[1] )
    f = open( 'test.xml', 'w' )
    f.write( '\n'.join( x.toXMLList( flags ) + [ '' ] ) )
    f.close( )
    if c: # covariances
        f = open( 'test-covar.xml', 'w' )
        f.write( '\n'.join( c.toXMLList( flags ) + [ '' ] ) )
        f.close()
