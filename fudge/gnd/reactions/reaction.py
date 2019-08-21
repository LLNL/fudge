# <<BEGIN-copyright>>
# <<END-copyright>>

"""
This module contains the reaction class.
"""

QNotApplicableToken = 'notApplicable'

from fudge.core.utilities import fudgeExceptions
from pqu import physicalQuantityWithUncertainty
from fudge.legacy.converting import endfFormats, endf_endl

from ..productData import distributions
from fudge.processing import processingInfo
import fudge
from . import base

__metaclass__ = type

class reaction( base.base_reaction ) :
    """This is the class for a normal gnd reaction"""

    def __init__( self, outputChannel, label, ENDF_MT, crossSection = None, URRProbabilityTable = None, documentation = None, attributes = {} ) :
        """Creates a new reaction object. Reaction is two-body or uncorrelated-body, depending on
        the outputChannel type. This class is only meant to be used for 'distinct' reactions (distinct reactions
        do not overlap any other reactions; altogether they sum up to total).
        To store a sum over these distinct reactions, use the summedReaction class instead. """

        base.base_reaction.__init__( self, base.reactionToken, label, ENDF_MT, crossSection, documentation, attributes )

        if( type( outputChannel ) == type( '' ) ) : outputChannel = fudge.gnd.channels.channel( 'output', outputChannel )
        if( not isinstance( outputChannel, fudge.gnd.channels.channel ) ) :
            raise fudgeExceptions.FUDGE_Exception( 'Input channel not instance of class channel.' )
        outputChannel.setParent( self )
        self.URRProbabilityTable = URRProbabilityTable
        self.data = {}
        self.outputChannel = outputChannel

    def __eq__(self, other):
        if (not base.isGNDReaction( other )): return False
        return( ( self.parent.projectile == other.parent.projectile ) and ( self.parent.target == other.parent.target ) 
            and ( self.outputChannel == other.outputChannel ) )
    
    def __cmp__( self, other ) :
        """Test if self is <, == or > other."""

        if( not base.isGNDReaction( other ) ) : raise fudgeExceptions.FUDGE_Exception( "Other not an reaction object." )
        if( self.parent.projectile < other.parent.projectile ) : return( -1 )
        if( self.parent.projectile > other.parent.projectile ) : return(  1 )
        if( self.parent.target < other.parent.target ) : return( -1 )
        if( self.parent.target > other.parent.target ) : return(  1 )
        if( self.outputChannel < other.outputChannel ) : return( -1 )
        if( self.outputChannel > other.outputChannel ) : return(  1 )
        return( 0 )

    def __str__( self ) :

        return( str(self.outputChannel) )

    def addData( self, data ) :

        if( data.genre not in self.data ) :
            self.data[data.genre] = data
        else :
            raise Exception( 'self.data has genre %s' % data.genre )

    def getQ( self, unit, final = True, groundStateQ = False ) :        # ????? unit for Q is needed, badly.
        """Returns the Q-value for this input/output channel. It will be converted to a float if possible, otherwise a string value is returned."""

        def getLevel( particle, unit ) :

            return( particle.getLevelAsFloat( unit = unit, default = 0. ) )

        def getFinalProductsLevelMultiplicity( particle, final ) :
            """Returns a list of [ productName, multiplicity ] for the final products of particle."""

            if( ( final == False ) or ( particle.decayChannel == None ) ) :
                if( particle.particle.getToken( ) == 'gamma' ) :
                    pms = [ [ particle, 0., 0. ] ]
                else :
                    pms = [ [ particle, particle.multiplicity.getConstant( ), getLevel( particle, unit ) ] ]
            else :
                pms = []
                for product in particle.decayChannel :
                    pms += getFinalProductsLevelMultiplicity( product, final )
            return( pms )

        if( self.getAttribute( 'genre' ) == 'production' ) :
            Q = QNotApplicableToken
        else :
            return( self.outputChannel.getConstantQAs( unit, final = final ) )
            Q = self.getAttribute( 'Q' )
            if( not( Q is None ) ) :
                if( isinstance( Q, physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty ) ) :
                    Q = Q.getValueAs( unit )
                else :
                    pass                                        # ????? This needs work. For example, check for energy dependent Q.
            else :
                pm = {}
                for particle in [ self.parent.projectile, self.parent.target ] :
                    if( particle.getToken( ) not in pm ) : pm[particle.getToken( )] = [ particle, 0, 0. ]
                    pm[particle.getToken( )][1] += 1
                    pm[particle.getToken( )][2] += getLevel( particle, unit )
                for particle in self.outputChannel :
                    pms = getFinalProductsLevelMultiplicity( particle, final = final )
                    for p, m, level in pms :
                        if( p.particle.getToken( ) not in pm ) : pm[p.particle.getToken( )] = [ p, 0, 0. ]
                        pm[p.particle.getToken( )][1] -= m
                        pm[p.particle.getToken( )][2] -= level
                Q = 0.
                levelQ = 0.
                try :
                    for p in pm :
                        Q += pm[p][1] * pm[p][0].getMass( 'MeV/c**2' )  # Unit 'MeV' should not be hard-wired.
                        levelQ += pm[p][2]
                except :
                    print p, pm[p]
                    raise
                if( groundStateQ ) : levelQ = 0.
                Q = Q + levelQ
        return( Q )

    def getThreshold( self, unit ) :

        Q = self.getQ( unit = unit, final = False )
        if( Q == QNotApplicableToken ) : return( Q )
        if( Q >= 0. ) : return( 0. )
        return( -Q * ( 1. + self.parent.projectile.getMass( 'amu' ) / self.parent.target.getMass( 'amu' ) ) )

    def getURRProbabilityTable( self ) :
        """Returns the channel's URR probability table present, else None. This is a kludge"""

        return( self.URRProbabilityTable )

    def heatCrossSection( self, temperature, EMin, lowerlimit = None, upperlimit = None, interpolationAccuracy = 0.002, heatAllPoints = False, 
        doNotThin = True, heatBelowThreshold = True, heatAllEDomain = True ) :

        return( self.crossSection.heat( temperature, EMin, lowerlimit, upperlimit, interpolationAccuracy, heatAllPoints, doNotThin, 
            heatBelowThreshold, heatAllEDomain ) )

    def calculateDepositionData( self, processInfo, verbosityIndent ) :

        if( isinstance( self.outputChannel, fudge.gnd.channels.sumOfRemainingOutputChannels ) ) : return
        tempInfo = processingInfo.tempInfo( )
        tempInfo['reactionSuite'] = self.getParent( )
        tempInfo['outputChannel'] = self
        tempInfo['EMin'], tempInfo['EMax'] = self.getDomain( )
        tempInfo['incidentEnergyUnit'] = self.crossSection.getIncidentEnergyUnit( )
        if( processInfo['verbosity'] >= 10 ) : print '%s%s' % ( verbosityIndent, self.outputChannel.toString( simpleString = True ) )
        for productIndex, product in enumerate( self.outputChannel ) : 
            tempInfo['productIndex'] = str( productIndex )
            tempInfo['productToken'] = product.getToken( )
            tempInfo['productLabel'] = product.getLabel( )
            product.calculateDepositionData( processInfo, tempInfo, verbosityIndent + '    ' )

    def check( self, info ) :

        warnings = self.__checkCrossSection__( info )
        if not info['crossSectionOnly']:
            warnings += self.__checkOutputChannel__( info )
        return warnings

    def isBasicReaction( self ) :

        return( True )

    def isCompleteReaction( self ) :

        return( True )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        from fudge.gnd import xParticle

#        if( isinstance( self.outputChannel, fudge.gnd.channels.sumOfRemainingOutputChannels ) ) : return
        crossSection = self.crossSection
        tempInfo['incidentEnergyUnit'] = self.crossSection.getIncidentEnergyUnit( )
        if( 'LLNL_Pn' in processInfo['styles'] ) :
            if( processInfo['particles'][tempInfo['reactionSuite'].projectile.getName( )].groups[-1] <= crossSection.domainMin( ) ) : return

        tempInfo['outputChannel'] = self
        if( processInfo['verbosity'] >= 10 ) : print '%s%s' % ( verbosityIndent, self.outputChannel.toString( simpleString = True ) )
        projectile, target = processInfo.getProjectileName( ), processInfo.getTargetName( )
        if( 'LLNL_Pn' in processInfo['styles'] ) :
            fluxl0 = processInfo['flux']['data'][0]
            tempInfo['groupedFlux'] = fluxl0.groupOneFunction( processInfo.getParticleGroups( projectile ) )
        tempInfo['crossSection'] = crossSection.process( processInfo, tempInfo, verbosityIndent + '    ' )

        if( 'LLNL_Pn' in processInfo['styles'] ) :
            try :
                availableEnergy = self.getData( fudge.gnd.reactionData.base.availableEnergyToken )
            except :
                availableEnergy = fudge.gnd.reactionData.availableEnergy.component( )
                self.addData( availableEnergy )
            availableEnergy.makeGrouped( processInfo, tempInfo )

            try :
                availableMomentum = self.getData( fudge.gnd.reactionData.base.availableMomentumToken )
            except :
                availableMomentum = fudge.gnd.reactionData.availableMomentum.component( )
                self.addData( availableMomentum )
            availableMomentum.makeGrouped( processInfo, tempInfo )

            if( self.outputChannel.Q.nativeData in [ fudge.gnd.tokens.pointwiseFormToken ] ) :
                self.getData( fudge.gnd.channelData.energyDependentQToken ).makeGrouped( processInfo, tempInfo )    # ????? This is not right.

            tempInfo['transferMatrixComment'] = tempInfo['reactionSuite'].inputParticlesToReactionString( suffix = " --> " ) +  \
                self.outputChannel.toString( simpleString = True )
        for productIndex, product in enumerate( self.outputChannel ) : 
            tempInfo['productIndex'] = str( productIndex )
            tempInfo['productToken'] = product.getToken( )
            tempInfo['productGroupToken'] = tempInfo['productToken']
            if( isinstance( product.particle, xParticle.nuclearLevel ) ) : tempInfo['productGroupToken'] = product.particle.groundState.getToken( )
            tempInfo['productLabel'] = product.getLabel( )
            product.process( processInfo, tempInfo, verbosityIndent + '    ' )

    def toLLNLSnCOut( self, toOtherData, moreOtherData, verbosityIndent ) :

        from fudge.processing.deterministic import gndToLLNLSnCOut
        if( toOtherData['verbosity'] >= 10 ) : print '%s%s' % ( verbosityIndent, self.outputChannel.toString( simpleString = True ) )
        if( moreOtherData.temperature is None ) : moreOtherData.temperature = self.attributes['temperature']
        crossSection = gndToLLNLSnCOut.LLNLSnCOutCrossSection( self )
        isProductionChannel = self.getAttribute( 'genre' ) == 'production'
        if( not isProductionChannel ) :
            moreOtherData.CNumbers.append( self.getENDL_CS_ENDF_MT( )['C'] )
            moreOtherData.coutData[5] = moreOtherData.coutData[5] + \
                self.data[fudge.gnd.reactionData.base.availableEnergyToken].getFormByToken( fudge.gnd.tokens.groupedWithCrossSectionFormToken ).data
            if( fudge.gnd.channelData.energyDependentQToken in self.data ) :
                moreOtherData.coutData[6] = moreOtherData.coutData[6] + \
                    self.getData( fudge.gnd.channelData.energyDependentQToken ).getFormByToken( fudge.gnd.tokens.groupedWithCrossSectionFormToken ).data
            else :
                Q = self.getQ( 'MeV' )
                moreOtherData.coutData[6] = moreOtherData.coutData[6] + Q * self.crossSection.getFormByToken( fudge.gnd.tokens.groupedFormToken ).data
        addToFissionMatrix = self.outputChannel.isFission( )
        reactionWithProduct = { }
        for product in moreOtherData.reactionWithProduct : reactionWithProduct[product] = 0
        for product in self.outputChannel : 
            product.toLLNLSnCOut( toOtherData, moreOtherData, crossSection, reactionWithProduct, addToFissionMatrix, verbosityIndent + '    ' )
        projectile = toOtherData.getProjectileName( )
        if( not isProductionChannel ) :
            moreOtherData.coutData[4].append( crossSection )
            if( projectile in reactionWithProduct ) : reactionWithProduct[projectile] = 1
        for product in moreOtherData.reactionWithProduct : moreOtherData.reactionWithProduct[product] += reactionWithProduct[product]

    def setURRProbabilityTable( self, URRProbabilityTable ) :
        """Sets channel's URR probability table to URRProbabilityTable. This is a kludge"""

        self.URRProbabilityTable = URRProbabilityTable

    def toENDF6( self, endfMFList, flags, targetInfo, verbosityIndent = '' ) :
        from fudge.legacy.converting import gndToENDF6

        def addDecayProducts( parent, products ) :

            if( parent.decayChannel is None ) : return( False )
            doMF4AsMF6 = False
            for product in parent.decayChannel :
                if( product.distributions.components and product.getName() != 'gamma' ) :
                    doMF4AsMF6 = True
                    products.append( product )
            return( doMF4AsMF6 )

        def divideIgnoring0DividedBy0( self, other ) :

            from fudge.core.math.xData import XYs
            d = self.union( other.xSlice( xMin = self.domainMin( ), xMax = self.domainMax( ) ) )
            result = []
            for p in d :
                vp = other.getValue( p[0] )
                if( vp == 0 ) :
                    if( p[1] != 0 ) : raise Exception( 'Divide non-zero number by zero at %e' % p[0] )
                else :
                    p[1] = p[1] / vp
                result.append( [ p[0], p[1] ] )
            return( XYs.pointwiseXY( data = result ) )

        def thinWeights( weights ) :

            i, n, thinnedWeights = 1, len( weights ), [ weights[0] ]
            if( n == 2 ) :
                thinnedWeights.append( weights[-1] )
            else :
                while( i < n ) :
                    y = thinnedWeights[-1][1]
                    while( True ) :
                        i += 1
                        if( i == n ) : break
                        if( abs( y - weights[i][1] ) > 1e-8 * y ) : break
                        if( abs( y - weights[i-1][1] ) > 1e-8 * y ) : break
                    thinnedWeights.append( weights[i-1] )
            return( thinnedWeights )

        outputChannel = self.outputChannel
        if( outputChannel.genre == fudge.gnd.channels.productionGenre ) :
            print '        toENDF6 does not support writing of "%s" channel' % outputChannel.genre
            return
        MT = int( self.attributes['ENDF_MT'] )
        if( flags['verbosity'] >= 10 ) : print '%s%s' % ( verbosityIndent, outputChannel.toString( simpleString = True ) )
        targetInfo['Q'] = self.getQ( 'eV', groundStateQ = True )
        targetInfo['QM'] = None

        LR, level, tryLR = 0, 0., False
        if( outputChannel.getGenre( ) == fudge.gnd.channels.twoBodyGenre ) :
            tryLR = True
        elif( outputChannel.getGenre( ) == fudge.gnd.channels.NBodyGenre ) :
            if( MT == 91 ) :
                tryLR = True
        if( tryLR ) :
            secondProduct = outputChannel[1]
            primaryResidualName, decayProducts, decayChannel, = secondProduct.getName( ).split( '_' )[0], [], secondProduct.decayChannel
            if( not( decayChannel is None ) ) :
                for decayProduct in decayChannel :
                    decayProductName = decayProduct.getName( )
                    if( decayProductName not in [ primaryResidualName, 'gamma' ] ) : decayProducts.append( decayProductName )
            if( len( decayProducts ) == 1 ) :   # Kludge for Carbon breakup into 3 alphas.
                if( ( primaryResidualName == 'C' ) and ( decayProducts == [ 'He4' ] ) ) : LR = 23
            elif( len( decayProducts ) > 1 ) :                                        # This must be a breakup reaction.
                MTProducts = endf_endl.endfMTtoC_ProductList( 0, '' )
                MTProducts.productCounts[outputChannel[0].getName( )] += 1
                for decayProduct in decayProducts[:-1] : MTProducts.productCounts[decayProduct] += 1
                for MT_LR in [ 22, 23, 24, 25, 28, 29, 30, 32, 33, 34, 35, 36 ] :   # 39 and 40 not allowed in ENDF6
                    if( endf_endl.endfMTtoC_ProductLists[MT_LR].productCounts == MTProducts.productCounts ) :
                        LR = MT_LR
                        break
                if( ( LR == 32 ) and ( primaryResidualName == 'B10' ) and ( decayProducts[-1] == 'He4' ) ) : LR = 35   # Kludge for bad data.
            if( LR != 0 ) : targetInfo['QM'] = decayChannel.Q.forms['constant'].getValue( 0, 'eV' )
        targetInfo['LRs'][MT] = LR

        for product in outputChannel :
            tmp = product.getLevelAsFloat( 'eV' )
            if tmp: level = tmp

        crossSection = self.getCrossSection( )
        targetInfo['EMin'], targetInfo['EMax'] = crossSection.getDomain( )
        crossSection.toENDF6( MT, endfMFList, targetInfo, level, LR )

        products = []
        doMF4AsMF6 = False
        for product in outputChannel :
            if( ( outputChannel.getGenre( ) == fudge.gnd.channels.twoBodyGenre ) and
                    ( distributions.base.angularComponentToken in product.distributions.components ) ) :
                nativeData = product.getDistributionNativeData( )
                if( nativeData != distributions.base.noneComponentToken ) :
                    component = product.getDistributionComponentByToken( nativeData )
                    if( component.nativeData in [ distributions.base.recoilFormToken ] ) : doMF4AsMF6 = True
            if( product.distributions.components ) : products.append( product )
            doMF4AsMF6 = doMF4AsMF6 or addDecayProducts( product, products )
        if( outputChannel.isFission( ) ) : 
            doMF4AsMF6 = False
            if( not( outputChannel.fissionEnergyReleased is None ) ) :
                outputChannel.fissionEnergyReleased.toENDF6( 1, endfMFList, flags, targetInfo )
            delayedNubar = None                                 # Special work to get delayed nubar weights.
            for product in outputChannel :
                if( 'emissionMode' in product.attributes ) :
                    if( product.getAttribute( 'emissionMode' ) == 'delayed' ) :
                        if( delayedNubar is None ) :
                            delayedNubar = product.multiplicity.getFormByToken( fudge.gnd.tokens.pointwiseFormToken )
                        else :
                            delayedNubar = delayedNubar + product.multiplicity.getFormByToken( fudge.gnd.tokens.pointwiseFormToken )
                        targetInfo['delayedRates'].append( product.getAttribute( 'decayRate' ).getValueAs('1/s') )
            if( delayedNubar is not None ) :
                targetInfo['totalDelayedNubar'] = delayedNubar
                for product in outputChannel :
                    if( 'emissionMode' in product.attributes ) :
                        if( product.getAttribute( 'emissionMode' ) == 'delayed' ) :
                            weight = divideIgnoring0DividedBy0( product.multiplicity.getFormByToken( fudge.gnd.tokens.pointwiseFormToken ), delayedNubar )
                            product.ENDF6_delayedNubarWeights = thinWeights( weight )
        targetInfo['doMF4AsMF6'] = doMF4AsMF6
        targetInfo['MF6LCTs'], targetInfo['gammas'] = [], []
        for productIndex, product in enumerate( outputChannel ) :
            if( ( product.getName( ) == 'gamma' ) and not( outputChannel.getGenre( ) == fudge.gnd.channels.twoBodyGenre ) ) :
                targetInfo['gammas'].append( product )
                continue
            targetInfo['productIndex'] = str( productIndex )
            targetInfo['productToken'] = product.getToken( )
            targetInfo['productLabel'] = product.getLabel( )
            product.toENDF6( MT, endfMFList, flags, targetInfo, verbosityIndent = verbosityIndent + '    ' )
        gammas = targetInfo['gammas']
        if( len( gammas ) ) :
            gamma = gammas[0]
            targetInfo['zapID'] = gamma.getName( )
            targetInfo['particleMass'] = gamma.getMass( 'eV/c**2' )
            multiplicity = gamma.multiplicity
            if( gamma.getAttribute( 'ENDFconversionFlag' ) == 'MF6' ) :
                targetInfo['multiplicity'] = gamma.multiplicity
                gndToENDF6.gammasToENDF6_MF6( MT, endfMFList, flags, targetInfo, gammas )
            elif( gamma.getAttribute( 'ENDFconversionFlag' ) == 'MF13' ) :
                try :
                    targetInfo['crossSection'] = self.crossSection.forms[fudge.gnd.tokens.pointwiseFormToken]
                except :
                    targetInfo['crossSection'] = self.crossSection.forms[fudge.gnd.tokens.linearFormToken]
                gndToENDF6.gammasToENDF6_MF12_13( MT, 13, endfMFList, flags, targetInfo, gammas )
            elif( multiplicity.nativeData == fudge.gnd.tokens.constantFormToken ) :      # This should probably be changed to unknown in some cases?????
                pass
            elif( multiplicity.nativeData == fudge.gnd.tokens.unknownFormToken ) :       # This should probably be changed to unknown in some cases?????
                pass
            else :
                gndToENDF6.gammasToENDF6_MF12_13( MT, 12, endfMFList, flags, targetInfo, gammas )
        if( len( targetInfo['MF6LCTs'] ) > 0 ) :
            LCT = targetInfo['MF6LCTs'][0]
            for i in targetInfo['MF6LCTs'] :
                if( not( LCT is None ) ) : break
            if( LCT is None ) : LCT = 2
            for i in targetInfo['MF6LCTs'] :
                if( i is None ) : continue
                if( i != LCT ) : LCT = 3
            endfMFList[6][MT].insert( 0, endfFormats.endfHeadLine( targetInfo['ZA'], targetInfo['mass'], 0, LCT, len( targetInfo['MF6LCTs'] ), 0 ) )
            endfMFList[6][MT].append( endfFormats.endfSENDLineNumber( ) )

    def toString( self, indent = '' ) :

        s, p = indent, ''
        for particle in self.outputChannel :
            s += p + particle.toString( )
            p = ' + '
        return( s + '\n' )

def parseXMLNode( reactionElement, linkData={} ):
    """Translate a <reaction> element from xml into a reaction class instance."""
    try:
        crossSection = fudge.gnd.reactionData.crossSection.parseXMLNode( reactionElement.find('crossSection'), linkData )
        outputChannel = fudge.gnd.channels.parseXMLNode( reactionElement.find('outputChannel'), linkData )
    except Exception as e:
        raise Exception, '/reactionSuite/reaction[@label="%s"]%s' % (reactionElement.get('label'),e)
    outputChannel.fissionGenre = reactionElement.get("fissionGenre")
    MT = int( reactionElement.get('ENDF_MT') )
    attributes = {'date': reactionElement.get('date')}
    reac = reaction( outputChannel, reactionElement.get('label'), MT, crossSection, attributes = attributes )

    if reactionElement.find('documentations'):
        for doc in reactionElement.find('documentations'):
            reac.addDocumentation( fudge.gnd.documentation.parseXMLNode(doc) )

    # extra data may be present in derived files:
    if reactionElement.find('availableEnergy') is not None:
        reac.data['availableEnergy'] = fudge.gnd.reactionData.availableEnergy.parseXMLNode(
                reactionElement.find('availableEnergy') )
    if reactionElement.find('availableMomentum') is not None:
        reac.data['availableMomentum'] = fudge.gnd.reactionData.availableMomentum.parseXMLNode(
                reactionElement.find('availableMomentum') )
    return reac
