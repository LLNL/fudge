# <<BEGIN-copyright>>
# <<END-copyright>>

import os, subprocess
from fudge.core.utilities import fudgeFileMisc, subprocessing, times
from pqu.physicalQuantityWithUncertainty import PhysicalQuantityWithUncertainty

from fudge.core.math.xData import axes, XYs, LegendreSeries
GND2ProcessingInterpolations = { axes.unitBaseToken : 'unit base', axes.flatToken : 'flat', 
    axes.linearToken + axes.linearToken : 'linear-linear', axes.linearToken + axes.logToken : 'linear-log', 
    axes.logToken    + axes.linearToken : 'log-linear',    axes.logToken    + axes.logToken : 'log-log' }

versionStr = "xndfgenTransferMatrix: version 1.0"
doubleFmt = "%19.12e"
lowerEps = 1e-8
upperEps = 1e-8

if( os.path.exists( 'get_transfer' ) ) :
    srcPath = os.path.abspath( './' )
else :
    srcPath = os.path.split( os.path.abspath( __file__ ) )[0]

def twoBodyTransferMatrix( processInfo, projectile, product, crossSection, angularData, mProj, mTarget, mProd, mRes, Q, comment = None ) :

    from fudge.gnd.productData import distributions

    if( isinstance( angularData, distributions.angular.mixedRanges ) ) :    # Legendre and tabulated data can still be piecewise.
        crossSectionLower = crossSection.xSlice( xMax = angularData.LegendreForm.domainMax() )
        crossSectionUpper = crossSection.xSlice( xMin = angularData.tabulatedForm.domainMin() )
        L_TM1, L_TME = twoBodyTransferMatrix( processInfo, projectile, product, crossSectionLower, angularData.LegendreForm, mProj, mTarget, mProd, mRes, 
            Q, comment = comment )
        t_TM1, t_TME = twoBodyTransferMatrix( processInfo, projectile, product, crossSectionUpper, angularData.tabulatedForm, mProj, mTarget, mProd, mRes, 
            Q, comment = comment )
        TM1 = addTMs( [ L_TM1, t_TM1 ] )
        TME = addTMs( [ L_TME, t_TME ] )
    elif( isinstance( angularData, distributions.angular.piecewise ) ) :
        raise Exception( 'distributions.angular.piecewise not supported' )
    elif( isinstance( angularData, distributions.angular.LegendrePiecewise ) ) :
        TM1s, TMEs = [], []
        for region in angularData :
            TM1, TME = twoBodyTransferMatrix2( processInfo, projectile, product, crossSection, region, mProj, mTarget, mProd, mRes, Q, comment = comment )
            TM1s.append( TM1 )
            TMEs.append( TME )
        TM1 = addTMs( TM1s )
        TME = addTMs( TMEs )
    else :
        TM1, TME = twoBodyTransferMatrix2( processInfo, projectile, product, crossSection, angularData, mProj, mTarget, mProd, mRes, Q, comment = comment )
    return( TM1, TME )
        
def twoBodyTransferMatrix2( processInfo, projectile, product, crossSection, angularData, mProj, mTarget, mProd, mRes, Q, comment = None ) :

    from fudge.gnd.productData import distributions

    logFile, lMax = processInfo['logFile'], processInfo.getParticle_lMax( product )
    projectileGroupBoundaries, productGroupBoundaries = processInfo.getParticleGroups( projectile ), processInfo.getParticleGroups( product )
    fluxes, workDir = processInfo['flux']['data'], processInfo['workDir']

    if( isinstance( angularData, distributions.angular.recoil ) ) : angularData = angularData.getNumericalDistribution( )

    s = versionStr + '\n'
    if( isinstance( angularData, LegendreSeries.W_XYs_LegendreSeries ) ) :
        s += "Process: 'Legendre two body transfer matrix'\n"
        s += "Quadrature method: adaptive\n"
    else :
        s += "Process: 'two body transfer matrix'\n"

    if( comment is not None ) : s += "Comment: %s\n" % comment
    s += "Projectile's mass: %s\n" % floatToString( mProj )
    s += "Target's mass: %s\n" % floatToString( mTarget )
    s += "Product's mass: %s\n" % floatToString( mProd )
    s += "Residual's mass: %s\n" % floatToString( mRes )

    s += "Reaction's Q value: %s\n" % floatToString( Q )

    s += commonDataToString( lMax, projectileGroupBoundaries, productGroupBoundaries, fluxes, crossSection )
    s += angularToString( angularData, crossSection )
    return( executeCommand( logFile, '%s/get_transfer' % srcPath, s, workDir ) )

def wholeAtomScattering( processInfo, projectile, product, formFactor, comment = None ) :

    logFile, lMax, projectileGroupBoundaries, productGroupBoundaries, fluxes, workDir = processInfo['logFile'], processInfo.getParticle_lMax( product ), \
        processInfo.getParticleGroups( projectile ), processInfo.getParticleGroups( product ), processInfo['flux']['data'], processInfo['workDir']

    s  = versionStr + '\n'
    s += "Process: 'whole atom scattering'\n"

    if( comment is not None ) : s += "Comment: %s\n" % comment
    s += commonDataToString( lMax, projectileGroupBoundaries, productGroupBoundaries, fluxes, None )
    s += '\n'.join( twoDToString( "FormFactorData", formFactor ) )
    return( executeCommand( logFile, '%s/get_transfer' % srcPath, s, workDir ) )

def comptonScattering( processInfo, projectile, product, formFactor, comment = None ) :

    logFile, lMax, projectileGroupBoundaries, productGroupBoundaries, fluxes, workDir = processInfo['logFile'], processInfo.getParticle_lMax( product ), \
        processInfo.getParticleGroups( projectile ), processInfo.getParticleGroups( product ), processInfo['flux']['data'], processInfo['workDir']

    s  = versionStr + '\n'
    s += "Process: 'Compton scattering'\n"

    if( comment is not None ) : s += "Comment: %s\n" % comment
    s += commonDataToString( lMax, projectileGroupBoundaries, productGroupBoundaries, fluxes, None )
    s += '\n'.join( twoDToString( "ScatteringFactorData", formFactor ) )
    return( executeCommand( logFile, '%s/get_transfer' % srcPath, s, workDir ) )

def ENDFEMuEpP_TransferMatrix( processInfo, projectile, product, crossSection, angularEnergyData, multiplicity, comment = None ) :
    """This is ENDF MF = 6, LAW = 7 type data."""

    logFile, lMax, projectileGroupBoundaries, productGroupBoundaries, fluxes, workDir = processInfo['logFile'], processInfo.getParticle_lMax( product ), \
        processInfo.getParticleGroups( projectile ), processInfo.getParticleGroups( product ), processInfo['flux']['data'], processInfo['workDir']

    s  = versionStr + '\n'
    s += "Process: ENDF Double differential EMuEpP data\n"

    if( comment is not None ) : s += "Comment: %s\n" % comment
    s += "Quadrature method: adaptive\n"
    s += commonDataToString( lMax, projectileGroupBoundaries, productGroupBoundaries, fluxes, crossSection, multiplicity = multiplicity )
    s += EMuEpPDataToString( angularEnergyData )
    return( executeCommand( logFile, '%s/get_transfer' % srcPath, s, workDir ) )

def EMuEpP_TransferMatrix( processInfo, projectile, product, crossSection, angularData, EMuEpPData, multiplicity, comment = None ) :
    """This is LLNL I = 1, 3 type data."""

    logFile, lMax, projectileGroupBoundaries, productGroupBoundaries, fluxes, workDir = processInfo['logFile'], processInfo.getParticle_lMax( product ), \
        processInfo.getParticleGroups( projectile ), processInfo.getParticleGroups( product ), processInfo['flux']['data'], processInfo['workDir']

    s  = versionStr + '\n'
    s += "Process: 'Double differential EMuEpP data transfer matrix'\n"

    if( comment is not None ) : s += "Comment: %s\n" % comment
    s += commonDataToString( lMax, projectileGroupBoundaries, productGroupBoundaries, fluxes, crossSection, multiplicity = multiplicity )
    s += angularToString( angularData, crossSection )
    s += EMuEpPDataToString( EMuEpPData )
    return( executeCommand( logFile, '%s/get_transfer' % srcPath, s, workDir ) )

def ELEpP_TransferMatrix( processInfo, projectile, product, crossSection, LEEpPData, multiplicity, comment = None ) :

    logFile, lMax, projectileGroupBoundaries, productGroupBoundaries, fluxes, workDir = processInfo['logFile'], processInfo.getParticle_lMax( product ), \
        processInfo.getParticleGroups( projectile ), processInfo.getParticleGroups( product ), processInfo['flux']['data'], processInfo['workDir']

    s  = versionStr + '\n'
    s += "Process: 'Legendre EEpP data transfer matrix'\n"
    if( comment is not None ) : s += "Comment: %s\n" % comment

    s += commonDataToString( lMax, projectileGroupBoundaries, productGroupBoundaries, fluxes, crossSection, multiplicity = multiplicity )
    s += LEEpPDataToString( LEEpPData )
    return( executeCommand( logFile, '%s/get_transfer' % srcPath, s, workDir ) )

def Legendre_TransferMatrix( processInfo, projectile, product, crossSection, LegendreData, multiplicity, comment = None ) :

    logFile, lMax, projectileGroupBoundaries, productGroupBoundaries, fluxes, workDir = processInfo['logFile'], processInfo.getParticle_lMax( product ), \
        processInfo.getParticleGroups( projectile ), processInfo.getParticleGroups( product ), processInfo['flux']['data'], processInfo['workDir']

    s  = versionStr + '\n'
    s += "Process: Legendre energy-angle data\n"
    if( comment is not None ) : s += "Comment: %s\n" % comment

    s += "Quadrature method: adaptive\n"
    s += commonDataToString( lMax, projectileGroupBoundaries, productGroupBoundaries, fluxes, crossSection, multiplicity = multiplicity )
    s += LegendreDataToString( LegendreData )

    return( executeCommand( logFile, '%s/get_transfer' % srcPath, s, workDir ) )

def uncorrelated_EMuP_EEpP_TransferMatrix( processInfo, projectile, product, crossSection, angularData, EEpPData, multiplicity, comment = None, weight = None ) :

    from fudge.gnd.productData import distributions

    logFile, lMax, projectileGroupBoundaries, productGroupBoundaries, fluxes, workDir = processInfo['logFile'], processInfo.getParticle_lMax( product ), \
        processInfo.getParticleGroups( projectile ), processInfo.getParticleGroups( product ), processInfo['flux']['data'], processInfo['workDir']

    sSpecific = ''
    sCommon = commonDataToString( lMax, projectileGroupBoundaries, productGroupBoundaries, fluxes, crossSection, multiplicity = multiplicity )
    if( isinstance( EEpPData, distributions.energy.semiPiecewise ) ) : EEpPData = EEpPData.toPointwiseLinear( 1e-6, lowerEps = 1e-12, upperEps = 1e-12 )
    if( isinstance( EEpPData, distributions.energy.functionalBase ) ) :
        sProcess, sSpecific, sData = energyFunctionToString( EEpPData, weight = weight )
    elif( isinstance( EEpPData, distributions.energy.weightedFunctionals ) ) :
        TM1s, TMEs = [], []
        for weight in EEpPData :
            TM1, TME = uncorrelated_EMuP_EEpP_TransferMatrix( processInfo, projectile, product, crossSection, angularData, weight.functional, multiplicity, 
                comment = comment, weight = weight )
        TM1 = addTMs( TM1s )
        TME = addTMs( TMEs )
        return( TM1, TME )
    else :
        sProcess = "Process: 'Uncorrelated energy-angle data transfer matrix'\n"
        sData = angularToString( angularData, crossSection )
        sData += EEpPDataToString( EEpPData )
    if( comment is not None ) : sProcess += "Comment: %s\n" % comment
    s  = versionStr + '\n' + sProcess + sSpecific + sCommon + sData
    return( executeCommand( logFile, '%s/get_transfer' % srcPath, s, workDir ) )

def discreteGammaAngularData( processInfo, projectile, product, gammaEnergy, crossSection, angularData, multiplicity = 1, comment = None ) :
    """Currently, only isotropic (i.e., l = 0) data are returned. That is, lMax and angularData are ignored. This routine is also used
    for pair-production which pass angularData as None."""

    logFile, lMax, projectileGroupBoundaries, productGroupBoundaries, fluxes, workDir = processInfo['logFile'], processInfo.getParticle_lMax( product ), \
        processInfo.getParticleGroups( projectile ), processInfo.getParticleGroups( product ), processInfo['flux']['data'], processInfo['workDir']

    nProj = len( projectileGroupBoundaries ) - 1
    nProd = len( productGroupBoundaries ) - 1
    TM_1, TM_E = [], []
    for i in xrange( nProj ) :
        rowTM = nProd * [ 0. ]
        TM_1.append( rowTM )
        rowTM = nProd * [ 0. ]
        TM_E.append( rowTM )
    indexEp = -1
    for Ep in productGroupBoundaries :
        if( gammaEnergy <= Ep ) : break
        indexEp += 1
    indexEp = max( indexEp, 0 )
    indexEp = min( indexEp, nProd - 1 )
    groupedFlux = fluxes[0].groupOneFunction( projectileGroupBoundaries )
    xsec = crossSection.groupTwoFunctions( projectileGroupBoundaries, fluxes[0], norm = groupedFlux )
    for indexE in xrange( nProj ) :
        x = multiplicity * xsec[indexE]
        TM_1[indexE][indexEp] = x
        TM_E[indexE][indexEp] = x * gammaEnergy

    return( [ [ TM_1 ], [ TM_E ] ] )

def primaryGammaAngularData( processInfo, projectile, product, bindingEnergy, massRatio, crossSection, angularData, multiplicity = 1, comment = None ) :
    """Currently, only isotropic (i.e., l = 0) data are returned. That is, lMax and angularData are ignored. massRatio is the
    target mass divide by the sum of the projectile and target masses.

    This function perform the integration 

        TM = \int_g dE \int_h dE' S(E) M(E) f(E) P(E -> E') / \int_g dE f(E)

    where \int_g is the integral of E from E_i to E_{i+1}, \int_g is the integral of E' from E'_j to E'_{j+1}, S(E) is the cross section,
    M(E) is the products multiplicity, f(E) is the flux weighting and P(E -> E') is the probability for a projectile of energy E producing
    a primary gamma of energy E' for binding energy bindingEnergy. This function assumes that M(E) is a constant. For primary gamma's captured
    into binding energy bindingEnergy P(E -> E') = deltaFunction( E' - ( bindingEnergy + massRatio E ) ) where massRatio = mt / ( mp + mt ),
    mp is the projectile's mass and mt is the target's mass. Note, this formula is from the ENDF manual an ignores the recoil of the residual
    nucleus."""

    logFile, lMax, projectileGroupBoundaries, productGroupBoundaries, fluxes, workDir = processInfo['logFile'], processInfo.getParticle_lMax( product ), \
        processInfo.getParticleGroups( projectile ), processInfo.getParticleGroups( product ), processInfo['flux']['data'], processInfo['workDir']

    nProj = len( projectileGroupBoundaries ) - 1
    nProd = len( productGroupBoundaries ) - 1
    TM_1, TM_E = [], []
    for i in xrange( nProj ) :
        TM_1.append( nProd * [ 0. ] )
        TM_E.append( nProd * [ 0. ] )
    Eg2 = bindingEnergy + massRatio * projectileGroupBoundaries[0]
    for indexEo, Eo in enumerate( productGroupBoundaries ) :
        if( Eg2 <= Eo ) : break
    indexEo = min( max( indexEo - 1, 0 ), nProd - 1 )
    groupedFlux = fluxes[0].groupOneFunction( projectileGroupBoundaries )
    EMin, EMax = crossSection.domainMin( ), crossSection.domainMax( )
    axes_ = axes.defaultAxes( labelsUnits = { 0 : crossSection.axes[0].getUnit( ), 1 : crossSection.axes[0].getUnit( ) } )
    Egp = XYs.XYs( axes_, [ [ EMin, bindingEnergy + massRatio * EMin ], [ EMax, bindingEnergy + massRatio * EMax ] ], accuracy = 1e-12 )
    for indexEi in xrange( nProj ) :
        Ei2 = projectileGroupBoundaries[indexEi + 1]
        Eg2 = bindingEnergy + massRatio * Ei2
        EiMin = projectileGroupBoundaries[indexEi]
        while( True ) :
            incrementIndexEo, EiMax = 0, Ei2
            if( indexEo < ( nProd - 1 ) ) :
                if( Eg2 > productGroupBoundaries[indexEo + 1] ) :
                    incrementIndexEo = 1
                    EiMax = ( productGroupBoundaries[indexEo + 1] - bindingEnergy ) / massRatio
            TM_1[indexEi][indexEo] = crossSection.integrateTwoFunctions( fluxes[0], xMin = EiMin, xMax = EiMax ) / groupedFlux[indexEi]
            TM_E[indexEi][indexEo] = crossSection.integrateThreeFunctions( fluxes[0], Egp, xMin = EiMin, xMax = EiMax ) / groupedFlux[indexEi]
            if( incrementIndexEo == 0 ) : break
            EiMin = EiMax
            if( indexEo < ( nProd - 1 ) ) : indexEo += 1
    return( [ [ TM_1 ], [ TM_E ] ] )

def NBodyPhaseSpace( processInfo, projectile, product, crossSection, numberOfProducts, mProj, mTarget, mProduct, mTotal, Q, multiplicity = 1, comment = None ) :

    logFile, lMax = processInfo['logFile'], processInfo.getParticle_lMax( product )
    projectileGroupBoundaries, productGroupBoundaries = processInfo.getParticleGroups( projectile ), processInfo.getParticleGroups( product )
    fluxes, workDir = processInfo['flux']['data'], processInfo['workDir']

    s  = versionStr + '\n'
    s += "Process: phase space spectrum\n"
    if( comment is not None ) : s += "Comment: %s\n" % comment
    s += "Number of particles: %s\n" % numberOfProducts
    s += "Projectile's mass: %s\n" % floatToString( mProj )
    s += "Target's mass: %s\n" % floatToString( mTarget )
    s += "Product's mass: %s\n" % floatToString( mProduct )
    s += "Total mass: %s\n" % mTotal
    s += "Q value: %s\n" % Q
    s += 'Quadrature method: Gauss6\n'
    s += commonDataToString( lMax, projectileGroupBoundaries, productGroupBoundaries, fluxes, crossSection, multiplicity = multiplicity )
    return( executeCommand( logFile, '%s/get_transfer' % srcPath, s, workDir ) )

def KalbachMann_TransferMatrix( processInfo, projectile, product, crossSection, particlesData, KalbachMannData, multiplicity = 1, comment = None ) :

    logFile, lMax = processInfo['logFile'], processInfo.getParticle_lMax( product )
    projectileGroupBoundaries, productGroupBoundaries = processInfo.getParticleGroups( projectile ), processInfo.getParticleGroups( product )
    fluxes, workDir = processInfo['flux']['data'], processInfo['workDir']

    energy_in_unit = 'MeV'
    s  = versionStr + '\n'
    s += "Process: Kalbach spectrum\n"
    if( comment is not None ) : s += "Comment: %s\n" % comment
    s += "Projectile's mass: %s\n" % particlesData['projectile']['mass']
    s += "Projectile's ZA: %s\n" % particlesData['projectile']['ZA']
    s += "Target's mass: %s\n" % particlesData['target']['mass']
    s += "Target's ZA: %s\n" % particlesData['target']['ZA']
    s += "Product's mass: %s\n" % particlesData['product']['mass']
    s += "Product's ZA: %s\n" % particlesData['product']['ZA']
    s += "Residual's mass: %s\n" % particlesData['residual']['mass']
    s += "Compound's mass: %s\n" % particlesData['compound']['mass']
    s += "Quadrature method: adaptive\n"
    s += commonDataToString( lMax, projectileGroupBoundaries, productGroupBoundaries, fluxes, crossSection, multiplicity = multiplicity, 
        energy_in_unit = energy_in_unit )
    s += KalbachMannDataToString( KalbachMannData, energy_in_unit )
    return( executeCommand( logFile, '%s/get_transfer' % srcPath, s, workDir ) )

def executeCommand( logFile, file, cmd, workDir ) :

    def checkNegative_l0( TM_l0, weight, f ) :

        negative_l0Counter, largestNegative = 0, 0.
        for ig in sorted( TM_l0.keys( ) ) :
            row = TM_l0[ig]
            for ih, cell in enumerate( row ) :
                if( cell < 0 ) :
                    negative_l0Counter += 1
                    largestNegative = min( largestNegative, cell )
                    f.write( 'negative l=0 for weight %s at row %3d column %3d: %s\n' % ( weight, ig, ih, cell ) )
        return( negative_l0Counter )
            
    dataFile = fudgeFileMisc.fudgeTempFile( dir = workDir, suffix = "_in" )
    dataFile.write( cmd )
    dataFile.close( )
    infoFile = '%s.info' %  dataFile.getName( )
    t0 = times.times( )
    try :
        status, stdout, stderr = subprocessing.executeCommand( [ file, '-output', '%s.out' % dataFile.getName( ), dataFile.getName( ) ], 
            stdout = infoFile, stderr = subprocess.STDOUT )
    except :
        f = open( "%s.err" % dataFile.getName( ), "w" )
        f.close( )
        raise
    f = open( infoFile, 'a' )
    f.write( str( t0 ) + "\n" )
    if( status != 0 ) :
        print "status = ", status
        print s
        raise Exception( 'Transfer failed for %s %s' % ( file, dataFile.getName( ) ) )
    TM1, TME = parseOutputFile( '%s.out' % dataFile.getName( ) )
    negative_l0Counter = checkNegative_l0( TM1[0], "0", f )
    negative_l0Counter += checkNegative_l0( TME[0], "E", f )
    if( negative_l0Counter > 0 ) : print 'WARNING: %d negative l=0 elements found in transfer matrix' % negative_l0Counter
    f.close( )
    if( workDir is None ) :
        os.remove( '%s.out' % dataFile.getName( ) )
        os.remove( '%s.info' % dataFile.getName( ) )
        dataFile.delete( )
    return( TM1, TME )

def parseOutputFile( file ) :

    def parseTransferMatrix( paras, index, ls, file ) :

        s = ls[index].split( ":" )                      # This line should start as 'Integrals, weight = ...'
        EBins = int( s[1].split( '=' )[1] )
        if( EBins != paras['numEinBins'] ) : raise Exception( 'Bad transfer numEinBins = %d != %d in file %s\n%s' % ( EBins, paras['numEinBins'], file, ls[index] ) )
        lMax = paras['lMax']
        lRange = xrange( lMax + 1 )
        TM = []
        for l in lRange : 
            TMl = {}
            for i in xrange( EBins ) : TMl[i] = paras['numEoutBins'] * [ 0. ]
            TM.append( TMl )
        for iE in xrange( EBins ) :
            index += 1
            s = ls[index].split( ":" )                  # This should start as 'EinBin = ...'
            EBin = int( s[0].split( '=' )[1] )
            if( not( 0 <= EBin < EBins ) ) : raise Exception( 'Bad transfer EinBin at index = %d in file %s\n%s' % ( index, file, ls[index] ) )
            EpBins = int( s[1].split( '=' )[1] )
            if( EpBins != paras['numEoutBins'] ) : raise Exception( 'Bad transfer numEoutBins at index = %d in file %s\n%s' % ( index, file, ls[index] ) )
            rowLTM = []
            for l in lRange : rowLTM.append( [] )
            for iEp in xrange( EpBins ) :
                index += 1
                rowLs = map( float, ls[index].split( ) )
                if( len( rowLs ) != ( lMax + 1 ) ) : raise Exception( 'Bad transfer lRows at index = %d in file %s\n%s' % ( index, file, ls[index] ) )
                for l in lRange : rowLTM[l].append( rowLs[l] )
            for l in lRange : TM[l][EBin] = rowLTM[l]
        return( TM, index )

    TM1 = None
    TME = None
    try :
        f = open( file )
    except :
        raise Exception( 'Could not open transfer file = %s' % file )
    ls = f.readlines( )
    f.close( )
    if( ls[0] != "xndfgenGetTransfer: version 1\n" ) : raise Exception( 'Bad transfer file version in file %s\n%s' % ( file, ls[0] ) )
    index = 1
    n = len( ls )
    paras = { 'lMax' : None, 'numEinBins' : None, 'numEoutBins' : None }
    while( index < n ) :
        s = ls[index].split( ":" )
        key = s[0]
        if( key in [ 'lMax', 'numEinBins', 'numEoutBins' ] ) :
            if( paras[key] is not None ) : raise Exception( "Bad transfer file, multiple key = %s lines in file %s" % ( key, file ) )
            try :
                paras[key] = int( s[1] )
            except :
                raise Exception( "Bad transfer file, could not convert data for key = %s to integer\n%s in file %s" % ( key, ls[index], file ) )
        elif( key == "Comment" ) : 
            pass
        elif( key == "Integrals, weight = 1" ) :
            if( TM1 is not None  ) : raise Exception( "Bad transfer file, multiple weight = 1 datasets in file %s" % file )
            TM1, index = parseTransferMatrix( paras, index, ls, file )
        elif( key == "Integrals, weight = E'" ) :
            if( TME is not None  ) : raise Exception( "Bad transfer file, multiple weight = E' datasets in file %s" % file )
            TME, index = parseTransferMatrix( paras, index, ls, file )
        else :
            raise Exception( 'Bad transfer file key in file %s\n%s' % ( file, ls[index] ) )
        index += 1
    return( TM1, TME )

def commonDataToString( lMax, projectileGroupBoundaries, productGroupBoundaries, fluxes, crossSection, multiplicity = None, energy_in_unit = None ) :

    crossSection_energy_in_unit = crossSection.axes[0].getUnit( )
    if( energy_in_unit is not None ) :
        factor = crossSection.axes[0].unitConversionFactor( energy_in_unit )
        crossSection = crossSection.convertAxisToUnit( 0, energy_in_unit )
        fluxes = [ f.convertAxisToUnit( 0, energy_in_unit ) for f in fluxes ]
        projectileGroupBoundaries = [ factor * b for b in projectileGroupBoundaries ]   # ????? Currently use endl1dmath projectileGroupBoundaries, need to
        productGroupBoundaries = [ factor * b for b in productGroupBoundaries ]         # have a new class that has convertAxisToUnit.
    s  = "lMax: %d\n" % lMax
    s += projGBToString( projectileGroupBoundaries )
    s += prodGBToString( productGroupBoundaries )
    s += fluxesToString( fluxes )
    if( not ( crossSection is None ) ) : s += crossSectionToString( crossSection )
    if( not ( multiplicity is None ) ) : s += multiplicityToString( multiplicity, crossSection[0][0], crossSection[-1][0], crossSection_energy_in_unit,
        energy_in_unit = energy_in_unit )
    return( s )

def projGBToString( gb ) :

    s  = "Projectile's group boundaries: n = %d\n" % len( gb )
    s += oneDToString( gb )
    return( s )

def prodGBToString( gb ) :

    s  = "Product's group boundaries: n = %d\n" % len( gb )
    s += oneDToString( gb )
    return( s )

def fluxesToString( fluxes ) :

    interpolationStr, s = twoDToString( "", fluxes[0], addHeader = False )
    s  = [ "Fluxes: n = %d" % len( fluxes ) ]
    s.append( interpolationStr )
    for l, ef in enumerate( fluxes ) :
        interpolationStr2, sp = twoDToString( "", ef, addHeader = False )
        if( interpolationStr != interpolationStr2 ) : raise Exception( 'interpolation = "%s" for l = %s is difference from l = 0 which is "%s"' %
            interpolationStr2, interpolationStr )
        s.append( '  n = %d' % len( sp  ) )
        s += sp
    s.append( '' )
    return( '\n'.join( s ) )

def crossSectionToString( crossSection ) :

    return( '\n'.join( twoDToString( "Cross section", crossSection ) ) )

def multiplicityToString( multiplicity_, EMin, EMax, crossSection_energy_in_unit, energy_in_unit = None ) :

    import fudge

    class template :
        """Dummy class needed for multiplicity.polynomial.toPointwiseLinear."""

        def __init__( self, EMin, EMax ) :

            self.EMin = EMin
            self.EMax = EMax

        def domainMin( self ) : return( self.EMin )

        def domainMax( self ) : return( self.EMax )

    if( energy_in_unit is not None ) : 
        factor = PhysicalQuantityWithUncertainty( 1, energy_in_unit ).getValueAs( crossSection_energy_in_unit )
        EMin, EMax = factor * EMin, factor * EMax

    try :
        multiplicity = XYs.XYs( fudge.gnd.productData.multiplicity.pointwise.defaultAxes( energyUnit = crossSection_energy_in_unit ), 
                [ [ EMin, multiplicity_.getConstant( ) ], [ EMax, multiplicity_.getConstant( ) ] ], 1e-6 )
    except :
        from fudge.gnd import tokens

        if( multiplicity_.getNativeDataToken( ) == tokens.referenceFormToken ) : 
            multiplicity_ = multiplicity_.forms[tokens.referenceFormToken].referenceInstance.multiplicity
        if( tokens.pointwiseFormToken in multiplicity_.forms ) :
            multiplicity = multiplicity_.forms[tokens.pointwiseFormToken]
        else :
            multiplicityForm = multiplicity_[multiplicity_.getNativeDataToken( )]
            if( isinstance( multiplicityForm, fudge.gnd.productData.multiplicity.piecewise ) ) :
                multiplicity = multiplicityForm.toPointwiseLinear( 1e-8, 1e-8 )
            elif( isinstance( multiplicityForm, fudge.gnd.productData.multiplicity.polynomial ) ) :
                multiplicity = multiplicityForm.toPointwiseLinear( 1e-8, 1e-8, template( EMin, EMax ) )
            else :
                raise Exception( 'multiplicity = %s is unsupported' % multiplicityForm.moniker )
    if( energy_in_unit is not None ) : multiplicity = multiplicity.convertAxisToUnit( 0, energy_in_unit )
    return( '\n'.join( twoDToString( "Multiplicity", multiplicity ) ) )

def angularToString( angularData, crossSection ) :

    from fudge.gnd.productData import distributions
    from fudge.core.utilities import brb

    if( isinstance( angularData, distributions.angular.isotropic ) ) :
        s  = "Angular data: n = 2\n"
        s += 'Incident energy interpolation: %s-%s\n' % ( axes.linearToken, axes.linearToken )
        s += 'Outgoing cosine interpolation: %s-%s\n' % ( axes.linearToken, axes.linearToken )
        for x in [ crossSection.domainMin( ), crossSection.domainMax( ) ] : s += " %s : n = 2\n   -1 0.5\n   1 0.5\n" % ( doubleFmt % x )

    elif( isinstance( angularData, distributions.angular.pointwise ) ) :
        s  = "Angular data: n = %s\n" % len( angularData )

        independent, dependent, qualifier = angularData.axes[0].interpolation.getInterpolationTokens( )
        if( qualifier is not None ) : raise Exception( 'interpolation qualifier = %s is not supported' % qualifier )
        if( ( independent not in [ axes.linearToken ] ) or ( dependent not in [ axes.linearToken, axes.flatToken ] ) ) :
            angularData = angularData.toPointwiseLinear( lowerEps = lowerEps, upperEps = upperEps )

        independent, dependent, qualifier = angularData.axes[1].interpolation.getInterpolationTokens( )
        if( qualifier is not None ) : raise Exception( 'interpolation qualifier = %s is not supported' % qualifier )
        if( ( independent not in [ axes.linearToken ] ) or ( dependent not in [ axes.linearToken, axes.flatToken ] ) ) :
            angularData = angularData.toPointwiseLinear( lowerEps = lowerEps, upperEps = upperEps )

        independent, dependent, qualifier = angularData.axes[0].interpolation.getInterpolationTokens( )
        s += 'Incident energy interpolation: %s-%s\n' % ( independent, dependent )
        independent, dependent, qualifier = angularData.axes[1].interpolation.getInterpolationTokens( )
        s += 'Outgoing cosine interpolation: %s-%s\n' % ( independent, dependent )

        s += threeDListToString( angularData )

    elif( isinstance( angularData, LegendreSeries.W_XYs_LegendreSeries ) ) :

        independent, dependent, qualifier = angularData.axes[0].interpolation.getInterpolationTokens( )
        if( qualifier is not None ) : raise Exception( 'interpolation qualifier = %s is not supported' % qualifier )
        if( ( independent not in [ axes.linearToken ] ) or ( dependent not in [ axes.linearToken, axes.flatToken ] ) ) :
            angularData = angularData.toLegendreLinear( 1e-6 )

        s = [ "Legendre coefficients: n = %s" % len( angularData ) ]

        independent, dependent, qualifier = angularData.axes[0].interpolation.getInterpolationTokens( )
        s.append( 'Interpolation: %s-%s' % ( independent, dependent ) )

        for energy_in in angularData :
            s.append( ' Ein: %s:  n = %d' % ( energy_in.value, len( energy_in ) ) )
            for coefficient in energy_in : s.append( "%s" % coefficient )
        s = '\n'.join( s ) + '\n'

    else :
        raise Exception( "angular data = %s not supported" % brb.getType( angularData ) )

    return( s )

def energyFunctionToString( energyData, weight = None ) :

    def getParameter( data, label ) :

        linlin = data.toPointwiseLinear( lowerEps = lowerEps, upperEps = upperEps )
        sData = [ "%s: n = %d" % ( label, len( linlin ) ) ]
        sData.append( 'Interpolation: linear-linear' )
        for x, y in linlin : sData.append( '%s %s' % ( x, y ) )
        return( sData )

    from fudge.gnd.productData.distributions import energy

    sData = []
    if( isinstance( energyData, ( energy.simpleMaxwellianFissionSpectrum, energy.evaporationSpectrum, energy.WattSpectrum ) ) ) :
        sData.append( 'U: %s' % energyData.U.getValueAs( energyData.parameter1.axes[0].getUnit( ) ) )

    parameter1 = energyData.parameter1.toPointwiseLinear( lowerEps = lowerEps, upperEps = upperEps )

    if( weight is None ) :
        weight = XYs.XYs( axes.defaultAxes( ), [ [ parameter1.domainMin( ), 1. ], [ parameter1.domainMax( ), 1. ] ], 1e-6 )
    else :
        weight = weight.weight
    sData += twoDToString( "weight", weight, addExtraBlankLine = False )

    if( isinstance( energyData, energy.generalEvaporationSpectrum ) ) :
        sProcess = 'Process: general evaporation\n'
        sSpecific = ''
        sData += getParameter( energyData.parameter1, 'Theta' )
        sData += getParameter( energyData.parameter2, 'g' )
    elif( isinstance( energyData, energy.simpleMaxwellianFissionSpectrum ) ) :
        sProcess = 'Process: Maxwell spectrum\n'
        sSpecific = ''
        sData += getParameter( energyData.parameter1, 'Theta' )
    elif( isinstance( energyData, energy.evaporationSpectrum ) ) :
        sProcess = 'Process: Evaporation spectrum\n'
        sSpecific = 'Interpolate Eout integrals: false\n'
        sSpecific += 'Quadrature method: Gauss6\n'
        sData += getParameter( energyData.parameter1, 'Theta' )
    elif( isinstance( energyData, energy.WattSpectrum ) ) :
        sProcess = 'Process: Watt spectrum\n'
        sSpecific = 'Interpolate Eout integrals: false\n'
        sSpecific += 'Conserve: number\n'
        sData += getParameter( energyData.parameter1, 'a' )
        sData += getParameter( energyData.parameter2, 'b' )
    elif( isinstance( energyData, energy.MadlandNix ) ) :
        sProcess = 'Process: Madland-Nix spectrum\n'
        sSpecific  = 'Quadrature method: adaptive\n'
        sSpecific += 'Interpolate Eout integrals: false\n'
        sSpecific += 'Conserve: number\n'
        sSpecific += 'EFL: %s\n' % energyData.EFL.getValueAs( energyData.parameter1.axes[0].getUnit( ) )
        sSpecific += 'EFH: %s\n' % energyData.EFH.getValueAs( energyData.parameter1.axes[0].getUnit( ) )
        sData += getParameter( energyData.parameter1, 'TM' )
    else :
        raise Exception( 'Unsupport energy functional = %s' % energyData.moniker )
    sData.append( '' )
    return( sProcess, sSpecific, '\n'.join( sData ) )

def EEpPDataToString( EEpPData ) :

    s = "EEpPData: n = %d\n" % len( EEpPData )

    independent, dependent, qualifier = EEpPData.axes[0].interpolation.getInterpolationTokens( )
    if( qualifier == axes.correspondingPointsToken ) :
        s += 'Incident energy interpolation: corresponding energies\n'
    elif( qualifier == axes.unitBaseToken ) :
        s += 'Incident energy interpolation: %s\n' % GND2ProcessingInterpolations[axes.unitBaseToken]
    else :
        s += 'Incident energy interpolation: %s\n' % GND2ProcessingInterpolations[axes.unitBaseToken]

    independent, dependent, qualifier = EEpPData.axes[1].interpolation.getInterpolationTokens( )
    if( dependent == axes.flatToken ) :
        s += 'Outgoing energy interpolation: flat\n'
    elif( independent == dependent == axes.linearToken ) :
        s += 'Outgoing energy interpolation: linear-linear\n'

    s += threeDToString( EEpPData )
    return( s )

def LEEpPDataToString( LEEpPData ) :

    s = "LEEpPData: n = %d\n" % len( LEEpPData )
    s += 'Incident energy interpolation: unit base\n'
    s += 'Outgoing energy interpolation: linear-linear\n'
    for l_EEpP in LEEpPData :
        EEpP = l_EEpP.EpP
        s += "  l = %d: n = %d\n" % ( l_EEpP.l, len( EEpP ) )
        s += threeDToString( EEpP )
    return( s )

def LegendreDataToString( LegendreData ) :

    s = [ 'Legendre data by incident energy:  n = %d' % len( LegendreData ) ]
    independent, dependent, qualifier = LegendreData.axes[0].interpolation.getInterpolationTokens( )
    if( qualifier == axes.correspondingPointsToken ) :
        s.append( 'Incident energy interpolation: corresponding energies' )
    elif( qualifier == axes.unitBaseToken ) :
        s.append( 'Incident energy interpolation: %s' % GND2ProcessingInterpolations[axes.unitBaseToken] )
    independent, dependent, qualifier = LegendreData.axes[1].interpolation.getInterpolationTokens( )
    if( dependent == axes.flatToken ) :
        s.append( 'Outgoing energy interpolation: flat' )
    elif( independent == dependent == axes.linearToken ) :
        s.append( 'Outgoing energy interpolation: linear-linear' )
    for EEpCl in LegendreData :
        s.append( " Ein: %s:  n = %d" % ( EEpCl.value, len( EEpCl ) ) )
        for EpCl in EEpCl :
            s.append( "  Eout: %s:  n = %d" % ( EpCl.value, len( EpCl ) ) )
            for Cl in EpCl : s.append( "   %s" % Cl )
    s.append( '' )
    return( '\n'.join( s ) )

def EMuEpPDataToString( EMuEpPData ) :

    s = "EMuEpPData: n = %d\n" % len( EMuEpPData )

    independent, dependent, qualifier = EMuEpPData.axes[0].interpolation.getInterpolationTokens( )
    if( qualifier not in [ None, axes.unitBaseToken ] ) : raise Exception( 'interpolation qualifier = %s is not supported' % qualifier )
    if( ( independent not in [ axes.linearToken ] ) or ( dependent not in [ axes.linearToken, axes.flatToken ] ) ) :
            raise Exception( 'interpolation independent = %s or dependent = %s is not supported' % ( independent, dependent ) )
    s += 'Incident energy interpolation: %s\n' % GND2ProcessingInterpolations[axes.unitBaseToken]

    s += 'Outgoing cosine interpolation: %s\n' % GND2ProcessingInterpolations[axes.unitBaseToken]

    s += 'Outgoing energy interpolation: %s-%s\n' % ( axes.linearToken, axes.linearToken )

    for muEpP in EMuEpPData :
        E = muEpP.value
        s += "  E = %s: n = %d\n" % ( doubleFmt % E, len( muEpP ) )
        s += threeDToString( muEpP )
    return( s )

def KalbachMannDataToString( KalbachMannData, energy_in_unit ) :

    nForm = 3
    if( KalbachMannData.form == 'fra' ) : nForm = 4
    s = "KalbachData: n = %s\n" % len( KalbachMannData )
    s += "Incident energy interpolation: %s\n" % GND2ProcessingInterpolations[axes.unitBaseToken]
    s += getOutgoingEnergyInterpolation( KalbachMannData.axes[1], [ axes.linearToken ], [ axes.linearToken, axes.flatToken ] )
    factor = KalbachMannData.axes[0].unitConversionFactor( energy_in_unit ) # Energies passed to get_transfer must be in energy_in_unit, traditionally MeV.
    factors = [ factor, 1. / factor ] + ( nForm - 2 ) * [ 1. ]
    for energies_in in KalbachMannData :
        s += 'Ein: %s:  n = %d\n' % ( doubleFmt % ( energies_in.value * factor ), len( energies_in ) / nForm )
        for iEp in xrange( 0, len( energies_in.coefficients ), nForm ) :
            cs = ''
            for ic in xrange( nForm ) : cs += ' %s' % ( doubleFmt % ( energies_in.coefficients[iEp + ic] * factors[ic] ) )
            s += cs + '\n'
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

def oneDToString( data ) :

    i = 1
    a = [ "" ]
    for d in data :
        a.append( doubleFmt % d )
        if( ( i % 10 ) == 0 ) : a.append( "\n" )
        i += 1
    if( ( i % 10 ) != 1 ) : a.append( "\n" )
    return( " ".join( a ) )

def twoDToString( label, data, addHeader = True, addExtraBlankLine = True ) :

    if( isinstance( data, XYs.XYs ) ) :
        independent, dependent, qualifier = data.axes[0].interpolation.getInterpolationTokens( )
        if( qualifier is not None ) : raise Exception( 'interpolation qualifier = %s is not supported' % qualifier )
        if( ( independent not in [ axes.linearToken ] ) or ( dependent not in [ axes.linearToken, axes.flatToken ] ) ) :
            data = data.toPointwiseLinear( lowerEps = lowerEps, upperEps = upperEps )
        independent, dependent, qualifier = data.axes[0].interpolation.getInterpolationTokens( )
        if( dependent == axes.flatToken ) :
            interpolationStr = 'Interpolation: %s' % dependent
        else :
            interpolationStr = 'Interpolation: %s-%s' % ( independent, dependent )
    else :
        from fudge.core.utilities import brb
        brb.objectoutline( data )
        raise 'hell'
    fmt = " %s %s" % ( doubleFmt, doubleFmt )
    a = [ fmt % ( x, y ) for x, y in data ]
    if( not( addHeader ) ) : return( interpolationStr, a )
    a.insert( 0, interpolationStr )
    a.insert( 0, "%s: n = %d" % ( label, len( data ) ) )
    if( addExtraBlankLine ) : a.append( '' )
    return( a )

def threeDToString( data ) :

    a = []
    wFmt = " %s : n = %%d" % ( doubleFmt )
    xyFmt = "   %s %s" % ( doubleFmt, doubleFmt )
    for w_xy in data :
        if( type( w_xy ) == type( [] ) ) :       # ???? This is a kludge until Legendre EEpP data is represented like angular EMuP data.
            raise Exception( 'need to update for new gnd forms' )
            x, yz = xyz
        else :
            w, xy = w_xy.value, w_xy
        a.append( wFmt % ( w, len( xy ) ) )
        for x, y in xy : a.append( xyFmt % ( x, y ) )
    a.append( "" )
    return( "\n".join( a ) )

def threeDListToString( data ) :

    a = []
    wFmt = " %s : n = %%d" % ( doubleFmt )
    xyFmt = "   %s %s" % ( doubleFmt, doubleFmt )
    for w_xy in data :
        w = w_xy.value
        a.append( wFmt % ( w, len( w_xy ) ) )
        for x, y in w_xy : a.append( xyFmt % ( x, y ) )
    a.append( "" )
    return( "\n".join( a ) )

def getInterpolationString( axis, allowedIndependent, allowedDependent ) :

    independent, dependent, qualifier = axis.interpolation.getInterpolationTokens( )
    if( qualifier == axes.unitBaseToken ) :
        return( GND2ProcessingInterpolations[axes.unitBaseToken] )
    elif( qualifier is None ) :
        if( independent not in allowedIndependent ) : raise Exception( 'independent interpolation = "%s" not in allowed "%s"' % 
            ( independent, allowedIndependent ) )
        if( dependent not in allowedDependent ) : raise Exception( 'dependent interpolation = "%s" not in allowed "%s"' % 
            ( dependent, allowedDependent ) )
        if( dependent == axes.flatToken ) :
            return( GND2ProcessingInterpolations[axes.flatToken] )
        elif( dependent not in [ axes.linearToken, axes.logToken ] ) :
            raise Exception ( 'Need to add support for dependent interpolation = %s' % dependent )
        else :
            if( independent not in [ axes.linearToken, axes.logToken ] ) :
                raise Exception ( 'Need to add support for independent interpolation = %s' % independent )
            return( GND2ProcessingInterpolations[independent+dependent] )
    else :
        raise Exception( 'Need to add support for interpolation qualifier = %s' % qualifier )
    
def getOutgoingEnergyInterpolation( axis, allowedIndependent, allowedDependent ) :

    return( "Outgoing energy interpolation: %s\n" % getInterpolationString( axis, allowedIndependent, allowedDependent ) )

def addTMs( TMs ) :

    def addTM_l( TM_l, tm_l ) :

        for iRow in TM_l :
            TM_lRow, tm_lRow = TM_l[iRow], tm_l[iRow]
            for iColumn in xrange( len( TM_lRow ) ) : TM_lRow[iColumn] += tm_lRow[iColumn]

    TM = []
    for i, tm in enumerate( TMs ) :     # Find the TM with the largest Legendre order.
        if( tm is not None ) :
            if( len( tm ) > len( TM ) ) : index, TM = i, tm
    if( len( TM ) == 0 ) : return( None )
    for i, tm in enumerate( TMs ) :
        if( i == index ) : continue
        if( tm is not None ) :
            for l, tm_l in enumerate( tm ) : addTM_l( TM[l], tm_l )
    return( TM )
