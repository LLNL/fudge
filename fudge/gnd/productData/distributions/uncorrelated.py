# <<BEGIN-copyright>>
# <<END-copyright>>

"""Uncorrelated double differential distribution classes"""

import base
import miscellaneous
import angular
import energy
import fudge
from fudge.core.math.xData import axes, LegendreSeries, W_XYs

__metaclass__ = type

class component( base.component ) :
    """ Contains uncorrelated distributions for outgoing angle and outgoing energy. Use when the correlations
    between angle and energy are unknown."""

    def __init__( self, angularComponent, energyComponent ) :

        nativeData = '%s=%s : %s=%s' % ( base.angularComponentToken, angularComponent.nativeData, base.energyComponentToken, energyComponent.nativeData )
        base.component.__init__( self, base.uncorrelatedComponentToken, nativeData )
        self.angularComponent = angularComponent
        self.angularComponent.setParent( self )
        self.energyComponent = energyComponent
        self.energyComponent.setParent( self )

    def calculateDepositionData( self, processInfo, tempInfo ) :

        def fixLimits( multiplicityLimits, Es ) :

            if( Es[0] < multiplicityLimits[0] ) :
                while( len( Es ) and ( Es[0] < multiplicityLimits[0] ) ) : del Es[0]
                if( len( Es ) and ( Es[0] > multiplicityLimits[0] ) ) : Es.insert( 0, multiplicityLimits[0] )
            if( Es[-1] > multiplicityLimits[1] ) :
                while( len( Es ) and ( Es[-1] > multiplicityLimits[1] ) ) : del Es[-1]
                if( len( Es ) and ( Es[-1] < multiplicityLimits[1] ) ) : Es.append( multiplicityLimits[1] )

        depData = []
        angularForm = self.angularComponent.forms[self.angularComponent.nativeData]
        energyForm = self.energyComponent.forms[self.energyComponent.nativeData]
        if( hasattr( energyForm, 'calculateDepositionData' ) ) : return( energyForm.calculateDepositionData( processInfo, tempInfo ) )
        multiplicity = tempInfo['multiplicity']
        energyAccuracy, momentumAccuracy = 1e-6, 1e-3
        energyUnit = tempInfo['incidentEnergyUnit']
        momentumDepositionUnit = energyUnit + '/c'

        Es = energyForm.getEnergyArray( tempInfo['EMin'], tempInfo['EMax'] )
        multiplicityLimits = multiplicity.getEnergyLimits( tempInfo['EMin'], tempInfo['EMax'] )
        fixLimits( multiplicityLimits, Es )
        depEnergy = [ [ E, multiplicity.getValue( E ) * energyForm.EpAverageAtE( E ) ] for E in Es ]
        axes_ = fudge.gnd.productData.energyDeposition.pointwise.defaultAxes( energyUnit = energyUnit, energyDepositionUnit = energyUnit )
        depData.append( fudge.gnd.productData.energyDeposition.pointwise( axes_, depEnergy, energyAccuracy ) )

        if( angularForm.moniker == base.isotropicFormToken ) :
            depMomentum = [ [ depEnergy[0][0], 0 ], [ depEnergy[-1][0], 0. ] ]
        else :
            for E in angularForm.getEnergyArray( tempInfo['EMin'], tempInfo['EMax'] ) :
                if( E not in Es ) : Es.append( E )
            Es.sort( )
            fixLimits( multiplicityLimits, Es )
            depMomentum = []
            for E in Es :
                if( tempInfo['product'].getName( ) == 'gamma' ) :
                    if( energyForm.moniker == base.constantFormToken ) :
                        pp = float( energyForm.data )
                    else :
                        pp = energyForm.EpAverageAtE( E )
                else :
                    if( hasattr( energyForm, 'sqrtEp_AverageAtE' ) ) :
                        pp = energyForm.sqrtEp_AverageAtE( E )
                    else :
                        pp = energyForm.getAtEnergy( E ).integrateWithWeight_sqrt_x( )
                muAverage = angularForm.muAverageAtEnergy( E )
                if( muAverage == 0. ) :                 # Handles isotropic cases where sqrtEp_AverageAtE is not fully implemented and returns None.
                    depMomentum.append( [ E, 0. ] )
                else :
                    if( pp is None ) : raise Exception( "sqrtEp_AverageAtE not fully implemented for %s" % energyForm.moniker )
                    depMomentum.append( [ E, multiplicity.getValue( E ) * pp * angularForm.muAverageAtEnergy( E ) ] )
        axes_ = fudge.gnd.productData.momentumDeposition.pointwise.defaultAxes( energyUnit = energyUnit, momentumDepositionUnit = momentumDepositionUnit )
        depData.append( fudge.gnd.productData.momentumDeposition.pointwise( axes_, depMomentum, momentumAccuracy ) )

        return( depData )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        from fudge.processing.deterministic import transferMatrices

        newComponents = []
        angularForm = self.angularComponent.forms[self.angularComponent.nativeData]
        energyForm = self.energyComponent.forms[self.energyComponent.nativeData]
        energyUnit = tempInfo['incidentEnergyUnit']
        massUnit = energyUnit + '/c**2'
        if( 'LLNL_MC' in processInfo['styles'] ) :

            if( not( base.isotropicFormToken in self.angularComponent.forms ) ) :
                if( angularForm.moniker != base.linearFormToken ) :
                    angularForm_ = angularForm
                    if( angularForm.moniker == base.pointwiseFormToken ) :
                        if( not( angularForm_.axes.isLinear( flatIsOk = True ) ) ) : self.angularComponent.addForm( angularForm_.toPointwiseLinear( 1e-3 ) )
                        angularForm_ = None
                    if( angularForm.moniker == base.LegendrePointwiseFormToken ) :
                        self.angularComponent.addForm( angularForm_.toPointwiseLinear( 1e-3 ) )
                        angularForm_ = None
                    if( angularForm_ is not None ) : raise Exception( 'Not implemented: %s' % angularForm.moniker )

            energyFormKeys = self.energyComponent.forms.keys( )
            if( not( base.linearFormToken in energyFormKeys ) ) :
                energyForm_ = energyForm
                if( base.pointwiseFormToken in energyFormKeys ) :
                    if( not( energyForm_.axes.isLinear( flatIsOk = True, qualifierOk = True ) ) ) : self.energyComponent.addForm( energyForm_.toPointwiseLinear( 1e-3 ) )
                    energyForm_ = None
                if( energyForm_ is not None ) :
                    if( energyForm_.moniker in [ base.weightedFunctionalsFormToken, base.MadlandNixFormToken, base.NBodyPhaseSpaceFormToken ] ) :
                        energyForm_ = None
                    elif( energyForm_.moniker in [ base.semiPiecewiseFormToken, base.simpleMaxwellianFissionFormToken, \
                                base.evaporationFormToken, base.WattFormToken ] ) :
                        self.energyComponent.addForm( energyForm_.toPointwiseLinear( 1e-3, lowerEps = 1e-7, upperEps = 1e-7  ) )
                        energyForm_ = None
                    elif( energyForm_.moniker in [ base.generalEvaporationFormToken ] ) :
                        if( energyForm_.isLinear( flatIsOk = True, qualifierOk = True ) ) : energyForm_ = None
                    if( energyForm_ is not None ) : raise Exception( 'Not implemented: %s' % energyForm_.moniker )
            
        if( 'LLNL_Pn' in processInfo['styles'] ) :
            if( processInfo['verbosity'] >= 30 ) : print '%sGrouping %s' % ( verbosityIndent, self.moniker )
            crossSection = tempInfo['crossSection']
            outputChannel = tempInfo['outputChannel'].outputChannel
            projectile, target, product = tempInfo['reactionSuite'].projectile, tempInfo['reactionSuite'].target, tempInfo['product']
            projectileName, productName = processInfo.getProjectileName( ), product.particle.getToken( )
            if( energyForm.moniker == base.constantFormToken ) :
                if( product.getName( ) == 'gamma' ) :
                    Ep = float( energyForm.data )
                    TM_1, TM_E = transferMatrices.discreteGammaAngularData( processInfo, projectileName, productName, Ep, crossSection,
                        angularForm, 1., comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )
                else :
                    raise Exception( 'See Bret' )
            else :
                if( isinstance( energyForm, energy.NBodyPhaseSpace ) ) :
                    mass1, mass2, mass3 = projectile.getMass( massUnit ), target.getMass( massUnit ), product.getMass( massUnit )
                    totalMass = energyForm.numberOfProductsMasses.getValueAs( massUnit )
                    Q = tempInfo['outputChannel'].getQ( energyUnit, final = False, groundStateQ = True )
                    TM_1, TM_E = transferMatrices.NBodyPhaseSpace( processInfo, projectileName, productName, crossSection, \
                        energyForm.numberOfProducts, mass1, mass2, mass3, totalMass, Q, tempInfo['multiplicity'], \
                        comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )
                else :
                    from fudge.gnd.productData.distributions import angular
                    TM_1, TM_E = transferMatrices.uncorrelated_EMuP_EEpP_TransferMatrix( processInfo, projectileName, productName, \
                        crossSection, angularForm, energyForm, tempInfo['multiplicity'], \
                        comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )
            fudge.gnd.miscellaneous.TMs2Form( processInfo, tempInfo, newComponents, TM_1, TM_E, crossSection.axes )

        return( newComponents )

    def toPointwiseLinear( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        import angularEnergy

        angularForm = self.angularComponent.toPointwiseLinear( accuracy = accuracy, lowerEps = lowerEps, upperEps = upperEps )
        energyForm = self.energyComponent.toPointwiseLinear( accuracy = accuracy, lowerEps = lowerEps, upperEps = upperEps )
        grid = angularForm.getDomainGrid( )
        for e in energyForm.getDomainGrid( ) :
            if( e not in grid ) : grid.append( e )
        grid.sort( )
        axes_ = angularEnergy.pointwise.defaultAxes( energyUnit = energyForm.axes[0].getUnit( ), energy_outUnit = energyForm.axes[1].getUnit( ),
            probabilityUnit = energyForm.axes[2].getUnit( ) )
        f_E_mu_Ep = angularEnergy.pointwise( axes_ )
        axesW_XY = axes.referenceAxes( f_E_mu_Ep, dimension = 3 )
        axesXY = axes.referenceAxes( f_E_mu_Ep, dimension = 2 )
        for e in grid :
            w_xys = W_XYs.W_XYs( axesW_XY, value = e )
            af = angularForm.interpolateAtW( e )
            ef = energyForm.interpolateAtW( e, unitBase = True )
            for mu, P in af :
                efp = P * ef
                efp.value = mu
                w_xys.append( efp )
            f_E_mu_Ep.append( w_xys )
        return( f_E_mu_Ep )

    def toXMLList( self, indent = "" ) :

        indent2 = indent + '  '
        endfFlag = ''
        xmlString = [ '%s<%s%s>' % ( indent, self.moniker, endfFlag ) ]
        xmlString += self.angularComponent.toXMLList( indent = indent2 )
        xmlString += self.energyComponent.toXMLList( indent = indent2 )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    def toENDF6( self, MT, endfMFList, flags, targetInfo ) :
        """In ENDF MF=6, some distributions should really be treated as uncorrelated: NBodyPhaseSpace, and also 
        Legendre expansions when only L=0 is listed.
        For GND we split these into uncorrelated angular (isotropic) and energy distributions.
        Must put back in original format when writing back to ENDF, however: """

        if( self.energyComponent.nativeData == base.NBodyPhaseSpaceFormToken ) :
            assert self.angularComponent.nativeData == base.isotropicFormToken
            form = self.energyComponent[ self.energyComponent.nativeData ]
            newComponent = fudge.gnd.productData.distributions.energyAngular.component( form.moniker )
            newComponent.addForm( form )
            newComponent.toENDF6( MT, endfMFList, flags, targetInfo )
        elif( targetInfo['product'].getAttribute( 'ENDFconversionFlag' ) == 'MF6' ) :
            import Legendre
            if( self.angularComponent.nativeData == base.isotropicFormToken ) :
                form = self.energyComponent.forms[ self.energyComponent.nativeData ]
                axis0, axis1, axis2 = form.axes[0], form.axes[1], form.axes[2]
                independent0, dependent0, qualifier0 = axis0.interpolation.getInterpolationTokens( )
                independent1, dependent1, qualifier1 = axis1.interpolation.getInterpolationTokens( )
                axes_ = Legendre.pointwise.defaultAxes( independent0, dependent0, energyInterpolationQualifier = qualifier0, frame = axis2.frame, 
                    energy_outInterpolation = independent1, C_lInterpolation = dependent1, energyUnit = axis0.getUnit( ), 
                    energy_outUnit = axis1.getUnit( ), C_lUnit = axis2.getUnit( ) )
                form2 = Legendre.pointwise( axes_ )          # change to LegendrePointwise.
                axes__ = axes.referenceAxes( form2 )

                for i, e_in in enumerate( form ) :
                    w_xys_LegendreSeries = LegendreSeries.W_XYs_LegendreSeries( axes__, index = i, value = e_in.value )
                    xys = e_in.copyDataToXYs( xUnit = 'eV', yUnit = '1/eV' )
                    for j, e_out_C_0 in enumerate( xys ) :
                        w_xys_LegendreSeries[j] = LegendreSeries.XYs_LegendreSeries( axes_[2].getUnit( ), e_out_C_0[1:], index = j, value = e_out_C_0[0] )
                    form2[i] = w_xys_LegendreSeries
                        # gamma is special case. Must gather all (primary, discrete and continuum) gamma info before writing back to ENDF.
                if( targetInfo.dict.get( "gammaToENDF6" ) ) : return( form2 )
                newComponent = Legendre.component( form2.moniker )
                newComponent.addForm( form2 )
                newComponent.toENDF6( MT, endfMFList, flags, targetInfo )
            elif( ( self.energyComponent.nativeData == base.pointwiseFormToken ) and ( self.angularComponent.nativeData == base.pointwiseFormToken ) ) :
                from fudge.legacy.converting import endfFormats, gndToENDF6

                energyForm = self.energyComponent.forms[ self.energyComponent.nativeData ]
                angularForm = self.angularComponent.forms[ self.angularComponent.nativeData ]
                LANG, LEP = 12, 2
                independent, dependent, qualifier = energyForm.axes[1].interpolation.getInterpolationTokens( )
                if( dependent == axes.flatToken ) : LEP = 1  # interpolation for E_out
                independent, dependent, qualifier = energyForm.axes[1].interpolation.getInterpolationTokens( )
                if( dependent == axes.flatToken ) :
                    LANG = 11
                elif( dependent == axes.logToken ) :
                    LANG = 14
                MF6 = [ endfFormats.endfContLine( 0, 0, LANG, LEP, 1, len( energyForm ) ) ]
                EInInterpolation = gndToENDF6.axisToEndfInterpolationFlag( energyForm.axes[0] )
                MF6 +=  endfFormats.endfInterpolationList( [ len( energyForm ), EInInterpolation ] )
                for indexE, EEpP in enumerate( energyForm ) :
                    EMuP = angularForm[indexE]
                    if( EEpP.value != EMuP.value ) : raise Exception( "EEpP.value = %s != EMuP.value = %s" % ( EEpP.value, EMuP.value ) )
                    NA, NEP = 2 * len( EMuP ), len( EEpP )
                    MF6.append( endfFormats.endfContLine( 0, EEpP.value, 0, NA, NEP * ( NA + 2 ), NEP ) )
                    data = []
                    for EpP in EEpP :
                        data = [ EpP[0], EpP[1] ]
                        for muP in EMuP : data += muP
                        MF6 += endfFormats.endfDataList( data )
                LAW, frame = 1, energyForm.axes[1].frame
                gndToENDF6.toENDF6_MF6( MT, endfMFList, flags, targetInfo, LAW, frame, MF6 )
            else :
                raise Exception( 'uncorrelated.toENDF6 not supported for energy form = %s and angular form = %s' %
                    ( self.energyComponent.nativeData, self.angularComponent.nativeData ) )
        else :                          # original data is in uncorrelated form
            self.angularComponent.toENDF6( MT, endfMFList, flags, targetInfo )
            self.energyComponent.toENDF6( MT, endfMFList, flags, targetInfo )

def parseXMLNode( uncorrelatedElement, linkData = {} ):
    """Translate <uncorreleted> element from xml. Must contain one <angular> and one <energy> component. """
    try:
        angular_ = angular.parseXMLNode( uncorrelatedElement.find(base.angularComponentToken) )
    except Exception as e:
        raise Exception, '/uncorrelated/angular%s' % e
    try:
        energy_ = energy.parseXMLNode( uncorrelatedElement.find(base.energyComponentToken) )
    except Exception as e:
        raise Exception, '/uncorrelated/energy%s' % e
    uncorrelated = component( angular_, energy_ )
    return uncorrelated
