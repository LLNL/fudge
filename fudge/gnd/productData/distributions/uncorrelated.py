# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# <<END-copyright>>

"""Uncorrelated double differential distribution classes"""

import math
import base, angular, energy, miscellaneous
import fudge
from fudge.core.math.xData import axes, LegendreSeries, W_XYs
from fudge.core.math import fudgemath

__metaclass__ = type

class component( base.component ) :
    """ Contains uncorrelated distributions for outgoing angle and outgoing energy. Use when the correlations
    between angle and energy are unknown."""

    def __init__( self, angularFormOrComponent, energyFormOrComponent ) :

        # if forms are supplied, wrap them in a component:
        if isinstance(angularFormOrComponent, angular.form):
            AC = angular.component( )
            AC.addForm( angularFormOrComponent, isNativeData=True )
            angularFormOrComponent = AC
        if isinstance(energyFormOrComponent, energy.form):
            EC = energy.component( )
            EC.addForm( energyFormOrComponent, isNativeData=True )
            energyFormOrComponent = EC

        if not isinstance(angularFormOrComponent, angular.component):
            raise TypeError("Needed angular distribution form or component, got %s" % type(angularFormOrComponent))
        if not isinstance(energyFormOrComponent, energy.component):
            raise TypeError("Needed energy distribution form or component, got %s" % type(energyFormOrComponent))

        angularComponent = angularFormOrComponent
        energyComponent = energyFormOrComponent

        nativeData = '%s=%s : %s=%s' % ( base.angularComponentToken, angularComponent.nativeData,
                base.energyComponentToken, energyComponent.nativeData )
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

        def calculateAverageEnergy( self, Ein ) :

            averageEp = self.energyForm.EpAverageAtE( Ein )
            if( self.productFrame == axes.centerOfMassToken ) :
                A_Ein = self.massRatio * Ein
                if( not( self.angularForm.isIsotropic( ) ) ) : averageEp += 2 * math.sqrt( A_Ein * averageEp ) * self.angularForm.muAverageAtEnergy( Ein )
                averageEp += A_Ein
            averageEp *= self.multiplicity.getValue( Ein )
            return( averageEp )

        def calculateAverageMomentum( self, Ein ) :

            pp, muAverageAtEnergy = 0., self.angularForm.muAverageAtEnergy( E )
            if( muAverageAtEnergy != 0. ) : pp = energyForm.sqrtEp_AverageAtE( E ) * muAverageAtEnergy
            if( self.productFrame == axes.centerOfMassToken ) : pp += math.sqrt( self.massRatio * Ein )
            return( multiplicity.getValue( E ) * math.sqrt( 2. * self.productMass ) * pp )

        def calculateAverageMomentumForPhoton( self, Ein ) :

            muAverage = angularForm.muAverageAtEnergy( E )
            if( muAverage == 0. ) : return( 0. )
            return( multiplicity.getValue( E ) * energyForm.EpAverageAtE( E ) * muAverage )

        class calculateDepositionInfo :

            def __init__( self, productFrame, productMass, massRatio, angularForm, energyForm, multiplicity ) :

                self.productFrame = productFrame
                self.productMass = productMass
                self.massRatio = massRatio
                self.angularForm = angularForm
                self.energyForm = energyForm
                self.multiplicity = multiplicity

            def evaluateAtX( self, E ) :

                return( self._evaluateAtX( self, E ) )

            def setEvaluateAtX( self, evaluateAtX ) :

                self._evaluateAtX = evaluateAtX

            def setTolerances( self, relativeTolerance, absoluteTolerance ) :

                self.relativeTolerance = relativeTolerance
                self.absoluteTolerance = absoluteTolerance

        class NBodyPhaseSpace :

            def __init__( self, energyForm, massUnit, projectileMass, targetMass, productMass, Q ) :

                self.energyForm = energyForm
                self.massUnit = massUnit
                self.projectileMass = projectileMass
                self.targetMass = targetMass
                self.productMass = productMass
                self.Q = Q

            def EpAverageAtE( self, Ein ) :

                return( self.energyForm.EpAverageAtE( Ein, self.massUnit, self.projectileMass, self.targetMass, self.productMass, self.Q ) )

            def getEnergyArray( self, EMin, EMax ) :

                return( [ EMin, EMax ] )

        depData = []
        energyAccuracy, momentumAccuracy = 1e-6, 1e-3
        energyUnit = tempInfo['incidentEnergyUnit']
        momentumDepositionUnit = energyUnit + '/c'
        massUnit = energyUnit + '/c**2'

        projectileMass = tempInfo['reactionSuite'].projectile.getMass( massUnit )
        targetMass = tempInfo['reactionSuite'].target.getMass( massUnit )
        productMass = tempInfo['product'].getMass( massUnit )
        massRatio = projectileMass * productMass / ( projectileMass + targetMass )**2

        angularForm = self.angularComponent.forms[self.angularComponent.nativeData]
        productFrame = angularForm.getProductFrame( )
        if( tempInfo['product'].getName( ) == 'gamma' ) : productFrame = axes.labToken  # All gamma data treated as in lab frame.
        energyForm = self.energyComponent.forms[self.energyComponent.nativeData]
        if( energyForm.moniker == base.NBodyPhaseSpaceFormToken ) :
            Q = tempInfo['outputChannel'].getQ( energyUnit, final = False, groundStateQ = True )
            energyForm = NBodyPhaseSpace( energyForm, massUnit, projectileMass, targetMass, productMass, Q )
        multiplicity = tempInfo['multiplicity']

        calculationData = calculateDepositionInfo( productFrame, productMass, massRatio, angularForm, energyForm, multiplicity )

        Es = energyForm.getEnergyArray( tempInfo['EMin'], tempInfo['EMax'] )
        multiplicityLimits = multiplicity.getEnergyLimits( tempInfo['EMin'], tempInfo['EMax'] )
        fixLimits( multiplicityLimits, Es )
        calculationData.setEvaluateAtX( calculateAverageEnergy )
        depEnergy = [ [ E, calculationData.evaluateAtX( E ) ] for E in Es ]
        absoluteTolerance = 1e-3 * energyAccuracy * max( [ Ep for E, Ep in depEnergy ] )
        calculationData.setTolerances( energyAccuracy, absoluteTolerance )
        depEnergy = fudgemath.thickenXYList( depEnergy, calculationData )
        axes_ = fudge.gnd.productData.energyDeposition.pointwise.defaultAxes( energyUnit = energyUnit, energyDepositionUnit = energyUnit )
        depData.append( fudge.gnd.productData.energyDeposition.pointwise( axes_, depEnergy, energyAccuracy ) )

        if( ( angularForm.moniker == base.isotropicFormToken ) and ( productFrame != axes.centerOfMassToken ) ) :
            depMomentum = [ [ depEnergy[0][0], 0. ], [ depEnergy[-1][0], 0. ] ]
        else :
            if( tempInfo['product'].getName( ) == 'gamma' ) :
                calculationData.setEvaluateAtX( calculateAverageMomentumForPhoton )
            else :
                calculationData.setEvaluateAtX( calculateAverageMomentum )
            for E in angularForm.getEnergyArray( tempInfo['EMin'], tempInfo['EMax'] ) :
                if( E not in Es ) : Es.append( E )
            Es.sort( )
            fixLimits( multiplicityLimits, Es )
            depMomentum = [ [ E, calculationData.evaluateAtX( E ) ] for E in Es ]
            absoluteTolerance = 1e-3 * momentumAccuracy * max( [ pp for E, pp in depMomentum ] )
            calculationData.setTolerances( momentumAccuracy, absoluteTolerance )
            depMomentum = fudgemath.thickenXYList( depMomentum, calculationData )
        axes_ = fudge.gnd.productData.momentumDeposition.pointwise.defaultAxes( energyUnit = energyUnit, momentumDepositionUnit = momentumDepositionUnit )
        depData.append( fudge.gnd.productData.momentumDeposition.pointwise( axes_, depMomentum, momentumAccuracy ) )

        return( depData )

    def checkProductFrame( self ) :
        """
        Checks that the productFrames for the angular and energy components are valid and the that their 'NativeData' forms 
        have the same frame. Returns None if all is okay. Otherwise, executes a raise from base.form.checkProductFrame
        or a ValueError if the angular and energy frames differ. Overrides the method in base.form.
        """

        self.angularComponent.checkProductFrame( )
        self.energyComponent.checkProductFrame( )
        if( self.angularComponent.getProductFrame( ) != self.energyComponent.getProductFrame( ) ) :
            raise ValueError( 'angular frame = "%s" != energy frame = "%s"' %
                ( self.angularComponent.getProductFrame( ) != self.energyComponent.getProductFrame( ) ) )

    def getProductFrame( self ) :
    
        return( self.angularComponent.getProductFrame( ) )

    def getSpectrumAtEnergy( self, energy ) :
        """Returns the energy spectrum for self at projectile energy using nativeData."""

        return( self.energyComponent.getSpectrumAtEnergy( energy ) )

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
                        if( not( angularForm_.axes.isLinear( flatIsOk = True ) ) ) : self.angularComponent.addForm( angularForm_.toPointwise_withLinearXYs( 1e-3 ) )
                        angularForm_ = None
                    if( angularForm.moniker == base.LegendrePointwiseFormToken ) :
                        self.angularComponent.addForm( angularForm_.toPointwise_withLinearXYs( 1e-3 ) )
                        angularForm_ = None
                    if( angularForm_ is not None ) : raise Exception( 'Not implemented: %s' % angularForm.moniker )

            energyFormKeys = self.energyComponent.forms.keys( )
            if( not( base.linearFormToken in energyFormKeys ) ) :
                energyForm_ = energyForm
                if( base.pointwiseFormToken in energyFormKeys ) :
                    if( not( energyForm_.axes.isLinear( flatIsOk = True, qualifierOk = True ) ) ) : self.energyComponent.addForm( energyForm_.toPointwise_withLinearXYs( 1e-3 ) )
                    energyForm_ = None
                if( energyForm_ is not None ) :
                    if( energyForm_.moniker in [ base.weightedFunctionalsFormToken, base.MadlandNixFormToken, base.NBodyPhaseSpaceFormToken ] ) :
                        energyForm_ = None
                    elif( energyForm_.moniker in [ base.semiPiecewiseFormToken, base.simpleMaxwellianFissionFormToken, \
                                base.evaporationFormToken, base.WattFormToken ] ) :
                        self.energyComponent.addForm( energyForm_.toPointwise_withLinearXYs( 1e-3, lowerEps = 1e-7, upperEps = 1e-7  ) )
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
                    totalMass = energyForm.numberOfProductsMasses.getValueAs( massUnit )
                    Q = tempInfo['outputChannel'].getQ( energyUnit, final = False, groundStateQ = True )
                    TM_1, TM_E = transferMatrices.NBodyPhaseSpace( processInfo, projectileName, productName, crossSection, \
                        energyForm.numberOfProducts, tempInfo['masses'], totalMass, Q, tempInfo['multiplicity'], \
                        comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )
                else :
                    from fudge.gnd.productData.distributions import angular
                    TM_1, TM_E = transferMatrices.uncorrelated_EMuP_EEpP_TransferMatrix( processInfo, projectileName, productName, tempInfo['masses'], \
                        crossSection, angularForm, energyForm, tempInfo['multiplicity'], \
                        comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )
            fudge.gnd.miscellaneous.TMs2Form( processInfo, tempInfo, newComponents, TM_1, TM_E, crossSection.axes )

        return( newComponents )

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        import angularEnergy

        angularForm = self.angularComponent.toPointwise_withLinearXYs( accuracy = accuracy, lowerEps = lowerEps, upperEps = upperEps )
        energyForm = self.energyComponent.toPointwise_withLinearXYs( accuracy = accuracy, lowerEps = lowerEps, upperEps = upperEps )
        grid = angularForm.getDomainGrid( )
        for e in energyForm.getDomainGrid( ) :
            if( e not in grid ) : grid.append( e )
        grid.sort( )
        axes_ = angularEnergy.pointwise.defaultAxes( energyUnit = energyForm.axes[0].getUnit( ), energy_outUnit = energyForm.axes[1].getUnit( ),
            probabilityUnit = energyForm.axes[2].getUnit( ) )
        f_E_mu_Ep = angularEnergy.pointwise( axes_, self.getProductFrame() )
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
                axes_ = Legendre.pointwise.defaultAxes( independent0, dependent0, energyInterpolationQualifier = qualifier0, 
                    energy_outInterpolation = independent1, C_lInterpolation = dependent1, energyUnit = axis0.getUnit( ), 
                    energy_outUnit = axis1.getUnit( ), C_lUnit = axis2.getUnit( ) )
                form2 = Legendre.pointwise( axes_, productFrame = form.getProductFrame( ) )          # change to LegendrePointwise.
                axes__ = axes.referenceAxes( form2 )

                for i, e_in in enumerate( form ) :
                    w_xys_LegendreSeries = LegendreSeries.W_XYs_LegendreSeries( axes__, index = i, value = e_in.value )
                    xys = e_in.copyDataToXYs( xUnit = 'eV', yUnit = '1/eV' )
                    for j, e_out_C_0 in enumerate( xys ) :
                        w_xys_LegendreSeries[j] = LegendreSeries.XYs_LegendreSeries( e_out_C_0[1:], index = j, value = e_out_C_0[0] )
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
                LAW, frame = 1, energyForm.productFrame
                gndToENDF6.toENDF6_MF6( MT, endfMFList, flags, targetInfo, LAW, frame, MF6 )
            else :
                raise Exception( 'uncorrelated.toENDF6 not supported for energy form = %s and angular form = %s' %
                    ( self.energyComponent.nativeData, self.angularComponent.nativeData ) )
        else :                          # original data is in uncorrelated form
            self.angularComponent.toENDF6( MT, endfMFList, flags, targetInfo )
            self.energyComponent.toENDF6( MT, endfMFList, flags, targetInfo )

def parseXMLNode( uncorrelatedElement, xPath=[], linkData = {} ):
    """Translate <uncorreleted> element from xml. Must contain one <angular> and one <energy> component."""

    xPath.append( uncorrelatedElement.tag )
    angular_ = angular.parseXMLNode( uncorrelatedElement.find(base.angularComponentToken), xPath )
    energy_ = energy.parseXMLNode( uncorrelatedElement.find(base.energyComponentToken), xPath )
    uncorrelated = component( angular_, energy_ )
    xPath.pop()
    return uncorrelated
