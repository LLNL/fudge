# <<BEGIN-copyright>>
# <<END-copyright>>

"""Angular distribution classes with components component, twoBodyComponent and CoulombElasticComponent."""

import math
import base, miscellaneous
import fudge
from fudge.core.utilities import brb
from fudge.core.math import fudgemath
from fudge.core.math.xData import axes, XYs, W_XYs, LegendreSeries, regions
from pqu import physicalQuantityWithUncertainty

__metaclass__ = type

def parseXMLNode( angularElement, linkData={} ):
    """ translate <angular> element from xml """
    angular = twoBodyComponent( nativeData = angularElement.get('nativeData') )
    for form in angularElement:
        formClass = {base.LegendrePointwiseFormToken: LegendrePointwise,
                base.LegendrePiecewiseFormToken: LegendrePiecewise,
                base.linearFormToken: linear,
                base.pointwiseFormToken: pointwise,
                base.isotropicFormToken: isotropic,
                base.recoilFormToken: recoil,
                base.mixedRangesFormToken: mixedRanges,
                }.get( form.tag )
        if formClass is None: raise Exception("encountered unknown angular form: %s" % form.tag)
        newForm = formClass.parseXMLNode( form, linkData )
        angular.addForm( newForm )
        newForm.parent = angular
    return angular

class component( base.component ) :

    def __init__( self, nativeData = base.noneFormToken ) :

        base.component.__init__( self, base.angularComponentToken, nativeData )

    def calculateDepositionData( self, processInfo, tempInfo ) :

        product = tempInfo['product']
        if( product.getName( ) != 'gamma' ) :
            raise Exception( 'For component %s, calculateDepositionData is only for gammas, not %s' % ( self.moniker, product.getName( ) ) )

        depData = []
        angularForm = self.forms[self.nativeData]
        energyUnit = tempInfo['incidentEnergyUnit']
        massUnit = energyUnit + '/c**2'
        momentumDepositionUnit = energyUnit + '/c'
        energyAccuracy, momentumAccuracy = 1e-6, 1e-3

        Es = angularForm.getEnergyArray( tempInfo['EMin'], tempInfo['EMax'] )
        projectileName, productName = processInfo.getProjectileName( ), product.getName( )
        if( ( 'discrete' in product.attributes ) or ( 'primary' in product.attributes ) ) :
            massRatio = 0.
            if( 'discrete' in product.attributes ) :
                Eg = product.attributes['discrete'].getValueAs( energyUnit )
            else :
                Eg = product.attributes['primary'].getValueAs( energyUnit )
                projectile, target = tempInfo['reactionSuite'].projectile, tempInfo['reactionSuite'].target
                mass1, mass2 = projectile.getMass( massUnit ), target.getMass( massUnit )
                massRatio = mass2 / ( mass1 + mass2 )
            depEnergy = [ [ E, Eg + massRatio * E ] for E in Es ]
            depMomentum = [ [ E, ( Eg + massRatio * E ) * angularForm.muAverageAtEnergy( E, accuracy = 0.1 * momentumAccuracy ) ] for E in Es ]
        else :
            raise Exception( 'Unsupported gamma; gamma must be "discrete" or "primary"' )

        axes_ = fudge.gnd.productData.energyDeposition.pointwise.defaultAxes( energyUnit = energyUnit, energyDepositionUnit = energyUnit )
        depData.append( fudge.gnd.productData.energyDeposition.pointwise( axes_, depEnergy, energyAccuracy ) )
        axes_ = fudge.gnd.productData.momentumDeposition.pointwise.defaultAxes( energyUnit = energyUnit, momentumDepositionUnit = momentumDepositionUnit )
        depData.append( fudge.gnd.productData.momentumDeposition.pointwise( axes_, depMomentum, momentumAccuracy ) )

        return( depData )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        from fudge.processing.deterministic import transferMatrices

        product = tempInfo['product']
        if( product.getName( ) != 'gamma' ) : raise Exception( 'For component %s, process is only for gammas, not %s' % product.getName( ) )

        newComponents = []
        angularForm = self.forms[self.nativeData]
        crossSection = tempInfo['crossSection']
        energyUnit = tempInfo['incidentEnergyUnit']
        massUnit = energyUnit + '/c**2'

        if( 'LLNL_MC' in processInfo['styles'] ) :
            raise Exception( 'Not implemented: %s' % self.moniker )

        if( 'LLNL_Pn' in processInfo['styles'] ) :
            if( processInfo['verbosity'] >= 30 ) : print '%sGrouping %s' % ( verbosityIndent, self.moniker )
            projectileName, productName = processInfo.getProjectileName( ), product.getName( )
            if( 'discrete' in product.attributes ) :
                Eg = product.attributes['discrete'].getValueAs( energyUnit )
                TM_1, TM_E = transferMatrices.discreteGammaAngularData( processInfo, projectileName, productName, Eg, crossSection,
                   angularForm, 1., comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )
            elif( 'primary' in product.attributes ) :
                projectile, target = tempInfo['reactionSuite'].projectile, tempInfo['reactionSuite'].target
                mass1, mass2 = projectile.getMass( massUnit ), target.getMass( massUnit )
                massRatio = mass2 / ( mass1 + mass2 )
                ELevel = product.attributes['primary'].getValueAs( energyUnit )
                TM_1, TM_E = transferMatrices.primaryGammaAngularData( processInfo, projectileName, productName, ELevel, massRatio, crossSection,
                   angularForm, 1., comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )
            else :
                raise Exception( 'Unsupported gamma; gamma must be "discrete" or "primary"' )
            fudge.gnd.miscellaneous.TMs2Form( processInfo, tempInfo, newComponents, TM_1, TM_E, crossSection.axes )

        return( newComponents )

    def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

        from fudge.legacy.converting import endfFormats, gndToENDF6
        if( hasattr( self.forms[self.nativeData], 'toENDF6' ) ) :
            if( self.nativeData == base.mixedRangesFormToken ) :
                LI, LTT, NM, frame, MF4 = self.forms[self.nativeData].toENDF6( flags, targetInfo )  # Note, LTT is returned as LAW = 2, 3 or 4
            else :
                NM = 0
                LI, LTT, frame, MF4 = self.forms[self.nativeData].toENDF6( flags, targetInfo )  # Note, LTT is returned as LAW = 2, 3 or 4
            if( targetInfo['doMF4AsMF6'] ) :
                if( LTT in [ 0 ] ) :
                    LAW = 3
                elif( LTT in [ 1, 2 ] ) :
                    LAW = 2
                elif( LTT in [ 4 ] ) :
                    LAW = 4
                else :
                    raise Exception( 'LTT = %s needs a LAW' % LTT )
                gndToENDF6.toENDF6_MF6( MT, endfMFList, flags, targetInfo, LAW, frame, MF4 )
            else :
                LCT = { axes.labToken : 1, axes.centerOfMassToken : 2 }[frame]
                endfMFList[4][MT] = [ endfFormats.endfHeadLine( targetInfo['ZA'], targetInfo['mass'],  0, LTT, 0, 0 ),
                                    endfFormats.endfHeadLine(                0, targetInfo['mass'], LI, LCT, 0, NM ) ] + MF4
        else :
            print 'WARNING: form %s does not have method toENDF6 for component %s' % ( self.nativeData, self.moniker )

class twoBodyComponent( component ) :

    def __init__( self, nativeData = base.noneFormToken ) :

        component.__init__( self, nativeData = nativeData )

    def calculateDepositionData( self, processInfo, tempInfo ) :

        def calculateDepositionEnergyAtE( angularData, E, parameters ) :

            a1x, a2y, Qp = parameters['a1x'], parameters['a2y'], parameters['Qp']
            dE = max( 0., E - Qp )
            return( a1x * E + a2y * dE + 2. * math.sqrt( a1x * a2y * E * dE ) * angularData.muAverageAtEnergy( E, accuracy = energyAccuracy ) )

        def calculateDepositionEnergyAtEForPhotoProjectile( angularData, E, parameters ) :

            mass2, massx, massy, Q = parameters['m2'], parameters['mx'], parameters['my'], parameters['Q']
            E__E_m2 = E / ( E + mass2 )
            dE = max( 0., ( E__E_m2 * mass2 ) + Q ) * massy / ( massx + massy )
            return( 0.5 * E__E_m2**2 * massx + dE + E__E_m2 * math.sqrt( 2. * dE * massx ) * angularData.muAverageAtEnergy( E, accuracy = energyAccuracy ) )

        class calculateDepositionEnergyThicken :

            def __init__( self, data, angular, func, parameters, relativeTolerance, absoluteTolerance ) :

                self.data = data
                self.angular = angular
                self.func = func
                self.parameters = parameters
                self.relativeTolerance = relativeTolerance
                self.absoluteTolerance = absoluteTolerance

            def evaluateAtX( self, E ) :

                return( self.func( self.angular, E, self.parameters ) )

        def calculateDepositionMomentumAtE( angularData, E, parameters ) :

            mass1, massx, b1x, a2y, Qp = parameters['m1'], parameters['mx'], parameters['b1x'], parameters['a2y'], parameters['Qp']
            dE = max( 0., E - Qp )
            return( b1x * math.sqrt( 2. * E / mass1 ) + math.sqrt( 2. * massx * a2y * dE ) * angularData.muAverageAtEnergy( E, accuracy = momentumAccuracy ) )

        def calculateDepositionMomentumAtEForPhotoProjectile( angularData, E, parameters ) :

            mass2, massx, massy, Q = parameters['m2'], parameters['mx'], parameters['my'], parameters['Q']
            E__E_m2 = E / ( E + mass2 )
            dE = max( 0., ( E__E_m2 * mass2 ) + Q ) * massy / ( massx + massy )
            return( E__E_m2 * massx + math.sqrt( 2. * dE * massx ) * angularData.muAverageAtEnergy( E, accuracy = momentumAccuracy ) )

        class calculateDepositionMomentumThicken :

            def __init__( self, data, angular, func, parameters, relativeTolerance, absoluteTolerance ) :

                self.data = data
                self.angular = angular
                self.func = func
                self.parameters = parameters
                self.relativeTolerance = relativeTolerance
                self.absoluteTolerance = absoluteTolerance

            def evaluateAtX( self, E ) :

                return( self.func( self.angular, E, self.parameters ) )

        angularForm = self.forms[self.nativeData]
        energyUnit = tempInfo['incidentEnergyUnit']
        massUnit = energyUnit + '/c**2'
        momentumDepositionUnit = energyUnit + '/c'

        depData = []
        projectile, target = tempInfo['reactionSuite'].projectile, tempInfo['reactionSuite'].target
        outputChannel = tempInfo['outputChannel'].outputChannel
        mass1 = projectile.getMass( massUnit )
        mass2 = target.getMass( massUnit )
        massx = outputChannel.particles[0].getMass( massUnit )
        massy = outputChannel.particles[1].getMass( massUnit )
        if( tempInfo['productIndex'] == '1' ) : massx, massy = massy, massx
        m12, mxy = mass1 + mass2, massx + massy
        b1x, a1x, a2y = mass1 * massx / m12, mass1 * massx / ( m12 * m12 ), mass2 * massy / ( m12 * mxy )
        Qm = m12 - mxy
        Q = tempInfo['outputChannel'].getQ( energyUnit, final = False, groundStateQ = True )
        Qp = -( Q ) * m12 / mass2    # This is the threshold in the COM frame.

        if( mass1 == 0. ) :         # Photo as projectile
            energyFunc = calculateDepositionEnergyAtEForPhotoProjectile
            momentumFunc = calculateDepositionMomentumAtEForPhotoProjectile
            parameters = { 'm2' : mass2, 'mx' : massx, 'my' : massy, 'Q' : Q }
        else :
            energyFunc = calculateDepositionEnergyAtE
            momentumFunc = calculateDepositionMomentumAtE
            parameters = { 'm1' : mass1, 'mx' : massx, 'a1x' : a1x, 'a2y' : a2y, 'b1x' : b1x, 'Qp' : Qp }

        energyAccuracy, momentumAccuracy = 1e-6, 1e-3
        Es = angularForm.getEnergyArray( tempInfo['EMin'], tempInfo['EMax'] )
        depEnergy = [ [ E, energyFunc( angularForm, E, parameters ) ] for E in Es ]
        depEnergy = fudgemath.thickenXYList( depEnergy, calculateDepositionEnergyThicken( depEnergy, angularForm, energyFunc, parameters, energyAccuracy, 1e-10 ) )
        axes_ = fudge.gnd.productData.energyDeposition.pointwise.defaultAxes( energyUnit = energyUnit, energyDepositionUnit = energyUnit )
        depData.append( fudge.gnd.productData.energyDeposition.pointwise( axes_, depEnergy, energyAccuracy ) )

        depMomentum = [ [ E, momentumFunc( angularForm, E, parameters ) ] for E in Es ]
        depMomentum = fudgemath.thickenXYList( depMomentum, calculateDepositionMomentumThicken( depMomentum, angularForm, momentumFunc, 
            parameters, momentumAccuracy, 1e-10 ) )
        axes_ = fudge.gnd.productData.momentumDeposition.pointwise.defaultAxes( energyUnit = energyUnit, momentumDepositionUnit = momentumDepositionUnit )
        depData.append( fudge.gnd.productData.momentumDeposition.pointwise( axes_, depMomentum, momentumAccuracy ) )

        return( depData )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        from fudge.processing.deterministic import transferMatrices

        energyUnit = tempInfo['incidentEnergyUnit']
        massUnit = energyUnit + '/c**2'
        newComponents = []
        crossSection = tempInfo['crossSection']
        angularForm = self.forms[self.nativeData]

        if( 'LLNL_MC' in processInfo['styles'] ) :
            formKeys = self.forms.keys( )
            if( ( base.linearFormToken not in formKeys ) and ( base.isotropicFormToken not in formKeys ) and ( base.recoilFormToken not in formKeys ) ) :
                form = self.getNativeData( )
                if( base.pointwiseFormToken in self.forms ) :
                    if( self.forms[base.pointwiseFormToken].axes.isLinear( qualifierOk = True ) ) : form = None
                if( form is not None ) : self.addForm( self.getNativeData( ).toPointwiseLinear( ) )

        if( 'LLNL_Pn' in processInfo['styles'] ) :
            outputChannel = tempInfo['outputChannel'].outputChannel
            if( processInfo['verbosity'] >= 30 ) : print '%sGrouping %s' % ( verbosityIndent, self.moniker )
            projectileName, productGroupName, productLabel = processInfo.getProjectileName( ), tempInfo['productGroupToken'], tempInfo['productLabel']
            projectile, target = tempInfo['reactionSuite'].projectile, tempInfo['reactionSuite'].target
            Q = tempInfo['outputChannel'].getQ( energyUnit, final = False, groundStateQ = True )
# The next lines do not work for non-zero excitation level as getMass does not includes it.
            mProj, eProj = fudge.gnd.miscellaneous.getMassLevel( projectile, massUnit, energyUnit )
            mTarget, eTarget = fudge.gnd.miscellaneous.getMassLevel( target, massUnit, energyUnit )
            mProd, eProd = fudge.gnd.miscellaneous.getMassLevel( outputChannel.particles[0], massUnit, energyUnit )
            mRes, eRes = fudge.gnd.miscellaneous.getMassLevel( outputChannel.particles[1], massUnit, energyUnit )
            if( tempInfo['productIndex'] == '1' ) : mProd, eProd, mRes, eRes = mRes, eRes, mProd, eProd
            TM_1, TM_E = transferMatrices.twoBodyTransferMatrix( processInfo, projectileName, productGroupName, crossSection,
                angularForm, mProj + eProj, mTarget + eTarget, mProd + eProd, mRes + eRes, Q,
                comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % productLabel )
            fudge.gnd.miscellaneous.TMs2Form( processInfo, tempInfo, newComponents, TM_1, TM_E, crossSection.axes )

        return( newComponents )

class CoulombElasticComponent( component ) :
    """Charged-particle coulomb scattering, similar to normal angular distribution with some additional information."""

    moniker = base.CoulombElasticComponentToken

    def __init__( self, nativeData = base.noneFormToken, identicalParticles = False, spin = 0 ) :

        component.__init__( self, nativeData = nativeData )
        self.identicalParticles = identicalParticles
        self.spin = spin
    
    def toXMLList( self, indent = "" ):

        xmlString = '%s<%s nativeData="%s"' % ( indent, self.moniker, self.nativeData )
        if( self.identicalParticles ) : xmlString += ' identicalParticles="true" spin="%s"' % self.spin
        xmlString += '>'
        xmlString = [ xmlString ] + self.forms[self.nativeData].toXMLList( indent=indent + '  ' )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

        from fudge.legacy.converting import endfFormats, gndToENDF6
        LI, LTT, frame, MF6 = self.forms[self.nativeData].toENDF6( flags, { 'doMF4AsMF6' : True, 'LIDP' : self.identicalParticles } )
        spin = 0            # first line of MF6 must change: need to know if we have identical particles, plus particle spin:
        if( self.identicalParticles ) : spin = self.spin.value
        NR, NE = [ int(a) for a in MF6[0].split()[-2:] ]
        MF6[0] = endfFormats.endfContLine( spin, 0, self.identicalParticles, 0, NR, NE )
        LAW = 5
        gndToENDF6.toENDF6_MF6( MT, endfMFList, flags, targetInfo, LAW, frame, MF6 )

    @staticmethod
    def parseXMLNode( element, linkData={} ):
        from fudge.gnd import xParticle
        identicalParticles = False
        if element.get('identicalParticles') in ('True','true'): identicalParticles = True
        spin = xParticle.spin( element.get('spin') )
        component = CoulombElasticComponent( nativeData = element.get('nativeData'),
                identicalParticles=identicalParticles, spin=spin )
        for form in element:
            formClass = {
                    base.nuclearPlusCoulombInterferenceFormToken: nuclearPlusCoulombInterference,
                    base.pointwiseFormToken: pointwise,
                    }.get( form.tag )
            if formClass is None: raise Exception("encountered unknown CoulombElastic form: %s" % form.tag)
            newForm = formClass.parseXMLNode( form, linkData )
            component.addForm( newForm )
            newForm.parent = component
        return component
    
class isotropic( base.form ) :

    def __init__( self, frame ) :

        base.form.__init__( self, base.isotropicFormToken )
        self.frame = frame

    def domainMin( self, unitTo = None, asPQU = False ) :

        from fudge.gnd import product
        return( self.findClassInAncestry( product.product ).domainMin( unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        from fudge.gnd import product
        return( self.findClassInAncestry( product.product ).domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getDomain( self, unitTo = None, asPQU = False ) :

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ), self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        return( [ EMin, EMax ] )

    def invert( self ) :

        return( self.toPointwiseLinear( ) )

    def isIsotropic( self ) :

        return( True )

    def muAverageAtEnergy( self, E, accuracy = None ) :

        return( 0. )

    def check( self, info ) :

        return []

    def toPointwiseLinear( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        axes_ = pointwise.defaultAxes( )
        p = linear( axes_ )
        axes_muP = axes.referenceAxes( p )
        p.append( XYs.XYs( axes_muP, [ [ -1, 0.5 ], [ 1, 0.5 ] ], value = self.domainMin( ), accuracy = 1e-6 ) )
        p.append( XYs.XYs( axes_muP, [ [ -1, 0.5 ], [ 1, 0.5 ] ], value = self.domainMax( ), accuracy = 1e-6 ) )
        return( p )

    def toXMLList( self, indent = ""  ) :

        return( [ '%s<%s frame="%s"/>' % ( indent, self.moniker, self.frame ) ] )

    def toENDF6( self, flags, targetInfo ) :

        from fudge.legacy.converting import endfFormats, gndToENDF6
        ENDFDataList = []
        if( not( targetInfo['doMF4AsMF6'] ) ) : ENDFDataList.append( endfFormats.endfSENDLineNumber( ) )
        return( 1, 0, self.frame, ENDFDataList )

    @staticmethod
    def parseXMLNode( element, linkData={} ): return isotropic( element.get('frame') )

class recoil( base.form ) :

    def __init__( self, product ) :

        base.form.__init__( self, base.recoilFormToken )
        self.setReference( product )

    def getAxes( self ):
        form = self.getPartnersForm( )
        return( form.axes )

    axes = property( getAxes )

    def getNumericalDistribution( self ) :

        form = self.getPartnersForm( ).invert( )
        return( form )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        return( self.getPartnersForm( ).getEnergyArray( EMin = EMin, EMax = EMax ) )

    def getPartnersForm( self ) :

        component = self.product.distributions.components[self.product.distributions.nativeData]
        return( component.forms[component.nativeData] )

    def setReference( self, product ) :

        self.product = product

    def isIsotropic( self ) :

        return( self.getPartnersForm( ).isIsotropic( ) )

    def muAverageAtEnergy( self, E, accuracy = None ) :

        return( -self.getPartnersForm( ).muAverageAtEnergy( E, accuracy = accuracy ) )

    def check( self, info ) :

        from fudge.gnd import warning

        warnings = []
        if self.getPartnersForm().moniker not in ( base.isotropicFormToken, base.pointwiseFormToken,
                base.LegendrePointwiseFormToken, base.LegendrePiecewiseFormToken ) :
            warnings.append( warning.missingRecoilDistribution( self ) )

        return warnings

    def toPointwiseLinear( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        return( self.getPartnersForm( ).invert( ).toPointwiseLinear( accuracy = accuracy, lowerEps = lowerEps, upperEps = upperEps ) )

    def toXMLList( self, indent = '' ) :
        """Returns the xml string representation of self."""

        return( [ '%s<%s xlink:type="simple" xlink:href="%s"/>' % ( indent, self.moniker, self.product.toXLink( ) ) ] )

    def toENDF6( self, flags, targetInfo ) :

        return( 1, 4, axes.centerOfMassToken, [] )

    @staticmethod
    def parseXMLNode( recoilElement, linkData={} ):
        """ translate <recoil> element from xml """
        from fudge.gnd import link
        xlink = link.parseXMLNode( recoilElement ).path
        ref = recoil( None )
        if 'unresolvedLinks' in linkData: linkData['unresolvedLinks'].append((xlink, ref))
        return ref

class pointwise( base.form, W_XYs.W_XYs ) :

    tag = base.pointwiseFormToken
    moniker = base.pointwiseFormToken

    def __init__( self, axes, **kwargs ) :

        base.form.__init__( self, self.tag )
        kwargs['isPrimaryXData'] = True
        W_XYs.W_XYs.__init__( self, axes, **kwargs )

    def getAtEnergy( self, energy ) :

        return( self.interpolateAtW( energy, unitBase = True, extrapolation = W_XYs.flatExtrapolationToken ) )

    def invert( self ) :

        c = pointwise( self.axes )
        axesMuP = axes.referenceAxes( c )
        for index, xys in enumerate( c ) :
            n_xys = xys.copyDataToXYs( ) 
            reverse( n_xys )
            for xy in n_xys  : xy[0] = -xy[0]
            c[index] = XYs.XYs( axesMuP, n_xys, accuracy = xys.getAccuracy( ), value = xys.value, index = index, parent = c )
        return( c )

    def isIsotropic( self ) :

        for energy_in in self :
            if( energy_in.yMin( ) != energy_in.yMax( ) ) : return( False )
        return( True )

    def muAverageAtEnergy( self, E, accuracy = None ) :

        muP = self.getAtEnergy( E )
        if( accuracy is None ) : accuracy = muP.getAccuracy( )
        muP = muP.changeInterpolationIfNeeded( [ [ axes.linearToken, axes.linearToken ], [ axes.linearToken, axes.flatToken ] ], accuracy = accuracy )
        n = muP.normalize( )
        return( muP.normalize( ).integrateWithWeight_x( ) )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        Es = [ data.value for data in self ]
        if( EMin is not None ) :
            if( EMin < Es[0] ) : Es.insert( 0, EMin )
        if( EMax is not None ) :
            if( EMax > Es[-1] ) : Es.append( EMax )
        return( Es )

    def check( self, info ) :

        from fudge.gnd import warning
        PQU = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty

        warnings = []
        for i in range(len(self)):
            integral = self[i].integrate(-1,1)
            if abs(integral - 1.0) > info['normTolerance']:
                warnings.append( warning.unnormalizedDistribution( PQU(self[i].value,self.axes[0].unit), i,
                    integral, self[i] ) )

            if self[i].yMin() < 0.0:
                warnings.append( warning.negativeProbability( PQU(self[i].value,self.axes[0].unit),
                    value = self[i].yMin(), obj=self ) )

        return warnings

    def toPointwiseLinear( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        return( W_XYs.W_XYs.toPointwiseLinear( self, accuracy = accuracy, lowerEps = lowerEps, upperEps = upperEps, cls = linear ) )

    def toENDF6( self, flags, targetInfo ) :

        from fudge.legacy.converting import endfFormats, gndToENDF6
        ENDFDataList = [ endfFormats.endfContLine( 0, 0, 0, 0, 1, len( self ) ) ]
        interpolation = gndToENDF6.axisToEndfInterpolationFlag( self.axes[0] )
        ENDFDataList += endfFormats.endfInterpolationList( [ len( self ), interpolation ] )
        for energy_in in self : ENDFDataList += gndToENDF6.angularPointwiseEnergy2ENDF6( energy_in, targetInfo )
        if( not( targetInfo['doMF4AsMF6'] ) ) : ENDFDataList.append( endfFormats.endfSENDLineNumber( ) )
        return( 0, 2, self.axes[1].frame, ENDFDataList )

    @staticmethod
    def defaultAxes( energyUnit = 'eV', energyInterpolation = axes.linearToken, energyFunctionInterpolation = axes.linearToken, energyInterpolationQualifier = None,
            energyFrame = axes.labToken, muInterpolation = axes.linearToken, muProbabilityFrame = axes.centerOfMassToken, probabilityInterpolation = axes.linearToken, probabilityUnit = '' ) :

        axes_ = axes.axes( dimension = 3 )
        axes_[0] = axes.axis( 'energy_in', 0, energyUnit, frame = axes.labToken, 
            interpolation = axes.interpolationXY( energyInterpolation, energyFunctionInterpolation, energyInterpolationQualifier ) )
        axes_[1] = axes.axis( 'mu', 1, '', frame = muProbabilityFrame, interpolation = axes.interpolationXY( muInterpolation, probabilityInterpolation ) )
        axes_[2] = axes.axis( 'P(mu|energy_in)', 2, probabilityUnit, frame = muProbabilityFrame )
        return( axes_ )

class linear( pointwise ) :

    tag = base.linearFormToken
    moniker = base.linearFormToken

    def __init__( self, axes, **kwargs ) :

        if( not axes.isLinear( qualifierOk = True ) ) : raise Exception( 'interpolation not linear: %s' % axes )
        pointwise.__init__( self, axes, **kwargs )

class piecewise( base.distributionFormWithRegions ) :
    """This class has never been tested and should not be used. Was in original ENDF to Gnd and that is why it is still here.
    Not a very good reason. Can be deleted if no database needs it, but this needs to be checked first."""

    def __init__( self, axes, IKnowWhatIAmDoing = False ) :

        if( not( IKnowWhatIAmDoing ) ) : raise Exception( "This class should not be used" )
        base.distributionFormWithRegions.__init__( self, base.piecewiseFormToken )
        self.axes = axes

    def isIsotropic( self ) :

        return( False )

    def toXMLList( self, indent = ""  ) :

        indent2 = indent + '  '
        xmlString = [ '%s<%s>' % ( indent, self.moniker ) ]
        xmlString += self.axes.toXMLList( indent = indent2 )
        xmlString += base.distributionFormWithRegions.toXMLList( self, indent2 )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    def toENDF6( self, flags, targetInfo, inMixedForm = False ) :

        from fudge.legacy.converting import endfFormats, gndToENDF6
        counter, lineData, interpolationFlatData, lastEnergyData, ENDFDataList = 0, [], [], None, []
        for region in self :
            nE, ENDFInterpolation, currentLineData, lastEnergyData = region.toENDF6( flags, targetInfo, lastEnergyData )
            counter += nE
            lineData += currentLineData
            interpolationFlatData.append( counter )
            interpolationFlatData.append( ENDFInterpolation )
        ENDFDataList += [ endfFormats.endfContLine( 0, 0, 0, 0, len( self ), counter ) ]
        ENDFDataList += endfFormats.endfInterpolationList( interpolationFlatData )
        ENDFDataList += lineData
        if( not( targetInfo['doMF4AsMF6'] ) ) : ENDFDataList.append( endfFormats.endfSENDLineNumber( ) )
        return( 0, 2, self.axes[1].frame, ENDFDataList )

class piecewiseRegion :

    def __init__( self, index, interpolation, yzInterpolation ) :

        self.index = index
        self.interpolationString = interpolation
        self.yzInterpolationString = yzInterpolation
        self.energies = []

    def __len__( self ) :

        return( len( self.energies ) )

    def append( self, energy, pdf_of_mu ) :

        if( len( self.energies ) > 0 ) :
            if( self.energies[-1].value >= energy ) : 
                raise Exception( 'Energy = %e must be greater than last energy = %e' % ( energy, self.energies[-1].value ) )
        self.energies.append( piecewiseRegionEnergy( len( self.energies ), energy, pdf_of_mu ) )

    def domain( self ) :

        return( self.energies[0].value, self.energies[-1].value )

    def toXMLList( self, indent = ""  ) :

        xmlString = [ '%s<region index="%d" type="3d.xyz" interpolation="%s" yzInterpolation="%s">' % 
            ( indent, self.index, self.interpolationString, self.yzInterpolationString ) ]
        for index, energy in enumerate( self.energies ) : xmlString += energy.toXMLList( indent + '  ', index = index )
        xmlString[-1] += '</region>'
        return( xmlString )

    def toENDF6( self, flags, targetInfo, lastEnergysFlatData ) :

        from fudge.legacy.converting import endfFormats, gndToENDF6
        counter, ENDFInterpolation, lineData = 0, endfFormats.XYStringToENDFInterpolation( self.interpolationString ), []
        yzInterpolation = endfFormats.XYStringToENDFInterpolation( self.yzInterpolationString )
        if( targetInfo['doMF4AsMF6'] ) : yzInterpolation += 10
        for energy in self.energies :
            counter += 1
            currentFlatData = energy.toENDF6( flags, targetInfo, interpolation = yzInterpolation )
            if( energy.index == 0 ) :
                if( not( lastEnergysFlatData is None ) ) :
                    if( lastEnergysFlatData == currentFlatData ) :
                        currentFlatData = []
                        counter  -= 1
            lineData += currentFlatData
        return( counter, ENDFInterpolation, lineData, currentFlatData )

class piecewiseRegionEnergy :

    moniker = 'energy_in'

    def __init__( self, index, value, xys ) :

        self.index = index
        self.value = value
        self.xys = xys 
    
    def __len__( self ) :
    
        return( len( self.xys ) )
    
    def __getitem__( self, i ) :
    
        return( self.xys[i] )

    def getValue( self ) :

        return( self.value )
    
    def toXMLList( self, indent = "", floatFormat = "%s", index = None ) :

        if( index is None ) :
            indexStr = ''
        else :
            indexStr = ' index="%s"' % index
        tag = '%s<%s value="%s"%s>' % ( indent, self.moniker, self.value, indexStr )
        return( [ tag + endl2dmathmisc.list2dToXMLPointwiseString( self.xys, indent = '', floatFormat = floatFormat ) + '</%s>' % self.moniker ] )

    def toENDF6( self, flags, targetInfo, interpolation = 2 ) :

        if( targetInfo['doMF4AsMF6'] ) :
            ENDFDataList = [ endfFormats.endfContLine( 0, self.value, interpolation, 0, 2 * len( self.xys ), len( self.xys ) ) ]
            ENDFDataList += endfFormats.endfNdDataList( self.xys )
        else :
            ENDFDataList = [ endfFormats.endfContLine( 0, self.value, 0, 0, 1, len( self.xys ) ) ]
            ENDFDataList += endfFormats.endfInterpolationList( [ len( self.xys ), interpolation ] )
            ENDFDataList += endfFormats.endfNdDataList( self.xys )
        return( ENDFDataList )

class LegendrePointwise( base.form, LegendreSeries.W_XYs_LegendreSeries ) :

    tag = base.LegendrePointwiseFormToken
    moniker = base.LegendrePointwiseFormToken

    def __init__( self, axes ) :

        base.form.__init__( self, base.LegendrePointwiseFormToken )
        LegendreSeries.W_XYs_LegendreSeries.__init__( self, axes, isPrimaryXData = True )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        Es = [ data.value for data in self ]
        if( EMin is not None ) :
            if( EMin < Es[0] ) : Es.insert( 0, EMin )
        if( EMax is not None ) :
            if( EMax > Es[-1] ) : Es.append( EMax )
        return( Es )

    def getAtEnergy( self, energy ) :

        return( self.getLegendreSeriesAtW( energy, extrapolation = W_XYs.flatExtrapolationToken ) )

    def invert( self ) :

        c = LegendrePointwise( self.axes )
        unit = self.axes[1].getUnit( )
        for index, coeffs in enumerate( self ) : c[index] = LegendreSeries.XYs_LegendreSeries( unit, coeffs.invert( ), value = coeffs.value, index = index )
        return( c )

    def muAverageAtEnergy( self, E, accuracy = None ) :

        LegendreSeries_ = self.getAtEnergy( E )
        if( len( LegendreSeries_ ) < 2 ) : return( 0. )
        return( LegendreSeries_[1] )

    def check( self, info ) :

        from fudge.gnd import warning
        PQU = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty

        warnings = []
        for Ein in self.LegendreSeries_s:
            if Ein.coefficients[0] != 1.0:
                warnings.append( warning.unnormalizedDistribution( PQU(Ein.value,self.axes[0].unit),
                    Ein.index, Ein.coefficients[0], Ein ) )
            distribution = Ein.toPointwiseLinear( 0.0001 )
            if distribution.yMin() < 0:
                warnings.append( warning.negativeProbability( PQU(Ein.value,self.axes[0].unit),
                    value = distribution.yMin(), obj = self ) )

        return warnings

    def toPointwiseLinear( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        if( accuracy is None ) : accuracy = 1e-3
        w_xys = LegendreSeries.W_XYs_LegendreSeries.toPointwiseLinear( self, accuracy )
        axes_ = pointwise.defaultAxes( )
        p = linear( axes_ )
        for xys in w_xys : p.append( xys )
        return( p )

    def toENDF6( self, flags, targetInfo ) :

        from fudge.legacy.converting import endfFormats, gndToENDF6
        ENDFDataList = [ endfFormats.endfContLine( 0, 0, 0, 0, 1, len( self ) ) ]
        ENDFDataList += endfFormats.endfInterpolationList( [ len( self ), 2 ] )
        for energy in self : 
            NW, NL = len( energy ) - 1, 0
            if( targetInfo['doMF4AsMF6'] ) : NL = NW
            ENDFDataList.append( endfFormats.endfContLine( 0, energy.value, 0, 0, NW, NL ) )
            ENDFDataList += endfFormats.endfDataList( energy.coefficients[1:] )
        if( not( targetInfo['doMF4AsMF6'] ) ) : ENDFDataList.append( endfFormats.endfSENDLineNumber( ) )
        return( 0, 1, self.axes[1].frame, ENDFDataList )

    @staticmethod
    def defaultAxes( xyFrame, energyInterpolation, dependentInterpolation, energyLabel = 'energy_in', energyUnit = 'eV', interpolationQualifier = None ) :

        axes_ = axes.axes( dimension = 2 )
        axes_[0] = axes.axis( energyLabel,    0, energyUnit, frame = axes.labToken, \
            interpolation = axes.interpolationXY( energyInterpolation, dependentInterpolation, interpolationQualifier ) )
        axes_[1] = axes.axis( 'C_l', 1, '', frame = xyFrame )
        return( axes_ )

    @staticmethod
    def parseXMLNode( LegendrePointwiseElement, linkData={} ):
        """ translate a <LegendrePointwise> element from xml """
        axes_ = axes.parseXMLNode( LegendrePointwiseElement[0] )
        pointwiseForm = LegendrePointwise( axes_ )
        for energy_in in LegendrePointwiseElement[1:]:
            value = float(energy_in.get("value"))
            coefs = map(float, energy_in.text.split())
            assert len(coefs)==int(energy_in.get("length"))
            pointwiseForm.append( LegendreSeries.XYs_LegendreSeries( axes_[-1].getUnit(), coefs, value=value ) )
        return pointwiseForm

class LegendrePiecewise( base.form, regions.regions ) :

    tag = base.LegendrePiecewiseFormToken
    xData = "regions" + LegendreSeries.W_XYs_LegendreSeries.xData

    def __init__( self, axes ) :

        base.form.__init__( self, base.LegendrePiecewiseFormToken )
        regions.regions.__init__( self, LegendreSeries.W_XYs_LegendreSeries, axes )

    def domainMin( self, unitTo = None, asPQU = False ) :

        return( self[0][0].value )

    def domainMax( self, unitTo = None, asPQU = False ) :

        return( self[-1][-1].value )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        Es = []
        for region in self : 
            for data in region :
                if( data.value not in Es ) : Es.append( data.value )
        if( EMin is not None ) :
            if( EMin < Es[0] ) : Es.insert( 0, EMin )
        if( EMax is not None ) :
            if( EMax > Es[-1] ) : Es.append( EMax )
        return( Es )

    def isIsotropic( self ) :

        for region in self :
            if( not( region.isIsotropic( ) ) ) : return( False )
        return( True )

    def maxLegendreOrder( self ) :

        NM = 0
        for region in self : NM = max( NM, region.maxLegendreOrder( ) )
        return( NM )

    def muAverageAtEnergy( self, E, accuracy = None ) :

        for i, region in enumerate( self ) :
            if( region.value > E ) : break
        if( i > 0 ) : i -= 1
        LegendreSeries_ = self[i].getLegendreSeriesAtW( E, extrapolation = W_XYs.flatExtrapolationToken )
        if( len( LegendreSeries_ ) < 2 ) : return( 0. )
        return( LegendreSeries_[1] )

    def check( self, info ) :
        from fudge.gnd import warning
        PQU = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty

        warnings = []
        for region in self:
            subwarnings = []
            for Ein in region.LegendreSeries_s:
                if Ein.coefficients[0] != 1.0:
                    subwarnings.append( warning.unnormalizedDistribution( PQU(Ein.value,self.axes[0].unit),
                        Ein.index, Ein.coefficients[0], Ein ) )

                distribution = Ein.toPointwiseLinear( 0.0001 )
                if distribution.yMin() < 0:
                    subwarnings.append( warning.negativeProbability( PQU(Ein.value,self.axes[0].unit),
                        value = distribution.yMin(), obj = self ) )
            if subwarnings: warnings.append( warning.context("Region %i:" % self.regions.index(region),
                subwarnings) )

        return warnings

    def toENDF6( self, flags, targetInfo, inMixedForm = False ) :

        from fudge.legacy.converting import endfFormats, gndToENDF6
        counter, lineData, interpolationFlatData, lastEnergyData, ENDFDataList = 0, [], [], None, []
        for region in self :
            nE, ENDFInterpolation, currentLineData, lastEnergyData = W_XYs_LegendreSeries_toENDF6( region, flags, targetInfo, lastEnergyData )
            counter += nE
            lineData += currentLineData
            interpolationFlatData.append( counter )
            interpolationFlatData.append( ENDFInterpolation )
        ENDFDataList += [ endfFormats.endfContLine( 0, 0, 0, 0, len( self ), counter ) ]
        ENDFDataList += endfFormats.endfInterpolationList( interpolationFlatData )
        ENDFDataList += lineData
        if( not( inMixedForm ) ) :
            if( not( targetInfo['doMF4AsMF6'] ) ) : ENDFDataList.append( endfFormats.endfSENDLineNumber( ) )
        return( 0, 1, self.axes[1].frame, ENDFDataList )

    def toPointwiseLinear( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        from fudge.core.utilities import mathMisc

        if( accuracy is None ) : accuracy = 1e-3
        axes_ = pointwise.defaultAxes( )
        pwl = linear( axes_ )
        w_xys = self[0].toPointwiseLinear( accuracy )
        for xys in w_xys : pwl.append( xys )
        for i, region in enumerate( self ) :
            if( i == 0 ) : continue
            regionForm = region.toPointwiseLinear( accuracy )
            for j, xys in enumerate( regionForm ) :
                if( j == 0 ) : xys.value = mathMisc.shiftFloatUpABit( xys.value, 1e-6 )
                pwl.append( xys )
        return( pwl )

    def toXMLList( self, indent = ""  ) :

        indent2 = indent + '  '
        xmlString = [ '%s<%s xData="%s">' % ( indent, self.tag, self.xData ) ]
        xmlString += self.axes.toXMLList( indent = indent2 )
        for region in self : xmlString += region.toXMLList( tag = 'region', indent = indent2 )
        xmlString[-1] += '</%s>' % self.tag
        return( xmlString )

    @staticmethod
    def defaultAxes( xyFrame, energyInterpolation = axes.byRegionToken, dependentInterpolation = axes.byRegionToken, energyLabel = 'energy_in', 
            energyUnit = 'eV', interpolationQualifier = None ) :

        axes_ = axes.axes( dimension = 2 )
        axes_[0] = axes.axis( energyLabel,    0, energyUnit, frame = axes.labToken, \
            interpolation = axes.interpolationXY( energyInterpolation, dependentInterpolation, interpolationQualifier ) )
        axes_[1] = axes.axis( 'C_l', 1, '', frame = xyFrame )
        return( axes_ )

    @staticmethod
    def parseXMLNode( LPiecewiseElement, linkData={} ):
        axes_ = axes.parseXMLNode( LPiecewiseElement[0] )
        LPW = LegendrePiecewise( axes_ )
        for reg in LPiecewiseElement[1:]:
            axes_ = axes.parseXMLNode( reg[0] )
            region_ = LegendreSeries.W_XYs_LegendreSeries( axes_, index=int(reg.get("index")), parent=LPW,
                    isPrimaryXData = False )
            for energy_in in reg[1:]:
                region_.append( LegendreSeries.XYs_LegendreSeries('', map(float, energy_in.text.split()),
                    index=int(energy_in.get("index")), value=float(energy_in.get("value")), parent=region_ ) )
            LPW.regions.append( region_ )
        return LPW

class nuclearPlusCoulombInterference( base.form ) :
    """For charged-particle elastic scattering, data may be in simple legendre or mu/P(mu) format,
    but it may also be given in terms of legendre coefficients for nuclear and interference terms."""

    def __init__( self, nuclear_term, interferenceReal_term, interferenceImaginary_term ) :

        if( isinstance( nuclear_term, LegendrePointwise ) ) :
            if( not( isinstance( interferenceReal_term, LegendrePointwise ) ) ) :  
                raise Exception( 'nuclear_term is type LegendrePointwise and so must interferenceReal_term  = %s' % brb.getType( interferenceReal_term ) )
            if( not( isinstance( interferenceImaginary_term, LegendrePointwise ) ) ) :
                raise Exception( 'nuclear_term is type LegendrePointwise and so must interferenceImaginary_term = %s' % brb.getType( interferenceReal_term ) )
        elif( isinstance( nuclear_term, LegendrePiecewise ) ) :
            if( not( isinstance( interferenceReal_term, LegendrePiecewise ) ) ) : 
                raise Exception( 'nuclear_term is type LegendrePiecewise and so must interferenceReal_term  = %s' % brb.getType( interferenceReal_term ) )
            if( not( isinstance( interferenceImaginary_term, LegendrePiecewise ) ) ) :
                raise Exception( 'nuclear_term is type LegendrePiecewise and so must interferenceImaginary_term = %s' % brb.getType( interferenceReal_term ) )
        else :
            raise Exception( 'Only LegendrePointwise and LegendrePiecewise are currently supported, not %s' % brb.getType( nuclear_term ) )
        base.form.__init__( self, base.nuclearPlusCoulombInterferenceFormToken )
        self.nuclear_term = nuclear_term
        self.interferenceReal_term = interferenceReal_term
        self.interferenceImaginary_term = interferenceImaginary_term

    def isIsotropic( self ) :

        return( False )

    def toXMLList( self, indent = "" ) :

        xmlString = [ '%s<%s>' % ( indent, self.moniker ) ]
        for term in ('nuclear_term','interferenceReal_term','interferenceImaginary_term'):
            xmlString += ['%s<%s>' % (indent+'  ', term)]
            xmlString += getattr(self, term).toXMLList( indent=indent+'    ' )
            xmlString[-1] += '</%s>' % term
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    def toENDF6( self, flags, targetInfo, inMixedForm = False ) :

        def LTP_oneSubParsing( LTP, nuclear, interferenceReal, interferenceImaginary, lineData ) :

            if( targetInfo['LIDP'] ):
                NL = len( nuclear ) - 1
                NW = 3 * NL + 3
            else :
                NL = ( len( nuclear ) - 1 ) // 2
                NW = 4 * NL + 3
            lineData.append( endfFormats.endfContLine( 0, nuclear.value, LTP, 0, NW, NL ) )
            legendreDat = nuclear.coefficients
            for j, r in enumerate( interferenceReal.coefficients ) :
                legendreDat.append( r )
                legendreDat.append( interferenceImaginary[j] )
            lineData += endfFormats.endfDataList( legendreDat )

        from fudge.legacy.converting import endfFormats, gndToENDF6
        counts, interpolationFlagsList, lineData = 0, [], []
        LTP = 1                     # indicates this is a nuclear + interference section
        if( isinstance( self.nuclear_term, LegendrePointwise ) ) :
            for ridx in xrange( len( self.nuclear_term ) ) :
                counts += 1
                nuclear, interferenceReal, interferenceImaginary = self.nuclear_term[ridx], self.interferenceReal_term[ridx], self.interferenceImaginary_term[ridx]
                LTP_oneSubParsing( LTP, nuclear, interferenceReal, interferenceImaginary, lineData )
            interpolationFlagsList += [ counts, gndToENDF6.axisToEndfInterpolationFlag( self.nuclear_term.axes[0] ) ]
        elif( isinstance( self.nuclear_term, LegendrePiecewise ) ) :
            for regionIndex, region in enumerate( self.nuclear_term ) :
                interferenceReal, interferenceImaginary = self.interferenceReal_term[regionIndex], self.interferenceImaginary_term[regionIndex]
                for energyIndex, nuclear in enumerate( region ) :
                    if( ( regionIndex != 0 ) and ( energyIndex == 0 ) ) : continue
                    counts += 1
                    LTP_oneSubParsing( LTP, nuclear, interferenceReal[energyIndex], interferenceImaginary[energyIndex], lineData )
                interpolationFlagsList += [ counts, gndToENDF6.axisToEndfInterpolationFlag( region.axes[0] ) ]
        else :
            raise Exception( 'Unsupported angular data type = %s' % brb.getType( self.nuclear_term ) )
        interpolationFlags = endfFormats.endfInterpolationList( interpolationFlagsList )
        ENDFDataList = [ endfFormats.endfContLine( 0, 0, 0, 0, len( interpolationFlagsList ) / 2, counts ) ] + interpolationFlags + lineData
        if( not( targetInfo['doMF4AsMF6'] ) ) : ENDFDataList.append( endfFormats.endfSENDLineNumber( ) )
        return( 0, 2, self.nuclear_term.axes[1].frame, ENDFDataList )

    @staticmethod
    def parseXMLNode( element, linkData={} ):
        # nuclearPlusCoulombInterference is divided into three sections: nuclear, real and imag. interference
        # each may be LegendrePointwise or LegendrePiecewise
        def getPointwiseOrPiecewise( element ):
            form = { base.LegendrePointwiseFormToken: LegendrePointwise,
                    base.LegendrePiecewiseFormToken: LegendrePiecewise, }[ element[0].tag ]
            return form.parseXMLNode( element[0], linkData )
        nuclear, realInterference, imagInterference = [getPointwiseOrPiecewise( term )
                for term in element]
        return nuclearPlusCoulombInterference( nuclear, realInterference, imagInterference )
    
class mixedRanges( base.form ) :

    def __init__( self, LegendreForm, tabulatedForm ) :

        if( not( isinstance( LegendreForm, ( LegendrePiecewise, LegendrePointwise ) ) ) ) :
            raise Exception( 'Invalid Legendre instance: %s' % brb.getType( LegendreForm ) )
        if( not( isinstance( tabulatedForm, ( linear, pointwise ) ) ) ) :
            raise Exception( 'Invalid tabulated form instance: %s' % brb.getType( LegendreForm ) )
        base.form.__init__( self, base.mixedRangesFormToken )
        self.LegendreForm = LegendreForm
        self.tabulatedForm = tabulatedForm

    def domainMin( self, unitTo = None, asPQU = False ) :

        return( self.LegendreForm.domainMin( unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        return( self.tabulatedForm.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getDomain( self, unitTo = None, asPQU = False ) :

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ), self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        Es = self.LegendreForm.getEnergyArray( EMin = EMin, EMax = EMax )
        for E in self.tabulatedForm.getEnergyArray( EMin = EMin, EMax = EMax ) :
            if( E not in Es ) : Es.append( E )
        Es.sort( )
        return( Es )

    def isIsotropic( self ) :

        return( self.LegendreForm.isIsotropic( ) and self.tabulatedForm.isIsotropic( ) )

    def muAverageAtEnergy( self, E, accuracy = None ) :

        if( E < self.LegendreForm.domainMax( ) ) : return( self.LegendreForm.muAverageAtEnergy( E, accuracy = accuracy ) )
        return( self.tabulatedForm.muAverageAtEnergy( E, accuracy = accuracy ) )

    def check( self, info ) :

        warnings = []
        warnings += self.LegendreForm.check( info )
        warnings += self.tabulatedForm.check( info )
        return warnings

    def toPointwiseLinear( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        from fudge.core.utilities import mathMisc

        LegendreForm = self.LegendreForm.toPointwiseLinear( accuracy, lowerEps = lowerEps, upperEps = upperEps )
        tabulatedForm = self.tabulatedForm.toPointwiseLinear( accuracy, lowerEps = lowerEps, upperEps = upperEps )
        for i, xys in enumerate( tabulatedForm ) :
            if( i == 0 ) : xys.value = mathMisc.shiftFloatUpABit( xys.value, 1e-6 )
            LegendreForm.append( xys )
        return( LegendreForm )

    def toXMLList( self, indent = ""  ) :

        indent2 = indent + '  '
        xmlString = [ '%s<%s Legendre="%s" tabulated="%s">' % ( indent, self.moniker, self.LegendreForm.moniker, self.tabulatedForm.moniker ) ]
        xmlString += self.LegendreForm.toXMLList( indent = indent2 )
        xmlString += self.tabulatedForm.toXMLList( indent = indent2 )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    def toENDF6( self, flags, targetInfo ) :

        from fudge.legacy.converting import endfFormats, gndToENDF6
        ENDFDataList = []

        LI, LTT, frame, d = self.LegendreForm.toENDF6( flags, targetInfo )
        ENDFDataList += d[:-1]

        LI, LTT, frame, d = self.tabulatedForm.toENDF6( flags, targetInfo )
        ENDFDataList += d[:-1]

        ENDFDataList.append( endfFormats.endfSENDLineNumber( ) )
        return( 0, 3, self.LegendreForm.maxLegendreOrder( ), frame, ENDFDataList )

    @staticmethod
    def parseXMLNode( form, linkData={} ):
        legendre, tabulated = form[0], form[1]
        # assume the two regions are 'LegendrePointwise' and 'pointwise':
        legendreClass = {base.LegendrePointwiseFormToken: LegendrePointwise,
                base.LegendrePiecewiseFormToken: LegendrePiecewise}.get( legendre.tag )
        if legendreClass is None: raise Exception("Can't handle %s legendre distribution yet" % legendre.tag)
        legendre = legendreClass.parseXMLNode( legendre, linkData )
        tabulatedClass = {base.pointwiseFormToken: pointwise}.get( tabulated.tag )
        if tabulatedClass is None: raise Exception("Can't handle %s tabulated distribution yet" % tabulated.tag)
        tabulated = tabulatedClass.parseXMLNode( tabulated, linkData )
        return mixedRanges( legendre, tabulated )

class equalProbableBins( base.equalProbableBinsFormBase ) :

    def __init__( self ) :

        base.equalProbableBinsFormBase.__init__( self )

def W_XYs_LegendreSeries_toENDF6( self, flags, targetInfo, lastEnergyData ) :

    from fudge.legacy.converting import endfFormats, gndToENDF6
    counter, ENDFInterpolation, ENDFDataList = 0, gndToENDF6.axesToEndfInterpolationFlag( self.axes ), []
    for energy in self :
        counter += 1
        NW, NL = len( energy ) - 1, 0
        if( targetInfo['doMF4AsMF6'] ) : NL = NW
        energyData = [ endfFormats.endfContLine( 0, energy.value, 0, 0, NW, NL ) ]
        energyData += endfFormats.endfDataList( energy.coefficients[1:] )
        if( counter == 1 ) :
            if( lastEnergyData is not None ) :
                if( lastEnergyData == energyData ) :
                    energyData = []
                    counter -= 1
        ENDFDataList += energyData
    return( counter, ENDFInterpolation, ENDFDataList, energyData )
