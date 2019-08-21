# <<BEGIN-copyright>>
# <<END-copyright>>

""" energy/angular double differential distribution classes """

import math
import base, miscellaneous
import fudge
from fudge.core.math.xData import axes, XYs, W_XYs, V_W_XYs

__metaclass__ = type

KalbachMann_Form_fr_Token = 'fr'
KalbachMann_Form_fra_Token = 'fra'

def parseXMLNode( energyAngularElement, linkData={} ):
    """ translate <energyAngular> element from xml """
    energyAngular = component( energyAngularElement.get("nativeData") )
    for form in energyAngularElement:
        if form.tag==base.KalbachMannFormToken:
            energyAngular.addForm( KalbachMann.parseXMLNode( form, linkData ) )
        elif form.tag==base.pointwiseFormToken:
            energyAngular.addForm( pointwise.parseXMLNode( form, linkData ) )
        else: raise Exception("encountered unknown energyAngular form: %s" % form.tag )
    return energyAngular

class component( base.component ) :

    def __init__( self, nativeData = base.noneFormToken ) :

        base.component.__init__( self, base.energyAngularComponentToken, nativeData = nativeData )

    def calculateDepositionData( self, processInfo, tempInfo ) :

        return( self.forms[self.nativeData].calculateDepositionData( processInfo, tempInfo ) )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        return( self.forms[self.nativeData].process( processInfo, tempInfo, verbosityIndent ) )

    def toENDF6( self, MT, endfMFList, flags, tempInfo ) :

        from fudge.legacy.converting import gndToENDF6
        if( hasattr( self.forms[self.nativeData], 'toENDF6' ) ) :
            nativeData = self.forms[self.nativeData]
            LAW, frame, MF6 = nativeData.toENDF6( flags, tempInfo )
            gndToENDF6.toENDF6_MF6( MT, endfMFList, flags, tempInfo, LAW, frame, MF6 )
        else :
            print 'WARNING: Component, no toENDF6 for nativeData = %s' % self.nativeData, self.forms[self.nativeData].__class__

class pointwise( base.form, V_W_XYs.V_W_XYs ) :

    tag = base.pointwiseFormToken
    moniker = base.pointwiseFormToken

    def __init__( self, axes ) :

        base.form.__init__( self, base.pointwiseFormToken )
        V_W_XYs.V_W_XYs.__init__( self, axes )

    def normalize( self, insitu = True ) :

        n = self
        if( not( insitu ) ) : n = self.copy( )
        for E_MuEpPs in n :
            sum = E_MuEpPs.integrate( )
            for muEpPs in E_MuEpPs : muEpPs.setData( muEpPs / sum )
        return( n )

    def toPointwiseLinear( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        return( V_W_XYs.V_W_XYs.toPointwiseLinear( self, accuracy, lowerEps = lowerEps, upperEps = upperEps, cls = linear ) )

    @classmethod
    def parseXMLNode( self, pwElement, linkData={} ):
        """ translate <pointwise> element from xml """
        from fudge.core.math.xData import axes, W_XYs
        axes_ = axes.parseXMLNode( pwElement[0] )
        pw = pointwise( axes_ )
        for energy_in in pwElement[1:]:
            w_xys = W_XYs.W_XYs( axes.referenceAxes(pw, 3), index=int(energy_in.get('index')),
                    value=float(energy_in.get('value')), parent=pw)
            for energy_out in energy_in:
                data = map(float, energy_out.text.split())
                data = zip( data[::2], data[1::2] )
                xys = XYs.XYs(
                        axes.referenceAxes(w_xys), data, float(energy_out.get('accuracy')),
                        value=float(energy_out.get('value')) )
                w_xys.append( xys )
            pw.append( w_xys )
        return pw

class linear( pointwise ) :

    tag = base.linearFormToken
    moniker = base.linearFormToken

    def __init__( self, axes ) :

        if( not axes.isLinear( qualifierOk = True ) ) : raise Exception( 'interpolation not linear: %s' % axes )
        pointwise.__init__( self, axes )

class KalbachMann( base.form ) :

    KalbachMann_Forms = [ KalbachMann_Form_fr_Token, KalbachMann_Form_fra_Token ]
    KalbachMann_a_parameters = { 'n'     : { 'M' : 1, 'm' : 0.5, 'I' : 0.0  }, 'H1'  : { 'M' : 1, 'm' : 1.0, 'I' : 0.0  },
                                 'H2'    : { 'M' : 1, 'm' : 1.0, 'I' : 2.22 }, 'H3'  : { 'M' : 0, 'm' : 1.0, 'I' : 8.48 },
                                 'He3'   : { 'M' : 0, 'm' : 1.0, 'I' : 7.72 }, 'He4' : { 'M' : 0, 'm' : 2.0, 'I' : 28.3 },
                                 'gamma' : { 'M' : 0, 'm' : 0.0, 'I' : 0.0 } }

    def __init__( self, form, axes ) :

        if( form not in self.KalbachMann_Forms ) : raise Exception( 'Invalid KalbachMann form = %s' % form )
        base.form.__init__( self, base.KalbachMannFormToken )
        self.form = form
        self.axes = axes.copy( parent = self )
        self.energies_in = []

    def __len__( self ) :

        return( len( self.energies_in ) )

    def __getitem__( self, index ) :

        return( self.energies_in[index] )

    def append( self, energy_coefficients ) :

        from fudge.core.utilities import brb
        if( not( isinstance( energy_coefficients, KalbachMannCoefficients ) ) ) : raise Exception( "Argument must be KalbachMannCoefficients: it is %s" % \
            brb.getType( energy_coefficients ) )
        if( len( self ) != 0 ) :
            if( energy_coefficients.value <= self.energies_in[-1].value ) : raise Exception( "value = %s <= last's value = %s" % \
                ( energy_coefficients.value, self.energies_in[-1].value ) )
        self.energies_in.append( energy_coefficients.copy( index = len( self ) ) )

    def check( self, info ):

        from fudge.gnd import warning
        from pqu.physicalQuantityWithUncertainty import PhysicalQuantityWithUncertainty as PQU
        warnings = []

        if self.form!='fr':
            raise Exception("Encountered unexpected Kalbach-Mann form: %s" % self.form)

        ax = axes.axes()
        ax[0] = self.axes[1]; ax[1] = self.axes[2]
        for energy_in in self:
            eprime = energy_in[::3]; f = energy_in[1::3]; r = energy_in[2::3]
            F = XYs.XYs( ax, zip( eprime, f ), 0.001 )
            integral = F.integrate()
            if abs(integral - 1.0) > info['normTolerance']:
                warnings.append( warning.unnormalizedKMDistribution( PQU(energy_in.value, self.axes[0].unit),
                    energy_in.index, integral, energy_in ) )
            if any( [f<0 for f in F] ):
                warnings.append( warning.negativeProbability( PQU(energy_in.value, self.axes[0].unit),
                    obj=energy_in ) )
            if not all( [0<= rv <=1 for rv in r] ):
                badR = min(r) if (0.5-min(r) > max(r)-0.5) else max(r)
                warnings.append( warning.valueOutOfRange("Invalid 'r' in KalbachMann distribution at incident energy %s"
                    % PQU(energy_in.value, self.axes[0].unit), badR, 0, 1, energy_in ) )

        return warnings

    def calculateDepositionData( self, processInfo, tempInfo ) :

        def sqrtE_MuComAverage( f, r, a, tolerance ) :          # This still needs to be tested.?????????

            def u_mu_func( energy_out, parameters ) :

                a = parameters[2].getValue( energy_out )
                if( abs( a ) < 1e-2 ) :
                    a2 = a * a
                    csa = a * ( 315. + a2 * ( -21. + 2. * a2 ) ) / 945.
                else :
                    csa = math.cosh( a ) / math.sinh( a ) - 1. / a
                return( math.sqrt( energy_out ) * parameters[0].getValue( energy_out ) * parameters[1].getValue( energy_out ) * csa )

            epMin = max( f.domainMin( ), r.domainMin( ), a.domainMin( ) )
            epMax = min( f.domainMax( ), r.domainMax( ), a.domainMax( ) )
            sqrtE_mu_com_, quadInfo = miscellaneous.GnG_adaptiveQuadrature( u_mu_func, epMin, epMax, [ f, r, a ], miscellaneous.GaussQuadrature2, tolerance )
            return( sqrtE_mu_com_ )

        energyUnit = tempInfo['incidentEnergyUnit']
        massUnit = energyUnit + '/c**2'
        momentumDepositionUnit = energyUnit + '/c'
        energyAccuracy, momentumAccuracy = 1e-6, 1e-3

        projectile, target, product = tempInfo['reactionSuite'].projectile, tempInfo['reactionSuite'].target, tempInfo['product']
        mass1 = projectile.getMass( massUnit )
        mass2 = target.getMass( massUnit )
        massx = product.getMass( massUnit )
        m1x = mass1 * massx / ( mass1 + mass2 )**2

        multiplicity = tempInfo['multiplicity']
        Es = [ coefficients.value for coefficients in self ]
        depEnergy, depMomentum = [], []
        for E in Es : 
            f, r, a = self.getFRAatEnergy_asLinearPointwise( E )
            E_com = m1x * E
            Ex_com = f.integrateWithWeight_x( )
            sqrtE_mu_com_ = sqrtE_MuComAverage( f, r, a, energyAccuracy )
            E_mu = 2. * math.sqrt( m1x * E ) * sqrtE_mu_com_
            depEnergy.append( [ E, multiplicity.getValue( E ) * ( E_com + Ex_com + E_mu ) ] )
            depMomentum.append( [ E, multiplicity.getValue( E ) * math.sqrt( 2. * massx ) * ( math.sqrt( m1x * E ) + sqrtE_mu_com_ ) ] )

        axes_ = fudge.gnd.productData.energyDeposition.pointwise.defaultAxes( energyUnit = energyUnit, energyDepositionUnit = energyUnit )
        depEnergy = fudge.gnd.productData.energyDeposition.pointwise( axes_, depEnergy, energyAccuracy )

        axes_ = fudge.gnd.productData.momentumDeposition.pointwise.defaultAxes( energyUnit = energyUnit, momentumDepositionUnit = momentumDepositionUnit )
        depMomentum = fudge.gnd.productData.momentumDeposition.pointwise( axes_, depMomentum, momentumAccuracy )
        return( [ depEnergy, depMomentum ] )

    @staticmethod
    def calculate_S_ab_MeV( Z_AB, N_AB, Z_C, N_C, I_ab ) :

        A_AB, A_C = float( Z_AB + N_AB ), float( Z_C + N_C )
        invA_AB_third = 1.0 / math.pow( A_AB, 1.0 / 3.0 )
        invA_C_third = 1.0 / math.pow( A_C, 1.0 / 3.0 )
        NZA_AB = ( N_AB - Z_AB ) * ( N_AB - Z_AB ) / A_AB
        NZA_C = ( N_C - Z_C ) * ( N_C - Z_C ) / A_C

        S =   15.68 * ( A_C - A_AB )                                             - 28.07 * ( NZA_C - NZA_AB ) \
            - 18.56 * ( A_C * invA_C_third - A_AB * invA_AB_third )              + 33.22 * ( NZA_C * invA_C_third - NZA_AB * invA_AB_third ) \
            - 0.717 * ( Z_C * Z_C * invA_C_third - Z_AB * Z_AB * invA_AB_third ) + 1.211 * ( Z_C * Z_C / A_C - Z_AB * Z_AB / A_AB ) \
            - I_ab
        return( S )

    def calculate_a( self, energy_in, energy_out_cmMin, energy_out_cmMax, accuracy = 1e-6 ) :

        from pqu import physicalQuantityWithUncertainty
        from fudge.core.math import fudgemath
        from fudge.core.utilities import fudgeZA

        def getParticleData( particle ) :

            Z, A, suffix, ZA = particle.getZ_A_SuffixAndZA( )
            return( particle.name, Z, max( 0, A - Z ), A, particle.getMass( 'MeV/c**2' ) )

        energyFactor = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( 1, 'MeV' ).getValueAs( self.axes[0].getUnit( ) )
        projectile = self.findAttributeInAncestry( 'projectile' )
        target = self.findAttributeInAncestry( 'target' )
        product = self.findClassInAncestry( fudge.gnd.product.product ).particle
        name_a, Z_a, N_a, A_a, AWRa = getParticleData( projectile )
        name_A, Z_A, N_A, A_A, AWRA = getParticleData( target )
        name_b, Z_b, N_b, A_b, AWRb = getParticleData( product )
        Z_C, N_C = Z_a + Z_A, N_a + N_A
        if( N_A == 0 ) : N_C = 0
        Z_B, N_B = Z_C - Z_b, max( 0, N_C - N_b )
        A_B = Z_B + N_B
        if( N_B == 0 ) : A_B = 0
        particles = self.findClassInAncestry( fudge.gnd.reactionSuite.reactionSuite ).particles
        residual = particles.getParticle( fudgeZA.ZAToGNDName( 1000 * Z_B + A_B ) )
        AWRB = residual.getMass( 'MeV/c**2' )
        Ma, Ia = self.KalbachMann_a_parameters[name_a]['M'], self.KalbachMann_a_parameters[name_a]['I']
        mb, Ib = self.KalbachMann_a_parameters[name_b]['m'], self.KalbachMann_a_parameters[name_b]['I']
        Sa, Sb = energyFactor * self.calculate_S_ab_MeV( Z_A, N_A, Z_C, N_C, Ia ), energyFactor * self.calculate_S_ab_MeV( Z_B, N_B, Z_C, N_C, Ib )

        C1, C2, C3 = 0.04 / energyFactor, 1.8e-6 / energyFactor**3, 6.7e-7 / energyFactor**4

        ea = energy_in * AWRA / ( AWRa + AWRA ) + Sa

        R1, R3 = 130 * energyFactor, 41 * energyFactor                        # MeV to self's energy units.
        if( ea < R1 ) : R1 = ea
        if( ea < R3 ) : R3 = ea

        def calculate_a2( energy_out_cm ) :

            eb = energy_out_cm * ( AWRb + AWRB ) / AWRB + Sb
            X1, X3 = R1 * eb / ea, R3 * eb / ea
            return( X1 * ( C1 + C2 * X1 * X1 ) + C3 * Ma * mb * X3**4 )

        class thicken_a :

            def __init__( self, calculate_a2, relativeTolerance, absoluteTolerance ) :

                self.calculate_a2 = calculate_a2
                self.relativeTolerance = relativeTolerance
                self.absoluteTolerance = absoluteTolerance

            def evaluateAtX( self, x ) :

                return( self.calculate_a2( x ) )

        a = [ [ energy_out_cmMin, calculate_a2( energy_out_cmMin ) ], [ energy_out_cmMax, calculate_a2( energy_out_cmMax ) ] ]
        a = fudgemath.thickenXYList( a, thicken_a( calculate_a2, accuracy, 1e-10 ) )

        axes_ = axes.axes( )
        axes_[0] = self.axes[1].copy( standAlone = True )
        axes_[0].setInterpolation( axes.interpolationXY( axes.linearToken, axes.linearToken ) )
        axes_[1] = axes.axis( 'a', 1, '', frame = self.axes[2].frame )
        return( XYs.XYs( axes_, a, accuracy = accuracy ) )

    def getPointwiseLinear( self, accuracy = 1e-3 ) :
        """Returns a pointwise represent of self in the center of mass frame. The returned data is as a W_XYs instance."""

        def fOfE( axesXY, accuracy, Ep, f, r, a ) :

            f_, r_, a_ = f.getValue( Ep ), r.getValue( Ep ), a.getValue( Ep )
            c = 0.5 * a_ * f_ / math.cosh( a_ )
            f_mu = []
            for imu_cm in xrange( -10, 11 ) :
                mu_cm = max( -1, min( 1, imu_cm / 10. ) )
                f_mu.append( [ mu_cm, c * ( math.cosh( a_ * mu_cm ) + r_ * math.sinh( a_ * mu_cm ) ) ] )
            return( XYs.XYs( axesXY, f_mu, accuracy, value = Ep ) )

        accuracy = min( max( accuracy, 1e-5 ), .2 )
        axes_ = axes.axes( dimension = 4 )
        axes_[0] = self.axes[0]
        axes_[1] = axes.axis( self.axes[1].getLabel( ), 1, self.axes[1].getUnit( ), frame = self.axes[1].getFrame( ), interpolation = axes.interpolationXY( axes.linearToken, axes.linearToken ) )
        axes_[2] = axes.axis( 'mu', 2, '', frame = axes.centerOfMassToken, interpolation = axes.interpolationXY( axes.linearToken, axes.linearToken ) )
        axes_[3] = axes.axis( 'f(%s,%s|%s)' % ( axes_[1].getLabel( ), axes_[2].getLabel( ), axes_[0].getLabel( ) ), 
            3, self.axes[2].getUnit( ), frame = self.axes[2].getFrame( ) )
        f_E_Ep_mu = linear( axes_ )
        axesW_XY = axes.referenceAxes( f_E_Ep_mu, dimension = 3 )
        axesXY = axes.referenceAxes( f_E_Ep_mu, dimension = 2 )
        for e_in in self :
            wxys = W_XYs.W_XYs( axesW_XY, value = e_in.value )
            f, r, a = self.getFRAatEnergy_asLinearPointwise( e_in.value )
            fp = f.union( r )
            fp = fp.union( a )
            for Ep, Fp in fp : wxys.append( fOfE( axesXY, accuracy, Ep, fp, r, a ) )
            f_E_Ep_mu.append( wxys )
        return( f_E_Ep_mu )

    def getFRAatEnergy_asLinearPointwise( self, E ) :

        def getAsLinearPointwise( form, coefficients ) :

            epsilon, fra_accuracy = 1e-8, 1e-6
            n, f, r, a = 4, [], [], []
            if( form == KalbachMann_Form_fr_Token ) : n = 3
            for i in xrange( 0, len( coefficients ), n ) :
                f.append( [ coefficients[i], coefficients[i+1] ] )
                r.append( [ coefficients[i], coefficients[i+2] ] ) 
                if( form == KalbachMann_Form_fra_Token ) : a.append( [ coefficients[i], coefficients[i+3] ] ) 
            axes_ = axes.axes( )
            axes_[0] = self.axes[1]
            axes_[1] = self.axes[2]
            f = XYs.XYs( axes_, f, accuracy = fra_accuracy )
            f = f.changeInterpolation( axes.linearToken, axes.linearToken, lowerEps = epsilon, upperEps = epsilon )

            axes_[1].label = 'r'
            r = XYs.XYs( axes_, r, accuracy = fra_accuracy )
            r = r.changeInterpolation( axes.linearToken, axes.linearToken, lowerEps = epsilon, upperEps = epsilon )
            if( r[-1][1] != 0. ) :
                x, y = r[-1]
                xl = ( 1. - epsilon ) * x
                if( r[-2][0] < ( xl - epsilon * x ) ) : r.setValue( xl, y )
                r.setValue( x, 0. )

            axes_[1].label = 'a'
            if( form == KalbachMann_Form_fra_Token ) :
                a = XYs.XYs( axes_, a, accuracy = fra_accuracy )
                a = a.changeInterpolation( axes.linearToken, axes.linearToken, lowerEps = epsilon, upperEps = epsilon )
            else :
                a = self.calculate_a( E, f.domainMin( ), f.domainMax( ), accuracy = fra_accuracy )

            return( f, r, a )

        coefficients1 = None
        for coefficients2 in self.energies_in :
            if( E <= coefficients2.value ) : break
        if( coefficients1 is None ) : 
            return( getAsLinearPointwise( self.form, coefficients2 ) )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        from fudge.processing.deterministic import transferMatrices

        energyUnit = tempInfo['incidentEnergyUnit']
        massUnit = energyUnit + '/c**2'
        newComponents = []

        if( 'LLNL_MC' in processInfo['styles'] ) :
            pass

        if( 'LLNL_Pn' in processInfo['styles'] ) :
            if( processInfo['verbosity'] >= 30 ) : print '%sGrouping %s' % ( verbosityIndent, self.moniker )
            outputChannel = tempInfo['outputChannel'].outputChannel
            projectile, target, product = tempInfo['reactionSuite'].projectile, tempInfo['reactionSuite'].target, tempInfo['product'].particle
            projectileZA, targetZA = projectile.getZ_A_SuffixAndZA( )[-1], target.getZ_A_SuffixAndZA( )[-1]
            productZA = product.getZ_A_SuffixAndZA( )[-1]
            compoundZA = projectileZA + targetZA
            residualZA = compoundZA - productZA
            particlesData = { 'projectile' : { 'ZA' : projectileZA, 'mass' : projectile.getMass( massUnit ) },
                              'target'     : { 'ZA' : targetZA,     'mass' : target.getMass( massUnit ) },
                              'product'    : { 'ZA' : productZA,    'mass' : product.getMass( massUnit ) },
# The next line is wrong.
                              'residual'   : { 'ZA' : residualZA,   'mass' : target.getMass( massUnit ) },      # ??????? This is wrong!
                              'compound'   : { 'ZA' : compoundZA, 'mass' : projectile.getMass( massUnit ) + target.getMass( massUnit ) } }
            TM_1, TM_E = transferMatrices.KalbachMann_TransferMatrix( processInfo, projectile.getName( ), product.getName( ), tempInfo['crossSection'], 
                particlesData, self, tempInfo['multiplicity'], 
                comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )
            fudge.gnd.miscellaneous.TMs2Form( processInfo, tempInfo, newComponents, TM_1, TM_E, self.axes )

        return( newComponents )

    @classmethod
    def parseXMLNode( self, KalbachMannElement, linkData={} ):
        """ translate <KalbachMann> element from xml """
        from fudge.core.math.xData import axes
        axes_ = axes.parseXMLNode( KalbachMannElement[0] )
        form = KalbachMannElement.get("form")
        KMForm = KalbachMann( form, axes_ )
        for energy_in in KalbachMannElement[1:]:
            KMForm.append( KalbachMannCoefficients( int(energy_in.get("index")), float(energy_in.get("value")),
                map(float, energy_in.text.split()) ) )
        return KMForm

    def toPointwiseLinear( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        return( self.getPointwiseLinear( accuracy ) )

    def toENDF6( self, flags, tempInfo ) :

        from fudge.core.math.xData import axes
        from fudge.legacy.converting import endfFormats
        LEP = 1
        independent, dependent, qualifier = self.axes[1].interpolation.getInterpolationTokens( )
        if( dependent != axes.flatToken ) : LEP = 2
        ENDFDataList = [ endfFormats.endfContLine( 0, 0, 2, LEP, 1, len( self ) ) ]
        ENDFDataList += endfFormats.endfInterpolationList( [ len( self ), 2 ] )
        for energy in self : ENDFDataList += energy.toENDF6( flags, tempInfo )
        return( 1, self.axes[2].frame, ENDFDataList )

    def toXMLList( self, indent = "" ) :

        indent2 = indent + '  '
        xmlString = [ '%s<%s form="%s">' % ( indent, self.moniker, self.form ) ]
        xmlString += self.axes.toXMLList( indent = indent2 )
        for index, data in enumerate( self ) : xmlString += data.toXMLList( indent = indent2 )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

class KalbachMannCoefficients :

    def __init__( self, index, value, coefficients ) :

        self.index = index
        self.value = value
        self.coefficients = coefficients

    def __len__( self ) :

        return( len( self.coefficients ) )

    def __getitem__( self, index ) :

        return( self.coefficients[index] )

    def getValue( self ) :

        return( self.value )

    def copy( self, index = None ) :

        if( index is None ) : index = self.index
        return( KalbachMannCoefficients( index, self.value, self.coefficients ) )

    def toENDF6( self, flags, tempInfo ) :

        from fudge.legacy.converting import endfFormats
        length = len( self.coefficients )
        ENDFDataList = [ endfFormats.endfContLine( 0, self.value, 0, 1, length, length / 3 ) ]
        ENDFDataList += endfFormats.endfDataList( self.coefficients )
        return( ENDFDataList )

    def toXMLList( self, indent = "" ) :

        from pqu import physicalQuantityWithUncertainty
        xmlString = '%s<energy_in value="%s" index="%d" length="%d">' % ( indent, physicalQuantityWithUncertainty.toShortestString( self.value ), self.index, len(self) )
        data = [ physicalQuantityWithUncertainty.toShortestString( datum ) for datum in self.coefficients ]
        xmlString += ' '.join( data ) + '</energy_in>'
        return( [ xmlString ] )
