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

"""Kalbach-Mann double differential distribution classes."""

import math
import base, miscellaneous
import fudge
from fudge.core.utilities import brb
from fudge.core.math import fudgemath as fudgemathModule

from pqu import PQU

import xData.base as xDataBaseModule
import xData.axes as axesModule
import xData.XYs as XYsModule
import xData.standards as standardsModule
import xData.multiD_XYs as multiD_XYsModule

from fudge.gnd.productData import energyDeposition as energyDepositionModule
from fudge.gnd.productData import momentumDeposition as momentumDepositionModule

from . import base as baseModule

__metaclass__ = type

KalbachMann_Form_fr_Token = 'fr'
KalbachMann_Form_fra_Token = 'fra'

class subform( baseModule.subform ) :
    """Abstract base class for energyAngular subforms."""

    def __init__( self, data ) :

        if( ( self.moniker == aSubform.moniker ) and ( data is None ) ) :
            pass
        else :
            if( not( isinstance( data, xDataBaseModule.xDataFunctional ) ) ) : raise TypeError( 'invalid data for KalbachMannCoefficient' )
            if( data.dimension != 2 ) : raise TypeError( 'invalid dimension = %s for KalbachMannCoefficient' % data.dimension )
        self.data = data

    def copy( self ) :

        return( self.__class__( self.data ) )

    def isEmptyASubform( self ) :

        return( ( self.moniker == aSubform.moniker ) and ( self.data is None ) )

    def toXMLList( self, indent = "", **kwargs ) :

        if( self.isEmptyASubform( ) ) : return( [] )
        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ "%s<%s>" % ( indent, self.moniker ) ]
        xmlStringList += self.data.toXMLList( indent = indent2, **kwargs )
        xmlStringList[-1] += "</%s>" % self.moniker
        return( xmlStringList )

    @staticmethod
    def defaultAxes( moniker ) :

        axes = axesModule.axes( rank = 3 )
        axes[2] = axesModule.axis( 'energy_in',  2, 'eV' )
        axes[1] = axesModule.axis( 'energy_out', 1, 'eV' )
        x0Unit = ''
        if( moniker == 'f' ) : x0Unit = '1/eV'
        axes[0] = axesModule.axis( moniker, 0, x0Unit )
        return( axes )

class fSubform( subform ) :
    # BRB, FIXME, data should be an energyModule.XYs2d instance.

    moniker = 'f'

    @staticmethod
    def defaultAxes( ) :

        return( subform.staticmethod( self.moniker ) )

class rSubform( subform ) :

    moniker = 'r'

    @staticmethod
    def defaultAxes( ) :

        return( subform.staticmethod( self.moniker ) )

class aSubform( subform ) :

    moniker = 'a'

    @staticmethod
    def defaultAxes( ) :

        return( subform.staticmethod( self.moniker ) )

class form( baseModule.form ) :

    moniker = 'KalbachMann'
    subformAttributes = ( 'fSubform', 'rSubform', 'aSubform' )
    KalbachMann_a_parameters = { 'n'     : { 'M' : 1, 'm' : 0.5, 'I' : 0.0  }, 'H1'  : { 'M' : 1, 'm' : 1.0, 'I' : 0.0  },
                                 'H2'    : { 'M' : 1, 'm' : 1.0, 'I' : 2.22 }, 'H3'  : { 'M' : 0, 'm' : 1.0, 'I' : 8.48 },
                                 'He3'   : { 'M' : 0, 'm' : 1.0, 'I' : 7.72 }, 'He4' : { 'M' : 0, 'm' : 2.0, 'I' : 28.3 },
                                 'gamma' : { 'M' : 0, 'm' : 0.0, 'I' : 0.0 } }

    def __init__( self, label, productFrame, _fSubform, _rSubform, _aSubform, makeCopy = True ) :

        if( not( isinstance( _fSubform, fSubform ) ) ) : raise TypeError( 'invalid Kalbach/Mann f data type' )
        if( not( isinstance( _rSubform, rSubform ) ) ) : raise TypeError( 'invalid Kalbach/Mann r data type' )
        if( not( isinstance( _aSubform, aSubform ) ) ) : raise TypeError( 'invalid Kalbach/Mann a data type' )
        if( makeCopy ) :
            _fSubform = _fSubform.copy( )
            _rSubform = _rSubform.copy( )
            _aSubform = _aSubform.copy( )
        baseModule.form.__init__( self, label, productFrame, ( _fSubform, _rSubform, _aSubform ) )

    def check( self, info ) :

        from fudge.gnd import warning
        warnings = []

        if( not( self.aSubform.isEmptyASubform( ) ) ) :
            raise NotImplementedError("Checking for Kalbach-Mann data with 'a' coefficients")

        for index, F in enumerate( self.fSubform.data ):    # F is like P(E' | E), must be normalized for each incident energy
            integral = F.integrate()
            if abs(integral - 1.0) > info['normTolerance']:
                warnings.append( warning.unnormalizedKMDistribution( PQU.PQU(F.value, F.axes[-1].unit),
                    index, integral, F ) )
            if F.rangeMin() < 0:
                warnings.append( warning.negativeProbability( PQU.PQU(F.value, F.axes[-1].unit),
                    obj=F ) )
        for R in self.rSubform.data:    # R = pre-compound fraction, must be between 0 and 1
            if R.rangeMin() < 0  or  R.rangeMax() > 1:
                badR = R.rangeMin() if (0.5-R.rangeMin() > R.rangeMax() - 0.5) else R.rangeMax()
                warnings.append( warning.valueOutOfRange("Invalid 'r' in KalbachMann distribution at incident energy %s"
                    % PQU.PQU(R.value, R.axes[-1].unit), badR, 0, 1, R ) )

        return warnings

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        def sqrtE_MuComAverage( f, r, a, tolerance ) :          # This still needs to be tested.?????????

            def u_mu_func( energy_out, parameters ) :

                a = parameters[2].evaluate( energy_out )
                if( abs( a ) < 1e-2 ) :
                    a2 = a * a
                    csa = a * ( 315. + a2 * ( -21. + 2. * a2 ) ) / 945.
                else :
                    csa = math.cosh( a ) / math.sinh( a ) - 1. / a
                return( math.sqrt( energy_out ) * parameters[0].evaluate( energy_out ) * parameters[1].evaluate( energy_out ) * csa )

            epMin = max( f.domainMin( ), r.domainMin( ), a.domainMin( ) )
            epMax = min( f.domainMax( ), r.domainMax( ), a.domainMax( ) )
            _sqrtE_mu_com, quadInfo = miscellaneous.GnG_adaptiveQuadrature( u_mu_func, epMin, epMax, [ f, r, a ], miscellaneous.GaussQuadrature2, tolerance )
            return( _sqrtE_mu_com )

        def calculateAverageProductDataAtEnergy( self, Ein ) :

            f, r, a = self.KalbackMannSelf.getFRAatEnergy_asLinearPointwise( Ein )
            E_com = self.m1x * Ein
            Ex_com = f.integrateWithWeight_x( )
            _sqrtE_mu_com = sqrtE_MuComAverage( f, r, a, energyAccuracy )
            E_mu = 2. * math.sqrt( self.m1x * Ein ) * _sqrtE_mu_com
            multi = self.multiplicity.getValue( Ein )
            return( multi * ( E_com + Ex_com + E_mu ), multi * math.sqrt( 2. * massx ) * ( math.sqrt( m1x * Ein ) + _sqrtE_mu_com ) )

        class calculateDepositionInfo :

            def __init__( self, KalbackMannSelf, multiplicity, massx, m1x ) :

                self.KalbackMannSelf = KalbackMannSelf
                self.multiplicity = multiplicity
                self.massx = massx
                self.m1x = m1x
                self.mode = 0

            def evaluateAtX( self, E ) :

                Eout, pout = calculateAverageProductDataAtEnergy( self, E )
                if( self.mode == 0 ) : return( Eout, pout )
                if( self.mode == 1 ) : return( Eout )
                return( pout )

            def setTolerances( self, relativeTolerance, absoluteTolerance ) :

                self.relativeTolerance = relativeTolerance
                self.absoluteTolerance = absoluteTolerance

        verbosity = kwargs.get( 'verbosity', 0 )
        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        energyUnit = kwargs['incidentEnergyUnit']
        momentumDepositionUnit = energyUnit + '/c'
        massUnit = energyUnit + '/c**2'
        multiplicity = kwargs['multiplicity']
        productMass = kwargs['product'].getMass( massUnit )
        energyAccuracy = kwargs['energyAccuracy']
        momentumAccuracy = kwargs['momentumAccuracy']
        reactionSuite = kwargs['reactionSuite']
        product = kwargs['product']
        mass1 = reactionSuite.projectile.getMass( massUnit )
        mass2 = reactionSuite.target.getMass( massUnit )
        massx = product.getMass( massUnit )
        EMin = kwargs['EMin']
        EMax = kwargs['EMax']
        
        m1x = mass1 * massx / ( mass1 + mass2 )**2

        calculationData = calculateDepositionInfo( self, multiplicity, massx, m1x )

        Es = [ coefficients.value for coefficients in self.fSubform.data ]
        if( EMin < Es[0] ) : EMin = Es[0]           # Fix some data issues.
        Es = sorted( set( Es + multiplicity.domainGrid( ) ) )
        while( Es[0] < EMin ) : del Es[0]
        aveEnergy = []
        aveMomentum = []
        for E in Es : 
            Eout, pout = calculateAverageProductDataAtEnergy( calculationData, E )
            aveEnergy.append( [ E, Eout ] )
            aveMomentum.append( [ E, pout ] )

        calculationData.mode = 1
        absoluteTolerance = 1e-3 * energyAccuracy * max( [ Ep for E, Ep in aveEnergy ] )
        calculationData.setTolerances( energyAccuracy, absoluteTolerance )
        aveEnergy = fudgemathModule.thickenXYList( aveEnergy, calculationData )

        calculationData.mode = 2
        absoluteTolerance = 1e-3 * momentumAccuracy * max( [ pp for E, pp in aveMomentum ] )
        calculationData.setTolerances( momentumAccuracy, absoluteTolerance )
        aveMomentum = fudgemathModule.thickenXYList( aveMomentum, calculationData )

        return( [ aveEnergy ], [ aveMomentum ] )

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

        from fudge.core.math import fudgemath
        from fudge.particles import nuclear

        def getParticleData( particle ) :

            Z, A, suffix, ZA = particle.getZ_A_SuffixAndZA( )
            return( particle.name, Z, max( 0, A - Z ), A, particle.getMass( 'MeV/c**2' ) )

        energyFactor = PQU.PQU( 1, 'MeV' ).getValueAs( self.fSubform.data.axes[-1].unit )
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
        residual = particles.getParticle( nuclear.nucleusNameFromZA( 1000 * Z_B + A_B ) )
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

        axes = axesModule.axes( )
        axes[0] = axesModule.axis( 'a', 0, '' )
        axes[1] = self.fSubform.data.axes[1].copy( )
        return( XYsModule.XYs1d( data = a, axes = axes, accuracy = accuracy ) )

    def copy( self ) :

        return( KalbachMann( self.fSubform, self.rSubform, self.aSubform, makeCopy = True ) )

    __copy__ = __deepcopy__ = copy

    def domainUnitConversionFactor( self, unitTo ) :

        if( unitTo is None ) : return( 1. )
        return( PQU.PQU( '1 ' + self.domainUnit( ) ).getValueAs( unitTo ) )

    def domainMin( self, unitTo = None, asPQU = False ) :

        from fudge.gnd import product
        return( self.findClassInAncestry( product.product ).domainMin( unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        from fudge.gnd import product
        return( self.findClassInAncestry( product.product ).domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def domain( self, unitTo = None, asPQU = False ) :

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ), self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def domainUnit( self ) :

        return( self.axes[-1].unit )

    def getPointwiseLinear( self, accuracy = 1e-3 ) :
        """Returns a pointwise represent of self in the center of mass frame. The returned data is as a XYs2d instance."""

        def fOfE( accuracy, Ep, f, r, a ) :

            f_, r_, a_ = f.evaluate( Ep ), r.evaluate( Ep ), a.evaluate( Ep )
            c = 0.5 * a_ * f_ / math.cosh( a_ )
            f_mu = []
            for imu_cm in xrange( -10, 11 ) :
                mu_cm = max( -1, min( 1, imu_cm / 10. ) )
                f_mu.append( [ mu_cm, c * ( math.cosh( a_ * mu_cm ) + r_ * math.sinh( a_ * mu_cm ) ) ] )
            return( pdfOfMu.pointwise( data = f_mu, accuracy = accuracy, value = Ep ) )

        accuracy = min( max( accuracy, 1e-5 ), .2 )
        # FIXME check for consistency between units on f,r and a parameters?
        axes_ = pointwise.defaultAxes(energyUnit='eV', energy_outUnit='eV', probabilityUnit='1/eV')
        f_E_Ep_mu = pointwise( axes=axes_ )
        # get list of incident energies
        egrid = self.fSubform.data.domainGrid() + self.rSubform.data.domainGrid()
        if( not( self.aSubform.isEmptyASubform( ) ) ) : egrid += self.aSubform.data.domainGrid()
        egrid = sorted( set( egrid ) )
        for e_in in egrid :
            wxys = pdfOfEpAndMu.pointwise( value = e_in )
            f, r, a = self.getFRAatEnergy_asLinearPointwise( e_in )
            fp = f.union( r )
            fp = fp.union( a )
            for Ep, Fp in fp : wxys.append( fOfE( accuracy, Ep, fp, r, a ) )
            f_E_Ep_mu.append( wxys )
        return( f_E_Ep_mu )

    def getFRAatEnergy_asLinearPointwise( self, E ) :

        epsilon, fra_accuracy = 1e-8, 1e-6
        f = self.fSubform.data.evaluate( E )
        f = f.changeInterpolation( standardsModule.interpolation.linlinToken, lowerEps = epsilon, upperEps = epsilon )

        r = self.rSubform.data.evaluate( E )
        r = r.changeInterpolation( standardsModule.interpolation.linlinToken, lowerEps = epsilon, upperEps = epsilon )
        if( r[-1][1] != 0. ) :
            x, y = r[-1]
            xl = ( 1. - epsilon ) * x
            if( r[-2][0] < ( xl - epsilon * x ) ) : r.setValue( xl, y )
            r.setValue( x, 0. )

        if( self.aSubform.isEmptyASubform( ) ) :
            a = self.calculate_a( E, f.domainMin( ), f.domainMax( ), accuracy = fra_accuracy )
        else :
            a = self.aSubform.data.evaluate( E )
            a = a.changeInterpolation( standardsModule.interpolation.linlinToken, lowerEps = epsilon, upperEps = epsilon )
        return( f, r, a )

    def processSnMultiGroup( self, style, tempInfo, indent ) :

        from fudge.processing import group as groupModule
        from fudge.processing.deterministic import transferMatrices as transferMatricesModule

        verbosity = tempInfo['verbosity']
        indent2 = indent + tempInfo['incrementalIndent']
        projectile = tempInfo['reactionSuite'].projectile
        target = tempInfo['reactionSuite'].target
        product = tempInfo['product'].particle
        productLabel = tempInfo['productLabel']

        energyUnit = tempInfo['incidentEnergyUnit']
        massUnit = energyUnit + '/c**2'

        if( verbosity > 2 ) : print '%sGrouping %s' % ( verbosityIndent, self.moniker )
        projectileZA = projectile.getZ_A_SuffixAndZA( )[-1]
        targetZA = target.getZ_A_SuffixAndZA( )[-1]
        productZA = product.getZ_A_SuffixAndZA( )[-1]
        compoundZA = projectileZA + targetZA
        residualZA = compoundZA - productZA
        particlesData = { 'projectile' : { 'ZA' : projectileZA },
                          'target'     : { 'ZA' : targetZA },
                          'product'    : { 'ZA' : productZA },
# The next line is wrong.
                          'residual'   : { 'ZA' : residualZA },             # ??????? This is wrong!
                          'compound'   : { 'ZA' : compoundZA, 'mass' : projectile.getMass( massUnit ) + target.getMass( massUnit ) } }
        residualMass = tempInfo['masses']['Residual']
        tempInfo['masses']['Residual'] = target.getMass( massUnit )         # ??????? This is wrong!
        TM_1, TM_E = transferMatricesModule.KalbachMann_TransferMatrix( style, tempInfo, tempInfo['crossSection'], 
                particlesData, self, tempInfo['multiplicity'], 
                comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )

        tempInfo['masses']['Residual'] = residualMass
        return( groupModule.TMs2Form( style, tempInfo, TM_1, TM_E ) )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):
        """Translate <KalbachMann> element from xml."""

        xPath.append( element.tag )
        childArrays = {}
        for child in element :
            _subformClass = None
            for subformClass in ( fSubform, rSubform, aSubform ) :
                if subformClass.moniker == child.tag:
                    _subformClass = subformClass
                    break
            if _subformClass is None:
                raise TypeError("Unexpected element '%s' encountered" % child.tag)
            xData = multiD_XYsModule.XYs2d.parseXMLNode( child[0], xPath, linkData )
            childArrays[ '_%sSubform' % child.tag ] = _subformClass( xData )
        if '_aSubform' not in childArrays: childArrays['_aSubform'] = aSubform( None )
        KM = form( element.get('label'), element.get('productFrame'), **childArrays )
        xPath.pop()
        return KM

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        return( self.getPointwiseLinear( accuracy ) )
