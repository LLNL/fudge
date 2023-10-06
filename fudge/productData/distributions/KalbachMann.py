# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""Kalbach-Mann double differential distribution classes."""

import math
from fudge.core.math import fudgemath as fudgemathModule

from pqu import PQU as PQUModule

from xData import enums as xDataEnumsModule
from xData import base as xDataBaseModule
from xData import axes as axesModule
from xData import values as valuesModule
from xData import XYs1d as XYs1dModule
from xData import multiD_XYs as multiD_XYsModule

from PoPs import specialNuclearParticleID as specialNuclearParticleIDPoPsModule
from PoPs import IDs as IDsPoPsModule
from PoPs.chemicalElements import misc as chemicalElementMiscPoPsModule

from . import base as baseModule
from . import energy as energyModule
from . import energyAngular as energyAngularModule
from . import miscellaneous as miscellaneousModule

KalbachMann_a_parameters = { IDsPoPsModule.neutron  : { 'M' : 1, 'm' : 0.5, 'I' : 0.0  }, 'H1'  : { 'M' : 1, 'm' : 1.0, 'I' : 0.0  },
                             'H2'                   : { 'M' : 1, 'm' : 1.0, 'I' : 2.22 }, 'H3'  : { 'M' : 0, 'm' : 1.0, 'I' : 8.48 },
                             'He3'                  : { 'M' : 0, 'm' : 1.0, 'I' : 7.72 }, 'He4' : { 'M' : 0, 'm' : 2.0, 'I' : 28.3 },
                             IDsPoPsModule.photon   : { 'M' : 0, 'm' : 0.0, 'I' : 0.0 } }

class XYs2d(multiD_XYsModule.XYs2d):
    '''Special XYs2d class for Kalbach/Mann functions r and a to make sure 'unitbase-unscaled' interpolation qualifier.'''

    def evaluate(self, domainValue, extrapolation=xDataEnumsModule.Extrapolation.none, interpolationQualifier=None, **kwargs):

        function = multiD_XYsModule.XYs2d.evaluate(self, domainValue, extrapolation=extrapolation, 
                interpolationQualifier=xDataEnumsModule.InterpolationQualifier.unitBaseUnscaled, **kwargs)

        return function

class Subform( baseModule.Subform ) :
    """Abstract base class for Kalback/Mann subforms."""

    ancestryMembers = ( 'data', )

    def __init__( self, data ) :

        baseModule.Subform.__init__( self )
        if( ( self.moniker == ASubform.moniker ) and ( data is None ) ) :
            pass
        else :
            if( not( isinstance( data, xDataBaseModule.XDataFunctional ) ) ) : raise TypeError( 'invalid data for KalbachMannCoefficient' )
            if( data.dimension != 2 ) : raise TypeError( 'invalid dimension = %s for KalbachMannCoefficient' % data.dimension )
        self.data = data
        if( data is not None ) : self.data.setAncestor( self )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        if( self.data is not None ) : self.data.convertUnits( unitMap )

    def copy( self ) :

        newData = None
        if self.data is not None: newData = self.data.copy( )
        return( self.__class__( newData ) )

    def fixDomains(self, energyMin, energyMax, domainToFix):
        """
        Calls the **fixDomains** for the **data** member.
        """

        if self.data is not None:
            return self.data.fixDomains(energyMin, energyMax, domainToFix, tweakLower = True)

        return 0

    def isEmptyASubform( self ) :

        return( ( self.moniker == ASubform.moniker ) and ( self.data is None ) )

    def toXML_strList( self, indent = "", **kwargs ) :

        if( self.isEmptyASubform( ) ) : return( [] )
        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ "%s<%s>" % ( indent, self.moniker ) ]
        xmlStringList += self.data.toXML_strList( indent = indent2, **kwargs )
        xmlStringList[-1] += "</%s>" % self.moniker
        return( xmlStringList )

    @staticmethod
    def defaultAxes( moniker, energyUnit ) :

        axes = axesModule.Axes(3)
        axes[2] = axesModule.Axis( 'energy_in',  2, energyUnit )
        axes[1] = axesModule.Axis( 'energy_out', 1, energyUnit )
        if( moniker == 'f' ) :
            axes[0] = axesModule.Axis( moniker, 0, '1/' + energyUnit )
        else :
            axes[0] = axesModule.Axis( moniker, 0, '' )
        return( axes )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """This is currently not called but should be."""

        pass

class FSubform( Subform ) :

    moniker = 'f'

    @staticmethod
    def defaultAxes( energyUnit ) :

        return( Subform.defaultAxes( FSubform.moniker, energyUnit ) )

class RSubform( Subform ) :

    moniker = 'r'

    @staticmethod
    def defaultAxes( energyUnit ) :

        return( Subform.defaultAxes( RSubform.moniker, energyUnit ) )

class ASubform( Subform ) :

    moniker = 'a'

    @staticmethod
    def defaultAxes( energyUnit ) :

        return( Subform.defaultAxes( ASubform.moniker, energyUnit ) )

class Form( baseModule.Form ) :

    moniker = 'KalbachMann'
    subformAttributes = ( 'fSubform', 'rSubform', 'aSubform' )
    ancestryMembers = subformAttributes

    def __init__( self, label, productFrame, _fSubform, _rSubform, _aSubform ) :

        if( not( isinstance( _fSubform, FSubform ) ) ) : raise TypeError( 'invalid Kalbach/Mann f data type' )
        if( not( isinstance( _rSubform, RSubform ) ) ) : raise TypeError( 'invalid Kalbach/Mann r data type' )
        if( not( isinstance( _aSubform, ASubform ) ) ) : raise TypeError( 'invalid Kalbach/Mann a data type' )
        baseModule.Form.__init__( self, label, productFrame, ( _fSubform, _rSubform, _aSubform ) )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.fSubform.convertUnits( unitMap )
        self.rSubform.convertUnits( unitMap )
        self.aSubform.convertUnits( unitMap )

    def check( self, info ) :

        from fudge import warning
        warnings = []

        if( not( self.aSubform.isEmptyASubform( ) ) ) :
            raise NotImplementedError("Checking for Kalbach-Mann data with 'a' coefficients")

        domainMins = set([term.data.domainMin for term in (self.rSubform, self.fSubform, self.aSubform)
                          if term.data is not None])
        domainMaxes = set([term.data.domainMax for term in (self.rSubform, self.fSubform, self.aSubform)
                          if term.data is not None])
        if len(domainMins) != 1 or len(domainMaxes) != 1:
            warnings.append( warning.KalbachMannDomainMismatch(self) )
        for index, F in enumerate(self.fSubform.data):    # F is like P(E' | E), must be normalized for each incident energy
            integral = F.integrate()
            if abs(integral - 1.0) > info['normTolerance']:
                warnings.append( warning.UnnormalizedKMDistribution( PQUModule.PQU( F.outerDomainValue, F.axes[-1].unit), index, integral, F ) )
            if F.rangeMin < 0:
                warnings.append(warning.NegativeProbability(F.rangeMin, PQUModule.PQU(F.outerDomainValue, F.axes[-1].unit), obj=F))
        for R in self.rSubform.data:    # R = pre-compound fraction, must be between 0 and 1
            if R.rangeMin < 0  or  R.rangeMax > 1:
                badR = R.rangeMin if (0.5-R.rangeMin > R.rangeMax - 0.5) else R.rangeMax
                warnings.append(warning.ValueOutOfRange("Invalid 'r' in KalbachMann distribution at incident energy %s"
                    % PQUModule.PQU(R.outerDomainValue, R.axes[-1].unit), badR, 0, 1, R))

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

            epMin = max( f.domainMin, r.domainMin, a.domainMin )
            epMax = min( f.domainMax, r.domainMax, a.domainMax )
            _sqrtE_mu_com, quadInfo = miscellaneousModule.GnG_adaptiveQuadrature( u_mu_func, epMin, epMax, [ f, r, a ],
                miscellaneousModule.GaussQuadrature2, tolerance )
            return( _sqrtE_mu_com )

        def calculateAverageProductDataAtEnergy( self, Ein, suppressNoResidualInPoPs = True ) :

            multiplicity = self.multiplicity.evaluate( Ein )
            if( multiplicity == 0 ) : return( 0.0, 0.0 )
            f, r, a = self.KalbackMannSelf.getFRAatEnergy_asLinearPointwise( Ein, suppressNoResidualInPoPs = suppressNoResidualInPoPs )
            E_com = self.m1x * Ein
            Ex_com = f.integrateWithWeight_x( )
            _sqrtE_mu_com = sqrtE_MuComAverage( f, r, a, energyAccuracy )
            E_mu = 2. * math.sqrt( self.m1x * Ein ) * _sqrtE_mu_com
            return( multiplicity * ( E_com + Ex_com + E_mu ), multiplicity * math.sqrt( 2. * massx ) * ( math.sqrt( m1x * Ein ) + _sqrtE_mu_com ) )

        class CalculateDepositionInfo :

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

        multiplicity = kwargs['multiplicity']
        energyAccuracy = kwargs['energyAccuracy']
        momentumAccuracy = kwargs['momentumAccuracy']

        mass1 = kwargs['projectileMass']
        mass2 = kwargs['targetMass']
        massx = kwargs['productMass']
        EMin = kwargs['EMin']
        EMax = kwargs['EMax']
        
        m1x = mass1 * massx / ( mass1 + mass2 )**2

        calculationData = CalculateDepositionInfo( self, multiplicity, massx, m1x )

        Es = [ coefficients.outerDomainValue for coefficients in self.fSubform.data ]
        if( EMin < Es[0] ) : EMin = Es[0]           # Fix some data issues.
        Es = valuesModule.Values.sortAndThin( set( Es + multiplicity.domainGrid ), 1e-6 )
        while( Es[0] < EMin ) : del Es[0]
        aveEnergy = []
        aveMomentum = []
        suppressNoResidualInPoPs = kwargs.get('suppressNoResidualInPoPs', False)
        for E in Es : 
            Eout, pout = calculateAverageProductDataAtEnergy( calculationData, E, suppressNoResidualInPoPs = suppressNoResidualInPoPs )
            aveEnergy.append( [ E, Eout ] )
            aveMomentum.append( [ E, pout ] )
            suppressNoResidualInPoPs = True

        calculationData.mode = 1
        absoluteTolerance = 1e-3 * energyAccuracy * max( [ Ep for E, Ep in aveEnergy ] )
        calculationData.setTolerances( energyAccuracy, absoluteTolerance )

        aveEnergy = fudgemathModule.thickenXYList( aveEnergy, calculationData )

        calculationData.mode = 2
        absoluteTolerance = 1e-3 * momentumAccuracy * max( [ pp for E, pp in aveMomentum ] )
        calculationData.setTolerances( momentumAccuracy, absoluteTolerance )
        aveMomentum = fudgemathModule.thickenXYList( aveMomentum, calculationData )

        return( [ aveEnergy ], [ aveMomentum ] )

    def calculate_a( self, energy_in, energy_out_cmMin, energy_out_cmMax, accuracy = 1e-6, **kwargs ) :

        suppressNoResidualInPoPs = kwargs.get('suppressNoResidualInPoPs')
        reactionSuite = self.rootAncestor

        def getParticleData( particleID ) :

            particle = reactionSuite.PoPs[particleID]
            if particleID in reactionSuite.PoPs.aliases:
                particle = reactionSuite.PoPs[ particle.pid ]
            Z, A, ZA, level = chemicalElementMiscPoPsModule.ZAInfo( particle )
            mass = particle.getMass( 'MeV/c**2' )
            return( particleID, Z, max( 0, A - Z ), A, mass )

        energyFactor = PQUModule.PQU( 1, 'MeV' ).getValueAs( self.fSubform.data.axes[-1].unit )

        projectileID = self.findAttributeInAncestry( 'projectile' )
        projectileIsPhoton = projectileID == IDsPoPsModule.photon
        if( projectileIsPhoton ) : projectileID = IDsPoPsModule.neutron

        targetID = self.findAttributeInAncestry( 'target' )
        productID = self.product.pid
        name_a, Z_a, N_a, A_a, AWRa = getParticleData( projectileID )
        name_A, Z_A, N_A, A_A, AWRA = getParticleData( targetID )
        name_b, Z_b, N_b, A_b, AWRb = getParticleData( productID )
        Z_C, N_C = Z_a + Z_A, N_a + N_A
        if( N_A == 0 ) : N_C = 0
        Z_B, N_B = Z_C - Z_b, max( 0, N_C - N_b )
        A_B = Z_B + N_B
        if( N_B == 0 ) : A_B = 0

        residualID = chemicalElementMiscPoPsModule.idFromZAndA( Z_B, A_B )
        try :
            residual = reactionSuite.PoPs[residualID]
            numberOfMasses = len( residual.mass )
        except :
            numberOfMasses = 0
        if( numberOfMasses == 0 ) :
            if A_B == 0:
                AWRB = AWRa + AWRA - AWRb
                if not suppressNoResidualInPoPs:
                    print(f'Residual "{residualID}" is not in PoPs. Its mass is estimated from the other masses: AWRB = {AWRB:.2f} MeV.')
            else:
                AWRB = PQUModule.PQU( A_B, 'amu' ).getValueAs( 'MeV/c**2' )
                if not suppressNoResidualInPoPs:
                    print(f'Residual "{residualID}" is not in PoPs. Its mass is estimated from its mass number: AWRB = {AWRB:.2f} MeV.')
        else :
            AWRB = residual.getMass( 'MeV/c**2' )

        name_a_ID = specialNuclearParticleIDPoPsModule.specialNuclearParticleID(name_a, specialNuclearParticleIDPoPsModule.Mode.nuclide)
        Ma, Ia = KalbachMann_a_parameters[name_a_ID]['M'], KalbachMann_a_parameters[name_a_ID]['I']
        name_b_ID = specialNuclearParticleIDPoPsModule.specialNuclearParticleID(name_b, specialNuclearParticleIDPoPsModule.Mode.nuclide)
        mb, Ib = KalbachMann_a_parameters[name_b_ID]['m'], KalbachMann_a_parameters[name_b_ID]['I']
        Sa = self.calculate_S_ab_MeV( Z_A, N_A, Z_C, N_C, Ia )
        Sb = self.calculate_S_ab_MeV( Z_B, N_B, Z_C, N_C, Ib )

        C1 = 0.04
        C2 = 1.8e-6
        C3 = 6.7e-7

        energy_in /= energyFactor
        ea = energy_in * AWRA / ( AWRa + AWRA ) + Sa

        R1 = 130
        R3 = 41
        if( ea < R1 ) : R1 = ea
        if( ea < R3 ) : R3 = ea

        def calculate_a2( energy_out_cm ) :

            eb = energy_out_cm * ( AWRb + AWRB ) / AWRB + Sb
            X1, X3 = R1 * eb / ea, R3 * eb / ea
            return( X1 * ( C1 + C2 * X1 * X1 ) + C3 * Ma * mb * X3**4 )

        class Thicken_a :

            def __init__( self, calculate_a2, relativeTolerance, absoluteTolerance ) :

                self.calculate_a2 = calculate_a2
                self.relativeTolerance = relativeTolerance
                self.absoluteTolerance = absoluteTolerance

            def evaluateAtX( self, x ) :

                return( self.calculate_a2( x ) )

        energy_out_cmMin /= energyFactor
        energy_out_cmMax /= energyFactor
        a = [ [ energy_out_cmMin, calculate_a2( energy_out_cmMin ) ], [ energy_out_cmMax, calculate_a2( energy_out_cmMax ) ] ]
        a = fudgemathModule.thickenXYList( a, Thicken_a( calculate_a2, accuracy, 1e-10 ) )

        if( projectileIsPhoton ) :
            factor = math.sqrt( energy_in / ( 2.0 * AWRa ) )
            for EpA in a :
                if( EpA[0] == 0.0 ) :
                    EpA[1] *= factor * 4.0
                else :
                    EpA[1] *= factor * min( 4.0, max( 1.0, 9.3 / math.sqrt( EpA[0] ) ) )

        for EpA in a : EpA[0] *= energyFactor
        axes = ASubform.defaultAxes( self.fSubform.data.axes[1].unit )
        return( XYs1dModule.XYs1d( data = a, axes = axes, outerDomainValue = energy_in * energyFactor ) )

    def copy( self ) :

        return( Form( self.label, self.productFrame, self.fSubform.copy( ), self.rSubform.copy( ), self.aSubform.copy( ) ) )

    __copy__ = copy

    def domainUnitConversionFactor( self, unitTo ) :

        if( unitTo is None ) : return( 1. )
        return( PQUModule.PQU( '1 ' + self.domainUnit ).getValueAs( unitTo ) )

    @property
    def domainMin(self):

        return self.fSubform.data.domainMin

    @property
    def domainMax(self):

        return self.fSubform.data.domainMax

    @property
    def domainUnit(self):

        return self.fSubform.data.axes[-1].unit

    def muProbability( self, mu, _r, _a ) :
        """
        Calculates the angular probabilty at mu given the Kalbach/Mann parameters r and a. Note, the caller is required to 
        have evaluated r(E,E') and a(E,E') at the required incident energy E and outgoing energy E' befure calling this function.
        """

        amu = _a * mu
        return( 0.5 * _a * ( math.cosh( amu ) + _r * math.sinh( amu ) ) / math.sinh( _a ) )

    def energySpectrumAtEnergy(self, energyIn, frame, **kwargs):

        muMin = kwargs.get('muMin', -1.0)
        muMax = kwargs.get('muMax',  1.0)

        if frame == xDataEnumsModule.Frame.centerOfMass:
            if muMin == -1 and muMax == -1:
                return self.fSubform.data.evaluate(energyIn)

            projectileIsPhoton = self.findAttributeInAncestry('projectile') == IDsPoPsModule.photon
            fAtEnergy, rAtEnergy, aAtEnergy = self.getFRAatEnergy_asLinearPointwise(energyIn, **kwargs)
            data = []
            for energyOut, PofEnergyOut in fAtEnergy:
                rValue = rAtEnergy.evaluate(energyOut)
                aValue = aAtEnergy.evaluate(energyOut)
                aMuMin = aValue * muMin
                aMuMax = aValue * muMax
                if projectileIsPhoton:
                    factor = 0.5 * ((1.0 - rValue) * (muMax - muMin) + rValue * aValue * (math.exp(aMuMax) - math.exp(aMuMin)) / math.sinh(aValue))
                else:
                    factor = 0.5 * ((math.sinh(aMuMax) - math.sinh(aMuMin)) +
                            rValue * ((math.cosh(aValue *muMax) - math.cosh(aValue * muMin))) / math.sinh(aValue))
                data.append([energyOut, factor * PofEnergyOut])
        else:
            xys2d = self.spectrumAtEnergy(energyIn, xDataEnumsModule.Frame.lab)
            data = [[xys1d.outerDomainValue, xys1d.integrate(domainMin=muMin, domainMax=muMax)] for xys1d in xys2d]

        return energyModule.XYs1d(data=data, axes=energyModule.defaultAxes(self.domainUnit))

    def fixDomains(self, energyMin, energyMax, domainToFix):
        """
        Calls the **fixDomains** for the **fSubform**, **rSubform** and **aSubform** members.
        """

        numberOfFixes  = self.fSubform.fixDomains(energyMin, energyMax, domainToFix)
        numberOfFixes += self.rSubform.fixDomains(energyMin, energyMax, domainToFix)
        numberOfFixes += self.aSubform.fixDomains(energyMin, energyMax, domainToFix)

        return numberOfFixes

    def spectrum( self, frame, **kwargs ) :

        if frame != xDataEnumsModule.Frame.lab:
            return( self.toPointwise_withLinearXYs( **kwargs ) )

        energiesIn = self.fSubform.data.domainGrid + self.rSubform.data.domainGrid
        if( not( self.aSubform.isEmptyASubform( ) ) ) : energiesIn += self.aSubform.data.domainGrid
        energiesIn = valuesModule.Values.sortAndThin( energiesIn, rel_tol = 1e-12 )
        xys3d = energyAngularModule.XYs3d( axes = energyAngularModule.defaultAxes( self.domainUnit ) )
        for energy in energiesIn : xys3d.append( self.spectrumAtEnergy( energy, frame ) )
        return( xys3d )

    def spectrumAtEnergy( self, energyIn, frame, **kwargs ) :

        class EnergyAngualarAtEnergyCOM :

            def __init__( self, KalbachMannSelf, rAtEnergy, aAtEnergy ) :

                self.KalbachMannSelf = KalbachMannSelf
                self.rAtEnergy = rAtEnergy
                self.aAtEnergy = aAtEnergy

            def probabilityCOM( self, energyPrimeCOM, muCOM ) :

                _r = rAtEnergy.evaluate( energyPrimeCOM )
                _a = aAtEnergy.evaluate( energyPrimeCOM )
                return( self.KalbachMannSelf.muProbability( muCOM, _r, _a ) )

        def fOfMu( projectileIsPhoton, EnergyOut, probability, rAtEnergy, aAtEnergy ) :

            _r = rAtEnergy.evaluate( EnergyOut )
            _a = aAtEnergy.evaluate( EnergyOut )
            if( projectileIsPhoton ) :
                norm1 = 0.5 * probability
                norm2 = _r * _a / math.sinh( _a )
            else :
                norm1 = 0.5 * _a * probability / math.sinh( _a )
            f_mu = []
            for imu_cm in range( -10, 11 ) :
                mu_cm = max( -1, min( 1, imu_cm / 10. ) )
                if( projectileIsPhoton ) :
                    f_mu.append( [ mu_cm, norm1 * ( 1.0 - _r + norm2 * math.exp( _a * mu_cm ) ) ] )
                else :
                    f_mu.append( [ mu_cm, norm1 * ( math.cosh( _a * mu_cm ) + _r * math.sinh( _a * mu_cm ) ) ] )
            return( energyAngularModule.XYs1d( data = f_mu, outerDomainValue = EnergyOut ) )

        projectileIsPhoton = self.findAttributeInAncestry( 'projectile' ) == IDsPoPsModule.photon
        fAtEnergy, rAtEnergy, aAtEnergy = self.getFRAatEnergy_asLinearPointwise( energyIn, **kwargs )

        if frame == xDataEnumsModule.Frame.centerOfMass:
            xys2d = energyAngularModule.XYs2d( outerDomainValue = energyIn )
            fAtEnergyUnion = fAtEnergy.union( rAtEnergy )
            fAtEnergyUnion = fAtEnergyUnion.union( aAtEnergy )
            for EnergyOut, probability in fAtEnergyUnion : xys2d.append( fOfMu( projectileIsPhoton, EnergyOut, probability, rAtEnergy, aAtEnergy ) )
            return( xys2d )

        return( miscellaneousModule.energyAngularSpectrumFromCOMSpectrumToLabAtEnergy( self, energyIn, fAtEnergy, EnergyAngualarAtEnergyCOM( self, rAtEnergy, aAtEnergy ) ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        Returns a pointwise represent of self in the center of mass frame. The returned data is as a XYs2d instance.
        This method will be deprecated, use spectrum method instead.
        """

        def fOfMu( accuracy, projectileIsPhoton, Ep, f, r, a ) :

            _f, _r, _a = f.evaluate( Ep ), r.evaluate( Ep ), a.evaluate( Ep )
            if( projectileIsPhoton ) :
                norm1 = 0.5 * _f
                norm2 = _r * _a / math.sinh( _a )
            else :
                norm1 = 0.5 * _a * _f / math.sinh( _a )
            f_mu = []
            for imu_cm in range( -10, 11 ) :
                mu_cm = max( -1, min( 1, imu_cm / 10. ) )
                if( projectileIsPhoton ) :
                    f_mu.append( [ mu_cm, norm1 * ( 1.0 - _r + norm2 * math.exp( _a * mu_cm ) ) ] )
                else :
                    f_mu.append( [ mu_cm, norm1 * ( math.cosh( _a * mu_cm ) + _r * math.sinh( _a * mu_cm ) ) ] )
            return( energyAngularModule.XYs1d( data = f_mu, outerDomainValue = Ep ) )

        accuracy = xDataBaseModule.getArguments( kwargs, { 'accuracy' : XYs1dModule.defaultAccuracy } )['accuracy']
        accuracy = min( max( accuracy, 1e-5 ), .2 )

        projectileIsPhoton = self.findAttributeInAncestry( 'projectile' ) == IDsPoPsModule.photon

        axes = energyAngularModule.defaultAxes( energyUnit = self.fSubform.data.axes[2].unit )
        xys3d = energyAngularModule.XYs3d( axes = axes )

        egrid = self.fSubform.data.domainGrid + self.rSubform.data.domainGrid       # get list of incident energies
        if( not( self.aSubform.isEmptyASubform( ) ) ) : egrid += self.aSubform.data.domainGrid
        egrid = sorted( set( egrid ) )
        for e_in in egrid :
            xys2d = energyAngularModule.XYs2d( outerDomainValue = e_in )
            f, r, a = self.getFRAatEnergy_asLinearPointwise( e_in, **kwargs )
            fp = f.union( r )
            fp = fp.union( a )
            for Ep, Fp in fp : xys2d.append( fOfMu( accuracy, projectileIsPhoton, Ep, fp, r, a ) )
            xys3d.append( xys2d )
        return( xys3d )

    def getFRAatEnergy_asLinearPointwise( self, E, **kwargs ) :

        epsilon, fra_accuracy = 1e-8, 1e-6
        _f = self.fSubform.data.evaluate( E )
        f = _f.changeInterpolation(xDataEnumsModule.Interpolation.linlin, XYs1dModule.defaultAccuracy, lowerEps=epsilon, upperEps=epsilon)
        if _f.interpolation == xDataEnumsModule.Interpolation.flat:
            f.setValue( f[-1][0], f[-2][1] )

        r = self.rSubform.data.evaluate(E, interpolationQualifier=xDataEnumsModule.InterpolationQualifier.unitBase)
        r = r.changeInterpolation(xDataEnumsModule.Interpolation.linlin, XYs1dModule.defaultAccuracy, lowerEps=epsilon, upperEps=epsilon)
# BRB removed 20/Dec/2017. Do not understand why this logic is here. It is okay for r[-1][1] to be 0.
#        if( r[-1][1] != 0. ) :
#            x, y = r[-1]
#            xl = ( 1. - epsilon ) * x
#            if( r[-2][0] < ( xl - epsilon * x ) ) : r.setValue( xl, y )
#            r.setValue( x, 0. )

        if( self.aSubform.isEmptyASubform( ) ) :
            a = self.calculate_a( E, f.domainMin, f.domainMax, accuracy = fra_accuracy, **kwargs )
        else :
            a = self.aSubform.data.evaluate( E )
            a = a.changeInterpolation(xDataEnumsModule.Interpolation.linlin, XYs1dModule.defaultAccuracy, lowerEps=epsilon, upperEps=epsilon)
        return( f, r, a )

    def processMC_cdf( self, style, tempInfo, indent, **kwargs ) :

        oldData = self.fSubform.data
        subform = energyModule.XYs2d( axes = oldData.axes, interpolation = oldData.interpolation,
                interpolationQualifier = oldData.interpolationQualifier )
        for xys in oldData :
            subform.append(energyModule.Xs_pdf_cdf1d.fromXYs(xys, xys.outerDomainValue, thinEpsilon=1e-14))
        _fSubform = FSubform(subform)

        if( self.aSubform.isEmptyASubform( ) ) :
            KalbachMann_a_Axes = ASubform.defaultAxes( oldData.domainUnit )
            _aSubform = XYs2d(axes=KalbachMann_a_Axes, interpolation=xDataEnumsModule.Interpolation.linlin,
                interpolationQualifier=xDataEnumsModule.InterpolationQualifier.unitBaseUnscaled)
            suppressNoResidualInPoPs = kwargs.get('suppressNoResidualInPoPs', False)
            for f_xys in oldData :
                _aSubform.append( self.calculate_a( f_xys.outerDomainValue, f_xys.domainMin, f_xys.domainMax, 
                        suppressNoResidualInPoPs = suppressNoResidualInPoPs, **kwargs ) )
                suppressNoResidualInPoPs = True
            _aSubform = ASubform( _aSubform )
        else :
            _aSubform = self.aSubform.copy( )               # BRB FIXME, f, a and r need to be on the same energy (E and E') grid.
        _form = Form( style.label, self.productFrame, _fSubform, self.rSubform.copy( ), _aSubform )
        return( _form )

    def processMultiGroup( self, style, tempInfo, indent ) :

        from fudge.processing import group as groupModule
        from fudge.processing.deterministic import transferMatrices as transferMatricesModule

        verbosity = tempInfo['verbosity']
        indent2 = indent + tempInfo['incrementalIndent']
        reactionSuite = tempInfo['reactionSuite']

        product = reactionSuite.PoPs[tempInfo['product'].pid]
        productLabel = tempInfo['productLabel']

        energyUnit = tempInfo['incidentEnergyUnit']
        massUnit = energyUnit + '/c**2'
        mass2MeVFactor = PQUModule.PQU( 1, energyUnit ).getValueAs( 'MeV' )

        if( verbosity > 2 ) : print('%sGrouping %s' % (indent, self.moniker))
        projectileZA = tempInfo['projectileZA']
        targetZA = tempInfo['targetZA']
        if( targetZA == 6000 ) :
            targetZA = 6012
            print('    Kludge for C_natural: changing targetZA from 6000 to %d' % targetZA)
        productZA = chemicalElementMiscPoPsModule.ZA( product )
        compoundZA = projectileZA + targetZA
        residualZA = compoundZA - productZA
        particlesData = { 'projectile' : { 'ZA' : projectileZA },
                          'target'     : { 'ZA' : targetZA },
                          'product'    : { 'ZA' : productZA },
                          'residual'   : { 'ZA' : residualZA },
                          'compound'   : { 'ZA' : compoundZA, 'mass' : tempInfo['projectileMass'] + tempInfo['targetMass'] } }

        residualMass = tempInfo['masses']['Residual']           # Save old value.
        tempInfo['masses']['Residual'] = self.residualMass( reactionSuite.PoPs, residualZA, massUnit, particlesData['compound']['mass'], product )

        masses = tempInfo['masses']
        tempInfo['masses'] = {}
        for particle in masses : tempInfo['masses'][particle] = masses[particle] * mass2MeVFactor
        particlesData['compound']['mass'] *= mass2MeVFactor

        try :
            TM_1, TM_E = transferMatricesModule.KalbachMann_TransferMatrix( style, tempInfo, tempInfo['crossSection'], 
                    particlesData, self, tempInfo['multiplicity'], 
                    comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )
        except :
            tempInfo['masses'] = masses
            raise
        tempInfo['masses'] = masses

        tempInfo['masses']['Residual'] = residualMass
        return( groupModule.TMs2Form( style, tempInfo, TM_1, TM_E ) )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        """Translate <KalbachMann> element from xml."""

        xPath.append( element.tag )
        childArrays = {}
        for child in element :
            _subformClass = None
            for subformClass in ( FSubform, RSubform, ASubform ) :
                if subformClass.moniker == child.tag:
                    _subformClass = subformClass
                    break
            if _subformClass is None:
                raise TypeError("Unexpected element '%s' encountered" % child.tag)
            if( child.tag == FSubform.moniker ) :
                xData = energyModule.XYs2d.parseNodeUsingClass(child[0], xPath, linkData, **kwargs)
            else :
                xData = XYs2d.parseNodeUsingClass(child[0], xPath, linkData, **kwargs)
            childArrays[ '_%sSubform' % child.tag ] = _subformClass( xData )
        if '_aSubform' not in childArrays: childArrays['_aSubform'] = ASubform( None )

        KM = cls( element.get( 'label' ), element.get( 'productFrame' ), **childArrays )

        xPath.pop()
        return KM

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
