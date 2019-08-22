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

import math
from fudge.gnd import baseClasses
from . import base
from fudge.gnd import miscellaneous, tokens
from fudge.core.math.xData import XYs, axes, regions

__metaclass__ = type

#: List of valid cross section forms 
crossSectionForms = [ tokens.pointwiseFormToken,                tokens.linearFormToken,             tokens.piecewiseFormToken, 
                      tokens.resonancesWithBackgroundFormToken, tokens.weightedPointwiseFormToken,  tokens.groupedFormToken, 
                      tokens.referenceFormToken ] 

lowerEps = 1e-8
upperEps = 1e-8

#
# crossSection genre and forms
#
class component( baseClasses.componentBase ) :

    genre = base.crossSectionToken

    def __init__( self, nativeData = tokens.unknownFormToken ) :

        baseClasses.componentBase.__init__( self, crossSectionForms, nativeData = nativeData )

    def domainMin( self, unitTo = None, asPQU = False ) :

        return( self.forms[self.nativeData].domainMin( unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        return( self.forms[self.nativeData].domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def domainUnitConversionFactor( self, unitTo ) :

        from pqu import physicalQuantityWithUncertainty

        if( unitTo is None ) : return( 1. )
        return( physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( '1 ' + self.getDomainUnit( ) ).getValueAs( unitTo ) )

    def getDomain( self, unitTo = None, asPQU = False ) :

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ), self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getDomainUnit( self ) :

        return( self.forms[self.nativeData].getDomainUnit( ) )

    def getIncidentEnergyUnit( self ) :

        return( self.forms[self.nativeData].getIncidentEnergyUnit( ) )

    def heat( self, temperature, EMin, lowerlimit = None, upperlimit = None, 
            interpolationAccuracy = 0.002, heatAllPoints = False, 
            doNotThin = True, heatBelowThreshold = True, heatAllEDomain = True ) : 
        """
        Returns the result of self.toPointwise_withLinearXYs( ).heat( ... ). See method pointwise.heat for more information.
        """

        return( self.toPointwise_withLinearXYs( ).heat( self, temperature, EMin, 
            lowerlimit, upperlimit, interpolationAccuracy, heatAllPoints, doNotThin,
            heatBelowThreshold, heatAllEDomain ) )

    def check( self, info ) :
        """
        Check cross section data for correct threshold, negative cross sections, etc.
        Returns a list of any warnings encountered during checking.
        
        :param dict info:  A Python dictionary containing the parameters that control the cross section checking.
        :keyword float Q: a parameter of `info`: if Q is positive, cross section must start at 1e-5 eV, otherwise, cross section threshold must agree with Q-value
        :keyword float kinematicFactor: a parameter of `info`:
        :keyword string dThreshold: a parameter of `info`: allowable threshold energy mismatch as a string suitable for the PQU class
        :keyword float crossSectionEnergyMax: a parameter of `info`: the cross section must extend up to limit (usually 20 MeV)
        """
        from fudge.gnd import warning
        from pqu.physicalQuantityWithUncertainty import PhysicalQuantityWithUncertainty as PQU
        warnings = []

        lower, upper = self.getDomain( asPQU = True )
        if 'Q' in info:
            # if Q is positive, cross section must start at 1e-5 eV
            # otherwise, cross section threshold must agree with Q-value:
            thresh = PQU( -info['Q'] * info['kinematicFactor'], 'eV' )
            if thresh.value>0:
                if abs(thresh-lower) > PQU( info['dThreshold'] ):
                    warnings.append( warning.threshold_mismatch( thresh, lower, self ) )
            elif lower != PQU(1e-5,'eV'):
                # ignore 2nd,3rd,4th-chance fission (they have Q>0 but are still threshold reactions):
                from fudge.gnd import channels
                if (not hasattr(self.parent, 'outputChannel') or (self.parent.outputChannel.getFissionGenre()
                        in (None, channels.fissionGenreTotal, channels.fissionGenreFirstChance) ) ):
                    warnings.append( warning.threshold_mismatch( lower, PQU(1e-5, 'eV'), self ) )

        # cross section must extend up to limit (usually 20 MeV):
        if upper < PQU( info['crossSectionEnergyMax'] ):
            warnings.append( warning.gapInCrossSection( upper, PQU( info['crossSectionEnergyMax'] ),
                self ) )

        if self.nativeData == tokens.resonancesWithBackgroundFormToken:
            # ensure no gaps between resonance and background:
            resDomain = info['reactionSuite'].resonances.getDomain( asPQU = True )
            bckDomain = self[ tokens.resonancesWithBackgroundFormToken ].tabulatedData.getDomain( asPQU = True )
            if bckDomain[0] > resDomain[1]:
                warnings.append( warning.gapInCrossSection(resDomain[1],bckDomain[0], self ) )

        # can it be converted to pointwise linear?
        if tokens.linearFormToken not in self.forms:
            try:
                self.forms[tokens.linearFormToken] = self.forms[ self.nativeData ].toPointwise_withLinearXYs(0,1e-8)
            except Exception as e:
                warnings.append( warning.ExceptionRaised( e, self ) )
                return warnings

        # following tests require a linearized cross section:
        # test for negative values, and for non-zero cross section at threshold
        pw = self.forms[tokens.linearFormToken]
        if pw.yMin() < 0:
            for i,(en,xsc) in enumerate(pw):
                if xsc < 0:
                    warnings.append( warning.negativeCrossSection( PQU(en, pw.axes[0].unit), i, self ) )

        if 'Q' in info and thresh.value>0:
            if pw[0][1] != 0:
                warnings.append( warning.nonZero_crossSection_at_threshold( pw[0][1], self ) )

        return warnings

    def process( self, processInfo, tempInfo, verbosityIndent ) :
        '''Process self into a form suitable for a specific application'''

        crossSection = self.toPointwise_withLinearXYs( )
        forms = crossSection.process( self, processInfo, tempInfo, verbosityIndent )
        for form in forms : self.addForm( form )
        return( crossSection )
        
    def computeSpectrumAverage( self, spectrum, Emin=None, Emax=None, covariance=None ):
        '''
        Computes the spectrum (i.e. weighted) average of self with the spectrum (i.e. weight) 
        specified by ``spectrum``  optionally from Emin to Emax.  If the covariance is 
        provided, the uncertainty on the spectrum average will be computed.
        
        :param spectrum: spectrum over which to average
        :type spectrum: XYs instance
        :param Emin: Lower integration limit.  
            If None (the default), then the lower limit of self's domain is used 
        :type Emin: PQU or None
        :param Emax: Upper integration limit.  
            If None (the default), then the upper limit of self's domain is used
        :type Emax: PQU or None
        :param covariance: covariance to use when computing uncertainty on the spectral average.  
            If None (default: None), no uncertainty is computed.
        :type covariance: covariance instance or None
        :rtype: PQU
        
        
        :How does it work?:
            Given a  weighting spectrum :math:`\phi(E)`, we intend to average the cross section in self as follows:
        
            .. math::
                < \sigma > = \int_{E_{min}}^{E_{max}} dE \phi(E) \sigma(E)
        
            To compute the uncertainty, we resort to the basis function expansion of the covariance
            as follows:
        
            .. math::
                \Delta^2\sigma(E,E') = \sum_{ij} \Delta^2\sigma_{ij} B_i(E) B_j(E')
            
            In the ENDF format, all covariances are assumed to be grouped in energy so the 
            basis functions are simple window functions.   With this,  
        
            .. math::
                \delta <\sigma> = \sqrt{ \int dE dE' \phi(E)\phi(E') \Delta^2\sigma(E,E') } 
                                = \sqrt{ \sum_{ij} \Delta^2\sigma_{ij} <B_i> <B_j> }
        
        .. warning:: Not implemented yet
        '''
        raise NotImplementedError()
        # Check that the covariance goes with the data
        # Check that the spectrum is given as an XYs instance
        # Convolve spectrum and self to get mean value
        # Compute uncertainty on convolution given mean value and covariance
    
    def computeRI( self, Ec='0.5 eV', covariance=None ):
        '''
        Compute the resonance integral (RI):
        
        .. math::
            RI = \int_{Ec}^{\infty} \sigma(E) dE/E
            
        Where the lower cut-off is usually taken to be the Cadmium cut-off energy of 
        Ec = 0.5 eV (see S. Mughabghab, Atlas of Neutron Resonances).  If the covariance 
        is provided, the uncertainty on the resonance integral will be computed.
        
        :param string Ec: lower bound of integration.  
            Usually taken to be the Cadmium cut-off energy (default: '0.5 eV')
        :param covariance: covariance to use when computing uncertainty on the spectral average.  
            If None (default: None), no uncertainty is computed.
        :type covariance: covariance instance or None
        :rtype: PQU
        
        .. warning:: Not implemented yet
        '''
        raise NotImplementedError()
        # Construct an XYs instance that spans the same domain as self.
        # The y values are all 1/E and the x values are all E.
        oneOverE = 0.0
        return self.spectrumAverage( oneOverE, Emin=Ec, covariance=covariance )
        
    def computeMACS( self, T, covariance=None ):
        '''
        Compute the Maxwellian average cross section (MACS):

        .. math::
            MACS(kT) = {\\frac{2}{\sqrt{\pi}}} {\\frac{a^2}{(kT)^2}} \int_0^{\infty}dE_n^{lab} \sigma(E_n^{lab})*E_n^{lab}*\exp{(-a E_n^{lab}/kT)}

        Here :math:`E_n^{lab}` is the incident neutron energy in the lab frame and :math:`a=m_2/(m_1+m_2)`.
        If the covariance is provided, the uncertainty on the MACS will be computed.

        :param PQU T: Temperature (as a physical quantity of the Maxwellian with which to average
        :param covariance: covariance to use when computing uncertainty on the spectral average.  
            If None (default: None), no uncertainty is computed.
        :type covariance: covariance instance or None
        :rtype: PQU
        
        .. warning:: Not implemented yet
        '''
        raise NotImplementedError()
        # Construct an XYs instance that spans the same domain as self.
        # The x values are the E values and the y values are a Maxwellian evaluated at x.
        maxwellian = 0.0
        return self.spectrumAverage( maxwellian, covariance=covariance )

    def computeCf252SpectrumAve( self, covariance=None ):
        '''
        Compute the Cf252 spontaneous fission neutron spectrum average:
        
        .. math::
            ave = \exp( -0.789 * E ) * \sinh( \sqrt{ 0.447 * E } )
        
        (where did this come from?  Found it in the sigma QA system I think. DAB 12/1/2012)
            
        If the covariance is provided, the uncertainty on the average will be computed.
        
        :param covariance: covariance to use when computing uncertainty on the spectral average.  
            If None (default: None), no uncertainty is computed.
        :type covariance: covariance instance or None
        :rtype: PQU
        
        .. warning:: Not implemented yet
        '''
        raise NotImplementedError()
        # Construct an XYs instance that spans the same domain as self.
        # The x values are the E values and the y values are a Cf252 spectrum evaluated at x.
        E=0 #Gotta write this
        Cf252SpectrumAve = math.exp( -0.789 * E ) * math.sinh( math.sqrt( 0.447 * E ) )
        return self.computeSpectrumAverage( Cf252SpectrumAve, covariance=covariance )

    def computeFourteenMeVPoint( self, covariance=None ):
        '''
        Compute the value of the cross section at 14.2 MeV.

        If the covariance is provided, the uncertainty on the 14.2 MeV point will be computed.
        
        :param covariance: covariance to use when computing uncertainty on the spectral average.  
            If None (default: None), no uncertainty is computed.
        :type covariance: covariance instance or None
        :rtype: PQU
        '''
        raise NotImplementedError()
        # Check that covariance goes with this data
        return self.getValue( '14.2 MeV' ), covariance.getUncertaintyVector( self, relative=False ).getValue( '14.2 MeV' )

    def computeRoomTempCS( self, covariance=None ): 
        '''
        The cross section exactly at room temperature E=0.0235339347844 eV or v=2200 m/s

        If the covariance is provided, the uncertainty on the 0.0235339347844 eV point will be computed.
        
        RoomTemp = 273.1 degK = 0.0235339347844 eV since kB = 8.6173324e-5 eV/degK

        :param covariance: covariance to use when computing uncertainty on the spectral average.  
            If None (default: None), no uncertainty is computed.
        :type covariance: covariance instance or None
        :rtype: PQU
        '''
        raise NotImplementedError()
        # Check that covariance goes with this data
        return self.getValue( '0.0235339347844 eV' ), covariance.getUncertaintyVector( self, relative=False ).getValue( '0.0235339347844 eV' )

    def computeWestcottFactor( self, covariance=None ):
        '''
        Wescott G-factor is the ratio of Maxwellian averaged cross section 
        (at room temparature) and the room temperature cross section.  
        Should be pretty close to 1 if cross section goes like 1/v.

        :param covariance: covariance to use when computing uncertainty on the spectral average.  
            If None (default: None), no uncertainty is computed.
        :type covariance: covariance instance or None
        :rtype: PQU
        '''
        macs,dmacs = self.computeMACS( "298.15 K", covariance )
        rt, drt = self.computeRoomTempCS( covariance ) 
        westcott = macs/rt
        return westcott, westcott * math.sqrt( pow( dmacs/macs, 2.0 ) + pow( drt/rt, 2.0 ) )

    def computeAstrophysicalReactionRate( self, T, covariance=None ):
        '''
        The astrophysical reaction rate R is defined as :math:`R = N_A <\sigma v>`, where 
        :math:`N_A` is Avagadro's number.  It is typically reported in units of [cm^3/mole s].

        :param PQU T: temperature as a physical quantity
        :param covariance: covariance to use when computing uncertainty on the spectral average.  
            If None (default: None), no uncertainty is computed.
        :type covariance: covariance instance or None
        :rtype: PQU
        
        .. warning:: Not implemented yet
        '''
        
        raise NotImplementedError()
    
    def computeStellerEnhancementFactor( self, covariance=None ):
        '''
        :param covariance: covariance to use when computing uncertainty on the spectral average.  
            If None (default: None), no uncertainty is computed.
        :type covariance: covariance instance or None
        :rtype: PQU
        
        .. warning:: Not implemented yet
        '''
        raise NotImplementedError()
      
    def computeScatteringRadius( self, covariance=None ):
        """
        R' the scattering radius
        
        :param covariance: covariance to use when computing uncertainty on the spectral average.  
            If None (default: None), no uncertainty is computed.
        :type covariance: covariance instance or None
        :rtype: PQU
        
        .. warning:: Not implemented yet
        """
        raise NotImplementedError()
    
    def computeSWaveScatteringLength( self, covariance=None ): 
        '''
        S0, the S-wave scattering length???
        
        :param covariance: covariance to use when computing uncertainty on the spectral average.  
            If None (default: None), no uncertainty is computed.
        :type covariance: covariance instance or None
        :rtype: PQU
        
        .. warning:: Not implemented yet

        '''
        raise NotImplemented()
    
    def computePWaveScatteringLength( self, covariance=None ): 
        '''
        S1, the P-wave scattering length???
        
        :param covariance: covariance to use when computing uncertainty on the spectral average.  
            If None (default: None), no uncertainty is computed.
        :type covariance: covariance instance or None     
        :rtype: PQU
        
        .. warning:: Not implemented yet
        '''
        raise NotImplemented()

    def toENDF6( self, MT, endfMFList, targetInfo, level, LR ) :
        '''
        Convert self into ENDF format
        
        :param int MT: The ENDF reaction designator, MT
        :param endfMFList:
        :param targetInfo:
        :param level:
        :param LR:
        '''

        from fudge.legacy.converting import endfFormats
        ZA, mass, QI = targetInfo['ZA'], targetInfo['mass'], targetInfo['Q']
        QM = targetInfo['QM']
        if( QM is None ) :
            if( MT in ( 2, 5 ) ) :
                QM = QI
            elif( MT == 4 ) :                               # Q should be 0 except for excited-state targets:
                if( targetInfo['reactionSuite'].target.getLevelIndex() > 0 ) :
                    QM = QI
                else :
                    QM = 0
            else :
                QM = QI + level
        interpolationFlatData, crossSectionFlatData = self.forms[self.nativeData].toENDF6Data( MT, endfMFList, targetInfo, level )
        endfMFList[3][MT] = [ endfFormats.endfHeadLine( ZA, mass, 0, 0, 0, 0 ) ]
        endfMFList[3][MT].append( endfFormats.endfContLine( QM, QI, 0, LR, len( interpolationFlatData ) / 2, len( crossSectionFlatData ) / 2 ) )
        endfMFList[3][MT] += endfFormats.endfInterpolationList( interpolationFlatData )
        endfMFList[3][MT] += endfFormats.endfDataList( crossSectionFlatData )
        endfMFList[3][MT].append( endfFormats.endfSENDLineNumber( ) )

class baseCrossSectionForm( baseClasses.formBase ) :

    genre = component.genre

    def __init__( self ) :

        baseClasses.formBase.__init__( self )

class pointwise( baseCrossSectionForm, XYs.XYs ) :

    tag = tokens.pointwiseFormToken
    form = tokens.pointwiseFormToken
    mutableYUnit = False

    def __init__( self, axes, data, accuracy, **kwargs ) :

        baseCrossSectionForm.__init__( self )
        kwargs['isPrimaryXData'] = True
        XYs.XYs.__init__( self, axes, data, accuracy, **kwargs )

    def getIncidentEnergyUnit( self ) :

        return( self.axes[0].getUnit( ) )

    def heat( self, ancestors, temperature, EMin, lowerlimit = None, upperlimit = None, interpolationAccuracy = 0.002, heatAllPoints = False, 
            doNotThin = True, heatBelowThreshold = True, heatAllEDomain = True ) :
        """
        Returns a linear version of the cross section heated to 'temperature'. If the current temperature of the cross section
        is greater than 'temperature' a raise is executed.
        If lowerlimit is None, it is set to 'oneOverV' except when the reaction is determined to be a threshold reaction, 
        then it is set to 'threshold'. Any cross section with domainMin greater than 2.5e-4 eV is determined to be a
        threshold reaction. If upperlimit is None it is set to 'constant'.
        If heatBelowThreshold is False, then EMin is set to the larger of EMin and self's domainMin.

        For more information on EMin, lowerlimit, upperlimit, interpolationAccuracy, heatAllPoints, doNotThin and heatAllEDomain
        see the module crossSectionAdjustForHeatedTarget.
        """

        from crossSectionAdjustForHeatedTarget import heat
        from pqu.physicalQuantityWithUncertainty import PhysicalQuantityWithUncertainty

        gotT = currentTemperature = ancestors.findAttributeInAncestry( 'getTemperature' )( )
        if( PhysicalQuantityWithUncertainty.isTemperature( gotT ) ) :
            gotT = PhysicalQuantityWithUncertainty( gotT.getValueAs( 'K' ), 'K' ) * PhysicalQuantityWithUncertainty( '1 k' )
        gotT = gotT.getValueAs( self.getIncidentEnergyUnit( ) )

        if( isinstance( temperature, str ) ) : temperature = PhysicalQuantityWithUncertainty( temperature )
        temperature_ = temperature
        if( PhysicalQuantityWithUncertainty.isTemperature( temperature ) ) :
            temperature_ = PhysicalQuantityWithUncertainty( temperature.getValueAs( 'K' ), 'K' ) * PhysicalQuantityWithUncertainty( '1 k' )
        wantedT = temperature_.getValueAs( self.getIncidentEnergyUnit( ) )

        dT = wantedT - gotT
        if( abs( dT ) <= 1e-2 * wantedT ) : dT = 0.
        if( dT < 0 ) : raise Exception( 'Current temperature "%s" (%.4e) higher than desired temperature "%s" (%.4e)' % ( currentTemperature, gotT, temperature, wantedT ) ) 

        projectile, target = ancestors.findAttributeInAncestry( 'projectile' ), ancestors.findAttributeInAncestry( 'target' )
        if( projectile.getMass( 'amu' ) == 0 ) : raise Exception( 'Heating with gamma as projectile not supported.' )
        massRatio = target.getMass( 'amu' ) / projectile.getMass( 'amu' )

        heated = unheated = self
        if( not( unheated.axes.isLinear( ) ) ) : unheated = self.toPointwise_withLinearXYs( lowerEps, upperEps )
        if( ( len( unheated ) > 0 ) and ( dT > 0. ) ) :
            if( isinstance( EMin, str ) ) : EMin = PhysicalQuantityWithUncertainty( EMin )
            EMin_ = EMin.getValueAs( self.getIncidentEnergyUnit( ) )
            if( not( heatBelowThreshold ) ) : EMin_ = max( unheated.domainMin( ), EMin_ )
            if( lowerlimit == None ) :
                lowerlimit = "oneOverV"
                domainMin = PhysicalQuantityWithUncertainty( 1e-5, 'eV' ).getValueAs( self.getIncidentEnergyUnit( ) )
                if( domainMin < 0.04 * unheated.domainMin( ) ) : lowerlimit = "threshold"
            if( upperlimit == None ) : upperlimit = 'constant'

            heated = heat.crossSectionAdjustForHeatedTarget( massRatio, dT, EMin_, unheated, lowerlimit = lowerlimit, \
                upperlimit = upperlimit, interpolationAccuracy = interpolationAccuracy, heatAllPoints = heatAllPoints, doNotThin = doNotThin,
                heatAllEDomain = heatAllEDomain )
        accuracy = max( interpolationAccuracy, self.getAccuracy( ) )
        return( linear( unheated.axes, heated, accuracy ) )

    def toPointwise_withLinearXYs( self, lowerEps, upperEps ) :

        return( XYs.XYs.toPointwise_withLinearXYs( self, lowerEps = lowerEps, upperEps = upperEps, cls = linear ) )

    def toENDF6Data( self, MT, endfMFList, targetInfo, level ) :

        from fudge.legacy.converting import gndToENDF6
        endfInterpolation = gndToENDF6.axesToEndfInterpolationFlag( self.axes )
        crossSectionFlatData = []
        for xy in self.copyDataToXYs( xUnit = 'eV', yUnit = 'b' ) : crossSectionFlatData += xy
        return( [ len( crossSectionFlatData ) / 2, endfInterpolation ], crossSectionFlatData )

    def toXMLList( self, indent = "" ) :

        return( XYs.XYs.toXMLList( self, indent = indent, incrementalIndent = '  ', pairsPerLine = 100 ) )

    @staticmethod
    def defaultAxes( energyUnit = 'eV', energyInterpolation = axes.linearToken, crossSectionUnit = 'b', crossSectionInterpolation = axes.linearToken ) :

        axes_ = axes.axes( dimension = 2 )
        axes_[0] = axes.axis( 'energy_in', 0, energyUnit, interpolation = axes.interpolationXY( energyInterpolation, crossSectionInterpolation ) )
        axes_[1] = axes.axis( base.crossSectionToken, 1, crossSectionUnit )
        return( axes_ )

class linear( pointwise ) :

    tag = tokens.linearFormToken
    form = tokens.linearFormToken

    def __init__( self, axes_, data, accuracy, **kwargs ) :

        if( not axes_.isLinear( qualifierOk = False ) ) :
            independent, dependent, qualifier = axes_[0].interpolation.getInterpolationTokens( )
            raise Exception( "invalid interpolation = (%s, %s, %s) for class %s: must be (%s, %s, %s)" %
                ( independent, dependent, qualifier, self.form, axes.linearToken, axes.linearToken, None ) )
        pointwise.__init__( self, axes_, data, accuracy, **kwargs )

    def process( self, component_, processInfo, tempInfo, verbosityIndent ) :

        from fudge.processing import miscellaneous

        if( processInfo['verbosity'] >= 20 ) : print '%sprocessing pointwise cross section' % verbosityIndent
        projectile, target = processInfo.getProjectileName( ), processInfo.getTargetName( )
        forms = []

        if( 'LLNL_MC' in processInfo['styles'] ) :
            try :               # Does linear cross section already exists in cross section component. If so, do nothing.
                crossSection = component_.getFormByToken( tokens.linearFormToken )
            except :
                try :           # Does pointwise cross section already exists in cross section component. If so, check if it is linear.
                    crossSection = component_.getFormByToken( tokens.pointwiseFormToken )
                    if( crossSection.axes.isLinear( qualifierOk = False ) ) : crossSection = None
                except :        # Else, add self as it is linear pointwise
                    crossSection = self
                if( crossSection is not None ) : forms.append( crossSection.toPointwise_withLinearXYs( 1e-8, 1e-8 ) )

        if( 'LLNL_Pn' in processInfo['styles'] ) :
            crossSection = miscellaneous.groupOneFunctionAndFlux( projectile, processInfo, self )
            tempInfo['groupedCrossSectionNorm'] = crossSection

            groupedFlux = tempInfo['groupedFlux']
            crossSectionGrouped_ = [ crossSection[i] / groupedFlux[i] for i in xrange( len( crossSection ) ) ]

            axes_ = self.axes.copy( standAlone = True )
            independent, dependent, qualifier = axes_[0].interpolation.getInterpolationTokens( )
            axes_[0].interpolation = axes.interpolationXY( independent, axes.flatToken )
            forms.append( grouped( axes_, crossSectionGrouped_ ) )

        return( forms )

class piecewise( baseCrossSectionForm, regions.regionsXYs ) :
    """This class stores a cross section in two or more regions, which may have different interpolations.
    Each region must contain at least two points. Each pair of adjacent regions must overlap at exactly one point."""

    moniker = tokens.piecewiseFormToken
    form = tokens.piecewiseFormToken
    tag = tokens.piecewiseFormToken

    def __init__( self, axes ) :

        baseCrossSectionForm.__init__( self )
        regions.regionsXYs.__init__( self, axes )

    def getIncidentEnergyUnit( self ) :

        return( self[0].axes[0].getUnit( ) )

    def toPointwise_withLinearXYs( self, lowerEps, upperEps ) :

        accuracy = 1e-3
        if( len( self ) > 0 ) : accuracy = self[0].getAccuracy( )
        xys = regions.regionsXYs.toPointwise_withLinearXYs( self, accuracy, lowerEps, upperEps )
        return( linear( xys.axes, xys, accuracy = xys.getAccuracy( ) ) )

    def toENDF6Data( self, MT, endfMFList, targetInfo, level ) :

        from fudge.legacy.converting import gndToENDF6
        interpolationFlatData, crossSectionFlatData = [], []
        counter = 0
        lastX, lastY = None, None
        for region in self :
            ENDFInterpolation = gndToENDF6.axesToEndfInterpolationFlag( region.axes )
            data = region.copyDataToXYs( xUnit = 'eV', yUnit = 'b' )
            if( lastX is not None ) :
                if( lastY == data[0][1] ) : data = data[1:]
            counter += len( data )
            interpolationFlatData.append( counter )
            interpolationFlatData.append( ENDFInterpolation )
            for xy in data : crossSectionFlatData += xy
            lastX, lastY = data[-1]
        return( interpolationFlatData, crossSectionFlatData )

    @staticmethod
    def defaultAxes( energyUnit = 'eV', energyInterpolation = axes.byRegionToken, crossSectionUnit = 'b', crossSectionInterpolation = axes.byRegionToken ) :

        axes_ = axes.axes( dimension = 2 )
        axes_[0] = axes.axis( 'energy_in', 0, energyUnit, interpolation = axes.interpolationXY( energyInterpolation, crossSectionInterpolation ) )
        axes_[1] = axes.axis( base.crossSectionToken, 1, crossSectionUnit )
        return( axes_ )

class grouped( baseClasses.groupedFormBase ) :

    genre = component.genre

    def __init__( self, axes_, groupData ) :

        baseClasses.groupedFormBase.__init__( self, axes_, groupData )


class resonancesWithBackground( baseCrossSectionForm ) :
    """This class stores cross sections that include resonances along with a background contribution.
    Contains a link to the resonances, and the 'background' which could be pointwise or piecewise data.
    The full pointwise cross section can be obtained by first reconstructing resonances
    and then adding the background contribution (users should use the reactionSuite.reconstructResonances method)."""

    form = tokens.resonancesWithBackgroundFormToken

    def __init__( self, tabulatedData, resonanceLink ) :
        
        baseCrossSectionForm.__init__( self )
        self.link = resonanceLink
        self.tabulatedData = tabulatedData

    def domainMin( self, unitTo = None, asPQU = False ) :

        return( self.tabulatedData.domainMin( unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        return( self.tabulatedData.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getDomain( self, unitTo = None, asPQU = False ) :

        return( self.tabulatedData.getDomain( unitTo = unitTo, asPQU = asPQU ) )

    def getDomainUnit( self ) :

        return( self.tabulatedData.getDomainUnit( ) )
    
    def getIncidentEnergyUnit( self ) :

        return( self.tabulatedData.getDomainUnit( ) )

    def toXMLList( self, indent = "" ) :

        xmlList = [ indent + '<%s>' % self.form, self.link.toXML( indent + '  ' )  ]
        xmlList += self.tabulatedData.toXMLList( indent = indent+'  ' )
        xmlList[-1] += '</%s>' % self.form
        return( xmlList )

    @staticmethod
    def parseXMLNode( element, xPath=[], linkData={} ):
        xPath.append( element.tag )
        from fudge.gnd import link
        xlink = link.parseXMLNode( element[0] )
        form = element[1]
        formClass = {tokens.pointwiseFormToken: pointwise,
                tokens.linearFormToken: linear,
                tokens.piecewiseFormToken: piecewise,
                tokens.weightedPointwiseFormToken: weightedPointwise,
                }.get( form.tag )
        if formClass is None: raise Exception("encountered unknown cross section form: %s" % form.tag)
        resWithBack = resonancesWithBackground( formClass.parseXMLNode(form,xPath,linkData), xlink )
        xPath.pop()
        return resWithBack
    
    def toENDF6Data( self, MT, endfMFList, targetInfo, level ) :
        
        return self.tabulatedData.toENDF6Data( MT, endfMFList, targetInfo, level )

    def toPointwise_withLinearXYs( self, lowerEps, upperEps ) :

        component_ = self.findClassInAncestry( component )
        if( tokens.linearFormToken in component_.forms ) : return( component_.forms[tokens.linearFormToken] )
        raise Exception( 'resonancesWithBackground cross section has not been reconstructed via reactionSuite.reconstructResonances' )

class reference( baseCrossSectionForm ) :
    """This cross section from contain a reference to the another cross section."""

    form = tokens.referenceFormToken

    def __init__( self, crossSection=None ) :

        baseCrossSectionForm.__init__( self )
        self.setReference( crossSection )

    def getIncidentEnergyUnit( self ) :

        return( self.crossSection.getIncidentEnergyUnit( ) )

    def domainMin( self, unitTo = None, asPQU = False ) :

        return( self.crossSection.domainMin( unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        return( self.crossSection.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getDomain( self, unitTo = None, asPQU = False ) :

        return( self.crossSection.getDomain( unitTo = unitTo, asPQU = asPQU ) )

    def getDomainUnit( self ) :

        return( self.crossSection.getDomainUnit( ) )

    def getReference( self ) :

        return( self.crossSection )

    def setReference( self, crossSection ) :

        if( not( isinstance( crossSection, (component, None.__class__) ) ) ) :
            raise Exception( 'crossSection argument must be a cross section component not type %s' % type( crossSection ) )
        self.crossSection = crossSection

    def toPointwise_withLinearXYs( self, lowerEps, upperEps ) :

        return( self.crossSection.toPointwise_withLinearXYs( lowerEps, upperEps ) )

    def check( self ) :

        from fudge.gnd import warning
        warnings = []
        if self.getRootParent() != self.getReference().getRootParent():
            warnings.append( warning.badCrossSectionReference() )
        return warnings

    def toXMLList( self, indent = "" ) :

        xmlList = [ '%s<%s xlink:type="simple" xlink:href="%s"/>' % ( indent, self.form, self.crossSection.toXLink( ) ) ]
        return( xmlList )

    @staticmethod
    def parseXMLNode( element, xPath = [], linkData = {} ):

        from fudge.gnd import link

        xlink = link.parseXMLNode( element ).path
        ref = reference(None)
        if( 'unresolvedLinks' in linkData ) : linkData['unresolvedLinks'].append( ( xlink, ref ) )
        return( ref )

class weightedPointwise( baseCrossSectionForm ) :
    """This class is used to store one cross section as a weighted multiple of another cross section.
    This is typically used to express the 'production' cross section for radioactive products of a reaction.
    The production cross section can be given as an energy-dependent multiple of the reaction cross section.

    The weightedPointwise class contains a link to another cross section, and also contains a set of
    energy-dependent weights that should be multiplied by the other cross section."""

    form = tokens.weightedPointwiseFormToken

    def __init__( self, crossSection, weights ) :

        baseCrossSectionForm.__init__( self )
        self.setReference( crossSection )
        self.weights = weights

    def getIncidentEnergyUnit( self ) :

        return( self.crossSection.getIncidentEnergyUnit( ) )

    def domainMin( self, unitTo = None, asPQU = False ) :

        return( self.weights.domainMin( unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        return( self.weights.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getDomain( self, unitTo = None, asPQU = False ) :

        return( self.weights.getDomain( unitTo = unitTo, asPQU = asPQU ) )

    def getDomainUnit( self ) :

        return( self.crossSection.getDomainUnit( ) )

    def setReference( self, crossSection ):

        self.crossSection = crossSection

    def toPointwise_withLinearXYs( self, lowerEps, upperEps ) :

        xys = self.crossSection.toPointwise_withLinearXYs( lowerEps, upperEps ) * self.weights
        return( linear( xys.axes, xys, accuracy = xys.getAccuracy( ) ) )

    def toXMLList( self, indent = "" ) :

        xmlList = [ '%s<%s xlink:type="simple" xlink:href="%s">' % ( indent, self.form, self.crossSection.toXLink( ) ) ]
        xmlList += self.weights.toXMLList( indent = indent + '  ', incrementalIndent = '  ', pairsPerLine = 100 )
        xmlList[-1] += '</%s>' % self.form
        return( xmlList )

    @staticmethod
    def parseXMLNode( WPelement, xPath=[], linkData={} ):
        from fudge.gnd import link
        xlink = link.parseXMLNode(WPelement).path
        weights = XYs.XYs.parseXMLNode( WPelement[0] )
        ref = weightedPointwise( None, weights )
        if 'unresolvedLinks' in linkData: linkData['unresolvedLinks'].append((xlink, ref))
        return ref

    def toENDF6Data( self, MT, endfMFList, targetInfo, level ) :

        from fudge.legacy.converting import gndToENDF6
        endfInterpolation = gndToENDF6.axesToEndfInterpolationFlag( self.weights.axes )
        crossSectionFlatData = []
        for xy in self.weights.copyDataToXYs() : crossSectionFlatData += xy
        return( [ len( crossSectionFlatData ) / 2, endfInterpolation ], crossSectionFlatData )

def chargeParticle_changeInterpolationSubFunction( self, xInterpolation, yInterpolation, accuracy, lowerEps = 0, upperEps = 0 ) :

    def chargeParticle_changeInterpolationSubFunction2( x1, y1, x2, y2, limit ) :

        if( limit > limitMax ) : return
        if( y1 == 0. ) : return
        B = math.log( x2 * y2 / ( x1 * y1 ) ) / ( 1. / math.sqrt( x1 - T ) - 1. / math.sqrt( x2 - T ) )
        A = x1 * y1 * math.exp( B / math.sqrt( x1 - T ) )
        x = 0.5 * ( x1 + x2 )
        y = A * math.exp( - B / math.sqrt( x - T ) ) / x
        ylinlin = ( y1 * ( x2 - x ) + y2 * ( x - x1 ) ) / ( x2 - x1 )
        if( abs( y - ylinlin ) > accuracy * ( abs( y1 ) + abs( y2 ) ) ) :
            chargeParticle_changeInterpolationSubFunction2( x1, y1, x, y, limit + 1 )
            data.append( [ x, y ] )
            chargeParticle_changeInterpolationSubFunction2( x, y, x2, y2, limit + 1 )

    if( xInterpolation != axes.linearToken ) : raise Exception( 'Only linear interpolation currently supported for x-axis, not %s' % xInterpolation )
    if( yInterpolation != axes.linearToken ) : raise Exception( 'Only linear interpolation currently supported for y-axis, not %s' % yInterpolation )
    limitMax = self.getBiSectionMax( )
    T = self.findAttributeInAncestry( 'getThreshold' )( self.axes[0].getUnit( ) )
    x1, y1, data = None, 0.0, []
    for x2, y2 in self :
        if( x1 is not None ) :
            chargeParticle_changeInterpolationSubFunction2( x1, y1, x2, y2, 0 )
        data.append( [ x2, y2 ] )
        x1, y1 = x2, y2
    return( data )

def parseXMLNode( crossSectionElement, xPath=[], linkData={} ):
    """Reads an xml <crossSection> element into fudge, including all cross section forms (pointwise, linear, piecewise etc)
    contained inside the crossSection. """
    xPath.append( crossSectionElement.tag )
    xsc = component( nativeData=crossSectionElement.get('nativeData') )
    for form in crossSectionElement:
        formClass = {tokens.pointwiseFormToken: pointwise,
                tokens.linearFormToken: linear,
                tokens.piecewiseFormToken: piecewise,
                tokens.resonancesWithBackgroundFormToken: resonancesWithBackground,
                tokens.weightedPointwiseFormToken: weightedPointwise,
                tokens.groupedFormToken: grouped,
                tokens.referenceFormToken: reference,
                }.get( form.tag )
        if formClass is None: raise Exception("encountered unknown cross section form: %s" % form.tag)
        newForm = formClass.parseXMLNode( form, xPath, linkData )
        xsc.addForm( newForm )
    xPath.pop()
    return xsc
