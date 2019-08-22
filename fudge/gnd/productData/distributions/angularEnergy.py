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

""" energy/angular double differential distribution classes """

import math
import base, miscellaneous
import fudge
from fudge.core.math.xData import axes, XYs, W_XYs, V_W_XYs

__metaclass__ = type

class component( base.component ) :

    def __init__( self, nativeData = base.noneFormToken ) :

        base.component.__init__( self, base.angularEnergyComponentToken, nativeData )

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

class form( base.form ) :
    """Abstract base class for angularEnergy forms."""

    pass

class pointwise( form, V_W_XYs.V_W_XYs ) :

    tag = base.pointwiseFormToken
    moniker = base.pointwiseFormToken

    def __init__( self, axes, productFrame ) :

        form.__init__( self, base.pointwiseFormToken, productFrame )
        V_W_XYs.V_W_XYs.__init__( self, axes )

    def check( self, info ) :
        """ check for incomplete angular distributions + any negative probabilities.
        Ignore normalization for this double-differential distribution """

        from fudge.gnd import warning
        from pqu.physicalQuantityWithUncertainty import PhysicalQuantityWithUncertainty as PQU
        warnings = []

        for energy_in in self:
            integral = energy_in.integrate()
            if abs(integral - 1.0) > info['normTolerance']:
                warnings.append( warning.unnormalizedDistribution( PQU(energy_in.value, self.axes[0].unit),
                    energy_in.index, integral, self.toXLink() ) )
            if energy_in.domainMin() != -1 or energy_in.domainMax() != 1 :
                warnings.append( warning.incompleteDistribution( PQU(energy_in.value, self.axis[0].unit),
                        energy_in.domainMin(), energy_in.domainMax(), energy_in.value, self.toXLink() ) )
            for mu in energy_in:
                if mu.domainMin() < 0:
                    warnings.append( warning.valueOutOfRange("Negative outgoing energy for energy_in=%s!"
                        % PQU(energy_in.value, self.axes[0].unit), mu.domainMin(), 0, 'inf', self.toXLink() ) )
                if mu.yMin() < 0:
                    warnings.append( warning.negativeProbability( PQU(energy_in.value, self.axes[0].unit),
                        mu=mu.value, obj=mu ) )

        return warnings

    def calculateDepositionData( self, processInfo, tempInfo ) :

        if( self.getProductFrame( ) == axes.centerOfMassToken ) : raise Exception( 'center of mass calculation not supported for %s' % self.moniker )

        energyUnit = self.axes[2].getUnit( )
        momentumDepositionUnit = energyUnit + '/c'
        massUnit = energyUnit + '/c**2'
        energyAccuracy, momentumAccuracy = 1e-6, 1e-3
        multiplicity = tempInfo['multiplicity']
        productMass = tempInfo['product'].getMass( massUnit )

        depEnergy = []
        for E_MuEpPs in self :
            muEp = [ [ mu_EpPs.value, mu_EpPs.integrateWithWeight_x( ) ] for mu_EpPs in E_MuEpPs ]
            Ep = multiplicity.getValue( E_MuEpPs.value ) * XYs.XYs( axes.defaultAxes( ), muEp, energyAccuracy ).integrate( )
            depEnergy.append( [ E_MuEpPs.value, Ep ] )

        const = math.sqrt( 2. * productMass )
        depMomentum = []
        for E_MuEpPs in self :
            muEp = [ [ mu_EpPs.value, mu_EpPs.integrateWithWeight_sqrt_x( ) ] for mu_EpPs in E_MuEpPs ]
            pp = const * multiplicity.getValue( E_MuEpPs.value ) * XYs.XYs( axes.defaultAxes( ), muEp, energyAccuracy ).integrate( )
            if( pp < 1e-12 ) : pp = 0.                          # This should be less arbitrary????????
            depMomentum.append( [ E_MuEpPs.value, pp ] )

        axes_ = fudge.gnd.productData.energyDeposition.pointwise.defaultAxes( energyUnit = energyUnit, energyDepositionUnit = energyUnit )
        depEnergy = fudge.gnd.productData.energyDeposition.pointwise( axes_, depEnergy, energyAccuracy )
        axes_ = fudge.gnd.productData.momentumDeposition.pointwise.defaultAxes( energyUnit = energyUnit, momentumDepositionUnit = momentumDepositionUnit )
        depMomentum = fudge.gnd.productData.momentumDeposition.pointwise( axes_, depMomentum, momentumAccuracy )

        return( [ depEnergy, depMomentum ] )

    def extraXMLAttributeString( self ) :

        return( 'productFrame="%s"' % self.productFrame )

    def normalize( self, insitu = True ) :

        n = self
        if( not( insitu ) ) : n = self.copy( )
        for E_MuEpPs in n :
            sum = E_MuEpPs.integrate( )
            for muEpPs in E_MuEpPs : muEpPs.setData( muEpPs / sum )
        return( n )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        from fudge.processing.deterministic import transferMatrices

        newComponents = []

        if( 'LLNL_MC' in processInfo['styles'] ) :
            if( not( self.axes.isLinear( qualifierOk = True ) ) ) : raise Exception( 'Not implemented: %s' % self.moniker )

        if( 'LLNL_Pn' in processInfo['styles'] ) :
            if( processInfo['verbosity'] >= 30 ) : print '%sGrouping %s' % ( verbosityIndent, self.moniker )
            outputChannel = tempInfo['outputChannel'].outputChannel
            projectile, product = processInfo.getProjectileName( ), tempInfo['product'].particle.getToken( )
            TM_1, TM_E = transferMatrices.ENDFEMuEpP_TransferMatrix( processInfo, projectile, product, tempInfo['masses'], tempInfo['crossSection'], self, 
                tempInfo['multiplicity'], comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )
            fudge.gnd.miscellaneous.TMs2Form( processInfo, tempInfo, newComponents, TM_1, TM_E, self.axes )

        return( newComponents )

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        return( V_W_XYs.V_W_XYs.toPointwise_withLinearXYs( self, accuracy, lowerEps = lowerEps, upperEps = upperEps, cls = pointwise ) )

    def toENDF6( self, flags, tempInfo ) :

        from fudge.legacy.converting import gndToENDF6, endfFormats
        MF6 = [ endfFormats.endfContLine( 0, 0, 0, 0, 1, len( self ) ) ]
        EinInterp = gndToENDF6.axisToEndfInterpolationFlag( self.axes[0] )
        MF6 += endfFormats.endfInterpolationList( [ len( self ), EinInterp ] )
        muInterpolation = gndToENDF6.axisToEndfInterpolationFlag( self.axes[1] )
        independent, dependent, qualifier = self.axes[2].interpolation.getInterpolationTokens( )
        pdf_of_EpInterpolation = gndToENDF6.gndToEndfInterpolationFlag( independent, dependent )
        for oneEin in self :
            Ein = oneEin.value
            numMu = len( oneEin )
            MF6 += [ endfFormats.endfContLine( 0, Ein, 0, 0, 1, numMu ) ]
            MF6 += endfFormats.endfInterpolationList( [ numMu, muInterpolation ] )
            for entries in oneEin :
                mu = entries.value
                numEout = len( entries )
                MF6 += [ endfFormats.endfContLine( 0, mu, 0, 0, 1, numEout ) ]
                MF6 += endfFormats.endfInterpolationList( [ numEout, pdf_of_EpInterpolation ] )
                xys = entries.copyDataToXYs( xUnit = 'eV', yUnit = '1/eV' )
                MF6 += endfFormats.endfNdDataList( xys )
        return( 7, axes.labToken, MF6 )

    @staticmethod
    def defaultAxes( energyUnit = 'eV', energyInterpolation = axes.linearToken, energyFunctionInterpolation = axes.linearToken, 
        energyInterpolationQualifier = None, muInterpolation = axes.linearToken, energy_outInterpolation = axes.linearToken, 
        energy_outUnit = 'eV', probabilityInterpolation = axes.linearToken, probabilityUnit = '1/eV' ) :

        axes_ = axes.axes( dimension = 4 )
        axes_[0] = axes.axis( 'energy_in', 0, energyUnit, interpolation = axes.interpolationXY( energyInterpolation, energyFunctionInterpolation, energyInterpolationQualifier ) )
        axes_[1] = axes.axis( 'mu', 1, '', interpolation = axes.interpolationXY( muInterpolation, probabilityInterpolation ) )
        axes_[2] = axes.axis( 'energy_out', 2, energy_outUnit, interpolation = axes.interpolationXY( energy_outInterpolation, probabilityInterpolation ) )
        axes_[3] = axes.axis( 'P(mu,energy_out|energy_in)', 3, probabilityUnit )
        return( axes_ )

    @staticmethod
    def parseXMLNode( element, xPath=[], linkData={} ):

        xPath.append( element.tag )
        axes_ = axes.parseXMLNode( element[0], xPath )
        pw = pointwise( axes_, element.get( 'productFrame' ) )
        for energy_in in element[1:]:
            w_xys = W_XYs.W_XYs(
                    axes.referenceAxes(pw, 3), index=int(energy_in.get('index')),
                    value=float(energy_in.get("value")), parent=pw)
            for mu in energy_in:
                data = map(float,mu.text.split())
                data = zip( data[::2], data[1::2] )
                xys = XYs.XYs(
                        axes.referenceAxes(w_xys), data, float(mu.get('accuracy')),
                        value=float(mu.get('value')) )
                w_xys.append( xys )
            pw.append( w_xys )
        xPath.pop()
        return pw

class LLNLComponent( base.component ) :

    def __init__( self, nativeData = base.noneFormToken ) :

        base.component.__init__( self, base.LLNLAngularEnergyComponentToken, nativeData )

    @staticmethod
    def parseXMLNode( element, xPath=[], linkData={} ):

        xPath.append( element.tag )
        if element.get('nativeData') != base.pointwiseFormToken:
            raise Exception, "Only pointwise distribution currently handled inside LLNLAngularEnergy"
        pointwise_ = pointwise.parseXMLNode( element[0], xPath )
        component = LLNLComponent( pointwise_.moniker )
        component.addForm( pointwise_ )
        xPath.pop()
        return component

class LLNLEqualProbableBins( form ) :

    def __init__( self, productFrame ) :

        form.__init__( self, base.equalProbableBinsFormToken, productFrame )
        self.numberOfBins = len( data[0][1][0][1] )

    def toXMLList( self, indent = '' ) :
        """Returns the xml string representation of self."""

        indent2 = indent + '  '
        indent3 = indent2 + '  '
        xmlString = [ '%s<%s>' % ( indent, self.moniker ) ]
        xmlString[-1] += '<equalProbableBins3d bins="%d">' % ( self.numberOfBins )
        xmlString += self.axes.toXMLList( indent = indent2 )
        for indexE, energyMuEp in enumerate( self.data ) :
            energy, muEps = energyMuEp
            xmlString.append( '%s<energy value="%s" index="%d">' % ( indent2, fudge.gnd.miscellaneous.floatToString( energy ), indexE ) )
            for indexMu, muEp in enumerate( muEps ) :
                mu, Ep = muEp
                EpString = base.list1dToXMLEqualProbableBins1dString( Ep )
                xmlString.append( '%s<mu value="%s" index="%d">%s</energy>' % ( indent3, fudge.gnd.miscellaneous.floatToString( mu ), indexMu, EpString ) )
            xmlString[-1] += '</energy>'
        xmlString[-1] += '</equalProbableBins3d></%s>' % self.moniker
        return( xmlString )

class LLNL_withAngularComponent( base.component ) :
    """This class is for legacy LLNL ENDL I = 1 and 3 data and is deprecated. Only use for the legacy ENDL data."""

    def __init__( self, angularComponent, angularEnergyComponent ) :

        nativeData = '%s=%s:%s=%s' % ( base.angularComponentToken, angularComponent.nativeData, base.LLNLAngularEnergyComponentToken, 
            angularEnergyComponent.nativeData )
        base.component.__init__( self, base.LLNL_withAngularComponentToken, nativeData )
        self.angularComponent = angularComponent
        self.angularEnergyComponent = angularEnergyComponent

    def calculateDepositionData( self, processInfo, tempInfo ) :

        from fudge.core.math import fudgemath

        def calculateDepositionMomentumAtMu( mu, parameters ) :

            f = ( mu - parameters.mu1 ) / ( parameters.mu2 - parameters.mu1 )
            P_mu = ( 1 - f ) * parameters.P1 + f * parameters.P2
            EpP = parameters.muEpPs.interpolateAtW( mu, unitBase = True, extrapolation = W_XYs.flatExtrapolationToken )
            Ep = EpP.integrateWithWeight_sqrt_x( )
            return( mu * P_mu * Ep )

        def calculateDepositionMomentum( muPs, muEpPs ) :

            class LLNLAngular_angularEnergy_MomentumParameters :

                def __init__( self, mu1, P1, mu2, P2, muEpPs ) :

                    self.mu1, self.P1, self.mu2, self.P2, self.muEpPs = mu1, P1, mu2, P2, muEpPs

                def nextMu( self, m2, P2 ) :

                    self.mu1, self.P1 = self.mu2, self.P2
                    self.mu2, self.P2 = mu2, P2

            parameters = LLNLAngular_angularEnergy_MomentumParameters( 0, 0, 0, 0, muEpPs )
            mu1, p = None, 0
            for mu2, P2 in muPs :
                parameters.nextMu( mu2, P2 )
                if( mu1 is not None ) :
                    p += miscellaneous.GnG_adaptiveQuadrature( calculateDepositionMomentumAtMu, mu1, mu2, parameters, miscellaneous.GaussQuadrature2,
                        tolerance = 1e-4 )[0]
                mu1 = mu2
            return( p )

        class calculateDepositionMomentumThicken :

            def __init__( self, angularForm, angularEnergyForm, relativeTolerance, absoluteTolerance ) :

                self.angularForm = angularForm
                self.angularEnergyForm = angularEnergyForm
                self.relativeTolerance = relativeTolerance
                self.absoluteTolerance = absoluteTolerance

            def evaluateAtX( self, E ) :

                muPs = self.angularForm.getAtEnergy( E )
                muEpPs = self.angularEnergyForm.interpolateAtV( E, unitBase = True )
                return( calculateDepositionMomentum( muPs, muEpPs ) )

        depData = []
        angularForm = self.angularComponent.forms[self.angularComponent.nativeData]
        angularEnergyForm = self.angularEnergyComponent.forms[self.angularEnergyComponent.nativeData]
        multiplicity = tempInfo['multiplicity']
        energyUnit = angularEnergyForm.axes[0].getUnit( )
        energyPrimeUnit = angularEnergyForm.axes[2].getUnit( )
        momentumDepositionUnit = energyPrimeUnit + '/c'
        energyAccuracy, momentumAccuracy = 1e-6, 1e-3

        depData.append( miscellaneous.calculateDepositionEnergyFromAngular_angularEnergy( angularForm, angularEnergyForm, multiplicity, 
            accuracy = energyAccuracy ) )

        if( tempInfo['product'].getName( ) == 'gamma' ) :
            depMomentum = miscellaneous.calculateDepositionEnergyFromAngular_angularEnergy( angularForm, angularEnergyForm, multiplicity, 
                True, accuracy = momentumAccuracy )
        else :
            relativeTolerance, absoluteTolerance = momentumAccuracy, momentumAccuracy
            depMomentum = []
            for indexE, muPs in enumerate( angularForm ) :
                depMomentum.append( [ muPs.value, calculateDepositionMomentum( muPs, angularEnergyForm[indexE] ) ] )
#            depMomentum = fudgemath.thickenXYList( depMomentum, calculateDepositionMomentumThicken( angularForm, angularEnergyForm, relativeTolerance, absoluteTolerance ) )
            const = math.sqrt( 2. * tempInfo['product'].getMass( energyPrimeUnit + '/c**2' ) )
            for EMomenutem in depMomentum : EMomenutem[1] *= const * multiplicity.getValue( EMomenutem[0] )
        axes_ = fudge.gnd.productData.momentumDeposition.pointwise.defaultAxes( energyUnit = energyUnit, momentumDepositionUnit = momentumDepositionUnit )
        depData.append( fudge.gnd.productData.momentumDeposition.pointwise( axes_, depMomentum, momentumAccuracy ) )

        return( depData )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        import angular 
        from fudge.processing.deterministic import transferMatrices

        newComponents = []
        angularForm = self.angularComponent.forms[self.angularComponent.nativeData]
        angularEnergyForm = self.angularEnergyComponent.forms[self.angularEnergyComponent.nativeData]
        if( 'LLNL_MC' in processInfo['styles'] ) :
            if( processInfo['verbosity'] >= 30 ) : print '%sEqual probably binning %s' % ( verbosityIndent, self.moniker )
            raise Exception( 'Not implemented: %s' % self.moniker )

        if( 'LLNL_Pn' in processInfo['styles'] ) :
            if( processInfo['verbosity'] >= 30 ) : print '%sGrouping %s' % ( verbosityIndent, self.moniker )
            outputChannel = tempInfo['outputChannel'].outputChannel
            projectile, product = processInfo.getProjectileName( ), tempInfo['product'].particle.getToken( )
            TM_1, TM_E = transferMatrices.EMuEpP_TransferMatrix( processInfo, projectile, product, tempInfo['masses'],
                tempInfo['crossSection'], angularForm, angularEnergyForm, tempInfo['multiplicity'],
                comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )
            fudge.gnd.miscellaneous.TMs2Form( processInfo, tempInfo, newComponents, TM_1, TM_E, angularForm.axes )

        return( newComponents )

    def toXMLList( self, indent = "" ) : 

        indent2 = indent + '  '
        xmlString = [ '%s<%s>' % ( indent, self.moniker ) ]
        xmlString += self.angularComponent.toXMLList( indent = indent2 )
        xmlString += self.angularEnergyComponent.toXMLList( indent = indent2 )
        xmlString[-1] += '</%s>' % self.moniker 
        return( xmlString )

    def toENDF6( self, MT, endfMFList, flags, tempInfo ) :

        from fudge.legacy.converting import gndToENDF6
        from pqu.physicalQuantityWithUncertainty import PhysicalQuantityWithUncertainty as PQU

        angularForm = self.angularComponent.forms[self.angularComponent.nativeData]
        angularEnergyForm = self.angularEnergyComponent.forms[self.angularEnergyComponent.nativeData]

        energy_inInterpolation, energy_inFunctionInterpolation, energy_inInterpolationQualifier = angularEnergyForm.axes[0].interpolation.getInterpolationTokens( )
        muInterpolation, muFunctionInterpolation, muQualifier = angularEnergyForm.axes[1].interpolation.getInterpolationTokens( )
        energy_outInterpolation, probabilityInterpolation, energy_outQualifier = angularEnergyForm.axes[2].interpolation.getInterpolationTokens( )
        frame = angularEnergyForm.getProductFrame( )
        axes_ = pointwise.defaultAxes( energyInterpolation = energy_inInterpolation, energyFunctionInterpolation = energy_inFunctionInterpolation, 
                energyInterpolationQualifier = energy_inInterpolationQualifier, muInterpolation = muInterpolation, 
                energy_outInterpolation = energy_outInterpolation, probabilityInterpolation = probabilityInterpolation )
        E_inRatio = PQU( 1, angularEnergyForm.axes[0].getUnit( ) ).getValueAs( 'eV' )
        E_outRatio = PQU( 1, angularEnergyForm.axes[2].getUnit( ) ).getValueAs( 'eV' )
        LAW7 = pointwise( axes_, self.angularComponent.getProductFrame( ) )
        axesW_XY = axes.referenceAxes( LAW7, dimension = 3 )
        axesXY = axes.referenceAxes( LAW7 )

        if( len( angularForm ) != len( angularEnergyForm ) ) :
            raise Exception( "len( angularForm ) = %s != len( angularEnergyForm ) = %s" % ( len( angularForm ), len( angularEnergyForm ) ) )
        for indexE, EMuP in enumerate( angularForm ) :
            EMuEpP = angularEnergyForm[indexE]
            if( EMuP.value != EMuEpP.value ) : raise Exception( "At indexE = %d, EMuP.value %s != EMuEpP.value = %s" % ( indexE, EMuP.value, EMuEpP.value ) )
            if( len( EMuP ) != len( EMuEpP ) ) :
                raise Exception( "At indexE = %d (E_in = %s), len( EMuP ) %s != len( EMuEpP ) = %s" % ( indexE, EMuP.value, len( EMuP ), len( EMuEpP ) ) )
            w_xys = W_XYs.W_XYs( axesW_XY, index = indexE, value = EMuP.value )
            for indexMu, muP in enumerate( EMuP ) :
                muEpP = EMuEpP[indexMu]
                if( muP[0] != muEpP.value ) : raise Exception( "At indexE = %d, mu = %s != muEpP.value = %s" % ( indexE, muP[0], muEpP.value ) )
                xys = [ [ E_outRatio * Ep, muP[1] * P / E_outRatio ] for Ep, P in muEpP ]
                xys = XYs.XYs( axesXY, xys, accuracy = muEpP.getAccuracy( ), value = muP[0], index = indexMu, parent = w_xys )
                w_xys.append( XYs.XYs( axesXY, muEpP * muP[1], accuracy = muEpP.getAccuracy( ), value = muP[0], index = indexMu, parent = w_xys ) )
            LAW7.append( w_xys )

        LAW, frame, MF6 = LAW7.toENDF6( flags, tempInfo )
        gndToENDF6.toENDF6_MF6( MT, endfMFList, flags, tempInfo, LAW, frame, MF6 )

    @staticmethod
    def parseXMLNode( element, xPath=[], linkData={} ):

        xPath.append( element.tag )
        angularComponent = fudge.gnd.productData.distributions.angular.parseXMLNode( element[0], xPath )
        angularEnergyComponent = LLNLComponent.parseXMLNode( element[1], xPath )
        LwAC = LLNL_withAngularComponent( angularComponent, angularEnergyComponent )
        xPath.pop()
        return LwAC

def parseXMLNode( angularEnergyElement, xPath=[], linkData={} ):

    xPath.append( angularEnergyElement.tag )
    angularEnergy = component( angularEnergyElement.get("nativeData") )
    for formElement in angularEnergyElement:
        if formElement.tag==base.pointwiseFormToken:
            angularEnergy.addForm( pointwise.parseXMLNode( formElement, xPath, linkData ) )
        else: raise Exception("encountered unknown angularEnergy form: %s" % formElement.tag )
    xPath.pop()
    return angularEnergy

