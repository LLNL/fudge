# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""Angular/energy double differential distribution classes."""

import math

from xData import enums as xDataEnumsModule
from xData import axes as axesModule
from xData import XYs1d as XYs1dModule
from xData import multiD_XYs as multiD_XYsModule

from fudge.processing import group as groupModule

from . import miscellaneous as miscellaneousModule

from . import base as baseModule
from . import angular as angularModule
from . import angularEnergyMC as angularEnergyMCModule


def defaultAxes( energyUnit, probabilityLabel = 'P(mu,energy_out|energy_in)' ) :

    axes = axesModule.Axes(4)
    axes[3] = axesModule.Axis( 'energy_in', 3, energyUnit )
    axes[2] = axesModule.Axis( 'mu', 2, '' )
    axes[1] = axesModule.Axis( 'energy_out', 1, energyUnit )
    axes[0] = axesModule.Axis( probabilityLabel, 0, '1/' + energyUnit )
    return( axes )

class XYs1d( XYs1dModule.XYs1d ) :

    def averageEnergy( self ) :

        return( self.integrateWithWeight_x( ) )
    
    def averageMomentum( self ) :

        return( self.integrateWithWeight_sqrt_x( ) )

    def toLinearXYsClass( self ) :

        return( XYs1d )
    
class XYs2d( multiD_XYsModule.XYs2d ) :

    def  __add__( self, other ) :

        functionNd = self.__class__( interpolation = self.interpolation, axes = self.axes.copy( ), index = self.index, outerDomainValue = self.outerDomainValue,
                label = self.label, interpolationQualifier = self.interpolationQualifier )

        if( isinstance( other, self.__class__ ) ) :
            x2s = list( set( [ function.outerDomainValue for function in self ] + [ function.outerDomainValue for function in other ] ) )
            for x2 in sorted( x2s ) :
                function = self.evaluate( x2 ) + other.evaluate( x2 )
                function.outerDomainValue = x2
                functionNd.append( function )
            raise NotImplementedError("Adding two angularEnergy.XYs2d instances")
        else :
            try :
                value = float( other )
            except :
                raise ValueError( 'Other must be a number or a subclass of "%s".' % other.__class__ )
            for function in self.functionals : functionNd.append( function + value )

        return( functionNd )

    __radd__ = __add__

    def averageEnergy( self ) :

        EpOfMu = [ [ pdfOfEpAtMu.outerDomainValue, pdfOfEpAtMu.averageEnergy( ) ] for pdfOfEpAtMu in self ]
        return XYs1dModule.XYs1d(EpOfMu).integrate()

    def averageMomentum( self ) :

        MpOfMu = [ [ pdfOfMpAtMu.outerDomainValue, pdfOfMpAtMu.averageMomentum( ) ] for pdfOfMpAtMu in self ]
        return XYs1dModule.XYs1d(MpOfMu).integrate()

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs1d, ) )

class Subform( baseModule.Subform ) :
    """Abstract base class for angularEnergy forms."""

    pass

class XYs3d( Subform, multiD_XYsModule.XYs3d ) :

    def __init__( self, **kwargs ) :

        multiD_XYsModule.XYs3d.__init__( self, **kwargs )
        Subform.__init__( self )

    def check( self, info ) :
        """
        Check for incomplete angular distributions + any negative probabilities.
        Ignore normalization for this double-differential distribution.
        """

        from fudge import warning
        from pqu import PQU
        warnings = []

        for index, energy_in in enumerate(self):
            integral = energy_in.integrate()
            if abs(integral - 1.0) > info['normTolerance']:
                warnings.append( warning.UnnormalizedDistribution( PQU.PQU( energy_in.outerDomainValue, self.axes[0].unit ), index, integral, self.toXLink() ) )
            if( energy_in.domainMin != -1 ) or ( energy_in.domainMax != 1 ) :
                warnings.append( warning.IncompleteDistribution( PQU.PQU( energy_in.outerDomainValue, self.axes[0].unit ), energy_in.domainMin, energy_in.domainMax, energy_in ) )
            for mu in energy_in:
                if( mu.domainMin < 0 ) :
                    warnings.append( warning.ValueOutOfRange("Negative outgoing energy for energy_in=%s!"
                        % PQU.PQU( energy_in.outerDomainValue, self.axes[0].unit ), mu.domainMin, 0, 'inf', self.toXLink() ) )
                if( mu.rangeMin < 0 ) :
                    warnings.append( warning.NegativeProbability( PQU.PQU( energy_in.outerDomainValue, self.axes[-1].unit ), mu=mu.outerDomainValue, obj=mu ) )

        return warnings

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        multiplicity = kwargs['multiplicity']
        productMass = kwargs['productMass']

        aveEnergy = []
        for pdfOfMuEpAtE in self :
            energy = pdfOfMuEpAtE.outerDomainValue
            aveEnergy.append( [ energy, multiplicity.evaluate( energy ) * pdfOfMuEpAtE.averageEnergy( ) ] )

        const = math.sqrt( 2. * productMass )
        aveMomentum = []
        for pdfOfMuEpAtE in self :
            energy = pdfOfMuEpAtE.outerDomainValue
            momemtum = const * multiplicity.evaluate( energy ) * pdfOfMuEpAtE.averageMomentum( )
            if( momemtum < 1e-12 ) : momemtum = 0.          # This should be less arbitrary????????
            aveMomentum.append( [ energy, momemtum ] )

        return( aveEnergy, aveMomentum )

    def evaluate( self, domainValue, epsilon = 0 ) :

        print( 'Hi' )
        position, function1, function2, frac, interpolation, interpolationQualifier = self.getBoundingSubFunctions( domainValue )
        print( type( function1 ), type( function2 ), position, frac, interpolationQualifier )
        if( position is None ) : raise Exception( "No data to interpolate" )

        if( frac <= epsilon ) :                 # If close to first point pick it.
            function = function1.copy( )
        elif( ( 1 - frac ) <= epsilon ) :       # If close to second point pick it.
            function = function2.copy( )
        else :
            if( position in ( '=', '<', '>' ) ) :
                if( position == '=' ) :
                    function = function1.copy( )
                elif( position == '>' ) :
                    function = function1.copy( )
                else :
                    index = { '<' : 0, '>' : -1 }[position]
                    raise Exception( "evaluation point = %s %s than %s" % ( domainValue, { '<' : 'less', '>' : 'greater' }[position], self[index].outerDomainValue ) )
            else :
                if( not( isinstance( function1, XYs1dModule.XYs1d ) ) ) :      # FIXME, accuracy, lowerEps and upperEps should not be hardwired.
                    function1 = function1.toPointwise_withLinearXYs( accuracy = 1e-4, lowerEps = 1e-6, upperEps = 1e-6 )
                if( not( isinstance( function2, XYs1dModule.XYs1d ) ) ) :
                    function2 = function2.toPointwise_withLinearXYs( accuracy = 1e-4, lowerEps = 1e-6, upperEps = 1e-6 )

                if( not( isinstance( function1, XYs2d ) ) ) : raise Exception( 'function1 must be an XYs2d instance' )
                if( not( isinstance( function2, XYs2d ) ) ) : raise Exception( 'function2 must be an XYs2d instance' )

                mus = sorted( set( [ function1d.outerDomainValue for function1d in function1 ] + [ function1d.outerDomainValue for function1d in function2 ] ) )
                function = XYs2d( axes = function1.axes )
                for mu in mus :
                    function1d1 = function1.evaluate( mu )
                    function1d2 = function2.evaluate( mu )
                    if( not( isinstance( function1d1, XYs1d ) ) or not( isinstance( function1d2, XYs1d ) ) ) :
                        raise Exception( 'function1d1 and function1d2 must be an XYs1d instance' )
                    xy = XYs1dModule.pointwiseXY_C.unitbaseInterpolate( domainValue, function1.outerDomainValue, function1d1, function2.outerDomainValue, function1d2, 1 )
                    xy = XYs1d( xy, outerDomainValue = mu )
                    function.append( xy )

        function.outerDomainValue = domainValue
        for ii in function : print( '   ii  ', ii.outerDomainValue, ii.integrate( ) )
        function.normalize( insitu = True )
        for ii in function : print( '   ii2 ', ii.outerDomainValue, ii.integrate( ) )
        return( function )

    def integrate( self, reaction_suite, energyIn, energyOut = None, muOut = None, phiOut = None, frame = xDataEnumsModule.Frame.product, LegendreOrder = 0 ) :

        if( energyIn < self.domainMin ) : return( 0.0 )
        if( energyIn > self.domainMax ) : energyIn = self.domainMax
        function2d = self.evaluate( energyIn )

        muMin, muMax = miscellaneousModule.domainLimits( muOut, function2d.domainMin, function2d.domainMax )
        if( muMax is None ) :
            function1d = function2d.evaluate( muMin )
            energyOutMin, energyOutMax = miscellaneousModule.domainLimits( energyOut, function1d.domainMin, function1d.domainMax )
            if( energyOutMax is None ) :
                muEnergyOutEvaluate = function1d.evaluate( energyOutMin )
            else :
                muEnergyOutEvaluate = function1d.integrate( energyOutMin, energyOutMax )
        else :
            xys = []
            for function1d in function2d :
                energyOutMin, energyOutMax = miscellaneousModule.domainLimits( energyOut, function1d.domainMin, function1d.domainMax )
                if( energyOutMax is None ) :
                    value = function1d.evaluate( energyOutMin )
                else :
                    value = function1d.integrate(energyOutMin, energyOutMax)
                xys.append( [ function1d.outerDomainValue, value ] )
            xys = XYs1d( xys )
            muEnergyOutEvaluate = xys.integrate( muMin, muMax )

        phiEvaluate = miscellaneousModule.muPhiEvaluate( None, phiOut )

        return( phiEvaluate * muEnergyOutEvaluate )

    def normalize( self, insitu = True ) :

        n = self
        if( not( insitu ) ) : n = self.copy( )
        for E_MuEpPs in n :
            sum = E_MuEpPs.integrate( )
            for muEpPs in E_MuEpPs : muEpPs.setData( muEpPs / sum )
        return( n )

    def processMultiGroup( self, style, tempInfo, indent ) :

        from fudge.processing.deterministic import transferMatrices as transferMatricesModule

        verbosity = tempInfo['verbosity']
        productFrame = tempInfo['productFrame']

        if( verbosity > 2 ) : print( '%sGrouping %s' % ( indent, self.moniker ) )
        TM_1, TM_E = transferMatricesModule.ENDFEMuEpP_TransferMatrix( style, tempInfo, productFrame, tempInfo['crossSection'], self,
            tempInfo['multiplicity'], comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )

        return( groupModule.TMs2Form( style, tempInfo, TM_1, TM_E ) )

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :

        angular = angularModule.XYs2d( axes = self.axes )
        for PofEpGivenMu in self :
            data = angularModule.XYs1d( [ [ PofEp.outerDomainValue, PofEp.integrate() ] for PofEp in PofEpGivenMu ] )
            angular.append( angularModule.Xs_pdf_cdf1d.fromXYs( angularModule.XYs1d( data ), outerDomainValue = PofEpGivenMu.outerDomainValue ) )

        xys3d = angularEnergyMCModule.XYs3d( axes = self.axes )
        for PofEpGivenMu in self :
            xys2d = angularEnergyMCModule.XYs2d( outerDomainValue = PofEpGivenMu.outerDomainValue )
            for PofEp in PofEpGivenMu :
                _PofEp = PofEp.toPointwise_withLinearXYs( accuracy = 1e-3, upperEps = 1e-8 )
                xys2d.append( angularEnergyMCModule.Xs_pdf_cdf1d.fromXYs( _PofEp, PofEp.outerDomainValue ) )
            xys3d.append( xys2d )

        return( angularEnergyMCModule.Angular( angular ), angularEnergyMCModule.AngularEnergy( xys3d ) )

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs2d, ) )

class Form(  baseModule.Form ) :

    moniker = 'angularEnergy'
    subformAttributes = ( 'angularEnergySubform', )
    ancestryMembers = subformAttributes

    def __init__( self, label, productFrame, angularEnergySubform ) :

        if( not( isinstance( angularEnergySubform, Subform ) ) ) : raise TypeError( 'instance is not an angularEnergy subform' )
        baseModule.Form.__init__( self, label, productFrame, ( angularEnergySubform, ) )

    @property
    def domainUnit( self ) :

        return( self.angularEnergySubform.domainUnit )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        if self.productFrame == xDataEnumsModule.Frame.centerOfMass:
            raise Exception( 'center of mass calculation not supported for %s' % self.moniker )

        kwargs['productFrame'] = self.productFrame
        aveEnergy, aveMomentum = self.angularEnergySubform.calculateAverageProductData( style, indent = indent, **kwargs )
        return( [ aveEnergy ], [ aveMomentum ] )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.angularEnergySubform.convertUnits( unitMap )

    def energySpectrumAtEnergy( self, energyIn, frame, **kwargs ) :

        if( self.productFrame == frame ) :
            self.angularEnergySubform.energySpectrumAtEnergy( energyIn )
            function2d = self.angularEnergySubform.evaluate( energyIn )
        else :
            raise TypeError( '%s to %s translation not supported.' % ( self.productFrame, frame ) )

    def fixDomains(self, energyMin, energyMax, domainToFix):
        """
        This method call \*\*fixDomains\* on the *angularEnergySubform* member.
        """

        return self.angularEnergySubform.fixDomains(energyMin, energyMax, domainToFix)
        
    def integrate( self, reaction_suite, energyIn, energyOut = None, muOut = None, phiOut = None, frame = xDataEnumsModule.Frame.product, LegendreOrder = 0 ) :

        return( self.angularEnergySubform.integrate( reaction_suite, energyIn, energyOut = energyOut, muOut = muOut, phiOut = phiOut, frame = frame, LegendreOrder = LegendreOrder ) )

    def processMC_cdf( self, style, tempInfo, indent ) :

        angular, angularEnergy = self.angularEnergySubform.to_xs_pdf_cdf1d( style, tempInfo, indent )
        return( angularEnergyMCModule.Form( style.label, self.productFrame, angular, angularEnergy ) )

    def processMultiGroup( self, style, tempInfo, indent ) :

        tempInfo['productFrame'] = self.productFrame
        return( self.angularEnergySubform.processMultiGroup( style, tempInfo, indent ) )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )
        formElement = element[0]
        subformClass = {
            XYs3d.moniker :  XYs3d,
        }[ formElement.tag ]
        if subformClass is None: raise Exception("encountered unknown angularEnergy subform: %s" % formElement.tag)
        subform = XYs3d.parseNodeUsingClass(formElement, xPath, linkData, **kwargs)
        AEC = cls( element.get( 'label' ), element.get('productFrame'), subform )
        xPath.pop()
        return AEC
