# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""Energy/angular double differential distribution classes."""

import math

import xData.axes as axesModule
import xData.standards as standardsModule
import xData.series1d as series1dModule
import xData.XYs as XYsModule
import xData.multiD_XYs as multiD_XYsModule
import xData.regions as regionsModule

from PoPs import IDs as PoPsIDsModule

from . import base as baseModule
from . import energy as energyModule
from . import miscellaneous as miscellaneousModule

__metaclass__ = type

def defaultAxes( energyUnit ) :

    axes = axesModule.axes( rank = 4 )
    axes[0] = axesModule.axis( 'P(energy_out,mu|energy_in)', 0, '1/' + energyUnit )
    axes[1] = axesModule.axis( 'mu', 1, '' )
    axes[2] = axesModule.axis( 'energy_out', 2, energyUnit )
    axes[3] = axesModule.axis( 'energy_in', 3, energyUnit )
    return( axes )

class XYs1d( XYsModule.XYs1d ) :    # FIXME: BRB, Should this class inherit from angular.XYs1d?

    def averageMu( self ) :

        allowedInterpolations = [ standardsModule.interpolation.linlinToken,
                                  standardsModule.interpolation.flatToken ]
        xys = self.changeInterpolationIfNeeded( allowedInterpolations, XYsModule.defaultAccuracy )
        return( xys.integrateWithWeight_x( ) )

    def toLinearXYsClass( self ) :

        return( XYs1d )

class Legendre( series1dModule.LegendreSeries ) :

    def averageMu( self ) :

        return( self.getCoefficientSafely( 1 ) )

    def toLinearXYsClass( self ) :

        return( XYs1d )

class XYs2d( multiD_XYsModule.XYs2d ) :

    def averageEnergy( self ) :

        integralOfMu = [ [ pdfOfMuAtEp.outerDomainValue, float( pdfOfMuAtEp.integrate( ) ) ] for pdfOfMuAtEp in self ]
        return( float( XYsModule.XYs1d( integralOfMu, interpolation = self.interpolation ).integrateWithWeight_x( ) ) )

    def averageForwardMomentum( self, sqrt_2_ProductMass ) :

        averageMu = [ [ pdfOfMuAtEp.outerDomainValue, pdfOfMuAtEp.averageMu( ) ] for pdfOfMuAtEp in self ]
        return( sqrt_2_ProductMass * float( XYsModule.XYs1d( averageMu, interpolation = self.interpolation ).integrateWithWeight_sqrt_x( ) ) )

    def toLinearXYsClass( self ) :

        return( XYs2d )

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs1d, Legendre ) )

class subform( baseModule.subform ) :
    """Abstract base class for energyAngular subforms."""

    pass

class XYs3d( subform, multiD_XYsModule.XYs3d ) :

    def __init__( self, **kwargs ) :

        multiD_XYsModule.XYs3d.__init__( self, **kwargs )
        subform.__init__( self )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        verbosity = kwargs.get( 'verbosity', 0 )
        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        multiplicity = kwargs['multiplicity']
        projectileMass = kwargs['projectileMass']
        targetMass = kwargs['targetMass']
        productMass = kwargs['productMass']
        energyAccuracy = kwargs['energyAccuracy']
        momentumAccuracy = kwargs['momentumAccuracy']
        productFrame = kwargs['productFrame']

        sqrt_2_ProductMass = math.sqrt( 2 * productMass )
#
# FIXME, still need to fill in between energy points.
#
        aveEnergy = []
        aveMomentum = []
        if( productFrame == standardsModule.frames.labToken ) :
            for pdfOfEpMuAtE in self :
                energy = pdfOfEpMuAtE.outerDomainValue
                aveEnergy.append( [ energy, multiplicity.evaluate( energy ) * pdfOfEpMuAtE.averageEnergy( ) ] )

            for pdfOfEpMuAtE in self :
                energy = pdfOfEpMuAtE.outerDomainValue
                momemtum = multiplicity.evaluate( energy ) * pdfOfEpMuAtE.averageForwardMomentum( sqrt_2_ProductMass )
                if( momemtum < 1e-12 ) : momemtum = 0.          # This should be less arbitrary????????
                aveMomentum.append( [ energy, momemtum ] )

        else :                              # Center-of-mass.
            const1 = math.sqrt( 2 * projectileMass ) / ( projectileMass + targetMass )
            for pdfOfEpMuAtE in self :
                energy = pdfOfEpMuAtE.outerDomainValue
                vCOM = const1 * math.sqrt( energy )
                EpCOM = pdfOfEpMuAtE.averageEnergy( )
                ppCOM = pdfOfEpMuAtE.averageForwardMomentum( sqrt_2_ProductMass )
                productLabEnergy = 0.5 * productMass * vCOM * vCOM + EpCOM + vCOM * ppCOM
                aveEnergy.append( [ energy, multiplicity.evaluate( energy ) * productLabEnergy ] )

            for pdfOfEpMuAtE in self :
                energy = pdfOfEpMuAtE.outerDomainValue
                vCOM = const1 * math.sqrt( energy )
                productLabForwardMomentum = productMass * vCOM + pdfOfEpMuAtE.averageForwardMomentum( sqrt_2_ProductMass )
                aveMomentum.append( [ energy, multiplicity.evaluate( energy ) * productLabForwardMomentum ] )

        return( [ aveEnergy ], [ aveMomentum ] )

    def check( self, info ) :
        from fudge import warning
        from pqu import PQU
        warnings = []
        axes = axesModule.axes()
        axes[0] = self.axes[-2]
        axes[1] = self.axes[-1]

        for index, energy_in in enumerate(self):
            if isinstance(energy_in[0], Legendre):
                integral = float( XYsModule.XYs1d( [ ( eout.outerDomainValue, eout.coefficients[0] ) for eout in energy_in ], axes = axes ).integrate( ) )
            else:
                integral = float( XYsModule.XYs1d( [ ( eout.outerDomainValue, float(eout.integrate()) ) for eout in energy_in ], axes = axes ).integrate( ) )
            if abs(integral-1.0) > info['normTolerance']:
                warnings.append( warning.unnormalizedDistribution( PQU.PQU( energy_in.outerDomainValue, axes[0].unit ), index, integral, self.toXLink() ) )
            energy_outs = {}
            for xy in energy_in :
                energy_outs[xy.toPointwise_withLinearXYs( accuracy = XYsModule.defaultAccuracy, upperEps = 1e-8 ).rangeMin] = xy.outerDomainValue
            minProb = min( energy_outs.keys( ) )
            if minProb < 0:
                warnings.append( warning.negativeProbability( PQU.PQU( energy_in.outerDomainValue, axes[0].unit ), value = minProb, 
                        energy_out = PQU.PQU( energy_outs[minProb], axes[1].unit ), obj = energy_in ) )
            if self.interpolationQualifier is standardsModule.interpolation.unitBaseToken:
                # check for more than one outgoing distribution integrating to 0 at high end of incident energies
                integrals = [eout.integrate() for eout in energy_in]
                for idx, integral in enumerate(integrals[::-1]):
                    if integral != 0: break
                if idx > 1:
                    warnings.append( warning.extraOutgoingEnergy( PQU.PQU( energy_in.outerDomainValue, axes[-1].unit ), obj = energy_in ) )
        return warnings

    def energySpectrumAtEnergy( self, energyIn ) :

        energyUnit = self.domainUnit
        flag, functional2d1, functional2d2, frac, interpolation, interpolationQualifier = self.getBoundingSubFunctions( energyIn )
        if( flag in [ '<', '>', None ] ) : return( energyModule.XYs1d( axes = energyModule.defaultAxes( energyUnit ) ) )
        data = [ [ functional1d.outerDomainValue, float( functional1d.integrate( ) ) ] for functional1d in functional2d1 ]
        functional1d = energyModule.XYs1d( data, interpolation = functional2d1.interpolation, axes = energyModule.defaultAxes( energyUnit ) )
        if( flag != '=' ) :
            data = [ [ functional1d.outerDomainValue, float( functional1d.integrate( ) ) ] for functional1d in functional2d2 ]
            functional1d2 = energyModule.XYs1d( data, interpolation = functional2d2.interpolation, axes = energyModule.defaultAxes( energyUnit ) )
            frac = ( energyIn - functional2d1.outerDomainValue ) / ( functional2d2.outerDomainValue - functional2d1.outerDomainValue )
            if( False ) :
                functional1d = ( 1.0 - frac ) * functional1d + frac * functional1d2
            else :
                functional1d = XYsModule.pointwiseXY_C.unitbaseInterpolate( energyIn, functional2d1.outerDomainValue, functional1d.nf_pointwiseXY,
                        functional2d2.outerDomainValue, functional1d2.nf_pointwiseXY, 1 )
                functional1d = energyModule.XYs1d( functional1d, interpolation = functional1d.getInterpolation( ), axes = energyModule.defaultAxes( energyUnit ) )

        functional1d.normalize( insitu = True )

        return( functional1d )

    def normalize( self, insitu = True ) :

        n = self
        if( not( insitu ) ) : n = self.copy( )
        for E_MuEpPs in n :
            sum = float( E_MuEpPs.integrate( ) )
            for muEpPs in E_MuEpPs : muEpPs.setData( muEpPs / sum )
        return( n )

    def processMultiGroup( self, style, tempInfo, indent ) :

        from fudge.processing import group as groupModule
        from fudge.processing.deterministic import transferMatrices as transferMatricesModule

        verbosity = tempInfo['verbosity']
        productFrame = tempInfo['productFrame']

        if( verbosity > 2 ) : print('%sGrouping %s' % (indent, self.moniker))
        TM_1, TM_E = transferMatricesModule.EEpMuP_TransferMatrix( style, tempInfo, productFrame, tempInfo['crossSection'], self,
            tempInfo['multiplicity'], comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )

        return( groupModule.TMs2Form( style, tempInfo, TM_1, TM_E ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        return( multiD_XYsModule.XYs3d.toPointwise_withLinearXYs( self, cls = XYs3d, **kwargs ) )

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :

        from . import energyAngularMC as energyAngularMCModule

        energy = energyModule.XYs2d( axes = self.axes, interpolation = self.interpolation, interpolationQualifier = self.interpolationQualifier )
        for PofMuGivenEp in self :
            data = energyModule.XYs1d( [ [ PofMu.outerDomainValue, float( PofMu.integrate( ) ) ] for PofMu in PofMuGivenEp ], interpolation = PofMuGivenEp.interpolation )
            energy.append( energyModule.xs_pdf_cdf1d.fromXYs( data, outerDomainValue = PofMuGivenEp.outerDomainValue ) )

        xys3d = energyAngularMCModule.XYs3d( axes = self.axes )
        for PofMuGivenEp in self :
            xys2d = energyAngularMCModule.XYs2d( outerDomainValue = PofMuGivenEp.outerDomainValue )
            for PofMu in PofMuGivenEp :
                _PofMu = PofMu
                if( isinstance( PofMu, Legendre ) ) :
                    if( PofMu[0] == 0 ) : _PofMu = XYs1d( [ [ -1, 0.5 ], [ 1, 0.5 ] ] )
                    _PofMu = _PofMu.toPointwise_withLinearXYs( accuracy = XYsModule.defaultAccuracy, upperEps = 1e-8 )
                xys2d.append( energyAngularMCModule.xs_pdf_cdf1d.fromXYs( _PofMu, PofMu.outerDomainValue ) )
            xys3d.append( xys2d )
        return( energyAngularMCModule.energy( energy ), energyAngularMCModule.energyAngular( xys3d ) )

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs2d, ) )

class regions3d( subform, regionsModule.regions3d ) :

    def __init__( self, **kwargs ) :

        regionsModule.regions3d.__init__( self, **kwargs )
        subform.__init__( self )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        raise Exception( 'Not implemented' )

    def check( self, info ) :

        raise Exception( 'Not implemented' )

    def normalize( self, insitu = True ) :

        n = self
        if( not( insitu ) ) : n = self.copy( )
        for region in n : region.normalize( insitu = True )
        return( n )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        raise Exception( 'Not implemented' )

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs3d, ) )

class form( baseModule.form ) :

    moniker = 'energyAngular'
    subformAttributes = ( 'energyAngularSubform', )
    ancestryMembers = subformAttributes

    def __init__( self, label, productFrame, energyAngularSubform ) :

        if( not( isinstance( energyAngularSubform, subform ) ) ) : raise TypeError( 'instance is not an energyAngular subform' )
        baseModule.form.__init__( self, label, productFrame, ( energyAngularSubform, ) )

    @property
    def domainMin( self ) :

        return( self.energyAngularSubform.domainMin )

    @property
    def domainMax( self ) :

        return( self.energyAngularSubform.domainMax )

    @property
    def domainUnit( self ) :

        return( self.energyAngularSubform.domainUnit )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        kwargs['productFrame'] = self.productFrame
        return( self.energyAngularSubform.calculateAverageProductData( style, indent = indent, **kwargs ) )

    def energySpectrumAtEnergy( self, energyIn, frame, **kwargs ) :

        if( self.product.id == PoPsIDsModule.photon ) : frame = self.productFrame       # Kludge for now.
        if( self.productFrame == frame ) :
            return( self.energyAngularSubform.energySpectrumAtEnergy( energyIn ) )
        else :
            if( self.productFrame == standardsModule.frames.labToken ) : raise TypeError( 'Lab to center-of-mass translation not supported.' )
            xys2d = self.spectrumAtEnergy( energyIn, standardsModule.frames.labToken )
            data = [ [ xys1d.outerDomainValue, float( xys1d.integrate( ) ) ] for xys1d in xys2d ]
            return( energyModule.XYs1d( data, axes = energyModule.defaultAxes( self.domainUnit ), interpolation = xys2d.interpolation ) )

    def integrate( self, reaction_suite, energyIn, energyOut = None, muOut = None, phiOut = None, frame = standardsModule.frames.productToken, LegendreOrder = 0 ) :

# Need to check frame.

        data = self.energyAngularSubform
        if( energyIn < data.domainMin ) : return( 0.0 )
        if( energyIn > data.domainMax ) : energyIn = self.domainMax

        XYs2d = data.evaluate( energyIn )
        points = []
        for function1d in XYs2d :
            muOutMin, muOutMax = miscellaneousModule.domainLimits( muOut, function1d.domainMin, function1d.domainMax )
            if( muOutMax is None ) :
                muPartialIntegral = function1d.evaluate( muOutMin )
            else :
                muPartialIntegral = float( function1d.integrate( muOutMin, muOutMax ) )
            points.append( [ function1d.outerDomainValue, float( muPartialIntegral ) ] )
        xys = XYs1d( data = points )

        energyOutMin, energyOutMax = miscellaneousModule.domainLimits( energyOut, xys.domainMin, xys.domainMax )
        if( energyOutMax is None ) :
            energyOutMuEvaluate = xys.evaluate( energyOutMin )
        else :
            energyOutMuEvaluate = xys.integrate( energyOutMin, energyOutMax )

        phiEvaluate = miscellaneousModule.muPhiEvaluate( None, phiOut )

        return( phiEvaluate * energyOutMuEvaluate )

    def processMC_cdf( self, style, tempInfo, indent ) :

        from . import energyAngularMC as energyAngularMCModule

        energy, energyAngular = self.energyAngularSubform.to_xs_pdf_cdf1d( style, tempInfo, indent )
        return( energyAngularMCModule.form( style.label, self.productFrame, energy, energyAngular ) )

    def processMultiGroup( self, style, tempInfo, indent ) :

        tempInfo['productFrame'] = self.productFrame
        return( self.energyAngularSubform.processMultiGroup( style, tempInfo, indent ) )

    def spectrum( self, frame, **kwargs ) :

        if( not( isinstance( self.energyAngularSubform, XYs3d ) ) ) : raise TypeError( 'Form "%s" is not supported' % self.energyAngularSubform.moniker )

        if( frame == self.productFrame ) : return( self.energyAngularSubform.toPointwise_withLinearXYs( accuracy = 1e-3 ) )
        if( frame == standardsModule.frames.centerOfMassToken ) : raise ValueError( 'Conversion from frame "%s" to "%s" not supported' % ( self.productFrame, frame ) )

        xys3d = XYs3d( axes = defaultAxes( self.domainUnit ) )
        for energy in self.energyAngularSubform.domainGrid : xys3d.append( self.spectrumAtEnergy( energy, frame ) )

        return( xys3d )

    def spectrumAtEnergy( self, energyIn, frame ) :

        class angualarAtEnergyCOM :

            def __init__( self, probability ) :

                self.probability = probability

            def probabilityCOM( self, energyPrimeCOM, muCOM ) :

                return( self.probability.evaluate( energyPrimeCOM ).evaluate( muCOM ) )

        _spectrumAtEnergy = self.energyAngularSubform.evaluate( energyIn ).toPointwiseLinear( lowerEps = 1e-7 ) # lowerEps should be user input and not hardwired.
        if( frame == self.productFrame ) : return( _spectrumAtEnergy )

        if( frame == standardsModule.frames.centerOfMassToken ) : raise ValueError( 'Conversion from frame "%s" to "%s" not supported' % ( self.productFrame, frame ) )

        energyProbability = [ [ probability.outerDomainValue, float( probability.integrate( ) ) ] for probability in _spectrumAtEnergy ]
        energyProbability = energyModule.XYs1d( energyProbability, axes = energyModule.defaultAxes( self.domainUnit ) )
        return( miscellaneousModule.energyAngularSpectrumFromCOMSpectrumToLabAtEnergy( self, energyIn, energyProbability, 
                angualarAtEnergyCOM( _spectrumAtEnergy ), angularIsNormalized = False ) )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):
        """Translate <energyAngular> element from xml."""

        xPath.append( element.tag )
        subformElement = element[0]
        subformClass = {
                XYs3d.moniker: XYs3d,
                }.get( subformElement.tag )
        if subformClass is None: raise Exception( "encountered unknown energyAngular subform: %s" % subformElement.tag )
        subForm = subformClass.parseXMLNode( subformElement, xPath, linkData )
        energyAngular = form( element.get( "label" ), element.get('productFrame'), subForm )
        xPath.pop()
        return energyAngular
