# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Photon coherent double difference cross section form and its supporting classes.
"""

import math

from pqu import PQU as PQUModule

import xData.ancestry as ancestryModule
import xData.standards as standardsModule
import xData.XYs as XYsModule

from fudge.core.math import fudgemath as fudgemathModule
from ... import crossSection as crossSectionModule
from .. import base as baseModule

__metaclass__ = type

class XYs1d( baseModule.XYs1d ) :

    pass

class regions1d( baseModule.regions1d ) :

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs1d, ) )

class coherentFunctionBase( ancestryModule.ancestry ) :

    def __init__( self, data ) :

        ancestryModule.ancestry.__init__( self )
        if( not( isinstance( data, ( baseModule.XYs1d, baseModule.regions1d ) ) ) ) :
            raise TypeError( "Needed %s baseModulesubform." % self.moniker )

        self.data = data
        self.data.setAncestor( self )

    def convertUnits( self, unitMap ):

        self.data.convertUnits( unitMap )

    def check( self, info ) :

        return( [] )

    def toXML( self, indent = "", **kwargs ) :

        return( '\n'.join( self.toXMLList( indent = indent, **kwargs ) ) )

    def toXMLList( self, indent = "", **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ "%s<%s>" % ( indent, self.moniker ) ]
        xmlStringList += self.data.toXMLList( indent = indent2, **kwargs )
        xmlStringList[-1] += "</%s>" % self.moniker
        return( xmlStringList )

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ) :

        xPath.append( element.tag )

        data = element[0]
        if( data.tag == baseModule.XYs1d.moniker ) :
            data = XYs1d.parseXMLNode( data, xPath, linkData )
        elif( data.tag == baseModule.regions1d.moniker ) :
            data = regions1d.parseXMLNode( data, xPath, linkData )
        else :
            raise TypeError( 'Invalid data "%s" for "%s"' % ( data.tag, cls.tag ) )

        xPath.pop( )
        return( cls( data ) )

class formFactor( coherentFunctionBase ) :

    moniker = 'formFactor'
    ENDFMT = 502

class imaginaryAnomalousFactor( coherentFunctionBase ) :

    moniker = 'imaginaryAnomalousFactor'
    ENDFMT = 505

class realAnomalousFactor( coherentFunctionBase ) :

    moniker = 'realAnomalousFactor'
    ENDFMT = 506

class form( baseModule.form ) :

    moniker = 'coherentPhotonScattering'
    subformAttributes = ( 'formFactor', 'anomalousScatteringFactor_realPart', 'anomalousScatteringFactor_imaginaryPart' )

    def __init__( self, pid, label, productFrame, _formFactor, realPart, imaginaryPart ) :

        if( not( isinstance( _formFactor, formFactor ) ) ) :
            raise Exception( 'Instance is class "%s" and not formFactor' % formFactor.__class__ )

        if( realPart is not None ) :
            if( not( isinstance( realPart, realAnomalousFactor ) ) ) :
                raise Exception( 'Instance is class "%s" and not realAnomalousFactor' % realPart.__class__ )

        if( imaginaryPart is not None ) :
            if( not( isinstance( imaginaryPart, imaginaryAnomalousFactor ) ) ) :
                raise Exception( 'Instance is class "%s" and not imaginaryAnomalousFactor' % imaginaryPart.__class__ )

        baseModule.form.__init__( self, pid, label, productFrame, ( _formFactor, realPart, imaginaryPart ) )

    def check( self, info ) :

        return( [] )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        multiplicity = kwargs['multiplicity']
        aveEnergy = [ [ multiplicity.domainMin, multiplicity.domainMin ], [ multiplicity.domainMax, multiplicity.domainMax ] ]
        aveMomentum = [ [ multiplicity.domainMin, 0 ], [ multiplicity.domainMax, 0 ] ]

        return( [ aveEnergy ], [ aveMomentum ] )

    def processMC_cdf( self, style, tempInfo, indent ) :

        return( None )

    def processMultiGroup( self, style, tempInfo, indent ) :

        from fudge.processing.deterministic import transferMatrices as transferMatricesModule
        from fudge.processing import group as groupModule

        verbosity = tempInfo['verbosity']
        if( verbosity > 2 ) : print('%sGrouping %s' % (indent, self.moniker))

        anomalousScatteringFactor = XYsModule.XYs1d( [ [ 0.0, 0.0 ], [ 21., 0.0 ] ] )
        anomalousScatteringFactor_realPart = anomalousScatteringFactor
        if( self.anomalousScatteringFactor_realPart is not None ) :
            anomalousScatteringFactor_realPart = self.anomalousScatteringFactor_realPart.data
        anomalousScatteringFactor_imaginaryPart = anomalousScatteringFactor
        if( self.anomalousScatteringFactor_imaginaryPart is not None ) :
            anomalousScatteringFactor_imaginaryPart = self.anomalousScatteringFactor_imaginaryPart.data

        TM_1, TM_E = transferMatricesModule.wholeAtomScattering( style, tempInfo, self.productFrame, self.formFactor.data,
                realAnomalousFactor = anomalousScatteringFactor_realPart, imaginaryAnomalousFactor = anomalousScatteringFactor_imaginaryPart,
                comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )
        return( groupModule.TMs2Form( style, tempInfo, TM_1, TM_E ) )

    def crossSection(  self, tolerance = 1e-3, interpolation = standardsModule.interpolation.linlinToken ) :
        """
        Integrates the double differential cross section to get the cross section :math:`\\sigma(E)` where :math:`E` is the projectile's energy.
        """

        class Parameters :

            pass

        class tester :

            def __init__( self, parameters, relativeTolerance, absoluteTolerance, anomalousScatteringFactor_realPart, anomalousScatteringFactor_imaginaryPart ) :

                self.parameters = parameters
                self.relativeTolerance = relativeTolerance
                self.absoluteTolerance = absoluteTolerance
                self.anomalousScatteringFactor_realPart = anomalousScatteringFactor_realPart
                self.anomalousScatteringFactor_imaginaryPart = anomalousScatteringFactor_imaginaryPart

            def evaluateAtX( self, energy ) :

                return( integratedCrossSection( energy, self.parameters, self.anomalousScatteringFactor_realPart, self.anomalousScatteringFactor_imaginaryPart ) )

        def integrateSub( _n, _a, logX, x1, _y1, x2, _y2 ) :

            epsilon = _a + _n + 1
            if( abs( epsilon ) < 1e-3 ) :
                epsilon_logX = epsilon * logX
                integral = _y1 * math.pow( x1, _n + 1 ) * logX * ( 1 + 0.5 * epsilon_logX * ( 1 + epsilon_logX / 3.0 * ( 1 + 0.25 * epsilon_logX ) ) )
            else :
                integral = ( _y2 * math.pow( x2, _n + 1 ) - _y1 * math.pow( x1, _n + 1 ) ) / epsilon
            return( integral )
            

        def integratedCrossSection( energy, parameters, anomalousScatteringFactor_realPart, anomalousScatteringFactor_imaginaryPart ) :

            parameters.anomalousScatteringFactor_realPart = anomalousScatteringFactor_realPart.evaluate( energy )
            if( parameters.anomalousScatteringFactor_realPart == None ) : parameters.anomalousScatteringFactor_realPart = 0

            parameters.anomalousScatteringFactor_imaginaryPart = anomalousScatteringFactor_imaginaryPart.evaluate( energy )
            if( parameters.anomalousScatteringFactor_imaginaryPart == None ) : parameters.anomalousScatteringFactor_imaginaryPart = 0.0

            formFactorLower = parameters.formFactor[0]

            integralPrime = 0.0
            integralUpper = 0.0
            if( energy <= formFactorLower.domainMax ) :
                muMin = -1.0
                integralLower = ( formFactorLower[0][1] + parameters.anomalousScatteringFactor_realPart )**2 + parameters.anomalousScatteringFactor_imaginaryPart**2
                integralLower *= 8.0 / 3.0
            else :
                ratio2 = 2.0 * ( formFactorLower.domainMax / energy )**2
                muMin = 1.0 - ratio2
                integralPrime = 8.0 / 3.0 * ( parameters.anomalousScatteringFactor_realPart**2 + parameters.anomalousScatteringFactor_imaginaryPart**2 )
                integralLower = ratio2 * ( 4 + muMin * ( 1 + muMin ) ) / 3.0 * \
                        ( formFactorLower[0][1] * ( formFactorLower[0][1] + 2 * parameters.anomalousScatteringFactor_realPart ) )

            i1 = 0
            if( muMin > -1.0 ) :
                formFactorLower = parameters.formFactor[1]
                for i1, xy in enumerate( formFactorLower ) :
                    _x2, y2 = xy
                    x2 = _x2
                    if( i1 > 0 ) :
                        xRatio = x2 / x1
                        logX = math.log( xRatio )
                        logY = math.log( y2 / y1 )
                        _a = logY / logX

                        if( x2 > energy ) :
                            x2 = energy
                            y2 = y1 * math.pow( energy / x1, _a )
                        x1_e = x1 / energy
                        x2_e = x2 / energy

                        integralSub  = 0.5 * integrateSub( 1, _a, logX, x1_e, y1, x2_e, y2 )
                        integralSub -=       integrateSub( 3, _a, logX, x1_e, y1, x2_e, y2 )
                        integralSub +=       integrateSub( 5, _a, logX, x1_e, y1, x2_e, y2 )
                        integralSub *= 2.0 * parameters.anomalousScatteringFactor_realPart

                        y1_2 = y1 * y1
                        y2_2 = y2 * y2
                        _2_a = 2 * _a
                        integralSub2  = 0.5 * integrateSub( 1, _2_a, logX, x1_e, y1_2, x2_e, y2_2 )
                        integralSub2 -=       integrateSub( 3, _2_a, logX, x1_e, y1_2, x2_e, y2_2 )
                        integralSub2 +=       integrateSub( 5, _2_a, logX, x1_e, y1_2, x2_e, y2_2 )

                        integralUpper += integralSub + integralSub2
                    x1 = x2
                    y1 = y2
                    if( _x2 >= energy ) : break
                integralUpper *= 16

            return( ( integralPrime + integralLower + integralUpper ) * parameters.factor_crossSection )

        parameters = Parameters( )

        product = self.ancestor.ancestor

        tolerance = max( 1e-6, min( 0.1, tolerance ) )
        parameters.tolerance = 1e-4 * tolerance

        factor_crossSection = PQUModule.PQU( 1.0, 'e**2 / ( 4. * pi * eps0 * me * c**2 )' )**2      # r_e^2 where r_e is the classical electron radius.
        factor_crossSection = math.pi * factor_crossSection.getValueAs( 'b' )
        parameters.factor_crossSection = factor_crossSection

        formFactor = self.formFactor.data
        if( isinstance( formFactor, regions1d ) ) :
            _formFactor = formFactor.copy( )
        else :
            raise Exception( 'Form factor must be a regions1d instance' )
        _formFactor.axes[1].unit = product.domainUnit
        factor_E2x = 1.0 / PQUModule.PQU( 1.0, '%s / hplanck / c' % product.domainUnit ).getValueAs( formFactor.axes[1].unit )
        for region in _formFactor : region.scaleOffsetXAndY( factor_E2x, 0.0, 1.0, 0.0, True )
        parameters.formFactor = _formFactor

        domainMin = product.domainMin
        domainMax = product.domainMax
        anomalousScatteringFactor_realPart = self.anomalousScatteringFactor_realPart
        if( anomalousScatteringFactor_realPart is None ) :
            anomalousScatteringFactor_realPart = XYs1d( [ [ domainMin, 0.0 ], [ domainMax, 0.0 ] ], axes = _formFactor.axes.copy( ) )
        else :
            anomalousScatteringFactor_realPart = anomalousScatteringFactor_realPart.data

        anomalousScatteringFactor_imaginaryPart = self.anomalousScatteringFactor_imaginaryPart
        if( anomalousScatteringFactor_imaginaryPart is None ) :
            anomalousScatteringFactor_imaginaryPart = XYs1d( [ [ domainMin, 0.0 ], [ domainMax, 0.0 ] ], axes = _formFactor.axes.copy( ) )
        else :
            anomalousScatteringFactor_imaginaryPart = anomalousScatteringFactor_imaginaryPart.data

        transitionEnergy = PQUModule.PQU( 1e3, 'eV' ).getValueAs( product.domainUnit )
        energy = domainMin
        initialEnergies = []
        while( energy < transitionEnergy ) :
            initialEnergies.append( energy )
            energy *= 1.01
        energy = transitionEnergy
        while( energy < domainMax ) :
            initialEnergies.append( energy )
            energy *= 10
        initialEnergies.append( domainMax )

        energy_crossSection = []
        for energy in initialEnergies :
            energy_crossSection.append( [ energy, integratedCrossSection( energy, parameters, anomalousScatteringFactor_realPart, anomalousScatteringFactor_imaginaryPart ) ] )

        _tester = tester( parameters, tolerance, tolerance * energy_crossSection[-1][1], 
                anomalousScatteringFactor_realPart, anomalousScatteringFactor_imaginaryPart )
        energy_crossSection = fudgemathModule.thickenXYList( energy_crossSection, _tester, biSectionMax = 10, interpolation = interpolation )
        return( crossSectionModule.XYs1d( data = energy_crossSection, axes = crossSectionModule.defaultAxes( product.domainUnit ),
                interpolation = interpolation ) )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )

        subForms = { formFactor.moniker : None, realAnomalousFactor.moniker : None, imaginaryAnomalousFactor.moniker : None }
        for child in element :
            for _class in ( formFactor, realAnomalousFactor, imaginaryAnomalousFactor ) :
                if( child.tag == _class.moniker ) : break

            if( child.tag != _class.moniker ) : raise TypeError( "Invalid element '%s' encountered '%s'" % ( child.tag, _class.moniker ) )
            subForms[_class.moniker] = _class.parseXMLNode( child, xPath, linkData )

        _form = form( element.get( 'pid' ), element.get( 'label' ), element.get( 'productFrame' ),
                        subForms[formFactor.moniker], subForms[realAnomalousFactor.moniker],
                        subForms[imaginaryAnomalousFactor.moniker] )

        xPath.pop( )
        return( _form )
