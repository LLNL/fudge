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

"""
Photon coherent double difference cross section form and its supporting classes.
"""

import math

from pqu import PQU as PQUModule

from numericalFunctions import integration as nf_integrationModule

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

    pass

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
            data = baseModule.XYs1d.parseXMLNode( data, xPath, linkData )
        elif( data.tag == baseModule.regions1d.moniker ) :
            data = baseModule.regions1d.parseXMLNode( data, xPath, linkData )
        else :
            raise TypeError( 'Invalid data "%s" for "%s"' % ( data.tag, cls.tag ) )

        xPath.pop( )
        return( cls( data ) )

class scatteringFunction( coherentFunctionBase ) :

    moniker = 'scatteringFactor'
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

    def __init__( self, pid, label, productFrame, formFactor, realPart, imaginaryPart ) :

        if( not( isinstance( formFactor, scatteringFunction ) ) ) :
            raise Exception( 'Instance is class "%s" and not scatteringFunction' % formFactor.__class__ )

        if( realPart is not None ) :
            if( not( isinstance( realPart, realAnomalousFactor ) ) ) :
                raise Exception( 'Instance is class "%s" and not realAnomalousFactor' % realPart.__class__ )

        if( imaginaryPart is not None ) :
            if( not( isinstance( imaginaryPart, imaginaryAnomalousFactor ) ) ) :
                raise Exception( 'Instance is class "%s" and not imaginaryAnomalousFactor' % imaginaryPart.__class__ )

        baseModule.form.__init__( self, pid, label, productFrame, ( formFactor, realPart, imaginaryPart ) )

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
        if( verbosity > 2 ) : print '%sGrouping %s' % ( indent, self.moniker )

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

    def crossSection( self, tolerance = 1e-3, interpolation = standardsModule.interpolation.linlinToken ) :
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

        def integrand( mu, parameters ) :

            _x = parameters.energy_hc * math.sqrt( 0.5 * ( 1.0 - mu ) )
            integral = ( parameters.formFactor.evaluate( _x ) + parameters.anomalousScatteringFactor_realPart )**2
            integral += parameters.anomalousScatteringFactor_imaginaryPart**2
            integral *= 1 + mu**2
            return( integral )

        def integratedCrossSection( energy, parameters, anomalousScatteringFactor_realPart, anomalousScatteringFactor_imaginaryPart ) :

            parameters.energy = energy
            parameters.energy_hc = energy * parameters.factor_E2x

            parameters.anomalousScatteringFactor_realPart = 0
            if( anomalousScatteringFactor_realPart is not None ) :
                parameters.anomalousScatteringFactor_realPart = anomalousScatteringFactor_realPart.evaluate( energy )
                if( parameters.anomalousScatteringFactor_realPart is None ) : parameters.anomalousScatteringFactor_realPart = 0

            parameters.anomalousScatteringFactor_imaginaryPart = 0
            if( anomalousScatteringFactor_imaginaryPart is not None ) :
                parameters.anomalousScatteringFactor_imaginaryPart = anomalousScatteringFactor_imaginaryPart.evaluate( energy )
                if( parameters.anomalousScatteringFactor_imaginaryPart is None ) : parameters.anomalousScatteringFactor_imaginaryPart = 0

            if( energy <= parameters.transitionEnergy ) :
                xSec, evaluations = nf_integrationModule.adaptiveQuadrature_GnG( 4, integrand, parameters, -1.0, 1.0, parameters.tolerance, 16 )
            else :
                mu = 1.0 - 0.1 * ( parameters.transitionEnergy / energy )**1.5
                xSec1, evaluations1 = nf_integrationModule.adaptiveQuadrature_GnG( 4, integrand, parameters, -1.0, mu, parameters.tolerance, 16 )
                xSec2, evaluations2 = nf_integrationModule.adaptiveQuadrature_GnG( 4, integrand, parameters, mu, 1.0, parameters.tolerance, 16 )
                evaluations = evaluations1 + evaluations2
                xSec = xSec1 + xSec2
            xSec *= parameters.factor_crossSection
            return( xSec )

        product = self.ancestor.ancestor
        formFactor = self.formFactor.data

        tolerance = max( 1e-6, min( 0.1, tolerance ) )
        parameters = Parameters( )
        parameters.tolerance = 1e-4 * tolerance
        parameters.formFactor = formFactor.toPointwise_withLinearXYs( lowerEps = 1e-12 )

        factor_crossSection = PQUModule.PQU( 1.0, 'e**2 / ( 4. * pi * eps0 * me * c**2 )' )**2
        factor_crossSection = math.pi * factor_crossSection.getValueAs( 'b' )
        parameters.factor_crossSection = factor_crossSection
        parameters.factor_E2x = PQUModule.PQU( 1.0, '%s / hplanck / c' % product.domainUnit ).getValueAs( formFactor.axes[1].unit )
        parameters.transitionEnergy = PQUModule.PQU( 2e5, 'eV' ).getValueAs( product.domainUnit )

        anomalousScatteringFactor_realPart = self.anomalousScatteringFactor_realPart
        if( anomalousScatteringFactor_realPart is not None ) : anomalousScatteringFactor_realPart = anomalousScatteringFactor_realPart.data
        anomalousScatteringFactor_imaginaryPart = self.anomalousScatteringFactor_imaginaryPart
        if( anomalousScatteringFactor_imaginaryPart is not None ) : anomalousScatteringFactor_imaginaryPart = anomalousScatteringFactor_imaginaryPart.data

        initialEnergies = [ product.domainMin, product.domainMax ]
        energy_crossSection = []
        for energy in initialEnergies :
            energy_crossSection.append( [ energy, integratedCrossSection( energy, parameters, anomalousScatteringFactor_realPart, anomalousScatteringFactor_imaginaryPart ) ] )
        _tester = tester( parameters, tolerance, tolerance * energy_crossSection[-1][1], 
                anomalousScatteringFactor_realPart, anomalousScatteringFactor_imaginaryPart )
        energy_crossSection = fudgemathModule.thickenXYList( energy_crossSection, _tester, biSectionMax = 20, interpolation = interpolation )

        return( crossSectionModule.XYs1d( data = energy_crossSection, axes = crossSectionModule.defaultAxes( product.domainUnit ),
                interpolation = interpolation ) )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )

        subForms = { scatteringFunction.moniker : None, realAnomalousFactor.moniker : None, imaginaryAnomalousFactor.moniker : None }
        for child in element :
            for _class in ( scatteringFunction, realAnomalousFactor, imaginaryAnomalousFactor ) :
                if( child.tag == _class.moniker ) : break

            if( child.tag != _class.moniker ) : raise TypeError( "Invalid element '%s' encountered '%s'" % ( child.tag, self.moniker ) )
            subForms[_class.moniker] = _class.parseXMLNode( child, xPath, linkData )

        _form = form( element.get( 'pid' ), element.get( 'label' ), element.get( 'productFrame' ),
                        subForms[scatteringFunction.moniker], subForms[realAnomalousFactor.moniker],
                        subForms[imaginaryAnomalousFactor.moniker] )

        xPath.pop( )
        return( _form )
