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

"""Angular/energy double differential distribution classes."""

import math
import fudge
from fudge.core.utilities import brb

from PoPs import IDs as IDsPoPsModule

from xData import standards as standardsModule
from xData import axes as axesModule
from xData import XYs as XYsModule
from xData import multiD_XYs as multiD_XYsModule

from fudge.processing import group as groupModule

from fudge.gnd.productData import momentumDeposition as momentumDepositionModule

from . import miscellaneous as miscellaneousModule

from . import base as baseModule
from . import angular as angularModule
from . import angularEnergyMC as angularEnergyMCModule

__metaclass__ = type

def defaultAxes( energyUnit, probabilityLabel = 'P(mu,energy_out|energy_in)' ) :

    axes = axesModule.axes( rank = 4 )
    axes[3] = axesModule.axis( 'energy_in', 3, energyUnit )
    axes[2] = axesModule.axis( 'mu', 2, '' )
    axes[1] = axesModule.axis( 'energy_out', 1, energyUnit )
    axes[0] = axesModule.axis( probabilityLabel, 0, '1/' + energyUnit )
    return( axes )

class XYs1d( XYsModule.XYs1d ) :

    def averageEnergy( self ) :

        return( self.integrateWithWeight_x( ) )
    
    def averageMomentum( self ) :

        return( self.integrateWithWeight_sqrt_x( ) )
    
class XYs2d( multiD_XYsModule.XYs2d ) :

    def averageEnergy( self ) :

        EpOfMu = [ [ pdfOfEpAtMu.value, pdfOfEpAtMu.averageEnergy( ) ] for pdfOfEpAtMu in self ]
        return( float( XYsModule.XYs1d( EpOfMu ).integrate( ) ) )

    def averageMomentum( self ) :

        MpOfMu = [ [ pdfOfMpAtMu.value, pdfOfMpAtMu.averageMomentum( ) ] for pdfOfMpAtMu in self ]
        return( float( XYsModule.XYs1d( MpOfMu ).integrate( ) ) )

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs1d, ) )

class subform( baseModule.subform ) :
    """Abstract base class for angularEnergy forms."""

    pass

class XYs3d( subform, multiD_XYsModule.XYs3d ) :

    def __init__( self, **kwargs ) :

        multiD_XYsModule.XYs3d.__init__( self, **kwargs )
        subform.__init__( self )

    def check( self, info ) :
        """
        Check for incomplete angular distributions + any negative probabilities.
        Ignore normalization for this double-differential distribution.
        """

        from fudge.gnd import warning
        from pqu import PQU
        warnings = []

        for index, energy_in in enumerate(self):
            integral = energy_in.integrate()
            if abs(integral - 1.0) > info['normTolerance']:
                warnings.append( warning.unnormalizedDistribution( PQU.PQU(energy_in.value, self.axes[0].unit),
                    index, integral, self.toXLink() ) )
            if( energy_in.domainMin != -1 ) or ( energy_in.domainMax != 1 ) :
                warnings.append( warning.incompleteDistribution( PQU.PQU(energy_in.value, self.axis[0].unit),
                        energy_in.domainMin, energy_in.domainMax, energy_in ) )
            for mu in energy_in:
                if( mu.domainMin < 0 ) :
                    warnings.append( warning.valueOutOfRange("Negative outgoing energy for energy_in=%s!"
                        % PQU.PQU(energy_in.value, self.axes[0].unit), mu.domainMin, 0, 'inf', self.toXLink() ) )
                if( mu.rangeMin < 0 ) :
                    warnings.append( warning.negativeProbability( PQU.PQU(energy_in.value, self.axes[-1].unit),
                        mu=mu.value, obj=mu ) )

        return warnings

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        multiplicity = kwargs['multiplicity']
        productMass = kwargs['productMass']

        aveEnergy = []
        for pdfOfMuEpAtE in self :
            energy = pdfOfMuEpAtE.value
            aveEnergy.append( [ energy, multiplicity.evaluate( energy ) * pdfOfMuEpAtE.averageEnergy( ) ] )

        const = math.sqrt( 2. * productMass )
        aveMomentum = []
        for pdfOfMuEpAtE in self :
            energy = pdfOfMuEpAtE.value
            momemtum = const * multiplicity.evaluate( energy ) * pdfOfMuEpAtE.averageMomentum( )
            if( momemtum < 1e-12 ) : momemtum = 0.          # This should be less arbitrary????????
            aveMomentum.append( [ energy, momemtum ] )

        return( aveEnergy, aveMomentum )

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

        if( verbosity > 2 ) : print '%sGrouping %s' % ( indent, self.moniker )
        TM_1, TM_E = transferMatricesModule.ENDFEMuEpP_TransferMatrix( style, tempInfo, productFrame, tempInfo['crossSection'], self,
            tempInfo['multiplicity'], comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )

        return( groupModule.TMs2Form( style, tempInfo, TM_1, TM_E ) )

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :

        angular = angularModule.XYs2d( axes = self.axes )
        for PofEpGivenMu in self :
            data = angularModule.XYs1d( [ [ PofEp.value, float( PofEp.integrate( ) ) ] for PofEp in PofEpGivenMu ] )
            angular.append( angularModule.xs_pdf_cdf1d.fromXYs( angularModule.XYs1d( data ),
                    value = PofEpGivenMu.value ) )

        xys3d = angularEnergyMCModule.XYs3d( axes = self.axes )
        for PofEpGivenMu in self :
            xys2d = angularEnergyMCModule.XYs2d( value = PofEpGivenMu.value )
            for PofEp in PofEpGivenMu :
                _PofEp = PofEp.toPointwise_withLinearXYs( accuracy = 1e-3, upperEps = 1e-8 )
                xys2d.append( angularEnergyMCModule.xs_pdf_cdf1d.fromXYs( _PofEp, PofEp.value ) )
            xys3d.append( xys2d )

        return( angularEnergyMCModule.angular( angular ), angularEnergyMCModule.angularEnergy( xys3d ) )

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs2d, ) )

class LLNLAngularOfAngularEnergySubform( baseModule.subform ) :

    moniker = 'LLNLAngularOfAngularEnergy'

    def __init__( self, data ) :

        if( not( isinstance( data, angularModule.XYs2d ) ) ) : raise TypeError( 'instance is not an angular.XYs2d' )
        baseModule.subform.__init__( self )
        self.data = data

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.data.convertUnits( unitMap )

    def copy( self ) :

        return( LLNLAngularOfAngularEnergySubform( self.data.copy ) )

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :

        return( angularEnergyMCModule.angular( self.data.to_xs_pdf_cdf1d( style, tempInfo, indent ) ) )

    def toXMLList( self, indent = "", **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        xmlStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        xmlStringList += self.data.toXMLList( indent = indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )
        data = angularModule.XYs2d.parseXMLNode( element[0], xPath, linkData )
        result = LLNLAngularOfAngularEnergySubform( data )
        xPath.pop( )
        return( result )

class LLNLAngularEnergyOfAngularEnergySubform( baseModule.subform ) :

    moniker = 'LLNLAngularEnergyOfAngularEnergy'

    def __init__( self, data ) :

        if( not( isinstance( data, XYs3d ) ) ) : raise TypeError( 'instance is not an angularEnergy.XYs3d' )
        baseModule.subform.__init__( self )
        self.data = data

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.data.convertUnits( unitMap )

    def copy( self ) :

        return( LLNLAngularEnergyOfAngularEnergySubform( self.data.copy( ) ) )

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :

        return( self.data.to_xs_pdf_cdf1d( style, tempInfo, indent ) )

    def toXMLList( self, indent = "", **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        xmlStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        xmlStringList += self.data.toXMLList( indent = indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )
        data = XYs3d.parseXMLNode( element[0], xPath, linkData )
        result = LLNLAngularEnergyOfAngularEnergySubform( data )
        xPath.pop( )
        return( result )

class LLNLAngularEnergyForm( baseModule.form ) :
    """This class is for legacy LLNL ENDL I = 1 and 3 data and is deprecated. Only use for the legacy ENDL data."""

    moniker = 'LLNLAngularEnergy'
    subformAttributes = ( 'angularSubform', 'angularEnergySubform' )

    def __init__( self, label, productFrame, angularSubform, angularEnergySubform ) :

        if( not( isinstance( angularSubform, LLNLAngularOfAngularEnergySubform ) ) ) :
            raise TypeError( 'instance is not an LLNLAngularOfAngularEnergySubform subform' )
        if( not( isinstance( angularEnergySubform, LLNLAngularEnergyOfAngularEnergySubform ) ) ) :
            raise TypeError( 'instance is not an LLNLAngularEnergyOfAngularEnergySubform subform' )
        baseModule.form.__init__( self, label, productFrame, ( angularSubform, angularEnergySubform ) )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        from fudge.core.math import fudgemath

        def calculateDepositionMomentumAtMu( mu, parameters ) :

            f = ( mu - parameters.mu1 ) / ( parameters.mu2 - parameters.mu1 )
            P_mu = ( 1 - f ) * parameters.P1 + f * parameters.P2
            EpP = parameters.muEpPs.interpolateAtValue( mu, unitBase = True, extrapolation = standardsModule.flatExtrapolationToken )
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
                    p += miscellaneousModule.GnG_adaptiveQuadrature( calculateDepositionMomentumAtMu, mu1, mu2, parameters, 
                            miscellaneousModule.GaussQuadrature2, tolerance = 1e-4 )[0]
                mu1 = mu2
            return( p )

        class calculateDepositionMomentumThicken :

            def __init__( self, angularSubform, angularEnergySubform, relativeTolerance, absoluteTolerance ) :

                self.angularSubform = angularSubform
                self.angularEnergySubform = angularEnergySubform
                self.relativeTolerance = relativeTolerance
                self.absoluteTolerance = absoluteTolerance

            def evaluateAtX( self, E ) :

                muPs = self.angularSubform.getAtEnergy( E )
                muEpPs = self.angularEnergySubform.interpolateAtV( E, unitBase = True )
                return( calculateDepositionMomentum( muPs, muEpPs ) )

        energyUnit = kwargs['incidentEnergyUnit']
        momentumDepositionUnit = kwargs['momentumDepositionUnit']
        multiplicity = kwargs['multiplicity']
        energyAccuracy = kwargs['energyAccuracy']
        momentumAccuracy = kwargs['momentumAccuracy']
        product = kwargs['product']
        productMass = kwargs['productMass']

        angularSubform = self.angularSubform.data
        angularEnergySubform = self.angularEnergySubform.data

        depEnergy = miscellaneousModule.calculateDepositionEnergyFromAngular_angularEnergy( style, 
                angularSubform, angularEnergySubform, multiplicity, accuracy = energyAccuracy )

        if( product.id == IDsPoPsModule.photon ) :
            depMomentum = miscellaneousModule.calculateDepositionEnergyFromAngular_angularEnergy( style, 
                    angularSubform, angularEnergySubform, multiplicity, True, accuracy = momentumAccuracy )
        else :
            relativeTolerance, absoluteTolerance = momentumAccuracy, momentumAccuracy
            depMomentum = []
            for indexE, muPs in enumerate( angularSubform ) :
                depMomentum.append( [ muPs.value, calculateDepositionMomentum( muPs, angularEnergySubform[indexE] ) ] )
#            depMomentum = fudgemath.thickenXYList( depMomentum, calculateDepositionMomentumThicken( angularSubform, angularEnergySubform, relativeTolerance, absoluteTolerance ) )
            const = math.sqrt( 2. * productMass )
            for EMomenutem in depMomentum : EMomenutem[1] *= const * multiplicity.evaluate( EMomenutem[0] )

        axes = momentumDepositionModule.defaultAxes( energyUnit = energyUnit, momentumDepositionUnit = momentumDepositionUnit )
        depMomentum = momentumDepositionModule.XYs1d( data = depMomentum, axes = axes,
                        label = style.label )

        return( [ depEnergy ], [ depMomentum ] )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.angularSubform.convertUnits( unitMap )
        self.angularEnergySubform.convertUnits( unitMap )

    def processMC( self, style, tempInfo, indent ) :

        angular = self.angularSubform.to_xs_pdf_cdf1d( style, tempInfo, indent )
        dummy, angularEnergy = self.angularEnergySubform.to_xs_pdf_cdf1d( style, tempInfo, indent )
        return( angularEnergyMCModule.form( style.label, self.productFrame, angular, angularEnergy ) )

    def processMultiGroup( self, style, tempInfo, indent ) :

        from fudge.processing.deterministic import transferMatrices as transferMatricesModule

        verbosity = tempInfo['verbosity']
        if( verbosity > 2 ) : print '%sGrouping %s' % ( indent, self.moniker )

        angularSubform = self.angularSubform.data
        angularEnergySubform = self.angularEnergySubform.data

        TM_1, TM_E = transferMatricesModule.ENDLEMuEpP_TransferMatrix( style, tempInfo, tempInfo['crossSection'], self.productFrame, 
            angularSubform, angularEnergySubform, tempInfo['multiplicity'],
            comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )

        return( groupModule.TMs2Form( style, tempInfo, TM_1, TM_E ) )

    def toXMLList( self, indent = "", **kwargs ) : 

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        indent3 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attributeStr = ''
        if( self.label is not None ) : attributeStr += ' label="%s"' % self.label
        if( self.productFrame is not None ) : attributeStr += ' productFrame="%s"' % self.productFrame
        xmlString = [ '%s<%s%s>' % ( indent, self.moniker, attributeStr ) ]

        xmlString += self.angularSubform.toXMLList( indent3, **kwargs )

        xmlString += self.angularEnergySubform.toXMLList( indent3, **kwargs )

        xmlString[-1] += '</%s>' % self.moniker 
        return( xmlString )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )
        for child in element :
            if( child.tag == LLNLAngularOfAngularEnergySubform.moniker ) :
                angularSubform = LLNLAngularOfAngularEnergySubform.parseXMLNode( child, xPath, linkData )
            elif( child.tag == LLNLAngularEnergyOfAngularEnergySubform.moniker ) :
                angularEnergySubform = LLNLAngularEnergyOfAngularEnergySubform.parseXMLNode( child, xPath, linkData )
            else :
                raise ValueError( 'Encountered unexpected child element "%s"' % child.tag )
        component = LLNLAngularEnergyForm( element.get( 'label' ), element.get( 'productFrame' ), 
                angularSubform, angularEnergySubform )

        xPath.pop( )
        return( component )

class form(  baseModule.form ) :

    moniker = 'angularEnergy'
    subformAttributes = ( 'angularEnergySubform', )
    ancestryMembers = subformAttributes

    def __init__( self, label, productFrame, angularEnergySubform ) :

        if( not( isinstance( angularEnergySubform, subform ) ) ) : raise TypeError( 'instance is not an angularEnergy subform' )
        baseModule.form.__init__( self, label, productFrame, ( angularEnergySubform, ) )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        if( self.productFrame == standardsModule.frames.centerOfMassToken ) :
            raise Exception( 'center of mass calculation not supported for %s' % self.moniker )

        kwargs['productFrame'] = self.productFrame
        aveEnergy, aveMomentum = self.angularEnergySubform.calculateAverageProductData( style, indent = indent, **kwargs )
        return( [ aveEnergy ], [ aveMomentum ] )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.angularEnergySubform.convertUnits( unitMap )

    def processMC( self, style, tempInfo, indent ) :

        angular, angularEnergy = self.angularEnergySubform.to_xs_pdf_cdf1d( style, tempInfo, indent )
        return( angularEnergyMCModule.form( style.label, self.productFrame, angular, angularEnergy ) )

    def processMultiGroup( self, style, tempInfo, indent ) :

        tempInfo['productFrame'] = self.productFrame
        return( self.angularEnergySubform.processMultiGroup( style, tempInfo, indent ) )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( element.tag )
        formElement = element[0]
        subformClass = {
            XYs3d.moniker :  XYs3d,
        }[ formElement.tag ]
        if subformClass is None: raise Exception("encountered unknown angularEnergy subform: %s" % formElement.tag)
        subform = XYs3d.parseXMLNode( formElement, xPath, linkData )
        AEC = form( element.get( 'label' ), element.get('productFrame'), subform )
        xPath.pop()
        return AEC
