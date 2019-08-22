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
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# <<END-copyright>>

""" energy/angular double differential distribution classes """

import math
import miscellaneous
import fudge
from fudge.core.utilities import brb

from xData import standards as standardsModule
from xData import axes as axesModule
from xData import XYs as XYsModule
from xData import multiD_XYs as multiD_XYsModule

from fudge.gnd.productData import energyDeposition as energyDepositionModule
from fudge.gnd.productData import momentumDeposition as momentumDepositionModule

from . import base as baseModule
from . import angular as angularModule

__metaclass__ = type

class pdfOfEp :

    class pointwise( XYsModule.XYs ) :

        def averageEnergy( self ) :

            return( self.integrateWithWeight_x( ) )
    
        def averageMomentum( self ) :

            return( self.integrateWithWeight_sqrt_x( ) )
    
class pdfOfMuAndEp :

    class pointwise( multiD_XYsModule.multiD_XYs ) :

        def averageEnergy( self ) :

            EpOfMu = [ [ pdfOfEpAtMu.value, pdfOfEpAtMu.averageEnergy( ) ] for pdfOfEpAtMu in self ]
            return( float( XYsModule.XYs( EpOfMu ).integrate( ) ) )

        def averageMomentum( self ) :

            MpOfMu = [ [ pdfOfMpAtMu.value, pdfOfMpAtMu.averageMomentum( ) ] for pdfOfMpAtMu in self ]
            return( float( XYsModule.XYs( MpOfMu ).integrate( ) ) )

        @staticmethod
        def allowedSubElements( ) :

            return( ( pdfOfEp.pointwise, ) )

class subform( baseModule.subform ) :
    """Abstract base class for angularEnergy forms."""

    pass

class pointwise( subform, multiD_XYsModule.multiD_XYs ) :

    def __init__( self, **kwargs ) :

        if( 'dimension' not in kwargs ) : kwargs['dimension'] = 3
        if( kwargs['dimension'] != 3 ) : raise ValueError( 'Dimension = %s != 3' % ( kwargs['dimension'] ) )
        multiD_XYsModule.multiD_XYs.__init__( self, **kwargs )
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
            if energy_in.domainMin() != -1 or energy_in.domainMax() != 1 :
                warnings.append( warning.incompleteDistribution( PQU.PQU(energy_in.value, self.axis[0].unit),
                        energy_in.domainMin(), energy_in.domainMax(), energy_in.value, self.toXLink() ) )
            for mu in energy_in:
                if mu.domainMin() < 0:
                    warnings.append( warning.valueOutOfRange("Negative outgoing energy for energy_in=%s!"
                        % PQU.PQU(energy_in.value, self.axes[0].unit), mu.domainMin(), 0, 'inf', self.toXLink() ) )
                if mu.rangeMin() < 0:
                    warnings.append( warning.negativeProbability( PQU.PQU(energy_in.value, self.axes[-1].unit),
                        mu=mu.value, obj=mu ) )

        return warnings

    def calculateDepositionData( self, processInfo, tempInfo, verbosityIndent ) :

        if( self.productFrame == standardsModule.frames.centerOfMassToken ) : raise Exception( 'center of mass calculation not supported for %s' % self.moniker )

        energyUnit = tempInfo['incidentEnergyUnit']
        momentumDepositionUnit = energyUnit + '/c'
        massUnit = energyUnit + '/c**2'
        energyAccuracy, momentumAccuracy = processInfo.energyAccuracy, processInfo.momentumAccuracy
        multiplicity = tempInfo['multiplicity']
        productMass = tempInfo['product'].getMass( massUnit )

        depEnergy = []

        for pdfOfMuEpAtE in self :
            energy = pdfOfMuEpAtE.value
            depEnergy.append( [ energy, multiplicity.getValue( energy ) * pdfOfMuEpAtE.averageEnergy( ) ] )

        const = math.sqrt( 2. * productMass )
        depMomentum = []
        for pdfOfMuEpAtE in self :
            energy = pdfOfMuEpAtE.value
            momemtum = const * multiplicity.getValue( energy ) * pdfOfMuEpAtE.averageMomentum( )
            if( momemtum < 1e-12 ) : momemtum = 0.          # This should be less arbitrary????????
            depMomentum.append( [ energy, momemtum ] )

        axes = energyDepositionModule.pointwise.defaultAxes( energyUnit = energyUnit, energyDepositionUnit = energyUnit )
        depEnergy = energyDepositionModule.pointwise( data = depEnergy, axes = axes,
                label = processInfo.style.label, accuracy = energyAccuracy )
        axes = momentumDepositionModule.pointwise.defaultAxes( energyUnit = energyUnit, momentumDepositionUnit = momentumDepositionUnit )
        depMomentum = momentumDepositionModule.pointwise( data = depMomentum, axes = axes,
                label = processInfo.style.label, accuracy = momentumAccuracy )

        return( [ depEnergy, depMomentum ] )

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

        if( 'LLNL_Pn' in processInfo['styles'] ) :
            if( processInfo.verbosity >= 30 ) : print '%sGrouping %s' % ( verbosityIndent, self.moniker )
            outputChannel = tempInfo['outputChannel'].outputChannel
            projectile, product = processInfo.getProjectileName( ), tempInfo['product'].particle.name
            TM_1, TM_E = transferMatrices.ENDFEMuEpP_TransferMatrix( processInfo, projectile, product, tempInfo['masses'], tempInfo['crossSection'], self, 
                tempInfo['multiplicity'], comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )
            fudge.gnd.miscellaneous.TMs2Form( processInfo, tempInfo, newComponents, TM_1, TM_E, self.axes )

        return( newComponents )

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        return( V_W_XYs.V_W_XYs.toPointwise_withLinearXYs( self, accuracy, lowerEps = lowerEps, upperEps = upperEps, cls = pointwise ) )

    @staticmethod
    def allowedSubElements( ) :

        return( ( pdfOfMuAndEp.pointwise, ) )

    @staticmethod
    def defaultAxes( energyUnit = 'eV', energy_outUnit = 'eV', probabilityUnit = '1/eV', probabilityLabel = 'P(mu,energy_out|energy_in)' ) :

        axes = axesModule.axes( rank = 4 )
        axes[3] = axesModule.axis( 'energy_in', 3, energyUnit )
        axes[2] = axesModule.axis( 'mu', 2, '' )
        axes[1] = axesModule.axis( 'energy_out', 1, energy_outUnit )
        axes[0] = axesModule.axis( probabilityLabel, 0, probabilityUnit )
        return( axes )

class LLNLAngularOfAngularEnergySubform( baseModule.subform ) :

    moniker = 'LLNLAngularOfAngularEnergy'

    def __init__( self, angularSubform, makeCopy = True ) :

        if( not( isinstance( angularSubform, angularModule.pointwise ) ) ) : raise TypeError( 'instance is not an angular.pointwise' )
        if( makeCopy ) : angularSubform = angularSubform.copy( )
        baseModule.subform.__init__( self )
        self.angularSubform = angularSubform

    def copy( self ) :

        return( LLNLAngularOfAngularEnergySubform( self.angularSubform, makeCopy = True ) )

    def toXMLList( self, indent = "", **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        xmlStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        xmlStringList += self.angularSubform.toXMLList( indent = indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

class LLNLAngularEnergyOfAngularEnergySubform( baseModule.subform ) :

    moniker = 'LLNLAngularEnergyOfAngularEnergy'

    def __init__( self, angularEnergySubform, makeCopy = True ) :

        if( not( isinstance( angularEnergySubform, pointwise ) ) ) : raise TypeError( 'instance is not an angular.pointwise' )
        if( makeCopy ) : angularEnergySubform = angularEnergySubform.copy( )
        baseModule.subform.__init__( self )
        self.angularEnergySubform = angularEnergySubform

    def copy( self ) :

        return( LLNLAngularEnergyOfAngularEnergySubform( self.angularEnergySubform, makeCopy = True ) )

    def toXMLList( self, indent = "", **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        xmlStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        xmlStringList += self.angularEnergySubform.toXMLList( indent = indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

class LLNLAngularEnergyForm( baseModule.form ) :
    """This class is for legacy LLNL ENDL I = 1 and 3 data and is deprecated. Only use for the legacy ENDL data."""

    moniker = 'LLNLAngularEnergy'
    subformAttributes = ( 'angularSubform', 'angularEnergySubform' )

    def __init__( self, label, productFrame, angularSubform, angularEnergySubform, makeCopy = True ) :

        if( not( isinstance( angularSubform, LLNLAngularOfAngularEnergySubform ) ) ) :
            raise TypeError( 'instance is not an LLNLAngularOfAngularEnergySubform subform' )
        if( not( isinstance( angularEnergySubform, LLNLAngularEnergyOfAngularEnergySubform ) ) ) :
            raise TypeError( 'instance is not an LLNLAngularEnergyOfAngularEnergySubform subform' )
        if( makeCopy ) : angularSubform = angularSubform.copy( )
        if( makeCopy ) : angularEnergySubform = angularEnergySubform.copy( )
        baseModule.form.__init__( self, label, productFrame, ( angularSubform, angularEnergySubform ) )

    def calculateDepositionData( self, processInfo, tempInfo, verbosityIndent ) :

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

            def __init__( self, angularSubform, angularEnergySubform, relativeTolerance, absoluteTolerance ) :

                self.angularSubform = angularSubform
                self.angularEnergySubform = angularEnergySubform
                self.relativeTolerance = relativeTolerance
                self.absoluteTolerance = absoluteTolerance

            def evaluateAtX( self, E ) :

                muPs = self.angularSubform.getAtEnergy( E )
                muEpPs = self.angularEnergySubform.interpolateAtV( E, unitBase = True )
                return( calculateDepositionMomentum( muPs, muEpPs ) )

        depData = []
        angularSubform = self.angularSubform
        angularEnergySubform = self.angularEnergySubform
        multiplicity = tempInfo['multiplicity']
        energyUnit = angularEnergySubform.axes[0].getUnit( )
        energyPrimeUnit = angularEnergySubform.axes[2].getUnit( )
        momentumDepositionUnit = energyPrimeUnit + '/c'
        energyAccuracy, momentumAccuracy = processInfo.energyAccuracy, processInfo.momentumAccuracy

        depData.append( miscellaneous.calculateDepositionEnergyFromAngular_angularEnergy( processInfo, 
                angularSubform, angularEnergySubform, multiplicity, accuracy = energyAccuracy ) )

        if( tempInfo['product'].name == 'gamma' ) :
            depMomentum = miscellaneous.calculateDepositionEnergyFromAngular_angularEnergy( processInfo, 
                    angularSubform, angularEnergySubform, multiplicity, True, accuracy = momentumAccuracy )
        else :
            relativeTolerance, absoluteTolerance = momentumAccuracy, momentumAccuracy
            depMomentum = []
            for indexE, muPs in enumerate( angularSubform ) :
                depMomentum.append( [ muPs.value, calculateDepositionMomentum( muPs, angularEnergySubform[indexE] ) ] )
#            depMomentum = fudgemath.thickenXYList( depMomentum, calculateDepositionMomentumThicken( angularSubform, angularEnergySubform, relativeTolerance, absoluteTolerance ) )
            const = math.sqrt( 2. * tempInfo['product'].getMass( energyPrimeUnit + '/c**2' ) )
            for EMomenutem in depMomentum : EMomenutem[1] *= const * multiplicity.getValue( EMomenutem[0] )

        axes = momentumDepositionModule.pointwise.defaultAxes( energyUnit = energyUnit, momentumDepositionUnit = momentumDepositionUnit )
        depData.append( momentumDepositionModule.pointwise( data = depMomentum, axes = axes,
                        label = processInfo.style.label, accuracy = momentumAccuracy ) )

        return( depData )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        import angular 
        from fudge.processing.deterministic import transferMatrices

        newForms = []
        angularSubform = self.angularSubform
        angularEnergySubform = self.angularEnergySubform

        if( 'LLNL_Pn' in processInfo['styles'] ) :
            if( processInfo.verbosity >= 30 ) : print '%sGrouping %s' % ( verbosityIndent, self.moniker )
            outputChannel = tempInfo['outputChannel'].outputChannel
            projectile, product = processInfo.getProjectileName( ), tempInfo['product'].particle.name
            TM_1, TM_E = transferMatrices.EMuEpP_TransferMatrix( processInfo, projectile, product, tempInfo['masses'],
                tempInfo['crossSection'], angularSubform, angularEnergySubform, tempInfo['multiplicity'],
                comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )
            fudge.gnd.miscellaneous.TMs2Form( processInfo, tempInfo, newForms, TM_1, TM_E, angularSubform.axes )

        return( newForms )

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

        import fudge.gnd.productData.distributions.angular as angularModule

        xPath.append( element.tag )
        for child in element :
            label = child.get( 'label' )
            if( label == 'angular' ) :                  # BRB - FIXME
                angularSubform = angularModule.form.parseXMLNode( child, xPath, linkData )
            elif( label == 'angularEnergy' ) :          # BRB - FIXME
                angularEnergySubform = pointwise.parseXMLNode( child, xPath, linkData )
            else :
                raise ValueError( 'unsupported label = "%s"' % label )
        component = LLNLAngularEnergyForm( element.get( 'label' ), element.get( 'productFrame' ), angularSubform, angularEnergySubform )
        xPath.pop( )
        return( component )

class form(  baseModule.form ) :

    moniker = 'angularEnergy'
    subformAttributes = ( 'angularEnergySubform', )

    def __init__( self, label, productFrame, angularEnergySubform, makeCopy = True ) :

        if( not( isinstance( angularEnergySubform, subform ) ) ) : raise TypeError( 'instance is not an angularEnergy subform' )
        if( makeCopy ) : angularEnergySubform = angularEnergySubform.copy( )
        baseModule.form.__init__( self, label, productFrame, ( angularEnergySubform, ) )

    def calculateDepositionData( self, processInfo, tempInfo, verbosityIndent ) :

        return( self.angularEnergySubform.calculateDepositionData( processInfo, tempInfo, verbosityIndent ) )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        return( self.angularEnergySubform.process( processInfo, tempInfo, verbosityIndent ) )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( element.tag )
        formElement = element[0]
        subformClass = {
            pointwise.moniker:  pointwise,
        }[ formElement.tag ]
        if subformClass is None: raise Exception("encountered unknown angularEnergy subform: %s" % formElement.tag)
        subform = pointwise.parseXMLNode( formElement, xPath, linkData )
        AEC = form( element.get( 'label' ), element.get('productFrame'), subform )
        xPath.pop()
        return AEC
