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
import base, miscellaneous
import fudge
from fudge.core.utilities import brb

from pqu import PQU

import xData.base as xDataBaseModule
import xData.axes as axesModule
import xData.XYs as XYsModule
import xData.standards as standardsModule
import xData.series1d as series1dModule
import xData.multiD_XYs as multiD_XYsModule
import xData.regions as regionsModule

from fudge.gnd.productData import energyDeposition as energyDepositionModule
from fudge.gnd.productData import momentumDeposition as momentumDepositionModule

from . import base as baseModule

__metaclass__ = type

class pdfOfMu :             # FIXME: BRB, Should this class inherit from angular.pdfOfMu?

    class pointwise( XYsModule.XYs ) :

        def averageMu( self ) :

            allowedInterpolations = [ standardsModule.interpolation.linlinToken,
                                      standardsModule.interpolation.flatToken ]
            xys = self.changeInterpolationIfNeeded( allowedInterpolations = allowedInterpolations )
            return( xys.integrateWithWeight_x( ) )

    class Legendre( series1dModule.LegendreSeries ) :

        def averageMu( self ) :

            return( self.getCoefficientSafely( 1 ) )

        def integrate( self ) :

            return( self.getCoefficientSafely( 0 ) )

class pdfOfEpAndMu :

    class pointwise( multiD_XYsModule.multiD_XYs ) :

        def averageEnergy( self ) :

            integralOfMu = [ [ pdfOfMuAtEp.value, pdfOfMuAtEp.integrate( ) ] for pdfOfMuAtEp in self ]
            return( float( XYsModule.XYs( integralOfMu, interpolation = self.interpolation ).integrateWithWeight_x( ) ) )

        def averageMomentum( self ) :

            integralOfMu = [ [ pdfOfMuAtEp.value, pdfOfMuAtEp.integrate( ) ] for pdfOfMuAtEp in self ]
            return( float( XYsModule.XYs( integralOfMu, interpolation = self.interpolation ).integrateWithWeight_sqrt_x( ) ) )

        @staticmethod
        def allowedSubElements( ) :

            return( ( pdfOfMu.pointwise, pdfOfMu.Legendre ) )

class subform( baseModule.subform ) :
    """Abstract base class for energyAngular subforms."""

    pass

class pointwise( subform, multiD_XYsModule.multiD_XYs ) :

    def __init__( self, **kwargs ) :

        if( 'dimension' not in kwargs ) : kwargs['dimension'] = 3
        if( kwargs['dimension'] != 3 ) : raise ValueError( 'Dimension = %s != 3' % ( kwargs['dimension'] ) )
        multiD_XYsModule.multiD_XYs.__init__( self, **kwargs )
        subform.__init__( self )

    def calculateDepositionData( self, processInfo, tempInfo, verbosityIndent ) :

        if( self.productFrame != standardsModule.frames.labToken ) :
            root = self.getRootAncestor( )
            f = open( 'energyAngular.COM', 'a' )
            f.write( '%s\n' % root.inputParticlesToReactionString( ) )
            f.close( )
            return( [] )

        energyUnit = tempInfo['incidentEnergyUnit']
        momentumDepositionUnit = energyUnit + '/c'
        massUnit = energyUnit + '/c**2'
        energyAccuracy, momentumAccuracy = processInfo.energyAccuracy, processInfo.momentumAccuracy
        multiplicity = tempInfo['multiplicity']
        productMass = tempInfo['product'].getMass( massUnit )

        depEnergy = []

        for pdfOfEpMuAtE in self :
            energy = pdfOfEpMuAtE.value
            depEnergy.append( [ energy, multiplicity.getValue( energy ) * pdfOfEpMuAtE.averageEnergy( ) ] )
        axes = energyDepositionModule.pointwise.defaultAxes( energyUnit = energyUnit, 
                energyDepositionUnit = energyUnit )
        depEnergy = energyDepositionModule.pointwise( data = depEnergy, axes = axes,
                label = processInfo.style.label, accuracy = energyAccuracy )

        const = math.sqrt( 2. * productMass )
        depMomentum = []
        for pdfOfEpMuAtE in self :
            energy = pdfOfEpMuAtE.value
            momemtum = const * multiplicity.getValue( energy ) * pdfOfEpMuAtE.averageMomentum( )
            if( momemtum < 1e-12 ) : momemtum = 0.          # This should be less arbitrary????????
            depMomentum.append( [ energy, momemtum ] )
        axes = momentumDepositionModule.pointwise.defaultAxes( energyUnit = energyUnit, 
                momentumDepositionUnit = momentumDepositionUnit )
        depMomentum = momentumDepositionModule.pointwise( data = depMomentum, axes = axes,
                label = processInfo.style.label, accuracy = momentumAccuracy )

        return( [ depEnergy, depMomentum ] )

    def check( self, info ) :
        from fudge.gnd import warning
        from pqu import PQU
        warnings = []
        axes = axesModule.axes()
        axes[0] = self.axes[-2]
        axes[1] = self.axes[-1]

        for index, energy_in in enumerate(self):
            integral = float( XYsModule.XYs( [ ( eout.value, eout.coefficients[0] ) for eout in energy_in ],
                    axes = axes, accuracy = 0.001 ).integrate( ) )
            if abs(integral-1.0) > info['normTolerance']:
                warnings.append( warning.unnormalizedDistribution( PQU.PQU(energy_in.value,axes[0].unit),
                    index, integral, self.toXLink() ) )
            minProb = min( [xy.toPointwise_withLinearXYs( 0.001 ).rangeMin() for xy in energy_in] )
            if minProb < 0:
                warnings.append( warning.negativeProbability( PQU.PQU(energy_in.value,axes[0].unit),
                    value = minProb, obj = energy_in ) )
            if self.interpolationQualifier is standardsModule.interpolation.unitBaseToken:
                # check for more than one outgoing distribution integrating to 0 at high end of incident energies
                integrals = [eout.integrate() for eout in energy_in]
                for idx, integral in enumerate(integrals[::-1]):
                    if integral != 0: break
                if idx > 1:
                    warnings.append( warning.extraOutgoingEnergy( PQU.PQU(energy_in.value,axes[-1].unit),
                        obj = energy_in ) )
        return warnings

    def normalize( self, insitu = True ) :

        n = self
        if( not( insitu ) ) : n = self.copy( )
        for E_MuEpPs in n :
            sum = E_MuEpPs.integrate( )
            for muEpPs in E_MuEpPs : muEpPs.setData( muEpPs / sum )
        return( n )

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        return( multiD_XYsModule.multiD_XYs.toPointwise_withLinearXYs( self, accuracy, lowerEps = lowerEps, upperEps = upperEps, cls = pointwise ) )

    @staticmethod
    def allowedSubElements( ) :

        return( ( pdfOfEpAndMu.pointwise, ) )

    @staticmethod
    def defaultAxes( energyUnit = 'eV', energy_outUnit = 'eV', probabilityUnit = '1/eV' ) :

        axes = axesModule.axes( rank = 4 )
        axes[0] = axesModule.axis( 'P(mu,energy_out|energy_in)', 0, probabilityUnit )
        axes[1] = axesModule.axis( 'mu', 1, '' )
        axes[2] = axesModule.axis( 'energy_out', 2, energy_outUnit )
        axes[3] = axesModule.axis( 'energy_in', 3, energyUnit )
        return( axes )

class piecewise( subform, regionsModule.regions ) :

    def __init__( self, **kwargs ) :

        if( 'dimension' not in kwargs ) : kwargs['dimension'] = 3
        if( kwargs['dimension'] != 3 ) : raise ValueError( 'Dimension = %s != 3' % ( kwargs['dimension'] ) )
        regionsModule.regions.__init__( self, **kwargs )
        subform.__init__( self )

    def calculateDepositionData( self, processInfo, tempInfo, verbosityIndent ) :

        raise Exception( 'Not implemented' )

    def check( self, info ) :

        raise Exception( 'Not implemented' )

    def normalize( self, insitu = True ) :

        n = self
        if( not( insitu ) ) : n = self.copy( )
        for region in n : region.normalize( insitu = True )
        return( n )

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        raise Exception( 'Not implemented' )

    @staticmethod
    def allowedSubElements( ) :

        return( ( pointwise, ) )

class form( baseModule.form ) :

    moniker = 'energyAngular'
    subformAttributes = ( 'energyAngularSubform', )

    def __init__( self, label, productFrame, energyAngularSubform, makeCopy = True ) :

        if( not( isinstance( energyAngularSubform, subform ) ) ) : raise TypeError( 'instance is not an energyAngular subform' )
        if( makeCopy ) : energyAngularSubform = energyAngularSubform.copy( )
        baseModule.form.__init__( self, label, productFrame, ( energyAngularSubform, ) )

    def calculateDepositionData( self, processInfo, tempInfo, verbosityIndent ) :

        return( self.energyAngularSubform.calculateDepositionData( processInfo, tempInfo, verbosityIndent ) )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        return( self.energyAngularSubform.process( processInfo, tempInfo, verbosityIndent ) )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):
        """ translate <energyAngular> element from xml """

        xPath.append( element.tag )
        subformElement = element[0]
        subformClass = {
                pointwise.moniker: pointwise,
                }.get( subformElement.tag )
        if subformClass is None: raise Exception( "encountered unknown energyAngular subform: %s" % subformElement.tag )
        subForm = subformClass.parseXMLNode( subformElement, xPath, linkData )
        energyAngular = form( element.get( "label" ), element.get('productFrame'), subForm )
        xPath.pop()
        return energyAngular
