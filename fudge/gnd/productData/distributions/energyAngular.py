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

"""Energy/angular double differential distribution classes."""

import math
import base, miscellaneous
import fudge
from fudge.core.utilities import brb

from pqu import PQU

import xData.base as xDataBaseModule
import xData.axes as axesModule
import xData.standards as standardsModule
import xData.series1d as series1dModule
import xData.XYs as XYsModule
import xData.multiD_XYs as multiD_XYsModule
import xData.regions as regionsModule

from . import base as baseModule

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

class Legendre( series1dModule.LegendreSeries ) :

    def averageMu( self ) :

        return( self.getCoefficientSafely( 1 ) )

    def integrate( self ) :

        return( self.getCoefficientSafely( 0 ) )

class XYs2d( multiD_XYsModule.XYs2d ) :

    def averageEnergy( self ) :

        integralOfMu = [ [ pdfOfMuAtEp.value, pdfOfMuAtEp.integrate( ) ] for pdfOfMuAtEp in self ]
        return( float( XYsModule.XYs1d( integralOfMu, interpolation = self.interpolation ).integrateWithWeight_x( ) ) )

    def averageForwardMomentum( self, sqrt_2_ProductMass ) :

        averageMu = [ [ pdfOfMuAtEp.value, pdfOfMuAtEp.averageMu( ) ] for pdfOfMuAtEp in self ]
        return( sqrt_2_ProductMass * float( XYsModule.XYs1d( averageMu, interpolation = self.interpolation ).integrateWithWeight_sqrt_x( ) ) )

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
                energy = pdfOfEpMuAtE.value
                aveEnergy.append( [ energy, multiplicity.evaluate( energy ) * pdfOfEpMuAtE.averageEnergy( ) ] )

            for pdfOfEpMuAtE in self :
                energy = pdfOfEpMuAtE.value
                momemtum = multiplicity.evaluate( energy ) * pdfOfEpMuAtE.averageForwardMomentum( sqrt_2_ProductMass )
                if( momemtum < 1e-12 ) : momemtum = 0.          # This should be less arbitrary????????
                aveMomentum.append( [ energy, momemtum ] )

        else :                              # Center-of-mass.
            const1 = math.sqrt( 2 * projectileMass ) / ( projectileMass + targetMass )
            for pdfOfEpMuAtE in self :
                energy = pdfOfEpMuAtE.value
                vCOM = const1 * math.sqrt( energy )
                EpCOM = pdfOfEpMuAtE.averageEnergy( )
                ppCOM = pdfOfEpMuAtE.averageForwardMomentum( sqrt_2_ProductMass )
                productLabEnergy = 0.5 * productMass * vCOM * vCOM + EpCOM + vCOM * ppCOM
                aveEnergy.append( [ energy, multiplicity.evaluate( energy ) * productLabEnergy ] )

            for pdfOfEpMuAtE in self :
                energy = pdfOfEpMuAtE.value
                vCOM = const1 * math.sqrt( energy )
                productLabForwardMomentum = productMass * vCOM + pdfOfEpMuAtE.averageForwardMomentum( sqrt_2_ProductMass )
                aveMomentum.append( [ energy, multiplicity.evaluate( energy ) * productLabForwardMomentum ] )

        return( [ aveEnergy ], [ aveMomentum ] )

    def check( self, info ) :
        from fudge.gnd import warning
        from pqu import PQU
        warnings = []
        axes = axesModule.axes()
        axes[0] = self.axes[-2]
        axes[1] = self.axes[-1]

        for index, energy_in in enumerate(self):
            integral = float( XYsModule.XYs1d( [ ( eout.value, eout.coefficients[0] ) for eout in energy_in ], axes = axes ).integrate( ) )
            if abs(integral-1.0) > info['normTolerance']:
                warnings.append( warning.unnormalizedDistribution( PQU.PQU(energy_in.value,axes[0].unit),
                    index, integral, self.toXLink() ) )
            minProb = min( [xy.toPointwise_withLinearXYs( accuracy = XYsModule.defaultAccuracy, upperEps = 1e-8 ).rangeMin for xy in energy_in] )
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

    def processMultiGroup( self, style, tempInfo, indent ) :

        from fudge.processing import group as groupModule
        from fudge.processing.deterministic import transferMatrices as transferMatricesModule

        verbosity = tempInfo['verbosity']
        productFrame = tempInfo['productFrame']

        if( verbosity > 2 ) : print '%sGrouping %s' % ( indent, self.moniker )
        TM_1, TM_E = transferMatricesModule.EEpMuP_TransferMatrix( style, tempInfo, productFrame, tempInfo['crossSection'], self,
            tempInfo['multiplicity'], comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )

        return( groupModule.TMs2Form( style, tempInfo, TM_1, TM_E ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        return( multiD_XYsModule.XYs3d.toPointwise_withLinearXYs( self, cls = XYs3d, **kwargs ) )

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :

        from . import energy as energyModule
        from . import energyAngularMC as energyAngularMCModule

        energy = energyModule.XYs2d( axes = self.axes )
        for PofMuGivenEp in self :
            data = energyModule.XYs1d( [ [ PofMu.value, PofMu.integrate( ) ] for PofMu in PofMuGivenEp ] )
            energy.append( energyModule.xs_pdf_cdf1d.fromXYs( energyModule.XYs1d( data ), 
                    value = PofMuGivenEp.value ) )

        xys3d = energyAngularMCModule.XYs3d( axes = self.axes )
        for PofMuGivenEp in self :
            xys2d = energyAngularMCModule.XYs2d( value = PofMuGivenEp.value )
            for PofMu in PofMuGivenEp :
                _PofMu = PofMu
                if( isinstance( PofMu, Legendre ) ) :
                    if( PofMu[0] == 0 ) : _PofMu = XYs1d( [ [ -1, 0.5 ], [ 1, 0.5 ] ] )
                    _PofMu = _PofMu.toPointwise_withLinearXYs( accuracy = XYsModule.defaultAccuracy, upperEps = 1e-8 )
                xys2d.append( energyAngularMCModule.xs_pdf_cdf1d.fromXYs( _PofMu, PofMu.value ) )
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

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        kwargs['productFrame'] = self.productFrame
        return( self.energyAngularSubform.calculateAverageProductData( style, indent = indent, **kwargs ) )

    def processMC( self, style, tempInfo, indent ) :

        from . import energyAngularMC as energyAngularMCModule

        energy, energyAngular = self.energyAngularSubform.to_xs_pdf_cdf1d( style, tempInfo, indent )
        return( energyAngularMCModule.form( style.label, self.productFrame, energy, energyAngular ) )

    def processMultiGroup( self, style, tempInfo, indent ) :

        tempInfo['productFrame'] = self.productFrame
        return( self.energyAngularSubform.processMultiGroup( style, tempInfo, indent ) )

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
