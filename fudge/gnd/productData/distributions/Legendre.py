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

""" angular distribution classes """

import math
import miscellaneous
from fudge.core.utilities import brb
from fudge.gnd import tokens as tokensModule

from xData import axes as axesModule
from xData import XYs as XYsModule
from xData import multiD_XYs as multiD_XYsModule

from fudge.gnd.productData import multiplicity as multiplicityModule
from fudge.gnd.productData import energyDeposition as energyDepositionModule
from fudge.gnd.productData import momentumDeposition as momentumDepositionModule

from . import base as baseModule

noExtrapolationToken = 'noExtrapolation'

__metaclass__ = type

class subform( baseModule.subform ) :
    """Abstract base class for angular forms."""

    def __init__( self ) :

        baseModule.subform.__init__( self )

class LLNLPointwise( subform ) :
    """
    This is for ENDL I = 4 data with l > 0.
    This is a temporary class, to be removed once testing is done and all coefficients in endl99, H2(n,2n)p data are filled in.
    """

    moniker = 'LLNLLegendrePointwise'

    def __init__( self, axes ) :

        subform.__init__( self )
        self.data = []
        self.axes = axes
        self.label = None

    def __getitem__( self, l ) :

        return( self.data[l] )

    def __len__( self ) :

        return( len( self.data ) )

    def append( self, EEpP ) :

        if( not( isinstance( EEpP, multiD_XYsModule.XYs2d ) ) ) : raise Exception( 'EEpP is an instance of %s' % brb.getType( EEpP ) )
        self.data.append( EEpP )

    def calculateDepositionData( self, style, indent = '', **kwargs ) :

        verbosity = kwargs.get( 'verbosity', 0 )
        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        energyUnit = kwargs['incidentEnergyUnit']
        momentumDepositionUnit = energyUnit + '/c'
        massUnit = energyUnit + '/c**2'
        multiplicity = kwargs['multiplicity']
        productMass = kwargs['product'].getMass( massUnit )
        energyAccuracy = kwargs['energyAccuracy']
        momentumAccuracy = kwargs['momentumAccuracy']
        EMin = kwargs['EMin']

        depData = []

        Legendre_l0, Legendre_l1 = self[0], None
        if( len( self ) > 1 ) : Legendre_l1 = self[1]

        if( isinstance( multiplicity, multiplicityModule.XYs1d ) ) :
            Es = [ EpP.value for EpP in Legendre_l0 ]
            doIt = False
            for E, m in multiplicity.getFormByToken( tokensModule.pointwiseFormToken ) :
                if( E < Es[0] ) : continue
                if( E not in Es ) :
                    doIt = True
                    Es.append( E )
            if( doIt ) :
                Es.sort( )
                Legendre_l0p, Legendre_l1p = LLNLPointwiseEEpP( 0 ), LLNLPointwiseEEpP( 1 )
                indexE, EpP1 = 0, None
                for EpP2 in Legendre_l0 :
                    while( ( indexE < len( Es ) ) and ( EpP2.value > Es[indexE] ) ) :
                        if( EpP1.value != Es[indexE] ) :
                            EpP = XYs.pointwiseXY_C.unitbaseInterpolate( Es[indexE], EpP1.value, EpP1, EpP2.value, EpP2 )
                            axes = EpP1.axes.copy( )
                            Legendre_l0p.append( XYs.XYs( axes, EpP, EpP1.getAccuracy( ), value = Es[indexE] ) )
                        indexE += 1
                    Legendre_l0p.append( EpP2 )
                    EpP1 = EpP2
                Legendre_l0 = Legendre_l0p
                if( Legendre_l1 is not None ) :
                    for EpP2 in Legendre_l1 :
                        while( ( indexE < len( Es ) ) and ( EpP2.value > Es[indexE] ) ) :
                            if( EpP1 is None ) : raise Exception( 'multiplicity energy = %s before Legendre l = 0 distribution energy = %s' 
                                % ( Es[indexE], EpP2.value ) )
                            if( EpP1.value != Es[indexE] ) :
                                EpP = XYs.pointwiseXY_C.unitbaseInterpolate( Es[indexE], EpP1.value, EpP1, EpP2.value, EpP2 )
                                axes = EpP1.axes.copy( )
                                Legendre_l1p.append( XYs.XYs( axes, EpP, EpP1.getAccuracy( ), value = Es[indexE] ) )
                            indexE += 1
                        Legendre_l1p.append( EpP2 )
                        EpP1 = EpP2
                    Legendre_l1 = Legendre_l1p

        depEnergy = [ [ EpP.value, multiplicity.getValue( EpP.value ) * miscellaneous.calculateDepositionEnergyFromEpP( EpP.value, EpP ) ] for EpP in Legendre_l0 ]
        if( depEnergy[0][0] > EMin ) : depEnergy.insert( 0, [ 'EMin', 0. ] ) # Special case for bad data

        axes = energyDepositionModule.XYs1d.defaultAxes( energyUnit = energyUnit, energyDepositionUnit = energyUnit )
        depData.append( energyDepositionModule.XYs1d( data = depEnergy, axes = axes,
                label = style.label, accuracy = energyAccuracy ) )

        if( Legendre_l1 is not None ) :
            const = math.sqrt( 2. * productMass )
            if( const == 0. ) : const = 1               # For gammas.
            depMomentum = []
            EpP = Legendre_l1[0]
            for EpP in Legendre_l1 :
                if( const == 0. ) :                     # For gammas.
                    depMomentum.append( [ EpP.value, const * multiplicity.getValue( EpP.value ) * EpP.integrateWithWeight_x( ) ] )
                else :
                    depMomentum.append( [ EpP.value, const * multiplicity.getValue( EpP.value ) * EpP.integrateWithWeight_sqrt_x( ) ] )
        else :
            depMomentum = [ [ Legendre_l0[0].value, 0. ], [ Legendre_l0[-1].value, 0. ] ]
        axes = momentumDepositionModule.XYs1d.defaultAxes( energyUnit = energyUnit, momentumDepositionUnit = momentumDepositionUnit )
        if( depMomentum[0][0] > EMin ) : depMomentum.insert( 0, [ EMin, 0. ] ) # Special case for bad data
        depData.append( momentumDepositionModule.XYs1d( data = depMomentum, axes = axes,
                label = style.label, accuracy = momentumAccuracy ) )
        return( depData )

    def processSnMultiGroup( self, style, tempInfo, indent ) :

        raise Exception( 'need to implement' )

#    def process( self, processInfo, tempInfo, verbosityIndent ) :

        import energy
        from fudge.processing.deterministic import transferMatrices

        newComponents = []
        Legendre = self.data

        if( 'LLNL_Pn' in processInfo['styles'] ) :
            if( processInfo.verbosity >= 30 ) : print '%sGrouping %s' % ( verbosityIndent, self.moniker )
            outputChannel = tempInfo['reaction'].outputChannel
            projectile, product = processInfo.getProjectileName( ), tempInfo['product'].particle.name
            TM_1, TM_E = transferMatrices.ELEpP_TransferMatrix( processInfo, projectile, product, tempInfo['crossSection'],
                Legendre, tempInfo['multiplicity'], comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )
            fudge.gnd.miscellaneous.TMs2Form( processInfo, tempInfo, newComponents, TM_1, TM_E, self.axes )

        return( newComponents )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ self.XMLStartTagString( indent = indent ) ]
        xmlString += self.axes.toXMLList( indent2, **kwargs )
        for l in self : xmlString += l.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )
        axes = axesModule.axes.parseXMLNode( element.find( axesModule.axes.moniker ), xPath, linkData )
        lpw = LLNLPointwise( axes )
        for lOrderData in element :
            if( lOrderData.tag == axesModule.axes.moniker ) : continue
            lpw.append( multiD_XYsModule.XYs2d.parseXMLNode( lOrderData, xPath, linkData ) )
        xPath.pop( )
        return( lpw )

    @staticmethod
    def defaultAxes( ) :

        axes = axesModule.axes( rank = 4 )
        axes[0] = axesModule.axis( 'c_l',        0, '1/MeV' )
        axes[1] = axesModule.axis( 'energy_out', 1, 'MeV' )
        axes[2] = axesModule.axis( 'energy_in',  2, 'MeV' )
        axes[3] = axesModule.axis( 'l',          3, '' )
        return( axes )

class LLNLPointwiseEEpP :

    def __init__( self, l ) :

        raise "This should no longer be used, but should not be removed until LLNLPointwise.calculateAverageProductData is changed as it uses it."
        self.l = l
        self.EpP = []

    def __getitem__( self, index ) :

        return( self.EpP[index] )

    def __len__( self ) :

        return( len( self.EpP ) )

    def append( self, EpP ) :

        if( not( isinstance( EpP, XYsModule.XYs1d ) ) ) : raise Exception( 'EpP is an instance of %s' % brb.getType( EpP ) )
        if( len( self ) > 0 ) :
            if( EpP.value is None ) : raise Exception( "EpP's value is None" )
            if( EpP.value <= self.EpP[-1].value ) : raise Exception( "EpP.value = %s <= prior's values = %s " % ( EpP.value, self.EpP[-1].value ) )
        self.EpP.append( EpP )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        kwargs['oneLine'] = True

        xmlString = [ '%s<l index="%s">' % ( indent, self.l ) ]
        for energy in self : xmlString += energy.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</l>'
        return( xmlString )

class form( baseModule.form ) :

    moniker = 'Legendre'
    subformAttributes = ( 'LegendreSubform', )

    def __init__( self, label, productFrame, LegendreSubform ) :

        if( not( isinstance( LegendreSubform, subform ) ) ) : raise TypeError( 'instance is not a Legendre subform' )
        baseModule.form.__init__( self, label, productFrame, ( LegendreSubform, ) )

    def calculateDepositionData( self, style, indent = '', **kwargs ) :

        raise 'hell - I do not think this is used anymore - BRB'
        return( self.forms[0].calculateDepositionData( style, indent = '', **kwargs ) )

    def processSnMultiGroup( self, style, tempInfo, indent ) :

        raise Exception( 'need to implement' )

#    def process( self, processInfo, tempInfo, verbosityIndent ) :

        return( self.forms[0].process( processInfo, tempInfo, verbosityIndent ) )

    @staticmethod
    def parseXMLNode( LegendreElement, xPath, linkData ) :

        xPath.append( LegendreElement.tag )
        subformElement = LegendreElement[0]
        subformClass = {
                LLNLPointwise.moniker : LLNLPointwise,
            }.get( subformElement.tag )
        if( subformClass is None ) : raise Exception( "unknown Legendre subform: %s" % subformElement.tag )
        LegendreSubform = subformClass.parseXMLNode( subformElement, xPath, linkData )
        form_ = form( LegendreElement.get( 'label' ), LegendreElement.get( 'productFrame' ), LegendreSubform )
        xPath.pop( )
        return( form_ )
