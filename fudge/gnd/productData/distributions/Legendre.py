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
This is for ENDL I = 4 data with l > 0.
This is a temporary module, to be removed once testing is done and all coefficients in endl99, H2(n,2n)p data are filled in.
"""

import math

from xData import standards as standardsModule
from xData import axes as axesModule
from xData import XYs as XYsModule
from xData import multiD_XYs as multiD_XYsModule

from fudge.core.utilities import brb

from fudge.gnd import tokens as tokensModule

from fudge.processing import group as groupModule

from fudge.gnd.productData import multiplicity as multiplicityModule
from fudge.gnd.productData import energyDeposition as energyDepositionModule
from fudge.gnd.productData import momentumDeposition as momentumDepositionModule

from . import base as baseModule
from . import miscellaneous as miscellaneousModule
from . import energyAngular as energyAngularModule

__metaclass__ = type

def defaultAxes( energyUnit ) :

    axes = axesModule.axes( rank = 4 )
    axes[3] = axesModule.axis( 'l',          3, '' )
    axes[2] = axesModule.axis( 'energy_in',  2, energyUnit )
    axes[1] = axesModule.axis( 'energy_out', 1, energyUnit )
    axes[0] = axesModule.axis( 'c_l',        0, '1/' + energyUnit )
    return( axes )

class subform( baseModule.subform ) :
    """Abstract base class for LLNLLegendre forms."""

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

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        verbosity = kwargs.get( 'verbosity', 0 )
        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        energyUnit = kwargs['incidentEnergyUnit']
        momentumDepositionUnit = kwargs['momentumDepositionUnit']
        multiplicity = kwargs['multiplicity']
        productMass = kwargs['productMass']
        energyAccuracy = kwargs['energyAccuracy']
        momentumAccuracy = kwargs['momentumAccuracy']
        EMin = kwargs['EMin']

        Legendre_l0, Legendre_l1 = self[0], None
        if( len( self ) > 1 ) : Legendre_l1 = self[1]

        if( isinstance( multiplicity, multiplicityModule.XYs1d ) ) :    # If multiplicity as points not in Legendre_l0 add them.
            Es = set( [ EpP.value for EpP in Legendre_l0 ] )
            for E, m in multiplicity : Es.add( E )
            if( len( Es ) > len( Legendre_l0 ) ) :
                Es = sorted( Es )

                Legendre_l0p = multiD_XYsModule.XYs2d( )
                for energy in Es :
                    Legendre_l0p.append( Legendre_l0.interpolateAtValue( energy, unitBase = True, 
                            extrapolation = standardsModule.flatExtrapolationToken ) )
                Legendre_l0 = Legendre_l0p

                if( Legendre_l1 is not None ) :
                    Legendre_l1p = multiD_XYsModule.XYs2d( )
                    for energy in Es :
                        Legendre_l1p.append( Legendre_l1.interpolateAtValue( energy, unitBase = True,
                                extrapolation = standardsModule.flatExtrapolationToken ) )
                    Legendre_l1 = Legendre_l1p

        calculateDepositionEnergyFromEpP = miscellaneousModule.calculateDepositionEnergyFromEpP
        depEnergy = [ [ EpP.value, multiplicity.evaluate( EpP.value ) * calculateDepositionEnergyFromEpP( EpP.value, EpP ) ] for EpP in Legendre_l0 ]
        if( depEnergy[0][0] > EMin ) : depEnergy.insert( 0, [ EMin, 0. ] ) # Special case for bad data
        axes = energyDepositionModule.defaultAxes( energyUnit )
        energyDep = energyDepositionModule.XYs1d( data = depEnergy, axes = axes, label = style.label )

        if( Legendre_l1 is not None ) :
            depMomentum = []
            const = math.sqrt( 2. * productMass )
            for EpP in Legendre_l1 :
                if( const == 0 ) :                          # For gammas.
                    depMomentum.append( [ EpP.value, multiplicity.evaluate( EpP.value ) * EpP.integrateWithWeight_x( ) ] )
                else :
                    depMomentum.append( [ EpP.value, const * multiplicity.evaluate( EpP.value ) * EpP.integrateWithWeight_sqrt_x( ) ] )
        else :
            depMomentum = [ [ Legendre_l0[0].value, 0. ], [ Legendre_l0[-1].value, 0. ] ]
        axes = momentumDepositionModule.defaultAxes( energyUnit, momentumDepositionUnit )
        if( depMomentum[0][0] > EMin ) : depMomentum.insert( 0, [ EMin, 0. ] ) # Special case for bad data
        momentumDep = momentumDepositionModule.XYs1d( data = depMomentum, axes = axes, label = style.label )

        return( [ energyDep ], [ momentumDep ] )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.axes.convertUnits( unitMap )
        for data in self.data : data.convertUnits( unitMap )

    def processMultiGroup( self, style, tempInfo, indent ) :

        from fudge.processing.deterministic import transferMatrices as transferMatricesModule

        verbosity = tempInfo['verbosity']
        if( verbosity > 2 ) : print '%sGrouping %s' % ( indent, self.moniker )

        TM_1, TM_E = transferMatricesModule.ELEpP_TransferMatrix( style, tempInfo, tempInfo['crossSection'], tempInfo['productFrame'],
            self.data, tempInfo['multiplicity'], comment = tempInfo['transferMatrixComment'] + ' outgoing data for %s' % tempInfo['productLabel'] )

        return( groupModule.TMs2Form( style, tempInfo, TM_1, TM_E ) )

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :

        linear = self.toPointwise_withLinearXYs( )
        return( linear.to_xs_pdf_cdf1d( style, tempInfo, indent ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        energy_ins = {}
        for order in self.data :
            for EEpCl in order :
                if( EEpCl.value not in energy_ins ) : energy_ins[EEpCl.value] = set( )
                E_outs = energy_ins[EEpCl.value]
                for Ep, Cl in EEpCl : E_outs.add( Ep )

        energy_ins_sorted = sorted( energy_ins )
        for energy_in in energy_ins_sorted : energy_ins[energy_in] = sorted( energy_ins[energy_in] )

        axes = energyAngularModule.defaultAxes( self.axes[-2].unit )

        XYs3d = energyAngularModule.XYs3d( axes = axes )
        for i1, energy_in in enumerate( energy_ins_sorted ) :
            XYs2d = energyAngularModule.XYs2d( axes = axes, value = energy_in )

            energy_outs = energy_ins[energy_in]
            for energy_out in energy_outs :
                coefficients = []
                for order in self.data :
                    value = order.evaluate( energy_in ).evaluate( energy_out )
                    if( value is None ) : value = 0                             # This is a Kludge for W data.
                    coefficients.append( value )
                XYs2d.append( energyAngularModule.Legendre( coefficients, axes = axes, value = energy_out ) )
            XYs3d.append( XYs2d )
        return( XYs3d )

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

class form( baseModule.form ) :
    """
    This is for ENDL I = 4 data with l > 0.
    This is a temporary class, to be removed once testing is done and all coefficients in endl99, H2(n,2n)p data are filled in.
    """

    moniker = 'LLNLLegendre'
    subformAttributes = ( 'Legendre', )

    def __init__( self, label, productFrame, LegendreSubform ) :

        if( not( isinstance( LegendreSubform, subform ) ) ) : raise TypeError( 'instance is not a Legendre subform' )
        baseModule.form.__init__( self, label, productFrame, ( LegendreSubform, ) )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        return( self.Legendre.calculateAverageProductData( style, indent = '', **kwargs ) )

    def processMultiGroup( self, style, tempInfo, indent ) :

        tempInfo['productFrame'] = self.productFrame
        return( self.Legendre.processMultiGroup( style, tempInfo, indent ) )

    def processMC( self, style, tempInfo, indent ) :

        from . import energyAngularMC as energyAngularMCModule

        energy, energyAngular = self.Legendre.to_xs_pdf_cdf1d( style, tempInfo, indent )
        return( energyAngularMCModule.form( style.label, self.productFrame, energy, energyAngular ) )        

    @staticmethod
    def parseXMLNode( LegendreElement, xPath, linkData ) :

        xPath.append( LegendreElement.tag )
        subformElement = LegendreElement[0]
        subformClass = {
                LLNLPointwise.moniker : LLNLPointwise,
            }.get( subformElement.tag )
        if( subformClass is None ) : raise Exception( "unknown Legendre subform: %s" % subformElement.tag )
        LegendreSubform = subformClass.parseXMLNode( subformElement, xPath, linkData )
        _form = form( LegendreElement.get( 'label' ), LegendreElement.get( 'productFrame' ), LegendreSubform )
        xPath.pop( )
        return( _form )
