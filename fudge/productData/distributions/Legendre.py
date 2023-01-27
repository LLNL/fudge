# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This is for ENDL I = 4 data with l > 0.
This is a temporary module, to be removed once testing is done and all coefficients in endl99, H2(n,2n)p data are filled in.
"""

import math

from xData import enums as xDataEnumsModule
from xData import axes as axesModule
from xData import multiD_XYs as multiD_XYsModule
from xData import series1d as series1dModule

from fudge.core.utilities import brb

from fudge.processing import group as groupModule

from .. import averageProductEnergy as averageProductEnergyModule
from fudge.productData import averageProductMomentum as averageProductMomentumModule

from . import base as baseModule
from . import miscellaneous as miscellaneousModule
from . import energyAngular as energyAngularModule


def defaultAxes( energyUnit ) :

    axes = axesModule.Axes(4)
    axes[3] = axesModule.Axis( 'l',          3, '' )
    axes[2] = axesModule.Axis( 'energy_in',  2, energyUnit )
    axes[1] = axesModule.Axis( 'energy_out', 1, energyUnit )
    axes[0] = axesModule.Axis( 'c_l',        0, '1/' + energyUnit )
    return( axes )

class Subform( baseModule.Subform ) :
    """Abstract base class for LLNLLegendre forms."""

    def __init__( self ) :

        baseModule.Subform.__init__( self )

class LLNLPointwise( Subform ) :
    """
    This is for ENDL I = 4 data with l > 0.
    This is a temporary class, to be removed once testing is done and all coefficients in endl99, H2(n,2n)p data are filled in.
    """

    moniker = 'LLNLLegendrePointwise'

    def __init__( self, axes ) :

        Subform.__init__( self )
        self.data = []
        self.axes = axes
        self.label = None

    def __getitem__( self, l ) :

        return( self.data[l] )

    def __len__( self ) :

        return( len( self.data ) )

    @property
    def domainUnit( self ) :

        return( self.data[0].domainUnit )

    def append( self, EEpP ) :

        if( not( isinstance( EEpP, multiD_XYsModule.XYs2d ) ) ) : raise Exception( 'EEpP is an instance of %s' % brb.getType( EEpP ) )
        self.data.append( EEpP )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        from .. import multiplicity as multiplicityModule

        verbosity = kwargs.get( 'verbosity', 0 )
        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        energyUnit = kwargs['incidentEnergyUnit']
        momentumUnit = kwargs['momentumUnit']
        multiplicity = kwargs['multiplicity']
        productMass = kwargs['productMass']
        EMin = kwargs['EMin']

        Legendre_l0, Legendre_l1 = self[0], None
        if( len( self ) > 1 ) : Legendre_l1 = self[1]

        if( isinstance( multiplicity, multiplicityModule.XYs1d ) ) :    # If multiplicity as points not in Legendre_l0 add them.
            Es = set( [ EpP.outerDomainValue for EpP in Legendre_l0 ] )
            for E, m in multiplicity : Es.add( E )
            if( len( Es ) > len( Legendre_l0 ) ) :
                Es = sorted( Es )

                Legendre_l0p = multiD_XYsModule.XYs2d( )
                for energy in Es :
                    Legendre_l0p.append(Legendre_l0.interpolateAtValue(energy, unitBase=True, extrapolation=xDataEnumsModule.Extrapolation.flat))
                Legendre_l0 = Legendre_l0p

                if( Legendre_l1 is not None ) :
                    Legendre_l1p = multiD_XYsModule.XYs2d( )
                    for energy in Es :
                        Legendre_l1p.append(Legendre_l1.interpolateAtValue(energy, unitBase=True, extrapolation=xDataEnumsModule.Extrapolation.flat))
                    Legendre_l1 = Legendre_l1p

        calculateDepositionEnergyFromEpP = miscellaneousModule.calculateDepositionEnergyFromEpP
        depEnergy = [ [ EpP.outerDomainValue, multiplicity.evaluate( EpP.outerDomainValue ) * calculateDepositionEnergyFromEpP( EpP.outerDomainValue, EpP ) ] for EpP in Legendre_l0 ]
        if( depEnergy[0][0] > EMin ) : depEnergy.insert( 0, [ EMin, 0. ] ) # Special case for bad data
        axes = averageProductEnergyModule.defaultAxes( energyUnit )
        energyDep = averageProductEnergyModule.XYs1d( data = depEnergy, axes = axes, label = style.label )

        if( Legendre_l1 is not None ) :
            depMomentum = []
            const = math.sqrt( 2. * productMass )
            for EpP in Legendre_l1 :
                if( const == 0 ) :                          # For gammas.
                    depMomentum.append( [ EpP.outerDomainValue, multiplicity.evaluate( EpP.outerDomainValue ) * EpP.integrateWithWeight_x( ) ] )
                else :
                    depMomentum.append( [ EpP.outerDomainValue, const * multiplicity.evaluate( EpP.outerDomainValue ) * EpP.integrateWithWeight_sqrt_x( ) ] )
        else :
            depMomentum = [ [ Legendre_l0[0].outerDomainValue, 0. ], [ Legendre_l0[-1].outerDomainValue, 0. ] ]
        axes = averageProductMomentumModule.defaultAxes( energyUnit, momentumUnit )
        if( depMomentum[0][0] > EMin ) : depMomentum.insert( 0, [ EMin, 0. ] ) # Special case for bad data
        momentumDep = averageProductMomentumModule.XYs1d( data = depMomentum, axes = axes, label = style.label )

        return( [ energyDep ], [ momentumDep ] )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.axes.convertUnits( unitMap )
        for data in self.data : data.convertUnits( unitMap )

    def energySpectrumAtEnergy(self, energyIn, **kwargs):

        muMin = kwargs.get('muMin', -1.0)
        muMax = kwargs.get('muMax',  1.0)

        muFactor = 1.0
        if muMin != -1.0 and muMax != 1.0:
            if len(self.data) > 1:
                raise Exception('Energy spectrum only supported for full mu domain (i.e., -1 to 1) and not (% to %s).' %(muMin, muMax))
            muFactor = muMax - muMin

        return muFactor * self.data[0].evaluate(energyIn)

    def fixDomains(self, domainMin, domainMax, fixToDomain):
        """
        Calls the **fixDomains** for each entry in *self.data*.
        """

        numberOfFixes = 0
        for data in self.data:
            numberOfFixes += data.fixDomains(domainMin, domainMax, fixToDomain)

        return numberOfFixes

    def integrate( self, reaction_suite, energyIn, energyOut = None, muOut = None, phiOut = None, frame = xDataEnumsModule.Frame.product, LegendreOrder = 0 ) :

        if( frame != xDataEnumsModule.Frame.lab) : raise TypeError( "For LLNL Legendre data (i.e., I=4, l>0), only lab frame is supported." )

        if( muOut is None ) :
            print( 'integrate of LegendreOrder not supported', type( Form ) )
            return( 0.0 )

        C_ls = []
        for data in self.data :
            C_l = data.evaluate( energyIn )
            domainMin, domainMax = miscellaneousModule.domainLimits( energyOut, C_l.domainMin, C_l.domainMax )
            C_ls.append(C_l.integrate(domainMin, domainMax))

        series = series1dModule.LegendreSeries( C_ls )
        domainMin, domainMax = miscellaneousModule.domainLimits( muOut, -1.0, 1.0 )

        phi_evaluate = miscellaneousModule.muPhiEvaluate( None, phiOut )

        return(phi_evaluate * series.integrate(domainMin, domainMax))

    def processMultiGroup( self, style, tempInfo, indent ) :

        from fudge.processing.deterministic import transferMatrices as transferMatricesModule

        verbosity = tempInfo['verbosity']
        if( verbosity > 2 ) : print('%sGrouping %s' % (indent, self.moniker))

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
                if( EEpCl.outerDomainValue not in energy_ins ) : energy_ins[EEpCl.outerDomainValue] = set( )
                E_outs = energy_ins[EEpCl.outerDomainValue]
                for Ep, Cl in EEpCl : E_outs.add( Ep )

        energy_ins_sorted = sorted( energy_ins )
        for energy_in in energy_ins_sorted : energy_ins[energy_in] = sorted( energy_ins[energy_in] )

        axes = energyAngularModule.defaultAxes( self.axes[-2].unit )

        XYs3d = energyAngularModule.XYs3d( axes = axes )
        for i1, energy_in in enumerate( energy_ins_sorted ) :
            XYs2d = energyAngularModule.XYs2d( axes = axes, outerDomainValue = energy_in )

            energy_outs = energy_ins[energy_in]
            for energy_out in energy_outs :
                coefficients = []
                for order in self.data :
                    value = order.evaluate( energy_in ).evaluate( energy_out )
                    if( value is None ) : value = 0                             # This is a Kludge for W data.
                    coefficients.append( value )
                XYs2d.append( energyAngularModule.Legendre( coefficients, axes = axes, outerDomainValue = energy_out ) )
            XYs3d.append( XYs2d )
        return( XYs3d )

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ self.XMLStartTagString( indent = indent ) ]
        xmlString += self.axes.toXML_strList( indent2, **kwargs )
        for l in self : xmlString += l.toXML_strList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )
        axes = axesModule.Axes.parseNodeUsingClass(element.find(axesModule.Axes.moniker), xPath, linkData, **kwargs)
        lpw = cls( axes )
        for lOrderData in element :
            if( lOrderData.tag == axesModule.Axes.moniker ) : continue
            lpw.append(multiD_XYsModule.XYs2d.parseNodeUsingClass(lOrderData, xPath, linkData, **kwargs))
        xPath.pop( )
        return( lpw )

class Form( baseModule.Form ) :
    """
    This is for ENDL I = 4 data with l > 0.
    This is a temporary class, to be removed once testing is done and all coefficients in endl99, H2(n,2n)p data are filled in.
    """

    moniker = 'LLNLLegendre'
    subformAttributes = ( 'Legendre', )

    def __init__( self, label, productFrame, LegendreSubform ) :

        if( not( isinstance( LegendreSubform, Subform ) ) ) : raise TypeError( 'instance is not a Legendre subform' )
        baseModule.Form.__init__( self, label, productFrame, ( LegendreSubform, ) )

        for data in LegendreSubform.data :
            if data.interpolationQualifier == xDataEnumsModule.InterpolationQualifier.none:
                data.interpolationQualifier = xDataEnumsModule.InterpolationQualifier.unitBase

    @property
    def domainUnit( self ) :

        return( self.Legendre.domainUnit )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        return( self.Legendre.calculateAverageProductData( style, indent = '', **kwargs ) )

    def energySpectrumAtEnergy(self, energyIn, frame, **kwargs):

        if frame == xDataEnumsModule.Frame.centerOfMass:
            raise TypeError('Lab to center-of-mass translation not supported.')

        return self.Legendre.energySpectrumAtEnergy(energyIn, **kwargs)

    def fixDomains(self, domainMin, domainMax, fixToDomain):
        """
        Calls the **fixDomains** for the **Legendre** members.
        """

        return self.Legendre.fixDomains(domainMin, domainMax, fixToDomain)

    def integrate( self, reaction_suite, energyIn, energyOut = None, muOut = None, phiOut = None, frame = xDataEnumsModule.Frame.product, LegendreOrder = 0 ) :

        return( self.Legendre.integrate( reaction_suite, energyIn, energyOut = energyOut, muOut = muOut, phiOut = phiOut, frame = frame, LegendreOrder = LegendreOrder ) )

    def processMultiGroup( self, style, tempInfo, indent ) :

        tempInfo['productFrame'] = self.productFrame
        return( self.Legendre.processMultiGroup( style, tempInfo, indent ) )

    def processMC_cdf( self, style, tempInfo, indent ) :

        from . import energyAngularMC as energyAngularMCModule

        energy, energyAngular = self.Legendre.to_xs_pdf_cdf1d( style, tempInfo, indent )
        return( energyAngularMCModule.Form( style.label, self.productFrame, energy, energyAngular ) )        

    @classmethod
    def parseNodeUsingClass(cls, LegendreElement, xPath, linkData, **kwargs):

        xPath.append( LegendreElement.tag )
        subformElement = LegendreElement[0]
        subformClass = {
                LLNLPointwise.moniker : LLNLPointwise,
            }.get( subformElement.tag )
        if( subformClass is None ) : raise Exception( "unknown Legendre subform: %s" % subformElement.tag )
        LegendreSubform = subformClass.parseNodeUsingClass(subformElement, xPath, linkData, **kwargs)

        _form = cls( LegendreElement.get( 'label' ), LegendreElement.get( 'productFrame' ), LegendreSubform )

        xPath.pop( )
        return( _form )
