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
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# <<END-copyright>>

""" energy distribution (spectra) classes """

import math
import base, miscellaneous
import fudge
from fudge.core.math.xData import axes, XYs, W_XYs, regions
from pqu.physicalQuantityWithUncertainty import PhysicalQuantityWithUncertainty
from fudge.core.utilities import brb
from fudge.core import ancestry

__metaclass__ = type

class component( base.component ) :

    def __init__( self, nativeData = base.noneFormToken ) :

        base.component.__init__( self, base.energyComponentToken, nativeData )

    def getSpectrumAtEnergy( self, energy ) :
        """Returns the energy spectrum for self at projectile energy using nativeData."""

        return( self.getNativeData( ).getSpectrumAtEnergy( energy ) )

    def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

        from fudge.legacy.converting import endfFormats, gndToENDF6
        weight = targetInfo['delayedNubarWeight']
        if( hasattr( self.forms[self.nativeData], 'toENDF6' ) ) :
            NK, MF5 = self.forms[self.nativeData].toENDF6( flags, targetInfo, weight = weight )
            if( MT == 455 ) :
                endfMFList[5][MT].append( MF5 ) 
            else :
                endfMFList[5][MT] = [ endfFormats.endfHeadLine( targetInfo['ZA'], targetInfo['mass'], 0, 0, NK, 0 ) ] + MF5 + \
                    [ endfFormats.endfSENDLineNumber( ) ]
        else :
            print 'WARNING: form %s does not have method toENDF6 for component %s' % ( self.nativeData, self.moniker )

class form( base.form ) :
    """Abstract base class for energy forms."""

    pass

class pointwise( form, W_XYs.W_XYs ) :

    tag = base.pointwiseFormToken
    moniker = base.pointwiseFormToken

    def __init__( self, axes, productFrame, **kwargs ):
        """
        @param: axes is an instance of xData.axes.axes
        >epointwise = pointwise( axes, productFrame )
        followed by:
        >epointwise[ 0 ] = XYs_data_1
        >epointwise[ 1 ] = XYs_data_2
        > ...
        >epointwise[ n-1 ] = XYs_data_n
        """

        form.__init__( self, base.pointwiseFormToken, productFrame )
        kwargs['isPrimaryXData'] = True
        W_XYs.W_XYs.__init__( self, axes, **kwargs )

    def extraXMLAttributeString( self ) :

        return( 'productFrame="%s"' % self.productFrame )

    def getAtEnergy( self, energy ) :
        "This method is deprecated, use getSpectrumAtEnergy."

        return( self.getSpectrumAtEnergy( energy ) )

    def getSpectrumAtEnergy( self, energy ) :
        """Returns the energy spectrum for self at projectile energy."""

        unitBase = self.axes[0].interpolation.isQualifierUnitBase( ) or self.axes[0].interpolation.isQualifierCorrespondingPoints( )
        unitBase = True
        return( self.interpolateAtW( energy, unitBase = unitBase, extrapolation = W_XYs.flatExtrapolationToken ) )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        Es = [ data.value for data in self ]
        if( EMin is not None ) :
            if( EMin < Es[0] ) : Es.insert( 0, EMin )
        if( EMax is not None ) :
            if( EMax > Es[-1] ) : Es.append( EMax )
        return( Es )

    def EpAverageAtE( self, E ) :

        return( self.getAtEnergy( E ).integrateWithWeight_x( ) )

    def check( self, info ) :

        from fudge.gnd import warning

        warnings = []
        for i in range(len(self)):
            integral = self[i].integrate()
            if abs(integral - 1.0) > info['normTolerance']:
                warnings.append( warning.unnormalizedDistribution( PhysicalQuantityWithUncertainty(self[i].value,self.axes[0].unit),
                    i, integral, self[i] ) )

            if self[i].yMin() < 0.0:
                warnings.append( warning.negativeProbability( PhysicalQuantityWithUncertainty(self[i].value,self.axes[0].unit),
                    value=self[i].yMin(), obj=self[i] ) )

        return warnings

    def sqrtEp_AverageAtE( self, E ) :

        return( self.getAtEnergy( E ).integrateWithWeight_sqrt_x( ) )

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        return( W_XYs.W_XYs.toPointwise_withLinearXYs( self, accuracy, lowerEps = lowerEps,
            upperEps = upperEps, cls = linear ) )

    def toENDF6( self, flags, targetInfo, weight = None ) :

        from fudge.legacy.converting import endfFormats, gndToENDF6
        NE = len( self )
        if( weight is None ) : weight = [ [ self[0].value, 1.0 ], [ self[-1].value, 1.0 ] ]
        EInInterpolation = gndToENDF6.axisToEndfInterpolationFlag( self.axes[0] )
        ENDFDataList = [ endfFormats.endfContLine( 0, 0, 0, 1, 1, len( weight ) ) ] + \
            endfFormats.endfInterpolationList( [ len( weight ), 2 ] ) + endfFormats.endfNdDataList( weight ) + \
            [ endfFormats.endfContLine( 0, 0, 0, 0, 1, NE ) ] + endfFormats.endfInterpolationList( [ NE, EInInterpolation ] )
        independent, dependent, qualifier = self.axes[1].interpolation.getInterpolationTokens( )
        interpolation = gndToENDF6.gndToEndfInterpolationFlag( independent, dependent )
        for energy in self :
            ENDFDataList.append( endfFormats.endfContLine( 0, energy.value, 0, 0, 1, len( energy ) ) )
            ENDFDataList += endfFormats.endfInterpolationList( [ len( energy ), interpolation ] )
            ENDFDataList += endfFormats.endfNdDataList( energy, xUnit = 'eV', yUnit = '1/eV' )
        return( 1, ENDFDataList )

    @staticmethod
    def defaultAxes( energyUnit = 'eV', energyInterpolation = axes.linearToken, energyFunctionInterpolation = axes.linearToken, 
            energyInterpolationQualifier = axes.unitBaseToken, energy_outInterpolation = axes.linearToken, probabilityInterpolation = axes.linearToken ) :

        axes_ = axes.axes( dimension = 3 )
        axes_[0] = axes.axis( 'energy_in',  0, energyUnit, 
            interpolation = axes.interpolationXY( energyInterpolation, energyFunctionInterpolation, energyInterpolationQualifier ) )
        axes_[1] = axes.axis( 'energy_out', 1, energyUnit, interpolation = axes.interpolationXY( energy_outInterpolation, probabilityInterpolation ) )
        axes_[2] = axes.axis( 'P(energy_out|energy_in)', 2, '1/' + energyUnit )
        return( axes_ )

class linear( pointwise ) :

    tag = base.linearFormToken
    moniker = base.linearFormToken

    def __init__( self, axes, productFrame, **kwargs ) :

        if( not axes.isLinear( qualifierOk = True ) ) : raise Exception( 'interpolation not linear: %s' % axes )
        pointwise.__init__( self, axes, productFrame = productFrame, **kwargs )

class semiPiecewise( form ) :

    def __init__( self, axes_, productFrame ) :

        form.__init__( self, base.semiPiecewiseFormToken, productFrame )
        self.axes = axes_.copy( parent = self )
        self.energies = []

    def __len__( self ) :

        return( len( self.energies ) )

    def __getitem__( self, index ) :

        return( self.energies[index] )

    def append( self, energy, region ) :

        if( len( self.energies ) > 0 ) :
            if( self.energies[-1].energy >= energy ) : raise Exception( 'energy = %e not greater than %e' % ( energy, self.energies[-1].energy ) )
        self.energies.append( semiPiecewiseEnergy( len( self.energies ), energy, region, self ) )

    def check( self, info ):

        from fudge.gnd import warning

        warnings = []
        for i in range(len(self)):
            pointwise = self[i].regions.toPointwise_withLinearXYs(0.001, 1e-8, 1e-8)
            integral = pointwise.integrate()
            if abs(integral - 1.0) > info['normTolerance']:
                warnings.append( warning.unnormalizedDistribution( PhysicalQuantityWithUncertainty(self[i].energy, self.axes[0].unit), i,
                    integral, self ) )

            if pointwise.yMin() < 0.0:
                warnings.append( warning.negativeProbability( PhysicalQuantityWithUncertainty(self[i].energy, self.axes[0].unit),
                    value=pointwise.yMin(), obj=self[i] ) )

        return warnings

    def getAtEnergy( self, energy ) :

        accuracy = 1e-3
        epsilon = 1e-8
        energyRegion1 = None
        for energyRegion2 in self :
            if( energy <= energyRegion2.energy ) : break
            energyRegion1 = energyRegion2
        if( energyRegion1 is None ) : energy = energyRegion2.energy
        if( energy >= energyRegion2.energy ) : return( energyRegion2.regions.toPointwise_withLinearXYs( accuracy, epsilon, epsilon, removeOverAdjustedPoints = True ) )
        f = ( energy - energyRegion1.energy ) / ( energyRegion2.energy - energyRegion1.energy )
        p1 = energyRegion1.regions.toPointwise_withLinearXYs( accuracy, epsilon, epsilon, removeOverAdjustedPoints = True )
        p2 = energyRegion2.regions.toPointwise_withLinearXYs( accuracy, epsilon, epsilon, removeOverAdjustedPoints = True )
        try :
            distribution = ( 1. - f ) * p1 + f * p2
        except :                        # This is a kludge for some bad data. If it fails a raise will be executed.
            m1, m2 = p1.mutualify( 0, epsilon, 1, p2, epsilon, 0, 1 )
            distribution = ( 1. - f ) * m1 + f * m2
        return( distribution )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        return( [ energy.energy for energy in self ] )

    def EpAverageAtE( self, E ) :

        return( self.getAtEnergy( E ).integrateWithWeight_x( ) )

    def sqrtEp_AverageAtE( self, E ) :

        return( self.getAtEnergy( E ).integrateWithWeight_sqrt_x( ) )

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        independent, dependent, qualifier = self.axes[0].interpolation.getInterpolationTokens( )
        if( independent not in [ axes.linearToken ] ) : raise Exception( 'independent = %s not supported' % independent )
        if( dependent not in [ axes.linearToken ] ) : raise Exception( 'dependent = %s not supported' % dependent )
        axes_ = pointwise.defaultAxes( )
        pwl = linear( axes_, self.productFrame )
        axesXY = axes.referenceAxes( pwl, dimension = 2 )
        for region in self :
            r = region.toPointwise_withLinearXYs( accuracy = accuracy, lowerEps = lowerEps, upperEps = upperEps )
            r.value = region.energy
            pwl.append( r )
        return( pwl )

    def toXMLList( self, indent = '' ) :
        """Returns the xml string representation of self."""

        indent2 = indent + '  '
        xmlString = [ self.XMLStartTagString( indent = indent ) ]
        xmlString += self.axes.toXMLList( indent = indent2 )
        for energy in self.energies : xmlString += energy.toXMLList( indent2 )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    def toENDF6( self, flags, targetInfo, weight = None ) :

        from fudge.legacy.converting import endfFormats, gndToENDF6
        NE = len( self.energies )
        if( weight is None ) : weight = [ [ self.energies[0].getValue( ), 1.0 ], [ self.energies[-1].getValue( ), 1.0 ] ]
        data = weight
        EInInterpolation = gndToENDF6.axisToEndfInterpolationFlag( self.axes[0] )
        ENDFDataList = [ endfFormats.endfContLine( 0, 0, 0, 1, 1, len( data ) ) ] + \
            endfFormats.endfInterpolationList( [ len( data ), 2 ] ) + endfFormats.endfNdDataList( data ) + \
            [ endfFormats.endfContLine( 0, 0, 0, 0, 1, NE ) ] + endfFormats.endfInterpolationList( [ NE, EInInterpolation ] )
        for energy in self.energies : ENDFDataList += energy.toENDF6( flags, targetInfo )
        return( 1, ENDFDataList )

    @staticmethod
    def defaultAxes( ) :

        axes_ = axes.axes( dimension = 3 )
        axes_[0] = axes.axis( 'energy_in',  0, 'eV',   interpolation = axes.interpolationXY( axes.linearToken, axes.linearToken ) )
        axes_[1] = axes.axis( 'energy_out', 1, 'eV',   interpolation = axes.interpolationXY( axes.byRegionToken, axes.byRegionToken ) )
        axes_[2] = axes.axis( 'P(energy_out|energy_in)', 2, '1/eV' )
        return( axes_ )

    @staticmethod
    def parseXMLNode( element, xPath=[], linkData={} ):

        xPath.append( element.tag )
        axes_ = axes.parseXMLNode( element[0], xPath )
        spw = semiPiecewise( axes_, element.get( 'productFrame' ) )
        for energy_in in element[1:]:
            regions_ = regions.regionsXYs( axes.referenceAxes(spw), parent=spw, isPrimaryXData = False )
            for region_ in energy_in[0]:
                interp = axes.interpolationXY.stringToInterpolationXY( region_[0].get("interpolation") )
                axes_ = axes.interpolationAxes( int(region_[0].get("index")), interp,
                        regions_ )
                data = map(float, region_[1].text.split() )
                data = zip( data[::2], data[1::2] )
                regions_.regions.append( XYs.XYs( axes_, data, float( region_.get( "accuracy" ) ), parent = regions_,
                    index = int( region_.get( "index" ) ), moniker = "region", isPrimaryXData = False ) )
            spw.append( float( energy_in.get( "value" ) ), regions_ )
        xPath.pop()
        return spw

class semiPiecewiseEnergy( ancestry.ancestry ) :

    def __init__( self, index, energy, regions_, parent ) :

        if( not( isinstance( regions_, regions.regions ) ) ) : raise Exception( 'regions_ not an instance of regions.regions: is an instance of %s' % \
            brb.getType( regions_ ) )
        ancestry.ancestry.__init__( self, 'semiPiecewiseEnergy', parent )
        self.index = index
        self.energy = energy
        self.regions = regions_

    def __len__( self ) :

        return( len( self.regions ) )

    def __getitem__( self, index ) :

        return( self.regions[index] )

    def getValue( self ) :

        return( self.energy )

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        return( self.regions.toPointwise_withLinearXYs( accuracy = accuracy, lowerEps = lowerEps, upperEps = upperEps ) )

    def toXMLList( self, indent = '' ) :
        """Returns the xml string representation of self."""

        xmlString = [ '%s<energy_in value="%s" index="%d">' % ( indent, fudge.gnd.miscellaneous.floatToString( self.energy ), self.index ) ]
        xmlString += self.regions.toXMLList( indent + '  ' )
        xmlString[-1] += '</energy_in>'
        return xmlString

    def toENDF6( self, flags, targetInfo ) :

        from fudge.legacy.converting import endfFormats, gndToENDF6
        interpolations, ENDFData = [], []
        for index, region in enumerate( self.regions ) :
            data = region.copyDataToXYs( xUnit = 'eV', yUnit = '1/eV' )
            if( index > 0 ) :
                if( ENDFData[-1] == data[0] ) : data = data[1:]
            ENDFData += data
            interpolations += [ len( ENDFData ), gndToENDF6.axesToEndfInterpolationFlag( region.axes ) ]
        return( [ endfFormats.endfContLine( 0., self.energy, 0, 0, len( interpolations ) / 2, len( ENDFData ) ) ] + \
            endfFormats.endfInterpolationList( interpolations ) + endfFormats.endfNdDataList( ENDFData ) )

class energyFunctionalData( XYs.XYs ) :

    mutableYUnit = False

    def __init__( self, axes, data, accuracy, dataForm = 'xys', **kwargs ) :

        moniker = kwargs['moniker']
        if( moniker not in [ 'theta', 'a', 'b', 'g', 'T_M' ] ) : raise Exception( 'Invalid moniker = "%s"' % moniker )
        self.moniker = moniker
        self.tag = moniker
        XYs.XYs.__init__( self, axes, data, accuracy, dataForm = dataForm, isPrimaryXData = True )

class functionalBase( form ) :

    def __init__( self, productFrame, moniker, LF, U, parameter1, parameter2 = None ) :

        form.__init__( self, moniker, productFrame )
        self.U = U
        self.LF = LF
        self.parameter1 = parameter1
        self.parameter2 = parameter2

    def check( self, info ):

        from fudge.gnd import warning

        warnings = []
        if( ( self.domainMin( unitTo = 'eV' ) - self.U.getValueAs( 'eV' ) ) < 0 ) :
            warnings.append( warning.energyDistributionBadU( self ) )
        return( warnings )

    def domainMin( self, unitTo = None, asPQU = False ) :

        return( self.parameter1.domainMin( unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        return( self.parameter1.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getDomain( self, unitTo = None, asPQU = False ):

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ), self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        Es = [ E for E, p in self.parameter1 ]
        if( EMin is not None ) :
            if( EMin < Es[0] ) : Es.insert( 0, EMin )
        if( EMax is not None ) :
            if( EMax > Es[-1] ) : Es.append( EMax )
        return( Es )

    def toXMLList( self, indent = '' ) :
        """Returns the xml string representation of self."""

        if( self.LF == 12 ) : 
            qualifiers = 'EFL="%s" EFH="%s"' % ( self.EFL, self.EFH )
        else :
            qualifiers = 'U="%s"' % self.U
        xmlString = [ self.XMLStartTagString( indent = indent, extraAttributesAsStrings = qualifiers ) ]
        xmlString += self.parameter1.toXMLList( indent = indent + '  ' )
        if( not( self.parameter2 is None ) ) : xmlString += self.parameter2.toXMLList( indent = indent + '  ' )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    def toENDF6( self, flags, targetInfo, weight = None ) :

        from fudge.legacy.converting import endfFormats, gndToENDF6
        U, EFL, EFH = 0, 0, 0
        if( self.LF == 12 ) :
            EFL, EFH = self.EFL.getValueAs( 'eV' ), self.EFH.getValueAs( 'eV' )
        else :
            U = self.U.getValueAs( 'eV' )
        energyFactor = PhysicalQuantityWithUncertainty( '1 eV' ) / PhysicalQuantityWithUncertainty( 1, self.parameter1.axes[0].getUnit( ) )
        if( weight is None ) :
            weight = [ [ energyFactor * self.parameter1[0][0], 1.0 ], [ energyFactor * self.parameter1[-1][0], 1.0 ] ]
            interp = 2
        elif( hasattr(weight, 'axes') ):
            interp = {axes.linearToken: 2, axes.flatToken: 1}.get( weight.axes[0].interpolation.dependent )
            if interp is None:
                raise Exception("Only linear and flat interpolations supported for functional energy distributions!")
        else:
            interp = 2
        ENDFDataList = [ endfFormats.endfContLine( U, 0, 0, self.LF, 1, len( weight ) ) ] + \
            endfFormats.endfInterpolationList( [ len( weight ), interp ] ) + endfFormats.endfNdDataList( weight )
        ENDFDataList += [ endfFormats.endfContLine( EFL, EFH, 0, 0, 1, len( self.parameter1 ) ) ] + \
            endfFormats.endfInterpolationList( [ len( self.parameter1 ), \
            gndToENDF6.axesToEndfInterpolationFlag( self.parameter1.axes ) ] )
        ENDFDataList += endfFormats.endfNdDataList( self.parameter1, xUnit = 'eV', yUnit = 'eV' )
        if( self.parameter2 is not None ) :
            yUnit = ''
            if( self.LF in [ 5, 11 ] ) : yUnit = '1/eV'
            ENDFDataList += [ endfFormats.endfContLine( 0, 0, 0, 0, 1, len( self.parameter2 ) ) ]
            ENDFDataList += endfFormats.endfInterpolationList( [ len( self.parameter2 ), 
                gndToENDF6.axesToEndfInterpolationFlag( self.parameter2.axes ) ] )
            ENDFDataList += endfFormats.endfNdDataList( self.parameter2, yUnit = yUnit )
        return( 1, ENDFDataList )

class generalEvaporationSpectrum( functionalBase ) :

    def __init__( self, productFrame, U, thetas, gs ) :

        functionalBase.__init__( self, productFrame, base.generalEvaporationFormToken, 5, U, thetas, gs )

    def EpAverageAtE( self, E ) :

        return( self.parameter1.getValue( E ) * self.parameter2.integrateWithWeight_x( ) )

    def sqrtEp_AverageAtE( self, E ) :

        return( math.sqrt( self.parameter1.getValue( E ) ) * self.parameter2 ).integrateWithWeight_sqrt_x( )

    def isLinear( self, qualifierOk = False, flatIsOk = False ) :
        """
        Returns the results of isLinear called on the axes of g(E'|E).
        """

        return( self.parameter2.axes.isLinear( qualifierOk = qualifierOk, flatIsOk = flatIsOk ) )

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        axes_ = pointwise.defaultAxes( )
        pwl = linear( axes_ )
        axesXY = axes.referenceAxes( pwl, dimension = 2 )
        thetas = self.parameter1.toPointwise_withLinearXYs( accuracy = None, lowerEps = lowerEps, upperEps = upperEps )
        gs = self.parameter2.toPointwise_withLinearXYs( accuracy = None, lowerEps = lowerEps, upperEps = upperEps )
        for E_in, theta in thetas :
            data = [ [ theta * x, y / theta ] for x, y in gs ]
            data = XYs.XYs( axesXY, data, accuracy, value = E_in )
            data.normalize( insitu = True )
            pwl.append( data )
        return( pwl )

    @staticmethod
    def parseXMLNode( element, xPath=[], linkData={} ):
        """Translate <generalEvaporation> element from xml."""

        xPath.append( element.tag )
        theta, g = [energyFunctionalData.parseXMLNode(node,xPath) for node in element]
        U = PhysicalQuantityWithUncertainty( element.get("U") )
        GES = generalEvaporationSpectrum( element.get( 'productFrame' ), U, theta, g )
        xPath.pop()
        return GES

class simpleMaxwellianFissionSpectrum( functionalBase ) :

    def __init__( self, productFrame, U, thetas ) :

        functionalBase.__init__( self, productFrame, base.simpleMaxwellianFissionFormToken, 7, U, thetas )

    def EpAverageAtE( self, E ) :

        theta = self.parameter1.getValue( E )
        a = ( E - self.U.getValue( ) ) / theta
        if( a < 1e-4 ) : return( theta * a * ( 1575. - a * ( 180. + 8 * a ) ) / 2625. )
        sqrt_a = math.sqrt( a )
        exp_a = math.exp( -a )
        erf_sqrt_a = math.sqrt( math.pi ) * math.erf( sqrt_a )
        return( theta * ( 0.75 * erf_sqrt_a - sqrt_a * ( 1.5 + a ) * exp_a ) / ( 0.5 * erf_sqrt_a - sqrt_a * exp_a ) )

    def sqrtEp_AverageAtE( self, E ) :

        theta = self.parameter1.getValue( E )
        a = ( E - self.U.getValue( ) ) / theta
        if( a < 1e-4 ) : return( math.sqrt( theta * a ) * ( 2100. - a * ( 140. + 9. * a ) ) / 2800. )
        sqrt_a = math.sqrt( a )
        exp_a = math.exp( -a )
        erf_sqrt_a = math.sqrt( math.pi ) * math.erf( sqrt_a )
        return( math.sqrt( theta ) * ( 1 - ( 1. + a ) * exp_a ) / ( 0.5 * erf_sqrt_a - sqrt_a * exp_a ) )

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        def evaluateAtX( self, x ) :

            return( math.sqrt( x ) * math.exp( -x / self.p1 ) )

        ef = energyFunctionalDataToPointwise( self, evaluateAtX )
        return( ef.toPointwise_withLinearXYs( accuracy, lowerEps = lowerEps, upperEps = upperEps ) )

    @staticmethod
    def parseXMLNode( MFelement, xPath=[], linkData={} ):

        xPath.append( MFelement.tag )
        theta = energyFunctionalData.parseXMLNode(MFelement[0],xPath)
        U = PhysicalQuantityWithUncertainty( MFelement.get("U") )
        SMF = simpleMaxwellianFissionSpectrum( MFelement.get( 'productFrame' ), U, theta )
        xPath.pop()
        return SMF

class evaporationSpectrum( functionalBase ) :

    def __init__( self, productFrame, U, thetas ) :

        functionalBase.__init__( self, productFrame, base.evaporationFormToken, 9, U, thetas )

    def EpAverageAtE( self, E ) :

        theta = self.parameter1.getValue( E )
        a = ( E - self.U.getValue( ) ) / theta
        if( a < 1e-4 ) : return( theta * a * ( 180. - a * ( 15. + a ) ) / 270. )
        exp_a = math.exp( -a )
        return( theta * ( 2. - a**2 * exp_a / ( 1. - ( 1. + a ) * exp_a ) ) )

    def sqrtEp_AverageAtE( self, E ) :

        theta = self.parameter1.getValue( E )
        a = ( E - self.U.getValue( ) ) / theta
        if( a < 1e-4 ) : return( math.sqrt( theta * a ) * ( 252. - a * ( 12. + a ) ) / 315. )
        sqrt_a = math.sqrt( a )
        exp_a = math.exp( -a )
        return( math.sqrt( theta ) * ( 1.32934038817913702 * math.erf( sqrt_a ) - sqrt_a * ( 1.5 + a ) * exp_a ) / ( 1. - ( 1. + a ) * exp_a ) )

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        def evaluateAtX( self, x ) :

            return( x * math.exp( -x / self.p1 ) )

        ef = energyFunctionalDataToPointwise( self, evaluateAtX )
        return( ef.toPointwise_withLinearXYs( accuracy, lowerEps = lowerEps, upperEps = upperEps ) )

    @staticmethod
    def parseXMLNode( evapElement, xPath=[], linkData={} ):

        xPath.append( evapElement.tag )
        theta = energyFunctionalData.parseXMLNode( evapElement[0], xPath )
        U = PhysicalQuantityWithUncertainty( evapElement.get("U") )
        ES = evaporationSpectrum( evapElement.get( 'productFrame' ), U, theta )
        xPath.pop()
        return ES

class WattSpectrum( functionalBase ) :

    def __init__( self, productFrame, U, a, b ) :

        functionalBase.__init__( self, productFrame, base.WattFormToken, 11, U, a, b )

    def EpAverageAtE( self, E ) :

        a, b = self.parameter1.getValue( E ), self.parameter2.getValue( E )
        xMax_a  = ( E - self.U.getValue( ) ) / a
        xMax_b  = math.sqrt( b * ( E - self.U.getValue( ) ) )
        ab = a * b
        sqrt_ab = math.sqrt( ab )
        I = 0.25 * math.sqrt( math.pi * ab ) * math.exp( 0.25 * ab ) * \
            ( math.erf( math.sqrt( xMax_a ) - 0.5 * sqrt_ab ) + math.erf( math.sqrt( xMax_a ) + 0.5 * sqrt_ab ) ) \
            - math.exp( -xMax_a ) * math.sinh( xMax_b )
        EI = a * math.sqrt( math.pi * ab ) * ( ab + 6 ) * math.exp( 0.25 * ab ) * \
            ( math.erf( math.sqrt( xMax_a ) - 0.5 * sqrt_ab ) + math.erf( math.sqrt( xMax_a ) + 0.5 * sqrt_ab ) ) \
            - 0.25 * a * math.exp( -xMax_a ) * math.sinh( xMax_b ) * ( 2. * xMax_b + ab + 4. + 4. * xMax_a )
        return( EI / ( 16 * I ) )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        aMin, aMax = self.parameter1.getDomain( )
        if( EMin is None ) : EMin = aMin
        if( EMax is None ) : EMax = aMax
        if( EMin < aMin ) : EMin = aMin
        if( EMax < aMax ) : EMax = aMax
        Es = [ EMin, EMax ]
        for E, a in self.parameter1 :
            if( E <= EMin ) : continue
            if( E >= EMax ) : continue
            if( E not in Es ) : Es.append( E )
        for E, b in self.parameter2 :
            if( E <= EMin ) : continue
            if( E >= EMax ) : continue
            if( E not in Es ) : Es.append( E )
        Es.sort( )
        return( Es )

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        def evaluateAtX( self, x ) :

            return( math.exp( -x / self.p1 ) * math.sinh( math.sqrt( self.p2 * x ) ) )

        ef = energyFunctionalDataToPointwise( self, evaluateAtX )
        return( ef.toPointwise_withLinearXYs( accuracy, lowerEps = lowerEps, upperEps = upperEps ) )

    @staticmethod
    def parseXMLNode( WattElement, xPath=[], linkData={} ):
        """ translate <Watt> element from xml """

        xPath.append( WattElement.tag )
        a, b = [energyFunctionalData.parseXMLNode(node,xPath) for node in WattElement]
        U = PhysicalQuantityWithUncertainty( WattElement.get("U") )
        WS = WattSpectrum( WattElement.get( 'productFrame' ), U, a, b )
        xPath.pop()
        return WS

class MadlandNix( functionalBase ) :

    def __init__( self, productFrame, EFL, EFH, Ts ) :

        functionalBase.__init__( self, productFrame, base.MadlandNixFormToken, 12, None, Ts )
        self.EFL = EFL
        self.EFH = EFH

    def check( self, info ) :

        from fudge.gnd import warning

        warnings = []
        if self.EFL.value <= 0 or self.EFH.value <= 0 or self.parameter1.yMin() <= 0:
            warnings.append( warning.MadlandNixBadParameters( self.EFL, self.EFH, self.parameter1.yMin(), self ) )

        return warnings

    def EpAverageAtE( self, E ) :

        unit = self.parameter1.axes[1].getUnit( )
        return( 0.5 * ( self.EFL.getValueAs( unit ) + self.EFH.getValueAs( unit ) ) + 4. * self.parameter1.getValue( E ) / 3. )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        return( [ x for x, y in self.parameter1 ] )

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        from numericalFunctions import specialFunctions

        def MadlandNixFunc( Ep, parameters ) :

            def g( Ep, E_F, T_M ) :

                u1 = ( math.sqrt( Ep ) - math.sqrt( E_F ) )
                u1 *= u1 / T_M
                u2 = ( math.sqrt( Ep ) + math.sqrt( E_F ) )
                u2 *= u2 / T_M
                E1 = 0                      # u1^3/2 * E1 is zero for u1 = 0. but E1 is infinity, whence, the next test.
                if( u1 != 0 ) : E1 = specialFunctions.exponentialIntegral( 1, u1 )
                E2 = specialFunctions.exponentialIntegral( 1, u2 )
                complementary = ( u1 > 2. )
                gamma1 = specialFunctions.incompleteGamma( 1.5, u1, complementary )
                gamma2 = specialFunctions.incompleteGamma( 1.5, u2, complementary )
                signG = 1
                if( complementary ) : signG = -1
                return( ( u2 * math.sqrt( u2 ) * E2 - u1 * math.sqrt( u1 ) * E1 + signG * ( gamma2 - gamma1 ) ) / ( 3 * math.sqrt( E_F * T_M ) ) )

            EFL, EFH, T_M = parameters
            return( 0.5 * ( g( Ep, EFL, T_M ) + g( Ep, EFH, T_M ) ) )

        axes_ = pointwise.defaultAxes( )
        pwl = linear( axes_ )
        axesXY = axes.referenceAxes( pwl, dimension = 2 )
        E_in_unit = self.parameter1.axes[0].getUnit( )
        EFL, EFH = self.EFL.getValueAs( E_in_unit ), self.EFH.getValueAs( E_in_unit )
        factor = PhysicalQuantityWithUncertainty( 1, 'eV' ).getValueAs( E_in_unit )
        xs_ = [ 1e-5, 1e-3, 1e-1, 1e1, 1e3, 1e5, 3e7 ]
        xs = [ factor * x for x in xs_ ]
        for E, T_M in self.parameter1 :             # This logic ignores the interpolation of parameter1 as the only two forms in ENDF/B-VII shows 
            parameters = [ EFL, EFH, T_M ]          # that linear-linear is better than the 'log-log' given in the ENDF/B-VII/data.
            g_Ep = XYs.XYs.createFromFunction( axesXY, xs, MadlandNixFunc, parameters, accuracy, biSectionMax = 12 )
            g_Ep.value = E                          # ????????? Class XYs does not the a proper setValue method. One should be added.
            g_Ep.normalize( insitu = True )
            pwl.append( g_Ep )
        return( pwl )

    @staticmethod
    def parseXMLNode( MNelement, xPath=[], linkData={} ):
        """ translate <MadlandNix> element from xml """

        xPath.append( MNelement.tag )
        T_M = energyFunctionalData.parseXMLNode( MNelement[0], xPath )
        EFL, EFH = [PhysicalQuantityWithUncertainty(tmp) for tmp in (MNelement.get("EFL"),MNelement.get("EFH") )]
        MN = MadlandNix( MNelement.get( 'productFrame' ), EFL, EFH, T_M )
        xPath.pop()
        return MN

class NBodyPhaseSpace( form ) :

    def __init__( self, numberOfProducts, numberOfProductsMasses ) :

        form.__init__( self, base.NBodyPhaseSpaceFormToken, axes.centerOfMassToken )
        self.numberOfProducts = numberOfProducts
        self.numberOfProductsMasses = numberOfProductsMasses

    def check( self, info ) :

        #if ( abs( self.numberOfProductsMasses - info['reactionSuite'].particles['n'].getMass('amu') * self.numberOfProducts ) >
        #        self.numberOfProductsMasses * 0.1 ) :    # return warning?

        return []

    def EpAverageAtE( self, E, massUnit, projectileMass, targetMass, productMass, Q ) :
        """
        Calculate the average energy of the product in the center-of-mass frame for projectile energy E. This method only works for
        a one-step reaction.
        """

        M = self.numberOfProductsMasses.getValueAs( massUnit )
        Ea = targetMass / ( targetMass + projectileMass ) * E + Q
        return( Ea * ( M - productMass ) / ( M * ( self.numberOfProducts - 1 ) ) )

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        from fudge.gnd import product
        from fudge.gnd import reactionSuite
        from fudge.gnd.reactions import reaction
        from fudge.gnd import channels
        from fudge.core.math import fudgemath

        class tester :

            def __init__( self, relativeTolerance, absoluteTolerance, n ) :

                self.relativeTolerance = relativeTolerance
                self.absoluteTolerance = absoluteTolerance
                self.n = n
                self.setEMax_i( 1 )

            def evaluateAtX( self, x ) :

                return( math.sqrt( x ) * math.pow( ( self.EMax_i - x ), 0.5 * ( 3. * self.n - 5. ) ) )

            def setEMax_i( self, EMax_i ) :

                self.EMax_i = EMax_i

        p = self.findClassInAncestry( product.product )
        numberOfProductsMasses, massUnit = self.numberOfProductsMasses.getValue( ), self.numberOfProductsMasses.getUnitName( )
        productMass = p.getMass( massUnit )

        r = self.findClassInAncestry( reaction.reaction )
        energyUnit = r.getDomainUnit( )
        EMin, EMax = r.getDomain( )

        rs = self.findClassInAncestry( reactionSuite.reactionSuite )
        projectileMass, targetMass = rs.projectile.getMass( massUnit ), rs.target.getMass( massUnit )

        c = self.findClassInAncestry( channels.channel )
        Q = c.Q.getConstantAs( energyUnit )

        axes_ = pointwise.defaultAxes( axes.centerOfMassToken, energyUnit = energyUnit )
        pwl = linear( axes_ )
        axesXY = axes.referenceAxes( pwl, dimension = 2 )

        t = tester( accuracy, 1e-10, self.numberOfProducts )
        n = 21
        f = math.pow( EMax / EMin, 1. / n )
        E_ins = [ EMin * f**i for i in xrange( n ) ]
        E_ins[-1] = EMax        # Fix possible round off issue.
        for i, E_in in enumerate( E_ins ) :
            Ea = targetMass / ( targetMass + projectileMass ) * E_in + Q
            EMax_i = Ea * ( numberOfProductsMasses - productMass ) / numberOfProductsMasses
            if( EMax_i < 0 ) : EMax_i = 1e-5                # This is a klude
            t.setEMax_i( EMax_i )
            t.absoluteTolerance = 1e-10 * t.evaluateAtX( 0.5 * EMax_i )
            data = fudgemath.thickenXYList( [ [ 0., 0. ], [ EMax_i, 0. ] ], t, biSectionMax = 10 )
            data = XYs.XYs( axesXY, data, accuracy, value = E_in )
            data.normalize( insitu = True )
            pwl.append( data )
        return( pwl )

    def toXMLList( self, indent = '' ) :
        """Returns the xml string representation of self."""

        xmlString = [ self.XMLStartTagString( indent = indent, extraAttributesAsStrings = 'numberOfProducts="%s" mass="%s"' %
            ( self.numberOfProducts, self.numberOfProductsMasses ), emptyTag = True ) ]
        return( xmlString )

    def toENDF6( self, flags, tempInfo ) :

        from fudge.legacy.converting import endfFormats, gndToENDF6
        mass = self.numberOfProductsMasses.getValueAs( 'eV/c**2' ) / tempInfo['neutronMass']
        ENDFDataList = [ endfFormats.endfContLine( mass, 0, 0, 0, 0, self.numberOfProducts ) ]
        return( 6, self.getProductFrame( ), ENDFDataList )

    @staticmethod
    def parseXMLNode( NBodyElement, xPath=[], linkData={} ):
        return NBodyPhaseSpace( int( NBodyElement.get( "numberOfProducts" ) ), PhysicalQuantityWithUncertainty(NBodyElement.get("mass")) )

class weighted( ancestry.ancestry ) :

    moniker = 'weighted'

    def __init__( self, weight, functional ) :

        self.weight = weight
        self.functional = functional
        self.W_index = 0    # applies to the weighted section, not to the weights

    def checkProductFrame( self ) :
        "Calls checkProductFrame on self's functional."

        self.functional.checkProductFrame( )

    def EpAverageAtE( self, E ) :

        return( self.weight.getValue( E ) * self.functional.EpAverageAtE( E ) )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        energyArray = self.functional.getEnergyArray( EMin, EMax )
        for x, y in self.weight :
            if( x not in energyArray ) : energyArray.append( x )
        energyArray.sort( )
        return( energyArray )

    def getPrductFrame( self ) :
        "Returns the results of calling getPrductFrame on self's functional."

        return( self.functional.getProductFrame( ) )

    def toXMLList( self, indent = '' ) :
        """Returns the xml string representation of self."""

        indent2 = indent + '  '
        xmlString = [ '%s<%s index="%i">' %(  indent, self.moniker, self.W_index ) ]
        xmlString += self.weight.toXMLList( tag = 'weight', indent = indent2 )
        xmlString += self.functional.toXMLList( indent = indent2 )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    def toENDF6( self, flags, targetInfo ) :

        return( self.functional.toENDF6( flags, targetInfo, weight = self.weight )[1] )

class weightedFunctionals( form ) :

    def __init__( self, productFrame ) :

        form.__init__( self, base.weightedFunctionalsFormToken, productFrame )
        self.weights = []

    def __len__( self ) :

        return( len( self.weights ) )

    def __getitem__( self, i ) :

        return( self.weights[i] )

    def addWeight( self, weight ) :

        if( not( isinstance( weight, weighted ) ) ) : raise Exception( 'Invalid weight of type %s' % brb.getType( weight ) )
        weight.W_index = len(self.weights)
        self.weights.append( weight )

    def check( self, info ) :

        from fudge.gnd import warning

        warnings = []
        totalWeight = sum( [w.weight for w in self.weights] )
        if (totalWeight.yMin() != 1.0) and (totalWeight.yMax() != 1.0):
            warnings.append( warning.weightsDontSumToOne( obj=self ) )

        for weight in self:
            warnings += weight.functional.check( info )

        return warnings

    def checkProductFrame( self ) :
        """
        Checks that the productFrames for the each weight is valid and that all frames are the same. Overrides the
        method in base.form. Returns None if all is okay. Otherwise, executes a raises from base.form.checkProductFrame
        or a ValueError if frames differ.
        """

        frames = set( )
        for weight in self :
            weight.checkProductFrame( )
            frames.add( weight.getPrductFrame( ) )
        if( len( frames ) != 1 ) : raise ValueError( 'multiple frames present :%s' % frames )

    def EpAverageAtE( self, E ) :

        Ep = 0
        for weight in self : Ep += weight.EpAverageAtE( E )
        return( Ep )

    def sqrtEp_AverageAtE( self, E ) :
        """This method has not been implemented. It returns None so the method uncorrelated.calculateDepositionData will still work 
        when calculating the momentum deposition for isotropic scattering in the lab frame, but will execute a raise otherwise."""

        return( None )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        energyArray = []
        for weight in self :
            for energy in weight.getEnergyArray( EMin, EMax ) :
                if( energy not in energyArray ) : energyArray.append( energy )
        energyArray.sort( )
        return( energyArray )

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        if( len( self ) > 2 ) : raise Exception( 'more than two weights currently not supported' )
        E_ins, data = [], []
        for weight in self :
            w = weight.toPointwise_withLinearXYs( accuracy = accuracy, lowerEps = lowerEps, upperEps = upperEps )
            e = weight.functional.toPointwise_withLinearXYs( accuracy = accuracy, lowerEps = lowerEps, upperEps = upperEps )
            data.append( [ w, e ] )
            for x, y in w :
                if( x not in E_ins ) : E_ins.append( x )
            for x in e :
                if( x.value not in E_ins ) : E_ins.append( x.value )
        E_ins.sort( )
        axes_ = pointwise.defaultAxes( )
        pwl = linear( axes_ )
        axesXY = axes.referenceAxes( pwl, dimension = 2 )
        for E_in in E_ins :
            wv1, ev1 = data[0]
            wv2, ev2 = data[1]
            w1, w2 = wv1.getValue( E_in ), wv2.getValue( E_in )
            if( w1 == 1 ) :
                e = ev1.interpolateAtW( E_in, unitBase = True )
            elif( w2 == 1 ) :
                e = ev2.interpolateAtW( E_in, unitBase = True )
            else :
                raise Exception( 'One weight must be zero: %s, %s' % ( w1, w2 ) )
            e.value = E_in
            e.normalize( insitu = True )
            pwl.append( e )
        return( pwl )

    def toXMLList( self, indent = '' ) :
        """Returns the xml string representation of self."""

        xmlString = [ self.XMLStartTagString( indent = indent ) ]
        for form in self.weights : xmlString += form.toXMLList( indent + '  ' )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    def toENDF6( self, flags, targetInfo, weight = None ) :                 # The weight is ignored in this method, only for compatibility with other toENDF6's.

        ENDFDataList = []
        for weight in self.weights : ENDFDataList += weight.toENDF6( flags, targetInfo )
        return( len( self.weights ), ENDFDataList )

    @staticmethod
    def parseXMLNode( WFelement, xPath=[], linkData={} ):

        xPath.append( WFelement.tag )
        WF = weightedFunctionals( WFelement.get( 'productFrame' ) )
        for weightSection in WFelement:
            weights, functional = weightSection
            weight_ = XYs.XYs.parseXMLNode( weights, xPath )
            # functional data:
            formClass = {
                    base.generalEvaporationFormToken: generalEvaporationSpectrum,
                    base.WattFormToken: WattSpectrum,
                    base.MadlandNixFormToken: MadlandNix,
                    base.evaporationFormToken: evaporationSpectrum,
                    }.get(functional.tag)
            functional = formClass.parseXMLNode( functional, xPath )
            WF.weights.append( weighted( weight_, functional ) )
            WF.weights[-1].W_index = len(WF.weights)-1
        xPath.pop()
        return WF

class equalProbableBins( base.equalProbableBinsFormBase ) :

    def __init__( self, productFrame ) :

        base.equalProbableBinsFormBase.__init__( self, productFrame )

class energyFunctionalDataToPointwise :

    def __init__( self, data, evaluateAtX ) :

        self.data = data
        self.evaluateAtX = evaluateAtX

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        from fudge.core.math import fudgemath

        class tester :

            def __init__( self, relativeTolerance, absoluteTolerance, evaluateAtX ) :

                self.relativeTolerance = relativeTolerance
                self.absoluteTolerance = absoluteTolerance
                self.evaluateAtX_ = evaluateAtX
                self.setParameter1( 1. )
                self.setParameter2( 1. )

            def evaluateAtX( self, x ) :

                return( self.evaluateAtX_( self, x ) )

            def setParameter1( self, p1 ) :

                self.p1 = p1

            def setParameter2( self, p2 ) :

                self.p2 = p2

        parameter1 = self.data.parameter1.toPointwise_withLinearXYs( accuracy = None, lowerEps = lowerEps, upperEps = upperEps )
        axes_ = pointwise.defaultAxes( )
        pwl = linear( axes_ )
        axesXY = axes.referenceAxes( pwl, dimension = 2 )
        one_eV = PhysicalQuantityWithUncertainty( '1 eV' ).getValueAs( axesXY[0].getUnit( ) )

        t = tester( accuracy, 1e-10, self.evaluateAtX )
        parameter2 = self.data.parameter2
        if( parameter2 is not None ) : parameter2 = parameter2.toPointwise_withLinearXYs( accuracy = None, lowerEps = lowerEps, upperEps = upperEps )
        for i, E_in_p1 in enumerate( parameter1 ) :
            E_in, p1 = E_in_p1
            EpMax = E_in - self.data.U.getValue( )
            EpMin = 0.              # Only used with debugging.
            if( EpMax == 0. ) :
                if( i != 0 ) : raise Exception( "i = %d, E_in = %s, U = %s" % ( E_in, self.data.U.getValue( ) ) )
                EpMax = E_in * 1e-6
                if( EpMax == 0. ) :
                    EpMax = parameter1[1][1] * 1e-6         # This and the next line are arbitary.
                if( EpMax > 1e-3 ) : EpMax = 2e-5
                data = [ [ EpMin, 1. ], [ EpMax, 1e-10 ] ]
            else :
                t.setParameter1( p1 )
                if( parameter2 is not None ) : t.setParameter2( parameter2.getValue( E_in ) )
                level, data = 0, [ [ EpMin, t.evaluateAtX( EpMin ) ], [ EpMax, t.evaluateAtX( EpMax ) ] ]
                while( data[1][1] < 1e-10 ) :               # Fill in some tail points while they are < 1e-10
                    level, EpMax = level + 1, 0.5 * EpMax
                    data.insert( 1, [ EpMax, t.evaluateAtX( EpMax ) ] )
                    if( level > 10 ) : break
                if( data[0][0] < one_eV ) :
                    if( data[1][0] > 1e6 * one_eV ) : data.insert( 1,  [ 1e6 * one_eV, t.evaluateAtX( 1e6 * one_eV ) ] )
                    if( data[1][0] > 1e3 * one_eV ) : data.insert( 1,  [ 1e3 * one_eV, t.evaluateAtX( 1e3 * one_eV ) ] )
                    if( data[1][0] > one_eV ) : data.insert( 1,  [ one_eV, t.evaluateAtX( one_eV ) ] )
                data = fudgemath.thickenXYList( data, t, biSectionMax = 10 )
            data = XYs.XYs( axesXY, data, accuracy, value = E_in )
            data.normalize( insitu = True )
            pwl.append( data )
        return( pwl )

def parseXMLNode( energyElement, xPath=[], linkData={} ):
    """ translate <energy> element from xml """

    xPath.append( energyElement.tag )
    energy = component( nativeData = energyElement.get('nativeData') )
    for formElement in energyElement:
        formClass = {base.pointwiseFormToken: pointwise,
                base.linearFormToken : linear,
                base.semiPiecewiseFormToken: semiPiecewise,
                base.generalEvaporationFormToken: generalEvaporationSpectrum,
                base.WattFormToken: WattSpectrum,
                base.MadlandNixFormToken: MadlandNix,
                base.simpleMaxwellianFissionFormToken: simpleMaxwellianFissionSpectrum,
                base.evaporationFormToken: evaporationSpectrum,
                base.weightedFunctionalsFormToken: weightedFunctionals,
                base.NBodyPhaseSpaceFormToken: NBodyPhaseSpace,
                }.get(formElement.tag)
        if formClass is None: raise Exception("encountered unknown energy form: %s" % formElement.tag)
        newForm = formClass.parseXMLNode( formElement, xPath, linkData )
        energy.addForm( newForm )
    xPath.pop()
    return energy
