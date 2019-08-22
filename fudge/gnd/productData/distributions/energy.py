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

"""Energy distribution classes."""

import math
import miscellaneous
import fudge
from pqu import PQU as PQUModule
from fudge.core.utilities import brb
import xData.ancestry as ancestryModule

from xData import standards as standardsModule
from xData import base as xDataBaseModule
from xData import axes as axesModule
from xData import XYs as XYsModule
from xData import multiD_XYs as multiD_XYsModule
from xData import regions as regionsModule

from fudge.gnd import physicalQuantity as physicalQuantityModule

from . import base as baseModule
from . import xs_pdf_cdf as xs_pdf_cdfModule

__metaclass__ = type

def defaultAxes( energyUnit ) :

    axes = axesModule.axes( rank = 3 )
    axes[2] = axesModule.axis( 'energy_in',  2, energyUnit )
    axes[1] = axesModule.axis( 'energy_out', 1, energyUnit )
    axes[0] = axesModule.axis( 'P(energy_out|energy_in)', 0, '1/' + energyUnit )
    return( axes )

class XYs1d( XYsModule.XYs1d ) :

    def averageEnergy( self ) :

        allowedInterpolations = [ standardsModule.interpolation.linlinToken,
                                      standardsModule.interpolation.flatToken ]
        xys = self.changeInterpolationIfNeeded( allowedInterpolations, XYsModule.defaultAccuracy )
        return( xys.integrateWithWeight_x( ) )

class regions1d( regionsModule.regions1d ) :

    def averageEnergy( self ) :

        averageEnergy = 0
        for region in self : averageEnergy += region.averageEnergy( )
        return( averageEnergy )

    def integrateWithWeight_x( self ) :

        sum = 0
        for region in self : sum += region.integrateWithWeight_x( )
        return( sum )

    def toLinearXYsClass( self ) :

        return( XYs1d )

    @staticmethod
    def allowedSubElements():

        return( XYs1d, )

class xs_pdf_cdf1d( xs_pdf_cdfModule.xs_pdf_cdf1d ) :

    pass

class subform( baseModule.subform ) :
    """Abstract base class for energy forms."""

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :

        return( None )

class discretePrimaryGamma( subform ) :

    dimension = 2
    ancestryMembers = ( '', )

    def __init__( self, value, domainMin, domainMax, axes = None ) :

        subform.__init__( self )

        if( isinstance( value, int ) ) : value = float( value )
        if( not( isinstance( value, float ) ) ) : raise TypeError( 'value must be a float.' )
        self.value = value

        if( isinstance( domainMin, int ) ) : domainMin = float( domainMin )
        if( not( isinstance( domainMin, float ) ) ) : raise TypeError( 'domainMin must be a float.' )
        self.domainMin = domainMin

        if( isinstance( domainMax, int ) ) : domainMax = float( domainMax )
        if( not( isinstance( domainMax, float ) ) ) : raise TypeError( 'domainMax must be a float.' )
        self.domainMax = domainMax

        if( axes is None ) :
            self.__axes = None
        else :
            if( not( isinstance( axes, axesModule.axes ) ) ) : raise TypeError( 'axes is not an axes instance' )
            if( len( axes ) <= self.dimension ) : raise Exception( 'len( axes ) = %d != ( self.dimension + 1 ) = %d' % ( len( axes ), ( self.dimension + 1 ) ) )
            self.__axes = axes.copy( )
            self.__axes.setAncestor( self )

    @property
    def axes( self ) :

        return( self.__axes )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        factors = self.axes.convertUnits( unitMap )
        self.value *= factors[1]
        self.domainMin *= factors[2]
        self.domainMax *= factors[2]

    def copy( self ):

        return self.__class__( self.value, self.domainMin, self.domainMax, self.axes )

    __copy__ = __deepcopy__ = copy

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s value="%s" domainMin="%s" domainMax="%s"' % ( indent, self.moniker, self.value, self.domainMin, self.domainMax ) ]
        if( self.axes is None ) : 
            XMLStringList[-1] += '/>'
        else :
            XMLStringList[-1] += '>'
            XMLStringList += self.axes.toXMLList( indent = indent2, **kwargs )
            XMLStringList[-1] += '</%s>' % self.moniker
        return( XMLStringList )

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ) :

        value = float( element.get( 'value' ) )
        domainMin = float( element.get( 'domainMin' ) )
        domainMax = float( element.get( 'domainMax' ) )

        axes = None
        for child in element :
            if( child.tag == axesModule.axes.moniker ) :
                axes = axesModule.axes.parseXMLNode( child, xPath, linkData )
            else :
                raise Exception( 'Invalid sub-element with tag = "%s"' % child.tag )
        
        return( cls( value, domainMin, domainMax, axes = axes ) )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        return( [ EMin, EMax ] )

class discreteGamma( discretePrimaryGamma ) :

    moniker = 'discreteGamma'

    def check( self, info ) :

        from fudge.gnd import warning

        warnings = []
        if( self.value <= 0 ) : warnings.append( warning.negativeDiscreteGammaEnergy() )
        return( warnings )

    def averageEp( self, E ) :

        return( self.value )

class primaryGamma( discretePrimaryGamma ) :

    moniker = 'primaryGamma'

    def __init__( self, value, domainMin, domainMax, axes = None ) :

        discretePrimaryGamma.__init__( self, value, domainMin, domainMax, axes = axes )
        self.__massRatio = None     # In ENDF lingo this is AWR / ( AWR + 1 ).

    @property
    def massRatio( self ) :

        if( self.__massRatio is None ) :
            self.__massRatio = self.findAttributeInAncestry( "getMassRatio" )( )
        return self.__massRatio

    def check( self, info ) :

        from fudge.gnd import warning

        warnings = []
# BRB6 hardwired
        Qvalue = self.findAttributeInAncestry('getQ')('eV')
        if isinstance( self.value, PQUModule.PQU ) :
            testValue = self.value.getValueAs( 'eV' )
        else:
            testValue = self.value
        if testValue > Qvalue:
            warnings.append( warning.primaryGammaEnergyTooLarge( self.value,
                            100 * testValue / Qvalue ) )
        return warnings

    def averageEp( self, E ) :

        return( float( self.value ) + self.massRatio * E )

class XYs2d( subform, multiD_XYsModule.XYs2d ) :

    def __init__( self, **kwargs ):
        """
        >pointwise = XYs2d( )
        followed by:
        >pointwise[ 0 ] = XYs_data_1
        >pointwise[ 1 ] = XYs_data_2
        > ...
        >pointwise[ n-1 ] = XYs_data_n
        """

        multiD_XYsModule.XYs2d.__init__( self, **kwargs )
        subform.__init__( self )

    def getAtEnergy( self, energy ) :
        "This method is deprecated, use getSpectrumAtEnergy."

        return( self.getSpectrumAtEnergy( energy ) )

    def getSpectrumAtEnergy( self, energy ) :
        """Returns the energy spectrum for self at projectile energy."""

#        unitBase = self.axes[0].interpolation.isQualifierUnitBase( ) or self.axes[0].interpolation.isQualifierCorrespondingPoints( )
        unitBase = True
        return( self.interpolateAtValue( energy, unitBase = unitBase, extrapolation = standardsModule.flatExtrapolationToken ) )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        Es = [ data.value for data in self ]
        if( EMin is not None ) :
            if( EMin < Es[0] ) : Es.insert( 0, EMin )
        if( EMax is not None ) :
            if( EMax > Es[-1] ) : Es.append( EMax )
        return( Es )

    def averageEp( self, energy ) :

        code, lower, upper, frac = self.getBoundingSubFunctions( energy )
        if( code is None ) :
            raise Exception( 'No distribution' )
        elif( code == '' ) :
            EpLowerMin, EpLowerMax = lower.domainMin, lower.domainMax
            EpUpperMin, EpUpperMax = upper.domainMin, upper.domainMax
            ELower = lower.value            # BRB FIXME after xData getValue is fixed.
            EUpper = upper.value
            EFraction = ( energy - ELower ) / ( EUpper - ELower )
            EpMidMin = ( 1 - EFraction ) * EpLowerMin + EFraction * EpUpperMin
            EpMidMax = ( 1 - EFraction ) * EpLowerMax + EFraction * EpUpperMax
            EpLower = ( lower.integrateWithWeight_x( ) - EpLowerMin ) / ( EpLowerMax - EpLowerMin ) * \
                    ( EpMidMax - EpMidMin ) + EpMidMin
            EpUpper = ( upper.integrateWithWeight_x( ) - EpUpperMin ) / ( EpUpperMax - EpUpperMin ) * \
                    ( EpMidMax - EpMidMin ) + EpMidMin
            return( ( 1 - EFraction ) * EpLower + EFraction * EpUpper )
        else :
            return( lower.integrateWithWeight_x( ) )

    def check( self, info ) :

        from fudge.gnd import warning

        warnings = []

        if self.interpolation == standardsModule.interpolation.flatToken:
            warnings.append( warning.flatIncidentEnergyInterpolation( ) )

        for idx in range(len(self)):
            integral = self[idx].integrate()
            if abs(integral - 1.0) > info['normTolerance']:
                warnings.append( warning.unnormalizedDistribution( PQUModule.PQU(self[idx].value,self.axes[-1].unit),
                    idx, integral, self[idx] ) )

            if( self[idx].rangeMin < 0.0 ) :
                warnings.append( warning.negativeProbability( PQUModule.PQU(self[idx].value,self.axes[-1].unit),
                    value=self[idx].rangeMin, obj=self[idx] ) )

        return warnings

    def sqrtEp_AverageAtE( self, E ) :

        return( self.getAtEnergy( E ).integrateWithWeight_sqrt_x( ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        return( multiD_XYsModule.XYs2d.toPointwise_withLinearXYs( self, cls = XYs2d, **kwargs ) )

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :

        linear = self
        for xys in self :
            if( isinstance( xys, XYs1d ) ) :
                if( xys.interpolation not in [ standardsModule.interpolation.linlinToken, standardsModule.interpolation.flatToken ] ) :
                    linear = self.toPointwise_withLinearXYs( accuracy = XYsModule.defaultAccuracy, upperEps = 1e-8 )
                    break
            else :
                linear = self.toPointwise_withLinearXYs( accuracy = XYsModule.defaultAccuracy, upperEps = 1e-8 )
                break
        subform = XYs2d( axes = self.axes, interpolation = self.interpolation, 
                interpolationQualifier = self.interpolationQualifier )
        for xys in linear : subform.append( xs_pdf_cdf1d.fromXYs( xys, xys.value ) )
        return( subform )

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs1d, regions1d, xs_pdf_cdf1d ) )

class regions2d( subform, regionsModule.regions2d ) :

    def __init__( self, **kwargs ):

        regionsModule.regions2d.__init__( self, **kwargs )
        subform.__init__( self )

    def check( self, info ) :

        from fudge.gnd import warning

        warnings = []
        for idx, region in enumerate( self ):
            regionWarnings = region.check( info )
            if regionWarnings:
                warnings.append( warning.context("Region %d:" % idx, regionWarnings) )
        return warnings

    def toPointwise_withLinearXYs( self, **kwargs ) :

        return( regionsModule.regions2d.toPointwise_withLinearXYs( self, cls = XYs2d, **kwargs ) )

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :

        _regions2d = regions2d( axes = self.axes )
        for region in self : _regions2d.append( region.to_xs_pdf_cdf1d( style, tempInfo, indent ) )
        return( _regions2d )

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs2d, ) )

class energyFunctionalData( ancestryModule.ancestry ) :

    ancestryMembers = ( 'data', )

    def __init__( self, data ) :

        ancestryModule.ancestry.__init__( self )
        self.data = data
        data.setAncestor( self )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.data.convertUnits( unitMap )

    def copy( self ):

        return self.__class__( self.data.copy( ) )

    __copy__ = __deepcopy__ = copy

    def toXMLList( self, indent = '', **kwargs ) :

        xml = ['%s<%s>' % (indent, self.moniker)]
        xml += self.data.toXMLList( indent + '  ', **kwargs )
        xml[-1] += '</%s>' % self.moniker
        return xml

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):

        xPath.append( element.tag )
        subClass = {
            XYsModule.XYs1d.moniker         : XYsModule.XYs1d,
            regionsModule.regions1d.moniker : regionsModule.regions1d,
            xs_pdf_cdfModule.xs_pdf_cdf1d.moniker : xs_pdf_cdfModule.xs_pdf_cdf1d
        }.get( element[0].tag )
        if( subClass is None ) : raise Exception( "encountered unknown energy functional subform: %s" % element[0].tag )
        EFD = cls( subClass.parseXMLNode( element[0], xPath, linkData ) )
        xPath.pop()
        return EFD

class a( energyFunctionalData ) :

    moniker = 'a'

class b( energyFunctionalData ) :

    moniker = 'b'

class theta( energyFunctionalData ) :

    moniker = 'theta'

class g( energyFunctionalData ) :

    moniker = 'g'

class T_M( energyFunctionalData ) :

    moniker = 'T_M'

class functionalBase( subform ) :

    ancestryMembers = ( 'parameter1', 'parameter2' )

    def __init__( self, LF, U, parameter1, parameter2 = None ) :

        subform.__init__( self )
        if( U is not None ) :
            if( not( isinstance( U, physicalQuantityModule.U ) ) ) : raise TypeError( 'Invalid U type' )
        self.U = U
        self.LF = LF
        self.parameter1 = parameter1
        parameter1.setAncestor( self )
        self.parameter2 = parameter2
        if( parameter2 is not None ) : parameter2.setAncestor( self )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        if( self.U is not None ) : self.U.convertUnits( unitMap )
        self.parameter1.convertUnits( unitMap )
        if( self.parameter2 is not None ) : self.parameter2.convertUnits( unitMap )

    def copy( self ):

        U = self.U
        if( U is not None ) : U = self.U.copy( )
        if( self.parameter2 is None ) :
            return self.__class__( U, self.parameter1.copy( ) )
        else :
            return self.__class__( U, self.parameter1.copy( ), self.parameter2.copy( ) )

    __copy__ = __deepcopy__ = copy

    def check( self, info ):

        from fudge.gnd import warning

        warnings = []
        if( ( self.domainMin - self.U.value ) < 0 ) :
            warnings.append( warning.energyDistributionBadU( self ) )
        return( warnings )

    @property
    def domainMin( self ) :

        return( self.parameter1.data.domainMin )

    @property
    def domainMax( self ) :

        return( self.parameter1.data.domainMax )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        if( isinstance( self.parameter1.data, regionsModule.regions1d ) ) :
            Es = []
            for region in self.parameter1.data :
                Es = Es[:-1] + [ E for E, p in region ]
        else :
            Es = [ E for E, p in self.parameter1.data ]
        if( EMin is not None ) :
            if( EMin < Es[0] ) : Es.insert( 0, EMin )
        if( EMax is not None ) :
            if( EMax > Es[-1] ) : Es.append( EMax )
        return( Es )

    def toXMLList( self, indent = '', **kwargs ) :
        """Returns the xml string representation of self."""

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ self.XMLStartTagString( indent = indent ) ]
        if( self.LF == 12 ) : 
            xmlString += self.EFL.toXMLList( indent2, **kwargs )
            xmlString += self.EFH.toXMLList( indent2, **kwargs )
        else :
            xmlString += self.U.toXMLList( indent2, **kwargs )
        xmlString += self.parameter1.toXMLList( indent2, **kwargs )
        if( not( self.parameter2 is None ) ) : xmlString += self.parameter2.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

class generalEvaporationSpectrum( functionalBase ) :

    moniker = 'generalEvaporation'

    def __init__( self, U, thetas, gs ) :

        functionalBase.__init__( self, 5, U, thetas, gs )

    def averageEp( self, E ) :

        return( self.parameter1.data.evaluate( E ) * self.parameter2.data.integrateWithWeight_x( ) )

    def sqrtEp_AverageAtE( self, E ) :

        return( math.sqrt( self.parameter1.data.evaluate( E ) ) * self.parameter2.data ).integrateWithWeight_sqrt_x( )

    def isLinear( self, qualifierOk = False, flatIsOk = False ) :
        """
        Returns the results of isLinear called on the axes of g(E'|E).
        """

        return( self.parameter2.axes.isLinear( qualifierOk = qualifierOk, flatIsOk = flatIsOk ) )

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :

        _form = generalEvaporationSpectrum( self.U, thetas = self.parameter1.copy( ),
            gs = g( xs_pdf_cdf1d.fromXYs( self.parameter2.data ) ) )
        return( _form )

    def toPointwise_withLinearXYs( self, **kwargs  ) :

        pwl = XYs2d( axes = defaultAxes( self.parameter1.data.domainUnit ) )
        thetas = self.parameter1.data.toPointwise_withLinearXYs( **kwargs )
        gs = self.parameter2.data.toPointwise_withLinearXYs( **kwargs )
        for E_in, theta in thetas :
            data = [ [ theta * x, y / theta ] for x, y in gs ]
            data = XYs1d( data, value = E_in )
            data.normalize( insitu = True )
            pwl.append( data )
        return( pwl )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :
        """Translate <generalEvaporation> element from xml."""

        xPath.append( element.tag )
        theta_ = theta.parseXMLNode( element.find(theta.moniker), xPath, linkData )
        g_ = g.parseXMLNode( element.find(g.moniker), xPath, linkData )
        U = physicalQuantityModule.U.parseXMLNode( element.find( 'U' ), xPath, linkData )
        GES = generalEvaporationSpectrum( U, theta_, g_ )
        xPath.pop()
        return GES

class simpleMaxwellianFissionSpectrum( functionalBase ) :

    moniker = 'simpleMaxwellianFission'

    def __init__( self, U, thetas ) :

        functionalBase.__init__( self, 7, U, thetas )

    def averageEp( self, E ) :

        theta = self.parameter1.data.evaluate( E )
        a = ( E - self.U.value ) / theta
        if( a < 1e-4 ) : return( theta * a * ( 1575. - a * ( 180. + 8 * a ) ) / 2625. )
        sqrt_a = math.sqrt( a )
        exp_a = math.exp( -a )
        erf_sqrt_a = math.sqrt( math.pi ) * math.erf( sqrt_a )
        return( theta * ( 0.75 * erf_sqrt_a - sqrt_a * ( 1.5 + a ) * exp_a ) / ( 0.5 * erf_sqrt_a - sqrt_a * exp_a ) )

    def sqrtEp_AverageAtE( self, E ) :

        theta = self.parameter1.data.evaluate( E )
        a = ( E - self.U.value ) / theta
        if( a < 1e-4 ) : return( math.sqrt( theta * a ) * ( 2100. - a * ( 140. + 9. * a ) ) / 2800. )
        sqrt_a = math.sqrt( a )
        exp_a = math.exp( -a )
        erf_sqrt_a = math.sqrt( math.pi ) * math.erf( sqrt_a )
        return( math.sqrt( theta ) * ( 1 - ( 1. + a ) * exp_a ) / ( 0.5 * erf_sqrt_a - sqrt_a * exp_a ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        def evaluateAtX( self, x ) :

            return( math.sqrt( x ) * math.exp( -x / self.p1 ) )

        ef = energyFunctionalDataToPointwise( self, evaluateAtX )
        return( ef.toPointwise_withLinearXYs( **kwargs ) )

    @staticmethod
    def parseXMLNode( MFelement, xPath, linkData ) :

        xPath.append( MFelement.tag )
        theta_ = theta.parseXMLNode( MFelement.find(theta.moniker), xPath, linkData )
        U = physicalQuantityModule.U.parseXMLNode( MFelement.find( 'U' ), xPath, linkData )
        SMF = simpleMaxwellianFissionSpectrum( U, theta_ )
        xPath.pop()
        return SMF

class evaporationSpectrum( functionalBase ) :

    moniker = 'evaporation'

    def __init__( self, U, thetas ) :

        functionalBase.__init__( self, 9, U, thetas )

    def averageEp( self, E ) :

        if( isinstance( self.parameter1.data, regionsModule.regions1d ) ) :
            for region in self.parameter1.data :
                if( E <= region[-1][0] ) : break
            theta = region.evaluate( E )
        else :
            theta = self.parameter1.data.evaluate( E )
        a = ( E - self.U.value ) / theta
        if( a < 1e-4 ) : return( theta * a * ( 180. - a * ( 15. + a ) ) / 270. )
        exp_a = math.exp( -a )
        return( theta * ( 2. - a**2 * exp_a / ( 1. - ( 1. + a ) * exp_a ) ) )

    def sqrtEp_AverageAtE( self, E ) :

        theta = self.parameter1.data.evaluate( E )
        a = ( E - self.U.value ) / theta
        if( a < 1e-4 ) : return( math.sqrt( theta * a ) * ( 252. - a * ( 12. + a ) ) / 315. )
        sqrt_a = math.sqrt( a )
        exp_a = math.exp( -a )
        return( math.sqrt( theta ) * ( 1.32934038817913702 * math.erf( sqrt_a ) - sqrt_a * ( 1.5 + a ) * exp_a ) / ( 1. - ( 1. + a ) * exp_a ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        def evaluateAtX( self, x ) :

            return( x * math.exp( -x / self.p1 ) )

        ef = energyFunctionalDataToPointwise( self, evaluateAtX )
        return( ef.toPointwise_withLinearXYs( **kwargs ) )

    @staticmethod
    def parseXMLNode( evapElement, xPath, linkData ) :

        xPath.append( evapElement.tag )
        theta_ = theta.parseXMLNode( evapElement.find(theta.moniker), xPath, linkData )
        U = physicalQuantityModule.U.parseXMLNode( evapElement.find( 'U' ), xPath, linkData )
        ES = evaporationSpectrum( U, theta_ )
        xPath.pop()
        return ES

class WattSpectrum( functionalBase ) :

    moniker = 'Watt'

    def __init__( self, U, a, b ) :

        functionalBase.__init__( self, 11, U, a, b )

    def averageEp( self, E ) :

        a, b = self.parameter1.data.evaluate( E ), self.parameter2.data.evaluate( E )
        domainMax_a  = ( E - self.U.value ) / a
        domainMax_b  = math.sqrt( b * ( E - self.U.value ) )
        ab = a * b
        sqrt_ab = math.sqrt( ab )
        I = 0.25 * math.sqrt( math.pi * ab ) * math.exp( 0.25 * ab ) * \
            ( math.erf( math.sqrt( domainMax_a ) - 0.5 * sqrt_ab ) + math.erf( math.sqrt( domainMax_a ) + 0.5 * sqrt_ab ) ) \
            - math.exp( -domainMax_a ) * math.sinh( domainMax_b )
        EI = a * math.sqrt( math.pi * ab ) * ( ab + 6 ) * math.exp( 0.25 * ab ) * \
            ( math.erf( math.sqrt( domainMax_a ) - 0.5 * sqrt_ab ) + math.erf( math.sqrt( domainMax_a ) + 0.5 * sqrt_ab ) ) \
            - 0.25 * a * math.exp( -domainMax_a ) * math.sinh( domainMax_b ) * ( 2. * domainMax_b + ab + 4. + 4. * domainMax_a )
        return( EI / ( 16 * I ) )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        aMin, aMax = self.parameter1.data.domainMin, self.parameter1.data.domainMax
        if( EMin is None ) : EMin = aMin
        if( EMax is None ) : EMax = aMax
        if( EMin < aMin ) : EMin = aMin
        if( EMax < aMax ) : EMax = aMax
        Es = [ EMin, EMax ]
        for E, a in self.parameter1.data :
            if( E <= EMin ) : continue
            if( E >= EMax ) : continue
            if( E not in Es ) : Es.append( E )
        for E, b in self.parameter2.data :
            if( E <= EMin ) : continue
            if( E >= EMax ) : continue
            if( E not in Es ) : Es.append( E )
        Es.sort( )
        return( Es )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        def evaluateAtX( self, x ) :

            return( math.exp( -x / self.p1 ) * math.sinh( math.sqrt( self.p2 * x ) ) )

        ef = energyFunctionalDataToPointwise( self, evaluateAtX )
        return( ef.toPointwise_withLinearXYs( **kwargs ) )

    @staticmethod
    def parseXMLNode( WattElement, xPath, linkData ):
        """Translate <Watt> element from xml."""

        xPath.append( WattElement.tag )
        _a = a.parseXMLNode( WattElement.find(a.moniker), xPath, linkData )
        _b = b.parseXMLNode( WattElement.find(b.moniker), xPath, linkData )
        U = physicalQuantityModule.U.parseXMLNode( WattElement.find( 'U' ), xPath, linkData )
        WS = WattSpectrum( U, _a, _b )
        xPath.pop()
        return WS

class MadlandNix( functionalBase ) :

    moniker = 'MadlandNix'

    def __init__( self, EFL, EFH, Ts ) :

        functionalBase.__init__( self, 12, None, Ts )
        self.EFL = EFL
        self.EFH = EFH

    def copy( self ):

        return MadlandNix( self.EFL.copy( ), self.EFH.copy( ), self.parameter1.copy( ) )

    __copy__ = __deepcopy__ = copy

    def check( self, info ) :

        from fudge.gnd import warning

        warnings = []
        if self.EFL.value <= 0 or self.EFH.value <= 0 or self.parameter1.data.rangeMin <= 0:
            warnings.append( warning.MadlandNixBadParameters( self.EFL, self.EFH, self.parameter1.data.rangeMin, self ) )

        return warnings

    def averageEp( self, E ) :

        unit = self.parameter1.data.axes[-1].unit
        return( 0.5 * ( self.EFL.getValueAs( unit ) + self.EFH.getValueAs( unit ) ) + 4. * self.parameter1.data.evaluate( E ) / 3. )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        functionalBase.convertUnits( self, unitMap )
        self.EFL.convertUnits( unitMap )
        self.EFH.convertUnits( unitMap )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        return( [ x for x, y in self.parameter1.data ] )

    def toPointwise_withLinearXYs( self, **kwargs ) :

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

        accuracy = kwargs.get( 'accuracy', XYsModule.defaultAccuracy )
        axes = defaultAxes( energyUnit = self.parameter1.data.axes[0].unit )
        pwl = XYs2d( axes = axes )
        E_in_unit = self.parameter1.data.axes[-1].unit
        EFL, EFH = self.EFL.getValueAs( E_in_unit ), self.EFH.getValueAs( E_in_unit )
        factor = PQUModule.PQU( 1, 'eV' ).getValueAs( E_in_unit )
        _xs = [ 1e-5, 1e-3, 1e-1, 1e1, 1e3, 1e5, 3e7 ]          # List of energies in eV, hence need for factor from eV to E_in_unit.
        xs = [ factor * x for x in _xs ]
        for E, T_M in self.parameter1.data :        # This logic ignores the interpolation of parameter1 as the only two subforms in ENDF/B-VII shows 
            parameters = [ EFL, EFH, T_M ]          # that linear-linear is better than the 'log-log' given in the ENDF/B-VII/data.
            g_Ep = XYs1d.createFromFunction( axes, xs, MadlandNixFunc, parameters, accuracy, biSectionMax = 12 )
            g_Ep.value = E                          # ????????? Class XYs1d does not the a proper setValue method. One should be added.
            g_Ep.normalize( insitu = True )
            pwl.append( g_Ep )
        return( pwl )

    def to_xs_pdf_cdf1d( self, style, tempInfo, indent ) :

        linear = self.toPointwise_withLinearXYs( )
        return( linear.to_xs_pdf_cdf1d( style, tempInfo, indent ) )

    @staticmethod
    def parseXMLNode( MNelement, xPath, linkData ):
        """Translate <MadlandNix> element from xml."""

        xPath.append( MNelement.tag )
        tm = T_M.parseXMLNode( MNelement.find(T_M.moniker), xPath, linkData )
        EFL = physicalQuantityModule.EFL.parseXMLNode( MNelement.find( "EFL" ), xPath, linkData )
        EFH = physicalQuantityModule.EFH.parseXMLNode( MNelement.find( "EFH" ), xPath, linkData )
        MN = MadlandNix( EFL, EFH, tm )
        xPath.pop()
        return MN

class NBodyPhaseSpace( subform ) :

    moniker = 'NBodyPhaseSpace'
    ancestryMembers = ( '', )

    def __init__( self, numberOfProducts, numberOfProductsMasses ) :

        subform.__init__( self )
        self.numberOfProducts = numberOfProducts
        self.numberOfProductsMasses = numberOfProductsMasses

    def convertUnits( self, unitMap ) :

        pass

    def copy( self ) :

        return NBodyPhaseSpace( self.numberOfProducts, self.numberOfProductsMasses.copy( ) )

    __copy__ = __deepcopy__ = copy

    def check( self, info ) :

        #if ( abs( self.numberOfProductsMasses - info['reactionSuite'].products['n'].getMass('amu') * self.numberOfProducts ) >
        #        self.numberOfProductsMasses * 0.1 ) :    # return warning?

        return []

    def averageEp( self, E, massUnit, projectileMass, targetMass, productMass, Q ) :
        """
        Calculate the average energy of the product in the center-of-mass frame for projectile energy E. This method only works for
        a one-step reaction.
        """

        M = self.numberOfProductsMasses.getValueAs( massUnit )
        Ea = targetMass / ( targetMass + projectileMass ) * E + Q
        return( Ea * ( M - productMass ) / ( M * ( self.numberOfProducts - 1 ) ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :

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
        numberOfProductsMasses, massUnit = self.numberOfProductsMasses.value, self.numberOfProductsMasses.unit
        productMass = p.getMass( massUnit )

        r = self.findClassInAncestry( reaction.reaction )
        EMin = r.domainMin
        EMax = r.domainMax
        energyUnit = r.domainUnit

        reactionSuite = self.findClassInAncestry( reactionSuite.reactionSuite )
        projectile = reactionSuite.PoPs[reactionSuite.projectile]
        projectileMass = projectile.mass[0].float( massUnit )
        target = reactionSuite.PoPs[reactionSuite.target]
        targetMass = target.mass[0].float( massUnit )

        c = self.findClassInAncestry( channels.channel )
        Q = c.Q.getConstant( )

        axes = defaultAxes( energyUnit )
        pwl = XYs2d( axes=axes )

        accuracy = kwargs.get( 'accuracy', XYsModule.defaultAccuracy )

        t = tester( accuracy, 1e-10, self.numberOfProducts )
        n = 21
        f = math.pow( EMax / EMin, 1. / n )
        E_ins = [ EMin * f**idx for idx in xrange( n ) ]
        E_ins[-1] = EMax        # Fix possible round off issue.
        for idx, E_in in enumerate( E_ins ) :
            Ea = targetMass / ( targetMass + projectileMass ) * E_in + Q
            EMax_i = Ea * ( numberOfProductsMasses - productMass ) / numberOfProductsMasses
            if( EMax_i < 0 ) : EMax_i = 1e-5                # This is a kludge
            t.setEMax_i( EMax_i )
            t.absoluteTolerance = 1e-10 * t.evaluateAtX( 0.5 * EMax_i )
            data = fudgemath.thickenXYList( [ [ 0., 0. ], [ EMax_i, 0. ] ], t, biSectionMax = 10 )
            data = XYs1d( data, value = E_in )
            data.normalize( insitu = True )
            pwl.append( data )
        return( pwl )

    def toXMLList( self, indent = '', **kwargs ) :
        """Returns the xml string representation of self."""

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ self.XMLStartTagString( indent = indent, 
                extraAttributesAsStrings = ' numberOfProducts="%s"' % ( self.numberOfProducts ), emptyTag = False ) ]
        xmlString += self.numberOfProductsMasses.toXMLList( indent2, **kwargs )
        xmlString[-1] += "</%s>" % self.moniker
        return( xmlString )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        mass = physicalQuantityModule.mass.parseXMLNode( element.find( 'mass' ), xPath, linkData )
        return NBodyPhaseSpace( int( element.get( "numberOfProducts" ) ), mass )

class weighted( ancestryModule.ancestry ) :

    moniker = 'weighted'
    ancestryMembers = ( 'weight', 'functional' )

    def __init__( self, weight, functional ) :

        ancestryModule.ancestry.__init__( self )

        self.weight = weight
        weight.setAncestor( self )
# BRB6 weight should be its own class so that label need not be set here.
        self.weight.label = 'weight'

        self.functional = functional
        functional.setAncestor( self )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        self.weight.convertUnits( unitMap )
        self.functional.convertUnits( unitMap )

    def copy( self ):

        return weighted( self.weight.copy( ), self.functional.copy( ) )

    __copy__ = __deepcopy__ = copy

    def checkProductFrame( self ) :
        "Calls checkProductFrame on self's functional."

        self.functional.checkProductFrame( )

    def averageEp( self, E ) :

        return( self.weight.evaluate( E ) * self.functional.averageEp( E ) )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        energyArray = self.functional.getEnergyArray( EMin, EMax )
        for x, y in self.weight :
            if( x not in energyArray ) : energyArray.append( x )
        energyArray.sort( )
        return( energyArray )

    def toXMLList( self, indent = '', **kwargs ) :
        """Returns the xml string representation of self."""

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ '%s<%s>' %(  indent, self.moniker ) ]
        xmlString += self.weight.toXMLList( indent2, **kwargs )
        xmlString += self.functional.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

class weightedFunctionals( subform ) :

    moniker = 'weightedFunctionals'
    ancestryMembers = ( '[weights', )

    def __init__( self ) :

        subform.__init__( self )
        self.weights = []

    def __len__( self ) :

        return( len( self.weights ) )

    def __getitem__( self, i ) :

        return( self.weights[i] )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        for weight in self.weights : weight.convertUnits( unitMap )

    def copy( self ) :

        newWF = weightedFunctionals( )
        for weight in self : newWF.addWeight( weight.copy( ) )
        return( newWF )

    __copy__ = __deepcopy__ = copy

    def addWeight( self, weight ) :

        if( not( isinstance( weight, weighted ) ) ) : raise Exception( 'Invalid weight of type %s' % brb.getType( weight ) )
        self.weights.append( weight )
        weight.setAncestor( self )

    def check( self, info ) :

        from fudge.gnd import warning

        warnings = []
        totalWeight = sum( [w.weight for w in self.weights] )
        if (totalWeight.rangeMin != 1.0) and (totalWeight.rangeMax != 1.0):
            warnings.append( warning.weightsDontSumToOne( obj=self ) )

        for weight in self:
            warnings += weight.functional.check( info )

        return warnings

    def averageEp( self, E ) :

        Ep = 0
        for weight in self : Ep += weight.averageEp( E )
        return( Ep )

    def sqrtEp_AverageAtE( self, E ) :
        """
        This method has not been implemented. It returns None so the method uncorrelated.calculateAverageProductData will still work
        when calculating the momentum deposition for isotropic scattering in the lab frame, but will execute a raise otherwise.
        """

        return( None )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        energyArray = []
        for weight in self :
            for energy in weight.getEnergyArray( EMin, EMax ) :
                if( energy not in energyArray ) : energyArray.append( energy )
        energyArray.sort( )
        return( energyArray )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        if( len( self ) > 2 ) : raise Exception( 'more than two weights currently not supported' )
        E_ins, data = [], []
        for weighted_ in self :
            w = weighted_.weight.toPointwise_withLinearXYs( **kwargs )
            e = weighted_.functional.toPointwise_withLinearXYs( **kwargs )
            data.append( [ w, e ] )
            for x, y in w :
                if( x not in E_ins ) : E_ins.append( x )
            for x in e :
                if( x.value not in E_ins ) : E_ins.append( x.value )
        E_ins.sort( )
        pwl = XYs2d( axes = defaultAxes( self[0].weight.domainUnit ) )
        for E_in in E_ins :
            wv1, ev1 = data[0]
            wv2, ev2 = data[1]
            w1, w2 = wv1.evaluate( E_in ), wv2.evaluate( E_in )
            if( w1 == 1 ) :
                e = ev1.getSpectrumAtEnergy( E_in )
            elif( w2 == 1 ) :
                e = ev2.getSpectrumAtEnergy( E_in )
            else :
                raise Exception( 'One weight must be zero: %s, %s' % ( w1, w2 ) )
            e.value = E_in
            e.normalize( insitu = True )
            pwl.append( e )
        return( pwl )

    def toXMLList( self, indent = '', **kwargs ) :
        """Returns the xml string representation of self."""

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ self.XMLStartTagString( indent = indent ) ]
        for _subform in self.weights : xmlString += _subform.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    @staticmethod
    def parseXMLNode( WFelement, xPath, linkData ) :

        xPath.append( WFelement.tag )
        WF = weightedFunctionals( )
        for weightSection in WFelement :
            weights, functional = weightSection
            _weight = XYsModule.XYs1d.parseXMLNode( weights, xPath, linkData )
            subformClass = {
                    XYs2d.moniker                           : XYs2d,
                    generalEvaporationSpectrum.moniker      : generalEvaporationSpectrum,
                    WattSpectrum.moniker                    : WattSpectrum,
                    MadlandNix.moniker                      : MadlandNix,
                    simpleMaxwellianFissionSpectrum.moniker : simpleMaxwellianFissionSpectrum,
                    evaporationSpectrum.moniker             : evaporationSpectrum,
                }.get( functional.tag )
            functional = subformClass.parseXMLNode( functional, xPath, linkData )
            WF.addWeight( weighted( _weight, functional ) )
        xPath.pop( )
        return WF

class energyFunctionalDataToPointwise :

    def __init__( self, data, evaluateAtX ) :

        self.data = data
        self.evaluateAtX = evaluateAtX

    def toPointwise_withLinearXYs( self, **kwargs ) :

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

        accuracy = xDataBaseModule.getArguments( kwargs, { 'accuracy' : XYsModule.defaultAccuracy } )['accuracy']

        parameter1 = self.data.parameter1.data.toPointwise_withLinearXYs( **kwargs )
        axes = defaultAxes( self.data.domainUnit )
        pwl = XYs2d( axes = axes )
# BRB6 hardwired
        one_eV = PQUModule.PQU( '1 eV' ).getValueAs( axes[-1].unit )

        t = tester( accuracy, 1e-10, self.evaluateAtX )
        parameter2 = self.data.parameter2
        if( parameter2 is not None ) : parameter2 = parameter2.data.toPointwise_withLinearXYs( **kwargs )
        for i, E_in_p1 in enumerate( parameter1 ) :
            E_in, p1 = E_in_p1
            EpMax = E_in - self.data.U.value
            EpMin = 0.              # Only used with debugging.
            if( EpMax == 0. ) :
                if( i != 0 ) : raise Exception( "i = %d, E_in = %s, U = %s" % ( E_in, self.data.U.value ) )
                EpMax = E_in * 1e-6
                if( EpMax == 0. ) :
                    EpMax = parameter1[1][1] * 1e-6         # This and the next line are arbitary.
                if( EpMax > 1e-3 ) : EpMax = 2e-5
                data = [ [ EpMin, 1. ], [ EpMax, 1e-10 ] ]
            else :
                t.setParameter1( p1 )
                if( parameter2 is not None ) : t.setParameter2( parameter2.evaluate( E_in ) )
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
            data = XYs1d( data = data, value = E_in )
            data.normalize( insitu = True )
            pwl.append( data )
        return( pwl )

class form( baseModule.form ) :
    "I think this is only used to convert to ENDF6. If so, should be moved to site_packages/legacy/toENDF6/productData/distributions/energy.py"

    moniker = 'energy'
    subformAttributes = ( 'energySubform', )

    def __init__( self, label, productFrame, energySubform ) :

        if( not( isinstance( energySubform, subform ) ) ) : raise TypeError( 'instance is not an energy subform' )
        baseModule.form.__init__( self, label, productFrame, ( energySubform, ) )

    @staticmethod
    def parseXMLNode( energyElement, xPath, linkData ) :
        """Translate energy component from xml."""

        subformClass = {
                discreteGamma.moniker :                     discreteGamma,
                primaryGamma.moniker :                      primaryGamma,
                XYs2d.moniker :                             XYs2d,
                regions2d.moniker :                         regions2d,
                generalEvaporationSpectrum.moniker :        generalEvaporationSpectrum,
                WattSpectrum.moniker :                      WattSpectrum,
                MadlandNix.moniker :                        MadlandNix,
                simpleMaxwellianFissionSpectrum.moniker :   simpleMaxwellianFissionSpectrum,
                evaporationSpectrum.moniker :               evaporationSpectrum,
                weightedFunctionals.moniker :               weightedFunctionals,
                NBodyPhaseSpace.moniker :                   NBodyPhaseSpace,
            }.get( energyElement.tag )
        if( subformClass is None ) : raise Exception( "encountered unknown energy subform: %s" % energyElement.tag )
        energyForm = subformClass.parseXMLNode( energyElement, xPath, linkData )
        return( energyForm )
