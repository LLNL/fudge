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
import base, miscellaneous
import fudge
from pqu import PQU
from fudge.core.utilities import brb
import xData.ancestry as ancestryModule

from xData import standards as standardsModule
from xData import axes as axesModule
from xData import XYs as XYsModule
from xData import multiD_XYs as multiD_XYsModule
from xData import regions as regionsModule

from . import base as baseModule

__metaclass__ = type


class XYs1d( XYsModule.XYs1d ) :

    def averageEnergy( self ) :

        allowedInterpolations = [ standardsModule.interpolation.linlinToken,
                                      standardsModule.interpolation.flatToken ]
        xys = self.changeInterpolationIfNeeded( allowedInterpolations = allowedInterpolations )
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

class subform( baseModule.subform ) :
    """Abstract base class for energy forms."""

    pass

class constant( subform ) :

    moniker = 'constant'

    def __init__( self, value ) :

        subform.__init__( self )
        self.value = PQU.PQU( value )

    def check( self, info ) :

        from fudge.gnd import warning

        warnings = []
        if self.value <= PQU.PQU(0,'eV'):
            warnings.append( warning.negativeDiscreteGammaEnergy() )
        return warnings

    def copy( self ):

        return constant( self.value )

    __copy__ = __deepcopy__ = copy

    def getEnergyArray( self, EMin = None, EMax = None ) :

        return( [ EMin, EMax ] )

    def averageEp( self, E ) :

        return( float( self.value ) )

    def toXMLList( self, indent = '', **kwargs ) :

        attributeStr = ''
        if( self.label is not None ) : attributeStr += ' label="%s"' % self.label
        return( [ '%s<%s%s value="%s"/>' % ( indent, self.moniker, attributeStr, self.value ) ] )

    @staticmethod
    def parseXMLNode( energyElement, xPath, linkData ) :

        value = energyElement.get( 'value' )
        return( constant( value ) )

class primaryGamma( subform ) :

    moniker = 'primaryGamma'

    def __init__( self, value ) :

        subform.__init__( self )
        self.value = PQU.PQU( value )
        self.__massRatio = None     # In ENDF lingo this is AWR / ( AWR + 1 ).

    @property
    def massRatio( self ) :

        if( self.__massRatio is None ) :
            self.__massRatio = self.findAttributeInAncestry( "getMassRatio" )( )
        return self.__massRatio

    def check( self, info ) :

        from fudge.gnd import warning

        warnings = []
        Qvalue = self.findAttributeInAncestry('getQ')('eV')
        if self.value.getValueAs('eV') > Qvalue:
            warnings.append( warning.primaryGammaEnergyTooLarge( self.value,
                            100 * self.value.getValueAs('eV') / Qvalue ) )
        return warnings

    def copy( self ):

        return primaryGamma( self.value )

    __copy__ = __deepcopy__ = copy

    def getEnergyArray( self, EMin = None, EMax = None ) :

        return( [ EMin, EMax ] )

    def averageEp( self, E ) :

        return( float( self.value ) + self.massRatio * E )

    def toXMLList( self, indent = '', **kwargs ) :

        attributeStr = ''
        if( self.label is not None ) : attributeStr += ' label="%s"' % self.label
        return( [ '%s<%s%s value="%s"/>' % ( indent, self.moniker, attributeStr, self.value ) ] )

    @staticmethod
    def parseXMLNode( energyElement, xPath, linkData ) :

        value = energyElement.get( 'value' )
        return primaryGamma( value )

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
            EpLowerMin, EpLowerMax = lower.domain( )
            EpUpperMin, EpUpperMax = upper.domain( )
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
        for idx in range(len(self)):
            integral = self[idx].integrate()
            if abs(integral - 1.0) > info['normTolerance']:
                warnings.append( warning.unnormalizedDistribution( PQU.PQU(self[idx].value,self.axes[-1].unit),
                    idx, integral, self[idx] ) )

            if self[idx].rangeMin() < 0.0:
                warnings.append( warning.negativeProbability( PQU.PQU(self[idx].value,self.axes[-1].unit),
                    value=self[idx].rangeMin(), obj=self[idx] ) )

        return warnings

    def sqrtEp_AverageAtE( self, E ) :

        return( self.getAtEnergy( E ).integrateWithWeight_sqrt_x( ) )

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        return( multiD_XYsModule.XYs2d.toPointwise_withLinearXYs( self, accuracy, lowerEps = lowerEps,
            upperEps = upperEps, cls = XYs2d ) )

    @staticmethod
    def defaultAxes( energyUnit = 'eV' ) :

        axes = axesModule.axes( rank = 3 )
        axes[0] = axesModule.axis( 'P(energy_out|energy_in)', 0, '1/' + energyUnit )
        axes[1] = axesModule.axis( 'energy_out', 1, energyUnit )
        axes[2] = axesModule.axis( 'energy_in',  2, energyUnit )
        return( axes )

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs1d, regions1d ) )

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

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        return( regionsModule.regions2d.toPointwise_withLinearXYs( self, accuracy, lowerEps = lowerEps,
            upperEps = upperEps, cls = XYs2d ) )

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs2d, ) )

    @staticmethod
    def defaultAxes( energyUnit = 'eV' ) :

        return( XYs2d.defaultAxes( energyUnit = energyUnit ) )

class energyFunctionalData( ancestryModule.ancestry ) :

    def __init__( self, data ) :

        ancestryModule.ancestry.__init__( self )
        self.data = data

    @classmethod
    def copy( cls, self ):

        return cls( self.data )

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

    def __init__( self, LF, U, parameter1, parameter2 = None ) :

        subform.__init__( self )
        self.U = U
        self.LF = LF
        self.parameter1 = parameter1
        self.parameter2 = parameter2

    def copy( self ):

        if( self.parameter2 is None ) :
            return self.__class__( self.U, self.parameter1 )
        else :
            return self.__class__( self.U, self.parameter1, self.parameter2 )

    __copy__ = __deepcopy__ = copy

    def check( self, info ):

        from fudge.gnd import warning

        warnings = []
        if( ( self.domainMin( unitTo = 'eV' ) - self.U.getValueAs( 'eV' ) ) < 0 ) :
            warnings.append( warning.energyDistributionBadU( self ) )
        return( warnings )

    def domainMin( self, unitTo = None, asPQU = False ) :

        return( self.parameter1.data.domainMin( unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        return( self.parameter1.data.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def domain( self, unitTo = None, asPQU = False ):

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ), self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

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

        if( self.LF == 12 ) : 
            EFL = self.EFL.toString( keepPeriod = False )
            EFH = self.EFH.toString( keepPeriod = False )
            qualifiers = ' EFL="%s" EFH="%s"' % ( EFL, EFH )
        else :
            qualifiers = ' U="%s"' % self.U.toString( keepPeriod = False )
        xmlString = [ self.XMLStartTagString( indent = indent, extraAttributesAsStrings = qualifiers ) ]
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

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        pwl = XYs2d( axes = XYs2d.defaultAxes( ) )
        thetas = self.parameter1.data.toPointwise_withLinearXYs( accuracy = None, lowerEps = lowerEps, upperEps = upperEps )
        gs = self.parameter2.data.toPointwise_withLinearXYs( accuracy = None, lowerEps = lowerEps, upperEps = upperEps )
        for E_in, theta in thetas :
            data = [ [ theta * x, y / theta ] for x, y in gs ]
            data = XYs1d( data, accuracy = accuracy, value = E_in )
            data.normalize( insitu = True )
            pwl.append( data )
        return( pwl )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :
        """Translate <generalEvaporation> element from xml."""

        xPath.append( element.tag )
        theta_ = theta.parseXMLNode( element.find(theta.moniker), xPath, linkData )
        g_ = g.parseXMLNode( element.find(g.moniker), xPath, linkData )
        U = PQU.PQU( element.get( "U" ) )
        GES = generalEvaporationSpectrum( U, theta_, g_ )
        xPath.pop()
        return GES

class simpleMaxwellianFissionSpectrum( functionalBase ) :

    moniker = 'simpleMaxwellianFission'

    def __init__( self, U, thetas ) :

        functionalBase.__init__( self, 7, U, thetas )

    def averageEp( self, E ) :

        theta = self.parameter1.data.evaluate( E )
        a = ( E - self.U.getValue( ) ) / theta
        if( a < 1e-4 ) : return( theta * a * ( 1575. - a * ( 180. + 8 * a ) ) / 2625. )
        sqrt_a = math.sqrt( a )
        exp_a = math.exp( -a )
        erf_sqrt_a = math.sqrt( math.pi ) * math.erf( sqrt_a )
        return( theta * ( 0.75 * erf_sqrt_a - sqrt_a * ( 1.5 + a ) * exp_a ) / ( 0.5 * erf_sqrt_a - sqrt_a * exp_a ) )

    def sqrtEp_AverageAtE( self, E ) :

        theta = self.parameter1.data.evaluate( E )
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
    def parseXMLNode( MFelement, xPath, linkData ) :

        xPath.append( MFelement.tag )
        theta_ = theta.parseXMLNode( MFelement.find(theta.moniker), xPath, linkData )
        U = PQU.PQU( MFelement.get("U") )
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
        a = ( E - self.U.getValue( ) ) / theta
        if( a < 1e-4 ) : return( theta * a * ( 180. - a * ( 15. + a ) ) / 270. )
        exp_a = math.exp( -a )
        return( theta * ( 2. - a**2 * exp_a / ( 1. - ( 1. + a ) * exp_a ) ) )

    def sqrtEp_AverageAtE( self, E ) :

        theta = self.parameter1.data.evaluate( E )
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
    def parseXMLNode( evapElement, xPath, linkData ) :

        xPath.append( evapElement.tag )
        theta_ = theta.parseXMLNode( evapElement.find(theta.moniker), xPath, linkData )
        U = PQU.PQU( evapElement.get( "U" ) )
        ES = evaporationSpectrum( U, theta_ )
        xPath.pop()
        return ES

class WattSpectrum( functionalBase ) :

    moniker = 'Watt'

    def __init__( self, U, a, b ) :

        functionalBase.__init__( self, 11, U, a, b )

    def averageEp( self, E ) :

        a, b = self.parameter1.data.evaluate( E ), self.parameter2.data.evaluate( E )
        domainMax_a  = ( E - self.U.getValue( ) ) / a
        domainMax_b  = math.sqrt( b * ( E - self.U.getValue( ) ) )
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

        aMin, aMax = self.parameter1.data.domain( )
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

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        def evaluateAtX( self, x ) :

            return( math.exp( -x / self.p1 ) * math.sinh( math.sqrt( self.p2 * x ) ) )

        ef = energyFunctionalDataToPointwise( self, evaluateAtX )
        return( ef.toPointwise_withLinearXYs( accuracy, lowerEps = lowerEps, upperEps = upperEps ) )

    @staticmethod
    def parseXMLNode( WattElement, xPath, linkData ):
        """Translate <Watt> element from xml."""

        xPath.append( WattElement.tag )
        _a = a.parseXMLNode( WattElement.find(a.moniker), xPath, linkData )
        _b = b.parseXMLNode( WattElement.find(b.moniker), xPath, linkData )
        U = PQU.PQU( WattElement.get( "U" ) )
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

        return MadlandNix( self.EFL, self.EFH, self.parameter1 )

    __copy__ = __deepcopy__ = copy

    def check( self, info ) :

        from fudge.gnd import warning

        warnings = []
        if self.EFL.value <= 0 or self.EFH.value <= 0 or self.parameter1.data.rangeMin() <= 0:
            warnings.append( warning.MadlandNixBadParameters( self.EFL, self.EFH, self.parameter1.data.rangeMin(), self ) )

        return warnings

    def averageEp( self, E ) :

        unit = self.parameter1.data.axes[-1].unit
        return( 0.5 * ( self.EFL.getValueAs( unit ) + self.EFH.getValueAs( unit ) ) + 4. * self.parameter1.data.evaluate( E ) / 3. )

    def getEnergyArray( self, EMin = None, EMax = None ) :

        return( [ x for x, y in self.parameter1.data ] )

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

        if( accuracy is None ) : accuracy = 1e-3
        axes = XYs2d.defaultAxes( energyUnit = self.parameter1.data.axes[0].unit )
        pwl = XYs2d( axes = axes )
        E_in_unit = self.parameter1.data.axes[-1].unit
        EFL, EFH = self.EFL.getValueAs( E_in_unit ), self.EFH.getValueAs( E_in_unit )
        factor = PQU.PQU( 1, 'eV' ).getValueAs( E_in_unit )
        xs_ = [ 1e-5, 1e-3, 1e-1, 1e1, 1e3, 1e5, 3e7 ]
        xs = [ factor * x for x in xs_ ]
        axes1d = axesModule.axes( )
        axes1d[0] = axes[0]
        axes1d[1] = axes[1]
        for E, T_M in self.parameter1.data :        # This logic ignores the interpolation of parameter1 as the only two subforms in ENDF/B-VII shows 
            parameters = [ EFL, EFH, T_M ]          # that linear-linear is better than the 'log-log' given in the ENDF/B-VII/data.
            g_Ep = XYs1d.createFromFunction( axes1d, xs, MadlandNixFunc, parameters, accuracy, biSectionMax = 12 )
            g_Ep.value = E                          # ????????? Class XYs1d does not the a proper setValue method. One should be added.
            g_Ep.normalize( insitu = True )
            pwl.append( g_Ep )
        return( pwl )

    @staticmethod
    def parseXMLNode( MNelement, xPath, linkData ):
        """Translate <MadlandNix> element from xml."""

        xPath.append( MNelement.tag )
        T_M_ = T_M.parseXMLNode( MNelement.find(T_M.moniker), xPath, linkData )
        EFL, EFH = [PQU.PQU(tmp) for tmp in (MNelement.get("EFL"),MNelement.get("EFH") )]
        MN = MadlandNix( EFL, EFH, T_M_ )
        xPath.pop()
        return MN

class NBodyPhaseSpace( subform ) :

    moniker = 'NBodyPhaseSpace'

    def __init__( self, numberOfProducts, numberOfProductsMasses ) :

        subform.__init__( self )
        self.numberOfProducts = numberOfProducts
        self.numberOfProductsMasses = numberOfProductsMasses

    def copy( self ):

        return NBodyPhaseSpace( self.numberOfProducts, self.numberOfProductsMasses )

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
        energyUnit = r.domainUnit( )
        EMin, EMax = r.domain( )

        rs = self.findClassInAncestry( reactionSuite.reactionSuite )
        projectileMass, targetMass = rs.projectile.getMass( massUnit ), rs.target.getMass( massUnit )

        c = self.findClassInAncestry( channels.channel )
        Q = c.Q.getConstantAs( energyUnit )

        axes = XYs2d.defaultAxes( standardsModule.frames.centerOfMassToken, energyUnit = energyUnit )
        pwl = XYs2d( axes )

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
            data = XYs1d( data, accuracy = accuracy, value = E_in )
            data.normalize( insitu = True )
            pwl.append( data )
        return( pwl )

    def toXMLList( self, indent = '', **kwargs ) :
        """Returns the xml string representation of self."""

        xmlString = [ self.XMLStartTagString( indent = indent, extraAttributesAsStrings = ' numberOfProducts="%s" mass="%s"' %
            ( self.numberOfProducts, self.numberOfProductsMasses ), emptyTag = True ) ]
        return( xmlString )

    @staticmethod
    def parseXMLNode( NBodyElement, xPath, linkData ):
        return NBodyPhaseSpace( int( NBodyElement.get( "numberOfProducts" ) ), PQU.PQU(NBodyElement.get("mass")) )

class weighted( ancestryModule.ancestry ) :

    moniker = 'weighted'

    def __init__( self, weight, functional ) :

        ancestryModule.ancestry.__init__( self )
        self.weight = weight
        self.weight.label = 'weight'
        self.functional = functional

    def copy( self ):

        return weighted( self.weight.copy(), self.functional.copy() )

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

    def __init__( self ) :

        subform.__init__( self )
        self.weights = []

    def __len__( self ) :

        return( len( self.weights ) )

    def __getitem__( self, i ) :

        return( self.weights[i] )

    def copy( self ):
        newWF = weightedFunctionals()
        for weight in self:
            newWF.addWeight( weight.copy() )
        return newWF

    __copy__ = __deepcopy__ = copy

    def addWeight( self, weight ) :

        if( not( isinstance( weight, weighted ) ) ) : raise Exception( 'Invalid weight of type %s' % brb.getType( weight ) )
        self.weights.append( weight )

    def check( self, info ) :

        from fudge.gnd import warning

        warnings = []
        totalWeight = sum( [w.weight for w in self.weights] )
        if (totalWeight.rangeMin() != 1.0) and (totalWeight.rangeMax() != 1.0):
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

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0 ) :

        if( len( self ) > 2 ) : raise Exception( 'more than two weights currently not supported' )
        E_ins, data = [], []
        for weighted_ in self :
            w = weighted_.weight.toPointwise_withLinearXYs( accuracy = accuracy, lowerEps = lowerEps, upperEps = upperEps )
            e = weighted_.functional.toPointwise_withLinearXYs( accuracy = accuracy, lowerEps = lowerEps, upperEps = upperEps )
            data.append( [ w, e ] )
            for x, y in w :
                if( x not in E_ins ) : E_ins.append( x )
            for x in e :
                if( x.value not in E_ins ) : E_ins.append( x.value )
        E_ins.sort( )
        pwl = XYs2d( axes = XYs2d.defaultAxes() )
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
            WF.weights.append( weighted( _weight, functional ) )
        xPath.pop( )
        return WF

class energyFunctionalDataToPointwise :

    def __init__( self, data, evaluateAtX ) :

        self.data = data
        self.evaluateAtX = evaluateAtX

    def toPointwise_withLinearXYs( self, accuracy = 0.001, lowerEps = 0, upperEps = 0 ) :

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

        parameter1 = self.data.parameter1.data.toPointwise_withLinearXYs( accuracy = None, lowerEps = lowerEps, upperEps = upperEps )
        axes = XYs2d.defaultAxes( )
        pwl = XYs2d( axes = axes )
        one_eV = PQU.PQU( '1 eV' ).getValueAs( axes[-1].unit )

        t = tester( accuracy, 1e-10, self.evaluateAtX )
        parameter2 = self.data.parameter2
        if( parameter2 is not None ) : parameter2 = parameter2.data.toPointwise_withLinearXYs( accuracy = None, lowerEps = lowerEps, upperEps = upperEps )
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

class energyLoss( subform, XYsModule.XYs1d ) :

    def __init__( self, **kwargs ) :

        XYsModule.XYs1d.__init__( self, **kwargs )
        subform.__init__( self )

    @staticmethod
    def defaultAxes( ) :

        axes = axesModule.axes( rank = 2 )
        axes[0] = axesModule.axis( 'energy_loss(energy_in)', 0, 'eV' )
        axes[1] = axesModule.axis( 'energy_in',  1, 'eV' )
        return( axes )

class form( baseModule.form ) :

    moniker = 'energy'
    subformAttributes = ( 'energySubform', )

    def __init__( self, label, productFrame, energySubform, makeCopy = True ) :

        if( not( isinstance( energySubform, subform ) ) ) : raise TypeError( 'instance is not an energy subform' )
        if( makeCopy ) : energySubform = energySubform.copy()
        baseModule.form.__init__( self, label, productFrame, ( energySubform, ) )

    def getSpectrumAtEnergy( self, energy ) :
        """Returns the energy spectrum for self at projectile energy."""

        return( self.evaluated.getSpectrumAtEnergy( energy ) )

    @staticmethod
    def parseXMLNode( energyElement, xPath, linkData ) :
        """Translate energy component from xml."""

        subformClass = {
                constant.moniker :                          constant,
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
