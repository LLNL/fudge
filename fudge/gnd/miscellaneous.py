# <<BEGIN-copyright>>
# <<END-copyright>>

import sys
import copy
import math

from fudge.core.ancestry import ancestry
from fudge.core.math import *
from fudge.legacy.converting import endfFormats
from fudge.core.math.xData import axes

__metaclass__ = type

class style :

    def __init__( self, name, attributes = {}, subStyle = None ) :

        self.name = name
        self.attributes = {}
        for attribute in attributes : self.attributes[attribute] = attributes[attribute]
        self.subStyle = subStyle

    def toXMLList( self, indent = '' ) :

        attributeString = ""
        for attribute in self.attributes : attributeString += ' %s="%s"' % ( attribute, self.attributes[attribute] )
        xmlString = [ '%s<style name="%s"%s>' % ( indent, self.name, attributeString ) ]
        if( not( self.subStyle is None ) ) : xmlString += self.subStyle.toXMLList( indent = indent + '  ' )
        xmlString[-1] += '</style>'
        return( xmlString )

    def getAttribute( self, name ) :
        """Returns value for attribute name if it exists; otherwise, returns None."""

        if( name in self.attributes ) : return( self.attributes[name] )
        return( None )

    def setAttribute( self, name, value ) :
        """Adds attribute name and its value to the list of attributes."""

        self.attributes[name] = value

    @staticmethod
    def parseXMLNode( styleElement, linkData={} ):
        attrs = dict( styleElement.items() )
        return style( attrs.pop("name"), attributes=attrs )

class subStyleLLNL_Pn :

    def __init__( self, particles ) :

        self.particles = particles

    def toXMLList( self, indent = '' ) :

        xmlString = []
        for particle in self.particles :
            particleInfo = self.particles[particle]
            xmlString.append( '%s<particle name="%s" group="%s" lMax="%d" conservationFlag="%s"/>' % \
                ( indent, particleInfo.name, particleInfo.groups.name, particleInfo.lMax, particleInfo.conservationFlag ) )
        return( xmlString )

def find_gndPath( reactionSuite, gndPath ) :

    def find_gndPath2( parent, paths, gndPath ) :

        if( len( paths ) == 0 ) : return( parent )
        current, others = paths[0], paths[1:]
        if( current[:6] == 'label:' ) :
            label = current[6:]
            for element in parent :
                if( element.getLabel( ) == label ) : return( find_gndPath2( element, others, gndPath ) )
            raise Exception( 'Could not find path fragment = %s in gndPath = %s' % ( '/'.join( paths ), gndPath ) )
        else :
            raise Exception( 'Finding element currently not supported for gndPath = %s' % gndPath )

    return( find_gndPath2( reactionSuite, gndPath.split( '/' )[1:], gndPath ) )

class undefinedLevel :

    def __init__( self, level ) :

        self.level = level

    def __str__( self ) :

        return( 'u:%s' % str( self.level ) )

    def getValueAs( self, unit ) :

        return( self.level.getValueAs( unit ) )

def floatToString( f, precision = 12 ) :
    """Returns the shortest string for the float from %f and %e conversion at precision."""

    s = ( '%%.%df' % precision ) % f
    s = s.rstrip( '0' )
    t = ( '%%.%de' % precision )  % f
    t2 = t.split( 'e' )
    p = t2[0].rstrip( '0' )
    t3 = p + 'e' + t2[1]
    e = int( t2[1] )
    if( ( e < 0 ) and ( ( len( p ) - 2 - e ) > precision ) ) : s = t3
    if( len( t3 ) <= len( s ) ) : s = t3
    return( s )

def copyVariables( variables, indices = [] ) :
    """Need to modify to handle axes instance for variables.???????"""

    if( indices == [] ) : indices = variables.keys( )
    variables_ = {}
    for index in indices :
        if( type( index ) == type( 1 ) ) :
            iFrom = iTo = index
        else :                              # Must be a list of length 2
            if( len( index ) != 2 ) : raise Exception( 'Index = "%s" not a list of integers of length 2' % str( index ) )
            iFrom, iTo = index
        variables_[iTo] = copy.copy( variables[iFrom] )
    return( variables_ )

def TMs2Form( processInfo, tempInfo, newComponents, TM_1, TM_E, axes_p ) :

    from productData import distributions
    from fudge.core.math import matrix
    crossSectionUnit = 'b'                          # ?????? 'b' should not be hardwired.
    axes_ = axes.axes( dimension = 3 )
    axes_[0] = axes.axis( axes_p[0].getLabel( ), 0, axes_p[0].getUnit( ), frame = axes.labToken, interpolation = axes.interpolationXY( axes.linearToken, axes.flatToken ) )
    axes_[1] = axes.axis( 'energy_out',          1, axes_p[0].getUnit( ), frame = axes.labToken, interpolation = axes.interpolationXY( axes.linearToken, axes.flatToken ) )
    axes_[2] = axes.axis( 'C_l(energy_in,energy_out)', 2, crossSectionUnit, frame = axes.labToken )
    # convert TM_1 and TM_E from dicts into matrix objects:
    if( not ( TM_1 is None ) ) :
        TM_1_new = []
        for i in range(len(TM_1)):
            TM_1_new.append( matrix.matrix( [TM_1[i][idx] for idx in range(len(TM_1[i]))], form="sparse_asymmetric" ) )
        component = distributions.Legendre.component( nativeData = distributions.base.groupedFormToken )
        component.addForm( distributions.Legendre.grouped( axes_, TM_1_new ) )
        newComponents.append( component )
    if( not ( TM_E is None ) ) :
        TM_E_new = []
        for i in range(len(TM_E)):
            TM_E_new.append( matrix.matrix( [TM_E[i][idx] for idx in range(len(TM_E[i]))], form="sparse_asymmetric" ) )
        component = distributions.Legendre.energyConservationComponent( nativeData = distributions.base.groupedFormToken )
        component.addForm( distributions.Legendre.energyConservationGrouped( axes_, TM_E_new ) )
        newComponents.append( component )

def getMassLevel( particle, massUnit, levelUnit ) :

    return( particle.getMass( massUnit ), particle.getLevelAsFloat( levelUnit ) )

def makeGrouped( self, processInfo, tempInfo, data, normType = 'groupedCrossSectionNorm' ) :

    from fudge.processing import miscellaneous

    projectile, target = processInfo.getProjectileName( ), processInfo.getTargetName( )
    norm = tempInfo[normType]
    crossSection = tempInfo['crossSection']
    label, unit = data.axes[1].getLabel( ), data.axes[1].getUnit( )
    axes_ = crossSection.axes.copy( )
    if( normType == 'groupedFlux' ) :
        axes_[1].label = "%s %s" % ( label,  axes_[1].label )
        axes_[1].unit = "%s %s" % ( unit,  axes_[1].getUnit( ) )
    else :
        axes_[1].label = "%s" % label
        axes_[1].unit = "%s" % unit
    independent, dependent, qualifier = axes_[0].interpolation.getInterpolationTokens( )
    axes_[0].interpolation = axes.interpolationXY( independent, axes.flatToken )
    grouped = miscellaneous.groupTwoFunctionsAndFlux( projectile, processInfo, crossSection, data, norm = norm )
    return( axes_, grouped )
