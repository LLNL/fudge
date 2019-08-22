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

import string
from fudge.core.ancestry import ancestry
from fudge.core.utilities import brb

__metaclass__ = type

linearToken = 'linear'
logToken = 'log'
flatToken = 'flat'
otherToken = 'other'
chargedParticleToken = 'charged-particle'

byRegionToken = 'byRegion'
unitBaseToken = 'unitBase'
correspondingPointsToken = 'correspondingPoints'

allowedIndependentInterpolations = [ linearToken, logToken, byRegionToken ]
allowedDependentInterpolations = allowedIndependentInterpolations + [ flatToken, chargedParticleToken ]
allowedInterpolationQualifiers = [ None, unitBaseToken, correspondingPointsToken ]

labToken = 'lab'
centerOfMassToken = 'centerOfMass'
allowedFrames = [ labToken, centerOfMassToken ]

allowedInterpolations = [ linearToken, logToken ]

class interpolationXY :

    def __init__( self, independent, dependent, qualifier = None ) :

        if( independent not in allowedIndependentInterpolations ) : raise Exception( 'Invalid independent interpolation = "%s"' % independent )
        if( dependent not in allowedDependentInterpolations ) : raise Exception( 'Invalid dependent interpolation = "%s": independent = %s' % ( dependent, independent ) )
        if( qualifier not in allowedInterpolationQualifiers ) : raise Exception( 'Invalid interpolation qualifier = "%s"' % qualifier )
        self.independent = independent
        self.dependent = dependent
        self.qualifier = qualifier

    def __eq__( self, other ) :

        if( self.independent != other.independent ) : return( False )
        if( self.dependent   != other.dependent ) :  return( False )
        if( self.qualifier   != other.qualifier ) :  return( False )
        return( True )

    def __ne__( self, other ) :

        return( not( self.__eq__( other ) ) )

    def __str__( self ) :

        if( self.independent is None ) : return( '' )
        q = ''
        if( self.qualifier is not None ) : q = '%s:' % self.qualifier
        return( '%s%s,%s' % ( q, self.independent, self.dependent ) )

    def copy( self ) :

        return( interpolationXY( self.independent, self.dependent, self.qualifier ) )

    def getInterpolationTokens( self ) :

        return( self.independent, self.dependent, self.qualifier )

    def isQualifierCorrespondingPoints( self ) :
        "Returns true if self's interpolation qualifier  is correspondingPoints."

        return( self.qualifier == correspondingPointsToken )

    def isQualifierUnitBase( self ) :
        "Returns true if self's interpolation qualifier  is unitbase."

        return( self.qualifier == unitBaseToken )

    @staticmethod
    def stringToTokens( s ) :

        qualifier = None
        independent, dependent = s.split( ',' )
        if( ':' in independent ) : qualifier, independent = independent.split( ':' )

        if(   independent == linearToken ) :
            independent = linearToken
        elif( independent == logToken ) :
            independent = logToken
        elif( independent == byRegionToken ) :
            independent = byRegionToken
        else :
            raise Exception( 'Invalid independent interpolation = "%s" from "%s"' % ( independent, s ) )

        if(   dependent == linearToken ) :
            dependent = linearToken
        elif( dependent == logToken ) :
            dependent = logToken
        elif( dependent == flatToken ) :
            dependent = flatToken
        elif( dependent == chargedParticleToken ) :
            dependent = chargedParticleToken
        elif( dependent == byRegionToken ) :
            dependent = byRegionToken
        else :
            raise Exception( 'Invalid dependent interpolation = "%s" from "%s"' % ( dependent, s ) )

        if(   qualifier == unitBaseToken ) :
            qualifier = unitBaseToken
        elif( qualifier == correspondingPointsToken ) :
            qualifier = correspondingPointsToken
        elif( qualifier is not None ) :
            raise Exception( 'Invalid qualifier interpolation = "%s" from "%s"' % ( qualifier, s ) )

        return( independent, dependent, qualifier )

    @staticmethod
    def stringToInterpolationXY( s ) :

        independent, dependent, qualifier = interpolationXY.stringToTokens( s )
        return( interpolationXY( independent, dependent, qualifier ) )

class axesBase( ancestry ) :
    """This is the base class for axes type classes."""

    def __init__( self, dimension, parent ) :

        self.dimension = dimension
        ancestry.__init__( self, self.moniker, parent )
        self.axes = []

    def __len__( self ) :

        return( self.dimension )

    def __getitem__( self, index ) :

        return( self.axes[index] )

    def __str__( self ) :

        l = [ str( axis ) for axis in self ]
        return( '\n'.join( l ) )

    def checkDimension( self, dimension ) :

        if( self.dimension != dimension ) : raise Exception( "self's dimension = %s != %s" % ( self.dimension, dimension ) )
        if( self.dimension != len( self.axes ) ) : raise Exception( 'Incomplete axes with dimension = %s and length = %s' % ( self.dimension, len( self.axes ) ) )

    def isLinear( self, qualifierOk = False, flatIsOk = False ) :
        """
        Returns True if all dependent axis of self are 'linear,linear' and False otherwise. In addition, if qualifierOk is False, then
        qualifier must also be None for returned value to be True; otherwise, qualifier is ignored.
        """

        for axis in self :
            if( not axis.isLinear( qualifierOk = qualifierOk, flatIsOk = flatIsOk ) ) : return( False )
        return( True )

    def toXML( self, indent = '' ) :

        return( '\n'.join( self.toXMLList( indent = indent ) ) )

class axis( ancestry ) :

    moniker = 'axis'

    def __init__( self, label, index, unit, interpolation = None, parent = None ) :
        "Returns a new instance of axis. Some checking of the input parameters is done."

        if( type( label ) != type( '' ) ) : raise Exception( 'label = "%s" is not a string' % label )
        if( type( index ) != type( 0 ) ) : raise Exception( 'index = "%s" is not an interger' % index )
        ancestry.__init__( self, self.moniker, parent, attribute = 'label' )

        self.setUnit( unit )
        self.label = label.strip( )
        self.index = index
        self.setInterpolation( interpolation )

    def __str__( self ) :

        interpolationStr = ''
        if( self.interpolation is not None ) : interpolationStr = ', interpolation="%s"' % self.interpolation
        return( 'label="%s", index="%s", unit="%s"%s' % ( self.label, self.index, self.unit, interpolationStr ) )

    def copy( self, index = None,  parent = None, standAlone = True ) :  # The standAlone agrument is not used, for compatibility with interpolationAxis.copy.
        "Returns a new instance that is a copy of self."

        if( index is None ) : index = self.index
        return( axis( self.label, index, self.unit, interpolation = self.interpolation, parent = parent ) )

    def getLabel( self ) :

        return( self.label )

    def getUnit( self ) :

        return( self.unit )

    def isLinear( self, qualifierOk = False, flatIsOk = False ) :

        if( self.interpolation is None ) : return( True )       # This is the dependent axis.
        independent, dependent, qualifier = self.interpolation.getInterpolationTokens( )
        if( flatIsOk ) : flatIsOk = ( independent == linearToken ) and ( dependent == flatToken )
        if( ( independent == dependent == linearToken ) or flatIsOk ) :
            if( qualifierOk or ( qualifier == None ) ) : return( True )
        return( False )

    def plotLabel( self, includeNonLinearInterpolation = True ) :

        label = self.label
        if( label == '' ) : label = 'unknown'
        if( self.unit != '' ) : label += ' (%s)' % self.unit
        if( includeNonLinearInterpolation and ( self.interpolation is not None ) ) :
            if( ( self.interpolation.independent != linearToken ) or ( self.interpolation.dependent != linearToken ) ) : label += ' (%s)' % self.interpolation
        return( label )

    def setInterpolation( self, interpolation ) :

        if( interpolation is not None ) :
            if( not( isinstance( interpolation, interpolationXY ) ) ) : raise Exception( 'Invalid interpolation: is instance of %s' % brb.getType( interpolation ) )
            interpolation = interpolation.copy( )
        self.interpolation = interpolation

    def setUnit( self, unit ) :
        "Sets self's unit. Only checks that unit is a string. If unit is None, it is set to an empty string (i.e., '')."

        if( unit is None ) : unit = ''
        if( type( unit ) != str ) : raise Exception( 'unit type "%s" is not a string' % type( unit ) )
        self.unit = unit.strip( )

    def toXML( self, indent = '' ) :

        interpolationStr = ''
        if( self.interpolation is not None ) : interpolationStr = ' interpolation="%s"' % self.interpolation
        return( '%s<axis index="%d" label="%s" unit="%s"%s/>' % ( indent, self.index, self.label, self.unit, interpolationStr ) )

    def unitConversionFactor( self, newUnit ) :
        "Returns as a float the factor needed to convert self's unit to newUnit. If units are not compatible, a raise is executed."

        from pqu import physicalQuantityWithUncertainty
        return( physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( 1., self.unit ).getValueAs( newUnit ) )

class axes( axesBase ) :

    moniker = 'axes'

    def __init__( self, dimension = 2, parent = None ) :

        axesBase.__init__( self, dimension, parent )

    def __setitem__( self, index_, axis_ ) :

        index = index_
        if( index < 0 ) : index += self.dimension
        if( index < 0 ) : raise IndexError( 'Index = %s' % index )
        if( index > len( self.axes ) ) : raise IndexError( 'index = %s is greater than length = %s' % ( index, len( self.axes ) ) )
        if( index >= self.dimension ) : raise IndexError( "index = %s is too big for self's dimension = %s" % ( index, self.dimension ) )
        if( not( isinstance( axis_, axis ) ) ) : raise Exception( 'axis_ is not an instance of axis: is an instance of %s', brb.getType( axis_ ) )
        if( index == ( self.dimension - 1 ) ) : 
            if( axis_.interpolation is not None ) : raise Exception( 'Last axis (i.e., dependent axis) cannot have interpolation' )
        else :
            if( axis_.interpolation is None ) : raise Exception( 'independent axis must have interpolation' )
        if( len( self.axes ) == index ) :
            self.axes.append( axis_.copy( index = index, parent = self ) )
        else :
            self.axes[index] = axis_.copy( index = index, parent = self )

    def copy( self, parent = None, standAlone = True ) :   # The standAlone and parent agruments are not used, for compatibility with interpolationAxes.copy.

        newAxisList = axes( dimension = self.dimension, parent = parent )
        for index, axis in enumerate( self ) : newAxisList[index] = axis        # __setitem__ makes a copy, so no need to here.
        return( newAxisList )

    def toXMLList( self, indent = '' ) :

        xmlString = [ '%s<axes>' % indent ]
        for axis in self : xmlString.append( axis.toXML( indent = indent + '  ' ) )
        xmlString[-1] += '</axes>'
        return( xmlString )

class referenceAxes( axesBase ) :

    moniker = 'referenceAxes'

    def __init__( self, parent, dimension = 2 ) :

        if( not( isinstance( parent.axes, axes ) or isinstance( parent.axes, referenceAxes ) ) ) : raise Exception( 'Parent does not have valid axes: parent type = %s' % brb.getType( parent ) )
        if( dimension > parent.axes.dimension ) : raise Exception( "dimension = %s greater than parent's dimension = %s" % ( dimension, parent.axes.dimension ) )
        axesBase.__init__( self, dimension, parent )

    def __getitem__( self, index ) :

        if( index < 0 ) : return( self.parent.axes[index] )
        return( self.parent.axes[self.parent.axes.dimension - self.dimension + index] )

    def checkDimension( self, dimension ) :

        if( self.dimension != dimension ) : raise Exception( "self's dimension = %s != %s" % ( self.dimension, dimension ) )

    def copy( self, parent = None, standAlone = False ) :   # The parent agruments is not used, for compatibility with interpolationAxes.copy.

        if( standAlone ) :
            axes_ = axes( dimension = self.dimension )
            for i in xrange( self.dimension ) : axes_[i]  = self[i]         # __setitem__ makes a copy, so no need to here.
            return( axes_ )
        return( referenceAxes( self.parent, dimension = self.dimension ) )

    def toXMLList( self, indent = '' ) :

        return( [] )

class interpolationAxis( ancestry ) :
    """This class is for internal use."""

    moniker = 'interpolationAxis'

    def __init__( self, index, interpolation, parent ) :

        ancestry.__init__( self, self.moniker, parent )
        self.index = index
        self.setInterpolation( interpolation )

    def __str__( self ) :

        return( 'index="%s" interpolation="%s"' % ( self.index, str( self.interpolation ) ) )

    def _getReferenceAxes( self ) :

        parent = self.parent
        while( parent is not None ) :
            if( hasattr( parent, 'axes' )  ) :
                if( isinstance( parent.axes, axes ) ) : break
            parent = parent.parent
        if( parent is None ) : raise Exception( 'No axes instance found in ancestry: %s' % self.toXLink( ) )
        return( parent.axes )

    def _getReferenceAxis( self ) :

        axes_ = self._getReferenceAxes( )
        return( axes_[self.index] )

    def copy( self, parent = None, standAlone = False ) :

        if( standAlone ) :
            index = 0
            if( self.index == -1 ) : index = 1
            axis_ = self._getReferenceAxis( ).copy( index = index )
            axis_.setInterpolation( self.interpolation )
        else :
            axis_ = interpolationAxis( self.index, self.interpolation, parent )
        return( axis_ )

    def getLabel( self ) :

        return( self._getReferenceAxis( ).getLabel( ) )

    def getUnit( self ) :

        return( self._getReferenceAxis( ).getUnit( ) )

    def isLinear( self, qualifierOk = False, flatIsOk = False ) :

        return( self._getReferenceAxis( ).isLinear( qualifierOk = qualifierOk, flatIsOk = flatIsOk ) )

    def plotLabel( self ) :

        return( self.copy( standAlone = True ).plotLabel( ) )

    def setInterpolation( self, interpolation ) :

        if( interpolation is not None ) :
            if( not( isinstance( interpolation, interpolationXY ) ) ) : raise Exception( 'Invalid interpolation: is instance of %s' % brb.getType( interpolation ) )
            interpolation = interpolation.copy( )
        self.interpolation = interpolation

    def unitConversionFactor( self, newUnit ) :
        "Returns as a float the factor needed to convert self's unit to newUnit. If units are not compatible, a raise is executed."

        return( self._getReferenceAxis( ).unitConversionFactor( newUnit ) )

class interpolationAxes( axesBase ) :
    """
    This axes type is designed for data with multiple regions and is used to describe the interpolation of a region.
    Each axis, an instance of interpolationAxis, only describes the interpolation of that axis and the index of the
    parent axes' axis. The axes type contains a link to the parent axes. The user should never instantiate an
    instance of interpolationAxis as the constructor for this class handles that task.
    """

    moniker = 'interpolationAxes'

    def __init__( self, independentIndex, interpolation, parent ) :

        axesBase.__init__( self, 2, parent )
        self[0] = interpolationAxis( independentIndex, interpolation, None )
        self[1] = interpolationAxis( -1, None, None )

    def __setitem__( self, index, axis_ ) :

        if( index < 0 ) : index += self.dimension
        if( index < 0 ) : raise Exception( 'Negative index = %s not allow' % index )
        if( index > len( self ) ) : raise Exception( 'index = %s is greater than length = %s' % ( index, len( self ) ) )
        if( index >= self.dimension ) : raise Exception( "index = %s is too big for self's dimension = %s" % ( index, self.dimension ) )
        if( not( isinstance( axis_, interpolationAxis ) ) ) : raise Exception( 'axis_ is not an instance of interpolationAxis' )
        if( index == len( self.axes ) ) :
            self.axes.append( axis_.copy( parent = self ) )
        else :
            self.axes[index] = axis_.copy( parent = self )

    def copy( self, parent = None, standAlone = False ) :

        if( standAlone and ( ( self.parent is not None ) or ( parent is not None ) ) ) :
            if( parent is None ) : parent = self.parent
            axesRef = self._getReferenceAxes( )
            axes_ = axes( 2, parent = parent )
            axes_[0] = axesRef[self[0].index]
            axes_[1] = axesRef[self[1].index]
            axes_[0].setInterpolation( self.axes[0].interpolation )
            return( axes_ )
        else :
            return( interpolationAxes( self.axes[0].index, self.axes[0].interpolation, parent ) )

    def _getReferenceAxes( self ) :

        parent = self.parent
        while( parent is not None ) :
            if( hasattr( parent, 'axes' )  ) :
                if( isinstance( parent.axes, axes ) ) : break
            parent = parent.parent
        if( parent is None ) : raise Exception( 'No axes instance found in ancestry: %s' % self.toXLink( ) )
        return( parent.axes )

    def toXMLList( self, indent = '' ) :

        return( [ '%s<%s index="%s" interpolation="%s"/>' % ( indent, self.moniker, self[0].index, self[0].interpolation ) ] )

def isValidAxes( axes_, doRaise = True ) :
    """Returns True if first argument is an instance one of the axes classes."""

    if( isinstance( axes_, axesBase ) ) : return( True )
    if( doRaise ) : raise Exception( 'axes_ argument is not a valid instance of an axes class: it is %s' % brb.getType( axes_ ) )
    return( False )

def defaultAxes( dimension = 2, independentInterpolation = linearToken, dependentInterpolation = linearToken, labelsUnits = {} ) :
    """
    Helper function to initialize ``axes``. Example, setting up ``axes`` for a pointwise lin-lin crossSection
    ::

        axes_ = defaultAxes( labelsUnits = { 0 : ( 'energy_in', 'eV' ), 1 : ( 'crossSection' , 'b' ) } ).
        
    """

    axes_ = axes( dimension = dimension )
    abcsOffset = string.ascii_lowercase.index( 'x' ) - ( dimension - 2 )
    interpolation = interpolationXY( independentInterpolation, dependentInterpolation )
    for i in xrange( dimension ) :
        if( i == ( dimension - 1 ) ) : interpolation = None
        label, unit = string.ascii_lowercase[abcsOffset+i], ''
        if( i in labelsUnits ) :
            label, unit = labelsUnits[i]
        axes_[i] = axis( label, i, unit, interpolation = interpolation )
    return( axes_ )

def parseXMLNode( axesElement, xPath=[] ):
    """Starting with '<axes>' or '<interpolationAxes>' element in xml-formatted GND document, read in axes information."""

    xPath.append( axesElement.tag )
    if( axesElement.tag == axes.moniker ) :
        ax = axes( dimension = len( axesElement ) )
        for axisElement in axesElement :
            interp = axisElement.get("interpolation")
            if interp is not None: interp = interpolationXY.stringToInterpolationXY( interp )
            ax.axes.append( axis( axisElement.get("label"), int(axisElement.get("index")),
                axisElement.get("unit"), interp, parent=ax) )
    elif( axesElement.tag == interpolationAxes.moniker ) :
        interp = interpolationXY.stringToInterpolationXY( axesElement.get("interpolation") )
        ax = interpolationAxes( int(axesElement.get("index")), interp, parent=None )
    else :
        raise Exception( 'Invalid axes type = %s' % ( axesElement.tag ) )
    xPath.pop()
    return ax

if( __name__ == "__main__" ) :

    class holder( ancestry ) :

        def __init__( self, axes ) :

            ancestry.__init__( self, 'holder', None )
            self.axes = axes
            axes.setParent( self )

    axes1 = axes( )
    axes1[0] = axis( 'energy_in', 0, 'eV', labToken, interpolation = interpolationXY( linearToken, linearToken ) )
    axes1[1] = axis( 'multiplicity', 1, ''  )
    print axes1
    print axes1.toXML( )

    print
#
# index does not matter as axes.__setitem__, copy axis an resets its index.
#
    axes2 = axes( dimension = 3 )
    axes2[0] = axis( 'energy_in',  0, 'eV',   labToken, interpolation = interpolationXY( linearToken, flatToken ) )
    axes2[1] = axis( 'energy_out', 0, 'eV',   centerOfMassToken, interpolation = interpolationXY( linearToken, logToken ) )
    axes2[2] = axis( "P(energy_in|energy_out)",    0, '1/eV', centerOfMassToken )
    print axes2
    print axes2.toXML( )
    print
    print 'Factor needed to convert "%s" to "MeV" is' % axes2[1].getUnit( ), axes2[1].unitConversionFactor( 'MeV' )
    print 'Factor needed to convert "%s" to "J" is' % axes2[1].getUnit( ), axes2[1].unitConversionFactor( 'J' )

    print
    axes3 = axes2.copy( )
    print axes3
    print
    axes3[2] = axes1[1].copy( )
    print axes3

    print axes3.toXML( )

    print '\nPlot labels:'
    for a in axes3 : print repr( a.plotLabel( ) )

    h1 = holder( axes2 )

    print
    iaxes1 = interpolationAxes( 1, interpolationXY( logToken, linearToken ), axes2 )
    print iaxes1

    print '\n%s\n' % iaxes1.toXML( )
    print '%s' % iaxes1[0].copy( standAlone = True )
    print '%s' % iaxes1[1].copy( standAlone = True )

    print '\n******* referenceAxes *******'
    print 'parent'
    print h1.axes
    print 'reference'
    ra1 = referenceAxes( h1 )
    print ra1
    print
    print ra1.copy( standAlone = True )
