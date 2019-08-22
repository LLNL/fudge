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

import abc
import datetime

from pqu import PQU as PQUModule

import xData.ancestry as ancestryModule
import xData.axes as axesModule

from fudge.gnds import physicalQuantity as physicalQuantityModule

from fudge.processing import flux as fluxModule
from fudge.processing import transportables as transportablesModule
from fudge.processing import inverseSpeed as inverseSpeedModule

__metaclass__ = type

def getClassDict():
    classDict = {}
    # Determine supported styles through module introspection
    import styles as tmpStyles
    for thing in filter(lambda x: not x.startswith('__'), dir(tmpStyles)):
        try:
            thingClass = tmpStyles.__dict__[thing]
            if issubclass(thingClass, style) and thingClass!=style:
                classDict[thingClass.moniker]=thingClass
        except TypeError: pass # catches imported modules that are irrelevant
    return classDict

class styles( ancestryModule.ancestry ) :
    """
    Stores the list of nuclear data styles that appear inside a file.
    The list generally includes one 'evaluated' style, plus 0 or more derived styles (heated, grouped, etc.)
    """

    moniker = 'styles'
    ancestryMembers = ( '', )

    def __init__( self ) :

        self.__styles = []
        ancestryModule.ancestry.__init__( self )

    def __contains__( self, label ) :

        for item in self :
            if( item.label == label ) : return( True )
        return( False )

    def __len__( self ) :

        return( len( self.__styles ) )

    def __iter__( self ) :

        n1 = len( self )
        for i1 in range( n1 ) : yield self.__styles[i1]

    def __getitem__( self, label ) :

        if( isinstance( label, int ) ) : return( self.__styles[label] )
        if( not( isinstance( label, str ) ) ) : raise TypeError( 'label must be a string' )
        for _style in self.__styles :
            if( _style.label == label ) : return( _style )
        raise IndexError( 'No style labelled == "%s"' % label )

    def add( self, _style ) :
        """
        Append a style to the list of styles.

        @param _style: style instance
        """

        if( not( isinstance( _style, style ) ) ) : raise TypeError( 'invalid style instance' )
        if( isinstance( _style, evaluated ) ) :
            for __style in self :
                if( isinstance( __style, evaluated ) ) : raise ValueError( 'multiple %s styles not allowed' % __style.moniker )
        for __style in self :
            if( __style.label == _style.label ) : raise ValueError( 'style labeled "%s" already exists' % _style.label )
        self.__styles.append( _style )
        _style.setAncestor( self, 'label' )

    def convertUnits( self, unitMap ) :

        for _style in self : _style.convertUnits( unitMap )

    def remove( self, _style ):
        """
        Remove the specified style. Accepts either a string or a style instance
        """

        index = None
        for idx,tmp in enumerate(self):
            if( ( tmp is _style ) or ( tmp.label == _style ) ) : index = idx
        if index is None:
            raise KeyError("style '%s' not found in styles" % _style)
        self.__styles.pop(index)

    def getStylesOfClass( self, cls ) :
        """
        Returns a list of all styles of class cls.
        """

        styleList = []
        for _style in self :
            if( isinstance( _style, cls ) ) : styleList.append( _style )
        return( styleList )

    def getStyleOfClass( self, cls ) :
        """
        Returns a list of all styles of class cls.
        """

        styleList = self.getStylesOfClass( cls )
        if( len( styleList ) == 0 ) : return( None )
        if( len( styleList )  > 1 ) : raise KeyError( '''multiple (%d) style's of type "%s" found''' % \
                ( len( styleList ), cls.moniker ) )
        return( styleList[0] )

    def getTempStyleNameOfClass( self, cls ):

        styleNames = [ _style.label for _style in self ]
        index = 0
        while( True ) :
            tmpLabel = 'tmp%d_%s' % ( index, cls.moniker )
            if( tmpLabel not in styleNames ) : return( tmpLabel )
            index += 1

    def getEvaluatedStyle( self ) :

        return( self.getStyleOfClass( evaluated ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = ['%s<%s>' % (indent, self.moniker)]
        for _style in self : xmlStringList += _style.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    def parseXMLNode( self, stylesElement, xPath, linkData ) :

        xPath.append( stylesElement.tag )

        classDict = {}
        for _style in ( evaluated, crossSectionReconstructed, angularDistributionReconstructed, averageProductData, MonteCarlo_cdf, multiGroup, 
                heated, heatedMultiGroup, griddedCrossSection ) :
            classDict[_style.moniker] = _style
        for styleElement in stylesElement :
            _class = classDict.get( styleElement.tag, None )
            if( _class is None ) :
                raise TypeError( 'encountered unknown style "%s"' % styleElement.tag )
            self.add( _class.parseXMLNode( styleElement, xPath, linkData ) )

        xPath.pop()

class style( ancestryModule.ancestry ) :

    __metaclass__ = abc.ABCMeta

    def __init__( self, label, derivedFrom, date = None ) :

        ancestryModule.ancestry.__init__( self )

        if( not( isinstance( label, str ) ) ) : raise TypeError( 'label must be a str instance.' )
        self.__label = label

        if( date is None ) : date = str( datetime.date.today( ) )
        self.__date = date

        if( not( isinstance( derivedFrom, str ) ) ) : raise TypeError( 'label must be a str instance.' )
        self.__derivedFrom = derivedFrom

    @property
    def date( self ) :

        return( self.__date )

    @property
    def derivedFrom( self ) :

        return( self.__derivedFrom )

    @property
    def derivedFromStyle( self ) :

        for _style in self.ancestor :
            if( _style.label == self.__derivedFrom ) : return( _style )
        return( None )

    @property
    def label( self ) :

        return( self.__label )

    @property
    def temperature( self ) :

        return( self.derivedFromStyle.temperature )

    def findFormMatchingDerivedStyle( self, component, styleFilter = None ) :
        """
        This method searches the link of derivedFroms, starting with self's derivedFrom,
        to find a form in component matching one of the derivedFroms. If a form is found
        matching one of the derivedFroms, that form is returned. If no match is found, None is returned.
        """

        def alwaysTrue( a_style ) : return( True )

        if( styleFilter is None ) : styleFilter = alwaysTrue

        parent = self.ancestor
        derivedFrom = self.derivedFrom
        while( derivedFrom != '' ) :
            for form in component :
                if( ( derivedFrom == form.label ) and styleFilter( parent[derivedFrom] ) ) : return( form )
            derivedFrom = parent[derivedFrom].derivedFrom
        return( None )
        
    def findDerivedFromStyle( self, cls ) :

        _style = self
        while( _style is not None ) :
            _style = _style.derivedFromStyle
            if( ( _style is not None ) and ( isinstance( _style, cls ) ) ) : return( _style )
        return( None )

    def sibling( self, label ) :
        """Returns the sibling of self's with label."""

        return( self.ancestor[label] )

    def convertUnits( self, unitMap ) :

        pass

    def XMLCommonAttributes( self ) :

        XMLCommon = 'label="%s"' % self.label
        if( self.derivedFrom != "" ) : XMLCommon += ' derivedFrom="%s"' % self.derivedFrom
        if( self.date is not None ) : XMLCommon += ' date="%s"' % self.date
        return( XMLCommon )

    @staticmethod
    def parseXMLNodeBase( element, xPath ) :

        label = element.get( 'label' )
        xPath.append( '%s[@label="%s"]' % ( element.tag, label ) )
        derivedFrom = element.get( 'derivedFrom', "" )
        date = element.get( 'date', None )
        return( label, derivedFrom, date )

class evaluated( style ) :

    moniker = 'evaluated'

    def __init__( self, label, derivedFrom, temperature, projectileEnergyDomain, library, version, date = None ) :

        style.__init__( self, label, derivedFrom, date = date )

        if( not( isinstance( library, str ) ) ) : raise TypeError( 'library must be a string' )
        self.__library = library

        if( not( isinstance( version, str ) ) ) : raise TypeError( 'version must be a string' )
        self.__version = version

        if( not( isinstance( temperature, physicalQuantityModule.temperature ) ) ) :
            raise TypeError( 'invalid temperature object' )
        self.__temperature = temperature

        self.projectileEnergyDomain = projectileEnergyDomain

    @property
    def library( self ) :

        return( self.__library )

    @library.setter
    def library( self, value ) :

        self.__library = value

    @property
    def projectileEnergyDomain(self): return self.__projectileEnergyDomain

    @projectileEnergyDomain.setter
    def projectileEnergyDomain(self, value):
        if not isinstance(value, projectileEnergyDomain):
            raise TypeError("Expected projectileEnergyDomain instance, got '%s' instead" % type(value))
        value.setAncestor(self)
        self.__projectileEnergyDomain = value

    @property
    def temperature( self ) :

        return( self.__temperature )

    @property
    def version( self ) :

        return( self.__version )

    @version.setter
    def version( self, value ) :

        self.__version = value

    def convertUnits( self, unitMap ) :

        self.__temperature.convertUnits( unitMap )
        self.__projectileEnergyDomain.convertUnits( unitMap )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s %s library="%s" version="%s">' %
                ( indent, self.moniker, self.XMLCommonAttributes( ), self.library, self.version ) ]
        xmlStringList += self.temperature.toXMLList( indent2, **kwargs )
        xmlStringList += self.projectileEnergyDomain.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        label, derivedFrom, date = style.parseXMLNodeBase( element, xPath )

        library = element.get( 'library' )
        version = element.get( 'version' )
        temperature = physicalQuantityModule.temperature.parseXMLNode( element.find( 'temperature' ), xPath, linkData )
        projectileEnergyDomain_ = projectileEnergyDomain.parseXMLNode(
            element.find( projectileEnergyDomain.moniker ), xPath, linkData )

        _evaluated = evaluated( label, derivedFrom, temperature, projectileEnergyDomain_, library, version, date = date )

        xPath.pop( )
        return( _evaluated )

class projectileEnergyDomain( ancestryModule.ancestry ) :

    moniker = 'projectileEnergyDomain'

    def __init__(self, min, max, unit):

        ancestryModule.ancestry.__init__(self)
        self.min = min
        self.max = max
        self.unit = unit

    def convertUnits( self, unitMap ) :

        if( self.unit in unitMap ) :
            factor = PQUModule.PQU( 1, self.unit ).getValueAs( unitMap[self.unit] )
            self.min *= factor
            self.max *= factor
            self.unit = unitMap[self.unit]

    def toXMLList( self, indent = '', **kwargs ) :

        return ['%s<%s min="%r" max="%r" unit="%s"/>'
                % (indent, self.moniker, self.min, self.max, self.unit)]

    @classmethod
    def parseXMLNode(cls, element, xPath, linkData):

        xPath.append(element.tag)
        PED = cls( float(element.get('min')), float(element.get('max')), element.get('unit') )
        xPath.pop()
        return PED

class crossSectionReconstructed( style ) :

    moniker = 'crossSectionReconstructed'

    def __init__( self, label, derivedFrom, date = None ) :

        style.__init__( self, label, derivedFrom, date = date )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s %s>' % ( indent, self.moniker, self.XMLCommonAttributes( ) ) ]
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        label, derivedFrom, date = style.parseXMLNodeBase( element, xPath )

        _crossSectionReconstructed = crossSectionReconstructed( label, derivedFrom, date = date )

        xPath.pop()
        return( _crossSectionReconstructed )

class angularDistributionReconstructed( style ) :

    moniker = 'angularDistributionReconstructed'

    def __init__( self, label, derivedFrom, date = None ) :

        style.__init__( self, label, derivedFrom, date = date )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s %s>' % ( indent, self.moniker, self.XMLCommonAttributes( ) ) ]
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        label, derivedFrom, date = style.parseXMLNodeBase( element, xPath )

        _angularDistributionReconstructed = angularDistributionReconstructed( label, derivedFrom, date = date )

        xPath.pop()
        return( _angularDistributionReconstructed )

class CoulombPlusNuclearElasticMuCutoff( style ) :

    moniker = 'CoulombPlusNuclearElasticMuCutoff'

    def __init__( self, label, derivedFrom, muCutoff, date = None ) :

        style.__init__( self, label, derivedFrom, date = date )

        if( not( isinstance( muCutoff, float ) ) ) : raise TypeError( 'invalid muCutoff object' )
        self.__muCutoff = muCutoff

    @property
    def muCutoff( self ) :

        return( self.__muCutoff )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s muCutoff="%s" %s>' % ( indent, self.moniker, self.__muCutoff, self.XMLCommonAttributes( ) ) ]
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        label, derivedFrom, date = style.parseXMLNodeBase( element, xPath )

        muCutoff = float( element.find( 'muCutoff' ), xPath, linkData )
        _CoulombPlusNuclearElasticMuCutoff = CoulombPlusNuclearElasticMuCutoff( label, derivedFrom, muCutoff, date = date )

        xPath.pop( )
        return( _CoulombPlusNuclearElasticMuCutoff )

class heated( style ) :

    moniker = 'heated'

    def __init__( self, label, derivedFrom, temperature, date = None ) :

        style.__init__( self, label, derivedFrom, date = date )

        if( not( isinstance( temperature, physicalQuantityModule.temperature ) ) ) :
            raise TypeError( 'invalid temperature object' )
        self.__temperature = temperature

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s %s>' % ( indent, self.moniker, self.XMLCommonAttributes( ) ) ]
        xmlStringList += self.temperature.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @property
    def temperature( self ) :

        return( self.__temperature )

    def convertUnits( self, unitMap ) :

        self.__temperature.convertUnits( unitMap )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        label, derivedFrom, date = style.parseXMLNodeBase( element, xPath )

        temperature = physicalQuantityModule.temperature.parseXMLNode( element.find( 'temperature' ), xPath, linkData )
        _heated = heated( label, derivedFrom, temperature, date = date )

        xPath.pop()
        return( _heated )

class averageProductData( style ) :

    moniker = 'averageProductData'

    def __init__( self, label, derivedFrom, date = None ) :

        style.__init__( self, label, derivedFrom, date = date )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s %s>' % ( indent, self.moniker, self.XMLCommonAttributes( ) ) ]
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        label, derivedFrom, date = style.parseXMLNodeBase( element, xPath )

        _averageProductData = averageProductData( label, derivedFrom, date = date )

        xPath.pop()
        return( _averageProductData )

class multiGroup( style ) :

    moniker = 'multiGroup'

    def __init__( self, label, lMax, date = None ) :

        style.__init__( self, label, "", date = date )
        self.__lMax = int( lMax )

        self.__transportables = transportablesModule.transportables( )
        self.__transportables.setAncestor( self )

    @property
    def lMax( self ) :

        return( self.__lMax )

    @property
    def transportables( self ) :

        return( self.__transportables )

    def convertUnits( self, unitMap ) :
        # BRB - FIXME Still need to fully implement.

        pass

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s %s lMax="%s">' % ( indent, self.moniker, self.XMLCommonAttributes( ), self.lMax ) ]
        xmlStringList += self.transportables.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        label, derivedFrom, date = style.parseXMLNodeBase( element, xPath )

        lMax = element.get( 'lMax' )

        _multiGroup = multiGroup( label, lMax, date = date )
        for child in element:
            if child.tag == transportablesModule.transportables.moniker:
                _multiGroup.transportables.parseXMLNode( child, xPath, linkData )
            else:
                raise TypeError("Encountered unexpected element '%s'" % child.tag)

        xPath.pop( )
        return( _multiGroup )

class heatedMultiGroup( style ) :

    moniker = 'heatedMultiGroup'

    def __init__( self, label, derivedFrom, parameters, flux, date = None ) :

        style.__init__( self, label, derivedFrom, date = date )

        if( not( isinstance( parameters, str ) ) ) : raise TypeError( 'parameters must be a str instance' )
        self.__parameters = parameters

        if( not( isinstance( flux, fluxModule.flux ) ) ) : raise TypeError( 'Invalid flux instance' )
        self.__flux = flux
        self.__flux.setAncestor( self )

        self.__inverseSpeed = None

        self.__multiGroupFlux = None                # Only for internal use, it not written to GNDS file.

    @property
    def flux( self ) :

        return( self.__flux )

    @property
    def inverseSpeed( self ) :

        return( self.__inverseSpeed )

    @inverseSpeed.setter
    def inverseSpeed( self, value ) :

        self.__inverseSpeed = value

    @property
    def lMax( self ) :

        _style = self.sibling( self.parameters )
        if( _style is None ) : raise TypeError( 'Cannot find a multiGroup in styles' )
        return( _style.lMax )

    @property
    def multiGroupFlux( self ) :

        return( self.__multiGroupFlux )

    @property
    def parameters( self ) :

        return( self.__parameters )

    @property
    def transportables( self ) :
 
        _style = self.sibling( self.parameters )
        if( _style is None ) : raise TypeError( 'Cannot find a multiGroup in styles' )
        return( _style.transportables )

    def convertUnits( self, unitMap ) :
        # BRB - FIXME Still need to fully implement.

        self.__inverseSpeed.convertUnits( unitMap )

    def processMultiGroup( self, style, tempInfo, indent ) :

        self.__multiGroupFlux = self.flux.processMultiGroup( self, tempInfo, indent )
        self.inverseSpeed = inverseSpeedModule.inverseSpeed( inverseSpeedModule.multiGroupInverseSpeed( self, tempInfo ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s %s parameters="%s">' % ( indent, self.moniker, self.XMLCommonAttributes( ), self.parameters ) ]
        xmlStringList += self.flux.toXMLList( indent2, **kwargs )
        if( self.inverseSpeed is not None ) : xmlStringList += self.inverseSpeed.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        label, derivedFrom, date = style.parseXMLNodeBase( element, xPath )
        parameters = element.get('parameters')

        flux = None
        inverseSpeed = None
        for child in element :
            if( child.tag == fluxModule.flux.moniker ) :
                flux = fluxModule.flux.parseXMLNode( child, xPath, linkData )
            elif( child.tag == inverseSpeedModule.inverseSpeed.moniker ) :
                inverseSpeed = inverseSpeedModule.inverseSpeed.parseXMLNode( child, xPath, linkData )

        _heatedMultiGroup = heatedMultiGroup( label, derivedFrom, parameters, flux, date = date )
        _heatedMultiGroup.inverseSpeed = inverseSpeed

        xPath.pop( )
        return( _heatedMultiGroup )

class MonteCarlo_cdf( style ) :

    moniker = 'MonteCarlo_cdf'

    def __init__( self, label, derivedFrom, date = None ) :

        style.__init__( self, label, derivedFrom, date = date )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s %s>' % ( indent, self.moniker, self.XMLCommonAttributes( ) ) ]
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        label, derivedFrom, date = style.parseXMLNodeBase( element, xPath )

        _MonteCarlo_cdf = MonteCarlo_cdf( label, derivedFrom, date = date )

        xPath.pop( )
        return( _MonteCarlo_cdf )

class griddedCrossSection( style ) :

    moniker = 'griddedCrossSection'

    def __init__( self, label, derivedFrom, date = None ) :

        style.__init__( self, label, derivedFrom, date = date )

        self.__grid = None

    @property
    def grid( self ) :

        return( self.__grid )

    @grid.setter
    def grid( self, grid ) :

        if( not( isinstance( grid, axesModule.grid ) ) ) : raise TypeError( 'grid not instance of axesModule.grid' )
        self.__grid = grid
        self.__grid.setAncestor( self )

    def convertUnits( self, unitMap ) :

        self.__grid.convertUnits( unitMap )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s %s>' % ( indent, self.moniker, self.XMLCommonAttributes( ) ) ]
        xmlStringList += self.__grid.toXMLList( indent = indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        label, derivedFrom, date = style.parseXMLNodeBase( element, xPath )

        grid = None
        for child in element :
            if( child.tag == axesModule.grid.moniker ) :
                grid = axesModule.grid.parseXMLNode( child, xPath, linkData )

        _griddedCrossSection = griddedCrossSection( label, derivedFrom, date = date )
        _griddedCrossSection.grid = grid

        xPath.pop( )
        return( _griddedCrossSection )

def findEvaluated( forms ) :

    return( findStyle( evaluated, forms ) )

def findAllOfStyle( cls, forms ) :

    formsFound = []
    for form in forms :
        if( isinstance( form, cls ) ) : formsFound.append( form )
    return( formsFound )

def findStyle( cls, forms ) :

    formsFound = findAllOfStyle( cls, forms )
    if( len( formsFound )  > 1 ) : raise Exception( 'multiple styles found' )
    if( len( formsFound ) == 0 ) : raise Exception( 'no style found' )
    return( formsFound[0] )
