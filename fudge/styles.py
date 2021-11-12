# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import abc
import datetime

from pqu import PQU as PQUModule

from xData import formatVersion as formatVersionModule
from xData import ancestry as ancestryModule
from xData import axes as axesModule
from xData.Documentation import documentation as documentationModule

from fudge import physicalQuantity as physicalQuantityModule

from fudge.processing import flux as fluxModule
from fudge.processing import transportables as transportablesModule
from fudge.processing import inverseSpeed as inverseSpeedModule

__metaclass__ = type

def getClassDict():
    """Determine supported styles through module introspection"""

    classDict = {}
    from . import styles as tmpStyles
    for thing in [ x for x in dir( tmpStyles ) if not x.startswith('__')  ] :
        try:
            thingClass = tmpStyles.__dict__[thing]
            if issubclass(thingClass, style) and thingClass!=style:
                classDict[thingClass.moniker]=thingClass
        except TypeError: pass # catches imported modules that are irrelevant
    return classDict

class styles( ancestryModule.ancestry ) :
    """
    Stores the list of nuclear data styles that appear inside a file.

    The list generally includes one 'evaluated' style, plus 0 or more derived styles (heated, heatedMultiGroup, etc.).
    """

    moniker = 'styles'
    """Label for use as the tag or identifier for this style"""

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
        for __style in self :
            if( __style.label == _style.label ) : raise ValueError( 'style labeled "%s" already exists' % _style.label )

        self.__styles.append( _style )
        _style.setAncestor( self, 'label' )

    def chains( self, ends = False, _styles = None ) :
        """For self, determine the derived from chains for each style. If ends is True, trains which are a part of a longer train are removed from the list."""

        if( _styles is None ) :
            _chains = [ _style.chain( ) for _style in self ]
        else :
            _chains = [ _style.chain( ) for _style in _styles ]

        if( ends ) :
            removeIndices = set( )
            for i1 in range( len( _chains ) ) :
                style1 = _chains[i1][0]
                for i2 in range( len( _chains ) ) :
                    if( style1 in _chains[i2][1:] ) : removeIndices.add( i1 )
            for index in reversed( sorted( removeIndices ) ) : _chains.pop( index )

        return( _chains )

    def chainWithLabel( self, label ) :

        chains = self.chains( ends = True )
        for chain in chains :
            for style in chain :
                if( style.label == label ) : return( chain )
        raise ValueError( 'No style with label = "%s" found.' % label )

    def preProcessingHeadInChainWithLabel( self, label ) :
        """
        Returns the pre-processing head of the chain containing the label *label*.
        """

        chain = self.chainWithLabel( label )
        preProcessingStyles = self.preProcessingStyles( )
        for style in chain :
            if( style.__class__ in preProcessingStyles ) : return( style )
        raise ValueError( 'No pre-processing style in chain with label = "%s" found.' % label )

    def preProcessingOnly( self ) :
        """
        Returns True if self only contains pre-processed styles and False otherwise.
        """

        preProcessingStyles = self.preProcessingStyles( )
        clean = True                            # True means no processed data.
        for style in self :
            if( not( isinstance( style, preProcessingStyles ) ) ) :
                clean = False
                break

        return( clean )

    def convertUnits( self, unitMap ) :

        for _style in self : _style.convertUnits( unitMap )

    def preProcessingChainHead( self, label = None ) :
        """
        Returns the pre-processing style that heads the chain containing the pre-processing label *label*. If label is None, there must be
        only one pre-processing chain and that head is returned.
        """

        chains = self.preProcessingChains( ends = True )
        if( len( chains ) == 0 ) : raise Exception( 'No pre-processing style exist.' )
        if( label is None ) :
            if( len( chains ) > 1 ) : raise Exception( 'Must have only one pre-processing chain, have %s.' % len( preProcessingChains ) )
            return( chains[0][0] )
        else :
            for chain in chains :
                for style in chain :
                    if( style.label == label ) : return( chain[0] )
            raise Exception( 'No pre-processing chain containing label "%s" found.' % label )

    def preProcessingChains( self, ends = True ) :
        """
        Returns chains for each style that can be processed (i.e., evaluated, crossSectionReconstructed, 
        angularDistributionReconstructed and Realization).
        """

        preProcessingStyles = self.preProcessingStyles( )
        preProcessingStyleInstances = []
        for style in self :
            if( isinstance( style, preProcessingStyles ) ) : preProcessingStyleInstances.append( style )

        return( self.chains( ends = ends, _styles = preProcessingStyleInstances ) )

    def preProcessingStyles( self ) :
        """Returns the list of style classes that can be processed. These are known as pre-processed styles."""

        return( evaluated, crossSectionReconstructed, angularDistributionReconstructed, Realization )

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
        Returns the style of class cls. If no style is found, None is returned. If more than 1 style is found, an exception is raised.
        """

        styleList = self.getStylesOfClass( cls )
        if( len( styleList ) == 0 ) : return( None )
        if( len( styleList )  > 1 ) : raise KeyError( '''multiple (%d) style's of type "%s" found''' % ( len( styleList ), cls.moniker ) )
        return( styleList[0] )

    def getTempStyleNameOfClass( self, cls ):

        styleNames = [ _style.label for _style in self ]
        index = 0
        while( True ) :
            tmpLabel = 'tmp%d_%s' % ( index, cls.moniker )
            if( tmpLabel not in styleNames ) : return( tmpLabel )
            index += 1

    def getEvaluatedStyle( self ) :

        _evaluateds = self.getStylesOfClass( evaluated )
        _styles = self.chains( ends = True, _styles = _evaluateds )
        if( len( _styles ) == 0 ) : raise Exception( '''No evaluated style found.''' )
        if( len( _styles ) > 1 ) : raise Exception( '''Multiple (%d) evaluated styles.''' % ( len( styleList ) ) )

        return( self[_styles[0][0].label] )

    def removeStyles( self, styleLabels ) :
        """
        Removes all forms whose label matches one of the labels in removeStyles.
        """

        for label in styleLabels :
            for index, style in enumerate( self.__styles ) :
                if( style.label == label ) :
                    del self.__styles[index]
                    break

    def temperatures( self, unit = 'MeV/k' ) :
        """
        Returns a listed of [ temperature, dictionary ] instances sorted by temperature. The directory keys are the
        monikers for the heated process styles. The associated values are the labels for that heated style. If the
        label is an empty string (i.e., ''), that moniker has not data.
        """

        heatedProcessedStyles = self.heatedProcessedStyles( )
        temperatureInfos = {}
        for instance in self :
            for style in heatedProcessedStyles :
                if( instance.moniker == style.moniker ) :
                    temperature = instance.temperature.getValueAs( unit )
                    if( temperature not in temperatureInfos ) :
                        temperatureInfos[temperature] = {}
                        for style2 in heatedProcessedStyles : temperatureInfos[temperature][style2.moniker] = ''
                    temperatureInfos[temperature][style.moniker] = instance.label
                    break

        return( [ [ temperaure, temperatureInfos[temperaure] ] for temperaure in sorted( temperatureInfos ) ] )

    def toXMLList( self, indent = '', **kwargs ) :
        """Returns an XML representation of self as a list of lines."""

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        do_multiGroup = kwargs['formatVersion'] == formatVersionModule.version_1_10

        parameters = ''
        if( do_multiGroup ) :
            heatedMultiGroups = [ _style for _style in self if isinstance( _style, heatedMultiGroup ) ]
            do_multiGroup = len( heatedMultiGroups ) > 0
            if( do_multiGroup ) :
                parameters = 'MultiGroup'
                multiGroup1 = multiGroup( parameters, 3 )
                for transportable in heatedMultiGroups[0].transportables : multiGroup1.transportables.add( transportable )

        xmlStringList = ['%s<%s>' % (indent, self.moniker)]
        for _style in self :
            if( do_multiGroup and isinstance( _style, heated ) ) :
                xmlStringList += multiGroup1.toXMLList( indent2, **kwargs )
                do_multiGroup = False
            xmlStringList += _style.toXMLList( indent2, **kwargs, parameters = parameters )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    def parseXMLNode( self, stylesElement, xPath, linkData ) :
        """
        Generate instances of allowed classes from the available XML nodes with the corresponding styles.

        Only the following styles are allowed
         - evaluated
         - crossSectionReconstructed
         - angularDistributionReconstructed
         - averageProductData
         - MonteCarlo_cdf
         - multiGroup       for GNDS 1.10 only
         - heated
         - heatedMultiGroup
         - griddedCrossSection
         - SnElasticUpScatter
         - CoulombPlusNuclearElasticMuCutoff
         - URR_probabilityTables
         - Realization

        """

        xPath.append( stylesElement.tag )

        classDict = {}
        for _style in ( evaluated, crossSectionReconstructed, angularDistributionReconstructed, Realization,
                        averageProductData, MonteCarlo_cdf, heated, heatedMultiGroup,
                        griddedCrossSection, SnElasticUpScatter, CoulombPlusNuclearElasticMuCutoff, URR_probabilityTables ) :
            classDict[_style.moniker] = _style

        multiGroupChildren = []
        for child in stylesElement :
            if( child.tag == multiGroup.moniker ) :       # Only in GNDS 1.10.
                multiGroupChildren.append( child )
            else :
                _class = classDict.get( child.tag, None )
                if( _class is None ) : raise TypeError( 'Encountered unknown style "%s".' % child.tag )
                self.add( _class.parseXMLNode( child, xPath, linkData ) )

        for child in multiGroupChildren :
            multiGroup1 = multiGroup.parseXMLNode( child, xPath, linkData )
            for style1 in self :
                if( isinstance( style1, heatedMultiGroup ) ) :
                    if( multiGroup1.label == style1.parameters ) :
                        for transportable in multiGroup1.transportables : style1.transportables.add( transportable )

        xPath.pop()

    @staticmethod
    def preProcessingStyles( ) :
        """Returns the list of style classes that can be processed. These are known as pre-processed styles."""

        return( evaluated, crossSectionReconstructed, angularDistributionReconstructed, Realization )

    @staticmethod
    def heatedProcessedStyles( ) :

        return( heated, griddedCrossSection, URR_probabilityTables, heatedMultiGroup, SnElasticUpScatter )

class style( ancestryModule.ancestry ) :
    """
    Abstract base class for the classes that represent actual GNDS style nodes.
    """

    __metaclass__ = abc.ABCMeta

    def __init__( self, label, derivedFrom, date = None ) :

        ancestryModule.ancestry.__init__( self )

        if( not( isinstance( label, str ) ) ) : raise TypeError( 'label must be a str instance.' )
        self.__label = label

        if( date is None ) : date = str( datetime.date.today( ) )
        self.__date = date

        if( not( isinstance( derivedFrom, str ) ) ) : raise TypeError( 'label must be a str instance.' )
        self.__derivedFrom = derivedFrom

        self.__documentation = documentationModule.Documentation( )
        self.__documentation.setAncestor( self )

    @property
    def date( self ) :
        """
        Property for the style release date.
        """

        return( self.__date )

    @property
    def derivedFrom( self ) :
        """
        Property indicating the parent style label from which the data in style is derived.
        """

        return( self.__derivedFrom )

    @property
    def derivedFromStyle( self ) :
        """
        Property for the style (not just the label) from which the data in this style is derived.
        """

        for _style in self.ancestor :
            if( _style.label == self.__derivedFrom ) : return( _style )
        return( None )

    @property
    def documentation( self ) :
        """Returns the documentation instance."""

        return( self.__documentation )

    @property
    def label( self ) :
        """
        Property for the string identifier of this style.
        """

        return( self.__label )

    @property
    def projectileEnergyDomain( self ) :
        """
        Property for the projectileEnergyDomain at which the data in this style is evaluated.
        """

        return( self.derivedFromStyle.projectileEnergyDomain )

    @property
    def temperature( self ) :
        """
        Property for the temperature at which the data in this style is evaluated.
        """

        return( self.derivedFromStyle.temperature )

    @property
    @abc.abstractmethod
    def sortOrderIndex( self ) :

        pass

    def chain( self ) :
        """Returns the list of derivedFrom style instances for self, ordered from self to root."""

        _chain = [ self ]
        while( True ) :
            _derivedFrom = _chain[-1].derivedFromStyle
            if( _derivedFrom is None ) : break
            _chain.append( _derivedFrom )

        return( _chain )

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
        """
        Method to find the derived style that corresponds to the input argument style.
        """

        _style = self
        while( _style is not None ) :
            _style = _style.derivedFromStyle
            if( ( _style is not None ) and ( isinstance( _style, cls ) ) ) : return( _style )
        return( None )

    def findInstancesOfClassInChildren( self, cls, level = 9999 ) :

        foundInstances = []
        level -= 1
        if( level < 0 ) : return( foundInstances )
        for item in self : foundInstances += item.findInstancesOfClassInChildren( cls, level )

        return( foundInstances )

    def hasTemperature( self ) :
        """
        Base method to indicate whether style has a temperature property.
        """

        return( False )

    def sibling( self, label ) :
        """Returns the sibling of the current style which has a label that correspond to the input argument."""

        return( self.ancestor[label] )

    def convertUnits( self, unitMap ) :
        """
        Base method to convert the style's units
        """

        pass

    def toXML( self, indent = "", **kwargs ) :
        """
        Convert an XML element from a list to a single newline delineated string.""
        """

        return( '\n'.join( self.toXMLList( **kwargs ) ) )

    def XMLCommonAttributes( self ) :
        """
        Return a string with the style label, and if available, the date and derivedFrom properties.
        """

        XMLCommon = 'label="%s"' % self.label
        if( self.derivedFrom != "" ) : XMLCommon += ' derivedFrom="%s"' % self.derivedFrom
        if( self.date is not None ) : XMLCommon += ' date="%s"' % self.date
        return( XMLCommon )

    @staticmethod
    def parseXMLNodeBase( element, xPath, linkData ) :
        """
        Update the XPath with the current XML element label and return the element label, derivedFrom and date attributes.

        For a given XML style element
         - the element tag is appended to the xPath; and
         - element label, derivedFrom and date attributes are returned as a tuple.
        """

        label = element.get( 'label' )
        xPath.append( '%s[@label="%s"]' % ( element.tag, label ) )

        derivedFrom = element.get( 'derivedFrom', "" )
        date = element.get( 'date', None )

        return( label, derivedFrom, date )

class styleWithTemperature( style ) :
    """
    Base class for the styles with a temperature attribute.
    """

    def __init__( self, label, derivedFrom, temperature, date = None ) :

        style.__init__( self, label, derivedFrom, date = date )

        if( not( isinstance( temperature, ( physicalQuantityModule.temperature, None.__class__ ) ) ) ) : raise TypeError( 'invalid temperature object' )
        self.__temperature = temperature

    @property
    def temperature( self ) :
        """
        Property for the style temperature.
        """

        if( self.__temperature is None ) : return( self.derivedFromStyle.temperature )
        return( self.__temperature )

    def convertUnits( self, unitMap ) :
        """
        Method to convert the units of the style temperature.
        """

        if( self.__temperature is not None ) : self.__temperature.convertUnits( unitMap )

    def hasTemperature( self ) :
        """
        Method indicating if the style has a temperature attribute.
        """

        return( self.__temperature is not None )

    @staticmethod
    def parseXMLNodeBase( element, xPath, linkData ) :
        """
        This methods reads the label, derivedFrom, date, and temperature attributes from the XML element and 
        return them as a tuple. The temperature is tranformed to an xData.physicalQuantity.physicalQantity
        object before it is returned in the tuple.
        """

        label, derivedFrom, date = style.parseXMLNodeBase( element, xPath, linkData )
        temperature = element.find( 'temperature' )
        if( temperature is not None ) : temperature = physicalQuantityModule.temperature.parseXMLNode( temperature, xPath, linkData )

        return( label, derivedFrom, date, temperature )

class evaluated( styleWithTemperature ) :
    """
    Style for data entered by an evaluator
    """

    moniker = 'evaluated'
    sortOrderIndex = 0

    def __init__( self, label, derivedFrom, temperature, projectileEnergyDomain, library, version, date = None ) :

        styleWithTemperature.__init__( self, label, derivedFrom, temperature, date = date )

        if( not( isinstance( library, str ) ) ) : raise TypeError( 'library must be a string' )
        self.__library = library

        if( not( isinstance( version, str ) ) ) : raise TypeError( 'version must be a string' )
        self.__version = version

        self.projectileEnergyDomain = projectileEnergyDomain

    @property
    def library( self ) :
        """
        Library property, i.e. the name of the library for this evaluation.
        """

        return( self.__library )

    @library.setter
    def library( self, value ) :
        """
        Setter for the library property, i.e. the name of the library for this evaluation.
        """

        self.__library = value

    @property
    def projectileEnergyDomain( self ) :
        """
        Projectile energy domain property.
        """

        if( self.__projectileEnergyDomain is None ) : return( self.derivedFromStyle.projectileEnergyDomain )
        return self.__projectileEnergyDomain

    @projectileEnergyDomain.setter
    def projectileEnergyDomain( self, value ) :
        """
        Setter for the projectile energy domain property.
        """

        if( not( isinstance( value, ( projectileEnergyDomain, None.__class__ ) ) ) ) :
            raise TypeError("Expected projectileEnergyDomain instance, got '%s' instead" % type(value))
        if( value is not None ) : value.setAncestor( self )
        self.__projectileEnergyDomain = value

    @property
    def version( self ) :
        """
        Version property, i.e. the version of this library's release.
        """

        return( self.__version )

    @version.setter
    def version( self, value ) :
        """
        Setter for the version property, i.e. the version of this library's release.
        """

        self.__version = value

    def convertUnits( self, unitMap ) :
        """
        Convert units for the style temperature and projectile energy domain.
        """

        styleWithTemperature.convertUnits( self, unitMap )
        if( self.__projectileEnergyDomain is not None ) : self.__projectileEnergyDomain.convertUnits( unitMap )

    def toXMLList( self, indent = '', **kwargs ) :
        """Return styles XML element and its child nodes as a list"""

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s %s library="%s" version="%s">' %
                ( indent, self.moniker, self.XMLCommonAttributes( ), self.library, self.version ) ]
        if( self.hasTemperature( ) ) : xmlStringList += self.temperature.toXMLList( indent2, **kwargs )
        if( self.__projectileEnergyDomain is not None ) : xmlStringList += self.projectileEnergyDomain.toXMLList( indent2, **kwargs )
        if( kwargs['formatVersion'] != formatVersionModule.version_1_10 ) :
            xmlStringList += self.documentation.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :
        """Generate an instance of this class from an XML node with the corresponding style"""

        label, derivedFrom, date, temperature = styleWithTemperature.parseXMLNodeBase( element, xPath, linkData )

        library = element.get( 'library' )
        version = element.get( 'version' )
        _projectileEnergyDomain = element.find( projectileEnergyDomain.moniker )
        if( _projectileEnergyDomain is not None ) : _projectileEnergyDomain = projectileEnergyDomain.parseXMLNode( _projectileEnergyDomain, xPath, linkData )

        _evaluated = evaluated( label, derivedFrom, temperature, _projectileEnergyDomain, library, version, date = date )

        _documentation = element.find( documentationModule.Documentation.moniker )
        if( _documentation is not None ) :
            _evaluated.documentation.parseNode( _documentation, xPath, linkData )

        xPath.pop( )
        return( _evaluated )

class projectileEnergyDomain( ancestryModule.ancestry ) :
    """
    Style for the storage of the energy domain and the unit of the projectile in the laboratory frame.
    """

    moniker = 'projectileEnergyDomain'

    def __init__(self, min, max, unit):

        ancestryModule.ancestry.__init__(self)
        self.min = min
        self.max = max
        self.unit = unit

    def convertUnits( self, unitMap ) :
        """
        Method for the conversion of the energy units.
        """

        if( self.unit in unitMap ) :
            factor = PQUModule.PQU( 1, self.unit ).getValueAs( unitMap[self.unit] )
            self.min *= factor
            self.max *= factor
            self.unit = unitMap[self.unit]

    def toXMLList( self, indent = '', **kwargs ) :
        """Returns an XML representation of self as a list of lines."""

        return ['%s<%s min="%s" max="%s" unit="%s"/>'
                % (indent, self.moniker, PQUModule.floatToShortestString( self.min ), PQUModule.floatToShortestString( self.max ), self.unit)]

    @classmethod
    def parseXMLNode(cls, element, xPath, linkData):
        """Generate an instance of this class from an XML node with the corresponding style"""

        xPath.append(element.tag)
        PED = cls( float(element.get('min')), float(element.get('max')), element.get('unit') )
        xPath.pop()
        return PED

class crossSectionReconstructed( style ) :
    """
    Style for cross section data reconstructed from resonance parameters.
    """

    moniker = 'crossSectionReconstructed'
    sortOrderIndex = 2

    def __init__( self, label, derivedFrom, date = None ) :

        style.__init__( self, label, derivedFrom, date = date )

    def toXMLList( self, indent = '', **kwargs ) :
        """Returns an XML representation of self as a list of lines."""

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s %s>' % ( indent, self.moniker, self.XMLCommonAttributes( ) ) ]
        xmlStringList += self.documentation.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :
        """Generate an instance of this class from an XML node with the corresponding style"""        

        label, derivedFrom, date = style.parseXMLNodeBase( element, xPath, linkData )

        _crossSectionReconstructed = crossSectionReconstructed( label, derivedFrom, date = date )

        _documentation = element.find( documentationModule.Documentation.moniker )
        if( _documentation is not None ) :
            _crossSectionReconstructed.documentation.parseNode( _documentation, xPath, linkData )

        xPath.pop()
        return( _crossSectionReconstructed )

class angularDistributionReconstructed( style ) :
    """
    Style for angular distribution data reconstructed from resonance parameters. 
    """

    moniker = 'angularDistributionReconstructed'
    sortOrderIndex = 3

    def __init__( self, label, derivedFrom, date = None ) :

        style.__init__( self, label, derivedFrom, date = date )

    def toXMLList( self, indent = '', **kwargs ) :
        """Returns an XML representation of self as a list of lines."""

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s %s>' % ( indent, self.moniker, self.XMLCommonAttributes( ) ) ]
        xmlStringList += self.documentation.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :
        """Generate an instance of this class from an XML node with the corresponding style"""

        label, derivedFrom, date = style.parseXMLNodeBase( element, xPath, linkData )

        _angularDistributionReconstructed = angularDistributionReconstructed( label, derivedFrom, date = date )

        _documentation = element.find(documentationModule.Documentation.moniker)
        if (_documentation is not None):
            _angularDistributionReconstructed.documentation.parseNode(_documentation, xPath, linkData)

        xPath.pop()
        return( _angularDistributionReconstructed )

class CoulombPlusNuclearElasticMuCutoff( style ) :
    """
    This is a derived data style for Coulomb interactions which have a singularity at 
    :math:`\\mu = 1.0`. The cross sections may consequently represented by the integral of 
    :math:`\\frac{d\\sigma \\left(E, \\mu\\right)}{d\\mu}` from :math:`\\mu = -1.0` to
    :math:`\\mu_{\\textrm{cutoff}}`, i.e. the cutoff scattering cosine which satisfies the condition 
    :math:`\\mu_{\\textrm{cutoff}} < 1.0`. 
    """

    moniker = 'CoulombPlusNuclearElasticMuCutoff'
    sortOrderIndex = 10

    def __init__( self, label, derivedFrom, muCutoff, date = None ) :

        style.__init__( self, label, derivedFrom, date = date )

        if( not( isinstance( muCutoff, float ) ) ) : raise TypeError( 'invalid muCutoff object' )
        self.__muCutoff = muCutoff

    @property
    def muCutoff( self ) :
        """Return the value for :math:`\\mu_{\\textrm{cutoff}}`"""

        return( self.__muCutoff )

    def toXMLList( self, indent = '', **kwargs ) :
        """Returns an XML representation of self as a list of lines."""

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s muCutoff="%s" %s>' % ( indent, self.moniker, self.__muCutoff, self.XMLCommonAttributes( ) ) ]
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :
        """Generate an instance of this class from an XML node with the corresponding style"""

        label, derivedFrom, date = style.parseXMLNodeBase( element, xPath, linkData )

        muCutoff = float( element.get( 'muCutoff' ) )
        _CoulombPlusNuclearElasticMuCutoff = CoulombPlusNuclearElasticMuCutoff( label, derivedFrom, muCutoff, date = date )

        xPath.pop( )
        return( _CoulombPlusNuclearElasticMuCutoff )

class heated( styleWithTemperature ) :
    """
    Style for data that has been heated to a specific temperature, e.g. Doppler broaded cross sections.
    """

    moniker = 'heated'
    sortOrderIndex = 100

    def __init__( self, label, derivedFrom, temperature, date = None ) :

        styleWithTemperature.__init__( self, label, derivedFrom, temperature, date = date )

    def toXMLList( self, indent = '', **kwargs ) :
        """Returns an XML representation of self as a list of lines."""

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s %s>' % ( indent, self.moniker, self.XMLCommonAttributes( ) ) ]
        xmlStringList += self.temperature.toXMLList( indent2, **kwargs )
        xmlStringList += self.documentation.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :
        """Generate an instance of this class from an XML node with the corresponding style"""

        label, derivedFrom, date, temperature = styleWithTemperature.parseXMLNodeBase( element, xPath, linkData )

        _heated = heated( label, derivedFrom, temperature, date = date )

        _documentation = element.find(documentationModule.Documentation.moniker)
        if (_documentation is not None):
            _heated.documentation.parseNode(_documentation, xPath, linkData)

        xPath.pop()
        return( _heated )

class averageProductData( style ) :
    """
    Style for data calculated for the averageProductEnergy and averageMomentum components.
    """

    moniker = 'averageProductData'
    sortOrderIndex = 11

    def __init__( self, label, derivedFrom, date = None ) :

        style.__init__( self, label, derivedFrom, date = date )

    def toXMLList( self, indent = '', **kwargs ) :
        """Returns an XML representation of self as a list of lines."""

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s %s>' % ( indent, self.moniker, self.XMLCommonAttributes( ) ) ]
        xmlStringList += self.documentation.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :
        """Generate an instance of this class from an XML node with the corresponding style"""        

        label, derivedFrom, date = style.parseXMLNodeBase( element, xPath, linkData )

        _averageProductData = averageProductData( label, derivedFrom, date = date )

        _documentation = element.find(documentationModule.Documentation.moniker)
        if (_documentation is not None):
            _averageProductData.documentation.parseNode(_documentation, xPath, linkData)

        xPath.pop()
        return( _averageProductData )

class multiGroup( style ) :
    """
    Style for the multigroup data. This style is deprecated (i.e., it does not exists it GNDS 2.0 or higher).
    """

    moniker = 'multiGroup'
    sortOrderIndex = -1

    def __init__( self, label, lMax, date = None ) :

        style.__init__( self, label, "", date = date )
        self.lMax = int( lMax )

        self.transportables = transportablesModule.transportables( )
        self.transportables.setAncestor( self )

    def toXMLList( self, indent = '', **kwargs ) :
        """Returns an XML representation of self as a list of lines."""

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s %s lMax="%s">' % ( indent, self.moniker, self.XMLCommonAttributes( ), self.lMax ) ]
        xmlStringList += self.transportables.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :
        """Generate an instance of this class from an XML node with the corresponding style"""

        label, derivedFrom, date = style.parseXMLNodeBase( element, xPath, linkData )

        lMax = element.get( 'lMax', 3 )

        _multiGroup = multiGroup( label, lMax, date = date )
        for child in element:
            if child.tag == transportablesModule.transportables.moniker:
                _multiGroup.transportables.parseXMLNode( child, xPath, linkData )
            else:
                raise TypeError("Encountered unexpected element '%s'" % child.tag)

        xPath.pop( )
        return( _multiGroup )

class heatedMultiGroup( style ) :
    """
    Style for data that has been heated and multi-grouped.
    """

    moniker = 'heatedMultiGroup'
    sortOrderIndex = 120

    def __init__( self, label, derivedFrom, flux, date = None, parameters = '' ) :

        style.__init__( self, label, derivedFrom, date = date )

        self.__transportables = transportablesModule.transportables( )
        self.__transportables.setAncestor( self )

        if( not( isinstance( flux, fluxModule.flux ) ) ) : raise TypeError( 'Invalid flux instance' )
        self.__flux = flux
        self.__flux.setAncestor( self )

        self.__inverseSpeed = None

        if( not( isinstance( parameters, str ) ) ) : raise TypeError( 'parameters must be a str instance' )
        self.parameters = parameters              # This is only used if GNDS 1.10.

        self.__multiGroupFlux = None                # Only for internal use, it is not written to GNDS file.

    @property
    def transportables( self ) :
        """
        Property for the transportables, i.e. the list of particles for which the multi-group processing was performed.
        """

        return( self.__transportables )

    @property
    def flux( self ) :
        """
        Property for the flux (used to collapse to the multi-group energy structure) as a fudge.processing.flux.flux object.
        """

        return( self.__flux )

    @property
    def inverseSpeed( self ) :
        """
        Property for the inverse speed.

        The group multi-group inverse speeds (:math:`\\left(\\upsilon_g\\right)`) are calculated from 
        :math:`\\upsilon_g = \\int_g dE \\frac{\\phi\\left(E\\right)}{\\upsilon\\left(E\\right)}`,
        where :math:`E` is the projectile energy, :math:`\\phi\\left(E\\right)` is the scalar flux,
        :math:`\\upsilon\\left(E\\right)` is the projectiles velocity, and the integral is over the
        energy range in group :math:`g`.
        """

        return( self.__inverseSpeed )

    @inverseSpeed.setter
    def inverseSpeed( self, value ) :
        """
        Setter for the inverseSpeed property.
        """

        self.__inverseSpeed = value

    @property
    def multiGroupFlux( self ) :
        """
        Property for the flux (used to collapse to the multi-group energy structure) as a fudge.processing.flux.gridded2d object.
        """

        return( self.__multiGroupFlux )

    def convertUnits( self, unitMap ) :
        # BRB - FIXME Still need to fully implement.

        self.__inverseSpeed.convertUnits( unitMap )

    def processMultiGroup( self, style, tempInfo, indent ) :
        """
        Method to set the multiGroupFlux and inverseSpeed properties.
        """

        self.__multiGroupFlux = self.flux.processMultiGroup( self, tempInfo, indent )
        self.inverseSpeed = inverseSpeedModule.inverseSpeed( inverseSpeedModule.multiGroupInverseSpeed( self, tempInfo ) )

    def toXMLList( self, indent = '', **kwargs ) :
        """Returns an XML representation of self as a list of lines."""

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        formatVersion = kwargs['formatVersion']
        parameters = kwargs.get( 'parameters', '' )

        attributes = self.XMLCommonAttributes( )
        if( parameters != '' ) : attributes += ' parameters="%s"' % parameters
        xmlStringList = [ '%s<%s %s>' % ( indent, self.moniker, attributes ) ]
        if( formatVersion != formatVersionModule.version_1_10 ) : xmlStringList += self.transportables.toXMLList( indent2, **kwargs )
        xmlStringList += self.flux.toXMLList( indent2, **kwargs )
        if( self.inverseSpeed is not None ) : xmlStringList += self.inverseSpeed.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :
        """Generate an instance of this class from an XML node with the corresponding style"""

        label, derivedFrom, date = style.parseXMLNodeBase( element, xPath, linkData )
        parameters = element.get( 'parameters', '' )

        transportables = None
        flux = None
        inverseSpeed = None
        for child in element :
            if( child.tag == transportablesModule.transportables.moniker ) :
                transportables = child
            elif( child.tag == fluxModule.flux.moniker ) :
                flux = fluxModule.flux.parseXMLNode( child, xPath, linkData )
            elif( child.tag == inverseSpeedModule.inverseSpeed.moniker ) :
                inverseSpeed = inverseSpeedModule.inverseSpeed.parseXMLNode( child, xPath, linkData )

        _heatedMultiGroup = heatedMultiGroup( label, derivedFrom, flux, date = date, parameters = parameters )
        if( transportables is not None ) : _heatedMultiGroup.transportables.parseXMLNode( transportables, xPath, linkData )
        _heatedMultiGroup.inverseSpeed = inverseSpeed

        xPath.pop( )
        return( _heatedMultiGroup )

class MonteCarlo_cdf( style ) :
    """
    This is a processed data style with distributions processed for use in Monte Carlo transport codes which
    requires:

        - Conversion of :math:`P\\left(\\mu, E ^{\\prime} \\mid E \\right)` to
          :math:`P\\left(\\mu \\mid E\\right) \\times P\\left(E ^{\\prime} \\mid E, \\mu\\right)`; and
        - Conversion of :math:`P\\left(E ^{\\prime}, \\mu \\mid E \\right)` to
          :math:`P\\left(E ^{\\prime} \\mid E \\right) \\times P\\left(\\mu \\mid E, E ^{\\prime} \\right)`

    """

    moniker = 'MonteCarlo_cdf'
    sortOrderIndex = 12

    def __init__( self, label, derivedFrom, date = None ) :

        style.__init__( self, label, derivedFrom, date = date )

    def toXMLList( self, indent = '', **kwargs ) :
        """Returns an XML representation of self as a list of lines."""


        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s %s>' % ( indent, self.moniker, self.XMLCommonAttributes( ) ) ]
        xmlStringList += self.documentation.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :
        """Generate an instance of this class from an XML node with the corresponding style"""

        label, derivedFrom, date = style.parseXMLNodeBase( element, xPath, linkData )

        _MonteCarlo_cdf = MonteCarlo_cdf( label, derivedFrom, date = date )

        _documentation = element.find(documentationModule.Documentation.moniker)
        if (_documentation is not None):
            _MonteCarlo_cdf.documentation.parseNode(_documentation, xPath, linkData)

        xPath.pop( )
        return( _MonteCarlo_cdf )

class SnElasticUpScatter( style ) :
    """
    This is a thermal elastic upscatter style for data that has been heated and multi-grouped.
    This is similar to the heatedMultigroup style but with an upscatter correction included for
    elastic scattering.
    """

    moniker = 'SnElasticUpScatter'
    sortOrderIndex = 121

    def __init__( self, label, derivedFrom, date = None ) :

        style.__init__( self, label, derivedFrom, date = date )
        
        self.__upperCalculatedGroup = 0 
        

    @property
    def upperCalculatedGroup( self ) :
        """Return or set  the index of the highest-energy group (counting from 0) impacted by the upscatter correction."""

        return( self.__upperCalculatedGroup )
    
    @upperCalculatedGroup.setter
    def upperCalculatedGroup( self, value ) :
        """Return or set  the index of the highest-energy group (counting from 0) impacted by the upscatter correction."""

        #if( not( isinstance( upperCalculatedGroup, int ) ) ) : raise TypeError( 'upperCalculatedGroup must be an integer' )
        self.__upperCalculatedGroup = value

    @property
    def flux( self ) :
        """Property for the flux used to collapse to group structure."""

        return( self.derivedFromStyle.flux )
    
    @property
    def transportables( self ) :
        """Property for the transportables, i.e. the list of particles for which the multi-group processing was performed."""

        return( self.derivedFromStyle.transportables )

    @property
    def inverseSpeed( self ) :
        """
        Property for the inverse speed.

        The group multi-group inverse speeds (:math:`\\left(\\upsilon_g\\right)`) are calculated from 
        :math:`\\upsilon_g = \\int_g dE \\frac{\\phi\\left(E\\right)}{\\upsilon\\left(E\\right)}`,
        where :math:`E` is the projectile energy, :math:`\\phi\\left(E\\right)` is the scalar flux,
        :math:`\\upsilon\\left(E\\right)` is the projectiles velocity, and the integral is over the
        energy range in group :math:`g`.
        """

        return( self.derivedFromStyle.inverseSpeed )

    def toXMLList( self, indent = '', **kwargs ) :
        """Returns an XML representation of self as a list of lines."""

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        xmlStringList = [ '%s<%s %s upperCalculatedGroup="%d">' % ( indent, self.moniker, self.XMLCommonAttributes( ), self.upperCalculatedGroup,  ) ]
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :
        """Generate an instance of this class from an XML node with the corresponding style"""

        label, derivedFrom, date = style.parseXMLNodeBase( element, xPath, linkData )

        styleOut = SnElasticUpScatter( label, derivedFrom, date = date )
        styleOut.upperCalculatedGroup = int(element.get('upperCalculatedGroup'))
        
        xPath.pop( )
        return( styleOut )

class griddedCrossSection( style ) :
    """
    Style for cross section data that has been converted to a common energy grid as used in Monte Carlo Transport.
    """

    moniker = 'griddedCrossSection'
    sortOrderIndex = 110

    def __init__( self, label, derivedFrom, date = None ) :

        style.__init__( self, label, derivedFrom, date = date )

        self.__grid = None

    @property
    def grid( self ) :
        """
        Common energy grid (as used in Monte Carlo transport) property.
        """

        return( self.__grid )

    @grid.setter
    def grid( self, grid ) :
        """
        Setter for the common energy grid (as used in Monte Carlo transport) property.
        """

        if( not( isinstance( grid, axesModule.grid ) ) ) : raise TypeError( 'grid not instance of axesModule.grid' )
        self.__grid = grid
        self.__grid.setAncestor( self )

    def convertUnits( self, unitMap ) :
        """
        Convert the style's energy grid units.
        """

        self.__grid.convertUnits( unitMap )

    def toXMLList( self, indent = '', **kwargs ) :
        """Returns an XML representation of self as a list of lines."""

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s %s>' % ( indent, self.moniker, self.XMLCommonAttributes( ) ) ]
        xmlStringList += self.__grid.toXMLList( indent = indent2, **kwargs )
        xmlStringList += self.documentation.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :
        """Generate an instance of this class from an XML node with the corresponding style"""

        label, derivedFrom, date = style.parseXMLNodeBase( element, xPath, linkData )

        grid = None
        for child in element :
            if( child.tag == axesModule.grid.moniker ) :
                grid = axesModule.grid.parseXMLNode( child, xPath, linkData )

        _griddedCrossSection = griddedCrossSection( label, derivedFrom, date = date )
        _griddedCrossSection.grid = grid

        _documentation = element.find(documentationModule.Documentation.moniker)
        if (_documentation is not None):
            _griddedCrossSection.documentation.parseNode(_documentation, xPath, linkData)

        xPath.pop( )
        return( _griddedCrossSection )

class URR_probabilityTables( style ) :
    """
    This is a precessed style for cross section probability tables :math:`P\\left(\\sigma \\mid E\\right)` 
    to better capture the possible rapid cross section fluctuations in the unresolved resonance region.
    """

    moniker = 'URR_probabilityTables'
    sortOrderIndex = 111

    def __init__( self, label, derivedFrom, date = None ) :

        style.__init__( self, label, derivedFrom, date = date )

    def convertUnits( self, unitMap ) :
        """
        unitMap is a dictionary with old/new unit pairs where the old unit is the key (e.g., { 'eV' : 'MeV', 'b' : 'mb' }).
        """

        pass

    def toXMLList( self, indent = '', **kwargs ) :
        """Returns an XML representation of self as a list of lines."""

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s %s/>' % ( indent, self.moniker, self.XMLCommonAttributes( ) ) ]
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :
        """Generate an instance of this class from an XML node with the corresponding style"""

        label, derivedFrom, date = style.parseXMLNodeBase( element, xPath, linkData )

        URR_probability_tables = URR_probabilityTables( label, derivedFrom, date = date )

        xPath.pop( )
        return( URR_probability_tables )

class Realization( style ) :
    """
    This style represents a change to the mean value of some of the data.
    """

    moniker = 'realization'
    sortOrderIndex = 1

    def __init__( self, label, derivedFrom, date = None ) :

        style.__init__( self, label, derivedFrom, date = date )

    def convertUnits( self, unitMap ) :
        """
        unitMap is a dictionary with old/new unit pairs where the old unit is the key (e.g., { 'eV' : 'MeV', 'b' : 'mb' }).
        """

        pass

    def toXMLList( self, indent = '', **kwargs ) :
        """Returns an XML representation of self as a list of lines."""

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        return( [ '%s<%s %s/>' % ( indent, self.moniker, self.XMLCommonAttributes( ) ) ] )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :
        """Returns a parsed instance of this class from an XML node."""

        label, derivedFrom, date = style.parseXMLNodeBase( element, xPath, linkData )

        realization = Realization( label, derivedFrom, date = date )

        xPath.pop( )
        return( realization )

def findEvaluated( forms ) :
    """
    Find the "evaluated" style class from a list of style classes. It calls a method that ensures that 
    the list of style classes contains exactly one occurence of the "evaluated" style.
    """

    return( findStyle( evaluated, forms ) )

def findAllOfStyle( cls, forms ) :
    """
    Find all the occurences of a guven style class from a list of style classes.
    """

    formsFound = []
    for form in forms :
        if( isinstance( form, cls ) ) : formsFound.append( form )
    return( formsFound )

def findStyle( cls, forms ) :
    """
    Find all occurences of a given style class and ensure that there is only a single occurence.
    """

    formsFound = findAllOfStyle( cls, forms )
    if( len( formsFound )  > 1 ) : raise Exception( 'multiple styles found' )
    if( len( formsFound ) == 0 ) : raise Exception( 'no style found' )
    return( formsFound[0] )
