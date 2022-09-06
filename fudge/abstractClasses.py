# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
import abc

from LUPY import ancestry as ancestryModule

from LUPY import times as timesModule

from . import suites as suitesModule
from . import styles as stylesModule


class Component( suitesModule.Suite, abc.ABC ) :

    def __init__( self, allowedClasses ) :

        suitesModule.Suite.__init__( self, allowedClasses )

    def __getitem__(self, key):
        """Returns entry with key value *key*."""

        instance = suitesModule.Suite.__getitem__(self, key)
        if isinstance( instance, LazyParsingHelper):
            instance = instance.parse()
            self.replace(instance)

        return instance

    def __iter__(self):

        hrefInstance = self.hrefInstance( )
        if( hrefInstance is not None ) :
            hrefInstance.__iter__( )
        else :
            n1 = len( self )
            for i1 in range( n1 ) :
                instance = self.__getitem__(i1)
                yield instance

    @property
    def evaluated( self ) :

        if not hasattr( self.rootAncestor, 'styles' ) :
            return self[0]                  # A hack to deal with orphaned components that are not part of a full reactionSuite.

        evaluated = self.rootAncestor.styles.getEvaluatedStyle( )
        try :
            return( self[evaluated.label] )
        except :
            return( self[0] )

    def keys( self ):

        return [form.label for form in self]

    def amendForPatch( self, fromLabel, toLabel ) :

        for key in self.keys( ) :
            if( key == fromLabel ) :
                self[key].label = toLabel
            else :
                self.pop( key )

    def cullStyles( self, styleList ) :

        keeperLabel = ''
        for style in styleList :                        # Determine form to keep.
            for label in self.labels( ) :
                if( label == style.label ) :
                    keeperLabel = label
                    break
            if( keeperLabel != '' ) : break

        for label in self.labels( ) :                   # Remove all but keeping form.
            if( label != keeperLabel ) : self.pop( label )

    def diff(self, other, diffResults):
        """
        Check the first form of *self* with that of *other*.
        """

        if len(self) == 0 and len(other) == 0: return

        if len(self) == 0:
            diffResults.append('%s missing - 1' % self.moniker, '', other[0].toXLink(), '')
        elif len(other) == 0:
            diffResults.append('%s missing - 2' % self.moniker, '', self[0].toXLink(), '')
        else:
            self[0].diff( other[0], diffResults )

    def getStylesOfClass( self, cls ) :

        styles = self.rootAncestor.styles.getStylesOfClass( cls )
        formList = []
        for _form in self :
            for style in styles :
                if( _form.label == style.label ) : formList.append( _form.label )
        return( formList )

    def getStyleOfClass( self, cls ) :

        style = self.rootAncestor.styles.getStyleOfClass( cls )
        if( style is not None ) : style = self[style.label]
        return( style )

    def processMC_cdf( self, style, tempInfo, indent = '' ) :

        indent2 = indent + tempInfo['incrementalIndent']

        addToComponent = tempInfo.get( 'addToComponent', True )

        _form = style.findFormMatchingDerivedStyle( self )
        try :
            MonteCarlo_cdf = _form.processMC_cdf( style, tempInfo, indent2 )
        except :
            print('==== form = %s failed' % _form.moniker)
            raise
        if( addToComponent and ( MonteCarlo_cdf is not None ) ) : self.add( MonteCarlo_cdf )

    def processMultiGroup( self, style, tempInfo, indent ) :

        def styleFilter( style ) :

            if( isinstance( style, ( stylesModule.MonteCarlo_cdf, stylesModule.GriddedCrossSection ) ) ) : return( False )
            return( True )

        indent2 = indent + tempInfo['incrementalIndent']
        verbosity = tempInfo['verbosity']
        addToComponent = tempInfo.get( 'addToComponent', True )

        if( verbosity > 3 ) : time = timesModule.Times( )

        _form = style.findFormMatchingDerivedStyle( self, styleFilter )
        if( _form is None ) :
#            if( verbosity ) : print('%sWARNING: no form found for component %s' % ( indent, self.moniker ))
            return( None )
        try :
            multiGroup = _form.processMultiGroup( style, tempInfo, indent2 )
        except :
            sys.stderr.write( 'form = %s: %s\n' % ( _form.moniker, self ) )
            raise

        if( verbosity > 3 ) : print('%s%s: %s' % ( indent2 + '    ', self.moniker, time.toString( current = False ) ))

        if( addToComponent and ( multiGroup is not None ) ) : self.add( multiGroup )
        return( multiGroup )

    def removeStyles( self, styleLabels ) :
        """
        Removes all forms whose label matches one of the labels in removeStyles.
        """

        for label in styleLabels : self.remove( label )

    def toLinear( self, label = None, **kwargs ) :

        if( not hasattr( self.rootAncestor, 'styles' ) ) :            # An orphaned component that is not part of a full reactionSuite.
            linear = self[0]
        else :
            linear = self[self.rootAncestor.styles.preProcessingHeadInChainWithLabel( label ).label]

        return( linear.toPointwise_withLinearXYs( **kwargs ) )              # This is not correct. Need toLinear for each class and use it.

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        This method calls the toPointwise_withLinearXYs method for the evaluated style and returns a "lin-lin" pointwise representation of it.

        @:param lowerEps: lower epsilon for converting discontinuities to monotonic function
        @:param upperEps: upper epsilon for converting discontinuities to monotonic function
        """

        evaluated = self.evaluated
        styleLabel = kwargs.get( 'styleLabel', evaluated.label )
        if( styleLabel != evaluated.label ) : evaluated = self[styleLabel]

        return( self.evaluated.toPointwise_withLinearXYs( **kwargs ) )

    def check( self, info ):

        from fudge import warning
        warnings = []

        for form in self:
            formWarnings = form.check( info )
            if formWarnings:
                warnings.append( warning.Context('form "%s":' % form.label, formWarnings) )

        return warnings

    def parseNode(self, node, xPath, linkData, **kwargs):

        if node is None: return

        xPath.append(node.tag)

        lazyParsing = kwargs.get('lazyParsing', False)

        for child in node:
            parseClass = None
            for class1 in self.allowedClasses:
                tag = self.legacyMemberNameMapping.get(child.tag, child.tag)
                if tag == class1.moniker:
                    parseClass = class1
                    break
            if parseClass is None: raise TypeError("Invalid node '%s' encountered in suite '%s'" % (child.tag, self.moniker))

            if lazyParsing:
                self.add(LazyParsingHelper('label', child, parseClass, linkData, **kwargs), addLazyParsingHelper=True)
            else:
                self.add(parseClass.parseNodeUsingClass(child, xPath, linkData, **kwargs))

        href = label = node.get('href')
        if href is not None: self.set_href(href)

        xPath.pop()

class LazyParsingHelper(ancestryModule.Ancestry):
    """
    This class stores the required data needed to parse a node. When parsing is needed, the *parse* method will returned the constructed instance.
    """

    moniker = 'LazyParsingHelper'

    def __init__(self, keyName, node, Class, linkData, **kwargs):
        """Constructor for LazyParsingHelper class."""

        self.__keyName = keyName
        self.__keyValue = node.get(keyName)

        self.node = node
        self.Class = Class
        self.linkData = linkData

        if 'LazyParsingHelperCounter' not in linkData: linkData['LazyParsingHelperCounter'] = 0
        linkData['LazyParsingHelperCounter'] += 1
        self.kwargs = kwargs

    @property
    def keyName(self):
        """Returns the *keyName* for *self*."""

        return self.__keyName

    @property
    def keyValue(self):
        """Returns the *keyValue* for *self*."""

        return self.__keyValue

    label = keyValue

    def parse(self):
        """Parses the stored node and returns the constructued instance."""

        self.linkData['LazyParsingHelperCounter'] -= 1
        return self.Class.parseNodeUsingClass(self.node, [], self.linkData, **self.kwargs)

class Form( ancestryModule.AncestryIO, abc.ABC ) :
    """
    This is the base class which is used as a form base class for outputChannelData, reactionData and productData form classes.
    """

    def getComponentsClass( self ) :

        return( sys.modules[self.__module__].Component )
