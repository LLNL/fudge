# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
import abc

from xData import ancestry as ancestryModule

from LUPY import times

from . import suites as suitesModule
from . import styles as stylesModule

__metaclass__ = type

class component( suitesModule.suite ) :

    __metaclass__ = abc.ABCMeta

    def __init__( self, allowedClasses ) :

        suitesModule.suite.__init__( self, allowedClasses )

    @property
    def evaluated( self ) :

        if not hasattr( self.getRootAncestor(), 'styles' ) :
            return self[0]                  # A hack to deal with orphaned components that are not part of a full reactionSuite.

        evaluated = self.getRootAncestor( ).styles.getEvaluatedStyle( )
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

    def getStylesOfClass( self, cls ) :

        styles = self.getRootAncestor( ).styles.getStylesOfClass( cls )
        formList = []
        for _form in self :
            for style in styles :
                if( _form.label == style.label ) : formList.append( _form.label )
        return( formList )

    def getStyleOfClass( self, cls ) :

        style = self.getRootAncestor( ).styles.getStyleOfClass( cls )
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

            if( isinstance( style, ( stylesModule.MonteCarlo_cdf, stylesModule.griddedCrossSection ) ) ) : return( False )
            return( True )

        indent2 = indent + tempInfo['incrementalIndent']
        verbosity = tempInfo['verbosity']
        addToComponent = tempInfo.get( 'addToComponent', True )

        if( verbosity > 3 ) : time = times.times( )

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

        if( not hasattr( self.getRootAncestor( ), 'styles' ) ) :            # An orphaned component that is not part of a full reactionSuite.
            linear = self[0]
        else :
            linear = self[self.getRootAncestor( ).styles.preProcessingHeadInChainWithLabel( label ).label]

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
                warnings.append( warning.context('form "%s":' % form.label, formWarnings) )

        return warnings

class form( ancestryModule.ancestry ) :
    """
    This is the base class which is used as a form base class for channelData, reactionData and productData form classes.
    """

    __metaclass__ = abc.ABCMeta

    def getComponentsClass( self ) :

        return( sys.modules[self.__module__].component )

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent = indent, **kwargs ) ) )
