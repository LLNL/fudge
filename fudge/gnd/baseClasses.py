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
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# <<END-copyright>>

from fudge.core.math import fudge2dGrouping

import xData.ancestry as ancestryModule

__metaclass__ = type

class componentBase( ancestryModule.ancestry ) :
    """This is the base class which is used as a component base class for channelData, reactionData and productData classes."""
    # FIXME should be an abstract base class
    # FIXME also switch 'genre' to 'moniker'?

    def __init__( self ) :

        self.moniker = self.genre
        ancestryModule.ancestry.__init__( self )
        self.forms = {}
        self.nativeData = None          # BRB, This will go away soon. Only a kludge for now.

    def __len__( self ) :
        "Returns the number of forms in self."

        return( len( self.forms ) )

    def __getitem__( self, label ) :
        "Returns self's form with label."

        return( self.forms[label] )

    def addForm( self, form ) :

        if( form.genre != self.genre ) : raise TypeError( "cannot add form with genre = '%s' to genre = '%s'" % ( form.genre, self.genre ) )
        if( form.label in self.forms ) : raise Exception( 'label = %s already in self' % form.label )
        self.forms[form.label] = form
        form.setAncestor( self )
        if( len( self.forms ) == 1 ) : self.nativeData = form.label

    def getStylesOfType( self, cls ) :

        styles = self.getRootAncestor( ).styles.getStylesOfType( cls )
        formList = []
        for formName in self.forms :
            for style in styles :
                if( formName == style.label ) : formList.append( self.label )
        return( formList )

    def getStyleOfType( self, cls ) :

        style = self.getRootAncestor( ).styles.getStyleOfType( cls )
        if( style is not None ) : style = self.forms[style.label]
        return( style )

    def getEvaluated( self ) :

        evaluated = self.getRootAncestor( ).styles.getEvaluatedStyle( )
        return( self.forms[evaluated.label] )

    def getFormTokens( self ) :
        """Returns a list of the name of each form currently in self."""

        return( self.forms.keys( ) )

# BRB, deprecated.
    def getFormByToken( self, formToken ) :
        "See __getitem__."

        return( self[formToken] )

    def getGenre( self ) :

        return( self.genre )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        raise Exception( "'%s' does not support processing" % self.genre )

    def removeStyle( self, form ) :
        """
        This method removes form from self. Form can be a label or a form instance.
        """

        if( isinstance( form, formBase ) ) :
            found = False
            for f in self.forms :
                if( self.forms[f] is form ) : found = True
            if( not( found ) ) : raise Exception( 'Form "%s" of genre "%s" is not a form of self with genre = "%s"' % ( form.label, form.genre, self.genre ) )
            form = form.label
        return( self.forms.pop( form, None ) )

    def toPointwise_withLinearXYs( self, lowerEps = 1.e-8, upperEps = 1.e-8 ) :
        """This method calls the toPointwise_withLinearXYs method for the evaluated style and returns a "lin,lin" pointwise representation of it."""

        evaluated = self.getEvaluated( )
        return( evaluated.toPointwise_withLinearXYs( lowerEps, upperEps ) )

    def toXML( self, indent = "", **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = "", **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        formXMLString = []
        for label in self.forms : formXMLString += self.forms[label].toXMLList( indent2, **kwargs )
        if( formXMLString == [] ) : return( [] )
        xmlString = [ '%s<%s>' % ( indent, self.genre ) ] + formXMLString
        xmlString[-1] += '</%s>' % self.genre
        return( xmlString )

class formBase( ancestryModule.ancestry ) :
    """This is the base class which is used as a form base class for channelData, reactionData and productData form classes."""
    # FIXME make it abstract?

    def getComponentsClass( self ) :

        import sys
        return( sys.modules[self.__module__].component )

    def getGenre( self ) :

        return( self.genre )

    def getForm( self ) :

        return( self.label )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        raise Exception( "processing of '%s' data from '%s' to '%s' is not supported" % ( self.genre, self.label, processInfo( 'to' ) ) )

class groupedFormBase( formBase, fudge2dGrouping.groupedData ) :
    """This is the base class which is used as a form grouped base class for channelData, reactionData and productData form classes."""

    moniker = 'grouped' 

    def __init__( self, label, axes_, groupData ) :

        fudge2dGrouping.groupedData.__init__( self, groupData, axes_ )

        if( label is not None ) : 
            if( not( isinstance( label, str ) ) ) : raise TypeError( 'label must be a string' )
        self.__label = label

    @property
    def label( self ) :

        return( self.__label )

class groupedWithCrossSectionFormBase( formBase, fudge2dGrouping.groupedData ) :
    """This is the base class which is used as a form grouped base class for channelData, reactionData and productData form classes."""

    moniker = 'groupedWithCrossSection' 

    def __init__( self, label, axes_, groupData, toForms = {} ) :

        fudge2dGrouping.groupedData.__init__( self, label, groupData, axes_ )
        self.toForms = toForms
