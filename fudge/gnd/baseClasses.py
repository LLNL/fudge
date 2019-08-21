# <<BEGIN-copyright>>
# <<END-copyright>>

from fudge.core.ancestry import ancestry
from fudge.core.math import fudge2dGrouping
from fudge.gnd import tokens

__metaclass__ = type

class componentBase( ancestry ) :
    """This is the base class which is used as a component base class for channelData, reactionData and productData classes."""

    def __init__( self, allowedForms, nativeData = tokens.unknownFormToken ) :

        self.allowedForms = allowedForms
        self.nativeData = nativeData
        self.forms = {}
        ancestry.__init__( self, self.genre, None )

    def __len__( self ) :
        "Returns the number of forms representing this data."

        return( len( self.forms ) )

    def __getitem__( self, formToken ) :
        "Returns the form with token formToken for self."

        return( self.forms[formToken] )

    def addForm( self, form ) :

        if( form.genre != self.genre ) : raise Exception( "cannot add form with genre = '%s' to genre = '%s'" % ( form.genre, self.genre ) )
        if( form.form not in self.allowedForms ) : raise Exception( "Unsupported form = '%s' for genre = '%s'" % ( form.form, self.genre ) )
        self.forms[form.form] = form
        if( len( self.forms ) == 1 and self.nativeData == tokens.unknownFormToken ) : self.nativeData = form.form
        form.setParent( self )

    def addFormAsNativeData( self, form ) :

        self.addForm( form )
        self.nativeData = form.form

    def getFormTokens( self ) :
        """Returns a list of the name of each form currently in self."""

        return( self.forms.keys( ) )

    def getFormByToken( self, formToken ) :
        "See __getitem__."

        return( self[formToken] )

    def getGenre( self ) :

        return( self.genre )

    def getNativeDataToken( self ) :

        return( self.nativeData )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        raise Exception( "'%s' does not support processing" % self.genre )

    def removeForm( self, form, force = False ) :
        """This method removes form from self. Form can be a token or a form instance. The nativeData of self can only be removed if force is True."""

        if( isinstance( form, formBase ) ) :
            found = False
            for f in self.forms :
                if( self.forms[f] is form ) : found = True
            if( not( found ) ) : raise Exception( 'Form "%s" of genre "%s" is not a form of self with genre = "%s"' % ( form.form, form.genre, self.genre ) )
            form = form.form
        if( self.getNativeDataToken( ) == form ) :
            if( not( force ) ) : raise Exception( 'Cannot remove nativeData unless force is True' )
            self.nativeData = tokens.unknownFormToken
        return( self.forms.pop( form, None ) )

    def toXMLList( self, indent = "" ) :

        indent2 = indent + '  '
        formXMLString = []
        for form in self.forms : formXMLString += self.forms[form].toXMLList( indent = indent2 )
        if( formXMLString == [] ) : return( [] )
        xmlString = [ '%s<%s nativeData="%s">' % ( indent, self.genre, self.nativeData ) ] + formXMLString
        xmlString[-1] += '</%s>' % self.genre
        return( xmlString )

class formBase( ancestry ) :
    """This is the base class which is used as a form base class for channelData, reactionData and productData form classes."""

    def __init__( self ) :

        ancestry.__init__( self, self.form, None )

    def getComponentsClass( self ) :

        import sys
        return( sys.modules[self.__module__].component )

    def getGenre( self ) :

        return( self.genre )

    def getForm( self ) :

        return( self.form )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        raise Exception( "processing of '%s' data from '%s' to '%s' is not supported" % ( self.genre, self.form, processInfo( 'to' ) ) )

class groupedFormBase( formBase, fudge2dGrouping.groupedData ) :
    """This is the base class which is used as a form grouped base class for channelData, reactionData and productData form classes."""

    form = tokens.groupedFormToken
    tag = tokens.groupedFormToken

    def __init__( self, axes_, groupData, toForms = {} ) :

        formBase.__init__( self )
        fudge2dGrouping.groupedData.__init__( self, groupData, axes_ )
        self.toForms = toForms

class groupedWithCrossSectionFormBase( formBase, fudge2dGrouping.groupedData ) :
    """This is the base class which is used as a form grouped base class for channelData, reactionData and productData form classes."""

    form = tokens.groupedWithCrossSectionFormToken
    tag = tokens.groupedWithCrossSectionFormToken

    def __init__( self, axes_, groupData, toForms = {} ) :

        formBase.__init__( self )
        fudge2dGrouping.groupedData.__init__( self, groupData, axes_ )
        self.toForms = toForms
