# <<BEGIN-copyright>>
# <<END-copyright>>

class fudgeDocumentationFile :

    def __init__( self, fileName ) :
        """The __init__ for fudgeDocumentationFile. This method does not read in any data, it only initializes some variables."""

        self.fileName = fileName
        self.text = None
        self.XML_ParserFriendlyText = None

    def read( self ) :
        """Read in the documentation for self. Note, this is not done by the __init__ method."""

        if( self.text == None ) :
            try :
                f = open( self.fileName )
                self.text = f.read( )
                f.close( )
            except :
                raise Exception( "\nError in fudgeDocumentationFile.read: could not open documentation file %s" % self.fileName )

    def getText( self ) :
        """Returns the documentation for self. Calls read if file has not been read."""

        if( self.text == None ) : self.read( )
        return( self.text )

    def getXML_ParserFriendlyText( self ) :
        """Converts &, < and > in the self.text to '&amp;', '&lt;' and '&gt;', respectively and stores it in self.XML_ParserFriendlyText"""

        if( self.XML_ParserFriendlyText == None ) :
            lines = self.getText( )
            if( lines != None ) :
                lines = lines.replace( '&', '&amp;' )
                lines = lines.replace( '<', '&lt;' )
                self.XML_ParserFriendlyText = lines.replace( '>', '&gt;' )
        return( self.XML_ParserFriendlyText )

    def setText( self, text ) :

        self.text = text
        self.XML_ParserFriendlyText = None

    def write( self, fileName = None ) :
        """If self.text is not equal to None, then its contents are written to file fileName. Self.text must be a python string. 
        If fileName is None then self.fileName is used."""

        if( self.text == None ) : return
        if( fileName == None ) : fileName = self.fileName
        try :
            f = open( fileName, 'w' )
            f.write( self.text )
            f.close( )
        except :
            raise Exception( "\nError in fudgeDocumentationFile.write: could not open documentation file %s" % self.fileName )
