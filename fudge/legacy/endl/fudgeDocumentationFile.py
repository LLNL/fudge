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
