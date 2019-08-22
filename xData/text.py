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

import base as baseModule

__metaclass__ = type

asciiEncodingToken = 'ascii'
utf8EncodingToken = 'utf8'
allowedEncodings = [ asciiEncodingToken, utf8EncodingToken ]

noneMarkupToken = 'none'
xmlMarkupToken = 'xml'
htmlMarkupToken = 'html'
latexMarkupToken = 'latex'

class text( baseModule.xDataCoreMembers ) :

    moniker = 'text'

    def __init__( self, textStr, encoding = asciiEncodingToken, markup = noneMarkupToken, allowedEncodings = allowedEncodings, 
            index = None, label = None ) :

        baseModule.xDataCoreMembers.__init__( self, self.moniker, index = index, label = label )

        if( not( isinstance( encoding, str ) ) ) : raise TypeError( 'encoding must a string' )
        if( encoding not in allowedEncodings ) : raise TypeError( 'Invalid encoding = "%s"' % encoding[:64] )
        self.__encoding = encoding

        if( not( isinstance( markup, str ) ) ) : raise TypeError( 'markup must a string' )
        self.__markup = markup

        if( not( isinstance( textStr, str ) ) ) : raise TypeError( 'text must a string' )
        self.__text = textStr

    @property
    def encoding( self ) :

        return( self.__encoding )

    @property
    def markup( self ) :

        return( self.__markup )

    @property
    def text( self ) :

        return( self.__text )

    def copy( self ) :

        return( text( self.text, self.encoding, self.markup, allowedEncodings = [ self.encoding ] ) )

    __copy__ = copy
    __deepcopy__ = __copy__

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent = indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

# BRB: Caleb, do we need to convert special XML characters?
        attributeStr = baseModule.xDataCoreMembers.attributesToXMLAttributeStr( self )
        if( self.encoding != asciiEncodingToken ) : attributeStr = ' encoding="%s"' % self.encoding
        if( self.markup != noneMarkupToken ) : attributeStr = ' markup="%s"' % self.markup
        XMLList = [ '%s<%s%s><![CDATA[%s]]></%s>' % ( indent, self.moniker, attributeStr, self.text, self.moniker ) ]
        return( XMLList )
