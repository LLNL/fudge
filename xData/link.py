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

__metaclass__ = type

import base as baseModule

class link( baseModule.xDataCoreMembers ) :
    """
    This class contains the path to another element, that could be stored in a different file.
    The 'follow' member function is used to get a pointer to the linked-to element,
    which will then be stored as self.link.
    
    Member data ::
    
        * label : the xml tag name (type: str)
        * link  : actual pointer to other data (type: str)
        * root  : file name where other data is stored, only needed if it's not the current file (type: str)
        * path  : self.path is an xPath expression that uniquely identifies the other data (type: str)

    sample links::
    <element xlink:type="simple" xlink:href='/reactionSuite/reaction[@label="102"]'>
        or
    <element xlink:type="simple" xlink:href='/covariances.xml#covId'>
    """

    moniker = 'link'

    def __init__( self, link = None, root = None, path = None, label = None, **attributes ) :

        baseModule.xDataCoreMembers.__init__( self, self.moniker, index = None, label = label )
        self.link = link    # actual pointer to other data
        self.root = root    # file name where other data is stored, only needed if it's not the current file
        self.path = path    # self.path is an xPath expression that uniquely identifies the other data

            # careful with nesting single/double quotes.
        if self.path is not None: self.path = self.path.replace('"',"'")
        self.attributes = attributes

    def __str__( self ) :

        return( "%s" % self.path )

    def __getitem__( self, key ):
        return self.attributes[key]

    def __setitem__( self, key, value ):
        self.attributes[key] = value

    def updateXPath( self ):
        '''ensure the xPath agrees with the linked-to data'''
        if hasattr(self.link, 'toXLink'): self.path = self.link.toXLink()

    def follow( self, startNode ):
        """
        :param startNode: instance corresponding to the beginning of self.path
        :return: class instance pointed to by self.path

        Uses ancestry.findEntity to find each element
        """
        import re
        # use regex on xPath expressions like "reaction[@label='2']"
        regex = re.compile("([a-zA-Z_]+)\[@([a-z]+)='([a-zA-Z0-9_]+|[a-zA-Z]*\([a-zA-Z0-9_,]+\))'\]")

        def follow2( xPathList, node ):
            # recursive helper function: descend the path to find the correct element
            if len(xPathList)==0: return node

            xPathNext = xPathList[0]
            match = regex.match(xPathNext)
            try:
                if match:
                    nodeNext = node.findEntity( *match.groups() )
                else:
                    nodeNext = node.findEntity( xPathNext )
            except:
                raise UnresolvableLink()
            return follow2(xPathList[1:], nodeNext)

        xPathList = self.path.split('/')
        while not xPathList[0]: # trim empty sections from the beginning
            xPathList = xPathList[1:]
        try:
            return follow2( xPathList, startNode )
        except UnresolvableLink:
            raise UnresolvableLink( "Cannot locate path '%s'" % self.path )


    def toXML( self, indent = '' , **kwargs ) :
        """Pointers show up in the attributes list on an xml element
        i.e., <element xlink:href="...xpath" anotherAttribute="foo" ...
        """

        self.updateXPath()
        attributesStr = self.attributesToXMLAttributeStr( )
        for key in sorted( self.attributes.keys( ) ) :
            attributesStr += ' %s="%s"' % ( key, self.attributes[key] )
        if( type( self.root ) == str ) :    # external link
            return '%s<%s%s xlink:href="%s#%s"/>' % ( indent, self.moniker, attributesStr, self.root, self.path )
        else:                               # link within this file
            return '%s<%s%s xlink:href="%s"/>' % ( indent, self.moniker, attributesStr, self.path )

    def toXMLList( self, indent = '', **kwargs ) :

        return [ self.toXML( indent, **kwargs ) ]

    @classmethod
    def parseXMLNode( cls, linkElement, xPath=[], linkData={} ):
        """
        Parse the xml-represented link back to python. The resulting link points to None,
        and must be resolved by the calling function.
        Link attributes can be converted to desired type by passing a 'typeConversion' dictionary in linkData
        """

        xPath.append( linkElement.tag )
        xlink_namespace = '{http://www.w3.org/1999/xlink}'
        path = linkElement.get(xlink_namespace+'href')
        if '#' in path: root, path = path.split('#')
        else: root = None
        # all optional (non-xlink) attributes:
        attributes = dict( [v for v in linkElement.items() if xlink_namespace not in v[0]] )
        for key in attributes:
            if key in linkData.get('typeConversion',{}):
                attributes[key] = linkData['typeConversion'][key](attributes[key])
        result = cls( link=None, root=root, path=path, **attributes )
        if( 'unresolvedLinks' in linkData ) : linkData['unresolvedLinks'].append( result )
        xPath.pop()

        return result


class UnresolvableLink( Exception ):
    pass
