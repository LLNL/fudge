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

    def __init__( self, link = None, root = None, path = None, label = None, relative = False, **attributes ) :

        baseModule.xDataCoreMembers.__init__( self, self.moniker, index = None, label = label )
        self.link = link    # actual pointer to other data
        self.root = root    # file name where other data is stored, only needed if it's not the current file
        self.path = path    # self.path is an xPath expression that uniquely identifies the other data
        self.__relative = relative      # whether to use relative link when writing out xPath

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
        """ensure the xPath agrees with the linked-to data"""
        if self.__relative and hasattr(self, 'toRelativeXLink'):
            self.path = self.toRelativeXLink( self.link )
        elif hasattr(self.link, 'toXLink'):
            self.path = self.link.toXLink()

    def follow( self, startNode ):
        """
        :param startNode: instance corresponding to the beginning of self.path (must inherit from ancestry)
        :return: class instance pointed to by self.path
        """
        return startNode.followXPath( self.path )

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
    def parseXMLNode( cls, linkElement, xPath, linkData ):
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
        relative = not path.startswith('/')
        # all optional (non-xlink) attributes:
        attributes = dict( [v for v in linkElement.items() if xlink_namespace not in v[0]] )
        for key in attributes:
            if key in linkData.get('typeConversion',{}):
                attributes[key] = linkData['typeConversion'][key](attributes[key])
        result = cls( link=None, root=root, path=path, relative=relative, **attributes )
        if( 'unresolvedLinks' in linkData ) : linkData['unresolvedLinks'].append( result )
        xPath.pop()

        return result
