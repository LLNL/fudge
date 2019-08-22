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
from . import ancestry as ancestryModule

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
    <element href='/reactionSuite/reaction[@label="102"]'>
        or
    <element href='/covariances.xml#covId'>
    """

    moniker = 'link'
    ancestryMembers = ( '', )

    def __init__( self, link = None, root = None, path = None, label = None, relative = False, **attributes ) :

        baseModule.xDataCoreMembers.__init__( self, self.moniker, index = None, label = label )
        self.link = link            # Pointer to another instance.
        self.root = root            # File name where other data is stored, only needed if it's not the current file.
        self.path = path            # Path is an xPath expression that uniquely identifies the other instance.
        self.__relative = relative  # Whether to use relative link when writing out xPath.

            # careful with nesting single/double quotes.
        if self.path is not None: self.path = self.path.replace('"',"'")
        self.attributes = attributes

    def __str__( self ) :

        return( "%s" % self.path )

    def __getitem__( self, key ):
        return self.attributes[key]

    def __setitem__( self, key, value ):
        self.attributes[key] = value

    def copy( self ):
        import copy

        if( self.path == None ) : self.updateXPath( )
        _link = link( link = None, root = self.root, path = self.path, label = self.label,
            relative = self.__relative, **copy.copy(self.attributes) )
        return( _link )

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
        i.e., <element href="...xpath" anotherAttribute="foo" ...
        """

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        self.updateXPath( )

        attributesStr = self.attributesToXMLAttributeStr( )
        for key in sorted( self.attributes.keys( ) ) :
            attributesStr += ' %s="%s"' % ( key, self.attributes[key] )

        if( type( self.root ) == str ) :    # external link
            XMLList = '%s<%s%s href="%s#%s"/>' % ( indent, self.moniker, attributesStr, self.root, self.path )
        else:                               # link within this file
            XMLList = '%s<%s%s href="%s"/>' % ( indent, self.moniker, attributesStr, self.path )

        return( [ XMLList ] )

    @classmethod
    def parseXMLNode( cls, linkElement, xPath, linkData ):
        """
        Parse the xml-represented link back to python. The resulting link points to None,
        and must be resolved by the calling function.
        Link attributes can be converted to desired type by passing a 'typeConversion' dictionary in linkData
        """

        xPath.append( linkElement.tag )
        path = linkElement.get('href')
        if '#' in path: root, path = path.split('#')
        else: root = None
        relative = not path.startswith('/')
        # all optional (non-xlink) attributes:
        attributes = dict( [v for v in linkElement.items() if not v[0].startswith('href')] )
        for key in attributes:
            if key in linkData.get('typeConversion',{}):
                attributes[key] = linkData['typeConversion'][key](attributes[key])
        result = cls( link=None, root=root, path=path, relative=relative, **attributes )
        if( 'unresolvedLinks' in linkData ) : linkData['unresolvedLinks'].append( result )
        xPath.pop()

        return result

class link2( ancestryModule.ancestry ) :
    """
    This class contains an href to another node as a xPath. The href could reference a node in a different file.
    The 'follow' member function is used to get a pointer to the linked-to node,
    which will then be stored as self.instance.
    
    Member data ::
    
        * moniker   : the moniker of the instance pointed to by instance (type: str)
        * instance  : python reference to the linked instance (type: reference)
        * root      : file name where other node exists, only needed if it's not the current file (type: str)
        * href      : an xPath expression that uniquely identifies the other instance (type: str)
        * relative  : Whether to use relative xlink when writing out xPath (type: bool)

    sample links::
    <element href='/reactionSuite/reaction[@label="102"]'/>
        or
    <element href='/covariances.xml#covId'/>
    """

    moniker = 'link'
    ancestryMembers = ( '', )

    def __init__( self, moniker, instance = None, root = None, href = None, relative = False, keyName = None, keyValue = None, label = None ) :

        self.__moniker = moniker
        ancestryModule.ancestry.__init__( self )

        self.__instance = instance      # Pointer to linked instance.
        self.__root = root              # File name where other data is stored, only needed if it's not the current file
        self.__href = href              # xPath expression to linked instance.
        self.__relative = relative      # Whether to use relative xlink when writing out xPath
        self.__keyName = None
        self.__keyValue = None
        self.updateKey( keyName = keyName, keyValue = keyValue )
        self.__label = label

# BRB is this what we want or something else.
        if( self.__href is not None ) : self.__href = self.__href.replace( '"', "'" )        # Careful with nesting single/double quotes.

    def __str__( self ) :

        return( "%s" % self.__href )

    def __getattr__( self, name ) :

        return( getattr( self.__instance, name ) )

    @property
    def moniker( self ) :

        return( self.__moniker )

    @property
    def instance( self ) :

        return( self.__instance )

    @instance.setter
    def instance( self, value ) :

        self.__instance = value

    @property
    def link( self ) :

        return( self.__instance )

    @link.setter
    def link( self, value ) :

        self.instance = value

    @property
    def keyName( self ) :

        return( self.__keyName )

    @property
    def keyValue( self ) :

        return( self.__keyValue )

    @property
    def label( self ) :

        return( self.__label )

    @property
    def root( self ) :

        return( self.__root )

    @property
    def href( self ) :

        return( self.__href )

    @property
    def path( self ) :

        return( self.__href )

    @property
    def relative( self ) :

        return( self.__relative )

    def copy( self, *args, **kwargs ) :

        try :
            if(self.__href is None) : self.updateXPath()
        except :
            pass
        return( link2( self.__moniker, instance = self.__instance, root = self.__root, href = self.__href, keyName = self.__keyName, 
                keyValue = self.__keyValue, label = self.__label ) )

    def updateKey( self, keyName = None, keyValue = None ) :

        if( self.__keyName is not None ) :
            if( self.__keyName != keyName ) : raise Exception( '''self's Keyname "%s" != "%s"''' % ( self.__keyName, keyName ) )

        if( self.__instance is not None ) :
            if( hasattr( self.__instance, 'keyName' ) ) :
                if( keyName is None ) :
                    keyName = self.instance.keyName
                else :
                    if( self.instance.keyName != keyName ) : raise Exception( '''Keyname "%s" != "%s"''' % ( self.instance.keyName, keyName ) )

        self.__keyName = keyName
        self.__keyValue = keyValue

    def updateXPath( self ) :
        """
        Ensure the xPath agrees with the linked instance.
        """

        if( self.__relative ) :
            self.__href = self.toRelativeXLink( self.__instance )
        else :
            if( self.__instance is not None ) :         # FIXME
                self.__href = self.__instance.toXLink( )

    def follow( self, startNode ) :
        """
        :param startNode: instance corresponding to the beginning of self.href (must inherit from ancestry)
        :return: instance pointed to by self.href
        """

        return( startNode.followXPath( self.__href ) )

    def toXML( self, indent = '' , **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :
        """
        Pointers show up in the attributes list on an xml element (i.e., <element href="...xpath" anotherAttribute="foo" .../>).
        """

        self.updateXPath( )

        attributes = ''
        if( ( self.__keyName is not None ) and ( self.__keyValue is not None ) ) : attributes += ' %s="%s"' % ( self.__keyName, self.__keyValue )
        if( self.__label is not None ) : attributes += ' label="%s"' % self.__label

        root = ''
        if( self.__root is not None ) : root += '%s#' % self.__root              # external file

        return( [ '%s<%s%s href="%s%s"/>' % ( indent, self.__moniker, attributes, root, self.__href ) ] )

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ) :
        """
        Parse the xml-represented link back to python. The resulting link points to None, and must be resolved by the calling function.
        Link attributes can be converted to desired type by passing a 'typeConversion' dictionary in linkData
        """

        xPath.append( element.tag )

        root = None
        href = None
        keyName = None
        keyValue = None
        label = None

        for key in element.keys( ) :
            if( key == 'href' ) :
                href = element.get( key )
            elif( key == 'label' ) :
                label = element.get( key )
            else :
                if( keyName is not None ) : raise Exception( 'Too many attributes for link: %s' % key )
                keyName = key
                keyValue = element.get( key )

        if( '#' in href ) : root, href = href.split( '#' )

        _link = cls( element.tag, instance = None, root = root, href = href, relative = href[0] == '/', keyName = keyName, keyValue = keyValue,
            label = label )

        if( 'unresolvedLinks' in linkData ) : linkData['unresolvedLinks'].append( _link )

        xPath.pop( )
        return( _link )
