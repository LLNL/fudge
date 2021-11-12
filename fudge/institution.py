# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
"""
This module defines the 'institution' class, used inside GNDS to store institution-specific data
that does not necessarily conform to GNDS specifications.

Within a reactionSuite (or possibly covarianceSuite), institution data is stored inside an <applicationData> section.
Each institution has a unique label and a list of data types. The only restriction on those data types is that
they be well-formed xml. Codes reading in GNDS files are free to ignore any unrecognized institution or data type.
"""

from xml.etree import ElementTree

from xData import ancestry as ancestryModule
from brownies.legacy.toENDF6 import ENDFconversionFlags as ENDFconversionFlagsModule
from fudge.processing import nuclearPlusCoulombInterference as nuclearPlusCoulombInterferenceModule

__metaclass__ = type

class institution(ancestryModule.ancestry):

    moniker = "institution"

    def __init__(self, label):

        ancestryModule.ancestry.__init__(self)
        self.__label = label
        self.__data = []

    def __getitem__(self, item):
        return self.__data[item]

    @property
    def label(self): return self.__label

    @property
    def data(self): return self.__data

    def append( self, item ) :

        self.__data.append( item )
        item.setAncestor( self )

    def findLinks( self, links ) :

        for item in self :
            if( hasattr( item, 'findLinks' ) ) : getattr( item, 'findLinks' )( links )

    def toXMLList( self, indent='', **kwargs ):

        indent2 = indent + kwargs.get('incrementalIndent','  ')
        xml = ['%s<%s label="%s">' % (indent,self.moniker,self.label)]
        for data in self:
            xml += data.toXMLList(indent2, **kwargs)
        xml[-1] += '</%s>' % self.moniker
        return xml

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):

        xPath.append('%s[@label="%s"]' % (element.tag, element.get('label')))

        if( element.get( 'label' ) != "LLNL" ) :
            _institution = UnknownInstitution( UnknownInstitutionXMLNode( element ) )
        else :
            _institution = cls(element.get('label'))
            for child in element:
                if child.tag == ENDFconversionFlagsModule.ENDFconversionFlags.moniker:
                    _institution.append( ENDFconversionFlagsModule.ENDFconversionFlags.parseXMLNode( child, xPath, linkData ) )
                elif child.tag == nuclearPlusCoulombInterferenceModule.NuclearPlusCoulombInterference.moniker:
                    _institution.append( nuclearPlusCoulombInterferenceModule.NuclearPlusCoulombInterference.parseNodeUsingClass( child, xPath, linkData ) )
                else:
                    print("WARNING: encountered unknown data type '%s' in %s" % (child.tag, element.tag))

        xPath.pop()

        return _institution

class UnknownInstitution( ancestryModule.ancestry ) :

    moniker = 'institution'

    def __init__( self, node ) :

        self.__node = node

    @property
    def label( self ) :
        """This method returng the label of self."""

        return( self.__node.label )

    @property
    def node( self ) :
        """This method returns the node instance of self."""

        return( self.__node )

    @node.setter
    def node( self, a_node ) :
        """This method replaces the current node with a_node. a_node must have methods label and toXMLList."""

        if( hasattr( a_node, 'label' ) and hasattr( a_node, 'toXMLList' ) ) :
            self.__node = a_node
        else :
            raise TypeError( 'Node is missing a label and/or toXMLList method.' )

    def toXMLList( self, indent = '', **kwargs ) :
        """This method returng the XML text that was not parsed."""

        return( self.__node.toXMLList( indent = indent, **kwargs ) )

class UnknownInstitutionXMLNode :

    def __init__( self, element ) :

        self.__label = element.get( 'label' )
        self.__stringList = ElementTree.tostring( element.data, encoding = 'unicode' ).rstrip( ).split( '\n' )

    @property
    def label( self ) :
        """This method returng the XML text that was not parsed."""

        return( self.__label )

    @property
    def stringList( self ) :
        """This method returng the XML text that was not parsed."""

        return( self.__stringList )

    def toXMLList( self, indent = '', **kwargs ) :
        """This method returng the XML text that was not parsed."""

        strings = self.stringList
        if( len( strings ) > 0 ) : strings[0] = indent + strings[0]
        return( self.stringList )

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ) :

        return( cls( element, xPath, linkData ) )
