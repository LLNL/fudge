# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the documentation class.
"""

__metaclass__ = type

import xData.ancestry as ancestryModule

class documentation( ancestryModule.ancestry ) :
    """For storing descriptive information, either about the reactionSuite or about an individual reaction. """

    moniker = 'documentation'

    def __init__( self, name, documentation ) :

        ancestryModule.ancestry.__init__( self )
        self.name = name
        self.documentation = documentation

    def __str__( self ) :

        return( self.documentation )

    def getLines( self ) :

        return( self.documentation.split( '\n' ) )
    
    def toXMLList( self, indent = '', **kwargs ) :
        
        xmlString = [ '%s<%s name="%s"><![CDATA[' % ( indent, self.moniker, self.name ) ]
        xmlString.append( self.documentation )
        xmlString[-1] += ']]></%s>' % self.moniker
        return( xmlString )

    @staticmethod
    def parseXMLNode(element, xPath, linkData):
        """ translate <documentation> element from XML: """
        xPath.append( element.tag )
        doc_ = documentation( element.get('name'), element.text.lstrip('\n') )
        xPath.pop()
        return doc_
