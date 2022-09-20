# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the documentation class.
"""

from LUPY import ancestry as ancestryModule

class Documentation( ancestryModule.AncestryIO ) :
    """For storing descriptive information, either about the reactionSuite or about an individual reaction. """

    moniker = 'documentation'

    def __init__( self, name, documentation ) :

        ancestryModule.AncestryIO.__init__( self )
        self.name = name
        self.documentation = documentation

    def __str__( self ) :

        return( self.documentation )

    def getLines( self ) :

        return( self.documentation.split( '\n' ) )
    
    def toXML_strList( self, indent = '', **kwargs ) :
        
        xmlString = [ '%s<%s name="%s"><![CDATA[' % ( indent, self.moniker, self.name ) ]
        xmlString.append( self.documentation )
        xmlString[-1] += ']]></%s>' % self.moniker
        return( xmlString )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        """ translate <documentation> element from XML: """
        xPath.append( element.tag )
        doc_ = Documentation( element.get('name'), element.text.lstrip('\n') )
        xPath.pop()
        return doc_
