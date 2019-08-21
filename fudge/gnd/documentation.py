# <<BEGIN-copyright>>
# <<END-copyright>>

"""
This module contains the documentation class.
"""

__metaclass__ = type

def parseXMLNode(docElement):
    """ translate XML <documentation> element: """
    return documentation( docElement.get('name'), docElement.text.lstrip('\n') )

class documentation :
    """For storing descriptive information, either about the reactionSuite or about an individual reaction. """

    def __init__( self, name, documentation ) :

        self.name = name
        self.documentation = documentation

    def __str__( self ) :

        return( self.documentation )

    def getLines( self ) :

        return( self.documentation.split( '\n' ) )
    
    def toXMLList( self, indent = '' ) :
        
        xmlString = [ '%s<documentation name="%s"><![CDATA[' % 
                ( indent, self.name ) ]
        xmlString.append( self.documentation )
        xmlString[-1] += ']]></documentation>'
        return( xmlString )
