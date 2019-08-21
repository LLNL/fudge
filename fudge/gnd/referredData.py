# <<BEGIN-copyright>>
# <<END-copyright>>

"""This module contains the referred class."""

from fudge.core.ancestry import ancestry

__metaclass__ = type

monikerReferredData = 'referredData'

def parseXMLNode( element, linkData={} ):
    """ translate <referredData> element from xml. Currently only crossSection elements allowed """
    import fudge
    referred = referredData()
    for dat in element:
        xsc = fudge.gnd.reactionData.crossSection.parseXMLNode( dat[0], linkData )
        referred.appendDatum( xsc )
    return referred

class referredData( ancestry ) :

    def __init__( self, parent = None ) :

        ancestry.__init__( self, monikerReferredData, parent )
        self.data = []

    def __len__( self ) :

        return( len( self.data ) )

    def __getitem__( self, index ) :

        return( self.data[index] )

    def appendDatum( self, data ) :

        label = str( len(self.data) )
        self.data.append( referredDatum( self, label, data ) )
        data.setParent( self.data[-1] )

    def toXMLList( self, indent = "" ) :

        if( len( self ) == 0 ) : return( [] )
        indent2, xmlString = indent + '  ', []
        xmlString.append( '%s<%s>' % ( indent, self.moniker ) )
        for datum in self : xmlString += datum.toXMLList( indent2 )
        xmlString[-1] += '</%s>' % ( self.moniker )
        return( xmlString )

class referredDatum( ancestry ) :

    def __init__( self, parent, label, data ) :

        ancestry.__init__( self, 'key', parent, attribute = "label" )
        self.label = label
        self.data = data
        data.setParent( self )

    def toXMLList( self, indent = "" ) :

        xmlString = [ '%s<key label="%s">' % ( indent, self.label ) ]
        xmlString += self.data.toXMLList( indent + '  ' )
        xmlString[-1] += "</%s>" % self.moniker
        return( xmlString )
