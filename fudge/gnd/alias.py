# <<BEGIN-copyright>>
# <<END-copyright>>

from fudge.core.ancestry import ancestry

__metaclass__ = type

class aliases( ancestry ) :

    moniker = 'aliases'

    def __init__( self, parent = None ) :

        ancestry.__init__( self, aliases.moniker, parent )
        self.aliases = {}

    def __len__( self ) :

        return( len( self.aliases ) )

    def __getitem__( self, key ) :

        return( self.aliases[key] )

    def __setitem__( self, key, value ) :

        if( type( key ) != type( '' ) ) : raise Exception( 'Key must be a string' )
        if( type( value ) != type( '' ) ) : raise Exception( 'Value must be a string' )
        if( key in self.aliases ) : raise Exception( 'Key already in aliases' )
        self.aliases[key] = value

    def toXMLList( self, indent = '' ) :

        indent2 = indent + '  '
        xmlString = []
        if( len( self ) ) :
            xmlString = [ '%s<%s>' % ( indent, self.moniker ) ]
            for key in sorted( self.aliases ) : xmlString.append( '%s<alias key="%s" value="%s"/>' % ( indent2, key, self.aliases[key] ) )
            xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

if( __name__ == '__main__' ) :

    a = aliases( )
    a['Co58_m1'] = 'Co58_e1'
    a['Ag110_m1'] = 'Ag110_e2'
    a['Cd115_m1'] = 'Cd115_e1'
    a['Te127_m1'] = 'Te127_e2'
    a['Te129_m1'] = 'Te129_e1'
    a['Pm148_m1'] = 'Pm148_e2'
    a['Ho166_m1'] = 'Ho166_e1'
    a['Am242_m1'] = 'Am242_e2'
    a['Am244_m1'] = 'Am244_e1'
    a['Es254_m1'] = 'Es254_e2'
    print '\n'.join( a.toXMLList( ) )
