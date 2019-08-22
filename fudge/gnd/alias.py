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
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# <<END-copyright>>

from fudge.core.ancestry import ancestry

__metaclass__ = type

class alias :

    moniker = 'alias'

    def __init__( self, key, value, attributes = {} ) :
        """
        Creates a new alias value and its attributes. The attributes argument must be of instance type dict and all keys and values in
        the dictionary must be an instance of str.
        """

        if( not( isinstance( key, str ) ) ) : raise Exception( 'key must be a string, not type "%s"' % type( key ) )
        if( not( isinstance( value, str ) ) ) : raise Exception( 'Value must be a string, not type "%s"' % type( value ) )
        if( not( isinstance( attributes, dict ) ) ) : raise TypeError( "attributes must be of type dict not '%s'" % type( attributes ) )
        self.key = key
        self.value = value
        self.attributes = {}
        for key, value_ in attributes.items( ) :
            if( not( isinstance( key, str ) ) ) : raise TypeError( "attribute's key must be an instance of str, not '%s'", type( key ) )
            if( not( isinstance( value_, str ) ) ) : raise TypeError( "The value for attribute '%s' must be an instance of str, not '%s'" %
                ( key, type( value_ ) ) )
            self.attributes[key] = value_
        self.attributes = attributes

    def __str__( self ) :

        return( 'alias -> key="%s" %s' % ( self.key, self.XMLData( ) ) )

    def getValue( self ) :
        """Returns the value for self."""

        return( self.value )

    def getAttribute( self, name ) :
        """Returns the value for the attribute name in self."""

        return( self.attributes[name] )

    def getAttributes( self ) :
        """Returns a copy of self's attributes."""

        attributes = {}
        for key, value in self.attributes.items( ) : attributes[key] = value_
        return( attributes )

    def hasAttribute( self, name ) :
        """Returns True if self has an attribute called name and False otherwise."""

        return( name in self.attributes )

    def XMLData( self ) :

        s1 = 'value="%s"' % self.value
        for key, value_ in self.attributes.items( ) :
            s1 += ' %s="%s"' % ( key, value_ )
        return( s1 )

class aliases( ancestry, dict ) :

    moniker = 'aliases'

    def __init__( self, parent = None ) :

        ancestry.__init__( self, aliases.moniker, parent )

    def __setitem__( self, key, value ) :
        "Sets the alias key to value. To set an alias' attributes, use add instead."

        self.add( key, value )

    def add( self, key, value, attributes = {} ) :
        """
        This is like __setitem__, but allows for attributes to be associated with the alias defined by key/value. The argument attributes
        must be a dictionary, and its keys and values must be an instance of str.
        """

        if( not( isinstance( key, str ) ) ) : raise Exception( 'Key must be a string' )
        if( key in self ) : raise Exception( "Key '%s' already in aliases" % key )
        dict.__setitem__( self, key, alias( key, value, attributes ) )

    def addNuclearMetaStable( self, isotopeName, nuclearLevelName, metaStableIndex ) :
        """
        Adds an alias to the list of aliases with the key 'isotopeName_m#' where # is metaStableIndex as an integer, a value 
        nuclearLevelName and attribute 'nuclearMetaStable' with value str( metaStableIndex ).
        """

        try :
            metaStableIndex = int( metaStableIndex )
        except :
            raise TypeError( 'metaStableIndex must be an instance of int or converible to int, it is of type "%s"' % type( metaStableIndex ) )
        self.add( self.nuclearMetaStableName( isotopeName, metaStableIndex ), nuclearLevelName, { 'nuclearMetaStable' : str( metaStableIndex ) } )

    def getAliasesFor( self, value ) :
        """Returns a list of all aliases that have value value."""

        aliases_ = []
        for key, alias in self.items( ) :
            if( alias.getValue( ) == value ) : aliases_.append( key )
        return( aliases_ )

    def toXMLList( self, indent = '' ) :
        "Returns a type list of XML strings representing self in XML."

        indent2 = indent + '  '
        xmlString = []
        if( len( self ) ) :
            xmlString = [ '%s<%s>' % ( indent, self.moniker ) ]
            for key in sorted( self ) : xmlString.append( '%s<alias key="%s" %s/>' % ( indent2, key, self[key].XMLData( ) ) )
            xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    @staticmethod
    def nuclearMetaStableName( isotopeName, metaStableIndex ) :

        return( "%s_m%d" % ( isotopeName, int( metaStableIndex ) ) )

if( __name__ == '__main__' ) :

    a = aliases( )
    a['e-'] = 'e'
    a['electron'] = 'e'
    a['e+'] = 'e_anti'
    a['h1'] = 'p'
    a.addNuclearMetaStable( 'Co58', 'Co58_e1', 1 )
    a.addNuclearMetaStable( 'Ag110', 'Ag110_e2', 1 )
    a.addNuclearMetaStable( 'Cd115', 'Cd115_e1', 1 )
    a.addNuclearMetaStable( 'Te127', 'Te127_e2', 1 )
    a.addNuclearMetaStable( 'Te129', 'Te129_e1', 1 )
    a.addNuclearMetaStable( 'Pm148', 'Pm148_e2', 1 )
    a.addNuclearMetaStable( 'Ho166', 'Ho166_e1', 1 )
    a.addNuclearMetaStable( 'Am242', 'Am242_e2', 1 )
    a.addNuclearMetaStable( 'Am244', 'Am244_e1', 1 )
    a.addNuclearMetaStable( 'Es254', 'Es254_e2', 1 )
    print '\n'.join( a.toXMLList( ) )
    print
    print ' ---- aliases for "e" ---- '
    print a.getAliasesFor( 'e' )
