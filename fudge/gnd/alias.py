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

import xData.ancestry as ancestryModule

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
        for key, value in self.attributes.items( ) : attributes[key] = value
        return( attributes )

    def hasAttribute( self, name ) :
        """Returns True if self has an attribute called name and False otherwise."""

        return( name in self.attributes )

    def XMLData( self ) :

        s1 = 'value="%s"' % self.value
        for key, value_ in self.attributes.items( ) :
            s1 += ' %s="%s"' % ( key, value_ )
        return( s1 )

class aliases( ancestryModule.ancestry, dict ) :

    moniker = 'aliases'

    def __init__( self ) :

        ancestryModule.ancestry.__init__( self )

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

    def toXMLList( self, indent = '', **kwargs ) :
        "Returns a type list of XML strings representing self in XML."

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

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
