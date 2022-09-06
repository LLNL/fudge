# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module defines the nuclides class, which stores the list of nuclides produced by fission.
Together with yields.
"""

from LUPY import ancestry as ancestryModule

class Nuclides(ancestryModule.AncestryIO):

    moniker = 'nuclides'
        
    def __init__( self, _nuclides ) :
    
        ancestryModule.AncestryIO.__init__(self)

        if( not( isinstance( _nuclides, Nuclides ) ) ) :
            if( not( isinstance( _nuclides, ( list, tuple ) ) ) ) : raise TypeError( 'Nuclides must be a list or tuple.' )
        self.__nuclides = tuple( nuclide for nuclide in _nuclides )

    def __len__( self ) :

        return( len( self.__nuclides ) )

    def __iter__( self ) :
        
        n1 = len( self.__nuclides )
        for i1 in range( n1 ) : yield self.__nuclides[i1]
        
    def __getitem__( self, index ) :

        return( self.__nuclides[index] )

    def __eq__( self, other ) :

        if( not( isinstance( other, ( list, tuple, Nuclides ) ) ) ) : return( False )
        if( len( self ) != len( other ) ) : return( False )
        for i1, nuclide in enumerate( self.__nuclides ) :
            if( nuclide != other[i1] ) : return( False )
        return( True )

    def __ne__( self, other ) :

        return( not( self == other ) )

    @property
    def data( self ) :

        return( self.__nuclides )

    def toXML_strList( self, indent = '', **kwargs ) :
    
        if len(self) == 0:
            return []

        XMLString = '%s<%s> ' % ( indent, self.moniker )
        XMLString += ' '.join( [ nuclide for nuclide in self.__nuclides ] )
        XMLString += '</%s>' % self.moniker

        return( [ XMLString ] )

    def parseNode(self, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)
        self.__nuclides = tuple( element.text.strip().split() )
        xPath.pop()

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)
        self = cls([])
        xPath.pop()
        self.parseNode(element, xPath, linkData, **kwargs)

        return self
