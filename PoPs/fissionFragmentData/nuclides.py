# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module defines the nuclides class, which stores the list of nuclides produced by fission.
Meant to be used along with the fissionFragmentData.yields.Yields class to store fission product yields.
"""

from LUPY import ancestry as ancestryModule

class Nuclides(ancestryModule.AncestryIO):
    """Stores a list of PoPs particle ids. The order (and length) must match the values in the
    corresponding fissionFragmentData.yields.Yields.values."""

    moniker = 'nuclides'
        
    def __init__(self, _nuclides):
        """Constructs a Nuclides container from a list of product ids.

        :param _nuclides: list or tuple of string PoPs particle ids."""
    
        ancestryModule.AncestryIO.__init__(self)

        # FIXME this test seems wrong/unnecessary, copy constructor should be different method...
        if( not( isinstance( _nuclides, Nuclides ) ) ) :
            if( not( isinstance( _nuclides, ( list, tuple ) ) ) ) : raise TypeError( 'Nuclides must be a list or tuple.' )
        self.__nuclides = tuple( nuclide for nuclide in _nuclides )

    def __len__(self):
        """Returns the number of nuclides."""

        return len(self.__nuclides)

    def __iter__(self):
        
        # CMM why not just 'for nuclide in self.__nuclides: yield nuclide'?
        n1 = len( self.__nuclides )
        for i1 in range( n1 ) : yield self.__nuclides[i1]
        
    def __getitem__(self, index):

        return self.__nuclides[index]

    def __eq__(self, other):
        """Returns True iff self and other contain the same list of nuclides in the same order."""

        if not isinstance(other, (list, tuple, Nuclides)): return False
        if len(self) != len(other): return False
        for i1, nuclide in enumerate(self.__nuclides):
            if nuclide != other[i1]: return False
        return True

    def __ne__( self, other ) :
        # FIXME CMM: I don't think this is necessary if __eq__ is defined

        return( not( self == other ) )

    @property
    def nuclides(self):
        """Accesses the list of nuclides."""

        return self.__nuclides

    @property
    def data(self):
        """Alternate method for accessing the list of nuclides."""

        return self.nuclides

    def toXML_strList(self, indent = '', **kwargs):
        """
        Returns a list of str instances representing the XML lines of self.

        :param indent:    The amount of indentation for each line. Child nodes and text may be indented more.
        :param kwargs:    A keyword list.

        :return:          List of str instances representing the XML lines of self.
        """
    
        if len(self) == 0:
            return []

        XMLString = '%s<%s> ' % ( indent, self.moniker )
        XMLString += ' '.join( [ nuclide for nuclide in self.__nuclides ] )
        XMLString += '</%s>' % self.moniker

        return [ XMLString ]

    def parseNode(self, node, xPath, linkData, **kwargs):

        xPath.append(node.tag)
        self.__nuclides = tuple( node.text.strip().split() )
        xPath.pop()

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """Create a new instance of class **cls** and parse contents of node into the instance.

        :param cls: FUDGE Python class to return.
        :param node: **nuclides** node to parse.
        :param xPath: List containing xPath to current node, useful mostly for debugging.
        :param linkData: dict that collects unresolved links.
        """

        xPath.append(node.tag)
        self = cls([])
        xPath.pop()
        self.parseNode(node, xPath, linkData, **kwargs)

        return self
