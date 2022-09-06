# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the class for the suite chemicalElements.
"""

from .. import suite as suiteModule

from ..chemicalElements import misc as chemicalElementMiscModule

from . import isotope as isotopeModule
from . import chemicalElement as chemicalElementModule

from ..families import nuclide as nuclideModule
from ..families import nucleus as nucleusModule

class ChemicalElements( suiteModule.SortedSuite ) :
    """
    This class stores a list of chemicalElement instances.
    """

    moniker = 'chemicalElements'

    def __init__( self ) :
        """
        Creates an empty suite of chemicalElements. Use the add() method to fill the suite.
        """

        suiteModule.SortedSuite.__init__( self, [ chemicalElementModule.ChemicalElement ], replace = True )

    def __contains__( self, key ) :
        """
        Recursively search for key inside self and all child particle groups.
        :param key: symbol or id of a particle or particle group
        :return: True if key found
        """

        baseID, chemicalElementSymbol, A, levelID, isNucleus, anti, qualifier = chemicalElementMiscModule.chemicalElementALevelIDsAndAnti( key )
        if( chemicalElementSymbol is None ) : return( False )
        if( suiteModule.SortedSuite.__contains__( self, chemicalElementSymbol ) ) :
            return( key in suiteModule.SortedSuite.__getitem__( self, chemicalElementSymbol ) )
        return( False )

    def __getitem__( self, key ) :
        """
        :param key: string or integer
        :return: chemicalElement, nuclide or nucleus.  Raises KeyError if key not found
        """
        try:
            baseID, chemicalElementSymbol, A, levelID, isNucleus, anti, qualifier = chemicalElementMiscModule.chemicalElementALevelIDsAndAnti( key )
        except:
            chemicalElementSymbol = None

        if( chemicalElementSymbol is not None ) :
            chemicalElement = suiteModule.SortedSuite.__getitem__( self, chemicalElementSymbol )
            if( baseID == chemicalElementSymbol ) : return( chemicalElement )
            return( chemicalElement[key] )

        return suiteModule.SortedSuite.__getitem__( self, key )

    def add( self, item ) :
        """
        Adds a particle or particle group to the proper location inside self
        :param item: instance of chemicalElement, isotope, nuclide or nucleus
        """

        if( isinstance( item, chemicalElementModule.ChemicalElement ) ) :
            suiteModule.SortedSuite.add( self, item )
        elif( isinstance( item, ( isotopeModule.Isotope, nuclideModule.Particle, nucleusModule.Particle ) ) ) :
            chemicalElementSymbol = item.chemicalElementSymbol
            if( self.hasSymbol( chemicalElementSymbol ) ) :
                suiteModule.SortedSuite.__getitem__( self, chemicalElementSymbol ).add( item )
            else :
                chemicalElement = chemicalElementModule.ChemicalElement( chemicalElementSymbol, item.Z, chemicalElementMiscModule.nameFromZ[item.Z] )
                suiteModule.SortedSuite.add( self, chemicalElement )
                chemicalElement.add( item )
        else :
            raise TypeError( "Particle not a nuclide or nucleus object." )

    def convertUnits( self, unitMap ) :
        """ See convertUnits documentation in PoPs.database """

        for chemicalElement in self : chemicalElement.convertUnits( unitMap )

    def getSymbol( self, key ) :
        """
        Get the chemicalElement with symbol == key, or raises KeyError if no such chemicalElement is found.
        :param key: chemicalElement symbol
        :return: chemicalElement instance (or raises KeyError)
        """

        return( suiteModule.SortedSuite.__getitem__( self, key ) )

    def hasSymbol( self, key ) :
        """
        Search for chemicalElement with symbol == key. Unlike __contains__, this method does not
        recursively search inside child particle groups
        :param key: chemicalElement symbol
        :return: True if key present
        """

        return( suiteModule.SortedSuite.__contains__( self, key ) )

    def parseNode(self, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )
        for child in element :
            self.add(chemicalElementModule.ChemicalElement.parseNodeUsingClass(child, xPath, linkData, **kwargs))
        xPath.pop()

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )
        self = cls( )
        xPath.pop()
        self.parseNode(element, xPath, linkData, **kwargs)

        return( self )
