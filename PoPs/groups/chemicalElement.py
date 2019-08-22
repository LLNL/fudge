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

"""
This module contains the class for the suite chemicalElement.
"""

from .. import misc as miscModule
from .. import suite as suiteModule

from ..families import nuclearLevel as nuclearLevelModule
from ..families import nucleus as nucleusModule

from . import isotope as isotopeModule

class suite( suiteModule.sortedSuite ) :
    """
    This class represents a chemical element (i.e., Oyxgen, Iron, Uranium). A chemical element
    stores all isotopes that contain the same number of protons (i.e., O16, O17, and O18). 
    The members of a chemicalElement are::

        +-------------------+-----------------------------------------------------------+
        | member            | type      | description                                   |
        +===================|===========+===============================================|
        | id                | string    | This is the symbol for the chemical element   |
        |                   |           | (e.g., 'O' for oxygen).                       |
        +-------------------+-----------------------------------------------------------+
        | Z                 | int       | Number of protons for this chemical element.  |
        |                   |           | (e.g., 8 for 'O').                            |
        +-------------------+-----------------------------------------------------------+
        | name              | string    | The name of the chemical element              |
        |                   |           | (i.e., 'Oxygen')                              |
        +-------------------+-----------------------------------------------------------+
        | list of isotopes  | suite     | The list of isotopes                          |
        +-------------------+-----------------------------------------------------------+
    """

    moniker = 'chemicalElement'

    def __init__( self, id, Z, name ) :

        suiteModule.sortedSuite.__init__( self, [ isotopeModule.suite ], key = 'id', replace = True )

        base, anti = miscModule.baseAntiFromID( id )

        self.__id = id

        self.__anti = anti == miscModule.antiSuffix

        checkZ( Z )
        self.__Z = Z
        if( symbolFromZ[Z] != base ) : raise ValueError( 'Z = %s and chemical element id = "%s" are not consistent' %
                ( Z, miscModule.toLimitedString( id ) ) )

# FIXME - What is the name if element is an '_anti'?
        if( not( isinstance( name, str ) ) ) : TypeError( 'name must be a string' )
        if( nameFromZ[Z] != name ) : 
            raise ValueError( 'Z = "%s" and chemical element name = "%s" are not consistent, should be "%s".' %
                    ( Z, miscModule.toLimitedString( name ), nameFromZ[Z] ) )
        self.__name = name

    def __contains__( self, key ) : 

        if( suiteModule.sortedSuite.__contains__( self, key ) ) : return( True )
        for chemicalElement in self :
            if( key in chemicalElement ) : return( True )
        return( False )

    def __getitem__( self, key ) :

        def getitem( key ) :

            for item in self :
                if( getattr( item, self.__key ) == key ) : return( item )
            return( None )

        if( isinstance( key, int ) ) : return( self.__items[key] )
        if( not( isinstance( key, str ) ) ) : raise TypeError( 'key must be a string' )

        item = getitem( key )
        if( item is not None ) : return( item )

        chemicalElementID = None
        try :
            isNucleus, chemicalElementID, A, levelIDs, anti, qualifier = nucleusModule.chemicalElementAAndLevelIDsFromNuclearLevelID( key )
        except :
            pass
        if( chemicalElementID is not None ) :
            isotopeID = isotopeModule.isotopeIDFromElementIDAndA( chemicalElementID, A )
            item = getitem( isotopeID )
            if( item is not None ) : return( item[key] )
        raise KeyError( 'key "%s" not found' % key )

    @property
    def Z( self ) :

        return( self.__Z )

    @property
    def id( self ) :

        return( self.__id )

    @property
    def isAnti( self ) :

        return( self.__anti )

    @property
    def key( self ) :

        return( self.__id )

    @key.setter
    def key( self, value ) :

        if( not( isinstance( value, str ) ) ) : raise TypeError( 'id must be a string instance.' )
        self.__id = value

    @property
    def name( self ) :

        return( self.__name )

    def add( self, item ) :

        if( isinstance( item, isotopeModule.suite ) ) :
            suiteModule.sortedSuite.add( self, item )
        elif( isinstance( item, nuclearLevelModule.particle ) ) :
            isotopeID = isotopeModule.isotopeIDFromElementIDAndA( item.chemicalElement, item.A )
            if( isotopeID not in self ) : suiteModule.sortedSuite.add( self, isotopeModule.suite( isotopeID, item.A ) )
            self[isotopeID].add( item )
        else :
            raise TypeError( "Object not an isotope or sub-isotope particle." )

    def copy( self ) :

        chemicalElement = suite( self.id, self.Z, self.name )
        for item in self : chemicalElement.add( item.copy( ) )
        return( chemicalElement )

    def sortCompare( self, other ) :

        if( not( isinstance( other, suite ) ) ) : raise TypeError( 'Invalid other.' )
        return( self.Z - other.Z )

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs )  ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s id="%s" Z="%s" name="%s">' % ( indent, self.moniker, self.id, self.Z, self.name ) ]
        for isotope in self : XMLStringList += isotope.toXMLList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker
        return( XMLStringList )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( element.tag )

        for child in element :
            isotope = isotopeModule.suite( child.attrib['id'], child.attrib['A'] )
            self.add( isotope.parseXMLNode( child,  xPath, linkData ) ) 

        xPath.pop( )
        return( self )

    @classmethod
    def parseXMLStringAsClass( cls, string ) :

        from xml.etree import cElementTree

        element = cElementTree.fromstring( string )

        attributes = {}
        for attributeName in element.attrib :
            if( attributeName == 'id' ) :
                id = element.attrib[attributeName]
            elif( attributeName == 'Z' ) :
                Z = int( element.attrib[attributeName] )
            elif( attributeName == 'name' ) :
                name = element.attrib[attributeName]
            else :
                raise ValueError( 'Unknown attribute = "%s"' % attributeName )

        return( cls( id, Z, name ).parseXMLNode( element, [], [] ) )

elementsZSymbolName = (
    (   1, "H",  "Hydrogen" ),      (   2, "He", "Helium" ),        (   3, "Li", "Lithium" ),
    (   4, "Be", "Beryllium" ),     (   5, "B",  "Boron" ),         (   6, "C",  "Carbon" ),
    (   7, "N",  "Nitrogen" ),      (   8, "O",  "Oxygen" ),        (   9, "F",  "Fluorine" ),
    (  10, "Ne", "Neon" ),          (  11, "Na", "Sodium" ),        (  12, "Mg", "Magnesium" ),
    (  13, "Al", "Aluminium" ),     (  14, "Si", "Silicon" ),       (  15, "P",  "Phosphorus" ),
    (  16, "S",  "Sulphur" ),       (  17, "Cl", "Chlorine" ),      (  18, "Ar", "Argon" ),
    (  19, "K",  "Potassium" ),     (  20, "Ca", "Calcium" ),       (  21, "Sc", "Scandium" ),
    (  22, "Ti", "Titanium" ),      (  23, "V",  "Vanadium" ),      (  24, "Cr", "Chromium" ),
    (  25, "Mn", "Manganese" ),     (  26, "Fe", "Iron" ),          (  27, "Co", "Cobalt" ),
    (  28, "Ni", "Nickel" ),        (  29, "Cu", "Copper" ),        (  30, "Zn", "Zinc" ),
    (  31, "Ga", "Gallium" ),       (  32, "Ge", "Germanium" ),     (  33, "As", "Arsenic" ),
    (  34, "Se", "Selenium" ),      (  35, "Br", "Bromine" ),       (  36, "Kr", "Krypton" ),
    (  37, "Rb", "Rubidium" ),      (  38, "Sr", "Strontium" ),     (  39, "Y",  "Yttrium" ),
    (  40, "Zr", "Zirconium" ),     (  41, "Nb", "Niobium" ),       (  42, "Mo", "Molybdenum" ),
    (  43, "Tc", "Technetium" ),    (  44, "Ru", "Ruthenium" ),     (  45, "Rh", "Rhodium" ),
    (  46, "Pd", "Palladium" ),     (  47, "Ag", "Silver" ),        (  48, "Cd", "Cadmium" ),
    (  49, "In", "Indium" ),        (  50, "Sn", "Tin" ),           (  51, "Sb", "Antimony" ),
    (  52, "Te", "Tellurium" ),     (  53, "I",  "Iodine" ),        (  54, "Xe", "Xenon" ),
    (  55, "Cs", "Cesium" ),        (  56, "Ba", "Barium" ),        (  57, "La", "Lanthanum" ),
    (  58, "Ce", "Cerium" ),        (  59, "Pr", "Praseodymium" ),  (  60, "Nd", "Neodymium" ),
    (  61, "Pm", "Promethium" ),    (  62, "Sm", "Samarium" ),      (  63, "Eu", "Europium" ),
    (  64, "Gd", "Gadolinium" ),    (  65, "Tb", "Terbium" ),       (  66, "Dy", "Dysprosium" ),
    (  67, "Ho", "Holmium" ),       (  68, "Er", "Erbium" ),        (  69, "Tm", "Thulium" ),
    (  70, "Yb", "Ytterbium" ),     (  71, "Lu", "Lutetium" ),      (  72, "Hf", "Hafnium" ),
    (  73, "Ta", "Tantalum" ),      (  74, "W",  "Tungsten" ),      (  75, "Re", "Rhenium" ),
    (  76, "Os", "Osmium" ),        (  77, "Ir", "Iridium" ),       (  78, "Pt", "Platinum" ),
    (  79, "Au", "Gold" ),          (  80, "Hg", "Mercury" ),       (  81, "Tl", "Thallium" ),
    (  82, "Pb", "Lead" ),          (  83, "Bi", "Bismuth" ),       (  84, "Po", "Polonium" ),
    (  85, "At", "Astatine" ),      (  86, "Rn", "Radon" ),         (  87, "Fr", "Francium" ),
    (  88, "Ra", "Radium" ),        (  89, "Ac", "Actinium" ),      (  90, "Th", "Thorium" ),
    (  91, "Pa", "Protactinium" ),  (  92, "U",  "Uranium" ),       (  93, "Np", "Neptunium" ),
    (  94, "Pu", "Plutonium" ),     (  95, "Am", "Americium" ),     (  96, "Cm", "Curium" ),
    (  97, "Bk", "Berkelium" ),     (  98, "Cf", "Californium" ),   (  99, "Es", "Einsteinium" ),
    ( 100, "Fm", "Fermium" ),       ( 101, "Md", "Mendelevium" ),   ( 102, "No", "Nobelium" ),
    ( 103, "Lr", "Lawrencium" ),    ( 104, "Rf", "Rutherfordium" ), ( 105, "Db", "Dubnium" ),
    ( 106, "Sg", "Seaborgium" ),    ( 107, "Bh", "Bohrium" ),       ( 108, "Hs", "Hassium" ),
    ( 109, "Mt", "Meitnerium" ),    ( 110, "Ds", "Darmstadtium" ),  ( 111, "Rg", "Roentgenium" ),
    ( 112, "Cn", "Copernicium" ),   ( 113, "Uut", "Ununtrium" ),    ( 114, "Fl", "Flerovium" ),
    ( 115, "Uup", "Ununpentium" ),  ( 116, "Lv", "Livermorium" ),   ( 117, "Uus", "Ununseptium" ),
    ( 118, "Uuo", "Ununoctium" ) )

symbolFromZ = {}
for Z, symbol, name in elementsZSymbolName : symbolFromZ[Z] = symbol

nameFromZ = {}
for Z, symbol, name in elementsZSymbolName : nameFromZ[Z] = name

ZFromSymbol = {}
for Z, symbol, name in elementsZSymbolName : ZFromSymbol[symbol] = Z

nameFromSymbol = {}
for Z, symbol, name in elementsZSymbolName : nameFromSymbol[symbol] = name

ZFromName = {}
for Z, symbol, name in elementsZSymbolName : ZFromName[name] = Z

symbolFromName = {}
for Z, symbol, name in elementsZSymbolName : symbolFromName[name] = symbol

def checkZ( Z ) :

    if( not( isinstance( Z, int ) ) ) : raise TypeError( 'Z not an int object.' )
    if( 0 < Z <= elementsZSymbolName[-1][0] ) : return( Z )
    raise ValueError( 'Z = "%s" out of range.' % Z )

def checkSymbol( symbol ) :

    if( not( isinstance( symbol, str ) ) ) : raise TypeError( 'symbol not a str object.' )
    if( symbol not in ZFromSymbol ) : raise ValueError( 'Invalid symbol: %s.' % miscModule.toLimitedString( symbol ) )
