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

from ..atomic import atomic as atomicDataModule

from . import misc as chemicalElementMiscModule
from . import isotope as isotopeModule

class chemicalElement( miscModule.classWithSymbolKey ) :
    """
    This class represents a chemical element (i.e., Oyxgen, Iron, Uranium). A chemical element
    stores all isotopes that contain the same number of protons (i.e., O16, O17, and O18). 
    The members of a chemicalElement are::

        +-------------------+---------------------------------------------------------------+
        | member            | type          | description                                   |
        +===================|===========+===================================================|
        | symbol            | string        | This is the symbol for the chemical element   |
        |                   |               | (e.g., 'O' for oxygen).                       |
        +-------------------+---------------------------------------------------------------+
        | Z                 | int           | Number of protons for this chemical element.  |
        |                   |               | (e.g., 8 for 'O').                            |
        +-------------------+---------------------------------------------------------------+
        | name              | string        | The name of the chemical element              |
        |                   |               | (i.e., 'Oxygen')                              |
        +-------------------+---------------------------------------------------------------+
        | atomicData        | atomicData    | Container for atomicData                      |
        +-------------------+---------------------------------------------------------------+
        | isotopes          | suite         | The list of isotopes                          |
        +-------------------+---------------------------------------------------------------+
    """

    moniker = 'chemicalElement'

    def __init__( self, symbol, Z, name ) :

        miscModule.classWithSymbolKey.__init__( self, symbol )

        base, anti = miscModule.baseAntiFromID( symbol )

        self.__anti = anti == miscModule.antiSuffix

        chemicalElementMiscModule.checkZ( Z )
        self.__Z = Z
        if( chemicalElementMiscModule.symbolFromZ[Z] != base ) : raise ValueError( 'Z = %s and chemical element symbol = "%s" are not consistent' %
                ( Z, miscModule.toLimitedString( symbol ) ) )

        if( not( isinstance( name, str ) ) ) : TypeError( 'name must be a string' )
        if( chemicalElementMiscModule.nameFromZ[Z] != name ) : 
            raise ValueError( 'Z = "%s" and chemical element name = "%s" are not consistent, should be "%s".' %
                    ( Z, miscModule.toLimitedString( name ), chemicalElementMiscModule.nameFromZ[Z] ) )
        self.__name = name

        self.__atomicData = None

        self.__isotopes = isotopeModule.isotopes( )
        self.__isotopes.setAncestor( self ) 

    def __contains__( self, key ) : 

        return( key in self.__isotopes )

    def __getitem__( self, key ) :

        return( self.__isotopes[key] )

    def __iter__( self ) :

        for isotope in self.__isotopes : yield isotope

    @property
    def Z( self ) :

        return( self.__Z )

    @property
    def isAnti( self ) :

        return( self.__anti )

    @property
    def name( self ) :

        return( self.__name )

    @property
    def atomicData( self ) :

        return( self.__atomicData )

    @atomicData.setter
    def atomicData( self, value ) :

        if( not( isinstance( value, atomicDataModule.atomic ) ) ) : raise TypeError( 'value must be an instance of atomicData.' )
        self.__atomicData = value
        self.__atomicData.setAncestor( self )

    @property
    def isotopes( self ) :

        return( self.__isotopes )

    def add( self, item ) :

        self.__isotopes.add( item )

    def convertUnits( self, unitMap ) :

        if( self.__atomicData is not None ) : self.__atomicData.convertUnits( unitMap )
        self.__isotopes.convertUnits( unitMap )

    def copy( self ) :

        _chemicalElement = chemicalElement( self.symbol, self.Z, self.name )
        if( self.__atomicData is not None ) : _chemicalElement.atomicData = self.__atomicData.copy( )
        for item in self.__isotopes : _chemicalElement.add( item.copy( ) )
        return( _chemicalElement )

    def sortCompare( self, other ) :

        if( not( isinstance( other, chemicalElement ) ) ) : raise TypeError( 'Invalid other.' )
        return( self.Z - other.Z )

    def check(self, info):

        from .. import warning as warningModule
        warnings = []

        for isotope in self:
            isotopeWarnings = isotope.check(info)
            if isotopeWarnings:
                warnings.append(warningModule.context('Isotope %s' % isotope.symbol, isotopeWarnings))

        return warnings

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs )  ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s symbol="%s" Z="%s" name="%s">' % ( indent, self.moniker, self.symbol, self.Z, self.name ) ]
        if( self.__atomicData is not None ) : XMLStringList += self.__atomicData.toXMLList( indent2, **kwargs )
        XMLStringList += self.__isotopes.toXMLList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker
        return( XMLStringList )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( '%s[@symbol="%s"]' % ( element.tag, element.get( 'symbol' ) ) )

        for child in element :
            if( child.tag == isotopeModule.isotopes.moniker ) :
                self.__isotopes.parseXMLNode( child, xPath, linkData )
            elif( child.tag == atomicDataModule.atomic.moniker ) :
                self.atomicData = atomicDataModule.atomic.parseXMLNodeAsClass( child, xPath, linkData )
            else :
                raise ValueError( 'Invalid child = "%s" for %s' % ( child.tag, self.moniker ) )

        xPath.pop( )
        return( self )

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ) :

        xPath.append( '%s[@symbol="%s"]' % ( element.tag, element.get( 'symbol' ) ) )

        self = cls( element.attrib['symbol'], int( element.attrib['Z'] ), element.attrib['name'] )
        xPath.pop( )

        self.parseXMLNode( element, xPath, linkData )

        return( self )

    @classmethod
    def parseXMLStringAsClass( cls, string ) :

        from xml.etree import cElementTree

        element = cElementTree.fromstring( string )

        attributes = {}
        for attributeName in element.attrib :
            if( attributeName == 'symbol' ) :
                symbol = element.attrib[attributeName]
            elif( attributeName == 'Z' ) :
                Z = int( element.attrib[attributeName] )
            elif( attributeName == 'name' ) :
                name = element.attrib[attributeName]
            else :
                raise ValueError( 'Unknown attribute = "%s"' % attributeName )

        return( cls( symbol, Z, name ).parseXMLNode( element, [], [] ) )
