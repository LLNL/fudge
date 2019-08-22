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
This module contains the isotope classes.
"""

from .. import suite as suiteModule
from .. import misc as miscModule

from ..groups import misc as chemicalElementMiscModule

from ..families import nuclide as nuclideModule
from ..families import nucleus as nucleusModule

class isotope( miscModule.classWithSymbolKey ) :
    """
    An isotope that contains a list of nuclides. The symbol is of the form  *SA[_anti]*
    where *S* is the istopes's chemical element symbol and *A* is the isotope's nucleon number.
    """

    moniker = 'isotope'

    def __init__( self, symbol, A ) :

        miscModule.classWithSymbolKey.__init__( self, symbol )

        baseID, chemicalElementSymbol, _A, levelID, isNucleus, anti, qualifier = chemicalElementMiscModule.chemicalElementALevelIDsAndAnti( symbol )
        isotopeSymbol = chemicalElementMiscModule.isotopeSymbolFromChemicalElementIDAndA( chemicalElementSymbol, A )
        if( symbol != isotopeSymbol ) : ValueError( 'invalid isotope symbol = "%s"' % symbol )

        if( A != _A ) : raise ValueError( 'symbol = "%s" does not agree with A = "%s"' % 
                ( miscModule.toLimitedString( symbol ), miscModule.toLimitedString( A ) ) )

        self.__chemicalElementSymbol = chemicalElementSymbol
        self.__A = A

        self.__Z = chemicalElementMiscModule.ZFromSymbol[chemicalElementSymbol]

        self.__atomicData = None

        self.__nuclides = nuclideModule.suite( )
        self.__nuclides.setAncestor( self )

    def __len__( self ) :

        return( len( self.__nuclides ) )

    def __contains__( self, key ) :

        baseID, chemicalElementSymbol, A, levelID, isNucleus, anti, qualifier = chemicalElementMiscModule.chemicalElementALevelIDsAndAnti( key )
        if( None in [ chemicalElementSymbol, A, levelID ] ) : return( False )
        isotopeSymbol = chemicalElementMiscModule.isotopeSymbolFromChemicalElementIDAndA( chemicalElementSymbol, A )
        nuclideID = chemicalElementMiscModule.nuclideIDFromIsotopeSymbolAndIndex( isotopeSymbol, levelID )
        if( nuclideID in self.__nuclides ) :
            if( not( isNucleus ) ) : return( True )
            nuclide = self.__nuclides[nuclideID]
            return( nuclide.nucleus is not None )
        return( False )

    def __getitem__( self, key ) :

        if( isinstance( key, int ) ) : return( self.__nuclides[0] )
        baseID, chemicalElementSymbol, A, levelID, isNucleus, anti, qualifier = chemicalElementMiscModule.chemicalElementALevelIDsAndAnti( key )
        if( None in [ chemicalElementSymbol, A, levelID ] ) : raise KeyError( 'key "%s" not found' % key )
        isotopeSymbol = chemicalElementMiscModule.isotopeSymbolFromChemicalElementIDAndA( chemicalElementSymbol, A )
        nuclideID = chemicalElementMiscModule.nuclideIDFromIsotopeSymbolAndIndex( isotopeSymbol, levelID )
        if( nuclideID in self.__nuclides ) :
            nuclide = self.__nuclides[nuclideID]
            if( not( isNucleus ) ) : return( nuclide )
            if( nuclide.nucleus is not None ) : return( nuclide.nucleus )
        raise KeyError( 'key "%s" not found' % key )

    def __iter__( self ) :

        for nuclide in self.__nuclides : yield nuclide

    @property
    def A( self ) :

        return( self.__A )

    @property
    def Z( self ) :

        return( self.__Z )

    @property
    def chemicalElementSymbol( self ) :

        return( self.__chemicalElementSymbol )

    @property
    def nuclides( self ) :

        return( self.__nuclides )

    def add( self, particle ) :

        self.__nuclides.add( particle )

    def check( self, info ):

        from .. import warning as warningModule
        warnings = []

        emax = -1
        for nuclide in self:
            nuclideWarnings = nuclide.check(info)
            if nuclideWarnings:
                warnings.append(warningModule.context('nuclide %s' % nuclide.id, nuclideWarnings))

            enow = nuclide.nucleus.energy.float('eV')
            if enow <= emax:
                warnings.append(warningModule.discreteLevelsOutOfOrder(nuclide.nucleus.index, nuclide))
            else:
                emax = enow

        return warnings

    def convertUnits( self, unitMap ) :

        if( self.__atomicData is not None ) : self.__atomicData.convertUnits( unitMap )
        self.__nuclides.convertUnits( unitMap )

    def copy( self ) :

        _isotope = isotope( self.symbol, self.A )
        for nuclide in self.__nuclides : _isotope.add( nuclide.copy( ) )
        return( _isotope )

    def sortCompare( self, other ) :

        if( not( isinstance( other, isotope ) ) ) : raise TypeError( 'Invalid other.' )
        return( self.A - other.A )

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs )  ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s symbol="%s" A="%s">' % ( indent, self.moniker, self.symbol, self.A ) ]
        if( self.__atomicData is not None ) : XMLStringList += self.__atomicData.toXMLList( indent = indent2, **kwargs )
        XMLStringList += self.__nuclides.toXMLList( indent = indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker
        return( XMLStringList )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( '%s[@symbol="%s"]' % ( element.tag, element.get( 'symbol' ) ) )
    
        for child in element :
            if( child.tag == nuclideModule.suite.moniker ) :
                self.nuclides.parseXMLNode( child,  xPath, linkData )
            else :
                raise Exception( 'Fix me: tag = %s' % child.tag )
    
        xPath.pop( )
        return( self )

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ) :

        xPath.append( '%s[@symbol="%s"]' % ( element.tag, element.get( 'symbol' ) ) )

        self = cls( element.get( 'symbol' ), int( element.get( 'A' ) ) )
        xPath.pop()
        self.parseXMLNode( element, xPath, linkData )

        return( self )

    @classmethod
    def parseXMLStringAsClass( cls, string ) :

        from xml.etree import cElementTree

        element = cElementTree.fromstring( string )
        self = cls.parseXMLNodeAsClass( element, [], [] )

        return( self )

class isotopes( suiteModule.sortedSuite ) :

    moniker = 'isotopes'

    def __init__( self, replace = True ) :

        suiteModule.sortedSuite.__init__( self, allowedClasses = ( isotope, ), replace = replace )

    def __contains__( self, key ) :

        baseID, chemicalElementSymbol, A, levelID, isNucleus, anti, qualifier = chemicalElementMiscModule.chemicalElementALevelIDsAndAnti( key )
        if( None in [ chemicalElementSymbol, A ] ) : return( False )
        isotopeSymbol = chemicalElementMiscModule.isotopeSymbolFromChemicalElementIDAndA( chemicalElementSymbol, A )
        if( suiteModule.sortedSuite.__contains__( self, isotopeSymbol ) ) :
            return( key in suiteModule.sortedSuite.__getitem__( self, isotopeSymbol ) )
        return( False )

    def __getitem__( self, key ) :

        if( isinstance( key, int ) ) : return( suiteModule.sortedSuite.__getitem__( self, key ) )
        baseID, chemicalElementSymbol, A, levelID, isNucleus, anti, qualifier = chemicalElementMiscModule.chemicalElementALevelIDsAndAnti( key )
        if( None in [ chemicalElementSymbol, A ] ) : raise KeyError( 'key "%s" not found' % key )
        isotopeSymbol = chemicalElementMiscModule.isotopeSymbolFromChemicalElementIDAndA( chemicalElementSymbol, A )
        return( suiteModule.sortedSuite.__getitem__( self, isotopeSymbol )[key] )

    def add( self, item ) :

        if( isinstance( item, isotope ) ) :
            suiteModule.sortedSuite.add( self, item )
        elif( isinstance( item, ( nuclideModule.particle, nucleusModule.particle ) ) ) :
            isotopeSymbol = chemicalElementMiscModule.isotopeSymbolFromChemicalElementIDAndA( item.chemicalElementSymbol, item.A )
            if( not( suiteModule.sortedSuite.__contains__( self, isotopeSymbol ) ) ) : suiteModule.sortedSuite.add( self, isotope( isotopeSymbol, item.A ) )
            suiteModule.sortedSuite.__getitem__( self, isotopeSymbol ).add( item )
        else :
            raise TypeError( "Object not an isotope, nuclide or nucleus particle." )

    def convertUnits( self, unitMap ) :

        for isotope in self : isotope.convertUnits( unitMap )

    def getSymbol( self, key ) :

        return( suiteModule.sortedSuite.__getitem__( self, key ) )

    def hasSymbol( self, key ) :

        return( suiteModule.sortedSuite.__contains__( self, key ) )
