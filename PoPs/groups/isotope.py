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
This module contains the class for the suite isotope.
"""

naturalAID = 1000000
MaxA = 400

from .. import suite as suiteModule
from .. import misc as miscModule

from ..families import nuclearLevel as nuclearLevelModule
from ..families import nucleus as nucleusModule

class suite( suiteModule.sortedSuite ) :
    """
    The class is and isotope that contains a list if nuclear levels. The ids for this class are SA[_anti]
    where S is the istopes's chemical element symbol, A is the isotope's nucleon number.
    """

    moniker = 'isotope'

    def __init__( self, id, A ) :

        from . import chemicalElement as chemicalElementModule

        chemicalElement, _A, anti, qualifier = chemicalElementAndAIDsFromIsotopeID( id )

        if( A != _A ) : raise ValueError( 'id = "%s" does not agree with A = "%s"' % 
                ( miscModule.toLimitedString( id ), miscModule.toLimitedString( A ) ) )

        suiteModule.sortedSuite.__init__( self, [ nuclearLevelModule.particle ], key = 'id', replace = True )

        self.__id = id

        self.__anti = anti == miscModule.antiSuffix

        self.__chemicalElement = chemicalElement
        self.__A = A
        self.__intA = checkA( A )

        self.__Z = chemicalElementModule.ZFromSymbol[chemicalElement]

    def __contains__( self, key ) :

        if( suiteModule.sortedSuite.__contains__( self, key ) ) : return( True )
        if( isinstance( key, str ) ) :
            keyUpper = key[:1].upper( ) + key[1:]
            if( suiteModule.sortedSuite.__contains__( self, keyUpper ) ) : return( True )   # all nuclearLevels must have a nucleus.
        return( False )

    def __getitem__( self, key ) :

        if( isinstance( key, int ) ) : return( suiteModule.sortedSuite.__getitem__( self, key ) )
        try :
            return( suiteModule.sortedSuite.__getitem__( self, key ) )
        except :
            keyUpper = key[:1].upper( ) + key[1:]
            if( suiteModule.sortedSuite.__contains__( self, keyUpper ) ) : return( suiteModule.sortedSuite.__getitem__( self, keyUpper ).nucleus )
        raise KeyError( 'key "%s" not found' % key )

    def __eq__( self, other ) :

        from ..families import nuclearLevel as nuclearLevelModule

        if(   isinstance( other, suite ) ) :
            return( self.id == other.id )
        elif( isinstance( other, nuclearLevelModule.particle ) ) :
            _particle = self.particle
            return( _particle.id == other.id )
        else :
            return( False )

    @property
    def A( self ) :

        return( self.__A )

    @property
    def intA( self ) :

        return( self.__intA )

    @property
    def Z( self ) :

        return( self.__Z )

    @property
    def chemicalElement( self ) :

        return( self.__chemicalElement )

    @property
    def groundState( self ) :

        if( len( self ) > 0 ) :
            groundState = self[0]
            if( groundState.intIndex == 0 ) : return( groundState )
        raise Exception( 'No ground state present for isotope id = "%s".' % self.id )

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
    def mass( self ) :

        return( self.groundState.mass )

    @property
    def spin( self ) :

        return( self.groundState.spin )

    @property
    def parity( self ) :

        return( self.groundState.parity )

    @property
    def charge( self ) :

        return( self.groundState.charge )

    @property
    def halflife( self ) :

        return( self.groundState.halflife )

    @property
    def energy( self ) :

        return( self.groundState.energy )

    @property
    def nucleus( self ) :

        return( self.groundState.nucleus )

    def add( self, particle ) :

        if( not( isinstance( particle, nuclearLevelModule.particle ) ) ) :
            raise TypeError( "Particle not a nuclearLevel object." )

        suiteModule.sortedSuite.add( self, particle )

    def check( self, info ):

        from .. import warning as warningModule
        warnings = []

        emax = -1
        for level in self:
            levelWarnings = level.check(info)
            if levelWarnings:
                warnings.append(warningModule.context('Level %s' % level.index, levelWarnings))

            if not level.index.isdigit(): continue  # 'c' and 'u' levels may be out of order
            enow = level.energy.float('eV')
            if enow <= emax:
                warnings.append(warningModule.discreteLevelsOutOfOrder(level.index, level))
            else:
                emax = enow

        return warnings

    def copy( self ) :

        isotope = suite( self.id, self.A )
        for item in self : isotope.add( item.copy( ) )
        return( isotope )

    def getMass( self, unit ) :

        return( self[0].getMass( unit ) )

    def sortCompare( self, other ) :

        if( not( isinstance( other, suite ) ) ) : raise TypeError( 'Invalid other.' )
        return( self.intA - other.intA )

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs )  ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s id="%s" A="%s">' % ( indent, self.moniker, self.id, self.A ) ]
        for level in self : XMLStringList += level.toXMLList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker
        return( XMLStringList )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( element.tag )
    
        for child in element :
            self.add( nuclearLevelModule.particle.parseXMLNodeAsClass( child,  xPath, linkData ) ) 
    
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
            elif( attributeName == 'A' ) :
                A = element.attrib[attributeName]
            else :
                raise ValueError( 'Unknown attribute = "%s"' % attributeName )

        return( cls( id, A ).parseXMLNode( element, [], [] ) )

def isotopeIDFromElementIDAndA( elementID, A ) :

    checkA( A )
    return( "%s%s" % ( elementID, A ) )

def chemicalElementAndAIDsFromIsotopeID( id, qualifierAllowed = False ) :
    """
    Returns the chemicalElement id and A from an isotope's id.
    """

    from . import chemicalElement as chemicalElementModule

    baseID, anti, qualifier = miscModule.baseAntiQualifierFromID( id, qualifierAllowed = qualifierAllowed )

    chemicalElementID, sep, A = baseID.rpartition( '_' )
    if( sep == '' ) :
        if(   baseID[1:].isdigit( ) ) :
            chemicalElementID, A = baseID[0], baseID[1:]
        elif( baseID[2:].isdigit( ) ) :
            chemicalElementID, A = baseID[:2], baseID[2:]
        elif( baseID[3:].isdigit( ) ) :
            chemicalElementID, A = baseID[:3], baseID[3:]
        else :
            raise ValueError( 'Invalid isotope id.' )
        try :
            Z = chemicalElementModule.ZFromSymbol[chemicalElementID]
        except :
            raise ValueError( 'Invalid chemical symbol "%s".' % chemicalElementID )
        intA = int( A )
        if( intA  != 0 ) :                  # Elementmental isotope.
            if( intA < Z ) : raise ValueError( 'A = %s must be greater than or equal to Z = %s: id = "%s"' % ( A, Z, id ) )
    else :                              # id is chemicalElement + '_natural'.
        if( A == '' ) : raise ValueError( 'Missing A in isotope id.' )
        if( A not in [ 'natural' ] ) : raise ValueError( 'Invalid A = "%s": id = "%s"' ( A, id ) )

    chemicalElementModule.checkSymbol( chemicalElementID )

    return( chemicalElementID, A, anti, qualifier )

def checkA( A ) :

    if( not( isinstance( A, str ) ) ) : raise TypeError( 'A attribute not str' )
    if( A not in ( '_natural', ) ) :
        try :
            intA = int( A )             # In the next line 300 is arbitary.
            if( not( 0 <= intA < MaxA ) ) : raise ValueError( 'integer A out of range: A = %s' % A )
            return( intA )
        except :
            raise ValueError( 'A attribute not str of an integer or "_natural"' )

    return( naturalAID )
