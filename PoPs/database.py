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
This module contains the class for the database class.
"""

import os

from xData import ancestry as ancestryModule

from .groups import chemicalElements as chemicalElementsModule
from .groups import chemicalElement as chemicalElementModule
from .groups import isotope as isotopeModule
from .families import particle as particleModule
from .families import gaugeBoson as gaugeBosonModule
from .families import lepton as leptonModule
from .families import baryon as baryonModule
from .families import nuclide as nuclideModule
from .families import unorthodox as unorthodoxModule

from . import styles as stylesModule
from . import alias as aliasModule

format = '0.1'

class database( ancestryModule.ancestry ) :

    moniker = 'PoPs'

    def __init__( self, name, version, documentation = '' ) :

        ancestryModule.ancestry.__init__( self )

        self.__format = format

        if( not( isinstance( name, str ) ) ) : raise TypeError( 'name must be a str' )
        self.name = name

        if( not( isinstance( version, str ) ) ) : raise TypeError( 'version must be a str' )
        self.version = version

        self.__styles = stylesModule.styles( )
        self.__styles.setAncestor( self )

        self.__gaugeBosons = gaugeBosonModule.suite( )
        self.__gaugeBosons.setAncestor( self )

        self.__leptons = leptonModule.suite( )
        self.__leptons.setAncestor( self )

        self.__baryons = baryonModule.suite( )
        self.__baryons.setAncestor( self )

        self.__chemicalElements = chemicalElementsModule.chemicalElements( )
        self.__chemicalElements.setAncestor( self )

        self.__unorthodoxes = unorthodoxModule.suite( )
        self.__unorthodoxes.setAncestor( self )

        self.__aliases = aliasModule.suite( )
        self.__aliases.setAncestor( self )

        self.documentation = documentation

    def __contains__( self, key ) :

        if( key in self.__gaugeBosons ) : return( True )
        if( key in self.__leptons ) : return( True )
        if( key in self.__baryons ) : return( True )
        if( key in self.__chemicalElements ) : return( True )
        if( key in self.__unorthodoxes ) : return( True )
        if( key in self.__aliases ) : return( True )
        return( False )

    def __getitem__( self, key ) :

        if( '{' in key ) : key = key.split( '{' )[0]
        return( self.__getitem( key, 0 ) )

    def __getitem( self, key, level ) :

        if( not( isinstance( key, str ) ) ) : raise TypeError( 'key must be a str' )
        if( key in self.__aliases ) :
            alias = self.__aliases[key]
            item = self.__getitem( alias.pid, level + 1 )
            if( level == 0 ) : item = item.alias( key, item )
            return( item )
        for suite in [ self.__gaugeBosons, self.__leptons, self.__baryons, self.__chemicalElements, self.__unorthodoxes ] :
            try :
                return( suite[key] )
            except KeyError :
                continue
            except :
                raise
        raise KeyError( 'id = "%s" not found' % key )

    def __iter__( self ) :

        for gaugeBoson in self.__gaugeBosons : yield gaugeBoson
        for lepton in self.__leptons : yield lepton
        for baryon in self.__baryons : yield baryon
        for chemicalElement in self.__chemicalElements :
            for isotope in chemicalElement :
                for nuclide in isotope :
                    yield nuclide
                    yield nuclide.nucleus
        for item in self.__unorthodoxes : yield item
        for item in self.__aliases : yield item

    def keys( self ) :

        return( [ particle.id for particle in self ] )

    @property
    def format( self ) :

        return( self.__format )

    @property
    def name( self ) :

        return( self.__name )

    @name.setter
    def name( self, value ) :

        if( not( isinstance( value, str ) ) ) : raise TypeError( 'name must be a string' )
        self.__name = value

    @property
    def version( self ) :

        return( self.__version )

    @version.setter
    def version( self, value ) :

        if( not( isinstance( value, str ) ) ) : raise TypeError( 'version must be a string' )
        self.__version = value

    @property
    def styles( self ) :

        return( self.__styles )

    @property
    def gaugeBosons( self ) :

        return( self.__gaugeBosons )

    @property
    def leptons( self ) :

        return( self.__leptons )

    @property
    def baryons( self ) :

        return( self.__baryons )

    @property
    def chemicalElements( self ) :

        return( self.__chemicalElements )

    @property
    def unorthodoxes( self ) :

        return( self.__unorthodoxes )

    @property
    def aliases( self ) :

        return( self.__aliases )

    @property
    def documentation( self ):

        return( self.__documentation )

    @documentation.setter
    def documentation( self, value ) :

        if( not( isinstance( value, str ) ) ) : raise TypeError( 'documentation must be a string' )
        self.__documentation = value

    def add( self, particle ) :

        if( isinstance( particle, gaugeBosonModule.particle ) ) :
            self.__gaugeBosons.add( particle )
        elif( isinstance( particle, leptonModule.particle ) ) :
            self.__leptons.add( particle )
        elif( isinstance( particle, baryonModule.particle ) ) :
            self.__baryons.add( particle )
        elif( isinstance( particle, ( chemicalElementModule.chemicalElement, isotopeModule.isotope, nuclideModule.particle ) ) ) :
            self.__chemicalElements.add( particle )
        elif( isinstance( particle, unorthodoxModule.particle ) ) :
            self.__unorthodoxes.add( particle )
        elif( isinstance( particle, aliasModule.alias ) ) :
            if( particle.id not in self.__aliases ) :      # Do not allow an alias id to overwrite a particle id.
                if( particle.id in self ) : raise ValueError( 'particle with id = "%s" already in database' % particle.id )
            self.__aliases.add( particle )
        else :
            raise TypeError( 'Unsupported particle family = "%s"' % particle.moniker )

    def check( self, **kwargs ) :
        from . import warning as warningModule
        warnings = []
        info = kwargs
        info['PoPs'] = self

        for suite in (
                self.__gaugeBosons, self.__leptons, self.__baryons,
                self.__chemicalElements, self.__unorthodoxes, self.__aliases):
            suiteWarnings = suite.check( info )
            if suiteWarnings:
                warnings.append( warningModule.context("%s" % suite.moniker, suiteWarnings) )

        result = warningModule.context('PoPs', warnings)
        result.info = info
        return result

    def convertUnits( self, unitMap ) :

        self.styles.convertUnits( unitMap )
        self.__gaugeBosons.convertUnits( unitMap )
        self.__leptons.convertUnits( unitMap )
        self.__baryons.convertUnits( unitMap )
        self.__chemicalElements.convertUnits( unitMap )
        self.__unorthodoxes.convertUnits( unitMap )

    def copy( self ) :

        _database = database( self.name, self.version, documentation = self.documentation )
        for item in self.__styles : _database.styles.add( item.copy( ) )
# FIXME, missing documentation.
        for item in self.__gaugeBosons : _database.add( item.copy( ) )
        for item in self.__leptons : _database.add( item.copy( ) )
        for item in self.__baryons : _database.add( item.copy( ) )
        for item in self.__chemicalElements : _database.add( item.copy( ) )
        for item in self.__unorthodoxes : _database.add( item.copy( ) )
        for item in self.__aliases : _database.add( item.copy( ) )

        return( _database )

    def hasAlias( self, id ) :

        return( id in self.__aliases )

    def saveToFile( self, fileName, **kwargs ) :

        dirname = os.path.dirname( fileName )
        if( dirname != '' ) :
            if( not( os.path.exists( dirname ) ) ) : os.makedirs( dirname )
        fOut = open( fileName, 'w' )
        fOut.write( '<?xml version="1.0"?>\n' )
        fOut.write( self.toXML( ) )
        fOut.write( '\n' )
        fOut.close( )

    def toXML( self, indent = "", **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ '%s<%s name="%s" version="%s" format="%s">' % ( indent, self.moniker, self.name, self.version, self.format ) ]
        xmlString += self.styles.toXMLList( indent2, **kwargs )
        if( self.documentation != '' ) :
            xmlString.append( '%s<documentations>\n%s<documentation name="endfDoc"><![CDATA[\n%s]]></documentation></documentations>' % ( indent2, indent2+indent2, self.documentation ) )

        xmlString += self.__aliases.toXMLList( indent2, **kwargs )
        xmlString += self.__gaugeBosons.toXMLList( indent2, **kwargs )
        xmlString += self.__leptons.toXMLList( indent2, **kwargs )
        xmlString += self.__baryons.toXMLList( indent2, **kwargs )
        xmlString += self.__chemicalElements.toXMLList( indent2, **kwargs )
        xmlString += self.__unorthodoxes.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( element.tag )

        children = {}
        for child in ( gaugeBosonModule, leptonModule, baryonModule, unorthodoxModule, aliasModule ) :
            children[child.suite.moniker] = child.suite

        children[chemicalElementsModule.chemicalElements.moniker] = chemicalElementsModule.chemicalElements

        for child in element :
            if( child.tag in children ) :
                children[child.tag].parseXMLNode( getattr( self, child.tag ), child, xPath, linkData )
            elif( child.tag == stylesModule.styles.moniker ) :
                self.styles.parseXMLNode( child, xPath, linkData )
            elif( child.tag == 'documentations' ) :
                self.documentation = child.find('documentation').text[1:]
            else :
                raise ValueError( 'sub-element = "%s" not allowed' % child.tag )

        xPath.pop( )
        return( self )

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ) :

        _format = element.attrib['format']
        if( format != _format ) : ValueError( 'Unsupported format = "%s"' % _format )

        self = cls( element.attrib['name'], element.attrib['version'] )
        self.parseXMLNode( element, xPath, linkData )
        return( self )

    @classmethod
    def parseXMLStringAsClass( cls, string ) :

        from xml.etree import cElementTree

        xPath = []
        try :
            element = cElementTree.fromstring( string )
            return( cls.parseXMLNodeAsClass( element, xPath, {} ) )
        except Exception :
            print( "Error encountered at xpath = /%s" % '/'.join( xPath ) )
            raise

    @staticmethod
    def readFile( fileName ) :

        from xml.etree import cElementTree

        xPath = []
        try :
            element = cElementTree.parse( fileName ).getroot( )
            return( database.parseXMLNodeAsClass( element, xPath, {} ) )
        except :
            print( "Error encountered at xpath = /%s" % '/'.join( xPath ) )
            raise
