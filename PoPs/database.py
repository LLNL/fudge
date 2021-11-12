# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module defines the PoPs database class, which is the top-level container in PoPs.
A database is a collection of particles, each with a unique id.
It may also store documentation and 'style' information.

Within the database, particles are organized by family (i.e. gauge bosons, leptons, baryons or chemical elements).
For nuclear data users, the largest part of the database will likely be inside the chemical elements section
where isotopes, nuclides and excited nuclear states are defined.

The database supports looking up particles by their id. For example, if 'pops' is a PoPs.database instance that
contains a particle with id 'Al27',

>>> 'Al27' in pops
True
>>> type( pops['Al27'] )
<class PoPs.families.nuclide.particle>
>>> pops.keys()   # list all particle ids in the database
['photon', 'n', ... 'Al28', 'Al26_m1', 'Na24_m1']

The particle database may also define aliases. These can be used to refer to a particle by another name.
For example, the id for the electron is normally "e-", but could also be aliased to 'e' or 'electron'.
In the example above, 'Al26_m1' and 'Na24_m1' are aliases pointing to excited nuclides.

>>> type( pops['Na24_m1'] )
<class 'PoPs.families.nuclide.alias'>

Each alias has a 'pid' field with the id of the actual particle it points to:

>>> pops['Na24_m1'].pid
'Na24_e1'

In this case, the first metastable state of Na24 is also the first excited state (as indicated by the '_e1' suffix).

A PoPs database has methods to write itself out to an XML-formatted file or string, to read back in from XML,
to convert units, etc.

One PoPs database may point to another parent database. In that case, if a particle is not found in the current
database it will also check in the parent database before throwing a KeyError.
This option allows an evaluator to override part of the parent database (i.e. to define a different mass
for one or more particles).
"""

import os

from xData import formatVersion as formatVersionModule
from xData import ancestry as ancestryModule
from xData.Documentation import documentation as documentationModule

from .groups import chemicalElements as chemicalElementsModule
from .groups import chemicalElement as chemicalElementModule
from .groups import isotope as isotopeModule
from .families import particle as particleModule
from .families import gaugeBoson as gaugeBosonModule
from .families import lepton as leptonModule
from .families import baryon as baryonModule
from .families import nuclide as nuclideModule
from .families import nucleus as nucleusModule
from .families import unorthodox as unorthodoxModule

from . import styles as stylesModule
from . import alias as aliasModule
from . import parentDatabase as parentDatabaseModule

class database( ancestryModule.ancestry ) :

    moniker = 'PoPs'
    ancestryMembers = ( 'styles', 'documentation', 'aliases', 'gaugeBosons', 'leptons', 'baryons', 'chemicalElements', 'unorthodoxes' )

    def __init__( self, name, version, formatVersion = formatVersionModule.default, parentDatabase = None ) :
        """
        :param name: string naming the database
        :param version: database version
        :param formatVersion: the GNDS format the data are represented in.
        :param parentDatabase: optional link to a parent PoPs database
        """

        ancestryModule.ancestry.__init__( self )

        self.__format = formatVersion

        self.name = name

        self.version = version

        self.__styles = stylesModule.styles( )
        self.__styles.setAncestor( self )

        self.__documentation = documentationModule.Documentation( )
        self.__documentation.setAncestor( self )

        self.__aliases = aliasModule.suite( )
        self.__aliases.setAncestor( self )

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

        self.__parentDatabase = parentDatabase

    def __contains__( self, key ) :

        if( key in self.__aliases ) : return( True )
        if( key in self.__gaugeBosons ) : return( True )
        if( key in self.__leptons ) : return( True )
        if( key in self.__baryons ) : return( True )
        if( key in self.__chemicalElements ) : return( True )
        if( key in self.__unorthodoxes ) : return( True )
        if( self.parentDatabase is not None ) : return( key in self.parentDatabase.link )
        return( False )

    def __getitem__( self, key ) :
        """
        :param key: id or symbol of a particle or particle group
        :return: particle with given key. KeyError is raised if the key is not found.
        """

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
        if( self.parentDatabase is not None ) :
            return( self.parentDatabase.link[ key ] )
        raise KeyError( 'id = "%s" not found' % key )

    def __iter__( self ) :
        """
        Iterates over all particles in the database. Particle groups like the
        chemicalElement and isotope classes are not included, although the particles they contain are.
        """

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
        if( self.parentDatabase is not None ) :
            for item in self.parentDatabase.link : yield item    # FIXME skip particles that are overridden in child database?

    def keys( self ) :
        """Return a list of ids for all particles in self"""

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
    def parentDatabase( self ) :

        return( self.__parentDatabase )

    @parentDatabase.setter
    def parentDatabase( self, value ) :

        if not isinstance( value, parentDatabaseModule.parentDatabase ) :
            raise TypeError("Expected parentDatabase instance, got '%s' instead" % type(value))
        self.__parentDatabase = value

    @property
    def documentation( self ) :
        """Returns the documentation instance."""

        return( self.__documentation )

    def add( self, particle ) :
        """
        Add a particle to the database, inside the appropriate particle family.

        :param particle: instance of families.particle.particle or derived class
        """

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

    def addFile( self, fileName, replace = False ) :
        """
        Opens the file *fileName* and adds its particles to self. If replace is False, particles already in self are
        not replaced.
        """

        pops = self.readFile( fileName )
        for particle in pops :
            if( not( replace ) and particle.id in self ) : continue
            if( isinstance( particle, nucleusModule.particle ) ) : continue
            self.add( particle )

    def check( self, **kwargs ) :
        """
        Check for physics problems in the database

        :param kwargs: options controlling i.e. how strict checks should be. Currently unused.
        :return: list of warnings along with context about where they appear in the database
        """
        from . import warning as warningModule
        warnings = []
        info = kwargs
        info['PoPs'] = self
        for suite in ( self.__gaugeBosons, self.__leptons, self.__baryons, self.__chemicalElements, self.__unorthodoxes, self.__aliases ) :
            suiteWarnings = suite.check( info )
            if suiteWarnings:
                warnings.append( warningModule.context("%s" % suite.moniker, suiteWarnings) )

        result = warningModule.context('PoPs', warnings)
        result.info = info
        return result

    def convertUnits( self, unitMap ) :
        """
        Recursively searches for units within all sections, converting units that appear in the unitMap.

        :param unitMap: dictionary with old/new unit pairs where the old unit is the key (e.g., { 'eV' : 'MeV', 'b' : 'mb' }).
        """

        self.__styles.convertUnits( unitMap )
        self.__gaugeBosons.convertUnits( unitMap )
        self.__leptons.convertUnits( unitMap )
        self.__baryons.convertUnits( unitMap )
        self.__chemicalElements.convertUnits( unitMap )
        self.__unorthodoxes.convertUnits( unitMap )

    def copy( self ) :
        """
        :return: deep copy of the entire database.
        """

        _database = database( self.name, self.version, parentDatabase = self.parentDatabase, formatVersion = self.format )
        for item in self.__styles : _database.styles.add( item.copy( ) )
            # FIXME, need to add documentation.
        for item in self.__aliases : _database.add( item.copy( ) )
        for item in self.__gaugeBosons : _database.add( item.copy( ) )
        for item in self.__leptons : _database.add( item.copy( ) )
        for item in self.__baryons : _database.add( item.copy( ) )
        for item in self.__chemicalElements : _database.add( item.copy( ) )
        for item in self.__unorthodoxes : _database.add( item.copy( ) )

        return( _database )

    def final( self, id ) :
        """
        Returns the particle matching id. If id is an alias, will follow the alias pid until a non-alias particle is found and that particle is returned.
        """

        particle = self[id]
        if( isinstance( particle, particleModule.alias ) ) : particle = particle.particle
        return( particle )

    def hasAlias( self, id ) :

        if id in self.__aliases:
            return True
        if self.parentDatabase is not None:
            return self.parentDatabase.link.hasAlias( id )
        return False

    def saveToFile( self, fileName, **kwargs ) :

        dirname = os.path.dirname( fileName )
        if( dirname != '' ) :
            if( not( os.path.exists( dirname ) ) ) : os.makedirs( dirname )
        fOut = open( fileName, 'w' )
        fOut.write( '<?xml version="1.0"?>\n' )
        fOut.write( self.toXML( **kwargs ) )
        fOut.write( '\n' )
        fOut.close( )

    def toXML( self, indent = "", **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        formatVersion = kwargs.get( 'formatVersion', self.format )
        if( formatVersion not in formatVersionModule.allowed ) : raise Exception( "Unsupported GNDS structure '%s'!" % str( formatVersion ) )
        if( formatVersion == formatVersionModule.version_1_10 ) : formatVersion = '0.1'       # There was no 1.10 for PoPs, it was 0.1. PoPs 0.1 and 2.0 are the same.

        xmlString = [ '%s<%s name="%s" version="%s" format="%s">' % ( indent, self.moniker, self.name, self.version, formatVersion ) ]
        if( self.parentDatabase is not None ) : xmlString += self.parentDatabase.toXMLList( indent2, **kwargs )
        xmlString += self.__styles.toXMLList( indent2, **kwargs )

        if formatVersion == '0.1' and len(self.styles) != 0:
            xmlString += self.styles.getEvaluatedStyle().documentation.toXMLList( indent2, **kwargs )
        # FIXME: currently both PoPs and PoPs/styles/evaluated have documentation nodes according to GNDS-2.0 specs
        xmlString += self.__documentation.toXMLList( indent2, **kwargs )

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

        formatVersion = element.get( 'format' )
        if( formatVersion == '0.1' ) :                             # PoPs 0.1 and 2.0 are the same so allow for now. There was no 1.10 for PoPs.
            formatVersion = formatVersionModule.version_1_10
        else :
            if( formatVersion not in formatVersionModule.allowed ) : raise Exception( "Unsupported GNDS structure '%s'!" % str( formatVersion ) )

        self.__format = formatVersion

        children = {}
        for child in ( gaugeBosonModule, leptonModule, baryonModule, unorthodoxModule, aliasModule ) :
            children[child.suite.moniker] = child.suite

        children[chemicalElementsModule.chemicalElements.moniker] = chemicalElementsModule.chemicalElements

        for child in element :
            if( child.tag in children ) :
                children[child.tag].parseXMLNode( getattr( self, child.tag ), child, xPath, linkData )
            elif( child.tag == self.__styles.moniker ) :
                self.__styles.parseXMLNode( child, xPath, linkData )
            elif( child.tag == 'parentDatabase' ) :
                self.parentDatabase = parentDatabaseModule.parentDatabase.parseXMLNode( child, xPath, linkData )
            elif( child.tag == self.documentation.moniker ) :
                if formatVersion == formatVersionModule.version_1_10 :
                    self.styles.getEvaluatedStyle().documentation.parseNode( child, xPath, linkData )
                else:
                    self.documentation.parseNode( child, xPath, linkData )
            else :
                raise ValueError( 'sub-element = "%s" not allowed' % child.tag )

        xPath.pop( )
        return( self )

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ) :

        formatVersion = element.get( 'format' )
        if( formatVersion == '0.1' ) :                             # PoPs 0.1 and 2.0 are the same so allow for now. There was no 1.10 for PoPs.
            formatVersion = formatVersionModule.version_1_10
        else :
            if( formatVersion not in formatVersionModule.allowed ) : raise Exception( "Unsupported GNDS structure '%s'!" % str( formatVersion ) )

        self = cls( element.get( 'name' ), element.get( 'version' ), formatVersion = formatVersion )
        self.parseXMLNode( element, xPath, linkData )
        return( self )

    @classmethod
    def parseXMLStringAsClass( cls, string ) :  # FIXME move this to PoPs/Test? That's the only place it is used

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
        """
        Load a PoPs database from file

        :param fileName: path to the file
        :return: PoPs.database instance containing all particles in the file
        """

        from xml.etree import cElementTree

        xPath = []
        try :
            element = cElementTree.parse( fileName ).getroot( )
            return( database.parseXMLNodeAsClass( element, xPath, {} ) )
        except :
            print( "Error encountered at xpath = /%s" % '/'.join( xPath ) )
            raise
