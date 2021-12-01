# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module defines the abstract base class for all types of particle appearing in PoPs.
"""

import abc

from .. import misc as miscModule
from .. import suite as suiteModule
from .. import warning as warningModule

from ..quantities import quantity as quantityModule
from ..quantities import mass as massModule
from ..quantities import spin as spinModule
from ..quantities import parity as parityModule
from ..quantities import charge as chargeModule
from ..quantities import halflife as halflifeModule

from ..decays import decayData as decayDataModule

class alias( miscModule.classWithIDKey ) :

    __metaclass__ = abc.ABCMeta

    def __init__( self, id, particle ) :

        miscModule.classWithIDKey.__init__( self, id )

        self.__particle = particle

    @property
    def particle( self ) :

        return( self.__particle )

    @property
    def pid( self ) :

        return( self.__particle.id )

    @property
    def isAnti( self ) :

        return( self.__particle.isAnti )

    @property
    def family( self ) :

        return( self.__particle.family )

    @property
    def mass( self ) :

        return( self.__particle.mass )

    @property
    def spin( self ) :

        return( self.__particle.spin )

    @property
    def parity( self ) :

        return( self.__particle.parity )

    @property
    def charge( self ) :

        return( self.__particle.charge )

    @property
    def halflife( self ) :

        return( self.__particle.halflife )

    def getMass( self, unit ) :

        return( self.__particle.getMass( unit ) )

    def isAlias( self ) :

        return( True )

class particle( miscModule.classWithIDKey ) :
    """
    This is the abstract base class for all particles.
    """

    __metaclass__ = abc.ABCMeta

    def __init__( self, id ) :
        """
        :param id: The particle id. Every particle in a PoPs database must have a unique id.
        """

        miscModule.classWithIDKey.__init__( self, id )

        base, anti, qualifier = miscModule.baseAntiQualifierFromID( id )

        self.__anti = anti == miscModule.antiSuffix

        self.__mass = self.addSuite( massModule )
        self.__spin = self.addSuite( spinModule )
        self.__parity = self.addSuite( parityModule )
        self.__charge = self.addSuite( chargeModule )
        self.__halflife = self.addSuite( halflifeModule )
        self.__decayData = decayDataModule.decayData( )
        self.__decayData.setAncestor( self )

    def __eq__( self, other ) :

        return( self.id == other.id )

    @property
    @abc.abstractmethod
    def familyOrder( self ) :

        pass

    @property
    def isAnti( self ) :
        """
        :return: whether this is an anti-particle
        """

        return( self.__anti )

    @property
    def family( self ) :
        """
        :return: family name of the particle
        """

        return( self.moniker )

    @property
    def mass( self ) :
        """
        :return: particle mass
        """

        return( self.__mass )

    @property
    def spin( self ) :
        """
        :return: particle spin
        """

        return( self.__spin )

    @property
    def parity( self ) :
        """
        :return: particle parity
        """

        return( self.__parity )

    @property
    def charge( self ) :
        """
        :return: particle charge
        """

        return( self.__charge )

    @property
    def halflife( self ) :
        """
        :return: particle halflife
        """

        return( self.__halflife )

    @property
    def decayData( self ) :
        """
        :return: particle decay information
        """

        return( self.__decayData )

    def addSuite( self, module ) :      # FIXME rename this __addSuite? Not intended for public use

        suite = module.suite( )
        suite.setAncestor( self )
        return( suite )

    def check( self, info ):
        """
        Basic physics checking for particles.

        Options & defaults:
            'branchingRatioSumTolerance'      1e-6
        """
        warnings = []

        BRSumAbsTol = info.get('branchingRatioSumTolerance', 1e-6)
# CMM FIXME
#        probabilitySum = sum( [decay.probability for decay in self.decays] )
#        if self.decays and abs(probabilitySum - 1.0) > BRSumAbsTol:
#            warnings.append(warningModule.unnormalizedDecayProbabilities(probabilitySum, self))
        return warnings

    def buildFromRawData( self, mass = None, spin = None, parity = None, charge = None, halflife = None, label = 'default' ) :
        """
        Helper method for quickly adding properties to a particle

        :param mass:  tuple(float value, string unit)
        :param spin:  tuple(int or fraction value, string unit)
        :param parity: tuple(int value, string unit)
        :param charge: tuple(int value, string unit)
        :param halflife: tuple(float or string value, string unit)
        :param label: style label (string) to apply to each quantity
        """

        def getUnit( unit ) :

            if( isinstance( unit, str ) ) : unit = quantityModule.stringToPhysicalUnit( unit )
            return( unit )

        if( mass is not None ) : self.mass.add( massModule.double( label, mass[0], getUnit( mass[1] ) ) )
        if( spin is not None ) :
            self.spin.add( spinModule.fraction( label, spin[0], spin[1] ) )
        if( parity is not None ) : self.parity.add( parityModule.integer( label, parity[0], getUnit( parity[1] ) ) )
        if( charge is not None ) : self.charge.add( chargeModule.integer( label, charge[0], getUnit( charge[1] ) ) )
        if( halflife is not None ) :
            if( isinstance( halflife[0], float ) ) :
                _halflife = halflifeModule.double( label, halflife[0], getUnit( halflife[1] ) )
            elif( isinstance( halflife[0], str ) ) :
                _halflife = halflifeModule.string( label, halflife[0], getUnit( halflife[1] ) )
            else :
                TypeError( 'halflife neither double or str instance' )
            self.halflife.add( _halflife )

    def convertUnits( self, unitMap ) :
        """ See convertUnits documentation in PoPs.database """

        self.__mass.convertUnits( unitMap )
        self.__spin.convertUnits( unitMap )
        self.__parity.convertUnits( unitMap )
        self.__charge.convertUnits( unitMap )
        self.__halflife.convertUnits( unitMap )
        self.__decayData.convertUnits( unitMap )

    def copy( self ) :
        """
        :return: deep copy of self
        """

        _particle = self.__class__( self.id )
        self.__copyStandardQuantities( _particle )
        return( _particle )

    def __copyStandardQuantities( self, other ) :

        if( other.moniker != self.moniker ) : TypeError( '''other's family = "%s" not same as self's family = "%s".''' %
            ( other.moniker, self.moniker ) )
        for item in self.__mass : other.mass.add( item.copy( ) )
        for item in self.__spin : other.spin.add( item.copy( ) )
        for item in self.__parity : other.parity.add( item.copy( ) )
        for item in self.__charge : other.charge.add( item.copy( ) )
        for item in self.__halflife : other.halflife.add( item.copy( ) )
        self.__decayData.copyItems( other.decayData )

    def familyOrderLessThan( self, other ) :

        if( not( isinstance( other, particle ) ) ) : raise TypeError( 'Other must be a particle: %s.' % other.__class__ )
        if( self.familyOrder < other.familyOrder ) : return( True )
        return( False )

    def replicate( self, other ) :
        """
        Copy data from other into self

        :param other: another particle instance
        """

        self.__mass.replicate( other.mass )
        self.__spin.replicate( other.spin )
        self.__parity.replicate( other.parity )
        self.__charge.replicate( other.charge )
        self.__halflife.replicate( other.halflife )
        self.__decayData = other.decayData.copy( )

    def getMass( self, unit ) :
        """
        Evaluate the particle mass in the desired unit.

        :param unit: desired unit (string)
        :return: mass (float)
        """

        return( self.mass[0].float( unit ) )

    def extraXMLAttributes( self ) :

        return( '' )

    def extraXMLElements( self, indent = '', **kwargs ) :

        return( [] )

    def isAlias( self ) :

        return( False )

    def sortCompare( self, other ) :
        """This should be defined in each sub-class."""

        if( not( isinstance( other, self.__class__ ) ) ) : raise TypeError( 'Invalid other.' )
        if( self.id > other.id ) : return( 1 )
        if( self.id == other.id ) : return( 0 )
        return( -1 )

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs )  ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s id="%s"%s>' % ( indent, self.moniker, self.id, self.extraXMLAttributes( ) ) ]

        XMLStringList += self.mass.toXMLList( indent = indent2, **kwargs )
        XMLStringList += self.spin.toXMLList( indent = indent2, **kwargs )
        XMLStringList += self.parity.toXMLList( indent = indent2, **kwargs )
        XMLStringList += self.charge.toXMLList( indent = indent2, **kwargs )
        XMLStringList += self.halflife.toXMLList( indent = indent2, **kwargs )

        XMLStringList += self.decayData.toXMLList( indent = indent2, **kwargs )

        XMLStringList += self.extraXMLElements( indent2, **kwargs )

        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    def parseExtraXMLElement( self, element, xPath, linkData ) :

        return( False )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( '%s[@id="%s"]' % ( element.tag, element.get( 'id' ) ) )

        children = { 'mass' : massModule,       'spin' : spinModule,            'parity' : parityModule, 
                     'charge' : chargeModule,   'halflife' : halflifeModule }
        for child in element :
            if( child.tag in children ) :
                children[child.tag].suite.parseXMLNode( getattr( self, child.tag ), child, xPath, linkData )
            elif( child.tag == decayDataModule.decayData.moniker ) :  # not a suite so can't be treated like other children
                self.decayData.parseXMLNode( child, xPath, linkData )
            else :
                if( not( self.parseExtraXMLElement( child, xPath, linkData ) ) ) :
                    raise ValueError( 'sub-element = "%s" not allowed' % child.tag )

        xPath.pop( )
        return( self )

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ) :

        xPath.append( '%s[@id="%s"]' % ( element.tag, element.get( 'id' ) ) )

        kwargs = {v[0]:v[1] for v in element.items()}
        del kwargs['id']
        self = cls( element.get('id'), **kwargs )
        xPath.pop()

        self.parseXMLNode( element, xPath, linkData )

        return( self )

    @classmethod
    def parseXMLStringAsClass( cls, string ) :

        from xml.etree import cElementTree

        element = cElementTree.fromstring( string )
        return( cls.parseXMLNodeAsClass( element, [], [] ) )

class suite( suiteModule.sortedSuite ) :

    def __init__( self, replace = True ) :

        suiteModule.sortedSuite.__init__( self, allowedClasses = ( self.particle, ), replace = replace )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( element.tag )

        for child in element :
            self.add( self.particle.parseXMLNodeAsClass( child,  xPath, linkData ) )

        xPath.pop( )

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ) :

        xPath.append( element.tag )

        self = cls( )
        for child in element :
            self.add( self.particle.parseXMLNodeAsClass( child,  xPath, linkData ) )

        xPath.pop( )
        return( self )

    @classmethod
    def parseXMLStringAsClass( cls, string ) :

        from xml.etree import cElementTree

        element = cElementTree.fromstring( string )
        return( cls.parseXMLNodeAsClass( element, [], [] ) )
