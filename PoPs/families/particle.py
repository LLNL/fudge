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
This module contains the parity classes.
"""

import abc

from xData import ancestry as ancestryModule

from .. import misc as miscModule
from .. import suite as suiteModule
from .. import warning as warningModule

from ..quantities import quantity as quantityModule
from ..quantities import mass as massModule
from ..quantities import spin as spinModule
from ..quantities import parity as parityModule
from ..quantities import charge as chargeModule
from ..quantities import halflife as halflifeModule

from ..decays import decay as decayModule
from ..decays import decayData as decayDataModule

class alias( ancestryModule.ancestry ) :

    __metaclass__ = abc.ABCMeta

    def __init__( self, id, particle ) :

        ancestryModule.ancestry.__init__( self )

        if( not( isinstance( id, str ) ) ) : TypeError( '''id must be a str instance.''' )
        self.__id = id

        self.__particle = particle

    @property
    def id( self ) :

        return( self.__id )

    @property
    def key( self ) :

        return( self.__id )

    @key.setter
    def key( self, value ) :

        if( not( isinstance( value, str ) ) ) : raise TypeError( 'id must be a string instance.' )
        self.__id = value

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

class particle( ancestryModule.ancestry ) :
    """
    This is the abstract base class for all particles.
    """

    __metaclass__ = abc.ABCMeta

    def __init__( self, id ) :

        base, anti, qualifier = miscModule.baseAntiQualifierFromID( id )

        ancestryModule.ancestry.__init__( self )

        self.__id = id

        self.__anti = anti == miscModule.antiSuffix

        self.__mass = self.addSuite( massModule )
        self.__spin = self.addSuite( spinModule )
        self.__parity = self.addSuite( parityModule )
        self.__charge = self.addSuite( chargeModule )
        self.__halflife = self.addSuite( halflifeModule )
        self.__decays = self.addSuite( decayModule )
        self.__decayData = decayDataModule.decayData( )
        self.__decayData.setAncestor( self )

    @property
    def id( self ) :

        return( self.__id )

    @property
    def key( self ) :

        return( self.__id )

    @key.setter
    def key( self, value ) :

        if( not( isinstance( value, str ) ) ) : raise TypeError( 'id must be a string instance.' )
        self.__id = value

    @property
    def isAnti( self ) :

        return( self.__anti )

    @property
    def family( self ) :

        return( self.moniker )

    @property
    def mass( self ) :

        return( self.__mass )

    @property
    def spin( self ) :

        return( self.__spin )

    @property
    def parity( self ) :

        return( self.__parity )

    @property
    def charge( self ) :

        return( self.__charge )

    @property
    def halflife( self ) :

        return( self.__halflife )

    @property
    def decays( self ) :

        return( self.__decays )

    @property
    def decayData( self ) :

        return( self.__decayData )

    def addSuite( self, module ) :

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
        probabilitySum = sum( [decay.probability for decay in self.decays] )
        if self.decays and abs(probabilitySum - 1.0) > BRSumAbsTol:
            warnings.append(warningModule.unnormalizedDecayProbabilities(probabilitySum, self))
        return warnings

    def buildFromRawData( self, mass = None, spin = None, parity = None, charge = None, halflife = None, label = 'default' ) :

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

        self.__mass.convertUnits( unitMap )
        self.__spin.convertUnits( unitMap )
        self.__parity.convertUnits( unitMap )
        self.__charge.convertUnits( unitMap )
        self.__halflife.convertUnits( unitMap )
        self.__decays.convertUnits( unitMap )
        self.__decayData.convertUnits( unitMap )

    def copy( self ) :

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
        for item in self.__decays : other.decays.add( item.copy( ) )
        self.__decayData.copyItems( other.decayData )

    def getMass( self, unit ) :

        return( self.mass[0].float( unit ) )

    def extraXMLAttributes( self ) :

        return( '' )

    def extraXMLElements( self, indent = '', **kwargs ) :

        return( [] )

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

        XMLStringList += self.decays.toXMLList( indent = indent2, **kwargs )
        XMLStringList += self.decayData.toXMLList( indent = indent2, **kwargs )

        XMLStringList += self.extraXMLElements( indent2, **kwargs )

        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    def parseExtraXMLElement( self, element, xPath, linkData ) :

        return( False )

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ) :

        xPath.append( element.tag )

        kwargs = element.attrib.copy( )
        del kwargs['id']
        self = cls( element.attrib['id'], **kwargs )

        children = { 'mass' : massModule,       'spin' : spinModule,            'parity' : parityModule, 
                     'charge' : chargeModule,   'halflife' : halflifeModule,    'decays' : decayModule }
        for child in element :
            if( child.tag in children ) :
                children[child.tag].suite.parseXMLNode( getattr( self, child.tag ), child, xPath, linkData )
            elif( child.tag == 'decayData' ) :  # not a suite so can't be treated like other children
                self.decayData.parseXMLNode( child, xPath, linkData )
            else :
                if( not( self.parseExtraXMLElement( child, xPath, linkData ) ) ) :
                    raise ValueError( 'sub-element = "%s" not allowed' % child.tag )

        xPath.pop( )
        return( self )

    @classmethod
    def parseXMLStringAsClass( cls, string ) :

        from xml.etree import cElementTree

        element = cElementTree.fromstring( string )
        return( cls.parseXMLNodeAsClass( element, [], [] ) )

class suite( suiteModule.sortedSuite ) :

    def __init__( self, replace = True ) :

        suiteModule.sortedSuite.__init__( self, allowedClasses = ( self.particle, ), key = 'id', replace = replace )

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
