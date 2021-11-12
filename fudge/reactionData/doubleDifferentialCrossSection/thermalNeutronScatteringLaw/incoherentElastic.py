# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from pqu import PQU as PQUModule

from xData import ancestry as ancestryModule
from xData import standards as standardsModule
from xData import physicalQuantity as physicalQuantityModule
from xData import XYs as XYs1dModule

from PoPs import IDs as IDsPoPsModule

from . import base as baseModule
from . import incoherentElasticMisc as incoherentElasticMiscModule

"""
Thermal neutron scattering law incoherent elastic double difference cross section form and its supporting classes.
"""

__metaclass__ = type

class characteristicCrossSection( physicalQuantityModule.physicalQuantity ):

    moniker = 'characteristicCrossSection'

class DebyeWaller( ancestryModule.ancestry ) :
    """
    For incoherent elastic sections, all we need is a characteristic cross section and a
    temperature-dependent list of the Debye-Waller integral.
    """

    moniker = "DebyeWaller"

    def __init__( self, function1d ) :

        ancestryModule.ancestry.__init__( self )

        self.__function1d = function1d
        self.__function1d.setAncestor( self )

    @property
    def function1d( self ) :

        return( self.__function1d )

    def convertUnits( self, unitMap ) :
        """
        unitMap is a dictionary of the for { 'eV' : 'MeV', 'b' : 'mb' }.
        """

        self.__function1d.convertUnits( unitMap )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        xmlStringList += self.__function1d.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker

        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )

        xys1d = XYs1dModule.XYs1d.parseXMLNode( element[0], xPath, linkData )
        DebyeWaller1 = DebyeWaller( xys1d )

        xPath.pop( )
        return( DebyeWaller1 )

class form( baseModule.form ) :

    moniker = 'thermalNeutronScatteringLaw_incoherentElastic'
    process = 'thermalNeutronScatteringLaw incoherent-elastic'
    subformAttributes = ( 'characteristicCrossSection', 'DebyeWaller' )

    def __init__( self, label, _characteristicCrossSection, _DebyeWaller ) :

        baseModule.form.__init__( self, IDsPoPsModule.neutron, label, standardsModule.frames.labToken, ( _characteristicCrossSection, _DebyeWaller, ) )

    @property
    def domainUnit( self ) :

        invUnit = 1 / PQUModule.PQU( 1, self.DebyeWaller.function1d.axes[0].unit )
        return( str( invUnit.unit ) )

    def processThermalNeutronScatteringLaw( self, style, kwargs ) :

        temperature = style.temperature.getValueAs( self.DebyeWaller.function1d.axes[1].unit )

        energyMin = PQUModule.PQU( 1e-11, 'MeV' ).getValueAs( kwargs['incidentEnergyUnit'] )
        energyMin = kwargs.get( 'energyMin', energyMin )
        energyMax = PQUModule.PQU(  5e-6, 'MeV' ).getValueAs( kwargs['incidentEnergyUnit'] )
        energyMax = kwargs.get( 'energyMax', energyMax )

        return( incoherentElasticMiscModule.process( self, style.label, energyMin, energyMax, temperature, kwargs ) )

    def temperatures( self, _temperatures ) :

        _temperatures['incoherent-elastic'] = [ self.DebyeWaller.function1d.axes[-1].unit, [ x for x, y in self.DebyeWaller.function1d ] ]

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )

        _characteristicCrossSection = characteristicCrossSection.parseXMLNode( element[0], xPath, linkData )
        _DebyeWaller = DebyeWaller.parseXMLNode( element[1], xPath, linkData )
        incoherentElastic = form( element.get( 'label' ), _characteristicCrossSection, _DebyeWaller )

        xPath.pop( )
        return( incoherentElastic )
