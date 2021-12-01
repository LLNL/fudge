# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Thermal neutron scattering law coherent elastic double difference cross section form and its supporting classes.
"""

__metaclass__ = type

from pqu import PQU as PQUModule

from xData import ancestry as ancestryModule
from xData import standards as standardsModule
from xData import gridded as griddedModule

from PoPs import IDs as IDsPoPsModule

from . import base as baseModule
from . import coherentElasticMisc as coherentElasticMiscModule

class S_table( ancestryModule.ancestry ) :
    """
    For elastic coherent, cumulative structure factor 'S' is given as function of incident energy and temperature.
    """

    moniker = 'S_table'

    def __init__(self, gridded2d):

        super().__init__()
        self.gridded2d = gridded2d

    @property
    def domainUnit( self ) :

        return( self.gridded2d.axes[1].unit )

    def convertUnits( self, unitMap ) :
        """
        unitMap is a dictionary of the for { 'eV' : 'MeV', 'b' : 'mb' }.
        """

        self.gridded2d.convertUnits( unitMap )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ indent + '<%s>' % self.moniker ]
        xmlStringList += self.gridded2d.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( element.tag )
        gridded2d = griddedModule.gridded2d.parseXMLNode( element[0], xPath, linkData )
        _S_table = S_table( gridded2d )

        xPath.pop()

        return( _S_table )

class form( baseModule.form ) :

    moniker = 'thermalNeutronScatteringLaw_coherentElastic'
    process = 'thermalNeutronScatteringLaw coherent-elastic'
    subformAttributes = ( 'S_table', )

    def __init__( self, label, _S_table ) :

        if( not( isinstance( _S_table, S_table ) ) ) : raise TypeError( "Invalid S_table for %s." % self.moniker )

        baseModule.form.__init__( self, IDsPoPsModule.neutron, label, standardsModule.frames.labToken, ( _S_table, ) )

    @property
    def domainUnit( self ) :

        return( self.S_table.domainUnit )

    def processThermalNeutronScatteringLaw( self, style, kwargs ) :

        temperature = style.temperature.getValueAs( kwargs['temperatureUnit'] )

        energyMin = PQUModule.PQU( 1e-11, 'MeV' ).getValueAs( kwargs['incidentEnergyUnit'] )
        energyMin = kwargs.get( 'energyMin', energyMin )
        energyMax = PQUModule.PQU(  5e-6, 'MeV' ).getValueAs( kwargs['incidentEnergyUnit'] )
        energyMax = kwargs.get( 'energyMax', energyMax )

        return( coherentElasticMiscModule.process( self, style.label, energyMin, energyMax, temperature, kwargs ) )

    def temperatures( self, _temperatures ) :

        _temperatures['coherent-elastic'] = [ self.S_table.gridded2d.axes[2].unit, self.S_table.gridded2d.axes[2].values.values ]

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )

        _S_table = S_table.parseXMLNode( element[0], xPath, linkData )
        coherentElastic = form( element.get( 'label' ), _S_table )

        xPath.pop( )
        return( coherentElastic )
