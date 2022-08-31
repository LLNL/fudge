# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Thermal neutron scattering law coherent elastic double difference cross section form and its supporting classes.
"""


from pqu import PQU as PQUModule

from LUPY import ancestry as ancestryModule

from xData import enums as xDataEnumsModule
from xData import gridded as griddedModule

from PoPs import IDs as IDsPoPsModule

from . import base as baseModule
from . import coherentElasticMisc as coherentElasticMiscModule

class S_table( ancestryModule.AncestryIO ) :
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

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ indent + '<%s>' % self.moniker ]
        xmlStringList += self.gridded2d.toXML_strList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )

        gridded2d = griddedModule.Gridded2d.parseNodeUsingClass(element[0], xPath, linkData, **kwargs)
        _S_table = cls( gridded2d )

        xPath.pop()

        return _S_table

class Form( baseModule.Form ) :

    moniker = 'thermalNeutronScatteringLaw_coherentElastic'
    keyName = 'label'

    process = 'thermalNeutronScatteringLaw coherent-elastic'
    subformAttributes = ( 'S_table', )

    def __init__( self, label, _S_table ) :

        if( not( isinstance( _S_table, S_table ) ) ) : raise TypeError( "Invalid S_table for %s." % self.moniker )

        baseModule.Form.__init__( self, IDsPoPsModule.neutron, label, xDataEnumsModule.Frame.lab, ( _S_table, ) )

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

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )

        _S_table = S_table.parseNodeUsingClass(element[0], xPath, linkData, **kwargs)
        coherentElastic = cls( element.get( 'label' ), _S_table )

        xPath.pop( )

        return coherentElastic
