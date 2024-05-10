# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the classes used to support a coherent elastic reaction for the thermal neutron scattering law.

This module contains the following classes: 
        
    +--------------+--------------------------------------------------------------------------------------------+
    | Class        | Description                                                                                |
    +==============+============================================================================================+
    | S_table      | This class represents the coherent elastic scattering :math:`S(E|T)` function.             |
    +--------------+--------------------------------------------------------------------------------------------+
    | Form         | This class represents pure-Coulomb scattering.                                             |
    +--------------+--------------------------------------------------------------------------------------------+
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
    This class represents the coherent elastic scattering :math:`S(E|T)` function.
    """

    moniker = 'S_table'

    def __init__(self, gridded2d):

        super().__init__()
        self.gridded2d = gridded2d

    @property
    def domainUnit( self ) :
        """
        This method rReturns the energy unit of the projectile.
        """

        return( self.gridded2d.axes[1].unit )

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.
        
        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        self.gridded2d.convertUnits( unitMap )

    def toXML_strList( self, indent = '', **kwargs ) :
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ indent + '<%s>' % self.moniker ]
        xmlStringList += self.gridded2d.toXML_strList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        """
        Parse *element* into an instance of *cls*.
        
        :param cls:         Form class to return.
        :param element:     Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.
        
        :return: an instance of *cls* representing *element*.
        """

        xPath.append( element.tag )

        gridded2d = griddedModule.Gridded2d.parseNodeUsingClass(element[0], xPath, linkData, **kwargs)
        _S_table = cls( gridded2d )

        xPath.pop()

        return _S_table

class Form( baseModule.Form ) :
    """
    This class is the form for the coherent elastic reaction for the thermal neutron scattering law.

    The following table list the primary members of this class:
    
    +-----------+---------------------------------------------------------------+
    | Member    | Description                                                   |
    +===========+===============================================================+
    | label     | This member is the form's label.                              |
    +-----------+---------------------------------------------------------------+
    | S_table   | This member is a :py:class:`S_table` that contains the data.  |
    +-----------+---------------------------------------------------------------+
    """

    moniker = 'thermalNeutronScatteringLaw_coherentElastic'
    keyName = 'label'

    process = 'thermalNeutronScatteringLaw coherent-elastic'
    subformAttributes = ( 'S_table', )

    def __init__( self, label, _S_table ) :
        """
        :param label:           The label for the form.
        :param _S_table:        An :py:class:`S_table` instance.
        """

        if( not( isinstance( _S_table, S_table ) ) ) : raise TypeError( "Invalid S_table for %s." % self.moniker )

        baseModule.Form.__init__( self, IDsPoPsModule.neutron, label, xDataEnumsModule.Frame.lab, ( _S_table, ) )

    @property
    def domainUnit( self ) :
        """
        This method rReturns the energy unit of the projectile.
        """

        return( self.S_table.domainUnit )

    def processThermalNeutronScatteringLaw( self, style, kwargs ) :
        """
        This method calculates the cross section, distribution, average product energy and average product momentum
        data for the double differential data and returns their forms. See :py:func:`coherentElasticMiscModule.process` for more details.

        :param style:           The style for the returned data.
        :param kwargs:          A dictionary containing data needed for processing.
        """

        temperature = style.temperature.getValueAs( kwargs['temperatureUnit'] )

        energyMin = PQUModule.PQU( 1e-11, 'MeV' ).getValueAs( kwargs['incidentEnergyUnit'] )
        energyMin = kwargs.get( 'energyMin', energyMin )
        energyMax = PQUModule.PQU(  5e-6, 'MeV' ).getValueAs( kwargs['incidentEnergyUnit'] )
        energyMax = kwargs.get( 'energyMax', energyMax )

        return( coherentElasticMiscModule.process( self, style.label, energyMin, energyMax, temperature, kwargs ) )

    def temperatures( self, _temperatures ) :
        """
        This method add to *_temperatures* the temperatures that *self* has been evaluated at.

        :param _temperatures:   A dictionary with the key 'coherent-elastic' added by this method.
        """

        _temperatures['coherent-elastic'] = [ self.S_table.gridded2d.axes[2].unit, self.S_table.gridded2d.axes[2].values.values ]

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        """
        Parse *element* into an instance of *cls*.
        
        :param cls:         Form class to return.
        :param element:     Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.
        
        :return: an instance of *cls* representing *element*.
        """

        xPath.append( element.tag )

        _S_table = S_table.parseNodeUsingClass(element[0], xPath, linkData, **kwargs)
        coherentElastic = cls( element.get( 'label' ), _S_table )

        xPath.pop( )

        return coherentElastic
