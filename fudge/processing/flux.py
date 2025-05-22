# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

r"""
This module contains classes for storing the data in a fluxes node.

This module contains the following classes:

    +---------------------------------------+-----------------------------------------------------------------------------------+
    | Class                                 | Description                                                                       |
    +=======================================+===================================================================================+
    | LegendreSeries                        | This class represents the angular part of the flux as a Legendre series.          |
    +---------------------------------------+-----------------------------------------------------------------------------------+
    | XYs2d                                 | This class represents the product energy and angular part of the flux.            |
    +---------------------------------------+-----------------------------------------------------------------------------------+
    | Gridded2d                             | This class represents a multi-group representation of the flux :math:`f(E,\mu)`.  |
    +---------------------------------------+-----------------------------------------------------------------------------------+
    | XYs3d                                 | This class represents the flux :math:`f(T,E,\mu)`.                                |
    +---------------------------------------+-----------------------------------------------------------------------------------+
    | Flux                                  | This class stores the flux as an XYs2d instance.                                  |
    +---------------------------------------+-----------------------------------------------------------------------------------+
    | Fluxes                                | This class stores a list of :py:class:Flux: instances.                            |
    +---------------------------------------+-----------------------------------------------------------------------------------+
"""

from LUPY import ancestry as ancestryModule

from xData import enums as xDataEnumsModule
from xData import axes as axesModule
from xData import XYs1d as XYs1dModule
from xData import series1d as series1dModule
from xData import gridded as griddedModule
from xData import xDataArray as arrayModule
from xData import values as valuesModule
from xData import multiD_XYs as multiD_XYsModule

from fudge import suites as suitesModule

def axes( energyUnit, fluxUnit = '1/s', temperatureUnit = 'MeV/k' ) :
    """
    The function returns an :py:class:`axesModule.Axes` instance that represents the data in a :py:class:`Flux` instance.

    :param energyUnit:          The unit of the incident energy.
    :param fluxUnit:            The unit of the flux.
    :param temperatureUnit:     The unit of the temperature.

    :returns:                   A :py:class:`axesModule.Axes` instance.
    """

    axes = axesModule.Axes(4)
    axes[0] = axesModule.Axis( 'flux', 0, fluxUnit )
    axes[1] = axesModule.Axis( 'mu', 0, '' )
    axes[2] = axesModule.Axis( 'energy_in', 0, energyUnit )
    axes[3] = axesModule.Axis( 'temperature', 0, temperatureUnit )
    return( axes )

class LegendreSeries( series1dModule.LegendreSeries ) :
    """
    This class represents the angular part of the flux as a Legendre series.
    """

    pass

class XYs2d( multiD_XYsModule.XYs2d ) :
    """
    This class represents the product energy and angular part of the flux.
    """

    def getFluxAtLegendreOrder( self, order ) :
        """
        This function returns an :py:class:`XYs1dModule.XYs1d` representing *self* at Legendre order *order*.

        :param order;           The desired order of *self*.

        :returns:               A :py:class:`XYs1dModule.XYs1d` instance.
        """

        axes = axesModule.Axes(2)
        axes[1] = self.axes[2].copy()
        axes[0] = self.axes[0].copy()
        for LS in self :
            if( len( LS ) > 1 ) : raise Exception( 'FIXME -- next line is wrong.' )
        return( XYs1dModule.XYs1d( [ [ LS.outerDomainValue, LS[0] ] for LS in self ], axes = axes, interpolation=self.interpolation ) )

    def processMultiGroup( self, style, tempInfo, indent ) :
        """
        The function returns a multi-group representation of *self*.

        :param style:           This is the multi-group style for the multi-group data.                                             
        :param tempInfo:        This is a dictionary with needed data.
        :param indent:          The amount of indentation if messages are printed.

        :returns:               ?
        """

        # BRB, Currently only processes LegendreSeries for 1d data and only its l=0 term. Need to convert units if needed.
        # The return instance is a list, this needs to be changed.

        from fudge.processing import miscellaneous as miscellaneousModule

        flux0 = self.getFluxAtLegendreOrder( 0 )
        projectileName = tempInfo['reactionSuite'].projectile
        return( flux0.groupOneFunction( style.transportables[projectileName].group.boundaries ) )

    @staticmethod
    def allowedSubElements( ) :
        """
        This function returns a python list of the allowed 1d-function classes that can be added to a :py:class:`XYs2d` instance.
        """

        return( ( LegendreSeries, ) )

class Gridded2d( griddedModule.Gridded2d ) :
    r"""
    This class represents a multi-group repreentation of the flux :math:`f(E,\mu)`.
    """

    pass

class XYs3d( multiD_XYsModule.XYs3d ) :
    r"""
    This class represents the flux :math:`f(T,E,\mu)` where :math:`T` is temperature, :math:`E` is the projectile energy
    and :math:`\mu` is the cosine an angle.
    """

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs2d, ) )

class Flux( ancestryModule.AncestryIO ) :
    """
    This class stores the flux as an XYs2d instance.
    """

    moniker = 'flux'

    def __init__( self, label, data ) :
        """
        :param label:       A unique label for a form.
        :param data:        The flux data.
        """

        ancestryModule.AncestryIO.__init__(self)

        if( not( isinstance( label, str ) ) ) : raise TypeError( 'label must be a string instance.' )
        self.label = label

        if( not( isinstance( data, ( XYs2d, ) ) ) ) : raise TypeError( 'Invalid flux data.' )
        self.data = data
        self.data.setAncestor( self )

    def convertUnits(self, unitMap):
        '''
        Converts unit per *unitMap*.

        :param unitMap:                 A dictionary that maps a existing unit into a new unit.
        '''

        self.data.convertUnits(unitMap)

    def getFluxAtLegendreOrder( self, order ) :
        """
        This function returns an :py:class:`XYs1dModule.XYs1d` representing *self* at Legendre order *order*.

        :param order;           The desired order of *self*.

        :returns:               A :py:class:`XYs1dModule.XYs1d` instance.
        """

        return( self.data.getFluxAtLegendreOrder( order ) )

    def processMultiGroup( self, style, tempInfo, indent ) :
        """
        This function returns a multi-group representation of *self* as a :py:class:`Gridded2d` instance.

        :param style:           This is the multi-group style for the multi-group data.                                             
        :param tempInfo:        This is a dictionary with needed data.
        :param indent:          The amount of indentation if messages are printed.

        :returns:               A :py:class:`Gridded2d` instance.
        """

        reactionSuite = tempInfo['reactionSuite']
        axes = self.data.axes.copy( )
        axes[1] = axesModule.Grid(axes[1].label, 1, axes[1].unit, style=xDataEnumsModule.GridStyle.parameters,
                values=valuesModule.Values([0], valueType=xDataEnumsModule.ValueType.integer32))
        axes[2] = style.transportables[reactionSuite.projectile].group.boundaries.copy()
        data = self.data.processMultiGroup( style, tempInfo, indent )
        starts = valuesModule.Values([0], valueType=xDataEnumsModule.ValueType.integer32)
        lengths = valuesModule.Values([len(data)], valueType=xDataEnumsModule.ValueType.integer32)
        array = arrayModule.Flattened( shape = ( len( data ), 1 ), data = data, starts = starts, lengths = lengths )
        data = Gridded2d( axes = axes, array = array )
        return( data )

    def toXML_strList( self, indent = '', **kwargs ) :
        """
        Returns a python list of str instances representing the XML lines of *self*.
        
        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.
        
        :return:                Python list of str instances.
        """ 

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s label="%s">' % ( indent, self.moniker, self.label ) ]
        xmlStringList += self.data.toXML_strList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance of *cls*.

        :param cls:         Form class to return.
        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    Dictionary that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return: an instance of *cls* representing *node*.
        """

        xPath.append( node.tag )

        label = node.get('label')
        for child in node:
            if( child.tag == XYs2d.moniker ) :
                fluxData = XYs2d.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            else :
                raise 'Unsupported tag = "%s"' % child.tag
        _flux = Flux( label, fluxData )

        xPath.pop( )
        return( _flux )

class Fluxes( suitesModule.Suite ) :
    """
    This class stores a list of :py:class:Flux: instances.
    """

    moniker = 'fluxes'

    def __init__( self ) :

        suitesModule.Suite.__init__( self, [ XYs3d ] )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance of *cls* which must be a sub-class of :py:class:`Fluxes`.

        :param cls:         Form class to return.
        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    Dictionary that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:            An instance of *cls* representing *node*.
        """

        instance = cls()
        instance.parseNode(node, xPath, linkData, **kwargs)

        return instance

    @staticmethod
    def read(fileName, **kwargs):
        """
        Reads in the file at path *fileName* and returns a **Fluxes** instance.

        :param fileName:    Path to the file to read.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is read.

        :returns:           A :py:class:`Fluxes` instance.
        """

        return Fluxes.readXML_file(fileName, **kwargs)

def read(fileName, **kwargs):
    """
    Reads in the file at path *fileName* and returns a **Fluxes** instance.

    :param fileName:    Path to the file to read.
    :param kwargs:      A dictionary of extra arguments that controls how *self* is read.

    :returns:           A :py:class:`Fluxes` instance.
    """

    return Fluxes.read(fileName, **kwargs)
