# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

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

    axes = axesModule.Axes(4)
    axes[0] = axesModule.Axis( 'flux', 0, fluxUnit )
    axes[1] = axesModule.Axis( 'mu', 0, '' )
    axes[2] = axesModule.Axis( 'energy_in', 0, energyUnit )
    axes[3] = axesModule.Axis( 'temperature', 0, temperatureUnit )
    return( axes )

class LegendreSeries( series1dModule.LegendreSeries ) :

    pass

class XYs2d( multiD_XYsModule.XYs2d ) :

    def getFluxAtLegendreOrder( self, order ) :

        axes = axesModule.Axes(2)
        axes[1] = self.axes[2].copy()
        axes[0] = self.axes[0].copy()
        for LS in self :
            if( len( LS ) > 1 ) : raise Exception( 'FIXME -- next line is wrong.' )
        return( XYs1dModule.XYs1d( [ [ LS.outerDomainValue, LS[0] ] for LS in self ], axes = axes ) )

    def processMultiGroup( self, style, tempInfo, indent ) :

        # BRB, Currently only processes LegendreSeries for 1d data and only its l=0 term. Need to convert units if needed.
        # The return instance is a list, this needs to be changed.

        from fudge.processing import miscellaneous as miscellaneousModule

        flux0 = self.getFluxAtLegendreOrder( 0 )
        projectileName = tempInfo['reactionSuite'].projectile
        return( flux0.groupOneFunction( style.transportables[projectileName].group.boundaries ) )

    @staticmethod
    def allowedSubElements( ) :

        return( ( LegendreSeries, ) )

class Gridded2d( griddedModule.Gridded2d ) :

    pass

class XYs3d( multiD_XYsModule.XYs3d ) :

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs2d, ) )

class Flux( ancestryModule.AncestryIO ) :
    """
    This class stores the flux as an XYs2d instance.
    """

    moniker = 'flux'

    def __init__( self, label, data ) :

        ancestryModule.AncestryIO.__init__(self)

        if( not( isinstance( label, str ) ) ) : raise TypeError( 'label must be a string instance.' )
        self.label = label

        if( not( isinstance( data, ( XYs2d, ) ) ) ) : raise TypeError( 'Invalid flux data.' )
        self.data = data
        self.data.setAncestor( self )

    def convertUnits(self, unitMap):
        '''
        Converts unit per *unitMap*.
        '''

        self.data.convertUnits(unitMap)

    def getFluxAtLegendreOrder( self, order ) :

        return( self.data.getFluxAtLegendreOrder( order ) )

    def processMultiGroup( self, style, tempInfo, indent ) :

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

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s label="%s">' % ( indent, self.moniker, self.label ) ]
        xmlStringList += self.data.toXML_strList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

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

    moniker = 'fluxes'

    def __init__( self ) :

        suitesModule.Suite.__init__( self, [ XYs3d ] )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        instance = cls()
        instance.parseNode(node, xPath, linkData, **kwargs)

        return instance

    @staticmethod
    def read(fileName, **kwargs):
        """
        Reads in the file name *fileName* and returns a **Fluxes** instance.
        """

        return Fluxes.readXML_file(fileName, **kwargs)

def read(fileName, **kwargs):
    """
    Reads in the file name *fileName* and returns a **Fluxes** instance.
    """

    return Fluxes.read(fileName, **kwargs)
