# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

__metaclass__ = type

from xData import standards as standardsModule
from xData import ancestry as ancestryModule
from xData import axes as axesModule
from xData import XYs as XYsModule
from xData import series1d as series1dModule
from xData import gridded as griddedModule
from xData import xDataArray as arrayModule
from xData import values as valuesModule
from xData import multiD_XYs as multiD_XYsModule

from fudge import suites as suitesModule

def axes( energyUnit, fluxUnit = '1/s', temperatureUnit = 'MeV/k' ) :

    axes = axesModule.axes( rank = 4 )
    axes[0] = axesModule.axis( 'flux', 0, fluxUnit )
    axes[1] = axesModule.axis( 'mu', 0, '' )
    axes[2] = axesModule.axis( 'energy_in', 0, energyUnit )
    axes[3] = axesModule.axis( 'temperature', 0, temperatureUnit )
    return( axes )

class LegendreSeries( series1dModule.LegendreSeries ) :

    pass

class XYs2d( multiD_XYsModule.XYs2d ) :

    def getFluxAtLegendreOrder( self, order ) :

        axes = axesModule.axes( )
        axes[1] = self.axes[2].copy( [] )
        axes[0] = self.axes[0].copy( [] )
        for LS in self :
            if( len( LS ) > 1 ) : raise Exception( 'FIXME -- next line is wrong.' )
        return( XYsModule.XYs1d( [ [ LS.outerDomainValue, LS[0] ] for LS in self ], axes = axes ) )

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

class gridded2d( griddedModule.gridded2d ) :

    pass

class XYs3d( multiD_XYsModule.XYs3d ) :

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs2d, ) )

class flux( ancestryModule.ancestry ) :
    """
    This class stores the flux as an XYs2d instance.
    """

    moniker = 'flux'

    def __init__( self, label, data ) :

        if( not( isinstance( label, str ) ) ) : raise TypeError( 'label must be a string instance.' )
        self.label = label

        if( not( isinstance( data, ( XYs2d, ) ) ) ) : raise TypeError( 'Invalid flux data.' )
        self.data = data
        self.data.setAncestor( self )

    def getFluxAtLegendreOrder( self, order ) :

        return( self.data.getFluxAtLegendreOrder( order ) )

    def processMultiGroup( self, style, tempInfo, indent ) :

        reactionSuite = tempInfo['reactionSuite']
        axes = self.data.axes.copy( )
        axes[1] = axesModule.grid( axes[1].label, 1, axes[1].unit, style = axesModule.parametersGridToken,
                values = valuesModule.values( [ 0 ], valueType = standardsModule.types.integer32Token ) )
        axes[2] = style.transportables[reactionSuite.projectile].group.boundaries.copy( [] )
        data = self.data.processMultiGroup( style, tempInfo, indent )
        starts = valuesModule.values( [ 0 ], valueType = standardsModule.types.integer32Token )
        lengths = valuesModule.values( [ len( data ) ], valueType = standardsModule.types.integer32Token )
        array = arrayModule.flattened( shape = ( len( data ), 1 ), data = data, starts = starts, lengths = lengths )
        data = gridded2d( axes = axes, array = array )
        return( data )

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent = indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s label="%s">' % ( indent, self.moniker, self.label ) ]
        xmlStringList += self.data.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )

        label = element.get( 'label' )
        for child in element :
            if( child.tag == XYs2d.moniker ) :
                fluxData = XYs2d.parseXMLNode( child, xPath, linkData )
            else :
                raise 'Unsupported tag = "%s"' % child.tag
        _flux = flux( label, fluxData )

        xPath.pop( )
        return( _flux )

class fluxes( suitesModule.suite ) :

    moniker = 'fluxes'

    def __init__( self ) :

        suitesModule.suite.__init__( self, [ XYs3d ] )
