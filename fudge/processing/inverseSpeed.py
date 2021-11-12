# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

__metaclass__ = type

import math

from pqu import PQU as PQUModule

from PoPs import IDs as IDsPoPsModule

from xData import ancestry as ancestryModule
from xData import axes as axesModule
from xData import gridded as griddedModule

from . import group as groupModule

from . import miscellaneous as miscellaneousModule

class gridded1d( griddedModule.gridded1d ) :

    def __init__( self, **kwargs ) :

        griddedModule.gridded1d.__init__( self, **kwargs )

    def convertUnits( self, unitMap ) :

        if( 'sh' in self.axes[0].unit ) : return
        griddedModule.gridded1d.convertUnits( self, unitMap )

class inverseSpeed( ancestryModule.ancestry ) :
    """
    This class stores the inverse speed.
    """

    moniker = 'inverseSpeed'

    def __init__( self, data ) :

        if( not( isinstance( data, ( gridded1d, ) ) ) ) : raise TypeError( 'Invalid flux data.' )
        self.data = data
        self.data.setAncestor( self )

    def convertUnits( self, unitMap ) :

        self.data.convertUnits( unitMap )

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent = indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        xmlStringList += self.data.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )

        for child in element :
            if( child.tag == gridded1d.moniker ) :
                data = gridded1d.parseXMLNode( child, xPath, linkData )
            else :
                raise 'Unsupported tag = "%s"' % child.tag
        _inverseSpeed = inverseSpeed( data )

        xPath.pop( )
        return( _inverseSpeed )

def multiGroupInverseSpeed( style, tempInfo ) :

    _groupBoundaries, flux = miscellaneousModule._groupFunctionsAndFluxInit( style, tempInfo, None )
    groupBoundaries = _groupBoundaries.values.values
    numberOfGroups = len( groupBoundaries ) - 1

    speedOfLight_m_mus = PQUModule.PQU( 1, 'c' ).getValueAs( 'm/mus' )

    if( tempInfo['projectile'].id == IDsPoPsModule.photon ) :
        inverseSpeeds = numberOfGroups * [ 1.0 / speedOfLight_m_mus ]
    else :
        inverseSpeeds = []
        projectileMass = tempInfo['projectileMass']
        relativisticTreatment = tempInfo.get('relativisticTreatment', 40 * groupBoundaries[-1] > projectileMass) # Treat legacy neutron boundaries as non-relativistic (i.e., 40 * 20 MeV is not greater than neutron mass energy).
        if relativisticTreatment:
            lowerRatio = groupBoundaries[0] / projectileMass
            lowerInverseSpeedValue = math.sqrt( lowerRatio * ( 2.0 + lowerRatio ) )
            for i1 in range( numberOfGroups ) :         # Note, currently the flux is not included in the group average.
                upperRatio = groupBoundaries[i1+1] / projectileMass
                upperInverseSpeedValue = math.sqrt( upperRatio * ( 2.0 + upperRatio ) )
                inverseSpeed = ( upperInverseSpeedValue - lowerInverseSpeedValue ) / ( upperRatio - lowerRatio ) / speedOfLight_m_mus
                if( i1 == 0 ) : print( '    ', groupBoundaries[i1+1], upperRatio, upperInverseSpeedValue, inverseSpeed, speedOfLight_m_mus )
                inverseSpeeds.append( inverseSpeed )
                lowerRatio = upperRatio
                lowerInverseSpeedValue = upperInverseSpeedValue
        else:
            const = 2. / ( 3. * math.sqrt( 2 ) * speedOfLight_m_mus )
            const *= math.sqrt( float( PQUModule.PQU( 1, 'MeV' ).convertToUnit( _groupBoundaries.unit ) ) )
            const *= math.sqrt( projectileMass )                    # Mass is in energy unit / c^2.

            E1 = groupBoundaries[0]

            for i1 in range( numberOfGroups ) :
                E2 = groupBoundaries[i1+1]
                fluxSlice = flux.domainSlice( domainMin = E1, domainMax = E2, fill = True )
                fluxIntegral = 0
                inverseSpeedIntegral = 0
                _E1 = None
                for _E2, _flux2 in fluxSlice :
                    if( _E1 is not None ) :
                        fluxIntegral += ( _E2 - _E1 ) * ( _flux2 + _flux1 )
                        x1 = math.sqrt( _E1 )
                        x2 = math.sqrt( _E2 )
                        inverseSpeedIntegral += ( x2 - x1 ) * ( _flux1 + _flux2 + ( _flux1 * x2 + _flux2 * x1 ) / ( x1 + x2 ) )
                    _E1 = _E2
                    _flux1 = _flux2
                inverseSpeed = 0
                if( fluxIntegral != 0 ) : inverseSpeed = const * inverseSpeedIntegral / ( fluxIntegral / 2 )
                inverseSpeeds.append( inverseSpeed )
                E1 = E2

    axes = axesModule.axes( rank = 2 )
    axes[0] = axesModule.axis( label = "inverse speed", unit = "mus/m", index = 0 )     # 'mus/m', which is a standard, is the same as 'sh/cm'.
    axes[1] = flux.axes[1]
    return( groupModule.toMultiGroup1d( gridded1d, style, tempInfo, axes, inverseSpeeds, addLabel = False, zeroPerTNSL = False ) )
