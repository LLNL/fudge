# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import math

from LUPY import ancestry as ancestryModule

from pqu import PQU as PQUModule

from PoPs import IDs as IDsPoPsModule

from xData import axes as axesModule
from xData import gridded as griddedModule

from . import group as groupModule

from . import miscellaneous as miscellaneousModule

class Gridded1d( griddedModule.Gridded1d ) :

    def convertUnits( self, unitMap ) :

        if( 'sh' in self.axes[0].unit ) : return
        griddedModule.Gridded1d.convertUnits( self, unitMap )

class InverseSpeed( ancestryModule.AncestryIO ) :
    """
    This class stores the inverse speed.
    """

    moniker = 'inverseSpeed'

    def __init__( self, data ) :

        ancestryModule.AncestryIO.__init__(self)

        if( not( isinstance( data, ( Gridded1d, ) ) ) ) : raise TypeError( 'Invalid flux data.' )
        self.data = data
        self.data.setAncestor( self )

    def convertUnits( self, unitMap ) :

        self.data.convertUnits( unitMap )

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        xmlStringList += self.data.toXML_strList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )

        for child in element :
            if( child.tag == Gridded1d.moniker ) :
                data = Gridded1d.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            else :
                raise 'Unsupported tag = "%s"' % child.tag
        _inverseSpeed = cls( data )

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
                inverseSpeeds.append( inverseSpeed )
                lowerRatio = upperRatio
                lowerInverseSpeedValue = upperInverseSpeedValue
        else:
            const = 2. / ( 3. * math.sqrt( 2 ) * speedOfLight_m_mus )
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

    axes = axesModule.Axes(2)
    axes[0] = axesModule.Axis( label = "inverse speed", unit = "mus/m", index = 0 )     # 'mus/m', which is a standard, is the same as 'sh/cm'.
    axes[1] = flux.axes[1]
    return( groupModule.toMultiGroup1d( Gridded1d, style, tempInfo, axes, inverseSpeeds, addLabel = False, zeroPerTNSL = False ) )
