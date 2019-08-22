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

__metaclass__ = type

from xData import standards as standardsModule
from xData import ancestry as ancestryModule
from xData import axes as axesModule
from xData import XYs as XYsModule
from xData import series1d as series1dModule
from xData import gridded as griddedModule
from xData import array as arrayModule
from xData import values as valuesModule
from xData import multiD_XYs as multiD_XYsModule

class LegendreSeries( series1dModule.LegendreSeries ) :

    pass

class XYs2d( multiD_XYsModule.XYs2d ) :

    def getFluxAtLegendreOrder( self, order ) :

        axes = axesModule.axes( )
        axes[1] = self.axes[2].copy( [] )
        axes[0] = self.axes[0].copy( [] )
        for LS in self :
            if( len( LS ) > 1 ) : raise Exception( 'FIXME -- next line is wrong.' )
        return( XYsModule.XYs1d( [ [ LS.value, LS[0] ] for LS in self ], axes = axes ) )

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
