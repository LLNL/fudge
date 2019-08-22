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

from xData import ancestry as ancestryModule
from xData import standards as standardsModule
from xData import axes as axesModule
from xData import values as valuesModule
from xData import array as arrayModule

class group( ancestryModule.ancestry ) :
    """
    This class stores the multi-group infomation.
    """

    moniker = 'group'

    def __init__( self, label, boundaries ) :

        ancestryModule.ancestry.__init__( self )
        if( not( isinstance( label, str ) ) ) : raise TypeError( 'label must only be a string instance.' )
        self.__label = label

        if( not( isinstance( boundaries, axesModule.grid ) ) ) : raise TypeError( 'Group boundaries must only be a grid instance.' )
        self.__boundaries = boundaries

    @property
    def label( self ) :

        return( self.__label )

    @property
    def boundaries( self ) :

        return( self.__boundaries )

    def copy( self ) :

        return( group( self.label, self.boundaries.copy( ) ) )

    __copy__ = copy
    __deepcopy__ = __copy__


    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s label="%s">' % ( indent, self.moniker, self.label ) ]
        xmlStringList += self.boundaries.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )
        boundaries = axesModule.grid.parseXMLNode( element.find('grid'), xPath, linkData )
        xPath.pop()
        return group( element.get( 'label' ), boundaries )

def toMultiGroup1d( cls, style, tempInfo, _axes, data, addLabel = True ) :
    """
    This function takes 1-d multi-group data as a list of floats and return to 1-d gridded instance
    containing the multi-group data. The list of data must contain n values where n is the
    number of groups.
    """

    reactionSuite = tempInfo['reactionSuite']

    axes = axesModule.axes( rank = 2 )
    axes[0] = axesModule.axis( _axes[0].label, 0, _axes[0].unit )
    energyInGrid = style.transportables[reactionSuite.projectile].group.boundaries.values.copy( )
    axes[1] = axesModule.grid( _axes[1].label, 1, _axes[1].unit, axesModule.boundariesGridToken, energyInGrid )
    shape = [ len( data ) ]
    data = valuesModule.values( data )
    for start, datum in enumerate( data ) :
        if( datum != 0 ) : break
    end = 0
    for i1, datum in enumerate( data ) :
        if( datum != 0 ) : end = i1
    if( start > end ) :                     # Only happens when all values are 0.
        start = 0
    else :
        end += 1
    data = data[start:end]
    starts = valuesModule.values( [ start ], valueType = standardsModule.types.integer32Token )
    lengths = valuesModule.values( [ len( data ) ], valueType = standardsModule.types.integer32Token )
    flattened = arrayModule.flattened( shape = shape, data = data, starts = starts, lengths = lengths )
    label = None
    if( addLabel ) : label = style.label
    return( cls( label = label, axes = axes, array = flattened ) )

def TMs2Form( style, tempInfo, TM_1, TM_E ) :

    from fudge.processing import transportables as transportablesModule
    from fudge.gnd.productData.distributions import multiGroup as multiGroupModule

    reactionSuite = tempInfo['reactionSuite']
    productName = tempInfo['productName']
    transportable = style.transportables[productName]
    conserve = transportable.conserve
    energyUnit = tempInfo['incidentEnergyUnit']

# BRB hardwired.
    crossSectionUnit = 'b'                          # ?????? 'b' should not be hardwired.
    axes = axesModule.axes( rank = 4 )
    axes[0] = axesModule.axis( 'C_l(energy_in,energy_out)', 0, crossSectionUnit )
    lMaxPlus1 = len( TM_1[0][0] )
    lGrid = valuesModule.values( [ i1 for i1 in range( lMaxPlus1 ) ], valueType = standardsModule.types.integer32Token )
    axes[1] = axesModule.grid( 'l',                         1,         '', axesModule.parametersGridToken, lGrid )
    energyOutGrid = style.transportables[productName].group.boundaries.values.copy( )
    axes[2] = axesModule.grid( 'energy_out',                2, energyUnit, axesModule.boundariesGridToken, energyOutGrid )
    energyInGrid = style.transportables[reactionSuite.projectile].group.boundaries.values.copy( )
    axes[3] = axesModule.grid( 'energy_in',                 3, energyUnit, axesModule.boundariesGridToken, energyInGrid )
    if( conserve == transportablesModule.conserve.number ) :
        TM = TM_1
    else :
        raise NotImplementedError( 'Need to implement' )

    n1 = len( TM )
    n2 = len( TM[0] )
    n3 = len( TM[0][0] )
    data, starts, lengths = [], [], []
    for i1 in sorted( TM.keys( ) ) :
        TM_i1 = TM[i1]
        start, length = None, 0
        for i2 in sorted( TM_i1.keys( ) ) :
            TM_i1_i2 = TM_i1[i2]
            cellMin = min( TM_i1_i2 )
            cellMax = max( TM_i1_i2 )
            if( ( cellMin != 0 ) or ( cellMax != 0 ) ) :
                if( start is None ) : start = n3 * ( n2 * i1 + i2 )
                length += n3
                data += TM_i1_i2
            else :
                if( start is not None ) : 
                    starts.append( start )
                    lengths.append( length )
                    start, length = None, 0
        if( start is not None ) : 
            starts.append( start )
            lengths.append( length )
    shape = [ n1, n2, n3 ]
    data = valuesModule.values( data )
    starts = valuesModule.values( starts, valueType = standardsModule.types.integer32Token )
    lengths = valuesModule.values( lengths, valueType = standardsModule.types.integer32Token )
    flattened = arrayModule.flattened( shape = shape, data = data, starts = starts, lengths = lengths, 
            dataToString = multiGroupModule.gridded3d.dataToString )
    gridded3d = multiGroupModule.gridded3d( axes = axes, array = flattened )
    return( multiGroupModule.form( style.label, standardsModule.frames.labToken, gridded3d ) )
