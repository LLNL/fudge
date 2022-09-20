# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import numpy

from LUPY import ancestry as ancestryModule

from xData import enums as xDataEnumsModule
from xData import axes as axesModule
from xData import values as valuesModule
from xData import xDataArray as arrayModule

from fudge import enums as enumsModule
from fudge import suites as suitesModule

class Group( ancestryModule.AncestryIO ) :
    """
    This class stores the multi-group information.
    """

    moniker = 'group'

    def __init__( self, label, boundaries ) :

        ancestryModule.AncestryIO.__init__( self )
        if( not( isinstance( label, str ) ) ) : raise TypeError( 'label must only be a string instance.' )
        self.__label = label

        if( not( isinstance( boundaries, axesModule.Grid ) ) ) : raise TypeError( 'Group boundaries must only be a grid instance.' )
        if list(boundaries.values) != sorted(boundaries.values):
            raise ValueError( 'Group boundaries must be sorted in increasing order!' )
        self.__boundaries = boundaries

    @property
    def label( self ) :

        return( self.__label )

    @property
    def boundaries( self ) :

        return( self.__boundaries )

    def copy( self ) :

        return Group(self.label, self.boundaries.copy())

    __copy__ = copy

    def convertUnits( self, unitMap ) :
        """
        unitMap is a dictionary of the form { 'eV' : 'MeV', 'b' : 'mb' }.
        """

        self.__boundaries.convertUnits( unitMap )

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s label="%s">' % ( indent, self.moniker, self.label ) ]
        xmlStringList += self.boundaries.toXML_strList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        xPath.append(node.tag)

        boundaries = axesModule.Grid.parseNodeUsingClass(node.find('grid'), xPath, linkData, **kwargs)

        xPath.pop()

        return cls(node.get( 'label' ), boundaries)

def toMultiGroup1d( cls, style, tempInfo, _axes, data, addLabel = True, zeroPerTNSL = True ) :
    """
    This function takes 1-d multi-group data as a list of floats and return to 1-d gridded instance
    containing the multi-group data. The list of data must contain n values where n is the
    number of groups.
    """

    reactionSuite = tempInfo['reactionSuite']

    axes = axesModule.Axes(2)
    axes[0] = axesModule.Axis( _axes[0].label, 0, _axes[0].unit )
    energyInGrid = style.transportables[reactionSuite.projectile].group.boundaries.values.copy( )
    axes[1] = axesModule.Grid( _axes[1].label, 1, _axes[1].unit, xDataEnumsModule.GridStyle.boundaries, energyInGrid )
    shape = [ len( data ) ]
    if( zeroPerTNSL ) : data[tempInfo['maximumProjectileGroupIndex']] = 0.0     # For TNSL data as data ends in middle of this group.
    data = valuesModule.Values( data )
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
    starts = valuesModule.Values([start], valueType=xDataEnumsModule.ValueType.integer32)
    lengths = valuesModule.Values([len(data)], valueType=xDataEnumsModule.ValueType.integer32)
    flattened = arrayModule.Flattened( shape = shape, data = data, starts = starts, lengths = lengths )
    label = None
    if( addLabel ) : label = style.label
    return( cls( label = label, axes = axes, array = flattened ) )

def TMs2Form( style, tempInfo, TM_1, TM_E, productName = None ) :

    from fudge.processing import transportables as transportablesModule
    from fudge.productData.distributions import multiGroup as multiGroupModule

    reactionSuite = tempInfo['reactionSuite']
    if( productName is None ) : productName = tempInfo['productName']
    transportable = style.transportables[productName]
    conserve = transportable.conserve
    energyUnit = tempInfo['incidentEnergyUnit']

# BRB hardwired.
    crossSectionUnit = 'b'                          # ?????? 'b' should not be hardwired.
    axes = axesModule.Axes(4)
    axes[0] = axesModule.Axis( 'C_l(energy_in,energy_out)', 0, crossSectionUnit )
    lMaxPlus1 = len( TM_1[0][0] )
    lGrid = valuesModule.Values([i1 for i1 in range(lMaxPlus1)], valueType=xDataEnumsModule.ValueType.integer32)
    axes[1] = axesModule.Grid( 'l',                         1,         '', xDataEnumsModule.GridStyle.parameters, lGrid )
    energyOutGrid = style.transportables[productName].group.boundaries.values.copy( )
    axes[2] = axesModule.Grid( 'energy_out',                2, energyUnit, xDataEnumsModule.GridStyle.boundaries, energyOutGrid )
    energyInGrid = style.transportables[reactionSuite.projectile].group.boundaries.values.copy( )
    axes[3] = axesModule.Grid( 'energy_in',                 3, energyUnit, xDataEnumsModule.GridStyle.boundaries, energyInGrid )
    if conserve == enumsModule.Conserve.number:
        TM = TM_1
    else:
        raise NotImplementedError('Need to implement')

    n1 = len(TM)
    n2 = len(TM[0])
    n3 = len(TM[0][0])

    if isinstance(TM, (list, tuple, numpy.ndarray)):
        TM_array3d = {}
        for i1 in range(n1):
            TM_array2d = {}
            for i2 in range(n2):
                TM_array2d[i2] = list(TM[i1,i2])
            TM_array3d[i1] = TM_array2d
        TM = TM_array3d

    data, starts, lengths = [], [], []
    if( tempInfo['zeroPerTNSL'] ) :
        maximumProjectileGroupIndex = tempInfo['maximumProjectileGroupIndex']
        for key in sorted( TM[maximumProjectileGroupIndex] ) : TM[maximumProjectileGroupIndex][key] = n3 * [ 0 ]
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
    data = valuesModule.Values( data )
    starts = valuesModule.Values(starts, valueType=xDataEnumsModule.ValueType.integer32)
    lengths = valuesModule.Values(lengths, valueType=xDataEnumsModule.ValueType.integer32)
    flattened = arrayModule.Flattened( shape = shape, data = data, starts = starts, lengths = lengths, 
            dataToString = multiGroupModule.Gridded3d.dataToString )
    gridded3d = multiGroupModule.Gridded3d( axes = axes, array = flattened )
    return( multiGroupModule.Form( style.label, xDataEnumsModule.Frame.lab, gridded3d ) )

class Groups( suitesModule.Suite ) :

    moniker = 'groups'

    def __init__( self ) :

        suitesModule.Suite.__init__( self, [ Group ] )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        instance = cls()
        instance.parseNode(node, xPath, linkData, **kwargs)

        return instance

    @staticmethod
    def read(fileName, **kwargs):
        """
        Reads in the file name *fileName* and returns a **Groups** instance.
        """

        return Groups.readXML_file(fileName, **kwargs)

def read(fileName, **kwargs):
    """
    Reads in the file name *fileName* and returns a **Groups** instance.
    """

    return Groups.read(fileName, **kwargs)
