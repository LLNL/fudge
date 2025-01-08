# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from LUPY import enums as LUPY_enumsModule
import numpy

from fudge.processing import group as groupModule
from fudge.processing import transportables as transportablesModule 

from xData import axes as axesModule
from xData import vector as vectorModule

class DelayedNeutrons(LUPY_enumsModule.Enum):
    """Custom enum class to indicated whether delayed neutrons are to be included in teh transport"""

    on  = LUPY_enumsModule.auto()
    off = LUPY_enumsModule.auto()


class Mode(LUPY_enumsModule.Enum):
    """Custom enum class to indicating a transport mode of either multiGroup, multiGroupWithSnElasticUpScatter or MonteCarloContinuousEnergy. """

    multiGroup = LUPY_enumsModule.auto()
    multiGroupWithSnElasticUpScatter = LUPY_enumsModule.auto()
    MonteCarloContinuousEnergy = LUPY_enumsModule.auto()


class TransportCorrectionType(LUPY_enumsModule.Enum):
    """ Custom enum class indicating the transport correction type of none or Pendlebury."""

    none = LUPY_enumsModule.auto()
    Pendlebury = LUPY_enumsModule.auto()


class MultiGroup():
    """ Class for multi-group data labels and boundaries. """

    aboveHighestBoundary = -1
    belowLowestBoundary  = -2
    noGroupBoundaries    = -3

    def __init__(self, labelOrMultiGroupInstance, boundaries=None):
        """
        Constructor for MultiGroup class.
        """

        if isinstance(labelOrMultiGroupInstance, MultiGroup):
            self.__label = labelOrMultiGroupInstance.label
            self.__boundaries = labelOrMultiGroupInstance.boundaries.copy()

        elif isinstance(labelOrMultiGroupInstance, groupModule.Group):
            self.__label = labelOrMultiGroupInstance.label
            self.__boundaries = vectorModule.Vector(values=labelOrMultiGroupInstance.boundaries.values.values)
        
        elif isinstance(labelOrMultiGroupInstance, str):
            self.__label = labelOrMultiGroupInstance

            if not isinstance(boundaries, (list, tuple, vectorModule.Vector)):
                raise TypeError(f'The first argument is expected to of the type {list}, {tuple} or {vectorModule.Vector}')

            self.__boundaries = vectorModule.Vector(values=boundaries)

        else:
            raise TypeError(f'First argument is expected to be of type {str}, {MultiGroup} or {groupModule.group}')

        if not self.boundaries.isSorted:
            raise ValueError( 'Group boundaries must be sorted in increasing order!' )

    def __getitem__(self, index):
        """
        Return value of multi-group boundary at location specified by the input argument.

        :index: Index of the self.boundaries values to return.
        :returns: scalar entry from self.boundaries.
        """        

        return self.boundaries[index]

    def __len__(self):
        """
        Returns the size of self.boundaries.
        """

        return len(self.boundaries)

    @property
    def label(self):
        """
        Return the multi-group label.
        """

        return self.__label

    @property
    def size(self):
        """
        Returns the size of self.boundaries.
        """

        return self.boundaries.size

    @property
    def numberOfGroups(self):
        """
        Returns the number of multi-groups.
        """

        return self.size - 1

    @property
    def boundaries(self):
        """
        Returns the multi-group boundaries.
        """
        
        return self.__boundaries

    def multiGroupIndexFromEnergy(self, energyValue, encloseOutOfRange=False):
        """
        Return the multi-group index whose boundaries enclose the energy argument.

        :param energyValue: The energy of the whose index is to be returned.
        :param encloseOutOfRange: Determines the action if energy is below or above the domain of the boundaries.
        """

        if self.size == 0:
            indexValue = self.noGroupBoundaries

        elif energyValue < self.boundaries[0]:
            indexValue =  0 if encloseOutOfRange else self.belowLowestBoundary

        elif energyValue >= self.boundaries[-1]:
            indexValue =  self.size-2 if encloseOutOfRange else self.aboveHighestBoundary

        else:
            indexValue = numpy.nonzero(energyValue < self.boundaries.vector)[0]
            if indexValue.shape[0] > 0:
                indexValue = indexValue[0] - 1
            else:
                indexValue = self.aboveHighestBoundary

        return indexValue





