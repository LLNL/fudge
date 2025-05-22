# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>


"""
This module contains the following classes:

    +---------------------------------------+-----------------------------------------------------------------------------------+
    | Class                                 | Description                                                                       |
    +=======================================+===================================================================================+
    | Particle                              |                                                                                   |
    +---------------------------------------+-----------------------------------------------------------------------------------+
    | Particles                             |                                                                                   |
    +---------------------------------------+-----------------------------------------------------------------------------------+
    | Flux                                  |                                                                                   |
    +---------------------------------------+-----------------------------------------------------------------------------------+
    | Fluxes                                |                                                                                   |
    +---------------------------------------+-----------------------------------------------------------------------------------+
    | ProcessedFlux                         |                                                                                   |
    +---------------------------------------+-----------------------------------------------------------------------------------+
"""


import re
import numpy

from fudge import styles as stylesModule
from fudge import enums as fudgeEnumsModule

from fudge.processing import flux as fluxModule
from fudge.processing import transporting as transportingModule
from fudge.processing.deterministic import multiGroupSettings as settingsModule

from fudge import suites as suitesModule
from fudge import reactionSuite as reactionSuiteModule

from LUPY import ancestry as ancestryModule

from xData import enums as xDataEnumsModule
from xData import axes as axesModule
from xData import values as valuesModule
from xData import vector as vectorModule
from xData import matrix as matrixModule
from xData import gridded as griddedModule
from xData import xDataArray as arrayModule
from xData import multiD_XYs as multiD_XYsModule

from PoPs import IDs as popsIDsModule

class Particle(ancestryModule.AncestryIO):
    """Specifies particle information as mainly needed for multi-group transport."""

    def __init__(self, particleInstanceOrJustID, multiGroup=None, mode=None, fluxes=None):
        """
        Constructor for Particle class.
        """

        self.conserve = None  # conserve flag is set by Particles.process()

        if isinstance(particleInstanceOrJustID, Particle):
            self.__particleID = particleInstanceOrJustID.particleID
            self.multiGroup = particleInstanceOrJustID.multiGroup
            self.fineMultiGroup = particleInstanceOrJustID.fineMultiGroup
            self.fluxes = particleInstanceOrJustID.fluxes
            self.mode = particleInstanceOrJustID.mode
            self.collapseIndices = particleInstanceOrJustID.collapseIndices
            self.processedFluxes = particleInstanceOrJustID.processedFluxes

        elif isinstance(particleInstanceOrJustID, str):
            if not settingsModule.Settings.isAllowedParticle(particleInstanceOrJustID):
                raise TypeError(f'First argument is expected to be one of the common particles in {popsIDsModule}.')

            self.__particleID = particleInstanceOrJustID

            if multiGroup is not None and not isinstance(multiGroup, transportingModule.MultiGroup):
                raise TypeError(f'The second argument expected to be of type {transportingModule.MultiGroup}.')

            self.multiGroup = multiGroup
            self.fineMultiGroup = None
            
            self.collapseIndices = []
            self.processedFluxes = []

            if mode is not None and not isinstance(mode, transportingModule.Mode):
                raise TypeError(f'The third argument is expected to be of type {transportingModule.Mode}')

            self.mode = mode

            if fluxes is not None and not isinstance(fluxes, Fluxes):
                raise TypeError(f'The fourth argument expected to be of type {Fluxes}.')

            self.fluxes = []
            if fluxes is not None:
                for flux in fluxes:
                    if not isinstance(flux, Flux):
                        raise TypeError(f'The entries in the third argument are expected to be of type {Flux}.')

                    self.appendFlux(flux)

    @property
    def particleID(self):
        """
        Returns the particle ID.
        """

        return self.__particleID

    @property
    def label(self):
        """
        Returns the particle ID.
        """

        return self.__particleID

    def moniker(self):
        pass

    def parseNodeUsingClass(self):
        pass

    def toXML_strList(self):
        pass

    def multiGroupIndexFromEnergy(self, energyValue, encloseOutOfRange=False):
        """
        Return the multi-group index whose boundaries enclose the energy argument.

        :param energyValue: The energy of the whose index is to be returned.
        :param encloseOutOfRange: Determines the action if energy is below or above the domain of the boundaries.
        """
        
        return self.multiGroup.multiGroupIndexFromEnergy(energyValue, encloseOutOfRange)

    def numberOfGroups(self):
        """
        Returns the number of multi-groups.
        """

        return self.multiGroup.numberOfGroups

    def appendFlux(self, flux):
        """
        Adds a flux at a specified temperature to the list of fluxes. Currently, the temperature of the flux must be greater than temperatures for the currently listed fluxes.

        :param flux: The flux to add.
        """

        if not isinstance(flux, Flux):
            raise TypeError(f'The argument to the appendFlux method is not of the type {Flux}')

        if len(self.fluxes) == 0:
            self.fluxes.append(flux)
        else:
            if any([flux.temperature<self.fluxes[i].temperature for i in range(len(self.fluxes))]):
                raise ValueError('The temperature of the flux provided as input argument must be greater than the temperatures of the currently listed fluxes.')

            self.fluxes.append(flux)

    def appendProcessedFlux(self, processedFlux):
        """
        Adds a processed flux at a specified temperature to the list of fluxes. Currently, the temperature of the flux must be greater than temperatures for the currently listed fluxes.

        :param processedFlux: The processed flux to add.
        """

        if not isinstance(processedFlux, ProcessedFlux):
            raise TypeError(f'The argument to the appendFlux method is not of the type {ProcessedFlux}')

        if len(self.processedFluxes) == 0:
            self.processedFluxes.append(processedFlux)
        else:
            if any([processedFlux.temperature<self.processedFluxes[i].temperature for i in range(len(self.processedFluxes))]):
                raise ValueError('The temperature of the processed flux provided as input argument must be greater than the temperatures of the currently listed processed fluxes.')

            self.processedFluxes.append(processedFlux)

    def nearestProcessedFluxToTemperature(self, temperature):
        """
        Returns the multi-group flux with its temperature closes to the input argument.

        :param temperature: The temperature of the desired flux.
        """

        availableTemperatures = numpy.array([processedFlux.temperature for processedFlux in self.processedFluxes])
        nearestIndex = numpy.abs(availableTemperatures - temperature).argmin()

        return self.processedFluxes[nearestIndex]

    def process(self, fineGroupBoundaries, epsilon=1e-6):
        """
        Determines the mapping of fine multi-group boundaries to coarse multi-group boundaries, and calculated multi-group fluxes from self.fluxes.

        :param fineGroupBoundaries: The fine group boundaries.
        :param epsilon: Specifies how close a coarse multi-group boundary must be to a fine multi-group boundary to be considered the same boundary.
        """

        if not isinstance(fineGroupBoundaries, axesModule.Grid):
            raise TypeError(f'The first argument is expected to of the type {axesModule.Grid}')

        if self.multiGroup is None:
            if self.mode != transportingModule.Mode.MonteCarloContinuousEnergy:
                raise RuntimeError(f'Multi-groundaries not set for particle "{self.particleID}".')

        else:
            fineGroupBoundariesVector = vectorModule.Vector(values=fineGroupBoundaries.values.values)
            if self.fineMultiGroup is not None:
                assert self.fineMultiGroup.size == fineGroupBoundariesVector.size, "Redefining particle's fine multi-group of different size not allowed."

                maxDifference = max([numpy.abs(self.fineMultiGroup[i] - fineGroupBoundariesVector[i]) for i in range(self.fineMultiGroup.size)])
                assert maxDifference < epsilon, "Redefining particle's fine multi-group not allowed."

                return

            else:
                def findMultiGroupInGroupBoundaries(_multiGroupIndex, _boundaryIndex):
                    while _boundaryIndex < fineGroupBoundariesVector.size:
                        _difference = numpy.abs(self.multiGroup[_multiGroupIndex] - fineGroupBoundariesVector[_boundaryIndex])
                        if _difference <= fineGroupBoundariesVector[_boundaryIndex] * epsilon:
                            break

                        _boundaryIndex += 1
                        if _boundaryIndex == fineGroupBoundariesVector.size:
                            raise RuntimeError(f'Multi-group {_multiGroupIndex + 1} is not compatible with the provided fine group boundaries for particle {self.particleID}.')

                    return _boundaryIndex

                boundaryIndex = 0
                for multiGroupIndex in range(self.multiGroup.size):
                    boundaryIndex = findMultiGroupInGroupBoundaries(multiGroupIndex, boundaryIndex)
                    self.collapseIndices.append(boundaryIndex)

                for flux in self.fluxes:
                    processedFlux = flux.process(fineGroupBoundaries)
                    self.appendProcessedFlux(processedFlux)

                self.fineMultiGroup = transportingModule.MultiGroup('', fineGroupBoundariesVector)

class Particles(suitesModule.Suite):
    """ This holds a dictionary of particles, whose keys are the particles' names (i.e., PoPs ID). """

    moniker = 'particles'

    def __init__(self):
        """
        Constructor for Particles class.
        """

        suitesModule.Suite.__init__( self, [ Particle ] )

    def hasParticle(self, particleID):
        """
        Return a boolean indicating if particle with given ID is present

        :param particleID: ID of the particle that is queried. 
        """

        return self.__contains__(particleID)

    def sortedIDs(self):
        """
        Returns a sorted list of particle ID present in suite.
        """

        return sorted(self.__particles.keys())

    def process(self, reactionSuite, styleLabel):
        """
        Process all the particles available in self.

        :param reactionSuite: The reactionSuite whose multi-group data are to accessed.
        :param styleLabel: The label of the multi-group data to process
        """

        if not isinstance(reactionSuite, reactionSuiteModule.ReactionSuite):
            raise TypeError(f'First argument is not of type {reactionSuiteModule.ReactionSuite}.')

        assert styleLabel in reactionSuite.styles, f'Second argument {styleLabel} is not a style label in the reactionSuite provided as the first argument'

        style = reactionSuite.styles[styleLabel]
        if isinstance(style, stylesModule.SnElasticUpScatter):
            style = style.derivedFromStyle

        assert isinstance(style, stylesModule.HeatedMultiGroup), f'The style label {styleLabel} does not yield a heatedMultiGroup style.'

        for particle in self:
            particle.conserve = style.transportables[particle.particleID].conserve
            groupBoundaries = style.transportables[particle.particleID].group.boundaries
            particle.process(groupBoundaries)
    

class Flux(fluxModule.Flux):
    """ The multi-group flux. """

    def __init__(self, labelOrFluxInstance, temperature=None, data=None):
        """
        The constructor for the Flux class.
        """

        if isinstance(labelOrFluxInstance, fluxModule.Flux):
            if temperature is None:
                raise RuntimeError('The flux temperature is required as the second argument.')

            label = labelOrFluxInstance.label
            data = labelOrFluxInstance.data.copy()

        elif isinstance(labelOrFluxInstance, Flux):
            label = labelOrFluxInstance.label
            data = labelOrFluxInstance.data.copy()
            temperature = labelOrFluxInstance.temperature

        else:
            if not isinstance(labelOrFluxInstance, str):
                raise TypeError('First argument is expected to be a string')

            label = labelOrFluxInstance

            if temperature is None:
                raise RuntimeError('The flux temperature is required as the second argument.')

            if data is None or not isinstance(data, multiD_XYsModule.XYs2d):
                raise TypeError(f'The third argument is requires to be of type {multiD_XYsModule.XYs2d}')

            if len(data.axes) != 4:
                raise TypeError(f'The third argument is expected to have an axes attribute of length 4')

        super().__init__(label, data)
        self.temperature = temperature

    def process(self, groupBoundaries):
        """
        Multi-group the flux and return as a ProcessedFlux instance.
        """

        if not isinstance(groupBoundaries, axesModule.Grid):
            raise TypeError(f'The first argument is expected to of the type {axesModule.Grid}')

        axes = self.data.axes.copy()
        axes.axes.pop(-1)
        axes[1] = axesModule.Grid(
            axes[1].label,
            1,
            axes[1].unit,
            style = xDataEnumsModule.GridStyle.parameters,
            values = valuesModule.Values([0], valueType=xDataEnumsModule.ValueType.integer32)
        )
        axes[2] = groupBoundaries.copy()

        legendre0Flux = self.getFluxAtLegendreOrder(0)
        data = legendre0Flux.groupOneFunction(groupBoundaries)
        dataSize = len(data)
        starts = valuesModule.Values([0], valueType=xDataEnumsModule.ValueType.integer32)
        lengths = valuesModule.Values([dataSize], valueType=xDataEnumsModule.ValueType.integer32)
        array = arrayModule.Flattened( shape=(dataSize, 1), data=data, starts=starts, lengths=lengths )
        data = griddedModule.Gridded2d( axes=axes, array=array )

        return ProcessedFlux(self.temperature, data.array.constructVector(0))


class Fluxes( suitesModule.Suite ) :
    """ Suite of Flux instances. """
    moniker = 'fluxes'

    def __init__( self ) :
        """
        Constructor for Fluxes class.
        """

        suitesModule.suite.__init__( self, [ Flux ] )


class ProcessedFlux():
    """ The multi-group collpased flux. """

    def __init__(self, temperatureOrFluxInstance, fluxVector=None):
        """
        Constructor for ProcessedFlux class.
        """

        if isinstance(temperatureOrFluxInstance, ProcessedFlux):
            self.temperature = temperatureOrFluxInstance.temperature
            self.multiGroupFlux = temperatureOrFluxInstance.multiGroupFlux

        elif isinstance(temperatureOrFluxInstance, float):
            if not isinstance(fluxVector, vectorModule.Vector):
                raise TypeError(f'The second argument is expected to be of type {vectorModule.Vector}')

            self.temperature = temperatureOrFluxInstance
            self.multiGroupFlux = fluxVector



def CollapseMatrix(preCollapsedMatrix, multiGroupSettings, particles, temperature, productID):
    """
    Collapse and return a multi-group matrix.

    :param preCollapsedMatrix: The matric to collapse.
    :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
    :param particles: The list of particles to be transported.
    :param temperatures: The temperature of the flux to use when collapsing.
    :param productID: Particle id of the outgoing particle.
    """

    if not isinstance(preCollapsedMatrix, matrixModule.Matrix):
        raise TypeError(f'The first argument is expected to of the type {matrixModule.Matrix}')

    if preCollapsedMatrix.shape == 0:
        return multiGroupSettings.multiGroupZeroMatrix(particles, productID, collapse=True)

    else:
        projectile = particles[multiGroupSettings.projectileID]
        processedFlux = projectile.nearestProcessedFluxToTemperature(temperature).multiGroupFlux
        projectileCollapsedIndices = projectile.collapseIndices

        product = particles[productID]
        productCollapsedIndices = product.collapseIndices
        productCollapsedIndices[0] = 0
        productCollapsedIndices[-1] = preCollapsedMatrix.numberOfColumns

        # weight by 1 to conserve number, or weight by mean E' to conserve energy:
        conserveEnergy = particles[productID].conserve is fudgeEnumsModule.Conserve.energyOut
        if conserveEnergy:
            gbs = particles[productID].fineMultiGroup.boundaries
            unitWeight = [0.5 * (gbs[idx] + gbs[idx+1]) for idx in range(len(gbs)-1)]
        else:
            unitWeight = [1] * preCollapsedMatrix.numberOfRows

        productCollapsed = matrixModule.Matrix()
        for rowIndex in range(preCollapsedMatrix.numberOfRows):
            productCollapsedRow = collapseRow(preCollapsedMatrix[rowIndex, :], productCollapsedIndices, unitWeight, returnNormalized=False)
            productCollapsed.append(productCollapsedRow)

        projectileCollapsed = matrixModule.Matrix()
        for columnIndex in range(productCollapsed.numberOfColumns):
            projectileCollapsedColumn = collapseRow(productCollapsed[:, columnIndex], projectileCollapsedIndices, processedFlux, returnNormalized=True)
            projectileCollapsed.append(projectileCollapsedColumn)

        if conserveEnergy:
            # re-weight by average outgoing energy in new bins:
            gbs = particles[productID].multiGroup.boundaries
            newWeight = [0.5 * (gbs[idx] + gbs[idx+1]) for idx in range(len(gbs)-1)]
            for idx, row in enumerate(projectileCollapsed.matrix):
                row /= newWeight[idx]

        projectileCollapsed = projectileCollapsed.transpose()

        return projectileCollapsed


def collapseVector(preCollapsedVector, multiGroupSettings, particles, temperature):
    """
    Collapse and return a multi-group vector.

    :param preCollapsedMatrix: The vector to collapse.
    :param multiGroupSettings: Object instance to instruct deterministic methods on what data are being requested.
    :param particles: The list of particles to be transported.
    :param temperature: The temperature of the flux to use when collapsing.
    """
    if not isinstance(preCollapsedVector, vectorModule.Vector):
        raise TypeError(f'The first argument is expected to of the type {vectorModule.Vector}')

    projectile = particles[multiGroupSettings.projectileID]
    processedFlux = projectile.nearestProcessedFluxToTemperature(temperature).multiGroupFlux
    projectileCollapsedIndices = projectile.collapseIndices

    return collapseRow(preCollapsedVector, projectileCollapsedIndices, processedFlux, returnNormalized=True)


def collapseRow(preCollapsedVector, collapseIndices, collapseWeight, returnNormalized):
    """
    Collapse and return a multi-group vector.

    :param preCollapsedVector: The vector to collapse.
    :param collapseIndices: Maps uncollapsed indices to collapsed indices.
    :param collapseWeight: The uncollapsed flux weighting.
    :param returnNormalized: Option to normalize the collapsed vector sum to unity.
    """

    collapseIndicesSize = len(collapseIndices)
    collapsedVector = vectorModule.Vector( collapseIndicesSize - 1 )

    if len(preCollapsedVector) > 0:
        collapseIndexStart = collapseIndices[0]
        for collapseVectorIndex in range(collapsedVector.size):
            fluxSum = 0
            valueSum = 0
            collapseIndexStop = collapseIndices[collapseVectorIndex+1]
            for collapseIndex in range(collapseIndexStart, collapseIndexStop):
                fluxSum += collapseWeight[collapseIndex]
                valueSum += collapseWeight[collapseIndex] * preCollapsedVector[collapseIndex]

            if returnNormalized and fluxSum > 0:
                collapsedVector[collapseVectorIndex] = valueSum / fluxSum

            else:
                collapsedVector[collapseVectorIndex] = valueSum

            collapseIndexStart = collapseIndexStop

    return collapsedVector



