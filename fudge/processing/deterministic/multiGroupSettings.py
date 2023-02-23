# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from typing import Type
from fudge import styles as stylesModule
from xData import vector as vectorModule
from xData import matrix as matrixModule

from PoPs import IDs as popsIDsModule

from fudge.processing import transporting as transportingModule
from fudge.processing.deterministic import particles as particlesModule


class Settings:
    """This class is used to instruct deterministic methods and the Monte Carlo API MCGIDI on what data are being requested."""

    def __init__(self, projectileID, delayedNeutrons=transportingModule.DelayedNeutrons.off, raiseOnError=True):
        """Constructor for Settings class."""

        if not isinstance(projectileID, str):
            raise TypeError('First argument is expected to be of type {str}')

#        if not self.isAllowedParticle(projectileID):
#            raise TypeError(f'First argument is expected to be one of the common particles in {popsIDsModule}.')

        self.projectileID = projectileID

        self.delayedNeutrons = transportingModule.DelayedNeutrons.checkEnumOrString(delayedNeutrons)

        self.nuclearPlusCoulombInterferenceOnly = False

        self.raiseOnError = raiseOnError

    @staticmethod
    def isAllowedParticle(particleID):
        """
        Return a boolean indicating whether particleID corresponds to one of the common particles in PoPs.IDs.
        """

        matchPoPsID = [x for x in dir(popsIDsModule) if isinstance(getattr(popsIDsModule, x), str) and particleID == getattr(popsIDsModule, x)]
        assert len(matchPoPsID) < 2

        return len(matchPoPsID) > 0

    def numberOfGroups(self, particles, particleID, collapsed):
        """
        Return the number of multi-groups.

        :param particles: Suite of particle objects.
        :param particleID: Label for the particle for which the number of groups is seeked.
        :param collapsed: Boolean indicating whether the number of collapsed or fine groups are to be returned.
        """

        if not isinstance(particles, particlesModule.Particles):
            raise TypeError(f'First argument is expected to be of type {particlesModule.Particles}')

        particle = particles[particleID]
        numberOfGroups = particle.numberOfGroups if collapsed else particle.fineMultiGroups.numberOfGroups

        return numberOfGroups

    def multiGroupZeroVector(self, particles, collapsed=True):
        """
        Returns a Vector of 0.0's of the proper length for the projectile's multi-group data.

        :param particles: Suite of particle objects.
        :param collapsed: Boolean indicating whether the number of collapsed or fine groups are to be used.
        """

        numberOfGroupsProjectile = self.numberOfGroups(particles, self.projectileID, collapsed)

        return vectorModule.Vector(numberOfGroupsProjectile)

    def multiGroupZeroMatrix(self, particles, productID, collapsed=True):
        """
        Returns a Matrix of 0.0's of the proper length for the projectile's and product's multi-group data.

        :param particles: Suite of particle objects.
        :param particleID: Label of the particle for which the zero matrix is seeked.
        :param collapsed: Boolean indicating whether the number of collapsed or fine groups are to be used.
        """

        numberOfGroupsProjectile = self.numberOfGroups(particles, self.projectileID, collapsed)
        numberOfGroupsProduct = self.numberOfGroups(particles, productID, collapsed)

        return matrixModule.Matrix(numberOfGroupsProjectile, numberOfGroupsProduct)


class MG(Settings):
    """This class is used to instruct deterministic methods on what data are being requested."""

    def __init__(self, projectileID, mode, delayedNeutrons=transportingModule.DelayedNeutrons.off, raiseOnError=True):
        """Constructor for MG class"""

        super().__init__(projectileID, delayedNeutrons, raiseOnError)

        self.mode = transportingModule.Mode.checkEnumOrString(mode)

    def form(self, suite, temperatureInfo):
        """
        Searches the suite for the form style specified matching one in the temperatureInfo argument.
        
        This only works for multi-group data (i.e., multiGroup or multiGroupWithSnElasticUpScatter type data).

        :param suite: The suite to search for the requested form.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :return: Entry from suite argument that correspond to the appropriate label in the temperatureInfo argument.
        """

        if self.mode == transportingModule.Mode.multiGroup:
            styleLabel = temperatureInfo.heatedMultiGroup
        elif self.mode == transportingModule.Mode.multiGroupWithSnElasticUpScatter:
            styleLabel = temperatureInfo.SnElasticUpScatter
        else:
            raise ValueError(f'Expected a mode of either {transportingModule.Mode.multiGroup} or {transportingModule.Mode.multiGroupWithSnElasticUpScatter}.')

        if styleLabel not in suite:
            if self.mode == transportingModule.Mode.multiGroupWithSnElasticUpScatter:
                styleLabel = temperatureInfo.heatedMultiGroup

        if styleLabel not in suite:
            if self.raiseOnError:
                raise KeyError(f'Style Label {styleLabel} not found in suite')
            else:
                return None

        else:
            return suite[styleLabel]
            
    def formAsVector(self, suite, temperatureInfo):
        """
        Convert the output from self.form to a xData.vector.Vector instance.

        :param suite: The suite to search for the requested form.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :return: xData.vector.Vector instance.
        """

        gridded1d = self.form(suite, temperatureInfo)
        if gridded1d is None:
            return vectorModule.Vector()
        else:
            return gridded1d.constructVector()

    def formAsMatrix(self, suite, temperatureInfo, thirdIndex):
        """
        Convert the output from self.form to a xData.matrix.Matrix instance.

        :param suite: The suite to search for the requested form.
        :param temperatureInfo: TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
        :param: Index for the third dimension for self.form 3-D array output
        :return: xData.matrix.Matrix instance.
        """

        gridded3dForm = self.form(suite, temperatureInfo)
        if gridded3dForm is None:
            return matrixModule.Matrix()

        else:
            if len(gridded3dForm.subforms) != 1:
                raise NotImplementedError('Method only implemented for subforms of length 1')

            return gridded3dForm.subforms[0].constructMatrix(thirdIndex)
