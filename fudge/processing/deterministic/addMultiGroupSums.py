# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
import numpy

from xData import enums as xDataEnumsModule
from xData import values as valuesModule
from xData import xDataArray as arrayModule

from PoPs import IDs as IDsModule

from fudge import enums as enumsModule
from fudge import institution as institutionModule
from fudge import outputChannel as outputChannelModule
from fudge import product as productModule
from fudge.reactions import reaction as reactionModule
from fudge.reactionData import crossSection as crossSectionModule
from fudge.reactionData import availableEnergy as availableEnergyModule 
from fudge.reactionData import availableMomentum as availableMomentumModule

from fudge.outputChannelData import Q as QModule
from fudge.outputChannelData.fissionFragmentData import delayedNeutron as delayedNeutronModule
from fudge.productData import multiplicity as multiplicityModule
from fudge.productData import averageProductEnergy as averageProductEnergyModule
from fudge.productData import averageProductMomentum as averageProductMomentumModule
from fudge.productData.distributions import multiGroup as multiGroupModule

from fudge.processing import transporting as transportingModule
from fudge.processing.deterministic import multiGroupSettings as multiGroupSettingsModule

from . import tokens as tokensModule

transportableIDs = (IDsModule.neutron, 'H1', 'H2', 'H3', 'He3', 'He4', IDsModule.photon)

def toGridded1d(module, axesTemplate1d, data, label):
    """
    Creates a Gridded1d instance from the input arguments.

    :param module:              The FUDGE module containing the Gridded1d instance to create.
    :param axesTemplate1d:      An axes templated used to aid in constructing the axes.
    :param data:                The data.
    :param label:               The style label for the data.
    """

    for index, value in enumerate(data):
        if value > 0:
            break
    if value == 0:
        index += 1

    starts = valuesModule.Values([index], valueType=xDataEnumsModule.ValueType.integer32)
    lengths = valuesModule.Values([len(data)-index], valueType=xDataEnumsModule.ValueType.integer32)
    array = arrayModule.Flattened(shape=[len(data)], starts=starts, lengths=lengths, data=valuesModule.Values(data[index:]))
    axes = module.defaultAxes(axesTemplate1d[1].unit)
    axes[1] = axesTemplate1d[1].copy()

    return module.Gridded1d(axes=axes, array=array, label=label)

def toGridded3d(axesTemplate3d, data, label):
    """Creates a Gridded3d instance from the input arguments.

    :param axesTemplate3d:      An axes templated used to aid in constructing the axes.
    :param data:                The data.
    :param label:               The style label for the data.
    """

    ndarray = data.array

    if ndarray.shape[2] > 1:
        if not numpy.any(ndarray[:,:,1:]):
            ndarray = ndarray[:,:,:1]

    starts = []
    lengths = []
    flattenedData = []
    if len(ndarray) > 0:
        n2 = len(ndarray[0])
        n3 = len(ndarray[0][0])
        for incidentEnergyIndex, incidentEnergyData in enumerate(ndarray):
            start = None
            length = 0
            for outgoingEnergyIndex, outgoingEnergyData in enumerate(incidentEnergyData):
                if min(outgoingEnergyData) != 0.0 or max(outgoingEnergyData) != 0.0:
                    if start is None: start = n3 * (n2 * incidentEnergyIndex + outgoingEnergyIndex)
                    length += n3
                    flattenedData += list(outgoingEnergyData)
                else:
                    if start is not None:
                        starts.append(start)
                        lengths.append(length)
                        start = None
                        length = 0
            if start is not None:
                starts.append(start)
                lengths.append(length)

    starts = valuesModule.Values(starts, valueType=xDataEnumsModule.ValueType.integer32)
    lengths = valuesModule.Values(lengths, valueType=xDataEnumsModule.ValueType.integer32)
    array = arrayModule.Flattened(shape=ndarray.shape, starts=starts, lengths=lengths, data=valuesModule.Values(flattenedData))
    axes = axesTemplate3d.copy()
    gridded3d = multiGroupModule.Gridded3d(axes=axes, array=array)

    return multiGroupModule.Form(label, xDataEnumsModule.Frame.lab, gridded3d)

def reactionMultiGrouping(protare, reactionOut, MG, MG_upscatter, temperatureInfo, axesTemplate1d, axesTemplate3d):
    """
    For the specified temperature, sums the multi-group data.

    :param protare:             The protare to process.
    :param reactionOut:         The reaction which the summed, multi-group data will be added to.
    :param MG:                  The multiGroupSettings.MG instance for non-upscatter data.
    :param MG_upscatter:        The multiGroupSettings.MG instance for upscatter data.
    :param temperatureInfo:     TemperatureInfo instance whose HeatedMultiGroup or SnElasticUpScatter label specifies the multi-group data to retrieve.
    :param axesTemplate1d:      The axes template for 1d multi-group data.
    :param axesTemplate3d:      The axes template for distribution multi-group data.
    """

    label = temperatureInfo.heatedMultiGroup

    data = protare.multiGroupCrossSection(MG, temperatureInfo)
    reactionOut.crossSection.add(toGridded1d(crossSectionModule, axesTemplate1d, data, label))

    data = protare.multiGroupAvailableEnergy(MG, temperatureInfo)
    reactionOut.availableEnergy.add(toGridded1d(availableEnergyModule, axesTemplate1d, data, label))

    data = protare.multiGroupAvailableMomentum(MG, temperatureInfo)
    reactionOut.availableMomentum.add(toGridded1d(availableMomentumModule, axesTemplate1d, data, label))

    data = protare.multiGroupQ(MG, temperatureInfo, True)
    reactionOut.outputChannel.Q.add(toGridded1d(QModule, axesTemplate1d, data, label))

    for transportableID in transportableIDs:
        try:
            data = protare.multiGroupMultiplicity(MG, temperatureInfo, transportableID)
        except:
            continue
        if len(data) > 0:
            try:
                product = reactionOut.outputChannel.products[transportableID]
            except:
                product = productModule.Product(transportableID, transportableID)
                reactionOut.outputChannel.products.add(product)

            product.multiplicity.add(toGridded1d(multiplicityModule, axesTemplate1d, data, label))

            data = protare.multiGroupAverageEnergy(MG, temperatureInfo, transportableID)
            product.averageProductEnergy.add(toGridded1d(averageProductEnergyModule, axesTemplate1d, data, label))

            data = protare.multiGroupAverageMomentum(MG, temperatureInfo, transportableID)
            product.averageProductMomentum.add(toGridded1d(averageProductMomentumModule, axesTemplate1d, data, label))

            data = protare.multiGroupProductArray(MG, temperatureInfo, None, transportableID)
            product.distribution.add(toGridded3d(axesTemplate3d, data, label))

            upscatterLabel = temperatureInfo.SnElasticUpScatter
            if transportableID == IDsModule.neutron and upscatterLabel != '':
                data = protare.multiGroupAverageEnergy(MG_upscatter, temperatureInfo, transportableID)
                product.averageProductEnergy.add(toGridded1d(averageProductEnergyModule, axesTemplate1d, data, upscatterLabel))

                data = protare.multiGroupProductArray(MG_upscatter, temperatureInfo, None, transportableID)
                product.distribution.add(toGridded3d(axesTemplate3d, data, upscatterLabel))

def addMultiGroupSums(self, replace=False):
    """
    Calculates total and delayedNeutron multi-group sums for *self* and adds the results to the applicationData node as
    separate institutions with labels "LLNL::multiGroupReactions" and "LLNL::multiGroupDelayedNeutrons", respectivley.

    :param replace:           If summed multi-group data already present, then a raise is executed if *replace* is **False**, otherwise the old data are replaced.
    """

    if tokensModule.multiGroupReactions in self.applicationData:
        if not replace:
            raise Exception('Summed "%s" data already present.' % tokensModule.multiGroupReactions)
        self.applicationData.remove(tokensModule.multiGroupReactions)

    if tokensModule.multiGroupDelayedNeutrons in self.applicationData:
        if not replace:
            raise Exception('Summed "%s" data already present.' % tokensModule.multiGroupDelayedNeutrons)
        self.applicationData.remove(tokensModule.multiGroupDelayedNeutrons)

    MG = multiGroupSettingsModule.MG(self.projectile, transportingModule.Mode.multiGroup, delayedNeutrons=transportingModule.DelayedNeutrons.off)
    MG_upscatter = multiGroupSettingsModule.MG(self.projectile, transportingModule.Mode.multiGroupWithSnElasticUpScatter, delayedNeutrons=transportingModule.DelayedNeutrons.off)
    MG_delayedNeutrons = multiGroupSettingsModule.MG(self.projectile, transportingModule.Mode.multiGroup, delayedNeutrons=transportingModule.DelayedNeutrons.on)

    totalReaction = reactionModule.Reaction('total', enumsModule.Genre.NBody, 1)

    fissionFragmentData = None
    fissionReaction = self.getReaction('fission')
    if fissionReaction is not None:
        fissionFragmentData = fissionReaction.outputChannel.fissionFragmentData
        outputChannelDelayedNeutrons = outputChannelModule.OutputChannel(enumsModule.Genre.NBody)
        productDelayedNeutrons = productModule.Product(IDsModule.neutron, IDsModule.neutron)
        outputChannelDelayedNeutrons.products.add(productDelayedNeutrons)

    elasticReaction = self.reactions[0]

    for temperatureInfo in self.styles.temperatures():

        label = temperatureInfo.heatedMultiGroup
        if label == '':
            continue

        axesTemplate1d = elasticReaction.crossSection[label].axes
        axesTemplate3d = elasticReaction.outputChannel.products[0].distribution[label].multiGroupSubform.axes

        reactionMultiGrouping(self, totalReaction, MG, MG_upscatter, temperatureInfo, axesTemplate1d, axesTemplate3d)

        if fissionFragmentData is not None:
            try:
                data = fissionFragmentData.multiGroupMultiplicity(MG_delayedNeutrons, temperatureInfo, IDsModule.neutron)
            except:
                print('WARNING: Delayed neutron data not complete and will not be included.')
                fissionFragmentData = None                      # Some delayed neutrons have unspecified multiplicities.
                data = []

            if len(data) > 0:
                numberOfProjectileGroups = len(data)
                productDelayedNeutrons.multiplicity.add(toGridded1d(multiplicityModule, axesTemplate1d, data, label))

                data = fissionFragmentData.multiGroupQ(MG_delayedNeutrons, temperatureInfo)
                if len(data) == 0:
                    data = numberOfProjectileGroups * [0]
                outputChannelDelayedNeutrons.Q.add(toGridded1d(QModule, axesTemplate1d, data, label))

                data = fissionFragmentData.multiGroupAverageEnergy(MG_delayedNeutrons, temperatureInfo, IDsModule.neutron)
                productDelayedNeutrons.averageProductEnergy.add(toGridded1d(averageProductEnergyModule, axesTemplate1d, data, label))

                data = fissionFragmentData.multiGroupAverageMomentum(MG_delayedNeutrons, temperatureInfo, IDsModule.neutron)
                productDelayedNeutrons.averageProductMomentum.add(toGridded1d(averageProductMomentumModule, axesTemplate1d, data, label))

                data = fissionFragmentData.multiGroupProductArray(MG_delayedNeutrons, temperatureInfo, None, IDsModule.neutron)
                productDelayedNeutrons.distribution.add(toGridded3d(axesTemplate3d, data, label))

    if len(totalReaction.crossSection) > 0:
        institution = institutionModule.Institution(tokensModule.multiGroupReactions)
        institution.append(totalReaction)
        self.applicationData.add(institution)

    if fissionFragmentData is not None:
        if len(outputChannelDelayedNeutrons.Q) > 0:
            institution = institutionModule.Institution(tokensModule.multiGroupDelayedNeutrons)
            institution.append(outputChannelDelayedNeutrons)
            self.applicationData.add(institution)
