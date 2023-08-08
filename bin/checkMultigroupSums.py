#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import pathlib
import argparse
import numpy

from fudge import reactionSuite as reactionSuiteModule
from fudge import styles as stylesModule
from fudge.processing.deterministic import tokens as tokensModule

description = '''
This script checks the multi-group sums data in the list of specified protares and prints the xlink of the data where issues are found.
If the option "--outputDir" is specfied and a protare has issues, then the protare with its issues fixed is written into the directory 
specified by this option and with the same file name as the original file.
'''

rtolDefault = 1e-8

parser = argparse.ArgumentParser(description=description)

parser.add_argument('files', nargs='+', type=pathlib.Path,      help = 'GNDS protares whose multi-group sums are checked.')
parser.add_argument('--outputDir', default=None,                help = 'Name of output directory to write fixed protares to.')
parser.add_argument('--rtol', type=float, default=rtolDefault,  help = 'Relative tolerance for comparing data. Default = %s.' % rtolDefault)

args = parser.parse_args()

labels = []
upscatterLabels = []
numberOfOIssues = 0

if args.outputDir is not None:
    args.outputDir = pathlib.Path(args.outputDir)

def message(msg):

    global numberOfOIssues

    numberOfOIssues += 1
    print(msg)

def whichLabels(suite1, suite2):
    '''
    Adds the contents of *upscatterLabels* to the list of labels to check for if required.
    '''

    allLabels = labels
    for label in upscatterLabels:
        if label in suite1:
            allLabels = labels + upscatterLabels
            break
        if label in suite2:
            allLabels = labels + upscatterLabels
            break

    return allLabels

def checkGridded1d(suite1, suite2):
    '''
    This functions compares the gridded1d multi-group data in *suite1* to that in *suite2* which are specified by *labels* 
    and maybe *upscatterLabels*.
    '''

    for label in whichLabels(suite1, suite2):
        if label in suite1:
            data1 = suite1[label]
        else:
            data1 = None
            message('    WARNING: %s data for label "%s" not in original.' % (suite1.moniker, label))

        if label in suite2:
            data2 = suite2[label]
        else:
            data2 = None
            message('    WARNING: %s data for label "%s" not in new.' % (suite2.moniker, label))

        if data1 is not None and data2 is not None:
            vector1 = data1.constructVector().vector
            vector2 = data2.constructVector().vector
            if not numpy.allclose(vector1, vector2, rtol=args.rtol, atol=0.0):
                message('    WARNING: %s in %s.' % (label, suite2.toXLink()))

def checkGridded3d(suite1, suite2):
    '''
    This function compares the griddede3d (i.e., distribution) multi-group data in *suite1* to that in *suite2* which are 
    specified by *labels* and maybe *upscatterLabels*.
    '''

    for label in whichLabels(suite1, suite2):
        array1 = suite1[label].multiGroupSubform.array.constructArray()
        array2 = suite2[label].multiGroupSubform.array.constructArray()
        if(array1.shape != array2.shape):
            message('    WARNING: %s in %s.' % (label, suite2.toXLink()))
        else:
            if not numpy.allclose(array1, array2, rtol=args.rtol, atol=0.0):
                message('    WARNGING: %s in %s.' % (label, suite2.toXLink()))

def checkProduct(productOrig, productNew):
    '''
    This function checks the multiplicity, distribution, averageProductEnergy and averageProductMomentum for the orginal and new product.
    '''

    if productOrig is None:
        message('    WARNING: product "%s" missing from original.' % productNew.pid)
    elif productNew is None:
        message('    WARNING: product "%s" missing from new.' % productOrig.pid)
    else:
        checkGridded1d(productOrig.multiplicity, productNew.multiplicity)
        checkGridded3d(productOrig.distribution, productNew.distribution)
        checkGridded1d(productOrig.averageProductEnergy, productNew.averageProductEnergy)
        checkGridded1d(productOrig.averageProductMomentum, productNew.averageProductMomentum)

def checkProducts(productsOrig, productsNew):
    '''
    This function checks each product in the list of products in *productsOrig* and *productsNew*.
    '''

    pids = set([product.pid for product in productsOrig] + [product.pid for product in productsNew])
    for pid in pids:
        if pid in productsOrig:
            product1 = productsOrig[pid]
        else:
            product1 = None

        if pid in productsNew:
            product2 = productsNew[pid]
        else:
            product2 = None

        checkProduct(product1, product2)

def checkOutputChannel(outputChannelOrig, outputChannelNew):
    '''
    This function checks the original and new outputChannels.
    '''

    checkGridded1d(outputChannelOrig.Q, outputChannelNew.Q)
    checkProducts(outputChannelOrig.products, outputChannelNew.products)

def checkReaction(totalReactionOrig, totalReactionNew):
    '''
    This function checks the original and new reactions.
    '''

    checkGridded1d(totalReactionOrig.crossSection, totalReactionNew.crossSection)
    checkOutputChannel(totalReactionOrig.outputChannel, totalReactionNew.outputChannel)
    checkGridded1d(totalReactionOrig.availableEnergy, totalReactionNew.availableEnergy)
    checkGridded1d(totalReactionOrig.availableMomentum, totalReactionNew.availableMomentum)

for file in args.files:
    numberOfOIssues = 0
    print(file)
    protare = reactionSuiteModule.read(file, lazyParsing=True)

    labels = [style.label for style in protare.styles if isinstance(style, stylesModule.HeatedMultiGroup)]
    upscatterLabels = [style.label for style in protare.styles if isinstance(style, stylesModule.SnElasticUpScatter)]
#
# Get current data in file.
#
    try:
        totalReactionOrig = protare.applicationData.pop(tokensModule.multiGroupReactions).data[0]
    except:
        totalReactionOrig = None

    try:
        delayedNeutronsOrig = protare.applicationData.pop(tokensModule.multiGroupDelayedNeutrons).data[0]
    except:
        delayedNeutronsOrig = None

    try:
        productsOrig = protare.applicationData.pop(tokensModule.multiGroupIncompleteProducts).data[0]
    except:
        productsOrig = None
#
# Recalculate multi-group summed data.
#
    protare.addMultiGroupSums()
#
# Check multi-group summed total reaction data.
#
    try:
        totalReactionNew = protare.applicationData[tokensModule.multiGroupReactions].data[0]
    except:
        totalReactionNew = None

    if totalReactionOrig is None:
        if totalReactionNew is not None:
            message('    WARNING: original file missing multi-group summed data.')
    else:
        if totalReactionNew is None:
            message('    WARNING: new file missing multi-group summed data.')
        else:
            checkReaction(totalReactionOrig, totalReactionNew)
#
# Check multi-group summed delayed neutron data.
#
    try:
        delayedNeutronsNew = protare.applicationData[tokensModule.multiGroupDelayedNeutrons].data[0]
    except:
        delayedNeutronsNew = None

    if delayedNeutronsOrig is None:
        if delayedNeutronsNew is not None:
            message('    WARNING: original file missing multi-group summed data.')
    else:
        if delayedNeutronsNew is None:
            message('    WARNING: new file missing multi-group summed data.')
        else:
            checkOutputChannel(delayedNeutronsOrig, delayedNeutronsNew)
#
# Check multi-group summed incomplete products data.
#
    try:
        productsNew = protare.applicationData[tokensModule.multiGroupIncompleteProducts].data[0]
    except:
        productsNew = None

    if productsOrig is None:
        if productsNew is not None:
            message('    WARNING: original file missing multi-group summed data for incomplete products.')
    else:
        if productsNew is None:
            message('    WARNING: new file missing multi-group summed data for incomplete products.')
        else:
            checkProducts(productsOrig, productsNew)

    if numberOfOIssues > 0 and args.outputDir is not None:
        protare.saveAllToFile(args.outputDir / file.name, hybrid=True)
