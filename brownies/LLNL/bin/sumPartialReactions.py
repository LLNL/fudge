#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

description = """
GNDS supports more fine-grained reactions than ENDF-6, such as capture to specific excited states
and reactions like (n,h+2a) that have no MT equivalent.

This utility takes a reactionSuite and sums together these
fine-grained reactions for compatibility with ENDF-6
"""
import argparse

from xData import enums as xDataEnumsModule

from fudge import enums as enumsModule
from fudge import reactionSuite as reactionSuiteModule
from fudge import product as productModule
from fudge import outputChannel as outputChannelModule
from fudge.reactions import reaction as reactionModule
from fudge.reactionData import crossSection as crossSectionModule
from fudge.productData import multiplicity as multiplicityModule
from fudge.productData import distributions as distributionsModule

from PoPs import IDs as PoPsIDsModule

parser = argparse.ArgumentParser(description)
parser.add_argument("input", help="reactionSuite to read and sum")
parser.add_argument("output", help="file name where summed reactionSuite will be saved")
parser.add_argument("-f", "--sumFission", action="store_true",
        help="Sum all partial fission reactions together into a single MT=18 reaction")
parser.add_argument("-c", "--sumCapture", action="store_true",
        help="Sum all partial capture reactions together into a single MT=102 reaction")
parser.add_argument("-5", "--sumMT5", action="store_true",
        help="Sum all reactions with no MT equivalent together into a single MT=5 reaction")


#<editor-fold desc="helper methods" defaultstate="collapsed">
def sumCrossSection( reactionList, info ):
    """
    Compute summed cross section for list of reactions,
    @return lin-lin crossSection.XYs1d instance
    """
    sum_ = None
    for reac in reactionList:
        xsc = reac.crossSection.toPointwise_withLinearXYs()
        if sum_ is None:
            sum_ = xsc
        else:
            sum_, xsc = sum_.mutualify(1e-8,1e-8,0, xsc, 1e-8,1e-8,0)
            sum_ += xsc

    # wrap result in correct class
    result = crossSectionModule.XYs1d( label=info['style'], axes=info['crossSectionAxes'],
            data=sum_.copyDataToXYs() )
    return result


def replacePartialsWithSum( ReactionSuite, partialReactions, newSummedReaction ):
    """
    Replace original list of reactions with new sum
    """
    # FIXME a bit awkward, have to remove all reactions, replace desired one and put them all back...
    fullReactionList = [reaction for reaction in ReactionSuite.reactions]

    index = fullReactionList.index(partialReactions[0])
    for reaction in partialReactions:
        fullReactionList.remove( reaction )
    fullReactionList.insert(index, newSummedReaction)

    ReactionSuite.reactions.clear()
    for reaction in fullReactionList:
        ReactionSuite.reactions.add( reaction )


def removeSums( ReactionSuite, reactionList, newTotal ):
    """
    Combining discrete reactions may mess up a summed cross section or multiplicity
    since it may sum over terms that no longer exist.
    Simplest way to deal with that: remove any sums that include removed terms UNLESS
    all of the removed terms are present and can be replaced with newTotal
    """
    from fudge import sums as sumsModule

    reactionListIds = list(map(id, reactionList))
    sumsToRemove = []
    for crossSectionSum in ReactionSuite.sums.CrossSectionSums:
        reactionsInSum = [summand.link.findClassInAncestry(reactionModule.Reaction)
                for summand in crossSectionSum.summands]
        reactionsInSumIds = list(map(id, reactionsInSum))

        intersection = set(reactionListIds).intersection( reactionsInSumIds )
        if not intersection:
            continue

        if intersection == set(reactionListIds):
            # can replace list of partial summands with newTotal
            for summand in crossSectionSum.summands.summandList[::-1]:
                if id(summand.link) in intersection:
                    crossSectionSum.summands.summandList.remove(summand)

            if len(crossSectionSum.summands) == 0:
                # this sum is no longer necessary, now present in ReactionSuite.reactions
                sumsToRemove.append( crossSectionSum )
            else:
                newSummand = sumsModule.Add( newTotal )
                crossSectionSum.summands.summandList.append( newSummand )

        elif intersection < set(reactionListIds):
            sumsToRemove.append( crossSectionSum )

        else:
            raise NotImplementedError("Shouldn't get here")

    for crossSectionSum in sumsToRemove:
        ReactionSuite.sums.crossSectionSums.remove( crossSectionSum.label )


def gammaCascade( PoPs, parentId ):
    """
    This method returns a complete gamma cascade starting with the parent particle
    down to the lowest possible state.

    :return: dictionary of {discreteEnergy: total_probability} for all gammas in the cascade
    """
    # FIXME method should be defined in PoPs, perhaps in PoPs.decays.misc.
    # Having svn trouble right now so defining it here for time being
    from PoPs import IDs as IDsModule
    parentNuclide = PoPs[parentId]

    def recursive_helper( product, probability, gammas ) :

        initialEnergy = product.nucleus.energy[0].value
        for decayMode in product.decayData.decayModes :
            branchingRatio = decayMode.probability[0].value * probability

            decayPath = decayMode.decayPath[0]
            residual = decayPath.products[0].pid
            if( residual == IDsModule.photon ) : residual = decayPath.products[1].pid
            residual = PoPs[residual]

            gammaEnergy = initialEnergy - residual.nucleus.energy[0].value
            if( gammaEnergy not in gammas ) : gammas[gammaEnergy] = 0.0
            gammas[gammaEnergy] += branchingRatio

            recursive_helper( residual, branchingRatio, gammas )

    gammas = {}
    recursive_helper( parentNuclide, 1.0, gammas )
    return gammas
#</editor-fold>


#<editor-fold desc="sum reactions" defaultstate="collapsed">
def sumCapture( ReactionSuite, info ):
    """
    Combine all capture reactions in ReactionSuite

    @param ReactionSuite: reactionSuite instance, will be modified in-place
    @param info: dictionary
    """

    photonId = PoPsIDsModule.photon
    captureReacs = [r for r in ReactionSuite.reactions if r.ENDF_MT == 102]

    if len(captureReacs) <= 1:  # nothing to do
        return ReactionSuite

    crossSection = sumCrossSection( captureReacs, info )
    Q = max([(reac.outputChannel.Q.getConstantAs(info['energyUnit']), reac.outputChannel.Q.evaluated)
        for reac in captureReacs])

    newCapture = reactionModule.Reaction(enumsModule.Genre.NBody, ENDF_MT=102)
    newCapture.crossSection.add( crossSection )
    newCapture.outputChannel.Q.add( Q[1] )

    # now combine all outgoing photons, recomputing multiplicities.
    # also need to add discrete photons if residual is left in a discrete excited state
    primaryGammas = []
    discreteGammas = []
    continuumGamma = None
    residualProduct = None

    for reaction in captureReacs:
        residual, = [product for product in reaction.outputChannel if product.pid != photonId]

        xsc = reaction.crossSection.toPointwise_withLinearXYs()
        multiplicityFactor = xsc / crossSection

        if reaction.outputChannel.genre == enumsModule.Genre.twoBody:
            # primary gamma, often followed by cascade of discrete gammas

            photon = reaction.outputChannel.getProductWithName(photonId)

            oldMult = photon.multiplicity.pop(info['style'])
            assert oldMult.rangeMax == 1

            newMult = multiplicityModule.XYs1d(axes=oldMult.axes, data=multiplicityFactor.copyDataToXYs())
            newMult.label = info['style']
            photon.multiplicity.add(newMult)

            # change distribution from 2-body to uncorrelated primaryGamma
            primaryGammaValue = reaction.outputChannel.Q.getConstantAs( info['energyUnit'] )
            angularDist = distributionsModule.angular.Isotropic2d()
            energyDist = distributionsModule.energy.PrimaryGamma( primaryGammaValue,
                    xsc.domainMin, xsc.domainMax,
                    axes = info['energySpectrumAxes'] )
            newDistribution = distributionsModule.uncorrelated.Form(
                    info['style'], xDataEnumsModule.Frame.lab,
                    distributionsModule.uncorrelated.AngularSubform( angularDist ),
                    distributionsModule.uncorrelated.EnergySubform( energyDist ) )
            photon.distribution.clear()
            photon.distribution.add( newDistribution )

            primaryGammas.append( photon )

            if residual.outputChannel:  # also have cascade of discrete gammas
                photon = residual.outputChannel.getProductWithName(photonId)
                distribution = photon.distribution.evaluated
                if isinstance(distribution, distributionsModule.branching3d.Form):

                    # FIXME: the same discrete gamma currently gets split into multiple products,
                    # for example it may appear both in the gamma cascade from 1st and 2nd excited states.
                    # Those should be combined together into a single product.

                    discretes = gammaCascade( ReactionSuite.PoPs,
                            photon.parentProduct.pid )

                    for energy, probability in discretes.items():
                        newDiscrete = productModule.Product(photonId)

                        newMultiplicity = multiplicityModule.XYs1d(
                                label=info['style'],
                                axes=info['multiplicityAxes'],
                                data=(probability * multiplicityFactor).copyDataToXYs()
                        )
                        newDiscrete.multiplicity.add(newMultiplicity)

                        angularDist = distributionsModule.angular.Isotropic2d()
                        energyDist = distributionsModule.energy.DiscreteGamma( energy,
                                xsc.domainMin, xsc.domainMax,
                                axes = info['energySpectrumAxes'] )
                        newDistribution = distributionsModule.uncorrelated.Form(
                                info['style'], xDataEnumsModule.Frame.lab,
                                distributionsModule.uncorrelated.AngularSubform( angularDist ),
                                distributionsModule.uncorrelated.EnergySubform( energyDist ) )
                        newDiscrete.distribution.add( newDistribution )

                        discreteGammas.append( newDiscrete )
                else:
                    print("Not yet implemented: %s distribution for residual photon emission"
                            % type(distribution))


        else:   # continuum gamma. Photon may be in reaction outputChannel or in 1st product outputChannel
            photons = reaction.outputChannel.getProductsWithName(photonId)
            if not photons:
                assert len(reaction.outputChannel) == 1
                photons = reaction.outputChannel[0].outputChannel.getProductsWithName(photonId)

            assert len(photons) == 1
            photon = photons[0]
            oldMult = photon.multiplicity.pop(info['style'])
            factor = multiplicityFactor.domainSlice( *oldMult.domain() )
            newMult = multiplicityModule.XYs1d(axes=oldMult.axes, data=oldMult * factor)
            newMult.label = info['style']
            photon.multiplicity.add(newMult)

            continuumGamma = photon

            if residual.outputChannel:
                residual, = [product for product in residual.outputChannel if product.pid != photonId]
            residualProduct = residual

    # add residual + all photons to outputChannel:
    for product in [ residualProduct, continuumGamma ] + primaryGammas + discreteGammas:
        if product is not None:
            newCapture.outputChannel.products.uniqueLabel( product )
            newCapture.outputChannel.products.add( product )

    replacePartialsWithSum( ReactionSuite, captureReacs, newCapture )

    removeSums( ReactionSuite, captureReacs, newCapture )


def sumFission( ReactionSuite, info ):
    """
    Combine all fission reactions in ReactionSuite

    @param ReactionSuite: reactionSuite instance, will be modified in-place
    @param info: dictionary
    """

    photonId = PoPsIDsModule.photon
    neutronId = PoPsIDsModule.neutron
    fissionReacs = [r for r in ReactionSuite.reactions if r.ENDF_MT == 18]

    if len(fissionReacs) <= 1:  # nothing to do
        return ReactionSuite

    crossSection = sumCrossSection( fissionReacs, info )
    Q = max([(reac.outputChannel.Q.getConstantAs(info['energyUnit']), reac.outputChannel.Q.evaluated)
        for reac in fissionReacs])

    newFission = reactionModule.Reaction(enumsModule.Genre.NBody, ENDF_MT = 18, fissionGenre=enumsModule.FissionGenre.total)
    newFission.crossSection.add( crossSection )
    newFission.outputChannel.Q.add( Q[1] )

    # add dummy neutron and gamma products. Hauser-Feshbach codes don't model fission output, so
    # need to obtain that information from somewhere else

    neutron = productModule.Product(neutronId)

    # place-holder for nubar
    nubar = multiplicityModule.XYs1d(
            label=info['style'],
            axes=info['multiplicityAxes'],
            data=[[crossSection.domainMin, 2.5], [crossSection.domainMax, 5.2]]
    )
    neutron.multiplicity.add(nubar)

    angularDist = distributionsModule.angular.Isotropic2d()
    energyDist = distributionsModule.energy.XYs2d(axes=info['energySpectrumAxes'])
    for energy in (crossSection.domainMin, crossSection.domainMax):
        energyDist.append( distributionsModule.energy.XYs1d(outerDomainValue=energy,
            data=[[1e-5, 2e+5], [2e-5, 0]] ) )

    newDistribution = distributionsModule.uncorrelated.Form(
            info['style'], xDataEnumsModule.Frame.lab,
            distributionsModule.uncorrelated.AngularSubform( angularDist ),
            distributionsModule.uncorrelated.EnergySubform( energyDist ) )
    neutron.distribution.add( newDistribution )

    newFission.outputChannel.products.uniqueLabel( neutron )
    newFission.outputChannel.products.add( neutron )

    replacePartialsWithSum( ReactionSuite, fissionReacs, newFission )

    removeSums( ReactionSuite, fissionReacs, newFission )


def sumMT5( ReactionSuite, info ):
    """
    Combine all reactions with no MT equivalent into MT=5

    @param ReactionSuite: reactionSuite instance, will be modified in-place
    @param info: dictionary
    """

    mt5Reacs = [r for r in ReactionSuite.reactions if r.ENDF_MT == -1]
    raise NotImplementedError("Haven't supported this step yet")
#</editor-fold>


if __name__ == '__main__':
    args = parser.parse_args()
    RS = reactionSuiteModule.ReactionSuite.readXML_file( args.input )
    info = {
            'style': RS.styles.getEvaluatedStyle().label,
            'energyUnit': RS.domainUnit,
            'crossSectionAxes': crossSectionModule.defaultAxes( RS.domainUnit ),
            'multiplicityAxes': multiplicityModule.defaultAxes( RS.domainUnit ),
            'energySpectrumAxes': distributionsModule.energy.defaultAxes( RS.domainUnit )
            }

    if args.sumCapture:
        sumCapture( RS, info )
    if args.sumFission:
        sumFission( RS, info )
    if args.sumMT5:
        sumMT5( RS, info )

    RS.saveToFile( args.output )

