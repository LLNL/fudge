#!/usr/bin/env python

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse
import sys
import os
import matplotlib.pyplot as plt

sys.path.append(os.path.split(__file__)[0] + os.sep + '..')

from brownies.BNL.utilities.io import *


# -------------------------------------------------------------------------------
# Comamnd line parsing
# -------------------------------------------------------------------------------
def parse_arguments():
    """Parse the command line"""
    parser = argparse.ArgumentParser(description="Extract the spectrum of an outgoing particle at a given energy")

    # Required command line options
    parser.add_argument('ENDF', type=str, help='the ENDF file(s) whose cross section you want to study.')

    # Set output
    parser.add_argument('-p', dest='plot', action='store_true', default=False, help="plot the spectrum")
    parser.add_argument('-t', dest='table', action='store_true', default=False, help="print the spectrum as a table")
    parser.add_argument('-o', dest='outPrefix', default=None, type=str, help='output to a file starting with OUTFILE, '
                                                                             'instead of printing to stdout or making '
                                                                             'an interactive plot')

    # Verbosity
    parser.add_argument('-v', dest='verbose', default=False, action='store_true', help="enable verbose output.")

    # Control what spectrum to extract
    parser.add_argument('-E', dest='energy', default=0.0235339347844, type=float,
                        help="Incident energy in eV (Default: 0.0235339347844 eV)")   # v=2200 m/s
    parser.add_argument('--MT', type=int, default=102, help='if given, work with this MT (Default: 102).')
    parser.add_argument('--product', default='photon', type=str,
                        choices=['n', 'p', 'H1', 'd', 'H2', 't', 'H3', 'He3', 'a', 'He4', 'g', 'photon'],
                        help='product whose spectrum we want to get (Default: "g")')
    parser.add_argument('-y', dest='yields', action='store_true', default=False, help="compute yield, not probability")

    return parser.parse_args()


# -------------------------------------------------------------------------------
# Main!
# -------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_arguments()

    # FIXME: only tests for particles named "photon" and "n", may need translators for others

    # Read in evaluation
    print(10*'*', 'Reading %s' % args.ENDF, 10*'*')
    theEvaluation = read_evaluation(args.ENDF)['reactionSuite']

    # Make sure the requested particle actually is in the reactionSuite,
    # otherwise there's no point in digging for its distributions.
    if theEvaluation.hasParticle(args.product):
        if args.verbose:
            print("Found particle %s in the evaluation" % args.product)
    else:
        raise ValueError("Cannot find particle %s in %s's reactionSuite" % (args.product, args.ENDF))

    # Get the correct reaction
    theReaction = theEvaluation.getReaction(args.MT)
    if theReaction is None:
        raise ValueError("Cannot find reaction MT=%i in %s's reactionSuite" % (args.MT, args.ENDF))
    else:
        if args.verbose:
            print("Found MT=%i in the evaluation" % args.MT)

    # Check that the energy in question is above threshold
    Ethreshold = theReaction.getThreshold('eV')
    if args.energy < Ethreshold:
        raise ValueError("Requested energy %s eV < Ethreshold = %s eV" % (args.energy, Ethreshold))
    if args.verbose:
        print("Reaction threshold for MT=%i is %s eV" % (args.MT, Ethreshold))

    # Get the product from the reaction outgoingChannel
    theProducts = theReaction.outputChannel.getProductsWithName(args.product)
    if theProducts is None or len(theProducts) == 0:
        raise ValueError("Cannot find product %s in reaction MT=%i" % (args.product, args.MT))
    else:
        if args.verbose:
            print("Found product %s in reaction MT=%i" % (args.product, args.MT))

    # FIXME: Add check for product the reaction **should** have, but doesn't

    # FIXME: Add logic for product that is there, but you have to get the info from PoPs (say a photon from (n,inel))

    print()
    print(10*'*', 'Processing product %s from MT=%i in file %s' % (args.product, args.MT, args.ENDF), 10*'*')
    # Loop over the product(s), do a yield weighted sum of the spectrum
    continuousSpectrum = None
    discreteLines = []
    totalYield = 0.0
    for iprod, prod in enumerate(theProducts):
        Y = prod.multiplicity.evaluate(args.energy)
        dist = prod.distribution['eval']
        print("For %s #%i at energy %s eV, got yield=%s" % (args.product, iprod, args.energy, Y))
        print("    Product frame:", dist.productFrame)
        print("    Distribution moniker:", dist.moniker)
        if dist.moniker != 'uncorrelated':
            raise NotImplementedError("FIXME: implement correlated distributions, found type %s" % dist.moniker)
        print("    Energy sub-distribution moniker:", dist.energySubform.data.moniker)
        if hasattr(dist.energySubform.data, 'getSpectrumAtEnergy'):
            enDist = dist.getSpectrumAtEnergy(args.energy)
            if continuousSpectrum is None:
                continuousSpectrum = enDist * Y
            else:
                continuousSpectrum += enDist * Y
            totalYield += Y
        elif dist.energySubform.data.moniker in ['discreteGamma']:
            discreteLines.append([dist.energySubform.data.value, Y])
        elif dist.energySubform.data.moniker in ['primaryGamma']:
            discreteLines.append([dist.energySubform.data.averageEp(args.energy), Y])
        else:
            msg = "FIXME: implement 'getSpectrumAtEnergy' for energy sub-distribution %s" % \
                  dist.energySubform.data.moniker
            if True:
                raise NotImplementedError(msg)
            else:
                print(msg)

    # Take care of normalization
    if continuousSpectrum is not None and not args.yields:
        continuousSpectrum = continuousSpectrum / totalYield

    # Plotting output
    if args.plot:
        figs, ax = plt.subplots()
        plt.ticklabel_format(axis='both', style='sci', scilimits=(-3, 3))
        plt.title("Spectrum for %s at %s eV from MT=%i from %s" % (args.product, args.energy, args.MT, args.ENDF))

        if continuousSpectrum is not None:
            xdata, ydata = continuousSpectrum.copyDataToXsAndYs()
            ax.plot(xdata, ydata, linewidth=2, label='continuum')
        if discreteLines:
            xdata = [x[0] for x in discreteLines]
            ydata = [x[1] for x in discreteLines]
            ax.stem(xdata, ydata, label='discrete', use_line_collection=True)

        ax.legend()
        ax.set_xlabel('$E_{%s}$ [eV]' % args.product)
        if args.yields:
            ax.set_ylabel('Y [# %s/eV]' % args.product)
        else:
            ax.set_ylabel('P [1/eV]')

        if args.outPrefix:
            plt.savefig(args.outPrefix+'.png')
        else:
            plt.show()

    # Table output
    if args.table:
        output = ''
        if continuousSpectrum is not None:
            output += "# Continuous spectrum\n"
            output += "# ===================\n"
            for pair in continuousSpectrum:
                output += ' '.join([str(x) for x in pair])+'\n'
        if discreteLines:
            output += '\n'
            output += "# Discrete lines\n"
            output += "# ==============\n"
            for pair in discreteLines:
                output += ' '.join([str(x) for x in pair]) + '\n'

        if args.outPrefix is None:
            print(output)
        else:
            with open(args.outPrefix+'.txt', mode='w') as outfile:
                outfile.write(output)
