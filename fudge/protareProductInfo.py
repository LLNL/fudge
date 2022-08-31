#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse

from PoPs import IDs as PoPsIDsModule
from PoPs.families import particle as PoPsParticleModule
from fudge import reactionSuite as reactionSuiteModule

description = '''For a protare (i.e., reactionSuite), list each reaction and its outgoing particles (i.e., products).'''

__doc__ = description

def protareProductInfo(fileName, indent):
    '''
    Opens the protare file *fileName* and calls **protareProductInfo2**.

    :param fileName:    The file name of the protare to read in.
    :param indent:      The amount of indentation for each line.

    :return:            The returned results from **protareProductInfo2**.
    '''

    return protareProductInfo2(reactionSuiteModule.read(fileName, lazyParsing=True), indent)

def protareProductInfo2(protare, indent):
    '''
    Constructs a line for a **ris** file for each reaction in *protare*. The list of lines constructed is returned.
    Each line contains four items separated by the colon (i.e., ":") character. The four items are:

        -) A formula for all final products.
        -) The effective threshold for the reaction.
        -) An intermediate product if one exists.
        -) The reaction process.

    :param protare:     A ReactionSuite instance whose reactions are analyzed.
    :param indent:      The amount of indentation for each line.

    :return:            The list of lines constructed.
    '''

    lines = []
    for reaction in protare.reactions:
        effectiveThreshold = reaction.crossSection.effectiveThreshold()

        reactionProducts = reaction.reactionProducts()
        particles2, intermediates, processes, dummy = reactionProducts.asSortedList(protare.PoPs, excludePhotons=False)
        particles = {}
        for pid, multiplicity in particles2:
            particles[pid] = multiplicity

        sortedParticles = PoPsParticleModule.ParticleSorter()
        for pid in particles:
            sortedParticles.add(protare.PoPs[pid])

        photons = ''
        label = ''
        sep = ''
        for particle in sortedParticles:
            multiplicity = ''
            if particles[particle.id] != 1:
                multiplicity = '%s' % particles[particle.id]
            if particle.id == PoPsIDsModule.photon:
                photons = '%s%s' % (multiplicity, particle.id)
            else:
                label += '%s%s%s' % (sep, multiplicity, particle.id)
                sep = ' + '
        if len(photons) > 0:
            label += '%s%s' % (sep, photons)

        intermediate = ''
        if len(intermediates) > 0:
            intermediate = intermediates[0][0]
        process = ''
        if len(processes) > 0:
            process = processes[0][0]

        lines.append('%s%-32s: %-12s : %-10s: %s' % (indent, label, effectiveThreshold, intermediate, process))

    return lines

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = description)
    parser.add_argument('input', type=str,                                      help='''The protare (i.e., reactionSuite) whose reactions are analyzed.''')

    args = parser.parse_args()
    print('\n'.join(protareProductInfo(args.input, '    ')))
