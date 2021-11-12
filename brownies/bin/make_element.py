#! /usr/bin/env python

import fudge.core.utilities.abundance as abund
import argparse
import json
import collections
import subprocess
import os.path

from PoPs.groups import misc as miscPoPsModule

__doc__="""
Construct a JSON file that collects all of the ENDF files for a specific element together, with the
correct abundances, for a specific sublibrary.  The resulting JSON file can be plotted with plot_evaluation.py
or certain metrics can be computed with inter.py.

Only transportable particles are allowed as projectiles.
"""

SUBLIBS=['n', 'p', 'd', 't', 'a', 'h', 'g'] # doesn't make sense for atomic, TSL, Decay, standards and FPY sublibraries
lineStyles = [':', "--", "-.", '-']

def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-v', dest='verbose', default=False, action='store_true', help="Enable verbose output.")
    parser.add_argument('-o', dest='outFile', default=None, help="Save to this file")
    parser.add_argument('-c', dest='color', default='blue', help='Color of the cross sections, etc., when plotted')
    parser.add_argument('-s', dest='sublibrary', default=None, choices=SUBLIBS, help='ENDF sublibrary prefix/projectile')
    parser.add_argument('-p', dest='path', default=None, help='Path to ENDF files for look for isotopic ENDF files')
    parser.add_argument('-h', dest='include_hash', default=False, action='store_true', help='Include the git hash in the generated JSON file')
    parser.add_argument('Z', type=int, help="Z of target to work with")
    return parser.parse_args()


def endf_file(_s, _Z, _sym, _A):
    return _s + '-' + str(_Z).zfill(3) + '_' + _sym + '_' + str(_A).zfill(3) + '.endf'


if __name__ == "__main__":
    args=parse_args()
    try:
        atable=abund.getElementsNaturalIsotopes(args.Z)
    except KeyError:
        print("Z=%i doesn't correspond to a naturally occuring isotope" % args.Z)
        exit()

    symbol = miscPoPsModule.symbolFromZ[args.Z]
    answer = collections.OrderedDict()
    answer["target"] = "%s0" % symbol
    if args.sublibrary is not None:
        answer["projectile"] = args.sublibrary
    answer["legend"] = "$^{nat}$%s" % miscPoPsModule.symbolFromZ[args.Z]
    answer["lineStyle"] = "-"
    answer["lineColor"] = args.color
    answer["lineWidth"] = 2
    answer["eMin"] = "1e-5 eV"
    answer["eMax"] = "20 MeV"
    answer['isotopes'] = collections.OrderedDict()
    for i, A in enumerate(atable):
        ZA = miscPoPsModule.idFromZAndA(args.Z, A)
        answer['isotopes'][ZA] = collections.OrderedDict()
        answer['isotopes'][ZA]["atomicFraction"] = atable[A][0]/100.0
        answer['isotopes'][ZA]["atomicFraction_unc"] = atable[A][1]/100.0
        if args.path is not None and os.path.exists(args.path):
            efile = endf_file(args.sublibrary, args.Z, symbol, A)
            fullpath = "%s/%s" % (args.path, efile)
            if os.path.exists(fullpath):
                answer['isotopes'][ZA]["pathToFile"] = fullpath
                if args.include_hash:
                    try:
                        fullpath_hash = subprocess.check_output(["git", "hash-object", fullpath])
                        answer['isotopes'][ZA]["githash"] = fullpath_hash.decode(encoding='utf-8').strip()
                    except CalledProcessError:
                        print("Could not hash?")
            else:
                raise FileNotFoundError(fullpath)
        answer['isotopes'][ZA]["legend"] = ZA
        answer['isotopes'][ZA]["lineStyle"] = lineStyles[i%len(lineStyles)]
        answer['isotopes'][ZA]["lineWidth"] = 1
        answer['isotopes'][ZA]["lineColor"] = args.color

    if args.outFile is None:
        print(json.dumps(answer, indent=4))
    else:
        with open(args.outFile, mode='w') as out:
            out.write(json.dumps(answer, indent=4))
