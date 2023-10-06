#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

description = """
This module creates ACE files for each GriddedCrossSection style in a processed GNDS file.
"""

__doc__ = description

import pathlib

from fudge import reactionSuite as reactionSuiteModule
from brownies.LANL.toACE import reactionSuite
from brownies.LANL.toACE import reaction
from brownies.LANL.toACE import production
from brownies.LANL.toACE import channels
from brownies.LANL.toACE import product
from brownies.LANL.toACE import multiplicity
from brownies.LANL.toACE import angularEnergy
from brownies.LANL.toACE import energy
from brownies.LANL.toACE import energyAngular
from brownies.LANL.toACE import KalbachMann

from brownies.LANL.build_xsdir import xsdir_util

from argparse import ArgumentParser

parser = ArgumentParser( description = description )
parser.add_argument( 'gnds', nargs='+',                                     help = 'gnds file(s) to convert to ACE.' )
parser.add_argument( '-o', '--output', default = '.',                       help = 'directory to save ACE files.' )
parser.add_argument( '-i', '--ID', type = int, default = 0,                 help = 'Evaluation ID for starting temperature.' )
parser.add_argument( '--NILm1', type = int, default = 20,                   help = 'Number of angular equal probable bins for TNSL inelastic scattering. Note this is NIL - 1.' )
parser.add_argument( '--NCL', type = int, default = 20,                     help = 'Number of angular equal probable bins for TNSL elastic scattering.' )
parser.add_argument( '-v', '--verbose', action = 'count', default = 0,      help = 'Verbose mode.' )
parser.add_argument( '--skipURR', action='store_true',                      help = 'Do not write URR probability tables even if they are present in GNDS.' )

args = parser.parse_args( )

xsdir = []
path = pathlib.Path(args.output)
if not path.exists():
    path.mkdir(parents=True)

for gndsFile in args.gnds:
    if args.verbose > 0:
        print(f'Reading {gndsFile}')
    gnds = reactionSuiteModule.read(gndsFile, lazyParsing=True)

    tempInfos = [t for t in gnds.styles.temperatures() if t.griddedCrossSection]

    basename = "%s-%s." % (gnds.projectile, gnds.target)
    idNow = args.ID
    for tempInfo in gnds.styles.temperatures():
        if args.verbose > 0:
            print(f'Generating ACE file for temperature {tempInfo.temperature:.3e} MeV/k, id {idNow}')

        output = path / (basename + str(idNow).zfill(2) + ".ace")
        gnds.toACE(args, tempInfo.griddedCrossSection, output, idNow, verbose=args.verbose, skipURR=args.skipURR)
        xsdir.append(xsdir_util.xsdir_entry(output))
        idNow += 1
        if idNow > 99:
            print(f'File {args.gnds} contains {len(tempInfos)} temperatures. Only wrote first {99 - args.ID} starting with {args.ID}')
            break

with open(path / 'xsdir', 'w') as fout:
    fout.write('\n'.join(xsdir))
