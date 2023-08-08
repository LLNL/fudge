#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

description = """
This module creates an ACE file from a GNDS file that has been processed for Monte Carlo transport.
"""

__doc__ = description

import pathlib

from fudge import reactionSuite as reactionSuiteModule
from fudge import styles as stylesModule

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

from argparse import ArgumentParser

parser = ArgumentParser( description = description )
parser.add_argument( '-a', '--annotate', action = 'store_true',                     help = 'If present, annotation is added to the ACE file.' )
parser.add_argument( '-i', '--ID', action = 'store', type = int, required = True,   help = 'The evaluation identification.' )
parser.add_argument( '-s', '--style', type = str, default = None,                   help = 'The griddedCrossSection style to convert to ACE.' )
parser.add_argument( '--NILm1', type = int, default = 20,                           help = 'Number of angular equal probable bins for TNSL inelastic scattering. Note this is NIL - 1.' )
parser.add_argument( '--NCL', type = int, default = 20,                             help = 'Number of angular equal probable bins for TNSL elastic scattering.' )
parser.add_argument( '-v', '--verbose', action = 'count', default = 0,              help = 'Verbose mode.' )
parser.add_argument( '--skipURR', action='store_true',                              help = 'Do not write URR probability tables even if they are present in GNDS.' )
parser.add_argument('--skipILF_logic', action='store_true',                         help = 'If present, the URR ILF and ILO flags are only set to -1. This is mainly for testing.')
parser.add_argument( 'gnds', type = str,                                            help = 'gnds file to convert to ACE.' )
parser.add_argument( 'output', type = str,                                          help = 'name of the outputted ACE file.' )

args = parser.parse_args( )

if args.verbose > 0:
    print('Reading GNDS file.')
gnds = reactionSuiteModule.read(args.gnds, lazyParsing=True)

if args.style is None:
    styleOptions = []
    for style in gnds.styles :
        if isinstance( style, stylesModule.GriddedCrossSection ): styleOptions.append( style.label )
    if len( styleOptions ) == 0: raise Exception( 'GNDS file does not contain Monte Carlo processed data.' )
    if len( styleOptions )  > 1:
        print( '    %16s | Temperature (%s)' % ("Style", gnds.styles[styleOptions[0]].temperature.unit))
        for style in styleOptions : print( '    %16s | %g' % ( style, gnds.styles[style].temperature ) )
        raise Exception( 'GNDS file contains multiple Monte Carlo processed data. Please select one of the above styles using option "-s"' )
    args.style = styleOptions[0]
if args.style not in gnds.styles: raise Exception( 'GNDS file does not contain style "%s".' % args.style )
if not isinstance(gnds.styles[args.style], stylesModule.GriddedCrossSection):
    raise Exception("Selected style must be an instance of 'GriddedCrossSection', not %s" % type(gnds.styles[args.style]))

path = pathlib.Path(args.output)
if not path.parent.exists():
    path.parent.mkdir(parents=True)

if args.verbose > 0:
    print('Calling toACE.')
gnds.toACE(args, args.style, args.output, args.ID, addAnnotation=args.annotate, verbose=args.verbose, skipURR=args.skipURR, skipILF_logic=args.skipILF_logic)
