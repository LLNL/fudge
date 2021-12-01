#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

description = """
This module creates an ACE file from a GNDS file that has been processed for Monte Carlo transport.
"""

__doc__ = description

from fudge import reactionSuite as reactionSuiteModule
from fudge import styles as stylesModule

from brownies.LANL.toACE import reactionSuite
from brownies.LANL.toACE import reaction
from brownies.LANL.toACE import production
from brownies.LANL.toACE import channels
from brownies.LANL.toACE import product
from brownies.LANL.toACE import multiplicity
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
parser.add_argument( 'gnds', type = str,                                            help = 'gnds file to convert to ACE.' )
parser.add_argument( 'output', type = str,                                          help = 'name of the outputted ACE file.' )

args = parser.parse_args( )

gnds = reactionSuiteModule.readXML( args.gnds )
if( args.style is None ) :
    styleOptions = []
    for style in gnds.styles :
        if( isinstance( style, stylesModule.griddedCrossSection ) ) : styleOptions.append( style.label )
    if( len( styleOptions ) == 0 ) : raise Exception( 'GNDS file does not contain Monte Carlo processed data.' )
    if( len( styleOptions )  > 2 ) :
        for style in styleOptions : print( '    %16s | %g' % ( style, gnds.styles[style].temperature ) )
        raise Exception( 'GNDS file does not contain multiple Monte Carlo processed data.' )
    args.style = styleOptions[0]
if( args.style not in gnds.styles ) : raise Exception( 'GNDS file does not contain style "%s".' % args.style )

gnds.toACE( args, args.style, args.output, args.ID, addAnnotation = args.annotate, verbose = args.verbose )
