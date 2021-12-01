#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

# This script takes a gnds/XML file, reads it in and rewrites it to a file as an ENDF6 file
# (in the currently directory) with the extension '.endf' added.

import os
import argparse

from LUPY import GNDSType as GNDSTypeModule
from fudge import styles as stylesModule

# add 'toENDF6' methods to GNDS classes:
import brownies.legacy.toENDF6.toENDF6

description = """Translate one or more GNDS files to ENDF-6.
Sample use: python gnds2endf.py n-001_H_001.xml n-001_H_002.xml ...
If file n-001_H_001-covar.xml exists, covariances will automatically be read from it.
Use option '-o' to change directory where resulting ENDF-6 files will be written."""

__doc__ = description

parser = argparse.ArgumentParser( description )
parser.add_argument( '-o', '--outputDir', default = '.',                    help = 'Directory to save translated ENDF-6 files' )
parser.add_argument( '-l', '--lineNumbers', action = 'store_true',          help = 'Add line numbers' )
parser.add_argument( '--writeReconstructed', action = 'store_true',         help = 'If present, the resonance reconstructed cross section data are written, if present, instead of the background cross section data.' )
parser.add_argument( '-v', '--verbose', action = 'count', default = 0,      help = 'Determines verbosity level, the more the merrier' )
parser.add_argument( 'gnds', nargs = "+",                                   help = "GNDS file(s) to translate" )

if( __name__ == '__main__' ) :
    args = parser.parse_args( )
    if( not os.path.exists( args.outputDir ) ) : os.makedirs( args.outputDir )

    for gndsFile in args.gnds :

        gndsCov = None
        name, dummy = GNDSTypeModule.type( gndsFile )

        gnds = GNDSTypeModule.read( gndsFile )
        if hasattr(gnds, 'loadCovariances'):
            covariances = gnds.loadCovariances()
            if len(covariances) == 1:
                gndsCov = covariances[0]
            elif len(covariances) > 1:
                raise NotImplementedError("Converting multiple covarianceSuites back to ENDF-6")

        styleLabel = gnds.styles[0].label
        if( hasattr( gnds.styles, 'preProcessingChains' ) ) :
            preProcessingChains = gnds.styles.preProcessingChains( ends = True )
            if( len( preProcessingChains ) != 1 ) : raise Exception( 'Currently, only one chain is supported.' )
            if( args.writeReconstructed ) :
                for style in preProcessingChains[0] :
                    if( isinstance( style, stylesModule.crossSectionReconstructed ) ) : styleLabel = style.label

        with open( os.path.join( args.outputDir, os.path.basename( gndsFile ) + '.endf' ), 'w' ) as fout :
            fout.write( gnds.toENDF6( styleLabel, flags = { 'verbosity' : args.verbose * 10 }, covarianceSuite = gndsCov, lineNumbers = args.lineNumbers ) )
