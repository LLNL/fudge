#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import argparse

from xData import formatVersion as formatVersionModule
from PoPs import database as databaseModule
from LUPY import GNDSType as GNDSTypeModule

extensionDefault = '.g2g'

description1 = """Read a GNDS file into Fudge, then write back to the GNDS/xml format.  Intent is to test 
for errors during reading or writing.  Sample usage: python gnds2gnds.py n-001_H_001.xml 
If file n-001_H_001-cov.xml (or -covar.xml) exists, covariances will automatically be read and re-written.
The output file's path and extension can be set via the -p and -e options respectively. If the 'output' argument
is specified, the path and extension arguments are ignored and it is used as the name of the output file.
"""

__doc__ = description1

parser = argparse.ArgumentParser( description1 )
parser.add_argument( 'input',                                                           help = 'GNDS and/or PoPs file to translate.' )
parser.add_argument( 'output', nargs = '?', default = None,                             help = 'The name of the output file.' )
parser.add_argument( '-e', '--extension', default = extensionDefault,                   help = 'The file extension to add to the output file. Default = "%s"' % extensionDefault )
parser.add_argument( '-o', '--outline', default = False, action = 'store_true',         help = 'The outputted GNDS files are written in outline mode.' )
parser.add_argument( '-p', '--path', default = None,                                    help = 'Path to write the file to. If absent, sent to same location as input.' )
parser.add_argument( '--energyUnit', type = str, default = None,                        help = 'Convert all energies in the gnds file to this unit.' )
parser.add_argument( '--hybrid', default = False, action = 'store_true',                help = 'Write out hybrid XML/HDF5' )
parser.add_argument( '--minLength', type = int, default = 50,                           help = 'Min length of datasets to store in HDF5' )
parser.add_argument( '--flatten', default = False, action = 'store_true',               help = 'Use flattened arrays in HDF5' )
parser.add_argument( '--compress', default = False, action = 'store_true',              help = 'Use gzip + shuffle compression in HDF5' )
parser.add_argument( '--formatVersion', default = formatVersionModule.default, choices = formatVersionModule.allowed,
                                                                                        help = 'Specifies the GNDS format for the outputted file.  Default = "%s".' % formatVersionModule.default )

if( __name__ == '__main__' ) :

    args = parser.parse_args( )

    fileName = args.input

    covariances = []
    name, dummy = GNDSTypeModule.type( fileName )
    if( name == databaseModule.database.moniker ) :
        gnds = GNDSTypeModule.read( fileName )
    else :
        gnds = GNDSTypeModule.read( fileName )
        if hasattr(gnds, 'loadCovariances'):
            covariances = gnds.loadCovariances()

    if( args.energyUnit is not None ) :
        gnds.convertUnits( { 'MeV' : args.energyUnit, 'eV' : args.energyUnit } )
        for covarianceSuite in covariances:
            covarianceSuite.convertUnits( { 'MeV' : args.energyUnit, 'eV' : args.energyUnit } )

    output = args.output
    path = args.path
    extension = args.extension

    if( output is None ) :
        if( path is None ) : path = os.path.dirname( fileName )
        output = os.path.join( path, os.path.basename( fileName ) ) + extension

    if args.hybrid:
        kwargs = {'minLength':args.minLength, 'flatten':args.flatten, 'compress':args.compress}
        gnds.saveToHybrid( output, formatVersion = args.formatVersion, **kwargs )
        # FIXME support hybrid storage for covariances
    else:
        gnds.saveToFile( output, outline = args.outline, xs_pdf_cdf1d_singleLine = True, formatVersion = args.formatVersion )

    for covarianceSuite in covariances:
        output = os.path.join( os.path.dirname( output ), os.path.basename( covarianceSuite.sourcePath ) ) + extension
        covarianceSuite.saveToFile( output, formatVersion = args.formatVersion )
