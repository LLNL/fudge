#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
import os
import traceback
import argparse

from fudge import GNDS_formatVersion as GNDS_formatVersionModule
from fudge import GNDS_file as GNDS_fileModule
from PoPs import database as databaseModule

extensionDefault = '.g2g'

description1 = """Read a GNDS file into Fudge, then write back to the GNDS/xml format.  Intent is to test 
for errors during reading or writing.  Sample usage: python gnds2gnds.py n-001_H_001.xml 
If file n-001_H_001-cov.xml (or -covar.xml) exists, covariances will automatically be read and re-written.
The output file's path and extension can be set via the -p and -e options respectively. If the 'output' argument
is specified, the path and extension arguments are ignored and it is used as the name of the output file.
"""

__doc__ = description1

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description1)
    parser.add_argument('input',                                                        help='GNDS and/or PoPs file to translate.')
    parser.add_argument('output', nargs='?', default=None,                              help='The name of the output file.')
    parser.add_argument('h5output', nargs='?', default=None,                            help='The name of the hdf5 output file, for use with the --hybrid option.')
    parser.add_argument('-e', '--extension', default=extensionDefault,                  help='The file extension to add to the output file. Default is "%s"' % extensionDefault)
    parser.add_argument('-o', '--outline', default=False, action='store_true',          help='The outputted GNDS files are written in outline mode.')
    parser.add_argument('-p', '--path', default=None,                                   help='Path to write the file to. If absent, sent to same location as input.')
    parser.add_argument('--energyUnit', type=str, default=None,                         help='Convert all energies in the gnds file to this unit.')
    parser.add_argument('--hybrid', default=False, action='store_true',                 help='Write out hybrid XML/HDF5')
    parser.add_argument('--skipCovariances', action='store_true',                       help='If present, any covariance files in are not written.')
    parser.add_argument('--formatVersion', default=GNDS_formatVersionModule.default, choices=GNDS_formatVersionModule.allowed,
                                                                                        help='Specifies the GNDS format for the outputted file.  Default is "%s".' % GNDS_formatVersionModule.default)
    parser.add_argument('--significantDigits', type=int, default=15,                    help='Sets the number of significantDigits used by the Values.toXML_strList method.')
    parser.add_argument('--traceback', action='store_true', default=False,              help='Print traceback on exception')

    args = parser.parse_args()

    fileName = args.input

    covariances = []
    name, dummy = GNDS_fileModule.type(fileName)
    if name == databaseModule.Database.moniker:
        gnds = GNDS_fileModule.read(fileName)
    else:
        gnds = GNDS_fileModule.read(fileName)
        if not args.skipCovariances:
            try:
                if hasattr(gnds, 'loadCovariances'): covariances = gnds.loadCovariances()
            except:
                print('WARNING: could not load covariance file(s).')
                if args.traceback:
                    exc_type, exc_value, exc_traceback = sys.exc_info()
                    traceback.print_tb(exc_traceback)

    if hasattr(gnds, 'convertUnits'):
        if args.energyUnit is not None:
            gnds.convertUnits({'MeV': args.energyUnit, 'eV': args.energyUnit})
            for covarianceSuite in covariances:
                covarianceSuite.convertUnits({'MeV': args.energyUnit, 'eV': args.energyUnit})

    output = args.output
    path = args.path
    extension = args.extension

    if output is None:
        if path is None: path = os.path.dirname(fileName)
        output = os.path.join(path, os.path.basename(fileName)) + extension

    if args.hybrid:
        gnds.saveToHybrid(output, hdfName=args.h5output, formatVersion=args.formatVersion, minLength=3, flatten=True)
        # FIXME support hybrid storage for covariances? Reenable compression?
    else:
        assert args.h5output is None, 'HDF5 output file listed without the --hybrid option!'
        gnds.saveToFile(output, outline=args.outline, xs_pdf_cdf1d_singleLine=True, formatVersion=args.formatVersion, 
                significantDigits=args.significantDigits)

    for covarianceSuite in covariances:
        output = os.path.join(os.path.dirname(output), os.path.basename(covarianceSuite.sourcePath)) + extension
        covarianceSuite.saveToFile(output, formatVersion=args.formatVersion)
