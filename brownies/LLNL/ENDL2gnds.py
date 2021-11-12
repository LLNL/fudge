#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
import os
from brownies.legacy.endl.endlProject import endlProject
import brownies.legacy.toENDF6.toENDF6 # adds 'toENDF6' methods to GNDS classes

from xData import formatVersion as formatVersionModule
import argparse

parser = argparse.ArgumentParser( description = "Translate ENDL evaluations into GNDS, and optionally also ENDF-6." )
parser.add_argument( 'ZA', type = str, help = "za of target to translate" )
parser.add_argument( '-l', '--library', default = 'endl2011.0',
            help = "ENDL library (default=endl2011.0). May be full path or identifier" )
parser.add_argument( '-p', '--projectile', default = 'n',
            type = lambda val: int(val) if val.isdigit() else val,
            help = "Choose projectile (default='n'). May be string (one of n,p,d,t,He3,a,g) or yi number" )
parser.add_argument( '--includeAverageProductData', action='store_true',
            help = "Include I=10 and 13 (average outgoing energy and momentum) data" )
parser.add_argument( '-6', '--toENDF6', action = 'store_true',
            help = "After creating GNDS, also translate to ENDF-6" )
parser.add_argument( '-o', '--output', default='endl2gnds', help='prefix for resulting .endf and .xml files' )
parser.add_argument( '-v', '--version', default = '', # '1.0.0',
            help = "Evaluation version number. Should be of form 'major.minor.patchlevel' (i.e. 2011.0.1)" )
parser.add_argument("--formatVersion", default=formatVersionModule.default, choices=formatVersionModule.allowed,
                    help="Specifies the format for the outputted GNDS file. ")

if( __name__ == '__main__' ) :

    args = parser.parse_args( )

    e = endlProject( args.library, projectile = args.projectile, readOnly = True )

    if args.ZA[-1].isalpha():   # isomer, e.g. 95242m
        za = e.readZA( args.ZA[:-1], suffix = args.ZA[-1] )
    else:
        za = e.readZA( args.ZA )
    za.read ()
    if not args.includeAverageProductData:
        za.removeFile(I=10)
        za.removeFile(I=13)

    rs, covars = za.toGNDS( evaluationLibrary = os.path.basename( args.library), evaluationVersion = args.version,
                            formatVersion = args.formatVersion )

    if covars is not None:
        from fudge import externalFile
        from LUPY import checksums
        cov_output = args.output + "-covar.xml"
        covars.externalFiles.add( externalFile.externalFile( "reactions", path=args.output+".xml" ) )
        covars.saveToFile( cov_output )

        sha1sum = checksums.sha1sum.from_file(cov_output)
        rs.externalFiles.add( externalFile.externalFile( "covariances", path=cov_output, checksum=sha1sum ) )

    rs.saveToFile( args.output + ".xml" )

    if( args.toENDF6 ) :
        r.convertUnits( {'MeV':'eV'} )
        with open( args.output + ".endf", "w" ) as fout :
            fout.write( r.toENDF6( style = "eval", flags = { 'verbosity' : 32 } ) )
