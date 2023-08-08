#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import argparse

from fudge import GNDS_formatVersion as GNDS_formatVersionModule
from fudge import externalFile as externalFileModule
from LUPY import checksums as checksumsModule
from PoPs import IDs as IDsPoPsModule

from brownies.legacy.endl.endlProject import endlProject
import brownies.legacy.toENDF6.toENDF6 # adds 'toENDF6' methods to GNDS classes

parser = argparse.ArgumentParser( description = "Translate ENDL evaluations into GNDS, and optionally also ENDF-6." )
parser.add_argument( 'ZA', type = str, help = "za of target to translate" )
parser.add_argument( '-l', '--library', default = 'endl2011.0',
            help = "ENDL library (default=endl2011.0). May be full path or identifier" )
parser.add_argument( '-p', '--projectile', default = IDsPoPsModule.neutron,
            type = lambda val: int(val) if val.isdigit() else val,
            help = "Choose projectile (default='n'). May be string (one of n,p,d,t,He3,a,g) or yi number" )
parser.add_argument( '--includeAverageProductData', action='store_true',
            help = "Include I=10 and 13 (average outgoing energy and momentum) data" )
parser.add_argument( '-6', '--toENDF6', action = 'store_true',
            help = "After creating GNDS, also translate to ENDF-6" )
parser.add_argument( '-o', '--output', default='endl2gnds', help='prefix for resulting .endf and .xml files' )
parser.add_argument( '-v', '--version', default = '', # '1.0.0',
            help = "Evaluation version number. Should be of form 'major.minor.patchlevel' (i.e. 2011.0.1)" )
parser.add_argument("--formatVersion", default=GNDS_formatVersionModule.default, choices=GNDS_formatVersionModule.allowed,
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

    outFile = args.output + ".xml"
    if covars is not None:
        cov_output = args.output + "-covar.xml"
        covars.externalFiles.add( externalFileModule.ExternalFile( "reactions", path=os.path.basename(outFile) ) )
        covars.saveToFile( cov_output )

        sha1sum = checksumsModule.Sha1sum.from_file(cov_output)
        rs.externalFiles.add( externalFileModule.ExternalFile( "covariances", path=os.path.basename(cov_output), checksum=sha1sum ) )

    rs.saveToFile( outFile )

    if( args.toENDF6 ) :
        r.convertUnits( {'MeV':'eV'} )
        with open( args.output + ".endf", "w" ) as fout :
            fout.write( r.toENDF6( style = "eval", flags = { 'verbosity' : 32 } ) )
