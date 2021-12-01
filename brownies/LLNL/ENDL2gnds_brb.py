#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import glob
from argparse import ArgumentParser

from brownies.legacy.endl.endlZA import endlZA as endlZAClass
from brownies.legacy.endl import bdfls as bdflsModule
from brownies.legacy.endl import fudgeParameters as fudgeParametersModule

from xData import formatVersion as formatVersionModule

fudgeParametersModule.VerboseMode = 0

defaultEvaluationName = "Test"
defaultEvaluationVersion = "1.0"

usage = """
Converts the data in an endlZA directory (e.g., './yi01/za008016') to GNDS. The path to the endlZA directory
must end with a projectile-directory (e.g., 'yi01') followed by a ZA-directory (e.g., 'za008016').
"""

parser = ArgumentParser( description = usage )
parser.add_argument( 'ZAPath', type = str,                                          help = 'path to an endl ZA.' )
parser.add_argument( '-i', '--includeAverageProductData', action = 'store_true',    help = 'Include ENDL average product data (i.e., I 10, 11, and 13 data) in conversion.' )
parser.add_argument( '-e', '--evaluationName', action = 'store', default = defaultEvaluationName,
                                                                                    help = 'Evaluation name. Default is "%s".' % defaultEvaluationName )
parser.add_argument( '-o', '--output', action = 'store', default = None,            help = 'Output file full path. If not entered, ZA-directory with "xml" appended.' )
parser.add_argument( '-v', '--verbose', action = 'count', default = 0,              help = 'Verbose mode.' )
parser.add_argument( '-V', '--version', action = 'store', default = defaultEvaluationVersion,
                                                                                    help = 'Evaluation version. Default is "%s".' % defaultEvaluationVersion )
parser.add_argument("--formatVersion", default=formatVersionModule.default, choices=formatVersionModule.allowed,
                    help="Specifies the format for the outputted GNDS file. ")

args = parser.parse_args( )

crossSectionFiles = glob.glob( args.ZAPath + os.sep + 'yo00c??i000s00?' )
if( len( crossSectionFiles ) == 0 ) : raise Exception( 'No ENDL cross section file in endl ZA directory given.' )

yiPath, ZA = os.path.split( args.ZAPath )
database, yi = os.path.split( yiPath )
if( database == '' ) : database = './'

asciiPath = os.path.dirname( os.path.dirname( yiPath ) )
bdflsFile = None
bdflsFileName = os.path.join( asciiPath, 'bdfls' )
if( os.path.exists( bdflsFileName ) ) : bdflsFile = bdflsModule.getBdflsFile( bdflsFileName )

endlZA = endlZAClass( ZA, yi, database, readOnly = True, bdflsFile = bdflsFile )
endlZA.read( )
GNDS, covars = endlZA.toGNDS( args.evaluationName, args.version, formatVersion = args.formatVersion,
                              excludeAverageProductData = not( args.includeAverageProductData ), verbose = args.verbose )

output = args.output
if( output is None ) : output = ZA + '.xml'

if covars is not None:
    from fudge import externalFile
    from LUPY import checksums
    cov_output = output.replace(".xml", "-covar.xml")
    covars.externalFiles.add( externalFile.externalFile( "reactions", path=output) )
    covars.saveToFile( cov_output )

    sha1sum = checksums.sha1sum.from_file(cov_output)
    GNDS.externalFiles.add( externalFile.externalFile( "covariances", path=cov_output, checksum=sha1sum ) )

GNDS.saveToFile( output )

