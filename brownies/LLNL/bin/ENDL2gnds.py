#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import glob
from argparse import ArgumentParser

from PoPs import database as databaseModule
from PoPs.chemicalElements import misc as PoPsChemicalElementsMiscModule
from fudge import GNDS_formatVersion as GNDS_formatVersionModule
from fudge import enums as enumsModule

from fudge import externalFile as externalFileModule
from LUPY import checksums as checksumsModule

from brownies.legacy.endl.endlZA import endlZA as endlZAClass
from brownies.legacy.endl import bdfls as bdflsModule
from brownies.legacy.endl import fudgeParameters as fudgeParametersModule

fudgeParametersModule.VerboseMode = 0

defaultEvaluationName = 'Test'
defaultEvaluationVersion = '1.0'

usage = '''
Converts the data in an endlZA directory (e.g., './yi01/za008016') to GNDS. The path to the endlZA directory
must end with a projectile-directory (e.g., 'yi01') followed by a ZA-directory (e.g., 'za008016').'''

parser = ArgumentParser(description=usage)
parser.add_argument('ZAPath', type=str,                                                help='path to an endl ZA.')
parser.add_argument('-i', '--includeAverageProductData', action='store_true',          help='Include ENDL average product data (i.e., I 10, 11, and 13 data) in conversion.')
parser.add_argument('-e', '--evaluationName', action='store', default=defaultEvaluationName,
                                                                                        help='Evaluation name. Default is "%s".' % defaultEvaluationName)
parser.add_argument('-o', '--output', action='store', default=None,                    help='Output file full path. If not entered, ZA-directory with "xml" appended.')
parser.add_argument('-r', '--rename', action='store_true',                             help='Apply standard name (proj-ZZZ_Sym_AAA) to output file(s).')
parser.add_argument('-v', '--verbose', action='count', default=0,                      help='Verbose mode.')
parser.add_argument('-V', '--version', action='store', default=defaultEvaluationVersion,
                                                                                        help='Evaluation version. Default is "%s".' % defaultEvaluationVersion)
parser.add_argument('--formatVersion', default=GNDS_formatVersionModule.default, choices=GNDS_formatVersionModule.allowed,
                                                                                        help='Specifies the format for the outputted GNDS file.')

args = parser.parse_args()

crossSectionFiles = glob.glob(args.ZAPath + os.sep + 'yo00c??i000s00?')
if len(crossSectionFiles) == 0:
    raise Exception('No ENDL cross section file in endl ZA directory given.')

yiPath, ZA = os.path.split(args.ZAPath)
database, yi = os.path.split(yiPath)
if database == '':
    database = './'

asciiPath = os.path.dirname(os.path.dirname(yiPath))
bdflsFile = None
bdflsFileName = os.path.join(asciiPath, 'bdfls')
if os.path.exists(bdflsFileName):
    bdflsFile = bdflsModule.getBdflsFile(bdflsFileName)

endlZA = endlZAClass(ZA, yi, database, readOnly=True, bdflsFile=bdflsFile)
endlZA.read()
GNDS, covars = endlZA.toGNDS(args.evaluationName, args.version, formatVersion=args.formatVersion,
                              excludeAverageProductData=not args.includeAverageProductData, verbose=args.verbose)

output = args.output
if output is None:
    output = ZA + '.xml'

if args.rename:
    if enumsModule.Interaction.LLNL_TNSL == GNDS.interaction:
        filename = 'LLNL_TNSL-%s.xml' % GNDS.target
    else:
        projectile = {1: 'n', 2: 'p', 3: 'd', 4: 't', 5: 'h', 6: 'a', 7: 'photoat', 9: 'e'}[endlZA.yi]
        Z, A = endlZA.Z, endlZA.A
        if Z != 0:
            if 99120 <= endlZA.ZA <= 99125:
                filename = '%s-ENDL_fissionProduct_%s.xml' % (projectile, endlZA.ZA)
            else:
                symbol = PoPsChemicalElementsMiscModule.symbolFromZ[Z]
                metaStable = ''
                if endlZA.suffix != '':
                    metaStable = 'm1'
                filename = '%s-%.3d_%s_%.3d%s.xml' % (projectile, Z, symbol, A, metaStable)
        else:
            if Z == 0 and A == 1:
                symbol = 'n'
            else:
                raise Exception('Do not know how to create file name for target symbol = "%s".' % endlZA.targetSymbol)
            filename = '%s-%.3d_%s_%.3d.xml' % (projectile, Z, symbol, A)

    GNDS.saveAllToFile(filename)
else:
    if covars is not None:
        cov_output = output.replace('.xml', '-covar.xml')
        covars.externalFiles.add(externalFileModule.ExternalFile('reactions', path=os.path.basename(output)))
        sha1sum = None
        if not args.rename:
            covars.saveToFile(cov_output)
            sha1sum = checksumsModule.Sha1sum.from_file(cov_output)

        GNDS.externalFiles.add(externalFileModule.ExternalFile('covariances', path=os.path.basename(cov_output), checksum=sha1sum))

    GNDS.saveToFile(output)
