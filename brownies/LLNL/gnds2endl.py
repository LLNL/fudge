#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from argparse import ArgumentParser

from brownies.legacy.toENDL import reactionSuite
from LUPY import argumentsForScripts as argumentsForScriptsModule

# from brownies.legacy.endl import fudgeParameters as fudgeParametersModule
# fudgeParametersModule.VerboseMode = 0

usage = """
Converts the data in a GNDS file to an endlZA directory (e.g., './yi01/za008016'). This will only work 
if the GNDS data came from ENDL data. Ergo, it will not work if the GNDS data are from an ENDF file.
"""

parser = ArgumentParser( description = usage )

singleProtareArguments = argumentsForScriptsModule.SingleProtareArguments( parser )

parser.add_argument( '-o', '--outputDir', action = 'store', default = None, required = True,    help = 'Output file full path.' )
parser.add_argument( '-v', '--verbose', action = 'count', default = 0,                          help = 'Sets verbosity.' )

args = parser.parse_args( )

protare = singleProtareArguments.protare( args )

ENDL = protare.toENDL( args.outputDir, args.verbose )
