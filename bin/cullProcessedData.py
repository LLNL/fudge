#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse
import pathlib

from LUPY import argumentsForScripts as argumentsForScriptsModule
import fudge.styles as stylesModule

summaryDocString__FUDGE = """This script removed all processed data (except resonance reconstructed data) for the specified GNDS protare file."""

description = summaryDocString__FUDGE + """ If option "--outputDir" is present, the directory is striped from the output path and replaced 
with the specified path."""

parser = argparse.ArgumentParser(description=description, allow_abbrev=False)
singleProtareArguments = argumentsForScriptsModule.SingleProtareArguments(parser)
parser.add_argument('outputPath', default=None, type=pathlib.Path,             help='Output file name.')
parser.add_argument('--outputDir', action='store', default=None, type=pathlib.Path,     help='The output directory to write the output file to.')

args = parser.parse_args()

protare = singleProtareArguments.protare(args)
protare.cullProcessedData()

outputPath = args.outputPath
if outputPath is None:
    outputPath = pathlib.Path(protare.sourcePath).with_suffix('.culled.xml')
if args.outputDir is not None:
    outputPath = args.outputDir / outputPath.name

protare.saveToFile(outputPath)
