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

description = """This script removed all processed data (except resonance reconstructed data) for the specified GNDS protare file."""

parser = argparse.ArgumentParser(description=description, allow_abbrev=False)
singleProtareArguments = argumentsForScriptsModule.SingleProtareArguments(parser)
parser.add_argument('output', nargs='?', default=None,                  help='Output file name.')

args = parser.parse_args()

protare = singleProtareArguments.protare(args)

stylesToRemove = []
preProcessingStyles = protare.styles.preProcessingStyles()
for style in protare.styles :
    if not isinstance(style, preProcessingStyles):
        stylesToRemove.append(style.label)
protare.removeStyles(stylesToRemove)

output = args.output
if output is None: output = pathlib.Path(protare.sourcePath).with_suffix('.culled.xml')
protare.saveToFile(output)
