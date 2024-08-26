#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse

from fudge import map as mapModule

summaryDocString__FUDGE = '''Converts a GNDS map file from one format to another.'''

description = """This script converts the specified map file into the requested GNDS format."""

parser = argparse.ArgumentParser(description=description)
parser.add_argument('input', type=str,                      help='The name map file to convert.')
parser.add_argument('output', type=str,                     help='The name of the output map file.')
parser.add_argument('--format', action='store', choices=mapModule.FormatVersion.allowed, default=mapModule.FormatVersion.default,
                                                            help='Specifies the GNDS format of the map file to be written.')

args = parser.parse_args()

map = mapModule.read(args.input)

map.saveToFile(args.output, format=args.format)
