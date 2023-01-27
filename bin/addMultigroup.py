#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os

from pqu import PQU as PQUModule
from xData import axes as axesModule
from xData import values as valuesModule

from fudge.processing import group as groupModule

energyUnitDefault = "MeV"

from argparse import ArgumentParser

summaryDocStringFUDGE = '''Adds a multi-group boundary definition (i.e., label and the multi-group boundaries) to a groups file (e.g., groups.xml).'''

description = """Read a list of energy boundaries from a specified file or the command line, and adds it to a multi-group file."""

parser = ArgumentParser( description = description )
parser.add_argument( "label",                                           help = """The label to assign to the energy boundaries.""" )
parser.add_argument( "boundaries",                                      help = """The name of a  file containing the list of energy boundaries or if the "--values" option is specified the list of boundaries in a quoted string (.e.g., "1e-10 1e-5 1 2 3 5 10 20").""" )
parser.add_argument( "output",                                          help = """Name of the outputted multi-group file. If the "input" argument is not present, then the output file, if present, is used for the input file. Otherwise, no input file is used.""" )
parser.add_argument( "input", nargs = "?", default = None,              help = """If present, the file to read existing multi-group data from. The new boundaries are appended to the boundaries in this file and all are written to the output file.""" )
parser.add_argument( "--override", action = "store_true",               help = """If option present and label exists in inputted multi-group file, replace multi-group with new energy boundaries; otherwise, execute a raise.""" )
parser.add_argument( "-u", "--unit", default = energyUnitDefault,       help = """The unit for the energy boundaries. Default is %s.""" % energyUnitDefault )
parser.add_argument( "--values", action = "store_true",                 help = """Qualifies how the "boundaries" argument is interpreted.""" )

args = parser.parse_args( )

unit = PQUModule.PQU( 1, args.unit )
if( not( unit.isEnergy( ) ) ) : raise TypeError( "Unit must be an energy unit." )

if( args.values ) :
    values = args.boundaries
else :
    fIn = open( args.boundaries )
    lines = fIn.readlines( )
    fIn.close( )
    values = " ".join( lines )

boundaries = list( map( float, values.split( ) ) )
if( sorted( boundaries ) != boundaries ) : raise ValueError( "Boundaries not in ascending order." )
boundaries = valuesModule.Values( boundaries )

grid = axesModule.Grid( "energy", 0, args.unit, "boundaries", boundaries )
group = groupModule.Group( args.label, grid )

input = args.input
if( input is None ) :
    if( os.path.exists( args.output ) ) : input = args.output

if( input is not None ) :
    multiGroups = groupModule.Groups.readXML_file(input)
else :
    multiGroups = groupModule.Groups( )

if( group.label in multiGroups ) :
    if( not( args.override ) ) : raise ValueError( """Label "%s" already in multi-group file.""" % group.label )
    multiGroups.replace( group )
else :
    multiGroups.add( group )

multiGroups.saveToFile( args.output )
