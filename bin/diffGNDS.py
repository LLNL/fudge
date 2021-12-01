#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

description = """Diffs two GNDS files. Prints the diff information to the terminal. 
Converts the energy unit of the second GNDS file to the first GNDS file if needed.
The GNDS files must be formatted in the same GNDS format version or a raise will be executed."""

__doc__ = description

import os
from argparse import ArgumentParser

from fudge import reactionSuite as reactionSuiteModule

parser = ArgumentParser( description = description )
parser.add_argument( "GNDSFile1",                                       help = "The name of the first GNDS file." )
parser.add_argument( "GNDSFile2",                                       help = "The name of the second GNDS file." )

args = parser.parse_args( )

protare1 = reactionSuiteModule.readXML( args.GNDSFile1 )
protare2 = reactionSuiteModule.readXML( args.GNDSFile2 )

if( protare1.format != protare2.format ) : raise Exception( 'GNDS formats not the same: "%s" vs. "%s".' % ( protare1.format, protare2.format ) )

if( protare1.domainUnit != protare2.domainUnit ) : protare2.convertUnits( { protare2.domainUnit : protare1.domainUnit } )

print( "FILE1: %s" % os.path.realpath( args.GNDSFile1 ) )
print( "FILE2: %s" % os.path.realpath( args.GNDSFile2 ) )
print( protare1.diff( protare2 ) )
