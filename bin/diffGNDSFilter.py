#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from argparse import ArgumentParser

description = """Removes the specified diff tags and their data from the specified diff file. The unfiltered tags are print to the terminal."""

parser = ArgumentParser( description = description )
parser.add_argument( "diffFile",                                        help = """The name of the diff file to filter.""" )
parser.add_argument( "tags", nargs = "+",                               help = """The list of tags. Each tag must be quoted (e.g., "Distribution unspecified - 2".""" )

args = parser.parse_args( )

fIn = open( args.diffFile )
lines = fIn.readlines( )
fIn.close( )

if( ( len( lines ) - 2 ) % 3 != 0 ) : raise Exception( "Input file does not have the correct number of lines to be a valid diff file." )

def checkFileLine( index ) :

    line = lines[index]
    if( ":" not in line ) : raise Exception( "Invalid file line at line %d is not valid." % ( index + 1 ) )
    if( "FILE%s" % ( index + 1 ) != line.split( ":" )[0] ) : raise Exception( "Invalid file line at line %d is not valid." % ( index + 1 ) )
    print( line[:-1] )

checkFileLine( 0 )
checkFileLine( 1 )

for index in range( 2, len( lines ), 3 ) :
    line = lines[index]
    if( ":" not in line ) : raise Exception( "Input tag at line %d is not valid." % ( 3 * index + 2 ) )
    tag = line.split( ":" )[0]

    if( tag in args.tags ) : continue
    print( line[:-1] )
    print( lines[index+1][:-1] )
    print( lines[index+2][:-1] )
