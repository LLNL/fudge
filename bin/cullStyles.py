#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse

from fudge import reactionSuite as reactionSuiteModule, GNDS_file as GNDS_fileModule

summaryDocStringFUDGE = """This module reads a GNDS file and removes (culls) all styles from each component but for a selected style."""

description1 = """This module reads a GNDS file and removes (culls) all styles from each component but for a selected style.
If the selected style is not present for a component, the nearest derived style for the selected style is kept.
"""

__doc__ = description1

parser = argparse.ArgumentParser( description1 )
group1 = parser.add_mutually_exclusive_group( required = True )
group1.add_argument( '-l', '--list', action = 'store_true',                                 help = 'List all styles and exit (i.e., does not cull any data).' )
group1.add_argument( '-s', '--style', action = 'store',                                     help = 'The label of the style to keep (i.e., the selected style).' )

parser.add_argument( 'gnds',                                                                help = 'Name of GNDS file to cull.' )
parser.add_argument( 'output',                                                              help = 'Name of output file.' )
parser.add_argument( '--removeDocumentation', action = 'store_true',                        help = 'Remove the documentation section.' )
parser.add_argument( '--removeApplicationData', action = 'store_true',                      help = 'Remove the applicationData section.' )

if __name__ == '__main__' :

    args = parser.parse_args( )

    name, dummy = GNDS_fileModule.type(args.gnds)

    if( name == reactionSuiteModule.ReactionSuite.moniker ) :
            gnds = GNDS_fileModule.read(args.gnds)
    else :
        raise Exception( 'Current only GNDS reactionSuited file is supported: name = %s.' % name )

    if( args.list ) :
        fmt = "%%-%ss -> %%s" % max( [ len( style.label ) for style in gnds.styles ] )
        for style in gnds.styles : print( fmt % ( style.label, style.moniker ) )
    else :
        gnds.cullStyles( args.style, removeDocumentation = args.removeDocumentation, removeApplicationData = args.removeApplicationData )
        gnds.saveToFile( args.output )
