#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import argparse

from brownies.legacy.endl import bdfls as bdflsModule

addMultigroupExec = os.path.join( os.path.dirname( os.path.dirname( os.path.dirname( os.path.dirname( os.path.realpath( __file__ ) ) ) ) ), "bin", "addMultigroup.py" )

parser = argparse.ArgumentParser( description = "Converts one or more bdfls multi-groups in a bdfls file into the new GNDS multi-group structure and outputs results to a file." )

parser.add_argument( 'bdfls',                                           help = "The name of the bdfls file to extract multi-group data from." )
parser.add_argument( "output",                                          help = "The name of the outputted multi-group file." )
parser.add_argument( "input", nargs = "?", default = None,              help = "The file to read existing multi-group data from." )
parser.add_argument( '-g', '--gids', action = 'append',                 help = "Append the gid to the bdfls multi-group id to convert. If absent, all are converted (e.g., -g 4)." )
parser.add_argument( "--override", action = "store_true",               help = "If label exists and option present, replace multi-group; otherwise, execute a raise." )

if( __name__ == '__main__' ) :

    args = parser.parse_args( )

    bdfls = bdflsModule.bdfls( template = args.bdfls )

    if( args.gids is None ) :
        gids = [ group.id for group in bdfls.g ]
    else :
        gids = [ int( gid ) for gid in args.gids ]

    override = ""
    if( args.override ) : override = " --override"

    input = ''
    if( args.input is not None ) : input = args.input

    for group in bdfls.g :
        if( group.id in gids ) :
            boundaries = ' '.join( [ str( value ) for value in group.gb ] )
            os.system( """%s %s --values %s "%s" %s %s""" % ( addMultigroupExec, override, group.name, boundaries, args.output, input ) )
            input = ''
