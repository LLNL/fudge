#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import argparse

from brownies.legacy.endl import bdfls as bdflsModule
from xData import XYs1d as XYs1dModule

addFluxExec = os.path.join( os.path.dirname( os.path.dirname( os.path.dirname( os.path.dirname( os.path.realpath( __file__ ) ) ) ) ), "bin", "addFlux.py" )

parser = argparse.ArgumentParser( description = "Converts one or more bdfls fluxesin a bdfls file into the new GNDS flux structure and outputs results to a file." )

parser.add_argument( 'bdfls',                                           help = "The name of the bdfls file to extract flux data from." )
parser.add_argument( "output",                                          help = "The name of the outputted flux file." )
parser.add_argument( "input", nargs = "?", default = None,              help = "The file to read existing flux data from." )
parser.add_argument( '-f', '--fids', action = 'append',                 help = "Append the fid to the bdfls flux id to convert. If absent, all are converted (e.g., -f 1)." )
parser.add_argument( "--override", action = "store_true",               help = "If label exists and option present, replace flux with new one; otherwise, execute a raise." )

if( __name__ == '__main__' ) :

    args = parser.parse_args( )

    bdfls = bdflsModule.bdfls( template = args.bdfls )

    if( args.fids is None ) :
        fids = [ flux.id for flux in bdfls.f ]
    else :
        fids = [ int( fid ) for fid in args.fids ]

    override = ""
    if( args.override ) : override = " --override"

    input = ''
    if( args.input is not None ) : input = args.input

    for bdflsFlux in bdfls.f :
        if( bdflsFlux.id in fids ) :
            orders = []
            grid = XYs1dModule.XYs1d( )                             # Only used to get a comment energy grid.
            for order in bdflsFlux.EF_l :
                flux = XYs1dModule.XYs1d( order )
                orders.append( flux )
                grid += flux
            flux = [ ]
            for energy, dummy in grid :
                energyLegendreCoefficients = '%s' % energy
                for order in orders : energyLegendreCoefficients += ' %s' % order.evaluate( energy )
                flux.append( energyLegendreCoefficients )
            os.system( """%s %s --values %s "%s" %s %s""" % ( addFluxExec, override, bdflsFlux.name, "; ".join( flux ), args.output, input ) )
            input = ''
