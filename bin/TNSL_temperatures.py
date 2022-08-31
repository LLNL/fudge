#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

description = """Prints the file name, the temperature unit and the list of evaluated temperatures for a thermal neutron scattering law (TNSL) GNDS file 
    as string following by a python list of floats.  If the GNDS file is not a TNSL file, an empty list of temperatrures is printed."""

from argparse import ArgumentParser

from pqu import PQU as PQUModule
from fudge import reactionSuite as reactionSuiteModule

temperatureUnitDefault = "K"

parser = ArgumentParser( description = description )
parser.add_argument( 'files', nargs = '*',                                                          help = 'List of GNDS files whose TNSL temperatures are printed.' )
parser.add_argument( '-a', '--all', action = 'store_true',                                          help = 'If present, temperatures for all double differential data are printed. Otherwise, just the "incoherent-inelastic" are generally printed.' )
parser.add_argument( '-u', '--temperatureUnit', action = 'store', default = temperatureUnitDefault, help = 'The unit to print temperatures. Default is "%s"' % temperatureUnitDefault )

args = parser.parse_args( )

def printTemperatures( temperatures ) :

    unit,  temperatures = temperatures
    if( unit != args.temperatureUnit ) :
        temperatures = [ float( "%.5g" % PQUModule.PQU( temperature, unit ).getValueAs( args.temperatureUnit ) ) for temperature in temperatures ]
        unit = args.temperatureUnit
    print( unit, temperatures )

dataOrder = [ 'incoherent-inelastic', 'coherent-elastic', 'incoherent-elastic' ]

for file in args.files :
    print( file, end = ' ' )
    temperatures = [ 'K', [] ]

    reactionSuite = reactionSuiteModule.ReactionSuite.readXML_file(file)
    allTemperatures = reactionSuite.thermalNeutronScatteringLawTemperatures( )
    if( args.all ) :
        print( )
        for name in dataOrder :
            if( name in allTemperatures ) : print( '    %s %s' % ( name, allTemperatures[name] ) )
    else :
        for name in dataOrder :
            if( name in allTemperatures ) :
                temperatures = allTemperatures[name]
                break

        printTemperatures( temperatures )
