#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

description = """Prints the file name, the temperature unit and the list of evaluated temperatures for a thermal neutron scattering law (TNSL) GNDS file 
    as string following by a python list of floats.  If the GNDS file is not a TNSL file, an empty list of temperatrures is printed."""

from argparse import ArgumentParser

from pqu import PQU as PQUModule
from fudge import reactionSuite as reactionSuiteModule
from fudge import styles as stylesModule
from LUPY import GNDSType as GNDSTypeModule

temperatureUnitDefault = "K"

parser = ArgumentParser( description = description )
parser.add_argument( 'files', nargs = '*',                                                          help = 'List of GNDS files whose TNSL temperatures are printed.' )
parser.add_argument( '-u', '--temperatureUnit', action = 'store', default = temperatureUnitDefault, help = 'The unit to print temperatures. Default is "%s"' % temperatureUnitDefault )

args = parser.parse_args( )

TNSL_dataOrder = [ 'incoherent-inelastic', 'coherent-elastic', 'incoherent-elastic' ]

heatedProcessedStyles = sorted( [ [ style.sortOrderIndex, style.moniker ] for style in stylesModule.styles.heatedProcessedStyles( ) ] )
heatedProcessedStyles = [ moniker for sortOrderIndex, moniker in heatedProcessedStyles ]

columnWidth = 24
labelFormat = '%%%ds' % columnWidth
temperatureFormat = '%.4g'
for file in args.files :
    try :
        reactionSuite = GNDSTypeModule.preview( file )
    except :
        print( '******** File is not a GNDS "reactionSuite" file: "%s"' % file )
        continue

    print( file )
    UnitAndTemperatures = [ 'K', [] ]

    if( reactionSuite.interaction == reactionSuiteModule.Interaction.TNSL ) :
        allTemperatures = reactionSuite.thermalNeutronScatteringLawTemperatures( )
        for name in TNSL_dataOrder :
            if( name in allTemperatures ) :
                UnitAndTemperatures = allTemperatures[name]
                break

        print( '  TNSL evaluated temperatures [%s] are:' % args.temperatureUnit, end = '' )
        unit, temperatures = UnitAndTemperatures
        sepChar = ' '
        for temperature in temperatures :
            print( '%s%s' % ( sepChar, float( "%.5g" % PQUModule.PQU( temperature, unit ).getValueAs( args.temperatureUnit ) ) ), end = '' )
            sepChar = ', '
        print( )

    else :
        preProcessingChains = reactionSuite.styles.preProcessingChains( ends = True )
        for preProcessingChain in preProcessingChains :
            temperature = preProcessingChain[0].temperature.getValueAs( args.temperatureUnit )
            print( '  temperature %.5s %s:' % ( temperature, args.temperatureUnit ), end = '' )
            sepChar = ' '
            for style in reversed( preProcessingChain ) :
                print( '%s%s' % ( sepChar, style.label ), end = '' )
                sepChar = ', '
            print( )

    temperatureInfos = reactionSuite.styles.temperatures( unit = args.temperatureUnit )
    if( len( temperatureInfos ) > 0 ) :
        print( labelFormat % ( 'temperature [%s]' % args.temperatureUnit ), end = '' )
        for heatedProcessedStyle in heatedProcessedStyles : print( heatedProcessedStyle.center( columnWidth ), end = '' )
        print( )
        print( '   ', ( ( len( heatedProcessedStyles ) + 1 ) * columnWidth - 2 ) * '-' )
        for temperature, labels in temperatureInfos  :
            print( labelFormat % ( temperatureFormat % temperature ), end = '' )
            for moniker in heatedProcessedStyles : print( labels[moniker].center( columnWidth ), end = '' )
            print( )
