#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

summaryDocStringFUDGE = '''Prints the list of temperatures in a GNDS reactionSuite and labels for each processed style for each temperature.'''

description = '''For each GNDS reactionSuite file, prints the list of temperatures and labels for each processed style for each temperature.'''

from argparse import ArgumentParser

from pqu import PQU as PQUModule
from fudge import enums as enumsModule
from fudge import GNDS_file as GNDS_fileModule
from fudge import styles as stylesModule

temperatureUnitDefault = "K"

parser = ArgumentParser( description = description )
parser.add_argument( 'files', nargs = '*',                                                          help = 'List of GNDS files whose TNSL temperatures are printed.' )
parser.add_argument( '-u', '--temperatureUnit', action = 'store', default = temperatureUnitDefault, help = 'The unit to print temperatures. Default is "%s"' % temperatureUnitDefault )

args = parser.parse_args( )

TNSL_dataOrder = [ 'incoherent-inelastic', 'coherent-elastic', 'incoherent-elastic' ]

heatedProcessedStyles = sorted( [ [ style.sortOrderIndex, style.moniker ] for style in stylesModule.Styles.heatedProcessedStyles( ) ] )
heatedProcessedStyles = [ moniker for sortOrderIndex, moniker in heatedProcessedStyles ]

columnWidth = 24
labelFormat = '%%%ds' % columnWidth
temperatureFormat = '%.4g'
for file in args.files :
    try :
        reactionSuite = GNDS_fileModule.preview(file)
    except :
        print( '******** File is not a GNDS "reactionSuite" file: "%s"' % file )
        continue

    print( file )
    UnitAndTemperatures = [ 'K', [] ]

    if( reactionSuite.interaction == enumsModule.Interaction.TNSL ) :
        reactionSuite = GNDS_fileModule.read(file, lazyParsing=True)
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

    temperatureInfos = reactionSuite.styles.temperatures(unit=args.temperatureUnit)
    if len(temperatureInfos) > 0:
        print(labelFormat % ('temperature [%s]' % args.temperatureUnit), end='')
        for heatedProcessedStyle in heatedProcessedStyles : print(heatedProcessedStyle.center(columnWidth), end = '')
        print()
        print('   ', ((len( heatedProcessedStyles ) + 1) * columnWidth - 2) * '-')
        for temperatureInfo in temperatureInfos:
            print(labelFormat % (temperatureFormat % temperatureInfo.temperature ), end='')
            for moniker in heatedProcessedStyles: print(getattr(temperatureInfo, moniker).center(columnWidth), end='')
            print()
