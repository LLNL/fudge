#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import glob
import argparse
import pathlib

from fudge import enums as enumsModule
from fudge import map as mapModule
from fudge import GNDS_file as GNDS_fileModule
from PoPs import database as PoPsDatabaseModule
from fudge import suites as suitesModule
from fudge import reactionSuite as reactionSuiteModule
from fudge.covariances import covarianceSuite as covarianceSuiteModule

summaryDocStringFUDGE = '''Checks a GNDS map file and its contents for consistency.'''

description = """
This script checks a map file for consistency. In the following description a protare is considered 
to be in the map file if it is directly specified in the map file via a **protare** or **TNSL** node,
or indirectly specified in the map file via an **import** node, including any nesting of **import** 
nodes. The protare directories are the directories that contain the protares listed in the map file 
(while all protares can be in the same directory as the map file this is not recommended).

This script prints out the following information:

    -) A list of all sub-directory inside the protare directories.
    -) A list of all protares in the protare directories but not specified in the map file.
    -) A list of all protares specified in the map file but not found in the protare directories.
    -) A list of PoPs files found in the protare directories.
    -) A list of map files found in the protare directories.
    -) A list of unknown files found in the protare directories.
    -) A list of external files referenced in the protares but not found on the disk system.
    -) A list of non GNDS, HDF5 files found in HDF5 directories.
    -) A list of map files with missing ris files.

The header for each item above is only printed if the item contains entries.

If a protare (reactionSuite) file references a non-existance HDF5 values files, then the parser
for reading a reactionSuite will raise an exception that will exit this script.

Note, here a protare and reactionSuite are synonymous.
"""

parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('mapFile', type=pathlib.Path,                   help='The path to a GNDS map file.')
parser.add_argument('--skipRIS', action='store_true',               help='If present, missing *ris* information is not printed.')

args = parser.parse_args( )

def simplePrint( _list ) :

    for item in _list : print( '    %s' % item )

def printMissingExternalFile( _list ) :

    for item in _list :
        externalPath, sourceFile, externalRealPath = item.split( ',' )
        print( '    "%s" referenced by "%s".' % ( externalPath, sourceFile ) )

def printSet( prefix, _set, printFunction = simplePrint ) :

    _list = sorted( list( _set ) )
    if( len( _list ) > 0 ) :
        print( '%s:' % prefix )
        printFunction( _list )

#
# Phase 1: Make a list (i.e., set) of all protare directories and protares.
#
protareDirectories = set()
protaresInMap = set()
hdf5ExternalFiles = set()
hdf5Directories = set()
covarianceExternalFiles = set()
missingExternalFiles = set()
unsupportedExternalFileTypes = set()
missingRIS_files = set()

masterMap = None
protareFormatMisMatch = set()
mapFormatMisMatch = set()

def checkProtare( protareFileName, map, entry ) :

    global masterMap

    protaresInMap.add( protareFileName )
    protareDirectories.add( os.path.dirname( protareFileName ) )
    if( os.path.exists( protareFileName ) ) :
        protare = GNDS_fileModule.preview(protareFileName, haltParsingMoniker = suitesModule.ExternalFiles.moniker)

        if protare.format != masterMap.format:
            protareFormatMisMatch.add(protareFileName)

        interaction = entry.interaction
        if isinstance(entry, mapModule.TNSL) and entry.interaction is None:
            interaction = enumsModule.Interaction.TNSL
        if protare.evaluation != entry.evaluation or protare.interaction != interaction:
            print('''For map %s''' % map.fileName)
            print('''and for path %s''' % entry.path)
            if protare.evaluation != entry.evaluation:
                print('''    map's evaluation "%s" does not match protare's evaluation "%s".''' % (entry.evaluation, protare.evaluation))

            if protare.interaction != interaction:
                print('''    map's interaction "%s" does not match protare's interaction "%s".''' % (interaction, protare.interaction))
            print()
        for externalFile in protare.externalFiles :
            realpath = externalFile.realpath( )
            if( os.path.exists( realpath ) ) :
                name, dummy = GNDS_fileModule.type(realpath)
                if( name == GNDS_fileModule.HDF5_values) :  
                    hdf5ExternalFiles.add( realpath )
                    hdf5Directories.add( os.path.dirname( realpath ) )
                elif( name == covarianceSuiteModule.CovarianceSuite.moniker ) :
                    covarianceExternalFiles.add( realpath )
                else :
                    unsupportedExternalFileTypes.add( realpath )
            else :
                missingExternalFiles.add( ','.join( [ externalFile.path, protareFileName, realpath ] ) )

mapFileInfo = []
missingTNSL_standardTarget = set()

def checkMap(mapFileName, path, priorMaps=[]):

    global masterMap, mapFileInfo

    risFile = pathlib.Path(mapFileName).with_suffix('.ris')
    if not risFile.exists():
        missingRIS_files.add(str(mapFileName))

    try:
        map = mapModule.Map.read(mapFileName)
    except FileNotFoundError:
        print('ERROR: map file %s does not exist.' % mapFileName)
        for priorMap in priorMaps:
            print('    From map file %s' % priorMap)
        return
    except:
        raise

    if masterMap is None:
        masterMap = map
    else:
        if map.format != masterMap.format:
            mapFormatMisMatch.add(mapFileName)

    mapFileIndex = len(mapFileInfo)
    mapFileInfo.append([len(priorMaps), path, 0])

    for entry in map:
        if isinstance(entry, mapModule.Import):
            map2 = pathlib.Path(entry.fileName)
            recursiveMaps = ''
            recursionFound = False
            for priorMap in priorMaps:
                recursiveMaps += ' > %s' % priorMap
                if pathlib.Path(mapFileName).samefile(priorMap):
                    recursionFound = True
                    print('Recursive mapping found: %s %s' % (mapFileName, recursiveMaps))
                    break
            if not recursionFound:
                checkMap(entry.fileName, entry.path, [mapFileName] + priorMaps)
        else:
            mapFileInfo[mapFileIndex][2] += 1
            checkProtare(entry.fileName, map, entry)
            if isinstance(entry, mapModule.TNSL):
                if masterMap.find(entry.projectile, entry.standardTarget) is None:
                    missingTNSL_standardTarget.add(entry.standardTarget)

checkMap(str(args.mapFile), str(args.mapFile))

#
# Phase 2: Analyze all files in the protare directories.
#
dirsInDirectories = set( )
mapsInDirectories = set( )
protaresInDirectories = set( )
covariancesInDirectories = set( )
popsInDirectories = set( )
unknownsInDirectories = set( )

for protareDirectory in protareDirectories :
    files = glob.glob( os.path.join( protareDirectory, '*' ) )
    for file in files :
        if( os.path.isdir( file ) ) :
            dirsInDirectories.add( file )
        else :
            try:
                name, dummy = GNDS_fileModule.type(file)
                if( name == mapModule.Map.moniker ) :
                    mapsInDirectories.add( file )
                elif( name == reactionSuiteModule.ReactionSuite.moniker ) :
                    protaresInDirectories.add( file )
                elif( name == covarianceSuiteModule.CovarianceSuite.moniker ) :
                    covariancesInDirectories.add( file )
                elif( name == PoPsDatabaseModule.Database.moniker ) :
                    popsInDirectories.add( file )
                elif(name == GNDS_fileModule.HDF5_values) :
                    pass
                else :
                    unknownsInDirectories.add( file )
            except :
                unknownsInDirectories.add( file )

unknownsInHDF5_directories = set( )
for hdf5Directory in hdf5Directories :
    if( hdf5Directory not in protareDirectories ) :
        files = glob.glob( os.path.join( hdf5Directory, '*' ) )
        for file in files :
            try :
                if( os.path.isdir( file ) ) :
                    mapsInDirectories.add( file )
                else :
                    name, dummy = GNDS_fileModule.type(file)
                    if( name != GNDS_fileModule.HDF5_values) : unknownsInHDF5_directories.add(file)
            except :
                unknownsInHDF5_directories.add( file )

print('A total of %s protares found in map file (directly or indirectly). Broken down as:' % len(protaresInMap))
for level, path, numberOfProtares in mapFileInfo:
    print('%-60s: %s' % (((level + 1) * '    ' + path), numberOfProtares))

printSet('Directories in protare directories', dirsInDirectories)
printSet('Protares in protare directories but not in map file', protaresInDirectories.difference(protaresInMap))
printSet('Protares in map file but not in protare directories', protaresInMap.difference(protaresInDirectories))
printSet('PoPs files in protare directories', popsInDirectories)
printSet('Map files in protare directories', mapsInDirectories)
printSet('Unknown files in protare directories', unknownsInDirectories)
printSet('Missing external files', missingExternalFiles, printFunction = printMissingExternalFile)
printSet('Unknown files in HDF5 directories', unknownsInHDF5_directories)
printSet('Missing TNSL standard targets', missingTNSL_standardTarget)
if not args.skipRIS:
    printSet('Map files with missing RIS files', missingRIS_files)
printSet('Maps with format that does not match that of the master map file', mapFormatMisMatch)
printSet('Protares with format that does not match that of the master map file', protareFormatMisMatch)

unreferenceCovarianceFiles = set()
for covarianceInDirectories in covariancesInDirectories:
    if covarianceInDirectories not in covarianceExternalFiles:
        unreferenceCovarianceFiles.add(covarianceInDirectories)
printSet('Unreferenced covariance files', unreferenceCovarianceFiles)
