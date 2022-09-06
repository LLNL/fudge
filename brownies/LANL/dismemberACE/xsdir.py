#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import sys
import argparse

description = """
Reads an MCNP xsdir file and prints requested information about its directory. An xsdir directory is a list of
"table"s. Each table starts with its name which is of the form '"ZA"."lib""cls"' (yes, the '"' are part of the 
symbolic names in this description) where "ZA", "lib" and "cls" are:

    +-------+---------------------------------------------------------------------+
    |     Z | The atomic number for an isotope.                                   |
    |     A | The atomic mass number for an isotope.                              |
    |     m | The meta-stable level for an isotope.                               |
    |  "ZA" | A + 1000 * Z + 1000000 * m or a special string for TNSL files.      |
    | "lib" | A number (1 to 3 digits) representing a library's identifier.       |
    | "cls" | To see the values for cls, run this script without any options.     |
    +-------+---------------------------------------------------------------------+

Examples include '100102.n' and '1027058.710nc'.

This script will only print out the legacy single character "cls" values. For example, it currently treats
'n' and 'nc' as the same. This is due in part because I could not find the complete list and have not a 
found documentation that describes the differences.

Two options ('--new' and '-q') do not print the requested information but exit early. They are mainly for debugging.
"""

parser = argparse.ArgumentParser( description = description, formatter_class = argparse.RawTextHelpFormatter )
parser.add_argument( '-v', '--verbose', action = 'count', default = 0,          help = 'Verbose mode.' )
parser.add_argument( 'xsdir', type = str,                                       help = 'xsdir file to read in.' )
parser.add_argument( '--ZA', type = str, default = None,                        help = 'The ZA of tables to find.' )
parser.add_argument( '--lib', type = int, default = -1,                         help = 'The library identifier of tables to find.' )
parser.add_argument( '--cls', type = str, default = None,                       help = 'The class of tables to find.' )
parser.add_argument( '--new', action = 'store_true',                            help = 'If present checks each ACE file found in the directory and prints its file name if it is an ACE 2.0 or greater file.' )
parser.add_argument( '-q', action = 'store_true',                               help = 'If present only reads the file and does not print anything. Mainly for debugging.' )
parser.add_argument( '--dismemberACE', action = 'store_true',                   help = 'If present, and the --cls and --lib are present then the dismemberACE.py command is written for each ZA.' )

priorFileName = ''
priorFileLines = []
badFiles = []
args = None

class Class :

    def __init__( self, classes, cls, label, supported ) :

        self.cls = cls
        self.label = label
        self.supported = supported
        self.libraries = {}

        classes[cls] = self

    def add( self, library, ZA, items ) :

        global priorFileName, priorFileLines

        ZA = ZA.strip( )

        if( library not in self.libraries ) : self.libraries[library] = {}
        if( ZA in self.libraries[library] ) :
            print( '    WARNING: ZA "%s" already in type "%s" and library "%s".' % ( ZA, self.type, library ) )
        else :
            self.libraries[library][ZA] = items

        if( args.new ) :
            fileName = items[2]
            if( fileName in badFiles ) : return
            if( fileName[0] != os.sep ) : fileName = os.path.join( os.path.dirname( args.xsdir ), fileName )
            if( priorFileName != fileName ) :
                priorFileName = fileName
                try :
                    with open( fileName ) as fIn : priorFileLines = fIn.readlines( )
                except :
                    badFiles.append( fileName )
                    print( 'Could not read file "%s".' % fileName )
                    return
            line = priorFileLines[int(items[5])-1]
            formatOrZA = line.split( )[0].strip( )
            if( formatOrZA[:4] == '2.0.' ) : print( '    %s' % fileName, items[5] )

if __name__ == '__main__':
    args = parser.parse_args( )

    classes = {}
    Class( classes, 'c', 'continuous-energy neutron', True )
    Class( classes, 'd', 'discrete-reaction neutron', False )
    Class( classes, 'y', 'dosimetry', False )
    Class( classes, 't', 'S(a,b) thermal', False )
    Class( classes, 'p', 'continuous-energy photoatomic', True )
    Class( classes, 'u', 'continuous-energy photonuclear', False )
    Class( classes, 'e', 'continuous-energy electron', False )
    Class( classes, 'm', 'multigroup neutron', False )
    Class( classes, 'g', 'multigroup photon', False )
    Class( classes, 'h', 'proton ?', False )
    Class( classes, 'o', 'deuteron ?', False )
    Class( classes, 'r', 'triton ?', False )
    Class( classes, 's', 'helium3 ?', False )
    Class( classes, 'a', 'alpha ?', False )

    classLabelMaximumLength = max( [ len( classes[cls].label ) for cls in classes ] )
    classLabelFormat = '%%-%ds' % ( classLabelMaximumLength + 2 )

    fIn = open( args.xsdir )
    lines = ''.join( fIn.readlines( ) )
    fIn.close( )

    lines = lines.split( '\ndire' )[1].split( '\n' )[1:-1]

    fileTypes = {}
    index = 0
    numberOfLines = len( lines )
    while( index < numberOfLines ) :
        lineIncrement = 1
        if( lines[index][0] == '#' ) :
            index += lineIncrement
            continue
        try :
            items = lines[index].split( )
            if( items[-1].strip( ) == '+' ) :
                items.pop( -1 )
                lineIncrement += 1
                items += lines[index+1].split( )
        except :
            print( index )
            print( lines[index] )
            raise

        try :
            ZA, suffix = items[0].split( '.' )
            cls = suffix
            library = ''
            while( cls[0].isdigit( ) ) :
                library += cls[0]
                cls = cls[1:]
            library = int( library )
            if( ( len( cls ) > 1 ) and ( cls != 'nc' ) ) : print( lines[index] )
            if( cls[-1] in classes ) :
                classes[cls[-1]].add( library, ZA, items )
            else :
                print( cls, lines[index] )
        except :
            print( lines[index] )
            raise

        index += lineIncrement

    if( args.q or args.new ) : sys.exit( 0 )

    if( ( args.cls is None ) and ( args.lib == -1 ) and ( args.ZA is None ) ) :
        for cls in classes :
            print( '    %2s %s: ' % ( cls, classLabelFormat % ( '[%s]' % classes[cls].label ) ) )
    else :
        if( args.cls is not None or True ) :
            clsIDs = [ args.cls ]
            if( args.cls is None ) : clsIDs = [ cls for cls in list( classes ) ]
            for clsID in clsIDs :
                cls = classes[clsID]
                printItems = []
                printItems.append( [ '    %2s %s: ' % ( clsID, classLabelFormat % ( '[%s]' % cls.label ) ), None ] )
                if( ( args.lib == -1 ) and ( args.ZA is None ) ) :
                    for lib in list( sorted( cls.libraries.keys( ) ) ) : print( '        %3s' % lib )
                else :
                    libs = [ args.lib ]
                    if( args.lib == -1 ) : libs = [ lib for lib in cls.libraries ]
                    for lib in libs :
                        if( lib in cls.libraries ) :
                            library = cls.libraries[lib]
                            if( args.ZA is None ) :
                                for ZA in library :
                                    if( args.dismemberACE ) :
                                        file = library[ZA][2]
                                        printItems.append(  [ '        ./dismemberACE.py --start %s %s   # %s' % ( library[ZA][5], file, ZA ), None ] )
                                    else :
                                        printItems.append( [ '        %-3s: %s' % ( ZA, library[ZA] ), None ] )
                            else :
                                if( args.ZA in library ) :
                                    file = library[args.ZA][2]
                                    if( file[0] != os.sep ) : file = os.path.join( os.path.dirname( args.xsdir), file )
                                    printItems.append( [ '        ./dismemberACE.py --start %s %s' % ( library[args.ZA][5], file ), library[args.ZA][0] ] )
                if( len( printItems ) > 1 ) :
                    print( printItems.pop( 0 )[0] )
                    if( printItems[0][1] is None ) :
                        for printItem, ZA in printItems : print( printItem )
                    else :
                        for printItem, ZA in printItems : print( "%-80s     # %s" % ( printItem, ZA ) )
