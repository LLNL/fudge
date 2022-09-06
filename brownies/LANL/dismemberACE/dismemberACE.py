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

outputDefault = 'dismemberACE.out'
ACE_format_1_0 = '1.0'
ACE_format_2_0 = '2.0'

from brownies.LANL.dismemberACE  import dismemberACE_misc as dismemberACE_miscModule
from brownies.LANL.dismemberACE import dismemberACE_c as dismemberACE_c_Module
from brownies.LANL.dismemberACE import dismemberACE_p as dismemberACE_p_Module
from brownies.LANL.dismemberACE import dismemberACE_t as dismemberACE_t_Module
from brownies.LANL.dismemberACE import dismemberACE_e as dismemberACE_e_Module

description = """
This module dismembers the MT data in an ACE file into separate sub-directories, one for each MT.
"""

parser = argparse.ArgumentParser( description = description )
parser.add_argument( '-a', '--addresses', action = 'store_true',                help = 'Add addresses (locators) to output (for debugging).' )
parser.add_argument( '-f', '--file', action = 'store_true',                     help = 'If present, the ACE file name is printed at startup.' )
parser.add_argument( '-k', '--keep', action = 'store_true',                     help = 'The contents of the output directory is removed unless this option is present.' )
parser.add_argument( '-q', '--quiet', action = 'store_true',                    help = 'Quiet mode. The citation data are not printed if all data have been cited.' )
parser.add_argument( '-O', '--Official', action = 'store_true',                 help = 'Only print information about official cites.' )
parser.add_argument( '-s', '--start', action = 'store', type = int, default = 1,    help = 'The starting offset of the data in the ace file. Default is 1.' )
parser.add_argument( '-v', '--verbose', action = 'count', default = 0,          help = 'Verbose mode.' )
parser.add_argument( 'ace', type = str,                                         help = 'ACE file name to dismember.' )
parser.add_argument( 'output', type = str, default = outputDefault,             help = 'Name of directory to put the MT sub-directories.' )

args = parser.parse_args( )

__doc__ = description

def dismemberACE( fileName ) :

    if( args.file ) : print( fileName )

    fIn = open( fileName )
    lines = fIn.readlines( )[args.start-1:]
    fIn.close( )

    print( lines[0] )
    formatOrZA = lines[0].split( )[0].strip( )
    if( formatOrZA[:len(ACE_format_2_0)] == ACE_format_2_0 ) : 
        ACE_format = ACE_format_2_0
        SZAID = lines[0].split( )[1].strip( )
        offset = int( lines[1].split( )[-1].strip( ) )
    else :
        ACE_format = ACE_format_1_0
        SZAID = lines[0].split( )[0].strip( )
        offset = 0
    ZA, suffix = SZAID.split( '.' )
    cls = suffix
    library = ''
    while( cls[0].isdigit( ) ) :
        library += cls[0]
        cls = cls[1:]
    library = int( library )

    fOut = dismemberACE_miscModule.openFile2( args, 0, 'header' )
    if( offset > 0 ) :
        fOut.write( ''.join( lines[:offset] ) )
        lines = lines[offset:]

    fOut.write( ''.join( lines[:12] ) )
    fOut.write( '\n' )

    NXS  = dismemberACE_miscModule.get8Integers( lines[6] )
    NXS += dismemberACE_miscModule.get8Integers( lines[7] )
    NXS = dismemberACE_miscModule.ListBase1( NXS )
    for i1 in range( 1, 16 ) : fOut.write( 'NXS[%2d] = %9d\n' % ( i1, NXS[i1] ) )

    fOut.write( '\n' )
    JXS  = dismemberACE_miscModule.get8Integers( lines[8] )
    JXS += dismemberACE_miscModule.get8Integers( lines[9] )
    JXS += dismemberACE_miscModule.get8Integers( lines[10] )
    JXS += dismemberACE_miscModule.get8Integers( lines[11] )
    JXS = dismemberACE_miscModule.ListBase1( JXS )
    for i1 in range( 1, 32 ) : fOut.write( 'JXS[%2d] = %9d\n' % ( i1, JXS[i1] ) )

    fOut.close( )

    XSS = []
    count = NXS[1]
    for i1, line in enumerate( lines[12:] ) :
        while( count == 0 ) : break
        count, lineXSS = dismemberACE_miscModule.get4Floats( line, count, i1 )
        XSS += lineXSS

    if( len( XSS ) != NXS[1] ) : print( 'WARNING: Number of stated XSS value (%d) differ than number found (%s)' % ( NXS[1], len( XSS ) ) )

    dismemberACE_miscModule.initialize( args, XSS )

    fOut = dismemberACE_miscModule.openFile( 0, 'XSS' )
    count = 0
    for XX in XSS :
        if( ( abs( XX ) < NXS[1] ) and ( int( XX ) == XX ) ) :
            count += 1
            fOut.write( '%9d\n' % XX )
        else :
            count += 1
            fOut.write( '%20.12e\n' % XX )
    fOut.close( )
    XSS = dismemberACE_miscModule.ListBase1( XSS )
    dismemberACE_miscModule.cites = len( XSS ) * [ 0 ]

    if( cls in [ 'c', 'nc' ] ) :
        dismemberACE_c_Module.dismemberACE( args, NXS, JXS, XSS )
    elif( cls in [ 'p' ] ) :
        dismemberACE_p_Module.dismemberACE( args, NXS, JXS, XSS )
    elif( cls in [ 't' ] ) :
        dismemberACE_t_Module.dismemberACE( args, NXS, JXS, XSS )
    elif( cls in [ 'e' ] ) :
        dismemberACE_e_Module.dismemberACE( args, NXS, JXS, XSS )
    else :
        print( 'Unsupported library class = "%s".' % cls )

    fOut = dismemberACE_miscModule.openFile( 0, 'cites' )
    uncited = 0
    multiplyCited = 0
    officialEnd = JXS[22]
    pastEndUncited = 0
    for i1, cite in enumerate( dismemberACE_miscModule.cites ) :
        if( cite == 0 ) : uncited += 1
        if( cite > 1 ) : multiplyCited += 1
        if( ( i1 >= officialEnd ) and ( cite == 0 ) ) : pastEndUncited += 1
        fOut.write( '%d\n' % cite )
    fOut.close( )

    fOut = dismemberACE_miscModule.openFile( 0, 'cites.compact' )
    priorCite = -1
    priorIndex = 0
    for index, cite in enumerate(dismemberACE_miscModule.cites):
        if cite != priorCite:
            fOut.write( '%d %8d %8d\n' % (cite, index + 1, index - priorIndex) )
            priorCite = cite
            priorIndex = index
    fOut.close( )

    if( ( uncited > 0 ) or not( args.quiet ) ) :
        if( not( args.Official and ( ( uncited - pastEndUncited ) == 0 ) ) ) :
            if( args.verbose > 0 ) : print( "    Official uncited data elements = %d, total uncited = %d, multiply cited = %d" % 
                ( ( uncited - pastEndUncited ), uncited, multiplyCited ) )

if( __name__ == '__main__' ) :

    if( not( args.keep ) ) : os.system( 'rm -rf %s/*' % args.output )
    dismemberACE( args.ace )
