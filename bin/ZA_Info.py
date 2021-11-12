#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse

from PoPs.groups import misc as PoPs_groupsMiscModule

description = """For each argument entered, which must be an isotope name specified by either its ZA (1000 * Z + A) or its PoPs id,
this script prints that argument, its ZA and PoPs id."""

parser = argparse.ArgumentParser( description = description )

parser.add_argument( 'names', nargs = '*',       help = 'The list of isotope names specified by either ZA and/or PoPs id.' )

args = parser.parse_args( )

width1 = 0
width3 = 0
data = []
for name in args.names :
    if( name[0] in '123456789' ) :
        ZA = int( name )
        popsId = PoPs_groupsMiscModule.idFromZA( ZA )
    else :
        popsId = name
        ZA = PoPs_groupsMiscModule.ZAInfo_fromString( popsId )[2]

    width1 = max( width1, len( name ) )
    width3 = max( width3, len( popsId ) )
    data.append( [ name, ZA, popsId ] )

format = '  %%-%ds |     ZA | %%-%ds' % ( width1, width3 )
header = format % ( 'name','PoPs id' )
print( header )
print( ' %s' % ( len( header ) * '-' ) )
format = '  %%-%ds | %%6d | %%-%ds' % ( width1, width3 )
for name, ZA, popsId in data :
    print( format % ( name, ZA, popsId ) )
