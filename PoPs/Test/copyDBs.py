#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from PoPs import database as databaseModule

fOut = open( 'copyDBs.py.copy.out', 'w' )

def _copy( fileName ) :

    fIn = open( fileName )
    fOut.writelines( fIn.readlines( ) )
    fIn.close( )

    database1 = databaseModule.database.readFile( fileName )
    database2 = database1.copy( )
    print( database2.toXML( ) )

_copy( 'Answers/database.py.out' )
_copy( 'Answers/database3.py.out' )
_copy( 'Answers/database4.py.out' )

fOut.close( )
