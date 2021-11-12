# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from PoPs.groups import misc as chemicalElementMiscPoPsModule

def getZandAFromName( name ) :

    if( name == 'n' ) : return( ( 0, 1 ) )
    for i, l in enumerate( name ) :
        if( l.isdigit( ) ) : break
    symbol, AStr = name[:i], name[i:]
    if '_' in AStr: AStr=AStr.split('_')[0]
    Z = chemicalElementMiscPoPsModule.ZFromSymbol[symbol]
    for i, l in enumerate( AStr ) :
        if( not( l.isdigit( ) ) ) : break
    A = int( AStr[:i+1] )
    return( Z, A )

def getZAFromName( name ):

    Z, A = getZandAFromName( name )
    return( 1000 * Z + int( A ) )

def getZ_A_suffix_andZAFromName( name ) :
    """
    Returns the tuple (Z, A, suffix, ZA) for an gnds isotope name (e.g., gnds name = 'Am242_m1'
    returns ( 95, 242, 'm1', 95242 ).
    """

    if( name == 'n' ) : return( 0, 1, '', 1 )
    if( name == 'gamma' ) : return( 0, 0, '', 0 )
    if( name == 'photon' ) : return( 0, 0, '', 0 )
    if( name[:18] == 'FissionProductENDL' ) :
        ZA = int( name[18:] )
        Z = ZA // 1000
        A = 1000 * Z - ZA
        return( Z, A, '', ZA )
    if( '__' in name ) : raise Exception ( "Name = %s" % name )
    naturalSuffix = ''
    if( '_' in name ) :         # Isotope names can have level designator (e.g., 'O16_e3') and naturals are of the form 'S_natural' or 'S_natural_l'
        s = name.split( '_' )   # where S is element's symbol and l is level designator (e.g., 'Xe_natural' or 'Xe_natural_c').
        sZA, suffix = s[:2]
        if( len( s ) > 2 ) :
            if( ( len( s ) > 3 ) or ( suffix != 'natural' ) ) : raise Exception( 'Invalid name for endl ZA particle = %s' % name )
            naturalSuffix = s[2]
    else :
        sZA = name
        suffix = ''
    for i, c in enumerate( sZA ) :
        if( c.isdigit( ) ) : break
    if( not c.isdigit( ) ) : i += 1
    sZ, sA = sZA[:i], sZA[i:]
    Z = chemicalElementMiscPoPsModule.ZFromSymbol[sZ]
    if( sA == '' ) :
        if( suffix == 'natural' ) : return( Z, 0, naturalSuffix, 1000 * Z )
        if( suffix == '' ) : return( Z, 0, '', 1000 * Z )
        raise Exception( 'No A for particle named %s' % name )
    elif( suffix == 'natural' ) :
        raise Exception( 'Natural element also has A defined for particle named %s' % name )
    else :
        try :
            A = int( sA )
        except :
            raise Exception( 'Could not convert A to an integer for particle named %s' % name )
    ZA = 1000 * Z + A
    return( Z, A, suffix, ZA )

# Crack the isotope name to get the A & symbol
def elementAFromName( name ) :

    sym = ''
    A = ''
    if( '_' in name ) :
        m = name.split('_')[1]
    else :
        m = None
    for c in name.split( '_' )[0] :
        if( c.isalpha( ) ) :
            sym += c
        else :
            A += c
    if( sym == 'n' ) : return( sym, 1, None )
    if( sym == 'g' ) : return( sym, 0, None )
    if( m == 'natural' ) : A = '0'
    return( sym, A, m )
