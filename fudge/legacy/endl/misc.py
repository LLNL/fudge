# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
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
        Z = ZA / 1000
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
