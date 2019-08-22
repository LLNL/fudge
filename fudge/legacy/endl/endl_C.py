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

"""
This module contains a set of routines that return information about the outgoing
particles for a given C-value as defined in ENDL.
"""
#  endl_CLabels is a structure of ( C, Lable, ( n, yo0, yo1, ... yon ) ).
#  Special cases for yo-tuple
#  ( )          Does not make sense to have an outgoing particle
#  ( -1, )      -1 = yoUnknown
#  ( -2, )      -2 = yo_eq_yi
#  ( -3, )      -3 = yoFission
#  ( -4, yo )   -4 = yoX

yoUnknown = -1
yo_eq_yi = -2
yoFission = -3
yoX = -4

endl_CLabels = (
    (   1,      "tot", ( ) ), \
    (   5,     "prod", ( ) ), \
    (   6,     "dest", ( ) ), \
    (   8,     "lacs", ( -2, ) ), \
    (   9,      "n+i", ( -2, ) ), \
    (  10,     "elas", ( -2, ) ), \
    (  11,        "n", ( 1, ) ), \
    (  12,       "2n", ( 1, 1 ) ), \
    (  13,       "3n", ( 1, 1, 1 ) ), \
    (  14,       "4n", ( 1, 1, 1, 1 ) ), \
    (  15,        "f", ( -3, ) ), \
    (  16,      "3np", ( 1, 1, 1, 1001 ) ), \
    (  17,      "n2p", ( 1, 1001, 1001 ) ), \
    (  18,       "2p", ( 1001, 1001 ) ), \
    (  19,       "pd", ( 1001, 1002 ) ), \
    (  20,      "n p", ( 1, 1001 ) ), \
    (  21,      "p n", ( 1001, 1 ) ), \
    (  22,      "n d", ( 1, 1002 ) ), \
    (  23,    "n d a", ( 1, 1002, 2004 ) ), \
    (  24,      "n t", ( 1, 1003 ) ), \
    (  25,    "n He3", ( 1, 2003 ) ), \
    (  26,      "n a", ( 1, 2004 ) ), \
    (  27,     "n 2a", ( 1, 2004, 2004 ) ), \
    (  28,    "n t a", ( 1, 1003, 2004 ) ), \
    (  29,     "2n p", ( 1, 1, 1001 ) ), \
    (  30,    "g n a", ( 7, 1, 2004 ) ), \
    (  31,   "2n p a", ( 1, 1, 1001, 2004 ) ), \
    (  32,     "2n d", ( 1, 1, 1002 ) ), \
    (  33,     "2n a", ( 1, 1, 2004 ) ), \
    (  34,    "n p a", ( 1, 1001, 2004 ) ), \
    (  35,      "d n", ( 1002, 1 ) ), \
    (  36,     "n 3a", ( 1, 2004, 2004, 2004 ) ), \
    (  37,       "2a", ( 2004, 2004 ) ), \
    (  38,    "He3 a", ( 2003, 2004 ) ), \
    (  39,      "p t", ( 1001, 1003 ) ), \
    (  40,        "p", ( 1001, ) ), \
    (  41,        "d", ( 1002, ) ), \
    (  42,        "t", ( 1003, ) ), \
    (  43,      "t a", ( 1003, 2004 ) ), \
    (  44,      "He3", ( 2003, ) ), \
    (  45,        "a", ( 2004, ) ), \
    (  46,        "g", ( 7, ) ), \
    (  47,      "d a", ( 1002, 2004 ) ), \
    (  48,      "p a", ( 1001, 2004 ) ), \
    (  49,     "2p a", ( 1001, 1001, 2004 ) ), \
    (  50,      "X p", ( -4, 1001 ) ), \
    (  51,      "X d", ( -4, 1002 ) ), \
    (  52,      "X t", ( -4, 1003 ) ), \
    (  53,    "X He3", ( -4, 2003 ) ), \
    (  54,      "X a", ( -4, 2004 ) ), \
    (  55,      "X g", ( -4, 7 ) ), \
    (  56,      "X n", ( -4, 1 ) ), \
    (  57,      "X e", ( -4, 8 ) ), \
    (  65,      "act", ( -1, ) ), \
    (  66,    "yield", ( -1, ) ), \
    (  70,     "totp", ( -1, ) ), \
    (  71,      "coh", ( -2, ) ), \
    (  72,    "incoh", ( -2, ) ), \
    (  73,    "photo", ( 9, ) ), \
    (  74,     "pair", ( 8, 9 ) ), \
    (  75,  "triplet", ( 8, 9, 9 ) ), \
    (  78,       "ic", ( -1, ) ), \
    (  81,      "ion", ( 9, 9 ) ), \
    (  82,     "brem", ( 9, 7 ) ), \
    (  83,    "excit", ( -2, ) ), \
    (  84,     "coll", ( -1, ) ), \
    (  91,    "shell", ( -1, ) ), \
    (  92,    "trans", ( -1, ) ), \
    (  93,    "whole", ( -1, ) ) )

def endl_nCMax( ) :
    """Returns the largest C-value for which information exist."""

    return( endl_CLabels[-1][0] )

def endl_CLabel( C ) :
    """Returns a short string (i.e., mnemonic) describing the outgoing particles - called yos in ENDL - 
or reaction type for the specified C-value. If C-value is not defined then None is returned.
"""

    for i in endl_CLabels :
        if ( i[0] == C ) : return i[1]
    return None

def endl_CLabelMaxWidth( ) :
    """Returns the maximum length of a string returned by endl_CLabel."""

    w = 0
    for i in endl_CLabels :
        if ( w < len( i[1] ) ) : w = len( i[1] )
    return w

def endl_C_yoInfo( C ) :
    """Returns a tuple of integers describing the outgoing particles - called yos in ENDL -
for the specified C-value, or a negative number (see "Special cases ..." below) if outgoing particle
types/numbers are not defined. Each particle is defined by its ZA = 1000 * Z + A value, with the
exceptions of gammas, positrons and electrons whose integer designations are 7, 8 and 9 respectively. For
example, neutron, Helium-3 and Pu-239 are designated as 1, 2003 and 94239 respectively.
If C-value is not defined then None is returned.

Special cases for first integer of returned tuple::
  -1  Unknown outgoing particles.
  -2  Outgoing particle is the same as the incident particle.
  -3  Unknown number of fission neutrons.
  -4  Unknown number of particles of type of the second integer of the tuple.
"""

    for i in endl_CLabels :
        if ( i[0] == C ) : return ( i[2] )
    return None

def endl_C_CList( ) :
    """Returns a list of C-values."""

    CList = []
    for C, symbols, yos in endl_CLabels : CList.append( C )
    return CList

yoDict = { 1:1, 1001:2, 1002:3, 1003:4, 2003:5, 2004:6 }

def endl_C_yosMultiplicity( C, yi ) :
    """Returns the list [ nn, np, nd, nt, nHe3, na, ng ] where nn is the multiplicity
    for neutrons, np is the multiplicity for protons, etc. If C-value is not defined
    then None is returned. If the reaction is fission then -1 is returned for neutron
    multiplity (i.e., C = 15 returns  [ -1, 0, 0, 0, 0, 0, 0 ] ). If multiplity for a
    particle is indeterminate then -2 is returned for that multiplity.  If multiplicities
    are undefined then nYos * [ None ] is returned. yi is only needed for elastic (C = 10)
    scattering."""

    nYos = 7
    yosMultiplicities = nYos * [ 0 ]
    yos = endl_C_yoInfo( C )
    if( yos == None ) :
        return( None )
    elif( yos[0] == yoUnknown ) :
        return( nYos * [ None ] )
    elif( yos[0] == yo_eq_yi ) :
        yo = yi
        if( yo in yoDict ) : yo = yoDict[yo]
        if( 1 <= yo <= 7 ) :
            yosMultiplicities[yo-1] = 1
        else :
            raise Exception( 'Error from endl_C.endl_C_yosMultiplicity: unsupported yi = %s' % `yi` )
    elif( yos[0] == yoFission ) :
        yosMultiplicities[0] = -1
    elif( yos[0] == yoX ) :
        if( yos[1] in yoDict ) :
            yosMultiplicities[yoDict[yos[1]]-1] = -2
        else :
            return( nYos * [ 0 ] )
    else :
        for yo in yos :
            if( yo not in [ 7, 8, 9 ] ) : yosMultiplicities[yoDict[yo]-1] += 1

    return( yosMultiplicities )
