# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains a set of routines that return information about ENDL I-values.
"""

endl_ILabels = (
    (  0,  2, "cross section" ), \
    (  1,  3, "angular dist." ), \
    (  3,  4, "E-angle dist." ), \
    (  4,  4, "l-order E-angle dist." ), \
    (  7,  2, "n_bar / fission" ), \
    (  8,  0, "histogram E dist." ), \
    (  9,  2, "multiplicity" ), \
    ( 10,  2, "secondary ave. E" ), \
    ( 11,  2, "residual ave. E" ), \
    ( 12,  2, "Q(E)" ), \
    ( 13,  2, "momentum" ), \
    ( 20,  4, "URR(E,T)" ), \
    ( 21,  3, "P(E|E_out)" ), \
    ( 22,  3, "angular dist. in (1-mu)" ), \
    ( 80,  0, "Maxwell reaction rate" ), \
    ( 81,  0, "Doppler x-sec." ), \
    ( 84,  0, "Maxwell E dist." ), \
    ( 89,  0, "Maxwell multiplicity" ), \
    ( 90,  0, "Maxwell secondary ave. E" ), \
    ( 91,  0, "Maxwell residual ave. E" ), \
    ( 92,  0, "Maxwell E of reaction part." ) )

def endl_ILabel( I ) :
    """Returns a mnemonic string for the specified I-value. If I-value is not defined than 'None' is returned."""

    for i in endl_ILabels :
        if ( i[0] == I ) : return i[2]
    return None

def endl_ILabelMaxWidth( ) :
    """Returns the maximum length of a string returned by endl_ILabel."""

    w = 0
    for i in endl_ILabels :
        if ( w < len( i[2] ) ) : w = len( i[2] )
    return w

def endl_IColumns( I ) :
    """Returns the number of columns of data for the specified I-value. If I-value is not defined than 'None' is returned."""

    for i in endl_ILabels :
        if ( i[0] == I ) : return i[1]
    return None
