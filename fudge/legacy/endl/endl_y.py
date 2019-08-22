# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# <<END-copyright>>

"""
This module contains a set of routines that return the symbol or name for an ENDL 
incident, yi, or outgoing, yo, particle designator.
"""

endl_yLabels = (
    (  1, "n", "neutron" ), \
    (  2, "p", "proton" ), \
    (  3, "d", "deuteron" ), \
    (  4, "t", "triton" ), \
    (  5, "He3", "helium3" ), \
    (  6, "He", "alpha" ), \
    (  7, "g", "gamma" ), \
    (  8, "e+", "positron" ), \
    (  9, "e-", "electron" ), \
    ( 10, "EC", "electron capture" ), \
    ( 11, "n as residual", "n as residual" ), \
    ( 12, "p as residual", "p as residual" ), \
    ( 13, "d as residual", "d as residual" ), \
    ( 14, "t as residual", "t as residual" ), \
    ( 15, "He3 as residual", "He3 as residual" ), \
    ( 16, "He as residual", "He as residual" ), \
    ( 18, "positron as residual", "positron as residual" ), \
    ( 19, "electron as residual", "electron as residual" ) )

def endl_yLabel( y ) :
    """Returns the symbol for the specified y or 'None' if y is out-of-bounds."""

    for i in endl_yLabels :
        if ( i[0] == y ) : return i[1]
    return None

def endl_yLongLabel( y ) :
    """Returns the label (i.e., name) for the specified y or 'None' if y is out-of-bounds."""

    for i in endl_yLabels :
        if ( i[0] == y ) : return i[2]
    return None

def endl_yLabelMaxWidth( ) :
    """Returns the maximum length of a string returned by endl_yLabel."""

    w = 0
    for i in endl_yLabels :
        if ( w < len( i[1] ) ) : w = len( i[1] )
    return w

def endl_yLongLabelMaxWidth( ) :
    """Returns the maximum length of a string returned by endl_yLongLabel."""

    w = 0
    for i in endl_yLabels :
        if ( w < len( i[2] ) ) : w = len( i[2] )
    return w
