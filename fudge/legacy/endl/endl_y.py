# <<BEGIN-copyright>>
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
