# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains a set of routines that return the symbol or name for an element.
"""

from PoPs.chemicalElements import misc as chemicalElementMiscModule

def endl_nZs( ) :
    """Returns the number of elements, starting at neutron with Z = 0, for which data is present.
    The largest element for which information can be obtained has Z = endl_nZs( ) - 1."""

    return( len( chemicalElementMiscModule.chemicalElementZSymbolNames ) )

def endl_ZSymbol( Z ) :
    """Returns the symbol for the specified Z or 'None' if Z is out-of-bounds."""
    try :
        return( chemicalElementMiscModule.symbolFromZ[Z] )
    except :
        return None

def endl_ZLabel( Z ) :
    """Returns the label (i.e., name) for the specified Z or 'None' if Z is out-of-bounds."""

    try :
        return( chemicalElementMiscModule.nameFromZ[Z] )
    except :
        return None

def endl_ZSymbolToZ( symbol ) :
    """Returns the Z for the specified symbol or 'None' if no match for symbol."""

    try :
        return( chemicalElementMiscModule.ZFromSymbol[symbol] )
    except :
        return None

def endl_ZLabelToZ( label ) :
    """Returns the Z for the specified label or 'None' if no match for label."""

    try :
        return( chemicalElementMiscModule.ZFromName[label] )
    except :
        return None
