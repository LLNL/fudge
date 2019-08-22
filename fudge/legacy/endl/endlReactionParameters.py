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

class endlReactionParameters :
    """This class contains the parameters C, S, X1, X2, X3, X4, Q for a reaction."""

    def __init__( self, C, S, X1, X2, X3, X4, Q ) :
        """Creates a new endlReactionParameters object."""

        self.C = C
        self.S = S
        self.X1 = X1
        self.X2 = X2
        self.X3 = X3
        self.X4 = X4
        self.Q = Q

    def __repr__( self ) :
        """Creates a string of self."""

        s = "C = %2d, S = %3d, X1 = %14.7e, X2 = %14.7e, X3 = %14.7e, X4 = %14.7e, Q = %14.7e" \
            % ( self.C, self.S, self.X1, self.X2, self.X3, self.X4, self.Q )
        return( s )

    def __getitem__( self, index ) :

        if( type( index ) == type( '' ) ) :
            if( index == 'C' ) : return( self.C )
            if( index == 'S' ) : return( self.S )
            if( index == 'X1' ) : return( self.X1 )
            if( index == 'X2' ) : return( self.X2 )
            if( index == 'X3' ) : return( self.X3 )
            if( index == 'X4' ) : return( self.X4 )
            if( index == 'Q' ) : return( self.Q )
            raise Exception( 'invalid string index = %s' % index )
        elif( type( index ) == type( 0 ) ) :
            if( index == 0 ) : return( self.C )
            if( index == 1 ) : return( self.S )
            if( index == 2 ) : return( self.X1 )
            if( index == 3 ) : return( self.X2 )
            if( index == 4 ) : return( self.X3 )
            if( index == 5 ) : return( self.X4 )
            if( index == 6 ) : return( self.Q )
            raise Exception( 'invalid integer index = %d' % index )
        raise Exception( 'invalid index type = %s; index must be string or integer' % type( index ) )

    def toPythonList( self ) :
        """Returns the [ self.C, self.S, self.X1, self.X2, self.X3, self.X4, self.Q ]"""

        return( [ self.C, self.S, self.X1, self.X2, self.X3, self.X4, self.Q ] )
