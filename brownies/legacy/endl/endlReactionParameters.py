# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
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
