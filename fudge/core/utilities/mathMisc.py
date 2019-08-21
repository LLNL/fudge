# <<BEGIN-copyright>>
# <<END-copyright>>

def shiftFloatABit( f, s, r_eps, a_eps, z_eps ) :
    "Only for internal use by shiftFloatDownABit and shiftFloatUpABit."

    if( a_eps is None ) : a_eps = r_eps
    if( z_eps is None ) : z_eps = r_eps
    if( r_eps < 0. ) : raise Exception( 'r_eps = %s < 0.' % r_eps )
    if( r_eps > 0.5 ) : raise Exception( 'r_eps = %s > 0.5' % r_eps )
    if( a_eps < 0. ) : raise Exception( 'a_eps = %s < 0.' % a_eps )
    if( a_eps > 0.5 ) : raise Exception( 'a_eps = %s > 0.5' % a_eps )
    if( z_eps < 0. ) : raise Exception( 'z_eps = %s < 0.' % z_eps )
    if( z_eps > 0.5 ) : raise Exception( 'z_eps = %s > 0.5' % z_eps )
    if( f < 0. ) : r_eps = -r_eps
    if( f == 0. ) : return( s * z_eps )
    return( ( 1 + s * r_eps ) * f + s * a_eps )

def shiftFloatDownABit( f, r_eps, a_eps = 0., z_eps = None ) :
    """
    Returns a float that is slightly less than f. The amount less depends on r_eps, a_eps and z_eps.
    If f is 0., returns -z_eps; otherwise returns ( 1 - r_eps ) * f - a_eps. r_eps, a_eps and z_eps
    all must be between 0. and 0.5 inclusive."""

    return( shiftFloatABit( f, -1, r_eps, a_eps, z_eps ) )

def shiftFloatUpABit( f, r_eps, a_eps = 0., z_eps = None ) :
    """
    Returns a float that is slightly more than f. The amount more depends on r_eps, a_eps and z_eps.
    If f is 0., returns z_eps; otherwise returns ( 1 + r_eps ) * f + a_eps. r_eps, a_eps and z_eps
    all must be between 0. and 0.5 inclusive."""

    return( shiftFloatABit( f, 1, r_eps, a_eps, z_eps ) )
