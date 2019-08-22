import fudge
import fudge.legacy.endl as endl

def deltafn( x, h=1.0 ):
    """
Returns an endl2dmath object of a fake delta function
with width from taken from last digit of endl4dmath format and area equal to h in size.
"""
    xx = float(x)
    hh = float(h)
    dx = nudgex( endl.endl4dmathmisc.endl4d_repr_xFormat, xx )
    minx = xx-dx
    midx = xx
    maxx = xx+dx
    return endl.endl2dmathmisc.triangle2d( minx, midx, maxx )*(hh/dx)

def nudgex( format, x ) :
    """Returns value that is equal to the smallest value that a given format can see a change in."""
    s = format % x
    p = int( format.split('.')[1].split('e')[0] )
    ls = len( s )
    expon = int( s[ls-3:ls] )
    nudge = 10.0**( -p+expon )
    return nudge
