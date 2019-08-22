import fudge, math, os
import fudge.legacy.endl.endl1dmathClasses 
import fudge.legacy.endl.endl2dmathClasses 
import fudge.legacy.endl.endl3dmathClasses 
from endlLikeDeltaFunc import nudgex
from fudge.core.math.xData import axes, XYs, W_XYs


class geft3dmath( fudge.legacy.endl.endl3dmathClasses.endl3dmath ) :
    def __init__( self, *args, **kw ) : fudge.legacy.endl.endl3dmathClasses.endl3dmath.__init__( self, *args, **kw )

    def __add__( self, v ) :
        if isinstance( v, fudge.legacy.endl.endl3dmathClasses.endl3dmath ) :
            d = self.copyData()
            s = self.copyData() # need separate copy so as to not affect interpolation as you add data
            r = v.copyData()
            x_ = s.xArray().union( r.xArray() )
            for x in x_ :
                yz_s = fudge.legacy.endl.endl2dmathClasses.endl2dmath( s.getAtX( x ) )
                yz_r = fudge.legacy.endl.endl2dmathClasses.endl2dmath( r.getAtX( x ) )
                d.setAtX( x, yz_s + yz_r )
        else :
            raise Exception( "\n(%s + %s) not supported" % (self.__class__, v.__class__) )
        return d

    def __sub__( self, v ) :
        if isinstance( v, fudge.legacy.endl.endl3dmathClasses.endl3dmath ) :
            d = self + ( -v )
        else :
            raise Exception( "\n(%s + %s) not supported" % (self.__class__, v.__class__) )
        return d

    def __neg__( self ) :
        d = self.copyData( )
        for x,yz in d.data :
            yz = fudge.legacy.endl.endl2dmathClasses.endl2dmath( yz )
            d.setAtX( x, -yz )
        return d

    def __mul__( self, v ) :
        if isinstance( v, ( int, float ) ) :
            v = float( v )
            d = self.copyData( )
            for x_yz in d.data :
                for yz in x_yz[1] :
                    yz[1] = yz[1] * v
        elif isinstance( v, fudge.legacy.endl.endl2dmathClasses.endl2dmath ) :
            d = self.copyData( )
            s = self.copyData( ) # need separate copy so as to not affect interpolation as you mul data
            for x in d.xArray( ) :
                yz_s = fudge.legacy.endl.endl2dmathClasses.endl2dmath( s.getAtX( x ) )
                z_v = v.getValue( x )
                d.setAtX( x, yz_s * z_v )
        else :
            raise Exception( "\n(%s * %s) not supported" % (self.__class__, v.__class__) )
        return d

    __rmul__ = __mul__

    def xArray( self ) :
        xa = []
        for xy in self.data : xa.append( xy[0] )
        return geft1dmath( xa )

    def copyData( self ) :
        d3 = []
        for x_yz in self.data :
            d2 = []
            for yz in x_yz[1] : d2.append( [ yz[0], yz[1] ] )
            d3.append( [ x_yz[0], d2 ] )
        return geft3dmath( d3, checkDataType = 0, template = self )

    def set( self, other ) :
        self.data = other.data

    def normalize( self ) :
        d = geft3dmath( [] )
        for x,yz in self.data :
            d.setAtX( x, fudge.legacy.endl.endl2dmathClasses.endl2dmath( yz ).normalize() )
        return d

    def slice( self, xMin = None, xMax = None ) :
        d = geft3dmath( [], checkDataType = 0, template = self )
        if ( len( self.data ) > 0 ) :
            if ( xMin == None ) : xMin = self.xArray()[0]
            if ( xMax == None ) : xMax = self.xArray()[-1]
            for x,yz in self.data :
                if xMin <= x <= xMax : d.setAtX( x, yz )
        return d

class geft2dmath( fudge.legacy.endl.endl2dmathClasses.endl2dmath ) :
    def __init__( self, *args, **kw ) : fudge.legacy.endl.endl2dmathClasses.endl2dmath.__init__( self, *args, **kw )

    def __mul__( self, v ) :
        if isinstance( v, geft3dmath ) :
            return geft3dmath.__mul__( v, self )
        else :
            return geft2dmath( fudge.legacy.endl.endl2dmathClasses.endl2dmath.__mul__( self, v ) )

class geft1dmath( fudge.legacy.endl.endl1dmathClasses.endl1dmath ) :
    def __init__( self, *args, **kw ) : fudge.legacy.endl.endl1dmathClasses.endl1dmath.__init__( self, *args, **kw )

    def union( self, other ) :
        d = self.copyData()
        for x in other.data :
            if x not in d.data : d.data.append( x )
        d.data.sort()
        return d

