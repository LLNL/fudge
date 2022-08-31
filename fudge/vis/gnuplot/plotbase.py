# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

def minMax( s, _min, _max ) :

    if( _min is not None ) : _min = float( _min )
    if( _max is not None ) : _max = float( _max )
    if( ( _min is not None ) and ( _max is not None ) ) :
        if( _min == _max ) :
            if( _min == 0 ) :
                _min, _max = -1, 1
            else :
                _min = ( 1 - 1e-3 ) * _min
                _max = ( 1 + 1e-3 ) * _max
                if( _max < _min ) : _min, _max = _max, _min
    options = [ ]
    if( _min is not None ) : options += [ '%sMin' % s, str( _min ) ]
    if( _max is not None ) : options += [ '%sMax' % s, str( _max ) ]
    return( options )

def parsePlotOptions( xMin, xMax, yMin, yMax, xLabel, yLabel, title, zMin = None, zMax = None, zLabel = None, tLabel = None, \
        tScaleLabel = None, xrot = None, zrot = None, delete = True, style = None ) :
    """For internal use only."""
    import math
    options = [ ]
    if( delete ) : options = [ "delete" ]
    options += minMax( 'x', xMin, xMax )
    options += minMax( 'y', yMin, yMax )
    options += minMax( 'z', zMin, zMax )
    if ( xLabel is not None ) : options += [ 'xLabel', str( xLabel ) ]
    if ( yLabel is not None ) : options += [ 'yLabel', str( yLabel ) ]
    if ( zLabel is not None ) : options += [ 'zLabel', str( zLabel ) ]
    if ( tLabel is not None ) : options += [ 'tLabel', str( tLabel ) ]
    if ( title is not None ) : options += [ 'title', str( title ) ]
    if ( style is not None ) : options += [ 'style', str( style ) ]
    if ( tScaleLabel is not None ) : options += [ 'tScaleLabel', str( tScaleLabel ) ]
    if ( xrot is not None ) :
        if( type( xrot ) == type( "" ) ) : xrot = float( xrot )
        xrot_ = math.fmod( float( xrot ), 180 )
        if ( xrot_ < 0. ) : xrot_ += 180
        options += [ 'xrot', str( xrot_ ) ]
    if ( zrot is not None ) :
        if( type( zrot ) == type( "" ) ) : zrot = float( zrot )
        zrot_ = math.fmod( float( zrot ), 360 )
        if ( zrot_ < 0. ) : zrot_ += 360
        options += [ 'zrot', str( zrot_ ) ]
    return( options )
