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
        tScaleLabel = None, xrot = None, zrot = None, delete = True ) :
    """For internal use only."""

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
