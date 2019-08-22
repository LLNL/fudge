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
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
