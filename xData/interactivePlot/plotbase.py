# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
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

    options = {'delete': delete}
    for minMaxOutput in [minMax('x', xMin, xMax), minMax('y', yMin, yMax)]:
        options[minMaxOutput[0]] = minMaxOutput[1]
        options[minMaxOutput[2]] = minMaxOutput[3]

    for key, value in (('xLabel', xLabel), ('yLabel', yLabel), ('title', title), ('style', style)):
        options[key] = '' if value is None else str(value)

    return options
