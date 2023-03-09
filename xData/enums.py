# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from LUPY import enums as enumsModule

class StorageOrder(enumsModule.Enum):
    '''
    StorageOrder is an enum class representing the way data can be stored in a table or an array. The two values
    for StorageOrder is 'row-major' or 'column-major' and are represented by the names rowMajor and columnMajor.
    '''

    rowMajor = 'row-major'
    columnMajor = 'column-major'

class ValueType(enumsModule.Enum):
    '''
    This enum spcifies the value types for the *Values* class.
    '''

    integer32 = 'Integer32'
    float32 = 'Float32'
    float64 = 'Float64'

class Interpolation(enumsModule.Enum):
    '''
    This enum specifies the interpolation values for functions.
    '''

    linlin = 'lin-lin'
    linlog = 'lin-log'
    loglin = 'log-lin'
    loglog = 'log-log'
    flat = 'flat'
    chargedParticle = 'charged-particle'

class InterpolationQualifier(enumsModule.Enum):
    '''
    This enum specifies the interpolation qualifier values for functions.
    '''
    # We also need direct, correspondingEnergies and correspondingEnergies-unscaled.

    none = ''
    unitBase = 'unitbase'
    unitBaseUnscaled = 'unitbase-unscaled'
    correspondingPoints = 'correspondingPoints'
    correspondingPointsUnscaled = 'correspondingPoints-unscaled'
    cumulativePoints = 'cumulativePoints'
    cumulativePointsUnscaled = 'cumulativePointsUnscaled'

class Frame(enumsModule.Enum):
    '''
    This enum specifies the allowed frame values.
    '''

    none = enumsModule.auto()
    lab = enumsModule.auto()
    centerOfMass = enumsModule.auto()
    product = enumsModule.auto()

class FixDomain(enumsModule.Enum):
    '''
    This enum specifies if the lower, upper or both domains are to be fixed for a function.
    '''

    upper = enumsModule.auto()
    lower = enumsModule.auto()
    both = enumsModule.auto()

class Extrapolation(enumsModule.Enum):
    '''
    This enum specifies supported extrapolations for some functions.
    '''

    none = enumsModule.auto()
    flat = enumsModule.auto()

class GridStyle(enumsModule.Enum):
    '''
    This enum specifies allowed styles for the Grid class.
    '''

    none = enumsModule.auto()
    points = enumsModule.auto()
    boundaries= enumsModule.auto()
    parameters= enumsModule.auto()
