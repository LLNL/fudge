# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from .endlProject import *
from .endl2 import *
from .endlZA import *
from .fudgemiscLegacy import *
from fudge.core.utilities.fudgeExceptions import *
from fudge.vis.gnuplot.fudgeMultiPlots import *
from .endl1dmathClasses import *
from .endl2dmathClasses import *
from .endl3dmathClasses import *
from .endl4dmathClasses import *
from fudge.core.math.fudgemath import isNumber
from .endlmath import fastSumOfManyAddends
from .fudgemisc import *

def readOnly( ) :
    """Sets Fudge's internal readonly flag which causes fudge to not create working directories
    for the isotopes read in. Note, you cannot unset this flag once readOnly has been called."""

    from . import fudgeParameters

    fudgeParameters.ReadOnly = 1

def verbose( mode = 1 ) :
    """Sets Fudge's internal verbose flag to mode. Mode = 0 suppress all informational messages."""

    from . import fudgeParameters

    fudgeParameters.VerboseMode = int( mode )
