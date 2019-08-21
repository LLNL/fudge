# <<BEGIN-copyright>>
# <<END-copyright>>

"""
top-level __init__ for the fudge project
cmattoon, 3/18/2011
"""
from __future__ import division, nested_scopes

__version__ = '4.0.0'

__docformat__ = 'epytext en'
from __doc__ import __doc__

import fudgeDefaults
import fudgeParameters
from core import *
import vis
import gnd
import particles

# if we want to export a smaller set of files with 'from fudge import *':
#__all__ = ['core','gnd',...]
