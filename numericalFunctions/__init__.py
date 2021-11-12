# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import glob
import pathlib

# Deployment via Makefiles places shared libraries in .../fudge/numericalFunctions/lib
if len(glob.glob(os.path.join(pathlib.Path(__file__).parent.absolute(), 'lib', '*pointwiseXY*'))) > 0:
    from .lib import pointwiseXY_C
    from .lib import Legendre
    from .lib.pointwiseXY_C import pointwiseXY_C as pointwiseXY
    from .lib import specialFunctions
    from .lib import angularMomentumCoupling
    from .lib import integration
    from .lib import listOfDoubles_C

# Deployment via `pip install` places shared libraries in .../site-packages/numericalFunctions
else:
    from . import pointwiseXY_C
    from . import Legendre
    from .pointwiseXY_C import pointwiseXY_C as pointwiseXY
    from . import specialFunctions
    from . import angularMomentumCoupling
    from . import integration
    from . import listOfDoubles_C
