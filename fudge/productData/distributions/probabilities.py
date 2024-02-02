# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains probability functions (e.g., P(x1|x2), P(x2,x1|x3)) as needed for distributions. 
Currently, only P(x1|x2) exists as PofX1GivenX2.
"""

from xData import multiD_XYs as multiD_XYsModule

from . import miscellaneous as miscellaneousModule

class PofX1GivenX2( multiD_XYsModule.XYs2d ) :
    """
    Abstract base class for a 2d probability function of the form P(x1|x2).
    """

    def integrate( self, x2, x1Domain ) :
        """
        This method returns the integral over x1 for *self* evaluated at *x2*.
        """

        try :
            PofX1AtX2 = self.evaluate( x2 )
        except :
            return( 0.0 )

        x1Min, x1Max = miscellaneousModule.domainLimits( x1Domain, PofX1AtX2.domainMin, PofX1AtX2.domainMax )

        if( x1Max is None ) : return( PofX1AtX2.evaluate( x1Min ) )

        return PofX1AtX2.integrate(x1Min, x1Max)
