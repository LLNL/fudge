# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from xData import axes as axesModule
from xData import multiD_XYs as multiD_XYsModule
from xData import xs_pdf_cdf as xs_pdf_cdfModule

class Axes( axesModule.Axes ) :

    def __init__( self, energyUnit ) :

        axesModule.Axes.__init__( self, 3, labelsUnits = { 0 : ( 'P(cross section|energy_in)', '1/b' ), 1 : ( 'cross section', 'b' ), 2 : ( 'energy_in', energyUnit ) } )

class Xs_pdf_cdf1d( xs_pdf_cdfModule.Xs_pdf_cdf1d ) :

    pass

class XYs2d( multiD_XYsModule.XYs2d ) :

    @staticmethod
    def allowedSubElements( ) :

        return( ( Xs_pdf_cdf1d, ) )
