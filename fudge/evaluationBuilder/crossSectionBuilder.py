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

from . import base
import fudge
from fudge.gnd.tokens import *
from fudge.gnd.reactionData import crossSection
from fudge.core.math.xData import XYs, axes

class crossSectionBuilder:
    def __init__( self, data, form = "Pointwise" ):
        self.component = crossSection.component()

        xs_types = {
                'Pointwise': Pointwise,
                'Piecewise': Piecewise,
                'ResonancesWithBackground': ResonancesWithBackground,
                }
        builder = xs_types.get( form, None )
        if builder is None:
            raise Exception( "'%s' is not a recognized form of cross section. Options are '%s'"
                    % ( form, xs_types.keys() ) )
        xs = builder( data )
        self.component.addForm( xs, asNativeData = True )

    def __call__(self): return self.component

    def setInterpolation( self, x_interpolation, y_interpolation, region=0 ):
        """
        By default each cross section region uses linear-linear interpolation.
        Use this method to change interpolation.
        """
        token = self.component.getNativeDataToken()
        form = self.component[ token ]

        if token == resonancesWithBackgroundFormToken:
            form = self.component[ token ].tabulatedData
            token = form.moniker

        if token == pointwiseFormToken:
            if region!=0: raise Exception("Only one region present in pointwise data!")
            form.axes[0].setInterpolation( axes.interpolationXY(x_interpolation, y_interpolation) )
        elif token == piecewiseFormToken:
            form[region].axes[0].setInterpolation( axes.interpolationXY(x_interpolation, y_interpolation) )
        else:
            raise Exception("Don't know how to set interpolation for %s cross section" % token)


def Pointwise( data ):
    """
    Create a pointwise cross section. By default use lin-lin interpolation, and units of 'eV' for incident energy
    and 'b' for the cross section.

    'data' should be a nested list of the form
    [(incident_energy_0, cross_section_0), (incident_energy_1, cross_section_1), ... ]
    """
    axes_ = crossSection.pointwise.defaultAxes()
    return crossSection.pointwise( axes_, data, accuracy = base.defaultAccuracy )


def Piecewise( data ):
    """
    Create a piecewise cross section, with multiple regions. Each region uses lin-lin interpolation, but
    the setInterpolation method of the crossSectionBuilder class can be used to change interpolation.

    'data' should be contain a list of regions:
    [
        [ # data for the first region
            (incident_energy_0, cross_section_0), (incident_energy_1, cross_section_1) ...
            ]
        [ # second region
            (incident_energy_0, cross_section_0), ...
            ]
        ...
    ]
    """
    axes_ = crossSection.piecewise.defaultAxes()
    pw = crossSection.piecewise( axes_ )
    for index, region in enumerate(data):
        ax1 = axes.interpolationAxes( 0, axes.interpolationXY(linearFormToken,linearFormToken), parent=None )
        xys = XYs.XYs( ax1, region, accuracy = base.defaultAccuracy )
        xys.parent = pw; xys.isPrimaryXData = False
        pw.append( xys )
    return pw


def ResonancesWithBackground( data ):
    def __init__(self): pass
