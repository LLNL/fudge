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

from . import base
import fudge
from fudge.gnd.tokens import *
from fudge.gnd.reactionData import crossSection
from fudge.core.math._xData import XYs, axes

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
        xys.setAncestor( pw )
        xys.isPrimaryXData = False
        pw.append( xys )
    return pw


def ResonancesWithBackground( data ):
    def __init__(self): pass
