# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
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
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
# <<END-copyright>>

from . import base
import fudge
from fudge.gnd.productData import distributions
from fudge.core.math._xData import XYs, axes
from fudge.core.math._xData.axes import labToken, centerOfMassToken

class distributionBuilder:
    def __init__( self, dataForm, data, frame=None ):
        dist_types = { 
                'AngularLegendre': AngularLegendre,
                'AngularPointwise': AngularPointwise,
                'AngularMixed': AngularMixed,
                'Isotropic': Isotropic,
                'Recoil': Recoil,
                'EnergyPointwise': EnergyPointwise,
                'Uncorrelated': Uncorrelated,
                'KalbachMann': KalbachMann,
                }
        componentBuilder = dist_types.get( dataForm, None )
        if componentBuilder is None:
            raise Exception( "'%s' is not a recognized form of distribution. Options are '%s'"
                    % (dataForm, dist_types.keys()) )
        opts = {}
        if frame: opts['frame'] = frame
        self.component = componentBuilder( data, **opts )

    def __call__(self): return self.component

# angular distributions:

def AngularLegendre( data, frame=centerOfMassToken ):
    """
    Create an angular distribution with Legendre coefficients. 'data' should be a list of tuples:
    [ (incident_energy_0, [list of coefficients]),
      (incident_energy_1, [list of coefficients]),
      ... ]

    A different number of Legendre coefficients may be supplied at each incident energy.
    The interpolation is assumed to be linear-linear, energies are assumed to be in eV.
    """

    axes_ = distributions.angular.LegendrePointwise.defaultAxes( energyInterpolation='lin', dependentInterpolation='lin' )
    lpw = distributions.angular.LegendrePointwise( axes_, frame )
    for energy, coefs in data:
        lpw.append( fudge.core.math._xData.LegendreSeries.XYs_LegendreSeries( coefficients=coefs, value=energy ) )
    component = distributions.angular.component()
    component.addForm( lpw )
    component.nativeData = lpw.name
    return component

def AngularPointwise( data, frame=centerOfMassToken ):
    """
    Create a pointwise angular distribution. 'data' should be a nested list:
    [ (incident_energy_0, [(mu0, P(mu0)), (mu1, P(mu1)), ... ]), 
      (incident_energy_1, [(mu0, P(mu0)), (mu1, P(mu1)), ... ]), 
      ...  ]

    A different mu grid may be used at each incident energy, but mu must always cover the range [-1:1]
    The interpolation is assumed to be linear-linear, energies are assumed to be in eV.
    """

    axes_ = distributions.angular.pointwise.defaultAxes( )
    pw = distributions.angular.pointwise( axes_, frame )
    for index, (energy, mu_list) in enumerate(data):
        axes_ = axes.referenceAxes( parent=pw )
        pw.append( XYs.XYs( axes_, mu_list, accuracy=base.defaultAccuracy, index=index, value=energy, parent=pw ) )
    component = distributions.angular.component()
    component.addForm( pw )
    component.nativeData = pw.name
    return component

def AngularMixed( data, frame=centerOfMassToken ):
    """
    Legendre coefficients are an efficient way to represent an angular distribution, but if the distribution
    is very forward-peaked (which is often true at high incident energies) we have trouble with negative probabilities
    due to truncating the Legendre series. One solution is to use the Mixed form, where Legendre coefficients are used
    at low incident energies, pointwise form at higher energies. This class builds a 'mixed' angular distribution.

    'data' should be a dictionary:
    { 'LegendreData': [ list of lists, just like the list used in the AngularLegendre function ],
      'PointwiseData': [ another list of lists, similar to the argument to AngularPointwise ] }

    The highest incident energy in 'LegendreData' should equal the lowest incident energy in 'PointwiseData'.
    In both cases, interpolation is lin-lin and energies are in eV.
    """

    if type(data) is not dict:
        raise TypeError("AngularMixed requires a dictionary with both Legendre and Pointwise data")
    LegendreComponent = AngularLegendre( data['LegendreData'], frame )
    PointwiseComponent = AngularPointwise( data['PointwiseData'], frame )
    Mixed = distributions.angular.MixedForm( LegendreComponent, PointwiseComponent )
    component = distributions.angular.component()
    component.addForm( Mixed )
    component.nativeData = Mixed.name
    return component

def Isotropic( data=None, frame=centerOfMassToken ):
    """
    Create an isotropic angular distribution. No arguments are required.
    """
    isotropic = distributions.angular.isotropic( frame )
    component = distributions.angular.component()
    component.addForm( isotropic )
    component.nativeData = isotropic.name
    return component

def Recoil( otherProduct, frame=None ):
    """
    For 2-body reactions, the distribution of one product can be calculated from the other.
    'otherProduct' should be a gnd.product or productBuilder instance, from which this product is recoiling
    """
    recoil = distributions.angular.recoil( otherProduct )
    component = distributions.angular.component()
    component.addForm( recoil )
    component.nativeData = recoil.name
    return component

# energy spectra:

def EnergyPointwise( data, frame=labToken ):
    """
    Create a pointwise energy distribution. 'data' should be a nested list:
    [ (incident_energy_0, [(E'_0, P(E'_0)), (E'_1, P(E'_1)), ... ]), 
      (incident_energy_1, [(E'_0, P(E'_0)), (E'_1, P(E'_1)), ... ]), 
      ...  ]

    A different outgoing energy grid may be used at each incident energy.
    The interpolation is assumed to be linear-linear, and energies are assumed to be in eV.
    """

    axes_ = distributions.energy.pointwise.defaultAxes( )
    pw = distributions.energy.pointwise( axes_, frame )
    for index, (energy, probability_list) in enumerate(data):
        axes_ = axes.referenceAxes( parent=pw )
        pw.append( XYs.XYs( axes_, probability_list, accuracy=base.defaultAccuracy, index=index, value=energy, parent=pw ) )
    component = distributions.energy.component()
    component.addForm( pw )
    component.nativeData = pw.name
    return component

# double-differential distributions:

def Uncorrelated( data, frame=labToken ):
    angular = distributionBuilder( *data['Angular'], frame=frame )( )
    energy = distributionBuilder( *data['Energy'], frame=frame )( )
    component = distributions.uncorrelated.component( angular, energy )
    return component

def KalbachMann( data, frame=centerOfMassToken ):
    """
    Create a Kalbach-Mann distribution. 'data' should be a dictionary with two keys:
    { 'form': a string, either 'f' or 'fr', depending on how many parameters are given per incident/outgoing energy combination,
      'data': a nested list, of form
        [ (incident_energy_0, [ outgoing_energy_0, f, (r), outgoing_energy_1, f, (r) ... ] ),
          (incident_energy_1, [ outgoing_energy_0, f, (r), outgoing_energy_1, f, (r) ... ] ),
          ...  ]
          
    The second parameter (r) is only required if the 'form' is 'fr'.
    Energies are assumed to be in eV.
    """
    # We still need a 'KalbachMann.defaultAxes' function:
    axes_ = axes.axes( dimension=3 )
    axes_.axes.append( axes.axis( 'energy_in', index=0, unit='eV', 
        interpolation = axes.interpolationXY('lin','lin') ) )   # should be able to just say interpolation='lin-lin'
    axes_.axes.append( axes.axis( 'energy_out', index=1, unit='eV', 
        interpolation = axes.interpolationXY('lin','flat') ) )
    axes_.axes.append( axes.axis( 'f', index=2, unit='1/eV' ) )

    form, data = data['form'], data['data']
    KM = distributions.energyAngular.KalbachMann( form, axes_ )
    # beware, this isn't checking for correct number of parameters (if form=='f', len(coefficients) must be multiple of 2. if form=='fr', must be multiple of 3)
    for i in range(len(data)):
        KM.append( distributions.energyAngular.KalbachMannCoefficients( index=i, value=data[i][0], coefficients=data[i][1] ) )
    component = distributions.energyAngular.component()
    component.addForm( KM )
    component.nativeData = KM.name
    return component
