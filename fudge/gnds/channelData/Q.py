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

from pqu import PQU as PQUModule

from xData import axes as axesModule
from xData import constant as constantModule
from xData import XYs as XYsModule
from xData import gridded as griddedModule

from fudge.gnds import abstractClasses as abstractClassesModule
from . import fissionEnergyReleased as fissionEnergyReleasedModule

__metaclass__ = type

def defaultAxes( energyUnit ) :

    axes = axesModule.axes( )
# BRB6 make 'energy_in' a token
    axes[1] = axesModule.axis( 'energy_in', 0, energyUnit )
    axes[0] = axesModule.axis( component.moniker, 1, energyUnit )
    return( axes )

#
# Q forms
#
class baseQForm( abstractClassesModule.form ) :

    pass

class constant1d( baseQForm, constantModule.constant1d ) :

    def __init__( self, Q, domainMin, domainMax, axes, label = None ) :

        baseQForm.__init__( self )
        constantModule.constant1d.__init__( self, Q, domainMin, domainMax, axes = axes, label = label )

    def processMultiGroup( self, style, tempInfo, indent ) :

        return( self.toPointwise_withLinearXYs( accuracy = 1e-5, upperEps = 1e-8 ).processMultiGroup( style, tempInfo, indent ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """This method returns the Q-value as linear-linear XYs1d data which spans self's domain."""

        return( XYs1d( data = [ [ self.domainMin, self.constant ], [ self.domainMax, self.constant ] ], axes = self.axes ) )

class XYs1d( baseQForm, XYsModule.XYs1d ) :

    mutableYUnit = False

    def __init__( self, **kwargs ) :

        baseQForm.__init__( self )
        XYsModule.XYs1d.__init__( self, **kwargs )

    def processMultiGroup( self, style, tempInfo, indent ) :

        from fudge.processing import miscellaneous as miscellaneousModule

        return( miscellaneousModule.groupFunctionCrossSectionAndFlux( gridded1d, style, tempInfo, self ) )

class gridded1d( baseQForm, griddedModule.gridded1d ) :

    def __init__( self, **kwargs ) :

        griddedModule.gridded1d.__init__( self, **kwargs )
#
# Q component
#
class component( abstractClassesModule.component ) :

    moniker = 'Q'

    def __init__( self ) :

        abstractClassesModule.component.__init__( self,
                ( constant1d, XYs1d, gridded1d, fissionEnergyReleasedModule.fissionEnergyReleased ) )

    def thresholdQAs( self, unit ) :

        for form in self :
            if( isinstance( form, constant1d ) ) :
                return( PQUModule.PQU( form.constant, form.axes[0].unit ).getValueAs( unit ) )
            elif( isinstance( form, fissionEnergyReleasedModule.fissionEnergyReleased ) ) :
                return( PQUModule.PQU( form.nonNeutrinoEnergy.data.evaluate( 0 ) ) )
        raise ValueError( 'Q-value is not a constant' )

    def getConstantAs( self, unit ) :

        return( self.thresholdQAs( unit ) )

    def getConstant( self ):

        for form in self :
            if( isinstance( form, constant1d ) ) :
                return form.constant
            elif( isinstance( form, fissionEnergyReleasedModule.fissionEnergyReleased ) ) :
                return form.nonNeutrinoEnergy.data.evaluate( 0 )
        raise ValueError( 'Q-value is not a constant' )

    def evaluate( self, E ) :

        return( self.getEvaluated( ).evaluate( E ) )

def parseXMLNode( QElement, xPath, linkData ):
    """
    Reads an xml <Q> element into fudge, including all child forms
    """
    from .fissionEnergyReleased import fissionEnergyReleased

    xPath.append( QElement.tag )
    kwargs = {}     # In case we need special interpolation rules for Q, see gnds/reactionData/crossSection.py

    Q = component( )
    for form in QElement :
        formClass = {
                constant1d.moniker     : constant1d,
                XYs1d.moniker          : XYs1d,
                gridded1d.moniker      : gridded1d,
                fissionEnergyReleased.moniker : fissionEnergyReleased
            }.get( form.tag )
        if( formClass is None ) : raise Exception( "unknown Q form: %s" % form.tag )
        newForm = formClass.parseXMLNode( form, xPath, linkData, **kwargs )
        Q.add( newForm )
    xPath.pop( )
    return( Q )
