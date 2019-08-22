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

import base as baseModule

from fudge.gnd import tokens as tokensModule
from fudge.gnd import abstractClasses as abstractClassesModule

import xData.axes as axesModule

__metaclass__ = type

class baseMomentumDepositionForm( abstractClassesModule.form ) :

    __genre = baseModule.momentumDepositionToken

class grouped( baseMomentumDepositionForm, abstractClassesModule.multiGroup ) :

    def __init__( self, label, axes, data ) :

        baseMomentumDepositionForm.__init__( self )
        abstractClassesModule.multiGroup.__init__( self, label, axes, data )

class pointwise( baseMomentumDepositionForm, baseModule.XYPointwiseFormBase ) :

    def __init__( self, **kwargs ) :

        baseMomentumDepositionForm.__init__( self )
        baseModule.XYPointwiseFormBase.__init__( self, **kwargs )

    @staticmethod
    def defaultAxes( energyUnit = 'eV', momentumDepositionName = 'momentumDeposition', momentumDepositionUnit = 'eV/c' ) :

        axes = axesModule.axes( rank = 2 )
        axes[0] = axesModule.axis( momentumDepositionName, 0, momentumDepositionUnit )
        axes[1] = axesModule.axis( 'energy_in', 1, energyUnit )
        return( axes )

class piecewise( baseMomentumDepositionForm, baseModule.XYPiecewiseFormBase ) :

    def __init__( self, **kwargs ) :

        baseMomentumDepositionForm.__init__( self )
        baseModule.XYPiecewiseFormBase.__init__( self, **kwargs )

    @staticmethod
    def allowedSubElements( ) :

        return( ( pointwise, ) )

    @staticmethod
    def defaultAxes( energyUnit = 'eV', momentumDepositionName = 'momentumDeposition', momentumDepositionUnit = 'eV/c' ) :

        return( pointwise.defaultAxes( energyUnit = energyUnit, momentumDepositionName = momentumDepositionName,
                momentumDepositionUnit = momentumDepositionUnit ) )

#
# momentumDeposition component
#
class component( abstractClassesModule.component ) :

    __genre = baseModule.momentumDepositionToken
    moniker = baseModule.momentumDepositionToken

    def __init__( self ) :

        abstractClassesModule.component.__init__( self, ( grouped, pointwise, piecewise ) )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        for form in self :
            if( form.label == tokensModule.pointwiseFormToken ) :
                ps = form.process( processInfo, tempInfo, verbosityIndent )
                for p in ps : self.add( p )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :
        """Reads an xml <depositionMomentum> component element into fudge, including all forms in the component."""

        xPath.append( element.tag )
        momentumDepositionComponent = component( )
        for form in element :
            formClass = {
                    pointwise.moniker : pointwise,
                    grouped.moniker :   grouped,
               }.get( form.tag )
            if( formClass is None ) : raise Exception( "encountered unknown depositionMomentum form: %s" % form.tag )
            try :
                newForm = formClass.parseXMLNode( form, xPath, linkData )
            except Exception as e :
                raise Exception, "/depositionMomentum/%s: %s" % ( form.tag, e )
            momentumDepositionComponent.add( newForm )
        xPath.pop()
        return( momentumDepositionComponent )
