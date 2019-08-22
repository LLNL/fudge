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

from fudge.gnd import baseClasses as baseClassesModule
from fudge.gnd import tokens

import xData.XYs as XYsModule
import xData.regions as regionsModule

__metaclass__ = type

energyDepositionToken = 'depositionEnergy'
momentumDepositionToken = 'depositionMomentum'
multiplicityToken = 'multiplicity'

class XYPointwiseFormBase( baseClassesModule.formBase, XYsModule.XYs ) :

    mutableYUnit = False

    def __init__( self, **kwargs ) :

        XYsModule.XYs.__init__( self, **kwargs )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        from fudge.processing import miscellaneous

        forms = []

        if( 'LLNL_Pn' in processInfo['styles'] ) :
            projectile, target = processInfo.getProjectileName( ), processInfo.getTargetName( )
            norm = tempInfo['groupedCrossSectionNorm']
            crossSection = tempInfo['crossSection']
            axes = self.axes.copy( )
            grouped = miscellaneous.groupTwoFunctionsAndFlux( projectile, processInfo, crossSection, self, norm = norm )
            forms = [ self.toForms[ tokens.groupedFormToken ]( axes, grouped ) ]
            if( tokens.groupedWithCrossSectionFormToken in self.toForms ) :
                norm = tempInfo['groupedFlux']
                axes = axes.copy( )
                axes[1].label = "%s %s" % ( axes[1].getLabel( ), crossSection.axes[1].getLabel( ) )
                if( axes[1].getUnit( ) == '' ) :
                    axes[1].unit = "%s" % ( crossSection.axes[1].getUnit( ) )
                else :
                    axes[1].unit = "%s * %s" % ( axes[1].getUnit( ), crossSection.axes[1].getUnit( ) )
                grouped = miscellaneous.groupTwoFunctionsAndFlux( projectile, processInfo, crossSection, self, norm = norm )
                forms.append( self.toForms[tokens.groupedWithCrossSectionFormToken]( axes, grouped ) )
        return( forms )

    def toXMLList( self, indent = "", **kwargs ) :

        from pqu import PQU

        def xyFormatter( x, y ) :

            return( "%s %s" % ( PQU.toShortestString( x, significantDigits = 9 ), PQU.toShortestString( y, significantDigits = 9 ) ) )

        if( 'xyFormatter' not in kwargs ) : kwargs['xyFormatter'] = xyFormatter
        return( XYsModule.XYs.toXMLList( self, indent, **kwargs ) )

class XYPiecewiseFormBase( baseClassesModule.formBase, regionsModule.regions ) :

    mutableYUnit = False

    def __init__( self, **kwargs ) :

        kwargs['dimension'] = 1
        regionsModule.regions.__init__( self, **kwargs )
