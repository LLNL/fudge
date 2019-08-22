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

from fudge.gnd import baseClasses
from fudge.core.math.xData import axes, XYs
from fudge.gnd import tokens

__metaclass__ = type

energyDepositionToken = 'depositionEnergy'
momentumDepositionToken = 'depositionMomentum'
multiplicityToken = 'multiplicity'

class XYPointwiseFormBase( baseClasses.formBase, XYs.XYs ) :

    form = tokens.pointwiseFormToken
    tag = tokens.pointwiseFormToken
    mutableYUnit = False

    def __init__( self, axes_, data, accuracy, **kwargs ) :

        baseClasses.formBase.__init__( self )
        kwargs['isPrimaryXData'] = True
        XYs.XYs.__init__( self, axes_, data, accuracy, **kwargs )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        from fudge.processing import miscellaneous

        forms = []

        if( 'LLNL_Pn' in processInfo['styles'] ) :
            projectile, target = processInfo.getProjectileName( ), processInfo.getTargetName( )
            norm = tempInfo['groupedCrossSectionNorm']
            crossSection = tempInfo['crossSection']
            axes_ = self.axes.copy( )
            axes_[0].interpolation = axes.interpolationXY( axes.linearToken, axes.flatToken )
            grouped = miscellaneous.groupTwoFunctionsAndFlux( projectile, processInfo, crossSection, self, norm = norm )
            forms = [ self.toForms[ tokens.groupedFormToken ]( axes_, grouped ) ]
            if( tokens.groupedWithCrossSectionFormToken in self.toForms ) :
                norm = tempInfo['groupedFlux']
                axes_ = axes_.copy( )
                axes_[1].label = "%s %s" % ( axes_[1].getLabel( ), crossSection.axes[1].getLabel( ) )
                if( axes_[1].getUnit( ) == '' ) :
                    axes_[1].unit = "%s" % ( crossSection.axes[1].getUnit( ) )
                else :
                    axes_[1].unit = "%s * %s" % ( axes_[1].getUnit( ), crossSection.axes[1].getUnit( ) )
                grouped = miscellaneous.groupTwoFunctionsAndFlux( projectile, processInfo, crossSection, self, norm = norm )
                forms.append( self.toForms[tokens.groupedWithCrossSectionFormToken]( axes_, grouped ) )
        return( forms )

    def toXMLList( self, indent = "" ) :

        from pqu.physicalQuantityWithUncertainty import toShortestString

        def xyFormatter( x, y ) :

            return( "%s %s" % ( toShortestString( x, precision = 9 ), toShortestString( y, precision = 9 ) ) )

        return( XYs.XYs.toXMLList( self, indent = indent, incrementalIndent = '  ', pairsPerLine = 100, xyFormatter = xyFormatter ) )
