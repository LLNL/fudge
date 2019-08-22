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
import base
from fudge.core.math.xData import axes, XYs
from fudge.gnd import tokens

__metaclass__ = type

momentumDepositionForms = [ tokens.constantFormToken, tokens.partialProductionFormToken, tokens.pointwiseFormToken, tokens.groupedFormToken ]

#
# momentumDeposition genre and forms
#
class component( baseClasses.componentBase ) :

    genre = base.momentumDepositionToken

    def __init__( self, nativeData ) :

        baseClasses.componentBase.__init__( self, momentumDepositionForms, nativeData )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        keys = self.forms.keys( )
        for form in keys :
            if( form == tokens.pointwiseFormToken ) :
                ps = self.forms[form].process( processInfo, tempInfo, verbosityIndent )
                for p in ps : self.addForm( p )

class grouped( baseClasses.groupedFormBase ) :

    genre = component.genre

    def __init__( self, axes, data ) :

        baseClasses.groupedFormBase.__init__( self, axes, data )

class pointwise( base.XYPointwiseFormBase ) :

    genre = component.genre

    def __init__( self, axes, data, accuracy, **kwargs ) :

        base.XYPointwiseFormBase.__init__( self, axes, data, accuracy, **kwargs )
        self.toForms = { tokens.groupedFormToken : grouped }

    @staticmethod
    def defaultAxes( energyUnit = 'eV', energyInterpolation = axes.linearToken, momentumDepositionName = 'momentumDeposition',
        momentumDepositionInterpolation = axes.linearToken, momentumDepositionUnit = 'eV/c' ) :

        axes_ = axes.axes( dimension = 2 )
        axes_[0] = axes.axis( 'energy_in', 0, energyUnit, interpolation = axes.interpolationXY( energyInterpolation, momentumDepositionInterpolation ) )
        axes_[1] = axes.axis( momentumDepositionName, 1, momentumDepositionUnit )
        return( axes_ )

def parseXMLNode( element, xPath=[], linkData={} ):
    """Reads an xml <depositionMomentum> component element into fudge, including all forms in the component."""

    xPath.append( element.tag )
    dmc = component( element.get('nativeData') )
    for form in element:
        formClass = {tokens.pointwiseFormToken: pointwise,
                tokens.groupedFormToken: grouped,
                }.get( form.tag )
        if formClass is None: raise Exception(" encountered unknown depositionMomentum form: %s" % form.tag)
        try:
            newForm = formClass.parseXMLNode( form, xPath, linkData )
        except Exception as e:
            raise Exception, "/depositionMomentum/%s: %s" % (form.tag, e)
        dmc.addForm( newForm )
    xPath.pop()
    return dmc
