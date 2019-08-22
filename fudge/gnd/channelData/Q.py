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
from fudge.gnd import miscellaneous
from fudge.core.math.xData import axes, XYs
from fudge.gnd import tokens

__metaclass__ = type

energyDependent = 'energyDependent'

QForms = [ tokens.notApplicableFormToken, tokens.constantFormToken, tokens.pointwiseFormToken, tokens.groupedFormToken, tokens.groupedWithCrossSectionFormToken ]

#
# Q genre and forms
#
class component( baseClasses.componentBase ) :

    genre = base.QToken

    def __init__( self, form = None ) :

        baseClasses.componentBase.__init__( self, QForms )
        if( form is not None ) : self.addForm( form )

    def getConstantAs( self, unit ) :

        if( self.nativeData in [ tokens.constantFormToken ]  ) : return( self.getValue( 0, unit ) )
        raise Exception( 'Q type = %s does not have a single value' % self.nativeData )

    def getValue( self, E, unit ) :

        return( self.forms[self.nativeData].getValue( E, unit ) )

    def getXMLAttribute( self ) :

        return( self.forms[self.nativeData].getXMLAttribute( ) )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        keys = self.forms.keys( )
        for form in keys :
            if( form == tokens.pointwiseFormToken ) :
                ps = self.forms[form].process( processInfo, tempInfo, verbosityIndent )
                for p in ps : self.addForm( p )

class constant( baseClasses.formBase ) :

    genre = component.genre
    form = tokens.constantFormToken

    def __init__( self, Q ) :
        """
        Q can be a string, convertible to a PhysicalQuantityWithUncertainty, or a PhysicalQuantityWithUncertainty. The
        Q must have units of energy.
        """

        from pqu.physicalQuantityWithUncertainty import PhysicalQuantityWithUncertainty
        from fudge.core.utilities import brb

        Q_ = Q
        if( type( Q_ ) == str ) :
            try :
                Q_ = PhysicalQuantityWithUncertainty( Q_ )
            except :
                raise Exception( "Could not convert '%s' to PhysicalQuantityWithUncertainty" % Q_ )
        if( isinstance( Q_, component ) ) : Q_ = Q_.getNativeData( )
        if( isinstance( Q_, constant ) ) : Q_ = Q_.Q
        if( not( isinstance( Q_, PhysicalQuantityWithUncertainty ) ) ) : raise Exception( "Invalid type for Q: type = '%s'" % brb.getType( Q ) )
        if( not( Q_.isEnergy( Q_ ) ) ) : raise Exception( "Q must have units of energy, entered Q was '%s'" % Q )
        baseClasses.formBase.__init__( self )
        self.Q = Q_

    def domainMin( self, unitTo = None, asPQU = False ) :

        from fudge.gnd.reactions import reaction
        return( self.findClassInAncestry( reaction.reaction ).domainMin( unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        from fudge.gnd.reactions import reaction
        return( self.findClassInAncestry( reaction.reaction ).domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getDomain( self, unitTo = None, asPQU = False ) :

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ), self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getDomainUnit( self ) :

        from fudge.gnd.reactions import reaction
        return( self.findClassInAncestry( reaction.reaction ).getDomainUnit( ) )

    def getValue( self, E, unit ) :

        return( self.Q.getValueAs( unit ) )

    def getXMLAttribute( self ) :

        return( self.Q )

    def toXMLList( self, indent ) :

        return( [] )

    def toPointwise_withLinearXYs( self, lowerEps, upperEps ) :
        """This method returns the Q-value as linear-linear pointwise data which spans self's domain."""

        axes_ = pointwise.defaultAxes( energyUnit = self.getDomainUnit( ), QUnit = self.Q.getUnitSymbol( ) )
        return( pointwise( axes_, [ [ self.domainMin( ), self.Q.getValue( ) ], [ self.domainMax( ), self.Q.getValue( ) ] ], 1e-12 ) )

class pointwise( baseClasses.formBase, XYs.XYs ) :

    genre = component.genre
    form = tokens.pointwiseFormToken
    tag = tokens.pointwiseFormToken
    mutableYUnit = False

    def __init__( self, axes, data, accuracy, **keyWords ) :

        baseClasses.formBase.__init__( self )
        keyWords['isPrimaryXData'] = True
        XYs.XYs.__init__( self, axes, data, accuracy, **keyWords )
        self.toForms = { tokens.groupedFormToken : grouped }

    def process( self, processInfo, tempInfo ) :

        projectile, target = processInfo.getProjectileName( ), processInfo.getTargetName( )
        groups = processInfo.getParticleGroups( projectile )
        Q_E = self.getFormByToken( tokens.pointwiseFormToken )
        grouped = miscellaneous.makeGrouped( self, processInfo, tempInfo, Q_E )
        self.addForm( grouped( grouped ) )
        grouped = miscellaneous.makeGrouped( self, processInfo, tempInfo, Q_E, normType = 'groupedFlux' )
        self.addForm( groupedWithCrossSection( grouped ) )

    def getXMLAttribute( self ) :

        return( energyDependent )

    @staticmethod
    def defaultAxes( energyUnit = 'eV', energyInterpolation = axes.linearToken, QUnit = 'eV', QInterpolation = axes.linearToken ) :

        axes_ = axes.axes( )
        axes_[0] = axes.axis( 'energy_in', 0, energyUnit, interpolation = axes.interpolationXY( energyInterpolation, QInterpolation ) )
        axes_[1] = axes.axis( base.QToken, 1, QUnit )
        return( axes_ )

class grouped( baseClasses.groupedFormBase ) :

    genre = component.genre

    def __init__( self, axes, data ) :

        baseClasses.groupedFormBase.__init__( self, axes, data )

class groupedWithCrossSection( baseClasses.groupedWithCrossSectionFormBase ) :

    genre = component.genre

    def __init__( self, axes, data ) :

        baseClasses.groupedWithCrossSectionFormBase.__init__( self, axes, data )

class notApplicable( baseClasses.formBase ) :

    genre = component.genre
    form = tokens.notApplicableFormToken

    def __init__( self ) :

        baseClasses.formBase.__init__( self )
        self.Q = 'N/A'

    def getValue( self, E, unit ) :

        return( 0. )

    def getXMLAttribute( self ) :

        return( self.Q )

    def toXMLList( self, indent ) :

        return( [] )
