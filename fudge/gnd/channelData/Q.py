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

from fudge.gnd import miscellaneous as miscellaneousModule
from fudge.gnd import abstractClasses as abstractClassesModule
from fudge.gnd import tokens as tokensModule

from pqu import PQU as PQUModule

import xData.axes as axesModule
import xData.XYs as XYsModule

__metaclass__ = type

#
# Q forms
#
class baseQForm( abstractClassesModule.form ) :

    __genre = baseModule.QToken

class constant( baseQForm ) :

    moniker = tokensModule.constantFormToken

    def __init__( self, label, Q ) :
        """
        Q can be a string, convertible to a PQU , or a PQU. The Q must have units of energy.
        """

        Q_ = Q
        if( isinstance( Q_, str ) ) :
            try :
                Q_ = PQUModule.PQU( Q_ )
            except :
                raise Exception( "Could not convert '%s' to PQU" % Q_ )
        if( isinstance( Q_, component ) ) : Q_ = Q_.getEvaluated( )
        if( isinstance( Q_, constant ) ) : Q_ = Q_.Q
        if( not( isinstance( Q_, PQUModule.PQU ) ) ) : raise TypeError( "Invalid type for Q: type = '%s'" % brb.getType( Q ) )
        if( not( Q_.isEnergy( ) ) ) : raise ValueError( "Q must have units of energy, entered Q was '%s'" % Q )
        self.Q = Q_

        if( label is not None ) :
            if( not( isinstance( label, str ) ) ) : raise TypeError( 'label must be a string' )
        self.__label = label

    def domainMin( self, unitTo = None, asPQU = False ) :

        from fudge.gnd.reactions import reaction
        return( self.findClassInAncestry( reaction.reaction ).domainMin( unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        from fudge.gnd.reactions import reaction
        return( self.findClassInAncestry( reaction.reaction ).domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def domain( self, unitTo = None, asPQU = False ) :

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ), self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def domainUnit( self ) :

        from fudge.gnd.reactions import reaction
        return( self.findClassInAncestry( reaction.reaction ).domainUnit( ) )

    def getValue( self, E, unit ) :

        return( self.Q.getValueAs( unit ) )

    @property
    def label( self ) :

        return( self.__label )

    def toXMLList( self, indent = '', **kwargs ) :

        attributeStr = ''
        if( self.label is not None ) : attributeStr += ' label="%s"' % self.label
        return( [ '%s<%s%s value="%s"/>' % ( indent, self.moniker, attributeStr, self.Q ) ] )

    def toPointwise_withLinearXYs( self, lowerEps, upperEps ) :
        """This method returns the Q-value as linear-linear pointwise data which spans self's domain."""

        axes = pointwise.defaultAxes( energyUnit = self.domainUnit( ), QUnit = self.Q.getUnitSymbol( ) )
        return( pointwise( [ [ self.domainMin( ), self.Q.getValue( ) ], [ self.domainMax( ), self.Q.getValue( ) ] ], axes = axes, accuracy = 1e-12 ) )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( element.tag )
        Q = constant( element.get( 'label' ), element.get( 'value' ) )
        xPath.pop( )
        return( Q )

class pointwise( baseQForm, XYsModule.XYs ) :

    mutableYUnit = False

    def __init__( self, **kwargs ) :

        baseQForm.__init__( self )
        XYsModule.XYs.__init__( self, **kwargs )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        projectile, target = processInfo.getProjectileName( ), processInfo.getTargetName( )
        groups = processInfo.getParticleGroups( projectile )
        grouped = miscellaneousModule.makeGrouped( self, processInfo, tempInfo, self )
        self.addForm( grouped( grouped ) )
        grouped = miscellaneousModule.makeGrouped( self, processInfo, tempInfo, Q_E, normType = 'groupedFlux' )
        self.addForm( groupedWithCrossSection( grouped ) )

    @staticmethod
    def defaultAxes( energyUnit = 'eV', QUnit = 'eV' ) :

        axes = axesModule.axes( )
        axes[0] = axesModule.axis( 'energy_in', 0, energyUnit )
        axes[1] = axesModule.axis( baseModule.QToken, 1, QUnit )
        return( axes )

class grouped( baseQForm, abstractClassesModule.multiGroup ) :

    moniker = tokensModule.groupedFormToken

    def __init__( self, axes, data ) :

        baseQForm.__init__( self )
        abstractClassesModule.multiGroup.__init__( self, axes, data )

class groupedWithCrossSection( baseQForm, abstractClassesModule.multiGroup ) :

    moniker = tokensModule.groupedFormToken

    def __init__( self, axes, data ) :

        baseQForm.__init__( self )
        abstractClassesModule.multiGroup.__init__( self, axes, data )

#
# Q component
#
class component( abstractClassesModule.component ) :

    __genre = baseModule.QToken
    moniker = baseModule.QToken

    def __init__( self ) :

        abstractClassesModule.component.__init__( self,
                ( constant, pointwise, grouped, groupedWithCrossSection ) )

    def getConstantAs( self, unit ) :

        for form in self :
            if( isinstance( form, constant ) ) : return( form.getValue( 0, unit ) )
        raise ValueError( 'Q-value is not a contant' )

    def getValue( self, E, unit ) :

        return( self.getEvaluated( ).getValue( E, unit ) )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        for form in self :
            if( isinstance( form, pointwise ) ) :
                ps = form.process( processInfo, tempInfo, verbosityIndent )
                for p in ps : self.addForm( p )

def parseXMLNode( QElement, xPath, linkData ):
    """
    Reads an xml <Q> element into fudge, including all child forms
    """

    xPath.append( QElement.tag )
    kwargs = {}     # In case we need special interpolation rules for Q, see gnd/reactionData/crossSection.py

    Q = component( )
    for form in QElement :
        formClass = {
                constant.moniker                : constant,
                pointwise.moniker               : pointwise,
                grouped.moniker                 : grouped,
                groupedWithCrossSection.moniker : groupedWithCrossSection
            }.get( form.tag )
        if( formClass is None ) : raise Exception( "unknown Q form: %s" % form.tag )
        newForm = formClass.parseXMLNode( form, xPath, linkData, **kwargs )
        Q.add( newForm )
    xPath.pop( )
    return( Q )
