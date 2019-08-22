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
from xData import XYs as XYsModule
from xData import gridded as griddedModule

from fudge.processing import group as groupModule

from fudge.gnd import miscellaneous as miscellaneousModule
from fudge.gnd import abstractClasses as abstractClassesModule
from fudge.gnd import tokens as tokensModule

__metaclass__ = type

#
# Q forms
#
class baseQForm( abstractClassesModule.form ) :

    pass

class constant( baseQForm ) :

    moniker = 'constant'

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

        from fudge.gnd.reactions import base as reactionsBaseModule

        return( self.findClassInAncestry( reactionsBaseModule.base_reaction ).domainMin( unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        from fudge.gnd.reactions import base as reactionsBaseModule

        return( self.findClassInAncestry( reactionsBaseModule.base_reaction ).domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def domain( self, unitTo = None, asPQU = False ) :

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ), self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def domainUnit( self ) :

        from fudge.gnd.reactions import base as reactionsBaseModule

        return( self.findClassInAncestry( reactionsBaseModule.base_reaction ).domainUnit( ) )

    def getValue( self, E, unit ) :

        return( self.Q.getValueAs( unit ) )

    @property
    def label( self ) :

        return( self.__label )

    def processSnMultiGroup( self, style, tempInfo, indent ) :

        return( self.toPointwise_withLinearXYs( 1e-8, 1e-8 ).processSnMultiGroup( style, tempInfo, indent ) )

    def toXMLList( self, indent = '', **kwargs ) :

        attributeStr = ''
        if( self.label is not None ) : attributeStr += ' label="%s"' % self.label
        return( [ '%s<%s%s value="%s"/>' % ( indent, self.moniker, attributeStr, self.Q ) ] )

    def toPointwise_withLinearXYs( self, lowerEps, upperEps ) :
        """This method returns the Q-value as linear-linear XYs1d data which spans self's domain."""

        axes = XYs1d.defaultAxes( energyUnit = self.domainUnit( ), QUnit = self.Q.getUnitSymbol( ) )
        return( XYs1d( data = [ [ self.domainMin( ), self.Q.getValue( ) ], [ self.domainMax( ), self.Q.getValue( ) ] ], axes = axes, accuracy = 1e-12 ) )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ):

        xPath.append( element.tag )
        Q = constant( element.get( 'label' ), element.get( 'value' ) )
        xPath.pop( )
        return( Q )

class XYs1d( baseQForm, XYsModule.XYs1d ) :

    mutableYUnit = False

    def __init__( self, **kwargs ) :

        baseQForm.__init__( self )
        XYsModule.XYs1d.__init__( self, **kwargs )

    def processSnMultiGroup( self, style, tempInfo, indent ) :

        from fudge.processing import miscellaneous as miscellaneousModule

        indent2 = indent + tempInfo['incrementalIndent']
        verbosity = tempInfo['verbosity']

        QGrouped = miscellaneousModule.groupOneFunctionAndFlux( style, tempInfo, self )
        return( groupModule.toMultiGroup1d( multiGroup, style, tempInfo, self.axes, QGrouped ) )

    @staticmethod
    def defaultAxes( energyUnit = 'eV', QUnit = 'eV' ) :

        axes = axesModule.axes( )
        axes[0] = axesModule.axis( 'energy_in', 0, energyUnit )
        axes[1] = axesModule.axis( component.moniker, 1, QUnit )
        return( axes )

class multiGroup( baseQForm, griddedModule.gridded ) :

    def __init__( self, **kwargs ) :

        griddedModule.gridded.__init__( self, **kwargs )
#
# Q component
#
class component( abstractClassesModule.component ) :

    moniker = 'Q'

    def __init__( self ) :

        abstractClassesModule.component.__init__( self,
                ( constant, XYs1d, multiGroup ) )

    def getConstantAs( self, unit ) :

        for form in self :
            if( isinstance( form, constant ) ) : return( form.getValue( 0, unit ) )
        raise ValueError( 'Q-value is not a contant' )

    def getValue( self, E, unit ) :

        return( self.getEvaluated( ).getValue( E, unit ) )

def parseXMLNode( QElement, xPath, linkData ):
    """
    Reads an xml <Q> element into fudge, including all child forms
    """

    xPath.append( QElement.tag )
    kwargs = {}     # In case we need special interpolation rules for Q, see gnd/reactionData/crossSection.py

    Q = component( )
    for form in QElement :
        formClass = {
                constant.moniker        : constant,
                XYs1d.moniker           : XYs1d,
                multiGroup.moniker      : multiGroup
            }.get( form.tag )
        if( formClass is None ) : raise Exception( "unknown Q form: %s" % form.tag )
        newForm = formClass.parseXMLNode( form, xPath, linkData, **kwargs )
        Q.add( newForm )
    xPath.pop( )
    return( Q )
