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

"""Base classes for distributions."""

from fudge.core.utilities import brb

import xData.ancestry as ancestryModule
import xData.standards as standardsModule

from fudge.gnds import abstractClasses as abstractClassesModule

__metaclass__ = type

#
# Standard genre can only be composited from one of the following standard forms.
#

angularTwoBodyGenre = 'angularTwoBody'
unknownGenre = 'unknown'
referenceGenre = 'reference'
NBodyGenre = 'NBody'

noneFormToken = 'none'
semiPiecewiseFormToken = 'semiPiecewise'
equalProbableBinsFormToken = 'equalProbableBins'
groupedFormToken = 'grouped'
LegendrePointwiseFormToken = 'LegendrePointwise'
LegendrePiecewiseFormToken = 'LegendrePiecewise'

class form( abstractClassesModule.form ) :

    def __init__( self, label, productFrame, subForms ) :

        abstractClassesModule.form.__init__( self )

        self.label = label

        if( productFrame is not None ) :
            if( productFrame not in standardsModule.frames.allowedFrames ) :
                raise TypeError( 'Invalid productFrame = "%s"' % brb.limitObjectToString( productFrame ) )
        self.__productFrame = productFrame

        for i1, subform_ in enumerate( subForms ) :
            setattr( self, self.subformAttributes[i1], subform_ )
            if( subform_ is not None ) : subform_.setAncestor( self )
        self.subforms = [ subform_ for subform_ in subForms ]

    @property
    def label( self ) :

        return( self.__label )

    @label.setter
    def label( self, value ) :

        if( value is not None ) :
            if( not( isinstance( value, str ) ) ) : raise TypeError( 'value must be a string' )
        self.__label = value

    @property
    def productFrame( self ) :

        return( self.__productFrame )

    def convertUnits( self, unitMap ) :
        """See documentation for reactionSuite.convertUnits."""

        if( hasattr( self, 'subforms' ) ) :
            for subform_ in self.subforms :
                if( subform_ is None ) : continue                    # FIXME, needed for photo-atomic coherentScattering data with no anomalousScattering subforms.
                subform_.convertUnits( unitMap )

    def findEntity( self, entityName, attribute = None, value = None ):
        """
        Find specific subform within the form. Overrides ancestry.findEntity.
        """

        for subform_ in self.subforms :
            if( subform_.moniker == entityName ) : return( subform_ )
        return( ancestryModule.ancestry.findEntity( self, entityName, attribute, value ) )

    def isTwoBody( self ) :

        return( False )

    def isSpecified( self ) :

        return( True )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        return( self.subforms[0].toPointwise_withLinearXYs( **kwargs ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attributeStr = ''
        if( self.label is not None ) : attributeStr += ' label="%s"' % self.label
        if( self.productFrame is not None ) : attributeStr += ' productFrame="%s"' % self.productFrame
        xmlString = [ '%s<%s%s>' % ( indent, self.moniker, attributeStr ) ]
        for subform_ in self.subforms : xmlString += subform_.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

class subform( ancestryModule.ancestry ) :

    def __init__( self ) :

        ancestryModule.ancestry.__init__( self )
        self.label = None

    def XMLStartTagString( self, indent = '', extraAttributesAsStrings = '', emptyTag = False ) :

        emptyTagStr = ''
        if( emptyTag ) : emptyTagStr = '/'
        attributeStr = ''
        if( self.label is not None ) : attributeStr += ' label="%s"' % self.label
        return( '%s<%s%s%s%s>' % ( indent, self.moniker, attributeStr, extraAttributesAsStrings, emptyTagStr ) )
