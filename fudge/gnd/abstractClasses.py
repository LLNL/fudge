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

from xData import ancestry as ancestryModule

from . import suites as suitesModule
from fudge.core.math import fudge2dGrouping as fudge2dGroupingModule

__metaclass__ = type

class component( suitesModule.suite ) :

    def __init__( self, allowedClasses ) :

        suitesModule.suite.__init__( self, allowedClasses )

    @property
    def evaluated( self ) :

        evaluated = self.getRootAncestor( ).styles.getEvaluatedStyle( )
        return( self[evaluated.label] )

    def genre( self ) :

        return( self.__genre )

    def getStylesOfClass( self, cls ) :

        styles = self.getRootAncestor( ).styles.getStylesOfClass( cls )
        formList = []
        for form in self :
            for style in styles :
                if( form.label == style.label ) : formList.append( form.label )
        return( formList )

    def getStyleOfClass( self, cls ) :

        style = self.getRootAncestor( ).styles.getStyleOfClass( cls )
        if( style is not None ) : style = self[style.label]
        return( style )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        raise Exception( "'%s' does not support processing" % self.genre )

    def toPointwise_withLinearXYs( self, lowerEps = 1.e-8, upperEps = 1.e-8 ) :
        """This method calls the toPointwise_withLinearXYs method for the evaluated style and returns a "lin,lin" pointwise representation of it."""

        return( self.evaluated.toPointwise_withLinearXYs( lowerEps, upperEps ) )

class form( ancestryModule.ancestry ) :
    """
    This is the base class which is used as a form base class for channelData, reactionData and productData form classes.
    """

    __genre = None

    @property
    def genre( self ) :

        return( self.__genre )

    def getComponentsClass( self ) :

        import sys
        return( sys.modules[self.__module__].component )

    def process( self, processInfo, tempInfo, verbosityIndent ) :

        raise Exception( "processing of '%s' data from '%s' to '%s' is not supported" % ( self.genre, self.style, processInfo( 'to' ) ) )

class multiGroup( form, fudge2dGroupingModule.groupedData ) :
    """This is the base class which is used as a form grouped base class for channelData, reactionData and productData form classes."""

    moniker = 'multiGroup'

    def __init__( self, style, axes_, groupData ) :

        fudge2dGroupingModule.groupedData.__init__( self, groupData, axes_ )
