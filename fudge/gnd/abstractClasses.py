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

import sys

from xData import ancestry as ancestryModule

from . import suites as suitesModule
from fudge.core.math import fudge2dGrouping as fudge2dGroupingModule

__metaclass__ = type

class component( suitesModule.suite ) :

    def __init__( self, allowedClasses ) :

        suitesModule.suite.__init__( self, allowedClasses )

    @property
    def evaluated( self ) :
        if not hasattr(self.getRootAncestor(),'styles'): return self['eval'] # a hack to deal with orphaned components that are not part of a full reactionSuite
        evaluated = self.getRootAncestor( ).styles.getEvaluatedStyle( )
        return( self[evaluated.label] )

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

    def processSnMultiGroup( self, style, tempInfo, indent ) :

        indent2 = indent + tempInfo['incrementalIndent']
        verbosity = tempInfo['verbosity']
        addToComponent = tempInfo.get( 'addToComponent', True )

        form = style.findFormMatchingDerivedStyles( self )
        if( form is None ) :
#            if( verbosity ) : print '%sWARNING: no form found for component %s' % ( indent, self.moniker )
            return( None )
        try :
            multiGroup = form.processSnMultiGroup( style, tempInfo, indent2 )
        except :
            sys.stderr.write( 'form = %s\n' % form.moniker )
            raise

        if( addToComponent and ( multiGroup is not None ) ) : self.add( multiGroup )
        return( multiGroup )

    def toPointwise_withLinearXYs( self, lowerEps = 1.e-8, upperEps = 1.e-8 ) :
        """
        This method calls the toPointwise_withLinearXYs method for the evaluated style and returns a "lin-lin" pointwise representation of it.

        @:param lowerEps: lower epsilon for converting discontinuities to monotonic function
        @:param upperEps: upper epsilon for converting discontinuities to monotonic function
        """

        return( self.evaluated.toPointwise_withLinearXYs( lowerEps, upperEps ) )

class form( ancestryModule.ancestry ) :
    """
    This is the base class which is used as a form base class for channelData, reactionData and productData form classes.
    """

    def getComponentsClass( self ) :

        return( sys.modules[self.__module__].component )

class multiGroup( form, fudge2dGroupingModule.groupedData ) :
    """This is the base class which is used as a form grouped base class for channelData, reactionData and productData form classes."""

    def __init__( self, label, axes_, groupData ) :

        fudge2dGroupingModule.groupedData.__init__( self, groupData, axes_ )
        self.__label = label

    @property
    def label( self ) :

        return( self.__label )
