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

"This module is deprecated."

import xData.ancestry as ancestryModule

__metaclass__ = type

class componentBase( ancestryModule.ancestry ) :
    """This class is deprecated and is only used for fissionEnergyReleased."""

    def __init__( self ) :

        ancestryModule.ancestry.__init__( self )
        self.forms = {}

    def __len__( self ) :
        "Returns the number of forms in self."

        return( len( self.forms ) )

    def __getitem__( self, label ) :
        "Returns self's form with label."

        return( self.forms[label] )

    def addForm( self, form ) :

        if( form.label in self.forms ) : raise Exception( 'label = %s already in self' % form.label )
        self.forms[form.label] = form
        form.setAncestor( self )

    def getEvaluated( self ) :

        evaluated = self.getRootAncestor( ).styles.getEvaluatedStyle( )
        return( self.forms[evaluated.label] )

    def removeStyle( self, form ) :
        """
        This method removes form from self. Form can be a label or a form instance.
        """

        if( isinstance( form, formBase ) ) :
            found = False
            for f in self.forms :
                if( self.forms[f] is form ) : found = True
            if( not( found ) ) : raise Exception( 'Form "%s" of "%s" is not a form of self of "%s"' % ( form.label, form.moniker, self.moniker ) )
            form = form.label
        return( self.forms.pop( form, None ) )

    def toPointwise_withLinearXYs( self, lowerEps = 1.e-8, upperEps = 1.e-8 ) :
        """This method calls the toPointwise_withLinearXYs method for the evaluated style and returns a "lin-lin" pointwise representation of it."""

        evaluated = self.getEvaluated( )
        return( evaluated.toPointwise_withLinearXYs( lowerEps, upperEps ) )

    def toXML( self, indent = "", **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = "", **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        formXMLString = []
        for label in self.forms : formXMLString += self.forms[label].toXMLList( indent2, **kwargs )
        if( formXMLString == [] ) : return( [] )
        xmlString = [ '%s<%s>' % ( indent, self.moniker ) ] + formXMLString
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

class formBase( ancestryModule.ancestry ) :
    """This class is deprecated."""

    def getComponentsClass( self ) :

        import sys
        return( sys.modules[self.__module__].component )

    def getForm( self ) :

        return( self.label )
