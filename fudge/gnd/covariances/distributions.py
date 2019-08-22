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

from xData.ancestry import ancestry

from . import mixed, covarianceMatrix

"""Special classes for storing covariances for product distributions."""

__metaclass__ = type

class LegendreOrderCovarianceForm( ancestry ):
    """ 
    Stores covariance between energy-dependent Legendre coefficients for a reaction.
    This class contains one or more LegendreLValue sections, each section containing the matrix
    between a pair of L-values 
    """
    
    moniker = 'LegendreOrderCovariance'

    def __init__(self, label = None, lvalues=None):
        ancestry.__init__( self )
        self.__label = label
        self.lvalues = lvalues or [] #: the l values of course

    def __getitem__(self, index):   return self.lvalues[index]

    @property
    def label( self ) :

        return( self.__label )

    def addLegendreOrder( self, LValue ):
        LValue.setAncestor( self )
        self.lvalues.append( LValue )

    def check( self, info ):

        from fudge.gnd import warning

        warnings = []
        for Lvalue in self:
            Lvalue_warnings = Lvalue.check( info )
            if Lvalue_warnings:
                warnings.append( warning.context( '%s L=%i vs %i' % (Lvalue.moniker,Lvalue.L1,Lvalue.L2), Lvalue_warnings ) )
        return warnings

    def fix( self, **kw ): return []

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ '%s<%s label="%s">' % ( indent, self.moniker, self.label ) ]
        for lvalue in self.lvalues : xmlString += lvalue.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):

        xPath.append( element.tag )
        form = LegendreOrderCovarianceForm( label = element.get( "label" ) )
        # add L's to each component:
        for lValue in element:
            form.addLegendreOrder( LegendreLValue.parseXMLNode( lValue, xPath, linkData ) )
        xPath.pop()
        return form


class LegendreLValue( ancestry ):
    """ 
    Represents one subsection of the Legendre coefficient covariance matrix:
    covariance between coefficients for two Legendre orders at various energies 
    """

    moniker = 'LegendreLValue'

    def __init__(self, L1, L2, frame, nativeData=None):
        ancestry.__init__( self )
        self.L1 = L1 #:
        self.L2 = L2 #:
        self.frame = frame #:
        self.nativeData = nativeData #:
        self.forms = {} #:

    def addForm( self, form ):
        form.setAncestor( self )
        self.forms[ form.moniker ] = form

    def check( self, info ):

        warnings = []
        if self.L1 == self.L2:
            warnings += self.forms[self.nativeData].check( info )
        # FIXME what about cross-terms between L-values?
        return warnings

    def fix( self, **kw ): return []

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ '%s<%s L1="%i" L2="%i" frame="%s" nativeData="%s">' % \
            ( indent, self.moniker, self.L1, self.L2, self.frame, self.nativeData ) ]
        for form in self.forms.values():
            xmlString += form.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):

        xPath.append( element.tag )
        component = LegendreLValue( int(element.get("L1")), int(element.get("L2")),
                element.get("frame"), element.get("nativeData") )
        for form in element:
            formCls = {
                    covarianceMatrix.moniker: covarianceMatrix,
                    mixed.mixedForm.moniker: mixed.mixedForm
                    }.get( form.tag )
            if formCls is None:
                raise TypeError("Encountered unknown covariance type '%s'" % form.tag)
            component.addForm( formCls.parseXMLNode( form, xPath, linkData ) )
        xPath.pop()
        return component

