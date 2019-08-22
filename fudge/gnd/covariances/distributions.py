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

from xData.ancestry import ancestry
from fudge.gnd import suites
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
        self.endfConversionFlag = ''

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

    def convertUnits( self, unitMap ):

        for lvalue in self: lvalue.convertUnits( unitMap )

    def fix( self, **kw ): return []

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ '%s<%s label="%s"' % ( indent, self.moniker, self.label ) ]
        if self.endfConversionFlag: xmlString[0] += ' endfConversionFlag="%s"' % self.endfConversionFlag
        xmlString[-1] += '>'
        for lvalue in self.lvalues : xmlString += lvalue.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):

        xPath.append( element.tag )
        form = LegendreOrderCovarianceForm( label = element.get( "label" ) )
        if element.get('endfConversionFlag') is not None:
            form.endfConversionFlag = element.get('endfConversionFlag')
        # add L's to each component:
        for lValue in element:
            form.addLegendreOrder( LegendreLValue.parseXMLNode( lValue, xPath, linkData ) )
        xPath.pop()
        return form


class LegendreLValue( suites.suite ):
    """ 
    Represents one subsection of the Legendre coefficient covariance matrix:
    covariance between coefficients for two Legendre orders at various energies 
    """

    moniker = 'LegendreLValue'

    def __init__(self, L1, L2, frame):
        suites.suite.__init__( self, [covarianceMatrix, mixed.mixedForm] )
        self.L1 = L1 #:
        self.L2 = L2 #:
        self.frame = frame #:

    def check( self, info ):

        warnings = []
        if self.L1 == self.L2:
            for form in self:
                warnings += form.check( info )
        # FIXME what about cross-terms between L-values?
        return warnings

    def fix( self, **kw ): return []

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlString = [ '%s<%s L1="%i" L2="%i" frame="%s">' %
            ( indent, self.moniker, self.L1, self.L2, self.frame ) ]
        for form in self:
            xmlString += form.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):

        xPath.append( element.tag )
        component = LegendreLValue( int(element.get("L1")), int(element.get("L2")),
                element.get("frame") )
        for form in element:
            formCls = {
                    covarianceMatrix.moniker: covarianceMatrix,
                    mixed.mixedForm.moniker: mixed.mixedForm
                    }.get( form.tag )
            if formCls is None:
                raise TypeError("Encountered unknown covariance type '%s'" % form.tag)
            component.add( formCls.parseXMLNode( form, xPath, linkData ) )
        xPath.pop()
        return component
