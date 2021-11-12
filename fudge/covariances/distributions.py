# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from xData.ancestry import ancestry
from fudge import suites as suitesModule
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

        from fudge import warning

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


class LegendreLValue( suitesModule.suite ):
    """ 
    Represents one subsection of the Legendre coefficient covariance matrix:
    covariance between coefficients for two Legendre orders at various energies 
    """

    moniker = 'LegendreLValue'

    def __init__(self, L1, L2, frame):
        suitesModule.suite.__init__( self, [covarianceMatrix, mixed.mixedForm] )
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
