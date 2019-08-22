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
from fudge.core.ancestry import ancestry
from pqu.physicalQuantityWithUncertainty import PhysicalQuantityWithUncertainty
from . import mixed, tokens, covarianceMatrix

"""Special classes for storing covariances for product distributions."""

class energyIntervalForm(mixed.mixedForm):
    """ 
    For distributions, the covariances may depend on incident energy. Similar to mixed,
    but each component has an associated energy 
    """

    moniker = tokens.energyIntervalFormToken

    def __init__(self, components=None):
        self.components = components or [] #: a list of components that are instances of ``mixedForm`` or ``covarianceMatrix``

    @classmethod
    def parseXMLNode( cls, element, xPath=[], linkData={} ):

        xPath.append( element.tag )
        form = super(energyIntervalForm, cls).parseXMLNode( element, xPath, linkData )
        # add energy bounds to each component:
        for i in range(len(element)):
            component = element[i]
            lower = PhysicalQuantityWithUncertainty( component.get("lowerBound") )
            upper = PhysicalQuantityWithUncertainty( component.get("upperBound") )
            form.components[i].energyBounds = (lower,upper)
        xPath.pop()
        return form

class LegendreOrderCovarianceForm( ancestry ):
    """ 
    Stores covariance between energy-dependent Legendre coefficients for a reaction.
    This class contains one or more LegendreLValue sections, each section containing the matrix
    between a pair of L-values 
    """
    
    moniker = tokens.legendreOrderCovarianceFormToken

    def __init__(self, lvalues=None):
        ancestry.__init__( self, tokens.legendreOrderCovarianceFormToken, None )
        self.lvalues = lvalues or [] #: the l values of course

    def __getitem__(self, index):   return self.lvalues[index]

    def addLegendreOrder( self, LValue ):
        LValue.setParent( self )
        self.lvalues.append( LValue )

    def check( self, info ): return []
    
    def fix( self, **kw ): return []

    def toXMLList(self, flags=None, indent=''):
        xmlString = [indent+'<%s>' % self.moniker]
        for lvalue in self.lvalues:
            xmlString += lvalue.toXMLList(flags, indent+'  ')
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    def toENDF6(self, flags, targetInfo):
        from fudge.legacy.converting import endfFormats
        endf = []
        for lVal in self.lvalues:
            LCT = {'frameOfMF4':0, 'lab':1, 'centerOfMass':2}[ lVal.frame ]
            form = lVal.forms[ lVal.nativeData ]
            NI = 1
            if form.moniker == tokens.mixedFormToken: NI = len(form)
            endf.append( endfFormats.endfHeadLine(0,0,lVal.L1, lVal.L2, LCT, NI) )
            endf += form.toENDF6(flags, targetInfo)
        return endf

    @classmethod
    def parseXMLNode( cls, element, xPath=[], linkData={} ):

        xPath.append( element.tag )
        form = LegendreOrderCovarianceForm()
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

    moniker = tokens.legendreLValueFormToken

    def __init__(self, L1, L2, frame, nativeData=None):
        ancestry.__init__( self, tokens.legendreLValueFormToken, None )
        self.L1 = L1 #:
        self.L2 = L2 #:
        self.frame = frame #:
        self.nativeData = nativeData #:
        self.forms = {} #:

    def addForm( self, form ):
        form.setParent( self )
        self.forms[ form.moniker ] = form

    def check( self, info ): return []
    
    def fix( self, **kw ): return []

    def toXMLList(self, flags=None, indent=''):
        xmlString = [indent+'<%s L1="%i" L2="%i" frame="%s" nativeData="%s">' % (self.moniker,
            self.L1, self.L2, self.frame, self.nativeData)]
        for form in self.forms.values():
            xmlString += form.toXMLList(flags, indent+'  ')
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseXMLNode( cls, element, xPath=[], linkData={} ):

        xPath.append( element.tag )
        component = LegendreLValue( int(element.get("L1")), int(element.get("L2")),
                element.get("frame"), element.get("nativeData") )
        for form in element:
            formCls = {
                    tokens.covarianceFormToken: covarianceMatrix,
                    tokens.mixedFormToken: mixed.mixedForm
                    }.get( form.tag )
            if formCls is None:
                raise BadCovariance
            component.addForm( formCls.parseXMLNode( form, xPath, linkData ) )
        xPath.pop()
        return component

