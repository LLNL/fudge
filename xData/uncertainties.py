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

"""
    This module contains the class ``uncertainties``, which can appear inside other xData classes.

    Multiple types of uncertainty are supported:
        pointwise (stored in the XYs class),
        links,
        covariances,
        ...

    May contain more than one uncertainty, e.g. to support asymmetric uncertainties
"""

__metaclass__ = type

import base as baseModule
import link as linkModule
import XYs as XYsModule

class uncertainties( baseModule.xDataCoreMembers ):

    moniker = 'uncertainties'

    def __init__(self, uncertainties):

        baseModule.xDataCoreMembers.__init__(self, self.moniker)
        self.__uncertainties = uncertainties or []

    def __getitem__(self, item):
        return self.__uncertainties[item]

    def __len__(self):
        return len( self.__uncertainties )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLList = ['%s<%s>' % (indent, self.moniker)]
        for uncertainty in self : XMLList += uncertainty.toXMLList( indent2, **kwargs )
        XMLList[-1] += '</%s>' % self.moniker
        return( XMLList )

    @classmethod
    def parseXMLNode( cls, element, xPath = [], linkData = {} ) :
        """
        Translate <uncertainties> element into uncertainties instance.
        """

        xPath.append( element.tag )
        uncertaintyList = [uncertainty.parseXMLNode( child, xPath, linkData ) for child in element ]
        uncertainties_ = uncertainties( uncertaintyList )
        xPath.pop()
        return uncertainties_


class uncertainty( baseModule.xDataCoreMembers ):

    moniker = 'uncertainty'

    defaultType='variance'
    defaultPdf='normal'
    defaultRelation='offset'

    def __init__(self, index=None, label=None, functional=None, type=None, pdf=None,
                 relation=None):

        baseModule.xDataCoreMembers.__init__(self, self.moniker, index=index, label=label )
        #FIXME check that each of the following are in the allowed list:
        self.__functional = functional
        self.__type = type or self.defaultType
        self.__pdf = pdf or self.defaultPdf
        self.__relation = relation or self.defaultRelation

    @property
    def type(self):

        return self.__type

    @property
    def pdf(self):

        return self.__pdf

    @property
    def relation(self):

        return self.__relation

    @property
    def data(self):

        return self.__functional

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attributes = ''
        for attr in ( 'index', 'label' ) :
            if getattr(self,attr) is not None:
                attributes += ' %s="%s"' % (attr, getattr(self,attr))
        if self.type != self.defaultType: attributes += ' type="%s"' % self.type
        if self.pdf != self.defaultPdf: attributes += ' pdf="%s"' % self.pdf
        if self.relation != self.defaultRelation: attributes += ' relation="%s"' % self.relation
        XMLList = ['%s<%s%s>' % (indent, self.moniker, attributes)]
        XMLList += self.data.toXMLList( indent2, **kwargs )
        XMLList[-1] += '</%s>' % self.moniker
        return( XMLList )

    @classmethod
    def parseXMLNode( cls, element, xPath = [], linkData = {} ) :
        """
        Translate <uncertainty> element into uncertainty instance.
        """

        xPath.append( element.tag )
        kwargs = {}
        for attribute in ( 'index', 'label', 'type', 'pdf', 'relation' ) :
            kwargs[attribute] = element.get(attribute,None)
        if len(element) != 1:
            raise TypeError("uncertainty element must contain exactly one functional")
        functionalClass = {
            linkModule.link.moniker: linkModule.link,
            XYsModule.XYs.moniker: XYsModule.XYs,
        }.get( element[0].tag )
        kwargs['functional'] = functionalClass.parseXMLNode( element[0] )
        uncertainty_ = uncertainty( **kwargs )
        xPath.pop()
        return uncertainty_
