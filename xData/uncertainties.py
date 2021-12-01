# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
    This module contains the class ``uncertainties``, which can appear inside other xData classes.

    Multiple types of uncertainty are supported:
        pointwise (stored in the XYs1d class),
        links,
        covariances,
        ...

    May contain more than one uncertainty, e.g. to support asymmetric uncertainties
"""

__metaclass__ = type

from . import ancestry as ancestryModule
from . import base as baseModule
from . import link as linkModule

class uncertainties( baseModule.xDataCoreMembers ):

    moniker = 'uncertainties'

    def __init__(self, uncertainties):

        baseModule.xDataCoreMembers.__init__(self, self.moniker)
        self.__uncertainties = uncertainties or []

    def __getitem__(self, item):
        return self.__uncertainties[item]

    def __len__(self):
        return len( self.__uncertainties )

    def convertUnits( self, unitMap ) :

        for uncertainty in self.__uncertainties : uncertainty.convertUnits( unitMap )

    def toXMLList( self, indent = '', **kwargs ) :

        if len(self) == 0: return( [] )

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLList = ['%s<%s>' % (indent, self.moniker)]
        for uncertainty in self : XMLList += uncertainty.toXMLList( indent2, **kwargs )
        XMLList[-1] += '</%s>' % self.moniker
        return( XMLList )

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ) :
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

    def convertUnits( self, unitMap ) :

        self.__functional.convertUnits( unitMap )

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
    def parseXMLNode( cls, element, xPath, linkData ) :
        """
        Translate <uncertainty> element into uncertainty instance.
        """

        from . import XYs as XYsModule
        from . import series1d as series1dModule

        xPath.append( element.tag )
        kwargs = {}
        for attribute in ( 'index', 'label', 'type', 'pdf', 'relation' ) :
            kwargs[attribute] = element.get(attribute,None)
        if len(element) != 1:
            raise TypeError("uncertainty element must contain exactly one functional")
        functionalClass = {
            linkModule.link.moniker: linkModule.link,
            XYsModule.XYs1d.moniker: XYsModule.XYs1d,
            series1dModule.polynomial1d.moniker: series1dModule.polynomial1d,
            covariance.moniker: covariance,
            listOfCovariances.moniker: listOfCovariances,
        }.get( element[0].tag )
        kwargs['functional'] = functionalClass.parseXMLNode( element[0], xPath, linkData )
        uncertainty_ = uncertainty( **kwargs )
        xPath.pop()
        return uncertainty_


class covariance(linkModule.link):

    moniker = 'covariance'

class listOfCovariances(ancestryModule.ancestry):

    moniker = 'listOfCovariances'

    def __init__(self):
        ancestryModule.ancestry.__init__(self)
        self.__items = []

    def __getitem__(self, item):
        return self.__items[item]

    def add(self, obj):
        if not isinstance(obj, covariance):
            raise TypeError("Expected covariance instance, got '%s'" % type(obj))
        obj.setAncestor(self)
        obj.label = 'cov%d' % len(self.__items)
        self.__items.append(obj)

    def toXMLList(self, indent='', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent','  ')
        xml = ['%s<%s>' % (indent, self.moniker)]
        for item in self:
            xml += item.toXMLList(indent2, **kwargs)
        xml[-1] += '</%s>' % self.moniker

        return xml

    @classmethod
    def parseXMLNode(cls, element, xPath, linkData):

        xPath.append(element.tag)
        CL = cls()
        for child in element:
            covClass = {
                covariance.moniker: covariance
            }.get( child.tag )
            CL.add( covClass.parseXMLNode(child, xPath, linkData) )

        xPath.pop()
        return CL
