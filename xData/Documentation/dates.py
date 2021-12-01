# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the classes representing the GNDS documentation nodes author and authors.
"""

from .. import ancestry as ancestryModule
from .. import suite as suiteModule
from .. import date as dateModule
from . import abstractClasses as abstractClassesModule

class DateType:

    allowed = ( 'accepted', 'available', 'copyrighted', 'collected', 'created', 'issued', 'submitted', 'updated', 'valid', 'withdrawn' )

class Date(ancestryModule.Ancestry2):

    moniker = 'date'
    keyName = 'dateType'

    def __init__(self, start, dateType):

        ancestryModule.Ancestry2.__init__(self)

        self.__start = dateModule.raiseIfNotDate(start)
        self.__dateType = abstractClassesModule.raiseIfNotInList(dateType, DateType.allowed, 'dateType')

    @property
    def start(self):

        return self.__start

    value = start

    @property
    def dateType(self):

        return self.__dateType

    def XML_extraAttributes( self, **kwargs ) :

        return ' %s dateType="%s"' % ( self.start.asXML_attribute(name = 'value'), self.dateType )

    def toXMLList(self, **kwargs):

        indent = kwargs.get('indent', '')

        return [ '%s<%s%s dateType="%s"/>' % ( indent, self.moniker, self.start.asXML_attribute(name = 'value'), self.dateType ) ]

    @staticmethod
    def parseConstructBareNodeInstance( node, xPath, linkData, **kwargs ) :

        start = dateModule.Date.parse(node.get('value'))
        dateType = node.get('dateType')

        return Date(start, dateType)

class Dates(suiteModule.Suite):

    moniker = 'dates'

    def __init__(self):

        suiteModule.Suite.__init__(self, [ Date ])

    def toXML(self, indent = '', **kwargs):

        return '\n'.join(self.toXMLList(**kwargs))
