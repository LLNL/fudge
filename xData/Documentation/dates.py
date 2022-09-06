# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the classes representing the GNDS documentation nodes author and authors.
"""

from LUPY import enums as enumsModule
from LUPY import ancestry as ancestryModule

from .. import suite as suiteModule
from .. import date as dateModule
from . import abstractClasses as abstractClassesModule

class DateType(enumsModule.Enum):

    accepted = enumsModule.auto()
    available = enumsModule.auto()
    copyrighted = enumsModule.auto()
    collected = enumsModule.auto()
    created = enumsModule.auto()
    issued = enumsModule.auto()
    submitted = enumsModule.auto()
    updated = enumsModule.auto()
    valid = enumsModule.auto()
    withdraw = enumsModule.auto()

class Date(ancestryModule.AncestryIO):

    moniker = 'date'
    keyName = 'dateType'

    def __init__(self, start, dateType):

        ancestryModule.AncestryIO.__init__(self)

        self.__start = dateModule.raiseIfNotDate(start)
        self.__dateType = DateType.checkEnumOrString(dateType)

    @property
    def start(self):

        return self.__start

    value = start

    @property
    def dateType(self):

        return self.__dateType

    def XML_extraAttributes( self, **kwargs ) :

        return ' %s dateType="%s"' % ( self.start.asXML_attribute(name = 'value'), self.dateType )

    def toXML_strList(self, **kwargs):

        indent = kwargs.get('indent', '')

        return [ '%s<%s%s dateType="%s"/>' % ( indent, self.moniker, self.start.asXML_attribute(name = 'value'), self.dateType ) ]

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        start = dateModule.Date.parse(node.get('value'))        # FIXME2, what is this? How is Date to have a date?
        dateType = node.get('dateType')

        return cls(start, dateType)

class Dates(suiteModule.Suite):

    moniker = 'dates'
    suiteName = 'dateType'

    def __init__(self):

        suiteModule.Suite.__init__(self, [ Date ])

    def toXML(self, indent = '', **kwargs):

        return '\n'.join(self.toXML_strList(**kwargs))
