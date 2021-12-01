# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
"""
Defines a class for storing special instructions for converting from GNDS back to ENDF-6.
These flags are sometimes required since ENDF-6 supports more than one way to store some types of data.
For example, particle distributions could be in MF=6, or in a combination of MF=4 and 5.

Internally this class simply stores a set of links to other data types within a GNDS file, along with
a comma-delimited list of ENDF conversion flags for each data type
"""

from xData import ancestry as ancestryModule
from xData import link as linkModule

__metaclass__ = type

class ENDFconversionFlags( ancestryModule.ancestry ) :

    moniker = "ENDFconversionFlags"

    def __init__(self):

        ancestryModule.ancestry.__init__( self )
        self.flags = []

    def add(self, item, flag):
        flag = conversion(item, flags=flag)
        self.flags.append(flag)

    def toXMLList(self, indent='', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent','  ')
        xml = ['%s<%s>' % (indent, self.moniker)]
        xml.append('%s<!-- flags to help with writing back to ENDF6 when data can be stored more than one way -->' % indent2)
        for link_ in self.flags:
            xml += link_.toXMLList(indent2, **kwargs)
        xml[-1] += '</%s>' % self.moniker
        return xml

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):

        xPath.append(element.tag)
        flags_ = cls()
        for child in element:
            flags_.flags.append( conversion.parseXMLNode(child, xPath, linkData) )
        xPath.pop()
        return flags_

class conversion(linkModule.link):

    moniker = "conversion"

    def __init__( self, link=None, root=None, path=None, label=None, relative=False,
                  flags=None ):
        linkModule.link.__init__(self, link=link, root=root, path=path, label=label, relative=relative)
        self.flags = flags

    def toXMLList( self, indent = '', **kwargs ) :

        attributesStr = ""
        for attr in ('flags',):
            if getattr(self, attr):
                attributesStr += ' %s="%s"' % (attr, getattr(self, attr))

        return [ '%s<%s%s href="%s"/>' % ( indent, self.moniker, attributesStr, self.build_href( **kwargs ) ) ]
