# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
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

from xml.sax import saxutils

from LUPY import ancestry as ancestryModule

from xData import link as linkModule


class ENDFconversionFlags( ancestryModule.AncestryIO ) :

    moniker = "ENDFconversionFlags"

    def __init__(self):

        ancestryModule.AncestryIO.__init__( self )
        self.flags = []

    def add(self, item, flag):
        flag = Conversion(item, flags=flag)
        self.flags.append(flag)

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent','  ')
        xml = ['%s<%s>' % (indent, self.moniker)]
        for link_ in self.flags:
            xml += link_.toXML_strList(indent2, **kwargs)
        xml[-1] += '</%s>' % self.moniker
        return xml

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        xPath.append(node.tag)

        flags_ = cls()
        for child in node:
            flags_.flags.append(Conversion.parseNodeUsingClass(child, xPath, linkData, **kwargs))

        xPath.pop()
        return flags_

class Conversion(linkModule.Link):

    moniker = "conversion"

    def __init__( self, link=None, root=None, path=None, label=None, relative=False,
                  flags=None ):
        linkModule.Link.__init__(self, link=link, root=root, path=path, label=label, relative=relative)
        self.flags = flags

    def __deepcopy__(self, memo = {}):

        if self.path is None: self.updateXPath()
        return self.__class__(link = self.link, root = self.root, path = self.path, label = self.label, relative = self.relative, flags = self.flags)

    def toXML_strList( self, indent = '', **kwargs ) :

        attributesStr = ""
        for attr in ('flags',):
            if getattr(self, attr):
                attributesStr += ' %s="%s"' % (attr, getattr(self, attr))

        return ['%s<%s%s href="%s" />' % (indent, self.moniker, attributesStr, self.build_href(**kwargs))]
