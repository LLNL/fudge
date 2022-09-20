# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Wrapper for xml parsers (Etree, DOM, etc), in case we need to support multiple parsers. For now only supporting
python's xml.etree.  This is only for reading (not writing) xml; creating new elements should happen within fudge
"""


class XML_node:

    etree = 'etree'
    dom = 'dom'
    allowedParsers = (etree, dom)

    def __init__(self, parsedXMLNode, parser):
        """Keep track of which parser is used ('etree','dom', etc). """

        if parser not in XML_node.allowedParsers:
            raise Exception("%s xml parser not supported" % parser)
        self.data = parsedXMLNode
        self.parser = parser
        if self.parser == XML_node.etree:
            self.text = self.data.text
            self.tail = self.data.tail
            self.tag = self.data.tag

    def __len__(self):
        """Returns the number of child elements."""
        return len(self.data)

    def __repr__(self):
        return "<Element '%s' at 0x%x>" % (self.tag, id(self))

    def __getitem__(self, index):
        """Return one or more child elements. The index could be a slice."""
        if self.parser == XML_node.etree:
            if isinstance(index, slice):
                # can't slice an iterator, must call getchildren() first:
                return [XML_node(a, self.parser) for a in list(self.data)[index]]
            return XML_node(self.data[index], self.parser)

    @property
    def attrib(self):
        """Access node attributes as dictionary"""
        if self.parser == XML_node.etree:
            return self.data.attrib

    def get(self, attribute, defaultValue=None):
        """Returns the value of the specified attribute. If the attribute is not present, returns defaultValue."""
        if self.parser == XML_node.etree:
            return self.data.get(attribute, defaultValue)

    def keys(self):
        """Returns the name of all attributes in this element, as a list."""
        if self.parser == XML_node.etree:
            return list(self.data.keys())

    def items(self):
        """Returns the name and value of all attributes in this element, as a list of tuples."""
        if self.parser == XML_node.etree:
            return list(self.data.items())

    def getchildren(self):
        """Returns a list of all child elements."""
        if self.parser == XML_node.etree:
            return [XML_node(a, self.parser) for a in self.data.getchildren()]

    def find(self, path):
        """Searches for child elements matching the path, and returns the first match."""
        if self.parser == XML_node.etree:
            findResult = self.data.find(path)
            if findResult is None:
                return
            return XML_node(findResult, self.parser)

    def findall(self, path):
        """Searches for child elements matching the path, and returns all matches."""
        if self.parser == XML_node.etree:
            return [XML_node(a, self.parser) for a in self.data.findall(path)]

    def xpath(self, xpath):
        """Searches for child elements matching the xpath, and returns all matches."""
        if self.parser == XML_node.etree:
            # etree doesn't currently have native support for xpath.
            # Make my own, using regex for xpath expressions like "reaction[@label='2']":
            import re
            import itertools
            regex = re.compile("""([a-zA-Z]+)[@([a-zA-Z0-9_]+)=['\"]([a-zA-Z0-9_]+|[a-zA-Z]*([a-zA-Z0-9_,]+))['\"]]""")

            def xpath2(nextInPath, node):
                match = regex.match(nextInPath)
                if match:
                    _element, attr, val = regex.match(nextInPath).groups()
                    children = [v for v in node.findall(_element) if v.get(attr) == val]
                else:
                    try:
                        children = node.findall(nextInPath)
                    except SyntaxError:
                        raise SyntaxError("%s is not a valid xPath expression!" % nextInPath)
                return children

            xPathList = xpath.split('/')
            results = [self]
            for element in xPathList:
                results = list(itertools.chain(*[xpath2(element, res) for res in results]))
            return results
