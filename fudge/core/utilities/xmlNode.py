# <<BEGIN-copyright>>
# <<END-copyright>>

"""
Wrapper for xml parsers (Etree, DOM, etc), in case we need to support multiple parsers. For now only supporting python's xml.etree.
This is only for reading (not writing) xml; creating new elements should happen within fudge
"""

__metaclass__ = type

class xmlNode:
    etree = 'etree'
    dom = 'dom'
    allowedParsers = (etree, dom)

    def __init__(self, parsedXMLNode, parser):
        """ keep track of which parser is used ('etree','dom', etc) """
        if parser not in xmlNode.allowedParsers:
            raise Exception("%s xml parser not supported" % parser)
        self.data = parsedXMLNode
        self.parser = parser
        if self.parser==xmlNode.etree:
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
        if self.parser==xmlNode.etree:
            if type(index) is slice:
                # can't slice an iterator, must call getchildren() first:
                return [xmlNode(a, self.parser) 
                        for a in self.data.getchildren()[ index ] ]
            return xmlNode( self.data[index], self.parser )

    def get(self, attribute, defaultValue=None):
        """Returns the value of the specified attribute. If the attribute is not present, returns defaultValue."""
        if self.parser==xmlNode.etree:
            return self.data.get(attribute, defaultValue)

    def keys(self):
        """Returns the name of all attributes in this element, as a list."""
        if self.parser==xmlNode.etree:
            return self.data.keys()

    def items(self):
        """Returns the name and value of all attributes in this element, as a list of tuples."""
        if self.parser==xmlNode.etree:
            return self.data.items()

    def getchildren(self):
        """Returns a list of all child elements."""
        if self.parser==xmlNode.etree:
            return [xmlNode(a) for a in self.data.getchildren()]

    def find(self, path):
        """Searches for child elements matching the path, and returns the first match."""
        if self.parser==xmlNode.etree:
            findResult = self.data.find(path)
            if findResult is None: return
            return xmlNode( findResult, self.parser )

    def findall(self, path):
        """Searches for child elements matching the path, and returns all matches."""
        if self.parser==xmlNode.etree:
            return [xmlNode(a,self.parser) for a in self.data.findall(path)]

    def xpath(self, xpath):
        """Searches for child elements matching the xpath, and returns all matches."""
        if self.parser==xmlNode.etree:
            # etree doesn't currently have native support for xpath.
            # Make my own, using regex for xpath expressions like "reaction[@label='2']":
            import itertools
            import re
            regex = re.compile("([a-zA-Z]+)\[@([a-zA-Z0-9_]+)=['\"]([a-zA-Z0-9_]+|[a-zA-Z]*\([a-zA-Z0-9_,]+\))['\"]\]")

            def xpath2( nextInPath, node ):
                match = regex.match( nextInPath )
                if match:
                    element, attr, val = regex.match( nextInPath ).groups()
                    children = [v for v in node.findall(element) if v.get(attr)==val]
                else:
                    try: children = node.findall( nextInPath )
                    except SyntaxError:
                        raise SyntaxError("%s is not a valid xPath expression!" % nextInPath)
                return children

            xPathList = xpath.split('/')
            results = [self]
            for element in xPathList:
                results = list( itertools.chain(*[xpath2( element, res ) for res in results]) )
            return results



