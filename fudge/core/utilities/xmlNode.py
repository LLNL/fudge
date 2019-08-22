# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
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
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
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
        """Keep track of which parser is used ('etree','dom', etc). """
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
            return [xmlNode(a,self.parser) for a in self.data.getchildren()]

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



