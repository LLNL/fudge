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

"""
Utility to handle linking in gnd.
xml representation of the link is implemented as xpath/xpointer

The 'link' class contains a link to a different location in the 
After reading in a link class instance from xml, use 'follow' to get the pointed-to element

sample links::
<element xlink:type="simple" xlink:href='/reactionSuite/reaction[@label="102"]'>
or
<element xlink:type="simple" xlink:href='/covariances.xml#covId'>

cmattoon, 2/4/2011
"""

__metaclass__ = type

def parseXMLNode( linkElement, linkData={} ):
    """ parse the xml-represented link back to python. The resulting link points to None,
    and must be resolved by the calling function """
    xlink_namespace = '{http://www.w3.org/1999/xlink}'
    label = linkElement.tag
    path = linkElement.get(xlink_namespace+'href')
    if '#' in path: root, path = path.split('#')
    else: root = None
    # all optional (non-xlink) attributes:
    attributes = dict( [v for v in linkElement.items() if xlink_namespace not in v[0]] )
    return link( label=label, root=root, path=path, **attributes )

class link:
    """
    the 'link' class contains the path to another element, that could be stored in a different file.
    The 'resolve' member function is used to get a pointer to the linked-to element,
    which will then be stored as self.link.
    
    Member data ::
    
        * label : the xml tag name (type: str)
        * link  : actual pointer to other data (type: str)
        * root  : file name where other data is stored, only needed if it's not the current file (type: str)
        * path  : self.path is an xPath expression that uniquely identifies the other data (type: str)
    
    """
    def __init__(self, label='link', link=None, root=None, path=None, **attributes):
        self.label = label  # becomes the xml tag name
        self.link = link    # actual pointer to other data
        self.root = root    # file name where other data is stored, only needed if it's not the current file
        self.path = path    # self.path is an xPath expression that uniquely identifies the other data
        # careful with nesting single/double quotes:
        if self.path is not None: self.path = self.path.replace('"',"'")
        self.attributes = attributes

    def __str__( self ) :

        return( "%s" % self.path )

    def __getitem__( self, key ):
        return self.attributes[key]

    def __setitem__( self, key, value ):
        self.attributes[key] = value

    def updateXPath( self ):
        '''ensure the xPath agrees with the linked-to data'''
        if hasattr(self.link, 'toXLink'): self.path = self.link.toXLink()

    def toXML(self, indent=''):
        """ pointers show up in the attributes list on an xml element
        ie, <element xlink:type="simple" ...
        """
        self.updateXPath()
        extraInfo = ''
        for key in self.attributes:
            extraInfo += ' %s="%s"' % (key, self.attributes[key])
        if type(self.root)==str:    # external link
            return indent+'<%s xlink:type="simple" xlink:href="%s#%s"%s/>' % (
                    self.label, self.root, self.path, extraInfo)
        else:                       # link within this file
            return indent+'<%s xlink:type="simple" xlink:href="%s"%s/>' % (
                    self.label, self.path, extraInfo)

def follow(xPath, node):
    """
    find the fudge data member pointed to by an xPath expression
    the 'node' argument should generally be either a reactionSuite or covarianceSuite instance
    """
    import re
    # use regex on xPath expressions like "reaction[@label='2']"
    regex = re.compile("([a-zA-Z]+)\[@([a-z]+)='([a-zA-Z0-9_]+|[a-zA-Z]*\([a-zA-Z0-9_,]+\))'\]")

    def follow2( xPathList, node ):
        # recursive helper function: descend the path to find the correct element
        if len(xPathList)==0: return node

        xPathNext = xPathList[0]
        try:
            nodeNext = getattr(node, xPathNext)
        except AttributeError:
            match = regex.match(xPathNext)
            if match:
                attribute, label, ID = match.groups()
                if attribute in ('summedReaction','fissionComponent','production','reactionSum',
                        'externalReaction'):
                    nodeNext = [e for e in getattr(node, attribute+'s') if getattr(e,label)==ID]
                else:
                    nodeNext = [e for e in node if getattr(e,label)==ID]
                assert len(nodeNext)==1
                nodeNext = nodeNext[0]

            # some special cases:
            elif xPathNext in ('','reactionSuite','covarianceSuite'):
                nodeNext = node
            elif xPathNext=='crossSection':
                nodeNext = getattr(node, 'crossSection', None)
                if nodeNext is None: nodeNext = node.data
            elif xPathNext=='energy':
                nodeNext = getattr(node,'energyComponent',None)
                if nodeNext is None: nodeNext = node['energy']
            elif xPathNext=='angular':
                nodeNext = getattr(node,'angularComponent',None)
                if nodeNext is None: nodeNext = node['angular']
            elif xPathNext=='uncorrelated':
                nodeNext = node['uncorrelated']

            else:
                raise Exception("Encountered unknown element in xPath: %s" % xPathNext)

        return follow2(xPathList[1:], nodeNext)

    return follow2( xPath.split('/'), node )


