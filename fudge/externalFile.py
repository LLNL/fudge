# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the externalFile class, used to keep track of
connections between GNDS files.
"""

__metaclass__ = type

import os

import xData.ancestry as ancestryModule
from LUPY import checksums

class externalFile( ancestryModule.ancestry ):
    """
    An externalFile contains a unique label,
    and the relative path to another file (unix path representation).
    """

    moniker = 'externalFile'

    def __init__(self, label, path, checksum=None, algorithm=checksums.sha1sum.algorithm):

        ancestryModule.ancestry.__init__(self)
        self.__label = label
        self.path = path
        self.__checksum = checksum
        self.__algorithm = algorithm

    @property
    def label(self):
        return self.__label

    @property
    def path(self):
        return self.__path

    @path.setter
    def path(self, value):
        if not isinstance(value,str):
            raise TypeError("Path must be a string!")
        self.__path = value

    @property
    def checksum(self):
        return self.__checksum

    @property
    def algorithm(self):
        return self.__algorithm

    def realpath( self ) :
        """Returns the realpath to the external file."""

        path = self.path
        if( path[0] != os.sep ) :
            path = os.path.join( os.path.dirname( self.rootAncestor.sourcePath ), path )

        return( os.path.realpath( path ) )


    def toXMLList(self, indent = '', **kwargs) :

        attrs = ""
        if self.checksum:
            attrs = ' checksum="%s" algorithm="%s"' % (self.checksum, self.algorithm)
        return [ '%s<%s label="%s" path="%s"%s/>' % ( indent, self.moniker, self.label, self.path, attrs ) ]

    @classmethod
    def parseXMLNode(cls, element, xPath, linkData):

        xPath.append('%s[@label="%s"]' % (element.tag, element.get('label')))
        ef = cls( element.get('label'), element.get('path'), element.get('checksum'), element.get('algorithm') )
        xPath.pop()
        return ef
