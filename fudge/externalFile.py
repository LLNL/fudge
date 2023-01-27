# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the externalFile class, used to keep track of
connections between GNDS files.
"""

import os

from LUPY import ancestry as ancestryModule
from LUPY import checksums as checksumsModule

class ExternalFile(ancestryModule.AncestryIO):
    """
    An externalFile contains a unique label,
    and the relative path to another file (unix path representation).
    """

    moniker = 'externalFile'

    def __init__(self, label, path, checksum=None, algorithm=checksumsModule.Sha1sum.algorithm):

        ancestryModule.AncestryIO.__init__(self)
        self.__label = label
        self.path = path
        self.__checksum = checksum
        self.__algorithm = algorithm
        self.__instance = None

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

    @property
    def instance(self):
        return self.__instance

    @instance.setter
    def instance(self, value):
        # For now instance can be set to None or any class instance. Could do more rigorous checking if we add a
        # 'fileType' attribute to the external file telling what it is (reactionSuite, covarianceSuite), etc.
        assert isinstance(value, (object, type(None)))
        self.__instance = value

    def updateChecksum(self, algorithm=checksumsModule.Sha1sum.algorithm):
        """
        Recompute checksum from external file, update checksum and algorithm fields.
        @param algorithm: 'sha1' or 'md5'
        """

        checker = checksumsModule.checkers[algorithm]
        self.__checksum = checker.from_file(self.realpath())
        self.__algorithm = algorithm

    def realpath( self ) :
        """Returns the realpath to the external file."""

        path = self.path
        if( path[0] != os.sep ) :
            path = os.path.join( os.path.dirname( self.rootAncestor.sourcePath ), path )

        return( os.path.realpath( path ) )

    def toXML_strList(self, indent='', **kwargs):

        attrs = ""
        if self.checksum:
            attrs = ' checksum="%s" algorithm="%s"' % (self.checksum, self.algorithm)
        return [ '%s<%s label="%s" path="%s"%s/>' % ( indent, self.moniker, self.label, self.path, attrs ) ]

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        xPath.append('%s[@label="%s"]' % (node.tag, node.get('label')))

        ef = cls(node.get('label'), node.get('path'), node.get('checksum'), node.get('algorithm'))

        xPath.pop()

        return ef
