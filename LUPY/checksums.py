# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Classes and static methods for computing checksums
"""

import hashlib

blocksize = 2**16

class NoneSum:
    """
    Dummy class to handle *None* as an algorithm.
    """

    algorithm = None

    @staticmethod
    def from_file(filename : str):
        """
        Always returns *None*.
        """

        return None

    @staticmethod
    def from_string(string : str):
        """
        Always returns *None*.
        """

        return None

class sha1sum:
    algorithm = "sha1"

    @staticmethod
    def from_file(filename : str):
        """
        Return SHA1 checksum for a file.
        @param filename:
        @return:
        """
        checksum = hashlib.sha1()
        with open(filename, 'rb') as source:
            block = source.read(blocksize)
            while len(block) != 0:
                checksum.update(block)
                block = source.read(blocksize)

        return str(checksum.hexdigest())

    @staticmethod
    def from_string(string : str):
        """
        Return SHA1 checksum for a string.
        @param string:
        @return:
        """
        checksum = hashlib.sha1()
        data = string.encode("utf8")
        block = data[:blocksize]
        while len(block) != 0:
            checksum.update(block)
            data = data[blocksize:]
            block = data[:blocksize]

        return str(checksum.hexdigest())


class md5sum:
    algorithm = "md5"

    @staticmethod
    def from_file(filename: str):
        """
        Return MD5 checksum for a file.
        @param filename:
        @return:
        """
        checksum = hashlib.md5()
        with open(filename, 'rb') as source:
            block = source.read(blocksize)
            while len(block) != 0:
                checksum.update(block)
                block = source.read(blocksize)

        return str(checksum.hexdigest())

    @staticmethod
    def from_string(string: str):
        """
        Return MD5 checksum for a string.
        @param string:
        @return:
        """
        checksum = hashlib.md5()
        data = string.encode("utf8")
        block = data[:blocksize]
        while len(block) != 0:
            checksum.update(block)
            data = data[blocksize:]
            block = data[:blocksize]

        return str(checksum.hexdigest())


checkers = {
    NoneSum.algorithm: NoneSum,
    sha1sum.algorithm: sha1sum,
    md5sum.algorithm: md5sum
}
supportedAlgorithms = checkers.keys()
