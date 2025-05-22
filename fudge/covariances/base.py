# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Abstract base class for 'standard' covariance matrices (not including parameter covariances).
All covariance classes need a label, ability to get row/column bounds, ability to convert to simple
numeric matrix, etc.
"""

import abc


class Covariance(abc.ABC):

    @property
    @abc.abstractmethod
    def label(self): pass

    @abc.abstractmethod
    def rowBounds(self, unit=None): pass

    @abc.abstractmethod
    def columnBounds(self, unit=None): pass

    @property
    @abc.abstractmethod
    def domainUnit(self): pass

    @property
    @abc.abstractmethod
    def isSymmetric(self): pass

    @abc.abstractmethod
    def toCovarianceMatrix(self, domain=None): pass

    @abc.abstractmethod
    def getUncertaintyVector(self, theData=None, relative=True): pass
