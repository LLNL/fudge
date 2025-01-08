#! /usr/bin/env python3
# encoding: utf-8

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
test fudge/covariances/
gert1 04/16/2020
"""

import unittest
import os

from fudge import reactionSuite as reactionSuiteModule
from fudge.covariances import covarianceSuite as covarianceSuiteModule
from fudge.covariances import enums as covarianceEnumsModule

TEST_DATA_PATH, this_filename = os.path.split(__file__)
NaEvaluation = reactionSuiteModule.ReactionSuite.readXML_file(
    os.path.join(TEST_DATA_PATH, 'n-011_Na_023.endf.gnds.xml'))
NaCovariance = covarianceSuiteModule.CovarianceSuite.readXML_file(
    os.path.join(TEST_DATA_PATH, 'n-011_Na_023.endf.gndsCov.xml'), reactionSuite=NaEvaluation)


class TesttoCovarianceMatrix(unittest.TestCase):
    def testToCovarianceMatrix(self):
        for section in NaCovariance.covarianceSections:
            self.assertIn(section['eval'].toCovarianceMatrix().type, covarianceEnumsModule.Type)


if __name__ == '__main__':
    unittest.main()
