# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import math
import unittest
import json
from pqu import PQU


class TestCaseWithIsClose(unittest.TestCase):

    def assertIsClose(self, a, b, percent=1.0, absTol=1e-14):
        """
        Uses isclose() function to do a fuzzy equality test of a & b

        :param a: first test value
        :param b: second test value
        :param percent: acceptable percent relative difference between a & b
        :param absTol:  acceptable absolute difference between a & b
        :return: None
        """
        self.assertTrue(math.isclose(a, b, rel_tol=percent / 100., abs_tol=absTol),
                        'a and b differ by %s' % str(a - b))

    def assertListAllAreClose(self, a, b, percent=1.0, absTol=1e-14):
        """
        Element-wise fuzzy equality test of a & b

        :param a: list of first test values
        :param b: list of second test values
        :param percent: acceptable percent relative difference between a & b values
        :param absTol:  acceptable absolute difference between a & b values
        :return: None
        """
        if len(a) != len(b):
            return False
        self.assertTrue(
            all([math.isclose(*x, rel_tol=percent / 100., abs_tol=absTol) for x in zip(a, b)]),
            'Items in list not equal')

    def assertDictAllAreClose(self, a, b, percent=1.0, absTol=1e-14, tryPQUonStrings=True, ignore_keys=None):
        """
        Value-wise fuzzy equality test of a & b.  We do some recursive magic to handle dict nesting.

        :param a: dict of first test values
        :param b: dict of first test values
        :param percent: acceptable percent relative difference between a & b values
        :param absTol:  acceptable absolute difference between a & b values
        :param tryPQUonStrings: if we encounter strings while traversing the tree, attempt to convert them to PQU's
                                for fuzzy equality testing
        :return: None
        """
        # Take care of keys we are ignoring
        _ignore_keys=[]
        if ignore_keys is not None:
            _ignore_keys=ignore_keys

        # Continue with diffs
        for k in a:
            if k in _ignore_keys:
                continue
            self.assertTrue(k in b, "Key in a '%s' not in argument b" % str(k))
        for k in b:
            if k in _ignore_keys:
                continue
            self.assertTrue(k in a, "Key in b '%s' not in argument a" % str(k))
        for k in a:
            if k in _ignore_keys:
                continue
            if isinstance(a[k], float) and isinstance(b[k], float):
                self.assertIsClose(a[k], b[k], percent=percent, absTol=absTol)
            elif isinstance(a[k], list) and isinstance(b[k], list):
                self.assertListAllAreClose(a[k], b[k], percent=percent, absTol=absTol)
            elif isinstance(a[k], tuple) and isinstance(b[k], tuple):
                self.assertListAllAreClose(a[k], b[k], percent=percent, absTol=absTol)
            elif isinstance(a[k], dict) and isinstance(b[k], dict):
                self.assertDictAllAreClose(a[k], b[k], percent=percent, absTol=absTol, tryPQUonStrings=tryPQUonStrings)
            elif isinstance(a[k], str) and isinstance(b[k], str) and tryPQUonStrings:
                try:
                    pak = PQU.PQU(a[k])
                    pbk = PQU.PQU(b[k])
                    self.assertIsClose(pak.getValue(), pbk.getValueAs(pak.unit), percent=percent, absTol=absTol)
                except TypeError:
                    self.assertEqual(a[k], b[k])
            else:
                self.assertEqual(a[k], b[k])

    def assertJSONEqual(self, a, b, ignore_keys=None):
        """
        Compares two strings that are JSON dumps of dictionaries by converting them back to dicts

        :param a: dict of first test values, in a JSON formatted string
        :param b: dict of first test values, in a JSON formatted string
        :return: None
        """
        aDict=json.loads(a)
        bDict=json.loads(b)
        # Continue with diffs
        self.assertDictEqual(aDict, bDict, ignore_keys=ignore_keys)

    def assertDictEqual(self, a, b, ignore_keys=None):
        if ignore_keys is not None:
            for k in ignore_keys:
                if k in a:
                    del a[k]
                if k in b:
                    del b[k]
        # Continue with diffs
        unittest.TestCase.assertDictEqual(self, a, b)


    def assertJSONAllAreClose(self, a, b, percent=1.0, absTol=1e-14, tryPQUonStrings=True, ignore_keys=None):
        """
        Fuzzy comparison of two strings that are JSON dumps of dictionaries by converting them back to dicts

        :param a: dict of first test values, in a JSON formatted string
        :param b: dict of first test values, in a JSON formatted string
        :param percent: acceptable percent relative difference between a & b values
        :param absTol:  acceptable absolute difference between a & b values
        :param tryPQUonStrings: if we encounter strings while traversing the tree, attempt to convert them to PQU's
                                for fuzzy equality testing
        :return: None
        """
        self.assertDictAllAreClose(json.loads(a), json.loads(b), percent=percent, absTol=absTol,
                                   tryPQUonStrings=tryPQUonStrings, ignore_keys=ignore_keys)

    def assertStringsEqual(self, a, b):
        """
        String equality assertion that should work for either Python2 or 3

        :param a: first test value
        :param b: second test value
        :return: None
        """
        listAssertEqual = getattr(self, 'assertCountEqual', None)
        if listAssertEqual is None:
            listAssertEqual = getattr(self, 'assertItemsEqual')

        aList = [x.rstrip() for x in a.split('\n')]
        bList = [x.rstrip() for x in b.split('\n')]
        listAssertEqual(aList, bList)
