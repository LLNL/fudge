"""
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

"""
testPhysicalQuantityWithUncertainty
N.R. Patel, infinidhi@llnl.gov, 7/31/2012
"""
import sys
import math
import unittest
from pqu.physicalQuantityWithUncertainty import PhysicalQuantityWithUncertainty as PQU

class testPhysicalQuantities(unittest.TestCase):
	def setUp(self):
		self.PQU = PQU

	def testInit(self):
		works = 'pass'
		re = 'raise exception'
		xs = [
			['', re],
			["'0.10'", works],
			["'-.12'", works],
			["'+.13'", works],
			["'1'", works],
			["'10'", works],
			["'10.'", works],
			["'1.11'", works],
			["'-1.12'", works],
			["'+1.13'", works],
			["'1.21e-2'", works],
			["'-1.22e2'", works],
			["'-1.23e-2'", works],
			["'-1.24e2 m**2 * eV'", works],
			["'1.31+/-0.31'", works],
			["'1.32e2+/-1.32%'", works],
			["'1.33e3 +/- 1.33e-1% MeV*b'", works],
			["'  1.41   +/-    1.41e0 '", works],
			["'  1.425   +/-    1.420  J   '", works],
			["'  1.51(5)  m **2  *eV'", works],
			["'1.52(5)e-3  kg'", works],
			["'1.61 kg', 'g'", re],
			["'1.62 K', 'K'", works],
			["1.63, 'W'", works],
			["1.64", works],
			["1.71, 'V', 0.11, 'absolute'", works],
			["1.72e3, 'mV', '0.12', 'absolute'", works],
			["1.73, 'kV', '0.13%'", re],		# FIXME doesn't raise... but doesn't work properly
			["1.81, 'A', significantDigits = 3", works],
			["'1.85','kb'", works],
			["1.91, 'Ang', '0.07', significantDigits = 4", re],
			["1.92, 'Ang', '0.7%', 'absolute'", re],
			["'2.0', 'Ang', '0.0'", re],
			["'2.0', 'Ang', '0.0003'", re],		# FIXME why should this raise?
			["'0.0', 'Ang', '0.0'", re],
			["0.0, 'Ang', '0.0'", re],
			["'0.0', 'Ang', 0.0", re],
			["'0.0 +/- 0.0 Ang'", re],
			["'0.0', 'Ang', '0.04'", works],
			["0.0, 'Ang', '0.045'", works],
			["'0.0', 'Ang', 0.050", works],
			["'0.0 +/- 0.06 Ang'", works],
			["'0 +/- 0.04 Ang'", works],
			["'0.0 +/- 0.050 Ang'", works],
			["'0.0 +/- 0.0555 Ang'", works]
		]

		for i, x in enumerate(xs):
			if x[1] == re:
				exec('self.assertRaises(TypeError, self.PQU, %s)' % x[0])
			else:
				xx = eval('self.PQU(%s)' % x[0])

	def testStrReprHash(self):
		xs = [self.PQU('-10.05e3'), self.PQU('-10.05e3 J')]
		rstr = ['-10050', '-10050 J']
		rrepr = ["PhysicalQuantityWithUncertainty(-10050.0, '')", "PhysicalQuantityWithUncertainty(-10050.0, 'J')"]
		rhash = [-7052590062163552704, 8978721800076804947]

		for i, x in enumerate(xs):
			self.assertEqual(str(x), rstr[i])
			self.assertEqual(repr(x), rrepr[i])
			self.assertEqual(hash(x), rhash[i])

	def testArithmetic(self):
		re = 'raise exception'
		xs = [self.PQU(10, ''), self.PQU('10 MeV')]
		ys = [100, self.PQU('100'), self.PQU('100 keV'), self.PQU('100 mb')]
		radd = [[self.PQU('110.0'), self.PQU('110.0'), re, re], [re, re, self.PQU('10.1 MeV'), re]]
		rsub = [[self.PQU('-90.0'), self.PQU('-90.0'), re, re], [re, re, self.PQU('9.9 MeV'), re]]
		rmul = [[self.PQU('1000'), self.PQU('1e3'), self.PQU('1 MeV'), self.PQU('1 b')], [self.PQU('1 GeV'), self.PQU('1 GeV'), self.PQU('1 MeV**2'), self.PQU('1 MeV*b')]]
		rdiv = [[self.PQU('0.1'), self.PQU('1e-1'), self.PQU('0.1e-3 1/eV'), self.PQU('0.1 1/mb')], [self.PQU('0.1 MeV'), self.PQU('100 keV'), self.PQU('100'), self.PQU('0.1 MeV/mb')]]
		rs = [radd, rsub, rmul, rdiv]
		zs = ['+', '-', '*', '/']
		ss = ['__add__', '__sub__', '__mul__', '__div__']

		for i, x in enumerate(xs):
			for j, y in enumerate(ys):
				for k, z in enumerate(zs): 
					if rs[k][i][j] == re:
						exec('self.assertRaises(TypeError, x.%s, y)' % ss[k])
					else:
						result = eval('x %s y' % z)
						self.assertEqual(result, rs[k][i][j])

	def testPower(self):
		re = 'raise exception'
		xs = [self.PQU('5 MeV'), self.PQU('5')]
		ys = [2, -2, 0.5, self.PQU('2'), self.PQU('2 s')]
		rs = [[self.PQU('25 MeV**2'), self.PQU('0.04 1/MeV**2'), re, re, re], [25, 0.04, math.sqrt(5), re, re]]

		for i, x in enumerate(xs):
			for j, y in enumerate(ys):
				if rs[i][j] == re:
					self.assertRaises(TypeError, x.__pow__, y)
				else:
					z = x**y
					self.assertEqual(z, rs[i][j])

	def testAbsolute(self):
		xs = [self.PQU('1.5'), self.PQU('-1.5'), self.PQU('-1.5 eV'), self.PQU('1.5 eV')]
		rs = [1.5, self.PQU('1.5'), self.PQU('1.5 eV'), self.PQU('1.5 eV')]

		for i,x in enumerate(xs):
			result = abs(x)
			self.assertEqual(result, rs[i])

	def testComparisonOperations(self):
		re = 'raise exception'
		xs = [self.PQU('10'), self.PQU('10 eV * Mb')]
		ys = [10, 10.000000000000001, 3.14, self.PQU('10.0'), self.PQU('10 MeV * b'), self.PQU('30 eV * Mb')]
		rsEq = [[True, True, False, True, False, False], [False, False, False, False, True, False]]
		rsNe = [[False, False, True, False, True, True], [True, True, True, True, False, True]]
		rsLe = [[True, True, False, True, re, re], [re, re, re, re, True, True]]
		rsGe = [[True, True, True, True, re, re], [re, re, re, re, True, False]]
		rsLt = [[False, False, False, False, re, re], [re, re, re, re, False, True]]
		rsGt = [[False, False, True, False, re, re], [re, re, re, re, False, False]]
		rs = [rsEq, rsNe, rsLe, rsGe, rsLt, rsGt]
		zs = ['==', '!=', '<=', '>=', '<', '>']
		ss = ['__eq__', '__ne__', '__le__', '__ge__', '__lt__', '__gt__']

		for i, pqu_ in enumerate(xs):
			for j, value in enumerate(ys):
				for k, operation in enumerate(zs):
					if rs[k][i][j] == re:
						exec('self.assertRaises(TypeError, pqu_.%s, value)' % ss[k])
					else:
						result = eval('pqu_ %s value' % operation)
						self.assertEqual(result, rs[k][i][j])

	def testConversions(self):
		re = 'raise exception'
		xs = [self.PQU('2.718'), self.PQU('-40 degC'), self.PQU('3.14 W')]
		ys = ['', 'K', 'degF', 'J/s', 'MeV/s', 'A*V']
		rGetValueAs = [[2.718, re, re, re, re, re], [re, 233.14999999999998, -40.000000000000057, re, re, re],
				[re, re, re, 3.14, 19598329980115.25, 3.14]]
		rInUnitsOf = [[self.PQU('2.718'), re, re, re, re, re],
				[re, self.PQU('233.14999999999998 K'), self.PQU('-40.000000000000057 degF'), re, re, re],
				[re, re, re, self.PQU('3.14 J/s'), self.PQU('19598329980115.25 MeV/s'), self.PQU('3.14 A*V')]]
		rIsCompatible = [[True, False, False, False, False, False], [False, True, True, False, False, False],
				[False, False, False, True, True, True]]
		rs = [rGetValueAs, rInUnitsOf, rIsCompatible]
		zs = ['getValueAs', 'inUnitsOf', 'isCompatible']

		for i, pqu_ in enumerate(xs):
			for j, toUnit in enumerate(ys):
				for k, method in enumerate(zs): 
					if rs[k][i][j] == re:
						exec('self.assertRaises(TypeError, pqu_.%s, toUnit)' % method)
					else:
						result = eval('pqu_.%s(toUnit)' % method)
						self.assertEqual(result, rs[k][i][j])

	def testConvertToUnit(self):
		re = 'raise exception'
		xs = [self.PQU('2.718'), self.PQU('-40 degC'), self.PQU('3.14 W')]
		ys = ['', 'K', 'degF', 'J/s', 'MeV/s', 'A*V']
		rs = [[self.PQU('2.718'), re, re, re, re, re],
				[re, self.PQU('233.14999999999998 K'), self.PQU('-40.000000000000057 degF'), re, re, re],
				[re, re, re, self.PQU('3.14 J/s'), self.PQU('19598329980115.25 MeV/s'), self.PQU('3.14 A*V')]]

		for i, x in enumerate(xs):
			for j, y in enumerate(ys):
				if rs[i][j] == re:
					self.assertRaises(TypeError, x.convertToUnit, y)

				else:
					x.convertToUnit(y)
					self.assertEqual(x, rs[i][j])

	def testInBaseUnits(self):
		re = 'raise exception'
		xs = [self.PQU('2.718'), self.PQU('-40 degC'), self.PQU('3.14 W')]
		rs = [re, self.PQU('233.14999999999998 K'), self.PQU('3.14 kg*m**2/s**3')]

		# FIXME: should inBaseUnits() raise TypeError for dimensionless numbers?

		for i, x in enumerate(xs):
			if rs[i] == re:
				self.assertRaises(TypeError, x.inBaseUnits)
			else:
				result = x.inBaseUnits()
				self.assertEqual(result, rs[i])

	def testTrigFunctions(self):
		xs = [self.PQU('30 deg'), self.PQU(math.pi/4, 'rad')]
		rs = [[0.5, math.sqrt(3 / 4.), 1 / math.sqrt(3.)], [1 / math.sqrt(2.), 1 / math.sqrt(2.), 1]]
		zs = ['sin', 'cos', 'tan']

		for i, x in enumerate(xs):
			for k, z in enumerate(zs): 
				result = eval('x.%s()' % z)
				self.assertAlmostEqual(result, rs[i][k])

if __name__ == '__main__':
	unittest.main()

