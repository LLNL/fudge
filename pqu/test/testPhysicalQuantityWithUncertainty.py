"""
# <<BEGIN-copyright>>
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

verbose = 2
doTestAll = True
doTestInit = False
doTestStrReprHash = False
doTestArithmetic = False
doTestPower = False
doTestAbsolute = False
doTestComparisonOperations = False
doTestConversions = False
doTestConvertToUnit = False
doTestInBaseUnits = False
doTestTrigFunctions = False
doTest = False
doTest = False
doTest = False
doTest = False

class testDimensionlessPhysicalQuantities(unittest.TestCase):
	def setUp(self):
		self.PQU = PQU

	if doTestInit or doTestAll:
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
["1.73, 'kV', '0.13%'", re],
["1.81, 'A', significantDigits = 3", works],
["'1.85','kb'", works],
["1.91, 'Ang', '0.07', significantDigits = 4", re],
["1.92, 'Ang', '0.7%', 'absolute'", re],
["'2.0', 'Ang', '0.0'", re],
["'2.0', 'Ang', '0.0003'", re],
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

			if verbose > 0:
				print '\ntesting __init__'

			if verbose > 1:
				rows, lenOfRows, columns = [], [], []
				lenOfCol0, lenOfCol1, lenOfCol2 = [len('input') + 1], [len('output') + 1], [len('value') + 1]
				lenOfCol3, lenOfCol4, lenOfCol5 = [len('unit') + 1], [len('unc') + 1], [len('uncType') + 1]
				lenOfCol6 = [len('sigDigs') + 1]

				for i, x in enumerate(xs):
					lenOfCol0.append(len(str(x[0])) + 6)

					if x[1] == re:
						lenOfCol1.append(len(re) + 1)
					else:
						xx = eval('self.PQU(%s)' % x[0])
						lenOfCol1.append(len(str(xx)) + 1)
						lenOfCol2.append(len(str(xx.value)) + 1)
						lenOfCol3.append(len(str(xx.unit)) + 1)

						try:
							lenOfCol4.append(len(str(xx.uncertainty)) + 1)
							lenOfCol5.append(len(str(xx.uncertaintyType)) + 1)
							lenOfCol6.append(len(str(xx.significantDigits)) + 1)
						except:
							pass

				maxLenOfCol0 = max(lenOfCol0)
				maxLenOfCol1 = max(lenOfCol1)
				maxLenOfCol2 = max(lenOfCol2)
				maxLenOfCol3 = max(lenOfCol3)
				maxLenOfCol4 = max(lenOfCol4)
				maxLenOfCol5 = max(lenOfCol5)
				maxLenOfCol6 = max(lenOfCol6) - 1

			for i, x in enumerate(xs):
				if verbose > 1:
					columns.__init__()
					columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol0) + 's')) % ('PQU(' + x[0] + ')'))

					if x[1] == re:
						exec('self.assertRaises(TypeError, self.PQU, %s)' % x[0])

						if verbose > 1:
							columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol1) + 's')) % 're')
					else:
						xx = eval('self.PQU(%s)' % x[0])

						if verbose > 1:
							columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol1) + 's')) % xx)
							columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol2) + 's')) % str(xx.value))
							columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol3) + 's')) % str(xx.unit))

							try:
								columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol4) + 's')) % str(xx.uncertainty))
								columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol5) + 's')) % str(xx.uncertaintyType))
								columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol6) + 's')) % str(xx.significantDigits))
							except:
								pass

					if verbose > 1:
						rows.append(''.join(columns))
						lenOfRows.append(len(''.join(columns)))

			if verbose > 1:
				rowHead = eval("'%s'" % ('%-' + str(maxLenOfCol0) + 's%-' + str(maxLenOfCol1) + 's%-' + str(maxLenOfCol2) + 's%-' + str(maxLenOfCol3) + 's%-' + str(maxLenOfCol4) + 's%-' + str(maxLenOfCol5) + 's%-' + str(maxLenOfCol6) + 's')) % ('input', 'output', 'value', 'unit', 'unc', 'uncType', 'sigDigs')
				lenOfRows.append(len(rowHead))
				rowLine = '%-s' % (max(lenOfRows) * '-')
				rows.insert(0, rowLine)
				rows.insert(0, rowHead)
				rows.append(rowLine)
				print '\n'.join(rows)


	if doTestStrReprHash or doTestAll:
		def testStrReprHash(self):
			xs = [self.PQU('-10.05e3'), self.PQU('-10.05e3 J')]
			rstr = ['-10050', '-10050 J']
			rrepr = ["PhysicalQuantityWithUncertainty(-10050.0, '')", "PhysicalQuantityWithUncertainty(-10050.0, 'J')"]
			rhash = [-7052590062163552704, 8978721800076804947]

			if verbose > 0:
				print '\ntesting str, repr, hash'

			if verbose > 1:
				rows, lenOfRows, columns = [], [], []
				lenOfCol0, lenOfCol1, lenOfCol2 = [], [], []

				for i, x in enumerate(xs):
					lenOfCol0.append(len(str(x)) + 6)
					lenOfCol1.append(len(str(x)) + 1)
					lenOfCol2.append(len(repr(x)) + 1)

				maxLenOfCol0 = max(lenOfCol0)
				maxLenOfCol1 = max(lenOfCol1)
				maxLenOfCol2 = max(lenOfCol2)

			for i, x in enumerate(xs):
				self.assertEqual(str(x), rstr[i])
				self.assertEqual(repr(x), rrepr[i])
				self.assertEqual(hash(x), rhash[i])

				if verbose > 1:
					if isinstance(x, self.PQU):
						columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol0) + 's')) % ('PQU(' + str(x) + ')'))
					else:
						columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol0) + 's')) % str(x))

					columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol1) + 's%-' + str(maxLenOfCol2) + 's%-s')) % (str(x), repr(x), hash(x)))
					rows.append(''.join(columns))
					lenOfRows.append(len(''.join(columns)))
					columns.__init__()

			if verbose > 1:
				rowHead = eval("'%s'" % ('%-' + str(maxLenOfCol0) + 's%-' + str(maxLenOfCol1) + 's%-' + str(maxLenOfCol2) + 's%-s')) % ('x', 'str(x)', 'repr(x)', 'hash(x)')
				lenOfRows.append(len(rowHead))
				rowLine = '%-s' % (max(lenOfRows) * '-')
				rows.insert(0, rowLine)
				rows.insert(0, rowHead)
				rows.insert(0, rowLine)
				rows.append(rowLine)
				print '\n'.join(rows)

	if doTestArithmetic or doTestAll:
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

			if verbose > 0:
				print '\ntesting arithmetic operations'

			if verbose > 1:
				rows, lenOfRows, columns = [], [], []
				lenOfCol0, lenOfCol1, lenOfCol2 = [], [], []
				lenOfCol3, lenOfCol4, lenOfCol5 = [], [], []

				for i, x in enumerate(xs):
					lenOfCol0.append(len(str(x)) + 6)

					for j, y in enumerate(ys):
						lenOfCol1.append(len(str(y)) + 6)

						for k, z in enumerate(zs): 
							if rs[k][i][j] == re:
								eval('lenOfCol%i.append(len(re) + 1)' % (k + 2))
							else:
								result = eval('x %s y' % z)
								eval('lenOfCol%i.append(len(str(result)) + 6)' % (k + 2))

				maxLenOfCol0 = max(lenOfCol0)
				maxLenOfCol1 = max(lenOfCol1)
				maxLenOfCol2 = max(lenOfCol2)
				maxLenOfCol3 = max(lenOfCol3)
				maxLenOfCol4 = max(lenOfCol4)
				maxLenOfCol5 = max(lenOfCol5) - 1

			for i, x in enumerate(xs):
				if verbose > 1:
					columns.__init__()
					if isinstance(x, self.PQU):
						columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol0) + 's')) % ('PQU(' + str(x) + ')'))
					else:
						columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol0) + 's')) % str(x))

				for j, y in enumerate(ys):
					if verbose > 1:
						if isinstance(y, self.PQU):
							columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol1) + 's')) % ('PQU(' + str(y) + ')'))
						else:
							columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol1) + 's')) % str(y))

					for k, z in enumerate(zs): 
						if rs[k][i][j] == re:
							exec('self.assertRaises(TypeError, x.%s, y)' % ss[k])

							if verbose > 1:
								columns.append(eval("'%s'" % ('%-' + str(eval('maxLenOfCol%i' % (k + 2))) + 's')) % re)
						else:
							result = eval('x %s y' % z)
							self.assertEqual(result, rs[k][i][j])

							if verbose > 1:
								if isinstance(result, self.PQU):
									columns.append(eval("'%s'" % ('%-' + str(eval('maxLenOfCol%i' % (k + 2))) + 's')) % ('PQU(' + str(result) + ')'))
								else:
									columns.append(eval("'%s'" % ('%-' + str(eval('maxLenOfCol%i' % (k + 2))) + 's')) % str(result))

					if verbose > 1:
						rows.append(''.join(columns))
						lenOfRows.append(len(''.join(columns)))
						del columns[1:]

			if verbose > 1:
				rowHead = eval("'%s'" % ('%-' + str(maxLenOfCol0) + 's%-' + str(maxLenOfCol1) + 's%-' + str(maxLenOfCol2) + 's%-' + str(maxLenOfCol3) + 's%-' + str(maxLenOfCol4) + 's%-' + str(maxLenOfCol5) + 's')) % ('x', 'y', 'x + y', 'x - y', 'x * y', 'x / y')
				lenOfRows.append(len(rowHead))
				rowLine = '%-s' % (max(lenOfRows) * '-')
				rows.insert(0, rowLine)
				rows.insert(0, rowHead)
				rows.insert(0, rowLine)
				rows.append(rowLine)
				print '\n'.join(rows)

	if doTestPower or doTestAll:
		def testPower(self):
			re = 'raise exception'
			xs = [self.PQU('5 MeV'), self.PQU('5')]
			ys = [2, -2, 0.5, self.PQU('2'), self.PQU('2 s')]
			rs = [[self.PQU('25 MeV**2'), self.PQU('0.04 1/MeV**2'), re, re, re], [25, 0.04, math.sqrt(5), re, re]]

			if verbose > 0:
				print '\ntesting __pow__'

			if verbose > 1:
				rows, lenOfRows, columns = [], [], []
				lenOfCol0, lenOfCol1, lenOfCol2 = [], [], []

				for i, x in enumerate(xs):
					lenOfCol0.append(len(str(x)) + 6)

					for j, y in enumerate(ys):
						lenOfCol1.append(len(str(y)) + 6)

				maxLenOfCol0 = max(lenOfCol0)
				maxLenOfCol1 = max(lenOfCol1)

			for i, x in enumerate(xs):
				if verbose > 1:
					columns.__init__()

					if isinstance(x, self.PQU):
						columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol0) + 's')) % ('PQU(' + str(x) + ')'))
					else:
						columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol0) + 's')) % str(x))

				for j, y in enumerate(ys):
					if verbose > 1:
						if isinstance(y, self.PQU):
							columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol1) + 's')) % ('PQU(' + str(y) + ')'))
						else:
							columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol1) + 's')) % str(y))

					if rs[i][j] == re:
						self.assertRaises(TypeError, x.__pow__, y)
						if verbose > 1:
							columns.append('%s' % re)
					else:
						z = x**y
						self.assertEqual(z, rs[i][j])

						if verbose > 1:
							if isinstance(z, self.PQU):
								columns.append('%s' % ('PQU(' + str(z) + ')'))
							else:
								columns.append('%s' % str(z))

					if verbose > 1: 
						rows.append(''.join(columns))
						lenOfRows.append(len(''.join(columns)))
						del columns[1:]

			if verbose > 1:
				rowHead = eval("'%s'" % ('%-' + str(maxLenOfCol0) + 's%-' + str(maxLenOfCol1) + 's%-s')) % ('x', 'y', 'x**y')
				lenOfRows.append(len(rowHead))
				rowLine = '%-s' % (max(lenOfRows) * '-')
				rows.insert(0, rowLine)
				rows.insert(0, rowHead)
				rows.insert(0, rowLine)
				rows.append(rowLine)
				print '\n'.join(rows)


	if doTestAbsolute or doTestAll:
		def testAbsolute(self):
			xs = [self.PQU('1.5'), self.PQU('-1.5'), self.PQU('-1.5 eV'), self.PQU('1.5 eV')]
			rs = [1.5, self.PQU('1.5'), self.PQU('1.5 eV'), self.PQU('1.5 eV')]

			if verbose > 0:
				print '\ntesting __abs__'

			if verbose > 1:
				rows, lenOfRows, columns, lenOfCol0 = [], [], [], []
				
				for i, x in enumerate(xs):
					lenOfCol0.append(len(str(x)) + 6)

				maxLenOfCol0 = max(lenOfCol0)

			for i,x in enumerate(xs):
				result = abs(x)
				self.assertEqual(result, rs[i])

				if verbose > 1:
					columns.__init__()

					if isinstance(x, self.PQU):
						columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol0) + 's')) % ('PQU(' + str(x) + ')'))
					else:
						columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol0) + 's')) % str(x))

					if isinstance(result, self.PQU):
						columns.append('%s' % ('PQU(' + str(result) + ')'))
					else:
						columns.append('%s' % str(result))

					if verbose > 1: 
						rows.append(''.join(columns))
						lenOfRows.append(len(''.join(columns)))

			if verbose > 1:
				rowHead = eval("'%s'" % ('%-' + str(maxLenOfCol0) + 's%-s')) % ('x', '|x|')
				lenOfRows.append(len(rowHead))
				rowLine = '%-s' % (max(lenOfRows) * '-')
				rows.insert(0, rowLine)
				rows.insert(0, rowHead)
				rows.insert(0, rowLine)
				rows.append(rowLine)
				print '\n'.join(rows)

	if doTestComparisonOperations or doTestAll:
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

			if verbose > 0:
				print '\ntesting comparison operations'

			if verbose > 1:
				rows, lenOfRows, columns = [], [], []
				lenOfCol0, lenOfCol1, lenOfCol2, lenOfCol3 = [], [], [], []
				lenOfCol4, lenOfCol5, lenOfCol6, lenOfCol7 = [], [], [], []

				for i, x in enumerate(xs):
					lenOfCol0.append(len(str(x)) + 6)

					for j, y in enumerate(ys):
						lenOfCol1.append(len(str(y)) + 6)

						for k, z in enumerate(zs):
							if rs[k][i][j] == re:
								eval('lenOfCol%i.append(len(re) + 1)' % (k + 2))
							else:
								result = eval('x %s y' % z)
								eval('lenOfCol%i.append(len(str(result)) + 2)' % (k + 2))

				maxLenOfCol0 = max(lenOfCol0)
				maxLenOfCol1 = max(lenOfCol1)
				maxLenOfCol2 = max(lenOfCol2)
				maxLenOfCol3 = max(lenOfCol3)
				maxLenOfCol4 = max(lenOfCol4)
				maxLenOfCol5 = max(lenOfCol5)
				maxLenOfCol6 = max(lenOfCol6)
				maxLenOfCol7 = max(lenOfCol7) - 1

			for i, x in enumerate(xs):
				if verbose > 1:
					columns.__init__()

					if isinstance(x, self.PQU):
						columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol0) + 's')) % ('PQU(' + str(x) + ')'))
					else:
						columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol0) + 's')) % str(x))

				for j, y in enumerate(ys):
					if verbose > 1:
						if isinstance(y, self.PQU):
							columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol1) + 's')) % ('PQU(' + str(y) + ')'))
						else:
							columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol1) + 's')) % str(y))

					for k, z in enumerate(zs):
						if rs[k][i][j] == re:
							exec('self.assertRaises(TypeError, x.%s, y)' % ss[k])

							if verbose > 1:
								columns.append(eval("'%s'" % ('%-' + str(eval('maxLenOfCol%i' % (k + 2))) + 's')) % re)
						else:
							result = eval('x %s y' % z)
							self.assertEqual(result, rs[k][i][j])

							if verbose > 1:
								if isinstance(result, self.PQU):
									columns.append(eval("'%s'" % ('%-' + str(eval('maxLenOfCol%i' % (k + 2))) + 's')) % ('PQU(' + str(result) + ')'))
								else:
									columns.append(eval("'%s'" % ('%-' + str(eval('maxLenOfCol%i' % (k + 2))) + 's')) % str(result))

					if verbose > 1:
						rows.append(''.join(columns))
						lenOfRows.append(len(''.join(columns)))
						del columns[1:]

			if verbose > 1:
				rowHead = eval("'%s'" % ('%-' + str(maxLenOfCol0) + 's%-' + str(maxLenOfCol1) + 's%-' + str(maxLenOfCol2) + 's%-' + str(maxLenOfCol3) + 's%-' + str(maxLenOfCol4) + 's%-' + str(maxLenOfCol5) + 's%-' + str(maxLenOfCol6) + 's%-' + str(maxLenOfCol7) + 's')) % ('x', 'y', 'x == y', 'x != y', 'x <= y', 'x >= y', 'x < y', 'x > y')
				lenOfRows.append(len(rowHead))
				rowLine = '%-s' % (max(lenOfRows) * '-')
				rows.insert(0, rowLine)
				rows.insert(0, rowHead)
				rows.insert(0, rowLine)
				rows.append(rowLine)
				print '\n'.join(rows)

	if doTestConversions or doTestAll:
		def testConversions(self):
			re = 'raise exception'
			xs = [self.PQU('2.718'), self.PQU('-40 degC'), self.PQU('3.14 W')]
			ys = ['', 'K', 'degF', 'J/s', 'MeV/s', 'A*V']
			rGetValueAs = [[2.718, re, re, re, re, re], [re, 233.14999999999998, -40.000000000000057, re, re, re], [re, re, re, 3.14, 19598329980115.25, 3.14]]
			rInUnitsOf = [[self.PQU('2.718'), re, re, re, re, re], [re, self.PQU('233.14999999999998 K'), self.PQU('-40.000000000000057 degF'), re, re, re], [re, re, re, self.PQU('3.14 J/s'), self.PQU('19598329980115.25 MeV/s'), self.PQU('3.14 A*V')]]
			rIsCompatible = [[True, False, False, False, False, False], [False, True, True, False, False, False], [False, False, False, True, True, True]]
			rs = [rGetValueAs, rInUnitsOf, rIsCompatible]
			zs = ['getValueAs', 'inUnitsOf', 'isCompatible']

			if verbose > 0:
				print '\ntesting conversion'

			if verbose > 1:
				rows, lenOfRows, columns = [], [], []
				lenOfCol0, lenOfCol1, lenOfCol2 = [], [], []
				lenOfCol3, lenOfCol4, lenOfCol5 = [], [], []

				for i, x in enumerate(xs):
					lenOfCol0.append(len(str(x)) + 6)

					for j, y in enumerate(ys):
						lenOfCol1.append(len(str(y)) + 1)

						for k, z in enumerate(zs): 
							if rs[k][i][j] == re:
								eval('lenOfCol%i.append(len(re) + 1)' % (k + 2))
							else:
								result = eval('x.%s(y)' % z)
								if isinstance(result, self.PQU):
									eval('lenOfCol%i.append(len(str(result)) + 6)' % (k + 2))
								else:
									eval('lenOfCol%i.append(len(str(result)) + 1)' % (k + 2))

				maxLenOfCol0 = max(lenOfCol0)
				maxLenOfCol1 = max(lenOfCol1)
				maxLenOfCol2 = max(lenOfCol2)
				maxLenOfCol3 = max(lenOfCol3)
				maxLenOfCol4 = max(lenOfCol4) - 1

			for i, x in enumerate(xs):
				if verbose > 1:
					columns.__init__()
					if isinstance(x, self.PQU):
						columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol0) + 's')) % ('PQU(' + str(x) + ')'))
					else:
						columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol0) + 's')) % str(x))

				for j, y in enumerate(ys):
					if verbose > 1:
						if isinstance(y, self.PQU):
							columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol1) + 's')) % ('PQU(' + str(y) + ')'))
						else:
							columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol1) + 's')) % str(y))

					for k, z in enumerate(zs): 
						if rs[k][i][j] == re:
							exec('self.assertRaises(TypeError, x.%s, y)' % z)

							if verbose > 1:
								columns.append(eval("'%s'" % ('%-' + str(eval('maxLenOfCol%i' % (k + 2))) + 's')) % re)
						else:
							result = eval('x.%s(y)' % z)
							self.assertEqual(result, rs[k][i][j])

							if verbose > 1:
								if isinstance(result, self.PQU):
									columns.append(eval("'%s'" % ('%-' + str(eval('maxLenOfCol%i' % (k + 2))) + 's')) % ('PQU(' + str(result) + ')'))
								else:
									columns.append(eval("'%s'" % ('%-' + str(eval('maxLenOfCol%i' % (k + 2))) + 's')) % str(result))

					if verbose > 1:
						rows.append(''.join(columns))
						lenOfRows.append(len(''.join(columns)))
						del columns[1:]

			if verbose > 1:
				rowHead = eval("'%s'" % ('%-' + str(maxLenOfCol0) + 's%-' + str(maxLenOfCol1) + 's%-' + str(maxLenOfCol2) + 's%-' + str(maxLenOfCol3) + 's%-' + str(maxLenOfCol4) + 's')) % ('x', 'y', 'getValueAs', 'inUnitsOf', 'isCompatible')
				lenOfRows.append(len(rowHead))
				rowLine = '%-s' % (max(lenOfRows) * '-')
				rows.insert(0, rowLine)
				rows.insert(0, rowHead)
				rows.insert(0, rowLine)
				rows.append(rowLine)
				print '\n'.join(rows)

	if doTestConvertToUnit or doTestAll:
		def testConvertToUnit(self):
			re = 'raise exception'
			xs = [self.PQU('2.718'), self.PQU('-40 degC'), self.PQU('3.14 W')]
			ys = ['', 'K', 'degF', 'J/s', 'MeV/s', 'A*V']
			rs = [[self.PQU('2.718'), re, re, re, re, re], [re, self.PQU('233.14999999999998 K'), self.PQU('-40.000000000000057 degF'), re, re, re], [re, re, re, self.PQU('3.14 J/s'), self.PQU('19598329980115.25 MeV/s'), self.PQU('3.14 A*V')]]

			if verbose > 0:
				print '\ntesting convertToUnit'

			if verbose > 1:
				rows, lenOfRows, columns = [], [], []
				lenOfCol0, lenOfCol1, lenOfCol2 = [], [], []

				for i, x in enumerate(xs):
					lenOfCol0.append(len(str(x)) + 6)

					for j, y in enumerate(ys):
						lenOfCol1.append(len(str(y)) + 1)

						if rs[i][j] == re:
							lenOfCol2.append(len(re) + 1)
						else:
							lenOfCol2.append(len(str(x.convertToUnit(y))) + 6)

				maxLenOfCol0 = max(lenOfCol0)
				maxLenOfCol1 = max(lenOfCol1)
				maxLenOfCol2 = max(lenOfCol2) - 1

			xs = [self.PQU('2.718'), self.PQU('-40 degC'), self.PQU('3.14 W')]

			for i, x in enumerate(xs):
				if verbose > 1:
					columns.__init__()
					if isinstance(x, self.PQU):
						columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol0) + 's')) % ('PQU(' + str(x) + ')'))
					else:
						columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol0) + 's')) % str(x))

				for j, y in enumerate(ys):
					if verbose > 1:
						if isinstance(y, self.PQU):
							columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol1) + 's')) % ('PQU(' + str(y) + ')'))
						else:
							columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol1) + 's')) % str(y))

					if rs[i][j] == re:
						self.assertRaises(TypeError, x.convertToUnit, y)

						if verbose > 1:
							columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol2) + 's')) % re)
					else:
						x.convertToUnit(y)
						self.assertEqual(x, rs[i][j])

						if verbose > 1:
							if isinstance(x, self.PQU):
								columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol2) + 's')) % ('PQU(' + str(x) + ')'))
							else:
								columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol2) + 's')) % str(x))

					if verbose > 1:
						rows.append(''.join(columns))
						lenOfRows.append(len(''.join(columns)))
						del columns[1:]

			if verbose > 1:
				rowHead = eval("'%s'" % ('%-' + str(maxLenOfCol0) + 's%-' + str(maxLenOfCol1) + 's%-' + str(maxLenOfCol2) + 's')) % ('x', 'y', 'convertToUnit')
				lenOfRows.append(len(rowHead))
				rowLine = '%-s' % (max(lenOfRows) * '-')
				rows.insert(0, rowLine)
				rows.insert(0, rowHead)
				rows.insert(0, rowLine)
				rows.append(rowLine)
				print '\n'.join(rows)


	if doTestInBaseUnits or doTestAll:
		def testInBaseUnits(self):
			re = 'raise exception'
			xs = [self.PQU('2.718'), self.PQU('-40 degC'), self.PQU('3.14 W')]
			rs = [re, self.PQU('233.14999999999998 K'), self.PQU('3.14 kg*m**2/s**3')]

			if verbose > 0:
				print '\ntesting inBaseUnits'

			if verbose > 1:
				rows, lenOfRows, columns = [], [], []
				lenOfCol0, lenOfCol1 = [], []

				for i, x in enumerate(xs):
					lenOfCol0.append(len(str(x)) + 6)

					if rs[i] == re:
						lenOfCol1.append(len(re) + 1)
					else:
						result = x.inBaseUnits()
						lenOfCol1.append(len(str(result)) + 6)

				maxLenOfCol0 = max(lenOfCol0)
				maxLenOfCol1 = max(lenOfCol1) - 1

			for i, x in enumerate(xs):
				if verbose > 1:
					columns.__init__()

					if isinstance(x, self.PQU):
						columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol0) + 's')) % ('PQU(' + str(x) + ')'))
					else:
						columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol0) + 's')) % str(x))

				if rs[i] == re:
					self.assertRaises(TypeError, x.inBaseUnits)

					if verbose > 1:
						columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol1) + 's')) % re)
				else:
					result = x.inBaseUnits()
					self.assertEqual(result, rs[i])

					if verbose > 1:
						if isinstance(result, self.PQU):
							columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol1) + 's')) % ('PQU(' + str(result) + ')'))
						else:
							columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol1) + 's')) % str(result))

				if verbose > 1:
					rows.append(''.join(columns))
					lenOfRows.append(len(''.join(columns)))
					del columns[1:]

			if verbose > 1:
				rowHead = eval("'%s'" % ('%-' + str(maxLenOfCol0) + 's%-' + str(maxLenOfCol1) + 's')) % ('x', 'inBaseUnits')
				lenOfRows.append(len(rowHead))
				rowLine = '%-s' % (max(lenOfRows) * '-')
				rows.insert(0, rowLine)
				rows.insert(0, rowHead)
				rows.insert(0, rowLine)
				rows.append(rowLine)
				print '\n'.join(rows)

	if doTestTrigFunctions or doTestAll:
		def testTrigFunctions(self):
			re = 'raise exception'
			xs = [self.PQU('30 deg'), self.PQU(math.pi/4, 'rad')]
			rs = [[0.5, math.sqrt(3 / 4.), 1 / math.sqrt(3.)], [1 / math.sqrt(2.), 1 / math.sqrt(2.), 1]]
			zs = ['sin', 'cos', 'tan']

			if verbose > 0:
				print '\ntesting Trig functions'

			if verbose > 1:
				rows, lenOfRows, columns = [], [], []
				lenOfCol0, lenOfCol1, lenOfCol2, lenOfCol3 = [], [], [], []

				for i, x in enumerate(xs):
					lenOfCol0.append(len(str(x)) + 6)

					for k, z in enumerate(zs): 
						if rs[i][k] == re:
							eval('lenOfCol%i.append(len(re) + 1)' % (k + 1))
						else:
							result = eval('x.%s()' % z)
							eval('lenOfCol%i.append(len(str(result)) + 1)' % (k + 1))

				maxLenOfCol0 = max(lenOfCol0)
				maxLenOfCol1 = max(lenOfCol1)
				maxLenOfCol2 = max(lenOfCol2)
				maxLenOfCol3 = max(lenOfCol3) - 1

			for i, x in enumerate(xs):
				if verbose > 1:
					columns.__init__()
					if isinstance(x, self.PQU):
						columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol0) + 's')) % ('PQU(' + str(x) + ')'))
					else:
						columns.append(eval("'%s'" % ('%-' + str(maxLenOfCol0) + 's')) % str(x))

					for k, z in enumerate(zs): 
						if rs[i][k] == re:
							exec('self.assertRaises(TypeError, x.%s)' % z)

							if verbose > 1:
								columns.append(eval("'%s'" % ('%-' + str(eval('maxLenOfCol%i' % (k + 1))) + 's')) % re)
						else:
							result = eval('x.%s()' % z)
							self.assertAlmostEqual(result, rs[i][k])

							if verbose > 1:
								if isinstance(result, self.PQU):
									columns.append(eval("'%s'" % ('%-' + str(eval('maxLenOfCol%i' % (k + 1))) + 's')) % ('PQU(' + str(result) + ')'))
								else:
									columns.append(eval("'%s'" % ('%-' + str(eval('maxLenOfCol%i' % (k + 1))) + 's')) % str(result))

					if verbose > 1:
						rows.append(''.join(columns))
						lenOfRows.append(len(''.join(columns)))
						del columns[1:]

			if verbose > 1:
				rowHead = eval("'%s'" % ('%-' + str(maxLenOfCol0) + 's%-' + str(maxLenOfCol1) + 's%-' + str(maxLenOfCol2) + 's%-' + str(maxLenOfCol3) + 's')) % ('x', 'sin', 'cos', 'tan')
				lenOfRows.append(len(rowHead))
				rowLine = '%-s' % (max(lenOfRows) * '-')
				rows.insert(0, rowLine)
				rows.insert(0, rowHead)
				rows.insert(0, rowLine)
				rows.append(rowLine)
				print '\n'.join(rows)

if __name__ == '__main__':
	unittest.main()

