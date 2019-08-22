# This is *now more than* a slight mod of the
# Scientific/Physics/PhysicalQuantities.py package, Physical quantities with 
# units. Written by Konrad Hinsen <hinsen@cnrs-orleans.fr> with contributions
# from Greg Ward. Last revision: 2007-5-25

# What is lacking?
#   -) Cannot handle unitless quantities 
#      (e.g., "a = PhysicalQuantityWithUncertainty('1')") and in
#      "a = PhysicalQuantityWithUncertainty('1 m'); b = a / a", b is returned as
#      a float not a PhysicalQuantityWithUncertainty.
#   -) Cannot handle numbers starting with ".", (e.g., "0.001 eV" works but
#      ".001 eV" does not.

# Modified by Nidhi R. Patel <infinidhi@llnl.gov> and Bret Beck <beck6@llnl.gov>
# Last revision: 2012-7-27

# What is new?
# 1. Added support for dimensionless quantity

# What is in the works?
# 1. Various ways to define physical quantity with uncertainty
# 2. Significant figures
# 3. '%' as unit
# 4. Linear propagation of errors
# 5. Fundamental constants read from NIST codata text file


"""
.. rubric:: PHYSICAL QUANTITY WITH UNCERTAINTY

This module provides a data type that represents a physical quantity together 
with its unit. It is possible to add and subtract these quantities if the units 
are compatible, and a quantity can be converted to another compatible unit. 
Multiplication, subtraction, and raising to integer powers is allowed without 
restriction, and the result will have the correct unit. A quantity can be raised
to a non-integer power only if the result can be represented by integer powers 
of the base units.

The values of physical constants are taken from the 1986 recommended values from
CODATA. Other conversion factors (e.g. for British units) come from various 
sources. I can't guarantee for the correctness of all entries in the unit table,
so use this at your own risk.
"""

import sys
import math
import re, string
from functools import reduce

__metaclass__ = type

def _getUnit(unit, default = None):
    """
    "unit" can only be a string (e.g., '', 'MeV', 'b * eV') or None.
    """

    if isinstance(default, PhysicalQuantityWithUncertainty):
        default = default.unit
    else:
        default = None


    if unit is None:
        if default is None:
            return None

        return _findUnit(default)

    if unit == '':
        return ''

    return _findUnit(unit)


def valueOrPQ(value, unitFrom = None, unitTo = None, asPQU = False):
    """
    This functions takes a float or a PhysicalQuantityWithUncertainty as input and 
    returns, base on asPQU, a float or PhysicalQuantityWithUncertainty. If the 
    unitTo is not None, the returned value is in units of unitTo. As 
    PhysicalQuantityWithUncertainty does not currently support a dimensionless 
    quantity, if the input is a float and unitFrom is not specified then a float is 
    returned, independent of asPQU.f value is a PhysicalQuantityWithUncertainty and 
    unitFrom is not None, they must be the same unit.

    Options for input are:
        +-----+--------+------+-----------------+-----------------------------------+
        |     |        |      |  from _getUnit  |       returned value asPQU        |
        +-----+--------+------+--------+--------+--------+--------------------------+
        |value|unitFrom|unitTo|unitFrom|unitTo  |False   |True                      |
        +-----+--------+------+--------+--------+--------+--------------------------+
        |float|None    |None  |None    |None    |float   |float                     |
        |     |        |      |        |        |        |(until PQU supports       |
        |     |        |      |        |        |        |dimensionless unit)       |
        +-----+--------+------+--------+--------+--------+--------------------------+
        |float|string  |None  |'' or PU|None    |float in|PQ in unitFrom            |
        |     |        |      |        |        |unitFrom|                          |
        +-----+--------+------+--------+--------+--------+--------------------------+
        |float|string  |string|'' or PU|'' or PU|float in|PQ in unitTo              |
        |     |        |      |        |        |unitTo  |                          |
        +-----+--------+------+--------+--------+--------+--------------------------+
        |PQU  |None    |None  |PU      |None    |float   |PQ                        |
        +-----+--------+------+--------+--------+--------+--------------------------+
        |PQU  |string  |None  |'' or PU|None    |float in|PQU (need to check that PQ|
        |     |        |      |        |        |unitFrom|and unitFrom are the same)|
        +-----+--------+------+--------+--------+--------+--------------------------+
        |PQU  |string  |string|'' or PU|'' or PU|float in|PQU in unitTo (ditto)     | 
        |     |        |      |        |        |unitTo  |                          |
        +-----+--------+------+--------+--------+--------+--------------------------+
        |PQU  |None    |string|PU      |'' or PU|float in|PQU in unitTo             |
        |     |        |      |        |        |unitTo  |                          |
        +-----+--------+------+--------+--------+--------+--------------------------+
    """
    unitFrom = _getUnit(unitFrom, value)

    if unitFrom is None:
        if unitTo is None:
            return value

        raise TypeError('Incompatible units: unitFrom is None while unitTo is "%s"' % unitTo)

    if isinstance(value, PhysicalQuantityWithUncertainty):
        if isinstance(unitFrom, PhysicalUnit):
            if unitFrom != value.unit:
                raise TypeError('Incompatible units: unitFrom is "%s" and value has "%s"' % (unitFrom, value.unit))
            else:
                raise TypeError('Incompatible units: value has = "%s" and unitFrom is ""' % value.unit)

        value = value.value

    # At this point value is a float and unitFrom is either '' or PU.
    if not (isinstance(unitFrom, PhysicalUnit)): # Is '' (e.g, dimensionless)
        if (unitTo is None) or (unitTo == ''):
            return(value)

        raise TypeError('Incompatible units: unitFrom is "%s" and unitTo is "%s"' % (unitFrom, unitTo))

    # At this point unitFrom has dimensions (i.e., is not '').
    if unitTo is None:
        unitTo = unitFrom
    else:
        if unitTo == '':
            raise TypeError('Incompatible units: unitFrom is "%s" and unitTo is ""' % unitFrom)

        unitTo = _findUnit(unitTo)

        if unitFrom != unitTo:
            value = PhysicalQuantityWithUncertainty(value, unitFrom).getValueAs(unitTo)

    if asPQU:
        return PhysicalQuantityWithUncertainty(value, unitTo)

    return value

def toShortestString(value, precision = 15):

    e = ('%%.%de' % (precision-1)) % value
    e1 = '%e' % value

    if float(e) == float(e1):
        e = e1

    s = e.split('e')
    m = s[0].rstrip('0')

    if m[-1] == '.':
        m = m[:-1]

    p = int( s[1] )
    e = m + 'e' + "%d" % p
    p = precision - p - 1

    if p > 25:
        return e

    p = max(p, 0)
    f = ('%%.%df' % p) % value

    if 'e' in f:
    # For very large numbers, the '%.?f' format returns 'e' format (e.g., .1234e300 --> "1e+299")
        return e

    if '.' in f:
        f = f.rstrip('0')

    if f[-1] == '.':
        f = f[:-1]

    if len(f) < len(e):
        return f

    return e


# Class definitions
class PatternRecognition:
    """
    .. rubric:: Regular Expression Matching
    """

    sign = '[+-]?' # optional sign
    num = '(%s\d*)(\.(\d*))?' % sign # value, can be decimal or integer
    exp = '(([eE])(%s\d+))?' % sign # optional exponent
    value = '^(%s%s)$' % (num, exp)
    valueWithPercent = '^((%s%s)(%s))$' % (num, exp, '\%')
    valueWithUnit = '^(%s%s)' % (num, exp)

    @staticmethod
    def number(self):
        """
        :param self: a number
        :type self: `str`
        :returns: components of the number
        :rtype: `tuple` of `str`
        :raises TypeError: if the number string is incompatible 
        """

        try:
            s0 = self.strip(' ')
            m = re.compile(PatternRecognition.value).match(s0)

            if isinstance(m, type(None)):
                raise TypeError('\n\n\tincompatible string\n\t(unable to parse value in %s)\n\n' % self)
            else:
                return m.groups()
        except:
            raise TypeError('\n\n\tincompatible string\n\t(unable to parse number in %s)\n\n' % self)

    @staticmethod
    def isNumber(self):
        """
        :param self: a number
        :type self: `str`
        :returns: `True` if self is a number
        :rtype: `bool`
        """

        try:
            m = PatternRecognition.number(self)
            return True
        except:
            return False

    @staticmethod
    def numberWithPercent(self):
        """
        :param self: a number with percent
        :type self: `str`
        :returns: components of the number with percent
        :rtype: `tuple` of `str`
        :raise TypeError: if the number with percent string is incompatible 
        """

        try:
            s0 = self.strip(' ')
            m = re.compile(PatternRecognition.valueWithPercent).match(s0)

            if isinstance(m, type(None)):
                raise TypeError("\n\n\tincompatible string\n\t(unable to parse valueWithPercent in '%s')\n\n" % self)
            else:
                return m.groups()
        except:
            raise TypeError("\n\n\tincompatible string\n\t(unable to parse numberWithPercent in '%s')\n\n" % self)

    @staticmethod
    def isNumberWithPercent(self):
        """
        :param self: a number with percent
        :type self: `str`
        :returns: `True` if self is a number with percent
        :rtype: `bool`
        """

        try:
            m = PatternRecognition.numberWithPercent(self)
            return True
        except:
            return False

    @staticmethod
    def numberWithUnit(self):
        """
        :param self: a number with units
        :type self: `str`
        :returns: components of the number with units
        :rtype: `tuple` of `str` and `PhysicalUnit`
        :raises TypeError: if the self or number is incompatible 
        """

        s0 = self.strip(' ').partition(' ')
        s00 = s0[0].strip(' ')
        s01 = s0[2].strip(' ')
        m = re.compile(PatternRecognition.valueWithUnit).match(self.strip(' '))
        v = re.compile(PatternRecognition.value).match(s00)
        u = _findUnit(s01)

        if isinstance(m, type(None)):
            raise TypeError('\n\n\tincompatible string\n\t(unable to parse numberWithUnit in %s)\n\n' % self)
        else:
            if isinstance(v, type(None)):
                raise TypeError('\n\n\tincompatible string\n\t(unable to parse number in %s)\n\n' % self)

            if isinstance(u, type(None)):
                raise TypeError('\n\n\tincompatible string\n\t(unable to parse number or unit in %s)\n\n' % self)
            return m.groups(), u

    @staticmethod
    def isNumberWithUnit(self):
        """
        :param self: a number with units
        :type self: `str`
        :returns: `True` if self is a number with units
        :rtype: `bool`
        """

        try:
            m = PatternRecognition.numberWithUnit(self)
            return True
        except:
            return False

    @staticmethod
    def numberWithUncertainty(self):
        """
        :param self: a number with uncertainty and/or units
        :type self: `str`
        :rtype: `dict`
        :raises TypeError: if the number with uncertainty and/or units string is incompatible 
        """

        if isinstance(self, str):
            s = self.strip(' ')
            x = {}

            if s.__contains__('+/-'):
                if s.count('+/-') == 1:
                    s0 = s.split('+/-')
                    m0 = PatternRecognition.number(s0[0])
                    x['value'] = m0[0]
                    m01 = PatternRecognition.numberWithUnit(s0[1].replace('%', ''))
                    x['uncertainty'] = m01[0][0]
                    x['unit'] = m01[1]

                    if s0[1].__contains__('%'):
                        x['uncertaintyType'] = '%'
                        x['absoluteUncertainty'] = toShortestString(float(x['uncertainty']) / 100 * float(x['value']))
                    else:
                        x['uncertaintyType'] = 'absolute'
                        x['absoluteUncertainty'] = x['uncertainty']

                    if float(x['absoluteUncertainty']) != 0.0:
                        if float(x['value']) != 0.0:
                            y = PatternRecognition.toLeastPreciseNumberWithUncertainty(x['value'], x['absoluteUncertainty'])
                            x['value'] = y['value']

                            if float(y['uncertainty']) != 0.0:
                                x['absoluteUncertainty'] = y['uncertainty']
                                x['significantDigits'] = PatternRecognition.findSignificantDigits(x['value'])

                                if x['uncertaintyType'] == 'absolute':
                                    x['uncertainty'] = x['absoluteUncertainty']
                            else:
                                pass
                                #raise TypeError('\n\n\tincompatible uncertainty value\n\t(absolute uncertainty cannot be equal to zero\n\twhen truncated to the same precision as the value)\n\n')
                        else:
                            sD = PatternRecognition.findSignificantDigits(x['absoluteUncertainty'])

                            if sD > 2:
                                x['significantDigits'] = 2
                            else:
                                x['significantDigits'] = sD

                            x['absoluteUncertainty'] = PatternRecognition.toSignificantDigits(x['absoluteUncertainty'], x['significantDigits'])
                            x['uncertainty'] = x['absoluteUncertainty']
                    else:
                        raise TypeError('\n\n\tincompatible uncertainty value\n\t(absolute uncertainty cannot be equal to zero)')
                else:
                    raise TypeError("\n\n\tincompatible string\n\t(too many '+/-' in %s\n\n)" %s)

            elif s.__contains__('(') and s.__contains__(')'):
                v = s.split('(')[0].strip(' ')
                u = s.split('(')[1].split(')')[0].strip(' ')
                e = s.split('(')[1].split(')')[1]
                ve = v + e
                m0 = PatternRecognition.numberWithUnit(ve)
                x['value'] = m0[0][0]
                x['unit'] = m0[1]
                x['uncertaintyType'] = 'absolute'
                ue = u + e
                pv = PatternRecognition.findPrecision(v)
                m1 = PatternRecognition.numberWithUnit(ue)

                if isinstance(m1[0][5], type(None)):
                    x['uncertainty'] = m1[0][0] + 'e' + str(pv)
                else:
                    pe = m1[0][0].find(m1[0][5])
                    ps = m1[0][0].find(' ')
                    x['uncertainty'] = m1[0][0][0:pe] + 'e' + str(pv + int(m1[0][6]))

                x['absoluteUncertainty'] = x['uncertainty']

                if float(x['absoluteUncertainty']) != 0.0:
                    if float(x['value']) != 0.0:
                        y = PatternRecognition.toLeastPreciseNumberWithUncertainty(x['value'], x['absoluteUncertainty'])
                        x['value'] = y['value']

                        if float(y['uncertainty']) != 0.0:
                            x['uncertainty'] = y['uncertainty']
                            x['absoluteUncertainty'] = y['uncertainty']
                            x['significantDigits'] = PatternRecognition.findSignificantDigits(x['value'])
                        else:
                            pass
                            #raise TypeError('\n\n\tincompatible uncertainty value\n\t(absolute uncertainty cannot be equal to zero\n\twhen truncated to the same precision as the value)\n\n')
                    else:
                        sD = PatternRecognition.findSignificantDigits(x['absoluteUncertainty'])

                        if sD > 2:
                            x['significantDigits'] = 2
                        else:
                            x['significantDigits'] = sD

                        x['absoluteUncertainty'] = PatternRecognition.toSignificantDigits(x['absoluteUncertainty'], x['significantDigits'])
                        x['uncertainty'] = x['absoluteUncertainty']

                else:
                    raise TypeError('\n\n\tincompatible uncertainty value\n\t(absolute uncertainty cannot be equal to zero)')

            else:
                m = PatternRecognition.numberWithUnit(s)
                x['value'] = m[0][0]
                x['unit'] = m[1]

        else:
            raise TypeError('\n\n\tincompatible string\n\t(argument is not string type)\n\n' %s)

        return x

    @staticmethod
    def findLeastPrecise(number1, number2):
        """
        The precision of the number1 and number2 are compared 

        :param number1: a number
        :type number1: `str`
        :param number2: a number
        :type number2: `str`
        :returns: 2 or 1 as key and precision as the value
        :rtype: `dict`
        """

        number1P = PatternRecognition.findPrecision(number1)
        number2P = PatternRecognition.findPrecision(number2)

        if number1P <= number2P:
            return {2: number2P}
        else:
            return {1: number1P}

    @staticmethod
    def toSamePrecision(number1, number2):
        """
        Number1 and number2 are adjusted so that they both have the same precision.

        :param number1: a number
        :type number1: `str`
        :param number2: a number
        :type number2: `str`
        :returns: adjusted number1 and number2
        :rtype: `list` of `str`
        """

        leastPrecise = PatternRecognition.findLeastPrecise(number1, number2)

        if leastPrecise.keys()[0] == 1:
            x = [number1, PatternRecognition.toPrecision(number2, leastPrecise.values()[0])]
        else:
            x = [PatternRecognition.toPrecision(number1, leastPrecise.values()[0]), number2]

        return x

    @staticmethod
    def toLeastPreciseNumberWithUncertainty(value, uncertainty, significantDigitOfUncertainty = 2):
        """
        The uncertainty is reduced to the specified significant digits (default is 2).
        Then the precision of the value and its uncertainty are compared. Both 
        quantities are adjusted to have the least precision.

        :param value: a number
        :type value: `str`
        :param uncertainty: a number
        :type uncertainty: `str`
        :param significantDigitOfUncertainty: a number
        :type significantDigitOfUncertainty: `int`
        :returns: new value and its uncertainty
        :rtype: `dict`
        """

        uncertaintySD = PatternRecognition.findSignificantDigits(uncertainty)
        leastSignificantDigit = min(uncertaintySD, significantDigitOfUncertainty)
        uncertainty0 = PatternRecognition.toSignificantDigits(uncertainty, leastSignificantDigit)
        leastPrecise = PatternRecognition.toSamePrecision(value, uncertainty0)

        return {'value': leastPrecise[0], 'uncertainty': leastPrecise[1]}

    @staticmethod
    def findSignificantDigitsWithUncertainty(value, uncertainty):
        """
        Using findLeastPrecise method with 2 significant digits set to the uncertainty,
        it finds significant digits of the value.

        :param value: a number
        :type value: `str`
        :param uncertainty: a number
        :type uncertainty: `str`
        :returns: significant digit of a number
        :rtype: `int`
        """

        return PatternRecognition.findSignificantDigits(PatternRecognition.findLeastPreciseNumberWithUncertainty(value, uncertainty)['value'])

    @staticmethod
    def findPrecision(self):
        """
        Finds precision of a number, the digit's place following a decimal point.

        :param self: a number
        :type self: `str`
        :raises TypeError: if the number string is incompatible
        :returns: a negative integer
        :rtype: `int`
        """

        try:
            m = PatternRecognition.number(self)

            if m[3] == None:
                rightOfDecimal = 0
            else:
                rightOfDecimal = -len(m[3])

            if m[6] == None:
                exponent = 0
            else:
                exponent = int(m[6])

            precision = rightOfDecimal + exponent

            if precision > 0:
                precision = 0

            return precision
        except:
            raise TypeError('\n\n\tincompatible string\n\t(unable to find precision in %s)\n\n' % self)

    @staticmethod
    def toPrecision(self, desiredPrecision = 0, forcedPrecision = False):
        """
        Transform a number to desired precision. The desired precision must be less 
        than the precision of the number itself.

        :param self: a number
        :type self: `str`
        :param desiredPrecision: a number
        :type desiredPrecision: `int`
        :returns: a number with desired precision
        :rtype: `str`
        :raises TypeError: if the number string is incompatible and/or if the precision of the number > desired precision > 0
        """

        m = PatternRecognition.number(self)
        precision = PatternRecognition.findPrecision(self)

        #if desiredPrecision >= precision and desiredPrecision <= 0:
        if desiredPrecision >= precision - 1 and desiredPrecision <= 0:
            newSelf = eval("'%s'" % ('%.' + str(abs(desiredPrecision)) + 'f')) % float(m[0])
            if desiredPrecision == 0:
                newSelf += '.'
        else:
            raise TypeError('\n\n\tincompatible desiredPrecision\n\t(precision of the value, %i, must be <= desiredPrecision, %i, which must be <= 0)\n\n' % (precision, desiredPrecision))

        return newSelf

    @staticmethod
    def findSignificantDigits(self):
        """
        :param self: a number
        :type self: `str`
        :returns: significant digits in a number
        :rtype: `int`
        :raises TypeError: if the number string is incompatible
        """

        try:
            m = PatternRecognition.number(self)
            l = 0

            if m[1] == '' or m[1] == '0' or m[1] == '+' or m[1] == '-':
                l = len(m[3].lstrip('0'))
            else:
                m1 = m[1].lstrip('+').lstrip('-')

                if m[2] == None:
                    l = len(m1.strip('0'))
                elif m[2] == '.':
                    l = len(m1.lstrip('0'))
                else:
                    l = len(m1.lstrip('0')) + len(m[3])

            return l
        except:
            raise TypeError('\n\n\tincompatible string\n\n')

    @staticmethod
    def toSignificantDigits(self, toSignificantDigits = 1):
        """
        :param self: a number
        :type self: `str`
        :param toSignificantDigits: a number
        :type toSignificantDigits: `int`
        :returns: a number with desired significant digits
        :rtype: `str`
        :raises TypeError: if the number string is incompatible
        """

        try:
            if float(self) != 0.0:
                m = PatternRecognition.number(self)
                es = eval("'%s'" % ('%+.' + str(toSignificantDigits - 1) + 'e')) % float(self)
                m = PatternRecognition.number(es)

                if m[2] == None:
                    if self.__contains__('.'):
                        rs = m[1] + '.'
                    else:
                        rs = m[1]
                else:
                    rs = m[1] + m[2]

                if int(m[6]) == 0:
                    pass
                else:
                    rs += m[5] + str(int(m[6]))

                return rs.replace('+', '')
            else:
                es = eval("'%s'" % ('%.' + str(toSignificantDigits) + 'f')) % float(self)
                return es
        except:
            raise TypeError('\n\n\tincompatible string\n\n')


class PhysicalQuantityWithUncertainty:
    """
    .. rubric:: PHYSICAL QUANTITY WITH UNCERTAINTY

    PhysicalQuantityWithUncertainty instances allow addition, subtraction, 
    multiplication, and division with each other as well as multiplication, 
    division, and exponentiation with numbers. Addition and subtraction check that 
    the units of the two operands are compatible and return the result in the units 
    of the first operand. A limited set of mathematical functions (from module 
    Numeric) is applicable as well:

    - sqrt: equivalent to exponentiation with 0.5.
    - sin, cos, tan: applicable only to objects whose unit is compatible with 'rad'.

    See the documentation of the PhysicalQuantityWithUncertainty module for a 
    list of the available units.
    
    Here is an example on usage:

    >>> import pqu
    >>> from pqu.physicalQuantityWithUncertainty \\
    ... import PhysicalQuantityWithUncertainty as PQU
    >>> distance1 = PQU('10 m')
    >>> distance2 = PQU('10 km')
    >>> total = distance1 + distance2
    >>> total
    PhysicalQuantityWithUncertainty(10010.0, 'm')
    >>> total.convertToUnit('km')
    >>> total.getValue()
    10.01
    >>> total.getUnitSymbol()
    'km'
    >>> total = total.inBaseUnits()
    >>> total
    PhysicalQuantityWithUncertainty(10010.0, 'm')
    >>> t = PQU(314159., 's')
    >>> # convert to days, hours, minutes, and second:
    >>> t2 = t.inUnitsOf('d', 'h', 'min', 's')
    >>> t2_print = ' '.join([str(i) for i in t2])
    >>> t2_print
    '3.0 d 15.0 h 15.0 min 59.0 s'
    >>> e = PQU('2.7 Hartree * Nav')
    >>> e.convertToUnit('kcal/mol')
    >>> e
    PhysicalQuantityWithUncertainty(1694.2757596034764,'kcal/mol')
    >>> e = e.inBaseUnits()
    >>> str(e)
    '7088849.77818 kg * m**2/s**2/mol'
    >>> freeze = PQU('0 degC')
    >>> freeze = freeze.inUnitsOf('degF')
    >>> str(freeze)
    '32.0 degF'


    .. rubric:: There are three constructor calling patterns:

    1.  PhysicalQuantityWithUncertainty(value, unit, uncertainty, uncertaintyType, significantDigits)

        Usage:
    
            A. value can be a float, int or string if by itself or a string if with unit, uncertainty, & uncertaintyType
        
            B. unit must be a string
        
        Optional:
    
            a. uncertainty can be a float, int or str if by itself
        
            b. uncertaintyType must be a string. Choices are 'absolute', or '%'
        
            c. significantDigits must be an int.
    
        Examples:
    
            >>> from pqu.physicalQuantityWithUncertainty import 
            ...    PhysicalQuantityWithUncertainty as PQU
            >>> x = PQU(3.0, 'm')
            >>> print x.value, x.unit
            3.0 m
            >>> PQU(3., 'm', uncertainty = 0.1)
            >>> print x.value, x.unit, x.uncertainty, x.uncertaintyType, 
            ...    x.significantDigits
            3.0 m 0.1 absolute 2
            >>> PQU(3, 'm', uncertainty = 1, uncertaintyType = '%')
            >>> print x.value, x.unit, x.uncertainty, x.uncertaintyType, 
            ...    x.significantDigits
            3.0 m 1.0 % 1
            >>> PQU(3, 'm', uncertainty = '1%')
            >>> print x.value, x.unit, x.uncertainty, x.uncertaintyType, 
            ...    x.significantDigits
            3.0 m 1.0 % 2

    2.  PhysicalQuantityWithUncertainty(value_with_unit, uncertainty, 
        uncertaintyType, significantDigits), where value_with_unit is a string that 
        contains both the value and the unit, i.e. '1.5 m/s'.

    3.  PhysicalQuantityWithUncertainty(value) where value is a string that 
        containts the value, the uncertainty, the uncertaintyType, and the unit, i.e 
        PQU('1. +/- 0.1 m/s') for absolute or PQU('1.0 +/- 10% m/s)' for percent 
        uncertaintyType. Other formats are PQU('1.0(1) m/s'), and PQU('1.0(1)e-2 m/s').
        All of these PQU have value 1.0 and 0.1 as uncertainty. Value can also be 
        another PhysicalQuantityWithUncertainty instance and unit, uncertainty, and
        uncertaintyType can be specified to overwrite the old values.

    See test/testPhysicalQuantityWithUncertainty.py for more examples.
    """

    def __init__(self, value, unit = None, uncertainty = None, uncertaintyType = None, significantDigits = None):

        if isinstance(value, (float, int)):
            self.value = value

            if isinstance(unit, type(None)):
                self.unit = _findUnit('')
            else:
                self.unit = _findUnit(unit)

        elif isinstance(value, PhysicalQuantityWithUncertainty):
            self.value = value.value

            if isinstance(unit, type(None)):
                self.unit = value.unit
            else:
                self.unit = _findUnit(unit)

            if hasattr(value, 'uncertainty'):
                self.uncertainty = value.uncertainty
                self.uncertaintyType = value.uncertaintyType
                self.significantDigits = value.significantDigits
                self.absoluteUncertainty = value.absoluteUncertainty

        elif isinstance(value, str):
            s = value.strip(' ')
            match = PatternRecognition.numberWithUncertainty(s)
            self.value, self.unit = match['value'], match['unit']

            if match.has_key('uncertainty'):
                self.uncertainty = match['uncertainty']
                self.uncertaintyType = match['uncertaintyType']
                self.absoluteUncertainty = match['absoluteUncertainty']
                self.significantDigits = match['significantDigits']

            if not isinstance(unit, type(None)):
                if not self.unit.isDimensionless():
                    if self.unit != _findUnit(unit):
                        raise TypeError('\n\n\tinconsistent units\n\t(2 different units specified)\n\n')
        else:
            raise TypeError('\n\n\tincompatible value\n\t(value must a str or a float or PhysicalQuantityWithUncertainty)\n\n')

        if not isinstance(uncertainty, type(None)):
            unc, uncType, absUnc = 0, 0, 0

            if isinstance(uncertainty, (float, int)):
                unc = uncertainty

                if not isinstance(uncertaintyType, type(None)):
                    uncType = uncertaintyType

                    if uncType == '%':
                        absUnc = toShortestString(unc / 100 * float(self.value))
                    elif uncType == 'absolute':
                        absUnc = toShortestString(unc)
                    else:
                        raise TypeError('\n\n\tincompatible uncertaintyType\n\t(unable to parse uncertaintyType in %s\n\n' % str(uncertaintyType))
                else:
                    uncType = 'absolute'
                    absUnc = toShortestString(unc)

            elif isinstance(uncertainty, str):
                try:
                    m = PatternRecognition.numberWithPercent(uncertainty)
                    unc = m[1]
                    uncType = '%'
                    absUnc = toShortestString(float(unc) / 100 * float(self.value))
                except:
                    m = PatternRecognition.number(uncertainty)
                    unc = m[0]
                    absUnc = unc

                    if isinstance(uncertaintyType, type(None)):
                        uncType = 'absolute'
                    else:
                        uncType = uncertaintyType

                if not isinstance(uncertaintyType, type(None)):
                    if uncType != uncertaintyType:
                        raise TypeError('\n\n\tinconsistent uncertaintyType\n\t(2 different uncertaintyTypedefined\n\n')
            else:
                raise TypeError('\n\n\tincompatible uncertainty\n\t(uncertainty must be a number or a number with percent\n\tin float, int, or str type)\n\n')

            if float(absUnc) != 0.0:
                if float(self.value) != 0.0:
                    y = PatternRecognition.toLeastPreciseNumberWithUncertainty(str(self.value), str(absUnc))
                    self.value = y['value']

                    if float(y['uncertainty']) != 0.0:
                        if uncType == '%':
                            self.uncertainty = unc
                            self.uncertaintyType = '%'
                            self.absoluteUncertainty = y['uncertainty']
                        else:
                            self.uncertainty = y['uncertainty']
                            self.uncertaintyType = 'absolute'
                            self.absoluteUncertainty = y['uncertainty']

                        self.significantDigits = PatternRecognition.findSignificantDigits(self.value)
                    else:
                        pass
                        #raise TypeError('\n\n\tincompatible uncertainty value\n\t(absolute uncertainty cannot be equal to zero\n\twhen truncated to the same precision as the value)\n\n')
                else:
                    self.uncertaintyType = 'absolute'
                    sD = PatternRecognition.findSignificantDigits(absUnc)

                    if sD > 2:
                        self.significantDigits = 2
                    else:
                        self.significantDigits = sD

                    self.absoluteUncertainty = PatternRecognition.toSignificantDigits(absUnc, self.significantDigits)
                    self.uncertainty = self.absoluteUncertainty

                if not isinstance(significantDigits, type(None)):
                    if significantDigits != self.significantDigits:
                        raise TypeError('\n\n\tinconsistent significantDigits\n\t(calculated significantDigits != significantDigits)\n\n')
            else:
                raise TypeError('\n\n\tincompatible uncertainty value\n\t(absolute uncertainty cannot be equal to zero)\n\n')

        if not isinstance(unit, type(None)):
            if self.unit.isDimensionless():
                self.unit = _findUnit(unit)
            else:
                if self.unit != _findUnit(unit):
                    raise TypeError('\n\n\tincompatible units\n\t(2 different units specified)')

        self.value = float(self.value)

        try:
            self.uncertainty = float(self.uncertainty)
            self.uncertaintyType = str(self.uncertaintyType)
            self.absoluteUncertainty = float(self.absoluteUncertainty)
            self.significantDigits = int(self.significantDigits)
        except:
            pass

    def __str__(self):
        s = toShortestString(self.value)

        if hasattr(self, 'uncertainty'):
            if self.uncertaintyType == '%':
                s += ' +/- ' + toShortestString(self.uncertainty) + '%'
            else:
                if self.value != 0.0:
                    s0 = PatternRecognition.toSignificantDigits(s, self.significantDigits)
                    s += ' +/- ' + PatternRecognition.toPrecision(toShortestString(self.absoluteUncertainty), PatternRecognition.findPrecision(s0))
                else:
                    s += ' +/- ' + PatternRecognition.toSignificantDigits(toShortestString(self.absoluteUncertainty), self.significantDigits)

        if not self.unit.isDimensionless():
            s += ' ' + self.unit.symbol()

        return s

    def __repr__(self):
        s = self.__class__.__name__ + '(' + repr(self.value)
        s += ', ' + repr(self.unit.symbol())

        if hasattr(self, 'uncertainty'):
            if self.uncertaintyType == '%':
                s += ', ' + repr(self.uncertainty)
            else:
                if self.value != 0.0:
                    s0 = PatternRecognition.toSignificantDigits(toShortestString(self.value), self.significantDigits)
                    s += ', ' + repr(PatternRecognition.toPrecision(toShortestString(self.absoluteUncertainty), PatternRecognition.findPrecision(s0)))
                else:
                    s += ', ' + repr(PatternRecognition.toSignificantDigits(toShortestString(self.absoluteUncertainty), self.significantDigits))

            s += ', ' + repr(self.uncertaintyType) + ', ' + repr(self.significantDigits) + ')'
        else:
            s += ')'

        return s

    def __hash__(self):
        return hash((self.value,self.unit.factor) + tuple(self.unit.powers))

    def _epsilonFactor(self, other):
        if isPhysicalQuantityWithUncertainty(other):
            return 2.0 * sys.float_info.epsilon * (abs(self.value) + abs(other.value))
        else:
            return 2.0 * sys.float_info.epsilon * (abs(self.value) + abs(float(other)))

    def _sum(self, other, sign1, sign2):
        if isPhysicalQuantityWithUncertainty(other):
            if other.unit.isDimensionless():
                if self.unit.isDimensionless():
                    new_value = sign1 * self.getValueAs('') + sign2 * other.getValueAs('')
                    result = self.__class__(new_value, self.unit)
                else:
                    raise TypeError('Incompatible units')
            else:
                new_value = sign1 * self.value + sign2 * other.value * other.unit.conversionFactorTo(self.unit)
                result = self.__class__(new_value, self.unit)
        else:
            if self.unit.isDimensionless():
                try:
                    new_value = float(sign1 * self.getValueAs('') + sign2 * other)
                    result = self.__class__(new_value, self.unit)
                except:
                    raise TypeError('Incompatible types')
            else:
                raise TypeError('Incompatible units')

        if abs(result.value) <= self._epsilonFactor(other):
            result = self.__class__(0, result.unit)

        if hasattr(self, 'uncertainty') and hasattr(other, 'uncertainty'):
            if self is other:
                # correlated
                if sign1 == sign2:
                    value = PatternRecognition.toSignificantDigits(toShortestString(result.value), self.significantDigits)
                    unc = PatternRecognition.toPrecision(toShortestString(self.absoluteUncertainty * 2), PatternRecognition.findPrecision(value), True)
                    result = PhysicalQuantityWithUncertainty(value, uncertainty = unc)
                else:
                    result = PhysicalQuantityWithUncertainty(result)
            else:
                # uncorrelated
                selfUnc = PhysicalQuantityWithUncertainty(self.absoluteUncertainty, self.unit)
                otherUnc = PhysicalQuantityWithUncertainty(other.absoluteUncertainty * other.unit.conversionFactorTo(self.unit), self.unit)

                selfSD = PatternRecognition.toSignificantDigits(toShortestString(self.value), self.significantDigits)
                otherSD = PatternRecognition.toSignificantDigits(toShortestString(other.value), other.significantDigits)
                selfP = PatternRecognition.findPrecision(selfSD)
                otherP = PatternRecognition.findPrecision(otherSD)
                leastPrecise = max(selfP, otherP)
                unc = math.sqrt(selfUnc.value**2 + otherUnc.value**2)
                result = PhysicalQuantityWithUncertainty(PatternRecognition.toPrecision(toShortestString(result.value), leastPrecise), uncertainty = PatternRecognition.toPrecision(toShortestString(unc), leastPrecise))

        return result

    def __add__(self, other):
        return self._sum(other, 1, 1)

    __radd__ = __add__

    def __sub__(self, other):
        return self._sum(other, 1, -1)

    def __rsub__(self, other):
        return self._sum(other, -1, 1)

    def __eq__(self, other):
        try:
            return self._sum(other, 1, -1).value == 0
        except TypeError:
            return False

    def __ne__(self, other):
        try:
            return self._sum(other, 1, -1).value != 0
        except TypeError:
            return True

    def __le__(self, other):
        return self._sum(other, 1, -1).value <= 0

    def __ge__(self, other):
        return self._sum(other, 1, -1).value >= 0

    def __lt__(self, other):
        return self._sum(other, 1, -1).value < 0

    def __gt__(self, other):
        return self._sum(other, 1, -1).value > 0

    """
    def __cmp__(self, other):
        diff = self._sum(other, 1, -1)
        return cmp(diff.value, 0)
    """

    def __mul__(self, other):
        if not isPhysicalQuantityWithUncertainty(other):
            return self.__class__(self.value*other, self.unit)

        value = self.value*other.value
        unit = self.unit*other.unit

        if unit.isDimensionless():
            return self.__class__(value*unit.factor, '')
        else:
            return self.__class__(value, unit)

    __rmul__ = __mul__

    def __div__(self, other):
        if not isPhysicalQuantityWithUncertainty(other):
            return self.__class__(self.value/other, self.unit)

        value = self.value/other.value
        unit = self.unit/other.unit

        if unit.isDimensionless():
            return self.__class__(value*unit.factor, '')
        else:
            return self.__class__(value, unit)

    def __rdiv__(self, other):
        if not isPhysicalQuantityWithUncertainty(other):
            return self.__class__(other/self.value, pow(self.unit, -1))

        value = other.value/self.value
        unit = other.unit/self.unit

        if unit.isDimensionless():
            return self.__class__(value*unit.factor, '')
        else:
            return self.__class__(value, unit)

    def __pow__(self, other):
        if isPhysicalQuantityWithUncertainty(other):
            raise TypeError('Exponents cannot be a PhysicalQuantityWithUncertainty')
        else:
            if self.unit.isDimensionless():
                return self.__class__(pow(self.value, other), '')
            else:
                return self.__class__(pow(self.value, other), pow(self.unit, other))

    def __rpow__(self, other):
        raise TypeError('Exponents cannot be a PhysicalQuantityWithUncertainty')

    def __abs__(self):
        return self.__class__(abs(self.value), self.unit)

    def __pos__(self):
        return self

    def __neg__(self):
        return self.__class__(-self.value, self.unit)

    def __bool__(self):
        return self.value != 0

    __nonzero__ = __bool__  # for python 2.x

    def convertToUnit(self, unit):
        """
        Change the unit and adjust the value such that the combination is equivalent to 
        the original one. The new unit must be compatible with the previous unit of the object.

        :param str unit: a unit
        :raises TypeError: if the unit string is not a know unit or a unit incompatible with the current one
        """

        unit = _findUnit(unit)
        self.value = _convertValue(self.value, self.unit, unit)
        self.unit = unit

    def inUnitsOf(self, *units):
        """
        Express the quantity in different units. If one unit is specified, a new 
        PhysicalQuantityWithUncertainty object is returned that expresses the quantity 
        in that unit. If several units are specified, the return value is a tuple of 
        PhysicalObject instances with with one element per unit such that the sum of all
        quantities in the tuple equals the the original quantity and all the values 
        except for the last one are integers. This is used to convert to irregular unit 
        systems like hour/minute/second.

        :param units: one or several units
        :type units: `str` or sequence of `str`
        :returns: one or more physical quantities
        :rtype: `PhysicalQuantityWithUncertainty` or `tuple` of `PhysicalQuantityWithUncertainty`
        :raises TypeError: if any of the specified units are not compatible with the original unit
        """

        units = list(map(_findUnit, units))

        if len(units) == 1:
            unit = units[0]
            value = _convertValue(self.value, self.unit, unit)
            return self.__class__(value, unit)
        else:
            units.sort()
            result = []
            value = self.value
            unit = self.unit

            for i in range(len(units) - 1, -1, -1):
                value = value*unit.conversionFactorTo(units[i])

                if i == 0:
                    rounded = value
                else:
                    rounded = _round(value)

                result.append(self.__class__(rounded, units[i]))
                value = value - rounded
                unit = units[i]

            return tuple(result)

    # Contributed by Berthold Hoellmann
    def inBaseUnits(self):
        """
        :returns: the same quantity converted to base units, i.e. SI units in most cases
        :rtype: `PhysicalQuantityWithUncertainty`
        """

        try:
            value = (self.value + self.unit.offset) * self.unit.factor
            num = ''
            denom = ''

            for i in range(9):
                unit = _base_symbols[i]
                power = self.unit.powers[i]

                if power < 0:
                    denom = denom + '/' + unit

                    if power < -1:
                        denom = denom + '**' + str(-power)

                elif power > 0:
                    num = num + '*' + unit

                    if power > 1:
                        num = num + '**' + str(power)

            if len(num) == 0:
                num = '1'
            else:
                num = num[1:]

            return self.__class__(value, num + denom)
        except:
            raise TypeError('Dimensionless unit has no base units')

    def isCompatible (self, unit):
        """
        :param unit: a unit
        :type unit: `str`
        :returns: `True` if the specified unit is compatible with the one of the quantity
        :rtype: `bool`
        """

        unit = _findUnit(unit)
        return self.unit.isCompatible(unit)

    def getValue(self):
        """
        Return value (float) of physical quantity (no unit).
        """

        return self.value

    def getUnitSymbol(self):
        """
        Return unit (string) of physical quantity.
        """

        return self.unit.symbol()

    def getValueAs(self, unit):
        """
        Returns value (float) in units of unit.
        """

        unit = _findUnit(unit)
        return _convertValue(self.value, self.unit, unit)

    def sqrt(self):
        return pow(self, 0.5)

    def sin(self):
        if self.unit.isAngle():
            return math.sin(self.value * self.unit.conversionFactorTo(_unit_table['rad']))
        else:
            raise TypeError('Argument of sin must be an angle')

    def cos(self):
        if self.unit.isAngle():
            return math.cos(self.value * self.unit.conversionFactorTo(_unit_table['rad']))
        else:
            raise TypeError('Argument of cos must be an angle')

    def tan(self):
        if self.unit.isAngle():
            return math.tan(self.value * self.unit.conversionFactorTo(_unit_table['rad']))
        else:
            raise TypeError('Argument of tan must be an angle')

    @staticmethod
    def isPhysicalUnitOf( self, unit ) :
        """
        Returns True if self has units as the argument unit and False otherwise. Self and unit can be an instance of string,
        PhysicalUnit or PhysicalQuantityWithUncertainty.
        """

        if( type( self ) == type( '' ) ) : self = _findUnit( self )
        selfUnit = _getPhysicalUnit( self )
        if( type( unit ) == type( '' ) ) : unit = _findUnit( unit )
        unitUnit = _getPhysicalUnit( unit )
        return( selfUnit.isCompatible( unitUnit ) )

    @staticmethod
    def isLength( self ) :
        "Returns True if self has units of length and False otherwise."

        return( PhysicalQuantityWithUncertainty.isPhysicalUnitOf( self, 'm' ) )

    @staticmethod
    def isMass( self ) :
        "Returns True if self has units of mass and False otherwise."

        return( PhysicalQuantityWithUncertainty.isPhysicalUnitOf( self, 'g' ) )

    @staticmethod
    def isTime( self ) :
        "Returns True if self has units of time and False otherwise."

        return( PhysicalQuantityWithUncertainty.isPhysicalUnitOf( self, 's' ) )

    @staticmethod
    def isCurrent( self ) :
        "Returns True if self has units of current and False otherwise."

        return( PhysicalQuantityWithUncertainty.isPhysicalUnitOf( self, 'A' ) )

    @staticmethod
    def isTemperature( self ) :
        "Returns True if self has units of temperature and False otherwise."

        return( PhysicalQuantityWithUncertainty.isPhysicalUnitOf( self, 'K' ) )

    @staticmethod
    def isMoles( self ) :
        "Returns True if self has units of mol and False otherwise."

        return( PhysicalQuantityWithUncertainty.isPhysicalUnitOf( self, 'mol' ) )

    @staticmethod
    def isLuminousIntensity( self ) :
        "Returns True if self has units of luminous intensity and False otherwise."

        return( PhysicalQuantityWithUncertainty.isPhysicalUnitOf( self, 'cd' ) )

    @staticmethod
    def isAngle( self ) :
        "Returns True if self has units of angle and False otherwise."

        return( PhysicalQuantityWithUncertainty.isPhysicalUnitOf( self, 'rad' ) )

    @staticmethod
    def isSolidAngle( self ) :
        "Returns True if self has units of solid angle and False otherwise."

        return( PhysicalQuantityWithUncertainty.isPhysicalUnitOf( self, 'sr' ) )

    @staticmethod
    def isEnergy( self ) :
        "Returns True if self has units of energy and False otherwise."

        return( PhysicalQuantityWithUncertainty.isPhysicalUnitOf( self, 'J' ) )

    @staticmethod
    def isSpin( self ) :
        "Returns True if self has units of spin and False otherwise."

        return( PhysicalQuantityWithUncertainty.isPhysicalUnitOf( self, 'hbar' ) )

    @staticmethod
    def isCharge( self ) :
        "Returns True if self has units of charge and False otherwise."

        return( PhysicalQuantityWithUncertainty.isPhysicalUnitOf( self, 'e' ) )

# Dictionary containing numbers
# These objects are meant to be used like arrays with generalized indices. Non-existent elements default to zero. Global operations are addition, subtraction, and multiplication/division by a scalar.
# Written by Konrad Hinsen <hinsen@cnrs-orleans.fr>
# Last revision: 2006-10-16
"""
DICTIONARY STORING NUMERICAL VALUES
"""


class NumberDict(dict):
    """
    .. rubric:: DICTIONARY STORING NUMERICAL VALUES

    Constructor: NumberDict()

    An instance of this class acts like an array of number with generalized 
    (non-integer) indices. A value of zero is assumed for undefined entries. 
    NumberDict instances support addition, and subtraction with other NumberDict 
    instances, and multiplication and division by scalars.
    """

    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            return 0

    def __coerce__(self, other):
        if type(other) == type({}):
            other = NumberDict(other)

        other.removeUnneededDimensionlessUnits()
        return self, other

    def __add__(self, other):
        sum_dict = NumberDict()

        for key in list(self.keys()):
            sum_dict[key] = self[key]

        for key in list(other.keys()):
            sum_dict[key] = sum_dict[key] + other[key]

        sum_dict.removeUnneededDimensionlessUnits()
        return sum_dict

    def __sub__(self, other):
        sum_dict = NumberDict()

        for key in list(self.keys()):
            sum_dict[key] = self[key]

        for key in list(other.keys()):
            sum_dict[key] = sum_dict[key] - other[key]

        sum_dict.removeUnneededDimensionlessUnits()
        return sum_dict

    def __mul__(self, other):
        new = NumberDict()

        for key in list(self.keys()):
            new[key] = other * self[key]

        new.removeUnneededDimensionlessUnits()
        return new

    __rmul__ = __mul__

    def __div__(self, other):
        new = NumberDict()

        for key in list(self.keys()):
            new[key] = self[key] / other

        new.removeUnneededDimensionlessUnits()
        return new

    def removeUnneededDimensionlessUnits(self):
        if (('' in self) and len(self)) > 1:
            del self['']


class PhysicalUnit:
    """
    .. rubric:: PHYSICAL UNIT

    A physical unit is defined by a symbol (possibly composite), a scaling factor, and
    the exponentials of each of the SI base units that enter into it. Units can be 
    multiplied, divided, and raised to integer powers.

    :param symbols: a dictionary mapping each symbol component to its associated integer
                power (e.g. ``{'m': 1, 's': -1}``) for `m/s`). As a shorthand, a string may be 
                passed which is assigned an implicit power 1.
    :type symbols: `dict` or `str`
    :param factor: a scaling factor
    :type factor: `float`
    :param powers: the integer powers for each of the nine base units
    :type powers: `list` of `int`
    :param offset: an additive offset to the base unit (used only for temperatures)
    :type offset: `float`
    """
    
    def __init__(self, symbols, factor, powers, offset=0):

        if symbols is not None:
            if '1' in symbols:
                rmOne = False

                for key in symbols:
                    if symbols[key] > 0:
                        rmOne = True

                if rmOne:
                    del symbols['1']

        if type(symbols) == type(''):
            self.symbols = NumberDict()
            self.symbols[symbols] = 1
        else:
            self.symbols = symbols

        self.factor = factor
        self.offset = offset
        self.powers = powers

    def __repr__(self):
        return self.__class__.__name__ + "('" + self.symbol() + "')"

    def __cmp__(self, other):
        raise DeprecationWarning("__cmp__ not supported in Python3, use __lt__ and other rich comparison operators")

        if self.powers != other.powers:
            raise TypeError('Incompatible units')

        return cmp(self.factor, other.factor)

    def __lt__(self,other):
        if self.powers != other.powers:
            raise TypeError('Incompatible units')

        return self.factor.__lt__(other.factor)

    def __eq__(self, other):
        if self.powers != other.powers:
            raise TypeError('Incompatible units')

        return self.factor.__eq__(other.factor)

    def __ne__(self, other):
        if self.powers != other.powers:
            raise TypeError('Incompatible units')

        return self.factor.__ne__(other.factor)

    def __mul__(self, other):
        if self.offset != 0 or (isPhysicalUnit (other) and other.offset != 0):
            raise TypeError("cannot multiply units with non-zero offset")

        if isPhysicalUnit(other):
            return PhysicalUnit(self.symbols + other.symbols, self.factor * other.factor, list(map(lambda a, b: a + b, self.powers, other.powers)))
        else:
            return PhysicalUnit(self.symbols + {str(other): 1}, self.factor * other, self.powers, self.offset * other)

    __rmul__ = __mul__

    def __truediv__(self, other):
        if self.offset != 0 or (isPhysicalUnit (other) and other.offset != 0):
            raise TypeError("cannot divide units with non-zero offset")

        if isPhysicalUnit(other):
            return PhysicalUnit(self.symbols - other.symbols, self.factor / other.factor, list(map(lambda a, b: a - b, self.powers, other.powers)))
        else:
            return PhysicalUnit(self.symbols + {str(other): -1}, self.factor / other, self.powers)

    def __rtruediv__(self, other):
        if self.offset != 0 or (isPhysicalUnit (other) and other.offset != 0):
            raise TypeError("cannot divide units with non-zero offset")

        if isPhysicalUnit(other):
            return PhysicalUnit(other.symbols - self.symbols, other.factor / self.factor, list(map(lambda a, b: a - b, other.powers, self.powers)))
        else:
            return PhysicalUnit(NumberDict({str(other): 1}) - self.symbols, other / self.factor, [-x for x in self.powers])

    __div__ = __truediv__   # for python 2.x
    __rdiv__ = __rtruediv__

    def __pow__(self, other):
        if self.offset != 0:
            raise TypeError("cannot exponentiate units with non-zero offset")

        if isinstance(other, int):
            return PhysicalUnit(other * self.symbols, pow(self.factor, other), list(map(lambda x, p = other: x * p, self.powers)))

        if isinstance(other, float):
            inv_exp = 1. / other
            rounded = int(math.floor(inv_exp + 0.5))

            if abs(inv_exp-rounded) < 1.e-10:
                if reduce(lambda a, b: a and b, list(map(lambda x, e = rounded: x % e == 0, self.powers))):
                    f = pow(self.factor, other)
                    p = list(map(lambda x, p = rounded: x / p, self.powers))

                    if reduce(lambda a, b: a and b, list(map(lambda x, e = rounded: x % e == 0, list(self.symbols.values())))):
                        symbols = self.symbols/rounded
                    else:
                        symbols = NumberDict()

                        if f != 1.:
                            symbols[str(f)] = 1

                        for i in range(len(p)):
                            symbols[_base_symbols[i]] = p[i]

                    return PhysicalUnit(symbols, f, p)
                else:
                    raise TypeError('Illegal exponent')

        raise TypeError('Only integer and inverse integer exponents allowed')

    def conversionFactorTo(self, other):
        """
        :param other: another unit
        :type other: `PhysicalUnit`
        :returns: the conversion factor from this unit to another unit
        :rtype: `float`
        :raises TypeError: if the units are not compatible
        """

        if self.powers != other.powers:
            raise TypeError('Incompatible units')

        if self.offset != other.offset and self.factor != other.factor:
            raise TypeError(('Unit conversion (%s to %s) cannot be expressed ' + 'as a simple multiplicative factor') % (self.symbol(), other.symbol()))

        return self.factor / other.factor

    def conversionTupleTo(self, other): # added 1998/09/29 GPW
        """
        :param other: another unit
        :type other: `PhysicalUnit`
        :returns: the conversion factor and offset from this unit to another unit
        :rtype: (`float`, `float`)
        :raises TypeError: if the units are not compatible
        """

        if self.powers != other.powers:
            raise TypeError('Incompatible units')

# let (s1,d1) be the conversion tuple from 'self' to base units (ie. (x+d1)*s1 converts a value x from 'self' to base units, and (x/s1)-d1 converts x from base to 'self' units) and (s2,d2) be the conversion tuple from 'other' to base units then we want to compute the conversion tuple (S,D) from 'self' to 'other' such that (x+D)*S converts x from 'self' units to 'other' units the formula to convert x from 'self' to 'other' units via the base units is (by definition of the conversion tuples):
#    ( ((x+d1)*s1) / s2 ) - d2
#  = ( (x+d1) * s1/s2) - d2
#  = ( (x+d1) * s1/s2 ) - (d2*s2/s1) * s1/s2
#  = ( (x+d1) - (d1*s2/s1) ) * s1/s2
#  = (x + d1 - d2*s2/s1) * s1/s2
# thus, D = d1 - d2*s2/s1 and S = s1/s2

        factor = self.factor / other.factor
        offset = self.offset - (other.offset * other.factor / self.factor)
        return (factor, offset)

    def isCompatible (self, other):     # added 1998/10/01 GPW
        """
        :param other: another unit
        :type other: `PhysicalUnit`
        :returns: `True` if the units are compatible, i.e. if the powers of the base units are the same
        :rtype: `bool`
        """
        return self.powers == other.powers

    def isDimensionless(self):
        return not reduce(lambda a, b: a or b, self.powers)

    def isAngle(self):
        return self.powers[7] == 1 and reduce(lambda a, b: a + b, self.powers) == 1

    def setSymbol(self, symbol):
        self.symbols = NumberDict()
        self.symbols[symbol] = 1

    def symbol(self):
        num = ''
        denom = ''

        for unit in list(self.symbols.keys()):
            power = self.symbols[unit]

            if power < 0:
                denom = denom + '/' + unit

                if power < -1:
                    denom = denom + '**' + str(-power)
            elif power > 0:
                num = num + '*' + unit

                if power > 1:
                    num = num + '**' + str(power)

        if len(num) == 0:
            num = '1'
        else:
            num = num[1:]

        return num + denom

    def __str__(self):
        return self.symbol()


# Type checks
def isPhysicalUnit(x):
    """
    :param x: an object
    :type x: any
    :returns: `True if x is a `PhysicalUnit`
    :rtype: `bool`
    """

    return hasattr(x, 'factor') and hasattr(x, 'powers')

def isPhysicalQuantityWithUncertainty(x):
    """
    :param x: an object
    :type x: any
    :returns: `True` if x is a `PhysicalQuantityWithUncertainty`
    :rtype: `bool`
    """

    return (hasattr(x, 'value') and hasattr(x, 'unit')) or (hasattr(x, 'value') and hasattr(x, 'uncertainty') and hasattr(x, 'uncertaintyType') and hasattr(x, 'significantDigits') and hasattr(x, 'unit'))


# Helper functions
def _findUnit(unit):
    if type(unit) == type(''):
        symbol = unit.strip()

        if symbol == '':
            unit = _unit_table[symbol]
        else:
            unit = eval(symbol, _unit_table)

        for cruft in ['__builtins__', '__args__']:
            try:
                del _unit_table[cruft]
            except:
                pass

    if not isPhysicalUnit(unit):
        raise TypeError(str(unit) + ' is not a unit')

    return unit

def _round(x):
    if x > 0.:
        return math.floor(x)
    else:
        return math.ceil(x)

def _convertValue(value, src_unit, target_unit):

    try :
        (factor, offset) = src_unit.conversionTupleTo(target_unit)
    except TypeError :
        raise TypeError( 'Unit "%s" not convertible to "%s"' % ( src_unit, target_unit ) )
    except :
        raise
    return (value + offset) * factor

def _getPhysicalUnit( self ) :

    if( isinstance( self, PhysicalUnit ) ) : return( self )
    if( isinstance( self, PhysicalQuantityWithUncertainty ) ) : return( self.unit )
    raise Exception( 'Valid instance for _getPhysicalUnit' )

# SI unit definitions

_base_symbols = ['m', 'kg', 's', 'A', 'K', 'mol', 'cd', 'rad', 'sr']

_base_units = [('m',   PhysicalUnit('m',   1.,    [1,0,0,0,0,0,0,0,0])),
               ('g',   PhysicalUnit('g',   0.001, [0,1,0,0,0,0,0,0,0])),
               ('s',   PhysicalUnit('s',   1.,    [0,0,1,0,0,0,0,0,0])),
               ('A',   PhysicalUnit('A',   1.,    [0,0,0,1,0,0,0,0,0])),
               ('K',   PhysicalUnit('K',   1.,    [0,0,0,0,1,0,0,0,0])),
               ('mol', PhysicalUnit('mol', 1.,    [0,0,0,0,0,1,0,0,0])),
               ('cd',  PhysicalUnit('cd',  1.,    [0,0,0,0,0,0,1,0,0])),
               ('rad', PhysicalUnit('rad', 1.,    [0,0,0,0,0,0,0,1,0])),
               ('sr',  PhysicalUnit('sr',  1.,    [0,0,0,0,0,0,0,0,1])),
               ]

_prefixes = [('Y',  1.e24, 'yotta'),
             ('Z',  1.e21, 'zetta'),
             ('E',  1.e18, 'exa'),
             ('P',  1.e15, 'peta'),
             ('T',  1.e12, 'tera'),
             ('G',  1.e9, 'giga'),
             ('M',  1.e6, 'mega'),
             ('k',  1.e3, 'kilo'),
             ('h',  1.e2, 'hecto'),
             ('da', 1.e1, 'deka'),
             ('d',  1.e-1, 'deci'),
             ('c',  1.e-2, 'centi'),
             ('m',  1.e-3, 'milli'),
             ('mu', 1.e-6, 'micro'),
             ('n',  1.e-9, 'nano'),
             ('p',  1.e-12, 'pico'),
             ('f',  1.e-15, 'femto'),
             ('a',  1.e-18, 'atto'),
             ('z',  1.e-21, 'zepto'),
             ('y',  1.e-24, 'yocto'),
             ]

_help = []
_unit_table = {}

for unit in _base_units:
    _unit_table[unit[0]] = unit[1]

def _addUnit(symbol, unit, name=''):
    if symbol in _unit_table:
        raise KeyError('Unit ' + symbol + ' already defined')

    if name:
        _help.append((name, unit, symbol))

    if type(unit) == type(''):
        unit = eval(unit, _unit_table)

        for cruft in ['__builtins__', '__args__']:
            try:
                del _unit_table[cruft]
            except:
                pass

    unit.setSymbol(symbol)
    _unit_table[symbol] = unit

def _addPrefixedUnit(unit):
    _prefixed_symbols = []

    for prefix in _prefixes:
        symbol = prefix[0] + unit
        _addUnit(symbol, prefix[1]*_unit_table[unit])
        _prefixed_symbols.append(symbol)

    for entry in _help:
        if isinstance(entry, tuple):
            if len(entry) == 3:
                name0, unit0, symbol0 = entry
                i = _help.index(entry)
                if unit == symbol0:
                    _help.__setitem__(i, (name0, unit0, symbol0, _prefixed_symbols))

# SI derived units; these automatically get prefixes
_help.append('.. rubric:: Prefixes: \n \n::')
_help.append('    ----  ------ ------\n' + '    name  factor symbol\n' + '    ----  ------ ------')
_help_prefix = []

#for i in range(20):
#   _help_prefix.append('%-5s %-6s %-.0e' % (_prefixes[i][2], _prefixes[i][0], _prefixes[i][1]))

#_help.append('\n'.join(_help_prefix))

for i in range(20):
    _help.append(('%-5s %-6.0e %-s' % (_prefixes[i][2], _prefixes[i][1], _prefixes[i][0]), '', ''))

_help.append('\n    %-26s %-32s %-s\n    %-26s %-32s %-s\n    %-26s %-32s %-s' % ('--------------------------', '--------------------------------', '------', 'name', 'value', 'symbol', '--------------------------', '--------------------------------', '------'))
_help.append('    SI derived units; these automatically get prefixes ............... ')

_unit_table['kg'] = PhysicalUnit('kg',   1., [0,1,0,0,0,0,0,0,0])
_addUnit('Hz', '1./s', 'Hertz')
_addUnit('N', 'm*kg/s**2', 'Newton')
_addUnit('Pa', 'N/m**2', 'Pascal')
_addUnit('J', 'N*m', 'Joule')
_addUnit('W', 'J/s', 'Watt')
_addUnit('C', 's*A', 'Coulomb')
_addUnit('V', 'W/A', 'Volt')
_addUnit('F', 'C/V', 'Farad')
_addUnit('ohm', 'V/A', 'Ohm')
_addUnit('S', 'A/V', 'Siemens')
_addUnit('Wb', 'V*s', 'Weber')
_addUnit('T', 'Wb/m**2', 'Tesla')
_addUnit('H', 'Wb/A', 'Henry')
_addUnit('lm', 'cd*sr', 'Lumen')
_addUnit('lx', 'lm/m**2', 'Lux')
_addUnit('Bq', '1./s', 'Becquerel')
_addUnit('Gy', 'J/kg', 'Gray')
_addUnit('Sv', 'J/kg', 'Sievert')

del _unit_table['kg']

for unit in list(_unit_table.keys()):
    _addPrefixedUnit(unit)

# NOTE everything below does not have prefixes added to units

# Dimensionless quantities
_unit_table[''] = PhysicalUnit('', 1., [0,0,0,0,0,0,0,0,0])

# Fundamental constants
_help.append('    Fundamental constants ............................................ ')

_unit_table['pi'] = math.pi
_addUnit('c', '299792458.*m/s', 'speed of light')
_addUnit('mu0', '4.e-7*pi*N/A**2', 'permeability of vacuum')
_addUnit('eps0', '1./mu0/c**2', 'permittivity of vacuum')
_addUnit('Grav', '6.67259e-11*m**3/kg/s**2', 'gravitational constant')
_addUnit('hplanck', '6.6260755e-34*J*s', 'Planck constant')
_addUnit('hbar', 'hplanck/(2*pi)', 'Planck constant / 2pi')
_addUnit('e', '1.60217733e-19*C', 'elementary charge')
_addUnit('me', '9.1093897e-31*kg', 'electron mass')
_addUnit('mp', '1.6726231e-27*kg', 'proton mass')
_addUnit('Nav', '6.0221367e23/mol', 'Avogadro number')
_addUnit('k', '1.380658e-23*J/K', 'Boltzmann constant')

# Time units
_help.append('    Time units ....................................................... ')

_addUnit('min', '60.*s', 'minute')
_addUnit('h', '60.*min', 'hour')
_addUnit('d', '24.*h', 'day')
_addUnit('wk', '7.*d', 'week')
_addUnit('yr', '365.25*d', 'year')

# Length units
_help.append('    Length units ..................................................... ')

_addUnit('inch', '2.54*cm', 'inch')
_addUnit('ft', '12.*inch', 'foot')
_addUnit('yd', '3.*ft', 'yard')
_addUnit('mi', '5280.*ft', '(British) mile')
_addUnit('nmi', '1852.*m', 'Nautical mile')
_addUnit('Ang', '1.e-10*m', 'Angstrom')
_addUnit('lyr', 'c*yr', 'light year')
_addUnit('Bohr', '4.*pi*eps0*hbar**2/me/e**2', 'Bohr radius')

# Area units
_help.append('    Area units ....................................................... ')

_addUnit('ha', '10000*m**2', 'hectare')
_addUnit('acres', 'mi**2/640', 'acre')
_addUnit('b', '1.e-28*m**2', 'barn')

_addPrefixedUnit('b')

# Volume units
_help.append('    Volume units ..................................................... ')

_addUnit('l', 'dm**3', 'liter')
_addUnit('dl', '0.1*l', 'deci liter')
_addUnit('cl', '0.01*l', 'centi liter')
_addUnit('ml', '0.001*l', 'milli liter')
_addUnit('tsp', '4.92892159375*ml', 'teaspoon')
_addUnit('tbsp', '3.*tsp', 'tablespoon')
_addUnit('floz', '2.*tbsp', 'fluid ounce')
_addUnit('cup', '8.*floz', 'cup')
_addUnit('pt', '16.*floz', 'pint')
_addUnit('qt', '2.*pt', 'quart')
_addUnit('galUS', '4.*qt', 'US gallon')
_addUnit('galUK', '4.54609*l', 'British gallon')

# Mass units
_help.append('    Mass units ....................................................... ')

_addUnit('amu', '1.6605402e-27*kg', 'atomic mass units')
_addUnit('oz', '28.349523125*g', 'ounce')
_addUnit('lb', '16.*oz', 'pound')
_addUnit('ton', '2000.*lb', 'ton')

# Force units
_help.append('    Force units ...................................................... ')

_addUnit('dyn', '1.e-5*N', 'dyne (cgs unit)')

# Energy units
_help.append('    Energy units ..................................................... ')

_addUnit('erg', '1.e-7*J', 'erg (cgs unit)')
_addUnit('eV', 'e*V', 'electron volt')
_addUnit('Hartree', 'me*e**4/16/pi**2/eps0**2/hbar**2', 'Wavenumbers/inverse cm')
_addUnit('Ken', 'k*K', 'Kelvin as energy unit')
_addUnit('cal', '4.184*J', 'thermochemical calorie')
_addUnit('kcal', '1000.*cal', 'thermochemical kilocalorie')
_addUnit('cali', '4.1868*J', 'international calorie')
_addUnit('kcali', '1000.*cali', 'international kilocalorie')
_addUnit('Btu', '1055.05585262*J', 'British thermal unit')

_addPrefixedUnit('eV')

# Power units
_help.append('    Power units ...................................................... ')

_addUnit('hp', '745.7*W', 'horsepower')

# Pressure units
_help.append('    Pressure units ................................................... ')

_addUnit('bar', '1.e5*Pa', 'bar (cgs unit)')
_addUnit('atm', '101325.*Pa', 'standard atmosphere')
_addUnit('torr', 'atm/760', 'torr = mm of mercury')
_addUnit('psi', '6894.75729317*Pa', 'pounds per square inch')

# Angle units
_help.append('    Angle units ...................................................... ')

_addUnit('deg', 'pi*rad/180.', 'degrees')

_help.append('    Temperature units ................................................ ')
# Temperature units -- can't use the 'eval' trick that _addUnit provides
# for degC and degF because you can't add units
Kelvin = _findUnit('K')
_addUnit('degR', '(5./9.)*K', 'degrees Rankine')
_addUnit('degC', PhysicalUnit(None, 1.0, Kelvin.powers, 273.15),
          'degrees Celcius')
_addUnit('degF', PhysicalUnit(None, 5./9., Kelvin.powers, 459.67),
          'degree Fahrenheit')
          
_help.append('\n')

del Kelvin

def description():
    """
    Return a string describing all available units.
    """

    s = [] # collector for description text

    for entry in _help:
        if isinstance(entry, str):
            if entry != None:
#               s.append('\n%-s\n' % ('-' * len(entry)) + entry + '\n%-s\n' % ('-' * len(entry)))
                s.append('\n'  + entry + '\n')
            else:
                pass
        elif isinstance(entry, tuple):
            if len(entry) == 4:
                name, unit, symbol, prefixes = entry
            elif len(entry) == 3:
                name, unit, symbol = entry

            s.append('    %-26s %-32s %-s\n' % (name, unit, symbol))
        else:
            # impossible
            raise TypeError('wrong construction of _help list')

    entry = '    ------------- ------ --------------\n    %-13s %-6s prefixed units\n    ------------- ------ --------------' % ('name', 'symbol')
#   s.append('\n%-s\n' % ('-' * len(entry)) + entry + '\n%-s\n' % ('-' * len(entry)))
    s.append('\n'  + entry + '\n')

    for entry in _help:
        if isinstance(entry, str):
            pass
        elif isinstance(entry, tuple):
            if len(entry) == 4:
                name, unit, symbol, prefixes = entry
                s.append('    %-13s %-6s %-s\n    %-13s %-6s %-s\n\n' % (name, symbol, ', '.join(prefixes[0:10]), ' ', ' ', ', '.join(prefixes[10:20])))

    return ''.join(s)

# add the description of the units to the module's doc string:
__doc__ += '\n' + description()

# Some demonstration code. Run with "python -i PhysicalQuantityWithUncertainty.py"
# to have this available.

if __name__ == '__main__':

    l = PhysicalQuantityWithUncertainty(10., 'm')
    big_l = PhysicalQuantityWithUncertainty(10., 'km')
    print big_l + l

    t = PhysicalQuantityWithUncertainty(314159., 's')
    print t.inUnitsOf('d', 'h', 'min', 's')

    pqu = PhysicalQuantityWithUncertainty # just a shorthand...

    e = pqu('2.7 Hartree * Nav')
    e.convertToUnit('kcal/mol')
    print e
    print e.inBaseUnits()

    freeze = pqu('0 degC')
    print freeze.inUnitsOf('degF')

    E = pqu(1000, 'MeV/c**2')
    print E.inUnitsOf('kg')

    # This is the similar to the last definition for getValueAs because BRB
    E = pqu(1000, 'MeV/c**2') 
    amu = E.getValueAs('amu') # coded it to handle masses and energy as the same thing when needed.

    E2 = pqu(pqu(amu, 'amu').getValueAs('MeV/c**2'), 'MeV/c**2')
    print E, amu, 'amu', E2, E.getValueAs('kg'), 'kg', pqu(1, 'kg').getValueAs('eV/c**2'), 'eV/c**2'
