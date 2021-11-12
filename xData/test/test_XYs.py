# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import unittest

import math
from numpy import logspace, flip, linspace

import xData.unittest
from fudge.reactionData.crossSection import upperEps
from pqu import PQU
from xData import XYs

DOPLOTS=False


def equal_lethargy_bins(numBins, domainMin=1e-5, domainMax=20.0e6, reverse=False):
    seq = logspace(start=math.log10(domainMin), stop=math.log10(domainMax), num=numBins)
    if reverse:
        seq = flip(seq)
    return seq


def function_to_XYs(func, fpars,
                    Egrid=equal_lethargy_bins(1000),
                    domainUnit='eV', domainName='energy_in', rangeUnit='1/eV', rangeName='flux',
                    accuracy=upperEps):
    """
    Helper function to convert a user-created function (OK, one of the spectra below) into an XYs instance
    that can be integrated, grouped, whatever.  We pre-defined a energy grid (in eV) that should work well
    even for pathological "spectra" like the problematic 1/E for the resonance integral.
    """
    return XYs.XYs1d.createFromFunction(
        XYs.XYs1d.defaultAxes(labelsUnits={
            XYs.yAxisIndex: (rangeName, rangeUnit),
            XYs.xAxisIndex: (domainName, domainUnit)}),
        Xs=Egrid,
        func=func,
        parameters=fpars,
        accuracy=accuracy,
        biSectionMax=20,
        checkForRoots=False,
        infill=1,
        safeDivide=1)


def simple_function_to_XYs(func, domain, nsteps=100):
    return function_to_XYs(lambda x, y: func(x), {},
                    Egrid=linspace(start=domain[0], stop=domain[1], num=nsteps),
                    domainUnit='', domainName='energy_in', rangeUnit='', rangeName='crossSection',
                    accuracy=upperEps)


def DEFAULTONEOVEREFUNCTION(E, *args):
    if E < float(CADMIUMCUTTOFF.value):
        return 0.0
    return 1.0 / E


CADMIUMCUTTOFF = PQU.PQU(0.5, 'eV')
DEFAULTONEOVEREGRID = [1e-5,
                       0.99999 * float(CADMIUMCUTTOFF.value)] + \
                      list(equal_lethargy_bins(2000, domainMin=float(CADMIUMCUTTOFF.value)))
DEFAULTONEOVEREXYs = function_to_XYs(DEFAULTONEOVEREFUNCTION, [], Egrid=DEFAULTONEOVEREGRID)
epsilon = 1e-6


class Test_XYs(xData.unittest.TestCaseWithIsClose):

    def test_function_to_XYs(self):
        def oneOverEFunc(E, *args):
            return 1.0 / E

        x = function_to_XYs(oneOverEFunc, [])
        for E in [1e-4, 2.3e-3, 33e-2, 0.1, 19.0]:
            self.assertIsClose(x.evaluate(E), 1.0 / E)

        def myExp(_E, *args):
            return math.exp(-_E / 1e3)

        y = function_to_XYs(myExp, [])
        for E in [33e-2, 0.1, 20.]:
            self.assertIsClose(y.evaluate(E), math.exp(-E / 1e3))

    def test_grouped_values_to_XYs(self):
        self.assertIsClose(0.0, 0.0)

    def test_integrateTwoFunctions(self):
        def myExp1(E, *args):
            return math.exp(-E)

        def myExp2(E, *args):
            return math.exp(-E / 2.0)

        a = function_to_XYs(myExp1, []).integrateTwoFunctions(function_to_XYs(myExp2, [])).inUnitsOf('1/eV')
        b = PQU.PQU(0.666657, '1/eV')
        self.assertIsClose(a, b, absTol=5e-7)

    def test_one_over_E(self):
        """
        Mathematica says the answer is \int_0.5^20e6 dE/E = 17.5044
        :return:
        """

        def oneFunc(x, *args): return 1.0

        one = function_to_XYs(oneFunc, [])
        self.assertIsClose(DEFAULTONEOVEREXYs.integrateTwoFunctions(one), PQU.PQU(17.5044, '1/eV'), percent=1e-4)

    def run_pdfOfY_test(self, testFunc, testFuncDomain, answer, ny=10, doPlot=False, verbose=False):
        """
        This is a common test runner for all of the pdfOfY tests.  In each case, we have an analytic result for the
        pdfOfY (`answer`) as well as the test function (`testFunc`) from which the analytic result is derived.

        Note:

            * y=testFunc(x), so domain in x variable is irrelevant, range in y is what is important

            * z=pdfOfY(y), range in y becomes the domain in y for the pdf, the range of the pdf is [0,1]

        :param testFunc: test function from which the pdfOfY is computed numnerically by FUDGE
        :param testFuncDomain: domain of the testFunc
        :param answer:  the analytic result for the pdfOfY
        :param ny:  number of steps in the y direction to loop over when testing values
        :param doPlot: do the plots
        :param verbose: turn on verbosity
        :return: None
        """

        # Contruct the numerical result for the pdfOfY from the test function testFunc
        testFuncAsXYs1d = XYs.XYs1d(data=simple_function_to_XYs(testFunc, domain=testFuncDomain, nsteps=2))
        testFuncPdfOfY = testFuncAsXYs1d.pdfOfY(epsilon=1e-8)

        # Plots, for debugging (if they currently work!)
        if doPlot:
            from fudge.vis.gnuplot import fudgeMultiPlots as fudgeMultiPlotsModule
            fudgeMultiPlotsModule.multiPlot([answer, testFuncPdfOfY])

        # Check the norm of the answer and the original.  Since the result is a pdf,
        # the norms had better both be 1.0, so this is a bit of a silly test.
        if verbose:
            print('norm:', answer.integrate(), testFuncPdfOfY.integrate())
        self.assertIsClose(answer.integrate(), testFuncPdfOfY.integrate())
        self.assertIsClose(1.0, testFuncPdfOfY.integrate())

        # Another silly test to make sure the final range in probability space is correct
        if verbose:
            print('range function:', testFuncPdfOfY.range(), answer.range())
            print('range min/max:', [testFuncPdfOfY.rangeMin, testFuncPdfOfY.rangeMax], [answer.rangeMin, answer.rangeMax])
        if False:
            self.assertListAllAreClose([answer.rangeMin, answer.rangeMax],
                                       [testFuncPdfOfY.rangeMin, testFuncPdfOfY.rangeMax])

        if verbose:
            print('domain:', testFuncPdfOfY.domain(), answer.domain())
        self.assertListAllAreClose(answer.domain(), testFuncPdfOfY.domain())

        for i in range(ny):
            dy = 1.0/ny
            y = dy*float(i)
            a = answer.evaluate(y)
            b = testFuncPdfOfY.evaluate(y)
            if a is None and b is None:
                continue
            if verbose:
                print(y, a, b)
            self.assertIsClose(a, b)

    @unittest.skip("Delta function test will never run correctly")
    def test_pdfOfY_test1(self):
        """
        Test 1
        ------
        A constant!  Should be easy.
            * f1(x) = 1
            * domain is x in [0,1]
            * The pdfOfY should be P1(y) = delta(y-1).
            * The plot of pdfOfY should be a very spiky triangle centered at 1.

        We skip this test because there is no way FUDGE can get this right except in the most approximate way
        with a very pointy triangle.
        """

        def f1(_x):
            return 1

        answer = XYs.XYs1d(data=[[1.0 - epsilon, 0.0], [1.0, 1.0/epsilon], [1.0 + epsilon, 0.0]])
        self.run_pdfOfY_test(testFunc=f1, answer=answer, testFuncDomain=(0.0, 1.0), ny=3, doPlot=DOPLOTS)

    def test_pdfOfY_test2(self):
        """
        Test 2
        ------
        Line with slope 1.  Should be easy!
            * f2(x) = x
            * domain is x in [0,1]
            * The pdfOfY should be P2(y) = 1.0 for y in domain [0,1]
        """

        def f2(_x):
            return _x

        answer = XYs.XYs1d(data=[[0.0, 1.0], [2.0, 1.0]])
        answer = answer.normalize()
        self.run_pdfOfY_test(testFunc=f2, answer=answer, testFuncDomain=(0.0, 2.0), doPlot=DOPLOTS, verbose=False)

    def test_pdfOfY_test3(self):
        """
        Test 3
        ------
        The hat function!  Might be easy, not sure.
            * f3(x) = x   if x < 1.0
                    = 2-x if x > 1.0
            * domain is x in [0,2]
            * The pdfOfY should be P3(y) = 1.0 for y in domain [0,1] as each
              segment independantly has a slope of +/-1, but otherwise below
              1.0 is basically Test 2 with weight 1/2 and above 1.0 is
              basically Test 2 with weight 1/2.  Should get `1`.
        """

        def f3(x):
            if x <= 1:
                return x
            return 2.0 - x

        testFuncDomain = (0.0, 2.0)
        answer = XYs.XYs1d(data=[[0.0, 1.0], [1.0, 1.0]])
        self.run_pdfOfY_test(testFunc=f3, answer=answer, testFuncDomain=testFuncDomain, doPlot=DOPLOTS)

    @unittest.skip("broken")
    def test_pdfOfY_test4(self):
        """
        Test 4
        ------
        The sine function.
            * f4(x) = sin( pi * x )
            * The pdfOfY should be P4(y) = (2/pi) / sqrt( 1 - y^2 )

        DAB 7 Sep 2018: FAILED!  looks like fluctuations around correct function.
        """

        def f4(x):
            return math.sin(math.pi * x)

        def f4pdfOfY(y):
            if abs(y) == 1.0:
                return 0.0
            return 1.0 / (math.pi * math.sqrt(1.0 - y ** 2))

        testFuncDomain = (-1.0, 1.0)
        answer = simple_function_to_XYs(f4pdfOfY, domain=(-1.0, 1.0), nsteps=100000)
        self.run_pdfOfY_test(testFunc=f4, answer=answer, testFuncDomain=testFuncDomain, doPlot=DOPLOTS, verbose=True)

    def test_pdfOfY_test5(self):
        """
        Test 5
        ------
        A truncated quadratic
            * f6(x) = 1-x^2 for |x| <= 1, otw. 0
            * The pdfOfY should be P6(y) = 0.5 / sqrt( 1 - y )

        """

        def f5(x):
            if abs(x) > 1.0:
                return 0.0
            return 1.0 - x * x

        def f5pdfOfY(y):
            if y == 1.0:
                return 0.0  # handle singularity
            return 0.5 / math.sqrt(1.0 - y)

        pdfDomain = (0.0, 1.0)
        testFuncDomain = (-1.0, 1.0)
        answer = simple_function_to_XYs(f5pdfOfY, domain=pdfDomain, nsteps=10000)
        answer.normalize()
        self.run_pdfOfY_test(testFunc=f5, answer=answer, testFuncDomain=testFuncDomain, doPlot=DOPLOTS, verbose=False)

    def test_pdfOfY_test6(self):
        """
        Test 6
        ------
        A Lorenzian centered at 1 with width 1 function.
            * f6(x) = 1 / ( ( x-1 )^2 + 1/4 )
            * The pdfOfY should be P6(y) = N / ( y^(3/2) * sqrt( 4 - y ) )

        DAB 7 Sep 2018: FAILED!  Again, funky steppy thing bracketing the correct function.
                                 but this time, my "answer" has a normalization factor that
                                 I didn't bother to work out.
        """

        def f6(x):
            x1 = x - 1.0
            return 1.0 / (x1 * x1 + 0.25)

        def f6pdfOfY(y):
            if (y > 3.999999) or (y < 1e-6):
                return 0.0
            return 2 * math.sqrt(1.0 / 1599.0) / (math.pow(y, 1.5) * math.sqrt(4.0 - y))

        testFuncDomain = (-19.0, 21.0)
        answer = simple_function_to_XYs(f6pdfOfY, domain=(1.0 / 400.25, 4.0), nsteps=100000)
        self.run_pdfOfY_test(testFunc=f6, answer=answer, testFuncDomain=testFuncDomain, doPlot=DOPLOTS, verbose=False)


# -------------------------------------------------------------------------------
# Run unit tests
# -------------------------------------------------------------------------------
if __name__ == "__main__":

    unittest.main()
