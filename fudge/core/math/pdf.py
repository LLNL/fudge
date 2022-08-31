# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import math
import numpy

from pqu.PQU import PhysicalQuantityWithUncertainty as PQU
from xData import XYs1d

SQRT2 = math.sqrt(2.0)
SQRTTWOOVERPI = math.sqrt(2.0 / math.pi)
SQRTTWOPI = math.sqrt(2.0 * math.pi)

# ---------------------------------------------------------------------------------
#
#    Probability Distribution Function helper classes
#
# ---------------------------------------------------------------------------------
import abc
import xData.base as baseModule


class XYs1dFunctional(baseModule.XDataFunctional):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def evaluate(self, x): pass

    @property
    @abc.abstractmethod
    def domainMin(self): pass

    @property
    @abc.abstractmethod
    def domainMax(self): pass

    @property
    @abc.abstractmethod
    def domainUnit(self): pass

    @abc.abstractmethod
    def toPointwise_withLinearXYs(self, **kwargs): pass


class PDF1d(XYs1dFunctional):
    """
    Abstract base class 1d probability density functions.
    """
    __metaclass__ = abc.ABCMeta

    @property
    @abc.abstractmethod
    def mean(self): pass

    @property
    @abc.abstractmethod
    def variance(self): pass

    @property
    def stddev(self): return math.sqrt(self.variance)

    @property
    @abc.abstractmethod
    def skew(self): pass

    @property
    @abc.abstractmethod
    def kurtosis(self): pass

    @property
    def median(self): raise NotImplementedError("WRITEME")

    @property
    def mode(self): raise NotImplementedError("WRITEME")

    @property
    def entropy(self): raise NotImplementedError("WRITEME")

    @property
    @abc.abstractmethod
    def confidenceInterval(self, conf): pass

    @property
    @abc.abstractmethod
    def cdf(self): pass

    @property
    @abc.abstractmethod
    def quantileFunction(self): pass

    @property
    def invcdf(self): return self.quantileFunction

    @property
    @abc.abstractmethod
    def drawSample(self, size=None): pass


class UnivariatePDF(XYs1d.XYs1d, PDF1d):
    """
    Simple implmentation of a generic probability distribution function.

    Relies on stuff build into fudge's XYs1d class.  In other words,  it takes the real PDF you
    define in the getValueRawPDF (and optionally getValueRawCDF) and uses it to populate
    the fudge XYs1d parent class info.  Then it uses that to fill out all the stuff needed to make
    this UnivariantePDF a valid PDF1d.
    """
    defaultDomainMin = -10.0
    defaultDomainMax = 10.0
    defaultSupportMin = 0.0
    defaultSupportMax = float("Inf")

    def __init__(self, xUnit='', domainMin=None, domainMax=None, supportMin=None, supportMax=None,
                 numXIntervals=100, xLabel="indep. variable", yLabel="PDF", accuracy=1e-8, params={}, **kwds):
        """
        Base class

        :param xUnit: units to use for the x axis
        :param domainMin: (float or PQU) pdf and cdf defined on interval (domainMin, domainMax), if domainMin a PQU,
                          will convert to xUnit, otherwise send in a float
        :param domainMax: (float or PQU) pdf and cdf defined on interval (domainMin, domainMax), if domainMax a PQU,
                          will convert to xUnit, otherwise send in a float
        :param numXIntervals: (int) number of intervals to subdivide the pdf into when turning it into an XYs1d linear
                              interpolation table (default=100)
        :param xLabel: (str) label for x axis (default "indep. variable")
        :param yLabel: (str) lable for y axis (default "PDF")
        :param accuracy: (float) relative accuracy to refine to when building the interpolation table from the
                         analytic form of the pdf
        :param params: (dict) other parameters needed by the analytic form of the pdf
        :param kwds: unused
        :return: the instance
        """

        # setup support interval for the PDF; this is the range on real axis for which the PDF is not zero
        self.supportMin, self.supportMax = self.fixSupportValues(supportMin, supportMax, xUnit)

        # Set up unit of PDF.  It is a density, so has units of 1/xUnit
        if xUnit == "":
            yUnit = ""
        else:
            yUnit = "1/(" + xUnit + ")"

        # Set up the grid in x, we don't need to save the domainMin and domainMax in self,
        # that's handled by the XYs1d constructor
        if 'xGrid' in kwds:
            xGrid = kwds['xGrid']
            domainMin = xGrid[0]
            domainMax = xGrid[-1]
        else:
            domainMin, domainMax = self.fixDomainValues(domainMin, domainMax, xUnit)
            xGrid = list(numpy.linspace(domainMin, domainMax, numXIntervals + 1))

        # Construct interpolation table of the PDF
        if 'pdfTable' in kwds:
            XYs1d.XYs1d.__init__(self,
                               axes=XYs1d.XYs1d.defaultAxes(labelsUnits={1: (xLabel, xUnit), 0: (yLabel, yUnit)}),
                               data=[xGrid, list(kwds['pdfTable']) + list(kwds['pdfTable'][-1:])],
                               interpolation='flat', dataForm='xsandys', **kwds)
        elif 'XYs1d' in kwds:
            XYs1d.XYs1d.__init__(self,
                               axes=kwds['XYs1d'].axes,
                               data=kwds['XYs1d'],
                               interpolation=kwds['XYs1d'].interpolation, index=kwds['XYs1d'].index)
        elif kwds:
            XYs1d.XYs1d.__init__(self, **kwds)
        else:
            tmp = XYs1d.XYs1d.createFromFunction(
                axes=XYs1d.XYs1d.defaultAxes(labelsUnits={1: (xLabel, xUnit), 0: (yLabel, yUnit)}),
                Xs=xGrid,
                func=self.__class__.getValueRawPDF,
                parameters=params,
                accuracy=accuracy,
                biSectionMax=50)
            XYs1d.XYs1d.__init__(self, axes=tmp.axes, data=tmp.copyDataToXsAndYs(), dataForm='xsandys', **kwds)

        # Compute the CDF
        self.__cdf = self.getCDF(params, xGrid, xLabel, xUnit)

    def fixDomainValues(self, domainMin, domainMax, xUnit):
        if isinstance(domainMin, PQU):
            domainMin = float(domainMin.inUnitsOf(xUnit).value)
        if isinstance(domainMax, PQU):
            domainMax = float(domainMax.inUnitsOf(xUnit).value)
        if domainMin is None:
            domainMin = self.__class__.defaultDomainMin
        if domainMax is None:
            domainMax = self.__class__.defaultDomainMax
        if domainMax <= domainMin:
            raise ValueError("domainMax=%g <= domainMin=%g" % (domainMax, domainMin))
        return domainMin, domainMax

    def fixSupportValues(self, supportMin, supportMax, xUnit):
        if isinstance(supportMin, PQU):
            supportMin = float(supportMin.inUnitsOf(xUnit).value)
        if isinstance(supportMax, PQU):
            supportMax = float(supportMax.inUnitsOf(xUnit).value)
        if supportMin is None:
            supportMin = self.__class__.defaultSupportMin
        if supportMax is None:
            supportMax = self.__class__.defaultSupportMax
        if supportMax <= supportMin:
            raise ValueError("supportMax=%g <= supportMin=%g" % (supportMax, supportMin))
        return supportMin, supportMax

    @staticmethod
    def getValueRawPDF(x, params):
        """
        Override this function to define the analytic form of the PDF.
        This is used to initialize self.

        :param x: dimensionless input parameter (float)
        :param params: other dimensionless parameters needed by the PDF
        :return: the value of the PDF evaluated at x (float)
        """
        raise NotImplementedError("Override this in derived classes")

    @staticmethod
    def getValueRawCDF(x, params):
        """
        Override this function to define the analytic form of the CDF.
        This is used to initialize self.cdf

        :param x: dimensionless input parameter (float)
        :param params: other dimensionless parameters needed by the CDF
        :return: the value of the CDF evaluated at x (float)
        """
        raise NotImplementedError("Override this in derived classes, or let getCDF() function handle it")

    def getValueRawInvCDF(self, x):
        """
        Override this function to define the analytic form of the inverse of the CDF.

        Is not static so it can access member data and can be broadcasted with numpy.

        This is used to initialize self.cdf
        :param x: dimensionless input parameter (float) in range [0,1] as it is a probability
        :return: the value of the CDF evaluated at x (float)
        """
        raise NotImplementedError("Override this in derived classes, or let getCDF() function handle it")

    @property
    def mean(self):
        return self.integrateWithFunction(lambda x, p: x, tolerance=1e-8)

    @property
    def variance(self):
        m = self.mean
        return self.integrateWithFunction(lambda x, p: x * x, tolerance=1e-8) - m * m

    @property
    def stddev(self):
        return math.sqrt(max(self.variance, 0.0))

    @property
    def skew(self):
        m = self.mean
        s = self.stddev
        if s == 0.0: return 0.0
        return self.integrateWithFunction(lambda x, p: pow((x - m) / s, 3.0), tolerance=1e-8)

    @property
    def kurtosis(self):
        m = self.mean
        s = self.stddev
        if s == 0.0: return 0.0
        return self.integrateWithFunction(lambda x, p: pow(x - m, 4.0), tolerance=1e-8) / pow(s, 4.0)

    def confidenceInterval(self, conf):
        if conf > 1.0 or conf < 0.0: return ValueError('must be in interval [0,1]')
        return self.invcdf.evaluate(1.0 - conf), self.invcdf.evaluate(conf)

    def getCDF(self, params, xGrid, xLabel, xUnit, yLabel="CDF", yUnit="", accuracy=1e-8):
        """
        Test whether getValueRawCDF is defined.  If it is, build the CDF interpolation table from it.

        The initialization is **much** faster if you provide an analytic CDF, integration is slow.  How slow?
        Profiling shows a 10x speed up with analytic CDF's.

        Note: The CDF's y-value is guaranteed to live on the domain [0,1).  However, for sampling, we need it on
        the domain [0,1].  The problem is, the PDF is defined on the interval [0,inf) and the large x values in the
        PDF map to CDF values near 1.  We can make the upper end of the  CDF's domain approach 1 by increasing the
        range of the PDF (which may result in a huge table of values with *very* small probabilities that we'll
        hardly ever sample) or we can cheat.   We're going to cheat below.

        :param params:
        :param xGrid:
        :param xLabel:
        :param xUnit:
        :param yLabel:
        :param yUnit:
        :param accuracy:
        :return:
        """
        have_getValueRawCDF = True
        try:
            self.getValueRawCDF(3.0, {})
        except NotImplementedError:
            have_getValueRawCDF = False
        except Exception:
            have_getValueRawCDF = True

        # OK, so there is a getValueRawCDF, let's use that to build the cdf
        if have_getValueRawCDF:
            cdf = XYs1d.XYs1d.createFromFunction(
                axes=XYs1d.XYs1d.defaultAxes(labelsUnits={1: (xLabel, xUnit), 0: (yLabel, yUnit)}),
                Xs=xGrid,
                func=self.__class__.getValueRawCDF,
                parameters=params,
                accuracy=accuracy,
                biSectionMax=50)

        # If a getValueRawCDF isn't written, we'll have to integrate the pdf
        else:
            cdf = self.indefiniteIntegral()
            cdf.axes[1].label = yLabel
            cdf.axes[1].unit = yUnit

        # The CDF's y-value is guaranteed to live on the domain [0,1).  However, for sampling, we need it on the
        # domain [0,1].  The problem is, the PDF is defined on the interval [0,inf) and the large x values in the
        # PDF map to CDF values near 1.  We can make the upper end of the  CDF's domain approach 1 by increasing the
        # range of the PDF (which may result in a huge table of values with *very* small probabilities that we'll
        # hardly ever sample) or we can cheat.   We're going to cheat.
        if cdf.evaluate(cdf.domainMax) < 1.0:
            cdf.setValue(cdf.domainMax, 1.0)

        return cdf

    def setCDFInverse(self, accuracy=1e-10):
        """
        Sets up the inverse CDF as an XYs1d interpolation table.

        aka Quantile function

        The inverse CDF is the trick used by MCNP to turn any PDF into samples:  you give a random number in the
        interval [0,1], this will give you the x value.  This routine also defines a broadcast ready function
        self.getValueFromInvCDF() which can be used to mass produce samples.


        Note: most of inverse CDF's are "cuspy" for probabilities approaching 1.  These probabilities correspond
        to independent variables in the PDF out on the upper tails.  Therefore, when the XYs1d that holds the
        self.invcdf could return a None (meaning it is off table).  Therefore, we have to be careful about the
        table ends.

        :param accuracy:
        :return: None
        """
        # Return the numpy broadcast capable version for speed!
        try:
            self.getValueRawInvCDF(
                0.5)  # If this doesn't raise an exception, we can make this into a broadcast ready version
            self.getValueFromInvCDF = numpy.frompyfunc(self.getValueRawInvCDF, 1, 1)
        except NotImplementedError:
            self.getValueFromInvCDF = numpy.frompyfunc(self.invcdf.evaluate, 1, 1)

    @property
    def quantileFunction(self):
        return self.cdf.domainSlice(
            domainMin=max(self.cdf.domainMin, self.supportMin),
            domainMax=min(self.cdf.domainMax, self.supportMax)).thin(1e-5).inverse()

    @property
    def cdf(self):
        return self.__cdf

    def drawSample(self, size=None):
        """
        Draw a sample from the pdf.  If size is given, you'll get a numpy.ndarray of size "size" of samples.

        :param size: if not None, the size of the ndarray you'll get as a result
        :return: the sample, either a float or an ndarray of floats
        """
        # Initialize the inverse CDF if it hasn't been done yet
        if not hasattr(self, 'invcdf'):
            self.setCDFInverse()

        # If we are doing many samples at once, use the broadcast ready version
        if type(size) == int:
            return self.getValueFromInvCDF(numpy.random.random(size))

        # Otherwise use the regular version
        return self.invcdf.evaluate(numpy.random.random())

    def isSampleConsistentWithPDF(self, sample, numDecimalPlaces=3):
        """
        Tests whether a list of samples are consistent with the given pdf.  Only checks mean and stddev currently.
        It would be better with skew, kurtosis tests added too

        :param sample: list of floats, these are the statistical samples of something to be tested against the pdf in self.
        :param numDecimalPlaces:
        :return: True if results consistent with pdf in self
        """
        return \
            (round(numpy.mean(sample) - self.mean(), numDecimalPlaces) == 0) and \
            (round(numpy.std(sample) - self.stddev(), numDecimalPlaces) == 0)


class WignerDistribution(UnivariatePDF):
    defaultDomainMin = -0.1
    defaultDomainMax = 5.0

    def __init__(self, xUnit='', xScale=1.0, domainMin=None, domainMax=None, xLabel="indep. variable", yLabel="PDF",
                 **kwds):
        """
        Famous distribution from Wigner's Surmise
        a.k.a. Rayleigh distribution (http://en.wikipedia.org/wiki/Rayleigh_distribution)

        Outside of the interval (domainMin,domainMax), getValue() evaluates to None since this pdf is undefined here

        :param xUnit: units to use for the x axis
        :param xScale: (float or PQU) scale the x values by this number, if domainMin a PQU, will convert to xUnit,
                       otherwise send in a float
        :param domainMin: (float or PQU) pdf and cdf defined on interval (domainMin, domainMax), if domainMin a PQU,
                          will convert to xUnit, otherwise send in a float
        :param domainMax: (float or PQU) pdf and cdf defined on interval (domainMin, domainMax), if domainMax a PQU,
                          will convert to xUnit, otherwise send in a float
        :param xLabel: (str) label for x axis (default "indep. variable")
        :param yLabel: (str) lable for y axis (default "PDF")
        :return: the instance
        """
        # Set up the x axis scale
        if isinstance(xScale, PQU):
            xScale = float(xScale.inUnitsOf(xUnit).value)
        params = {'xScale': xScale}
        self.xScale = xScale

        # Set up domain to define pdf on if not defined
        if domainMin is None:
            domainMin = -xScale / 10.0
        if domainMax is None:
            domainMax = xScale * 6.0

        # Make the pdf
        UnivariatePDF.__init__(self, xUnit=xUnit, domainMin=domainMin, domainMax=domainMax,
                               params=params, xLabel=xLabel, yLabel=yLabel, **kwds)

    @staticmethod
    def getValueRawPDF(x, params):
        """
        Famous distribution from Wigner's Surmise

        :param x: dimensionless input parameter (float)
        :param params:
        :return: the value of the PDF evaluated at x (float)
        """
        if x < 0.0: return 0.0
        xScale = params.get('xScale', 1.0)
        y = x / xScale
        c = math.pi / 2.0
        return c * y * math.exp(-c * y * y / 2.0) / xScale

    @staticmethod
    def getValueRawCDF(x, params):
        """
        CDF of Famous distribution from Wigner's Surmise

        :param x: dimensionless input parameter (float)
        :param params:
        :return: the value of the CDF evaluated at x (float)
        """
        if x < 0.0: return 0.0
        xScale = params.get('xScale', 1.0)
        y = x / xScale
        c = math.pi / 2.0
        return 1.0 - math.exp(-c * y * y / 2.0)

    def drawSample(self, size=None):
        """
        Draw a sample from the pdf.  If size is given, you'll get a numpy.ndarray of size "size" of samples.

        Since we use numpy's random.rayleigh function, samples are draw from

            (x/xScale^2)*exp(-x^2/(2*xScale^2)) #Wikipedia

        so results need to be multiplied a factor of 2 and xScale.  This, however, eliminates troubles
        sampling near y->0 which drives the interpolation in XYs1d nuts.

        :param size: if not None, the size of the ndarray you'll get as a result
        :return: the sample, either a float or an ndarray of floats
        """
        return numpy.random.rayleigh(self.xScale * SQRTTWOOVERPI, size=size)


class GOEDistribution(UnivariatePDF):
    defaultDomainMin = -5.0
    defaultDomainMax = 5.0
    defaultSupportMin = float("-Inf")
    defaultSupportMax = float("Inf")

    def __init__(self, a=1.0, domainMin=None, domainMax=None, xUnit='', xLabel="indep. variable", yLabel="PDF", **kwds):
        r"""
        Simple implementation of the large N limit of the GOE eigenvalue distribution
        (Wigner's semicircle distribution):

            ..math::

                P(x) = \frac{2}{\pi R^2} \sqrt{R^2-x^2}

        on the domain :math:`|x| < R`.  See e.g. https://en.wikipedia.org/wiki/Wigner_semicircle_distribution

        Outside of the interval (domainMin,domainMax), getValue() evaluates to None since this pdf is undefined here

        :param a: (float or PQU) radius of Wigner semi-circle
        :param xUnit: units to use for the x axis
        :param domainMin: (float or PQU) pdf and cdf defined on interval (domainMin, domainMax), if domainMin a PQU,
                          will convert to xUnit, otherwise send in a float
        :param domainMax: (float or PQU) pdf and cdf defined on interval (domainMin, domainMax), if domainMax a PQU,
                          will convert to xUnit, otherwise send in a float
        :param xLabel: (str) label for x axis (default "indep. variable")
        :param yLabel: (str) lable for y axis (default "PDF")
        :return: the instance
        """
        # Get "a" in common units
        if isinstance(a, PQU):
            a = float(a.inUnitsOf(xUnit).value)

        # Set up domain to define pdf on
        domainMin, domainMax = self.fixDomainValues(domainMin, domainMax, xUnit)

        # Set interval of the uniform dist
        renormalize = domainMin > -a or domainMax < a
        params = {'a': a}

        # Make the pdf
        UnivariatePDF.__init__(self, xUnit=xUnit, domainMin=domainMin, domainMax=domainMax,
                               supportMin=-a, supportMax=a, params=params, xLabel=xLabel, yLabel=yLabel, **kwds)
        if renormalize:
            self.normalize(insitu=True)

    @staticmethod
    def getValueRawPDF(x, params):
        """
        Override this function to define the analytic form of the PDF.
        This is used to initialize self.

        :param x: dimensionless input parameter (float)
        :param params: other dimensionless parameters needed by the PDF
        :return: the value of the PDF evaluated at x (float)
        """
        a = params['a']
        if abs(x) > a:
            return 0.0
        a2 = a * a
        return 2.0 * math.sqrt(a2 - x * x) / math.pi / a2

    @staticmethod
    def getValueRawCDF(x, params):
        """
        Override this function to define the analytic form of the CDF.
        This is used to initialize self.cdf.  The actual mathematical
        expression is

            ..math::

                C(x) = \\frac{1}{2} +\\frac{1}{\\pi}\\left[ \\frac{x}{R^2}\\sqrt{R^2-x^2} + \\sin^{-1}(x/r)\\right]

        See e.g. https://en.wikipedia.org/wiki/Wigner_semicircle_distribution

        :param x: dimensionless input parameter (float)
        :param params: other dimensionless parameters needed by the CDF
        :return: the value of the CDF evaluated at x (float)
        """
        a = params['a']
        if x <= -a:
            return 0.0
        if x >= a:
            return 1.0
        a2 = a * a
        x2 = x * x
        y = math.sqrt(a2 - x2)
        return 0.5 + (x * y / a2 + math.asin(x / a)) / math.pi

    @property
    def mean(self):
        return 0.0


class BrodyDistribution(UnivariatePDF):
    def __init__(self, xUnit='', xScale=1.0, w=0.5, domainMin=None, domainMax=None, xLabel="indep. variable",
                 yLabel="PDF", **kwds):
        """
        Brody distribution vaguely interpolates between Poisson and Wigner distributions

        Outside of the interval (domainMin,domainMax), getValue() evaluates to None since this pdf is undefined here

        :param xUnit: units to use for the x axis
        :param w: (float) exponent to use in the equation
        :param xScale: (float or PQU) scale the x values by this number, if domainMin a PQU, will convert to xUnit,
                       otherwise send in a float
        :param domainMin: (float or PQU) pdf and cdf defined on interval (domainMin, domainMax), if domainMin a PQU,
                          will convert to xUnit, otherwise send in a float
        :param domainMax: (float or PQU) pdf and cdf defined on interval (domainMin, domainMax), if domainMax a PQU,
                          will convert to xUnit, otherwise send in a float
        :param xLabel: (str) label for x axis (default "indep. variable")
        :param yLabel: (str) lable for y axis (default "PDF")
        :return: the instance
        """
        # Set up the x axis scale
        if isinstance(xScale, PQU):
            xScale = float(xScale.inUnitsOf(xUnit).value)

        # Save the omega parameter
        self.w = w

        # Load the params dict
        params = {'xScale': xScale, 'w': w}

        # Set up domain to define pdf on if not defined
        if domainMin is None:
            domainMin = -xScale / 10.0
        if domainMax is None:
            domainMax = xScale * 6.0

        # Make the pdf
        UnivariatePDF.__init__(self, xUnit=xUnit, domainMin=domainMin, domainMax=domainMax,
                               params=params, xLabel=xLabel, yLabel=yLabel, **kwds)

    @staticmethod
    def getValueRawPDF(x, params):
        """
        Brody distribution vaguely interpolates between Poisson and Wigner distributions

        :param x: dimensionless input parameter (float)
        :param params: the only other dimensionless parameters needed by the PDF is w,
                       a dimensionless parameter that defines the exponent(s) of x (float),
                       defaults to 0.5, midway between Poisson and Wigner distributions
        :return: the value of the PDF evaluated at x (float)
        """
        if x < 0.0: return 0.0
        xScale = params.get('xScale', 1.0)
        y = x / xScale
        w = params.get('w', 0.5)
        alpha = pow(math.gamma((w + 2.0) / (w + 1.0)), w + 1.0)
        return alpha * (w + 1.0) * pow(y, w) * math.exp(-alpha * pow(y, w + 1.0)) / xScale


class PoissonDistribution(UnivariatePDF):
    def __init__(self, xUnit='', xScale=1.0, domainMin=None, domainMax=None, xLabel="indep. variable", yLabel="PDF",
                 **kwds):
        """
        Poisson distribution (http://en.wikipedia.org/wiki/Poisson_distribution)
        since only k=0 implemented, this is also an exponential distribution
        (http://en.wikipedia.org/wiki/Exponential_distribution)

        Outside of the interval (domainMin,domainMax), getValue() evaluates to None since this pdf is undefined here

        :param xUnit: units to use for the x axis
        :param xScale: (float or PQU) scale the x values by this number, if domainMin a PQU, will convert to xUnit,
                       otherwise send in a float
        :param domainMin: (float or PQU) pdf and cdf defined on interval (domainMin, domainMax), if domainMin a PQU,
                          will convert to xUnit, otherwise send in a float
        :param domainMax: (float or PQU) pdf and cdf defined on interval (domainMin, domainMax), if domainMax a PQU,
                          will convert to xUnit, otherwise send in a float
        :param xLabel: (str) label for x axis (default "indep. variable")
        :param yLabel: (str) lable for y axis (default "PDF")
        :return: the instance
        """
        # Set up the x axis scale
        if isinstance(xScale, PQU):
            xScale = float(xScale.inUnitsOf(xUnit).value)
        params = {'xScale': xScale}

        # Set up domain to define pdf on if not defined
        if domainMin is None:
            domainMin = -xScale / 10.0
        if domainMax is None:
            domainMax = xScale * 6.0

        # Make the pdf
        UnivariatePDF.__init__(self, xUnit=xUnit, domainMin=domainMin, domainMax=domainMax, params=params,
                               xLabel=xLabel, yLabel=yLabel, **kwds)

    @staticmethod
    def getValueRawPDF(x, params):
        """
        Poisson distribution

        :param x: dimensionless input parameter (float)
        :param params:
        :return: the value of the PDF evaluated at x (float)
        """
        if x < 0.0:
            return 0.0
        xScale = params.get('xScale', 1.0)
        y = x / xScale
        return math.exp(-y) / xScale

    @staticmethod
    def getValueRawCDF(x, params):
        """
        Poisson distribution

        :param x: dimensionless input parameter (float)
        :param params:
        :return: the value of the PDF evaluated at x (float)
        """
        if x < 0.0:
            return 0.0
        xScale = params.get('xScale', 1.0)
        y = x / xScale
        return 1.0 - math.exp(-y)


class UniformDistribution(UnivariatePDF):
    def __init__(self, xUnit='', supportMin=0.0, supportMax=1.0, domainMin=None, domainMax=None,
                 xLabel="indep. variable", yLabel="PDF", **kwds):
        """
        Simple implementation of a uniform pdf (http://en.wikipedia.org/wiki/Uniform_distribution_(continuous))

        Outside of the interval (domainMin,domainMax), getValue() evaluates to None since this pdf is undefined here

        :param supportMin: (float or PQU) pdf is constant on interval (supportMin, supportMax), if supportMin a PQU,
                           will convert to xUnit, otherwise send in a float
        :param supportMax: (float or PQU) pdf is constant on interval (supportMin, supportMax), if supportMax a PQU,
                           will convert to xUnit, otherwise send in a float
        :param xUnit: units to use for the x axis
        :param domainMin: (float or PQU) pdf and cdf defined on interval (domainMin, domainMax), if domainMin a PQU,
                          will convert to xUnit, otherwise send in a float
        :param domainMax: (float or PQU) pdf and cdf defined on interval (domainMin, domainMax), if domainMax a PQU,
                          will convert to xUnit, otherwise send in a float
        :param xLabel: (str) label for x axis (default "indep. variable")
        :param yLabel: (str) lable for y axis (default "PDF")
        :return: the instance
        """
        # Set up domain to define pdf on
        domainMin, domainMax = self.fixDomainValues(domainMin, domainMax, xUnit)

        # Set interval of the uniform dist
        if isinstance(supportMin, PQU):
            xMin = float(supportMin.inUnitsOf(xUnit).value)
        if isinstance(supportMax, PQU):
            xMax = float(supportMax.inUnitsOf(xUnit).value)
        if domainMax < supportMax:
            raise ValueError("domainMax=%g < supportMax=%g" % (domainMax, supportMax))
        if supportMin < domainMin:
            raise ValueError("supportMin=%g <= domainMin=%g" % (supportMin, domainMin))
        if supportMax <= supportMin:
            raise ValueError("supportMin=%g >= supportMax=%g" % (supportMin, supportMax))
        params = {}
        params['supportMin'] = supportMin
        params['supportMax'] = supportMax

        # Make the pdf
        UnivariatePDF.__init__(self, xUnit=xUnit, domainMin=domainMin, domainMax=domainMax, supportMin=supportMin,
                               supportMax=supportMax, params=params, xLabel=xLabel, yLabel=yLabel, **kwds)

    @staticmethod
    def getValueRawPDF(x, params):
        """
        Override this function to define the analytic form of the PDF.
        This is used to initialize self.

        :param x: dimensionless input parameter (float)
        :param params:
        :return: the value of the PDF evaluated at x (float)
        """
        supportMin = params.get('supportMin', 0.0)
        supportMax = params.get('supportMax', 1.0)
        p = 1.0 / (supportMax - supportMin)
        if supportMin <= x <= supportMax:
            return p
        return 0.0

    @staticmethod
    def getValueRawCDF(x, params):
        """
        Override this function to define the analytic form of the CDF.
        This is used to initialize self.cdf

        :param x: dimensionless input parameter (float)
        :param params:
        :return: the value of the CDF evaluated at x (float)
        """
        supportMin = params.get('supportMin', 0.0)
        supportMax = params.get('supportMax', 1.0)
        if x < supportMin:
            return 0.0
        if x >= supportMax:
            return 1.0
        return (x - supportMin) / (supportMax - supportMin)


class NormalDistribution(UnivariatePDF):
    """(http://en.wikipedia.org/wiki/Normal_distribution)"""
    defaultSupportMin = -float("Inf")
    defaultSupportMax = float("Inf")

    def __init__(self, xUnit='', mean=0.0, stddev=0.0, domainMin=None, domainMax=None, xLabel="indep. variable",
                 yLabel="PDF", **kwds):
        """
        Normal a.k.a. Gaussian distribution (http://en.wikipedia.org/wiki/Normal_distribution)

        Outside of the interval (domainMin,domainMax), evaluate() evaluates to None since this pdf is undefined here

        :param xUnit: units to use for the x axis
        :param mean: (float or PQU) XXXXXXX, if domainMin a PQU, will convert to xUnit, otherwise send in a float
        :param variance: (float or PQU) XXXXXXX, if domainMin a PQU, will convert to xUnit, otherwise send in a float
        :param domainMin: (float or PQU) pdf and cdf defined on interval (domainMin, domainMax), if domainMin a PQU,
                          will convert to xUnit, otherwise send in a float
        :param domainMax: (float or PQU) pdf and cdf defined on interval (domainMin, domainMax), if domainMax a PQU,
                          will convert to xUnit, otherwise send in a float
        :param xLabel: (str) label for x axis (default "indep. variable")
        :param yLabel: (str) lable for y axis (default "PDF")
        :return: the instance
        """
        # Set up the x axis scale
        if isinstance(mean, PQU):
            mean = float(mean.inUnitsOf(xUnit).value)
        if isinstance(stddev, PQU):
            stddev = float(stddev.inUnitsOf(xUnit).value)
        params = {'mean': mean, 'stddev': stddev}
        self.__mean = mean
        self.__stddev = stddev

        # Set up domain to define pdf on if not defined
        if domainMin is None:
            domainMin = -6.0 * stddev + mean
        if domainMax is None:
            domainMax = 6.0 * stddev + mean

        # Make the pdf
        UnivariatePDF.__init__(self, xUnit=xUnit, domainMin=domainMin, domainMax=domainMax,
                               params=params, xLabel=xLabel, yLabel=yLabel, **kwds)

    @property
    def mean(self):
        return self.__mean

    @property
    def variance(self):
        return self.stddev * self.stddev

    @property
    def stddev(self):
        return self.__stddev

    @staticmethod
    def getValueRawPDF(x, params):
        """
        :param x: dimensionless input parameter (float)
        :param params: other dimensionless parameters needed by the PDF
        :return: the value of the PDF evaluated at x (float)
        """
        dx = params['stddev']
        x0 = params['mean']
        xx = (x - x0) / dx
        return numpy.exp(-xx * xx / 2.) / SQRTTWOPI / dx

    @staticmethod
    def getValueRawCDF(x, params):
        """
        :param x: dimensionless input parameter (float)
        :param params:
        :return: the value of the CDF evaluated at x (float)
        """
        import scipy.special as sp
        dx = params['stddev']
        x0 = params['mean']
        return 0.5 * (1. + sp.erf((x - x0) / SQRT2 / dx))

    def drawSample(self, size=None):
        """
        Draw a sample from the pdf.  If size is given, you'll get a numpy.ndarray of size "size" of samples.

        :param size: if not None, the size of the ndarray you'll get as a result
        :return: the sample, either a float or an ndarray of floats
        """
        return numpy.random.normal(loc=self.mean, scale=self.stddev, size=size)


class LogNormalDistribution(UnivariatePDF):
    """(http://en.wikipedia.org/wiki/Log-normal_distribution)"""
    defaultSupportMin = 0.0
    defaultSupportMax = float("Inf")

    def __init__(self, xUnit='', sigma=0.0, mu=0.0, domainMin=None, domainMax=None, xLabel="indep. variable",
                 yLabel="PDF", **kwds):
        """
        Log normal distribution (http://en.wikipedia.org/wiki/Log-normal_distribution)

        Outside of the interval (domainMin,domainMax), evaluate() evaluates to None since this pdf is undefined here

        :param xUnit: units to use for the x axis
        :param mean: (float or PQU) XXXXXXX, if domainMin a PQU, will convert to xUnit, otherwise send in a float
        :param variance: (float or PQU) XXXXXXX, if domainMin a PQU, will convert to xUnit, otherwise send in a float
        :param domainMin: (float or PQU) pdf and cdf defined on interval (domainMin, domainMax),
                          if domainMin a PQU, will convert to xUnit, otherwise send in a float
        :param domainMax: (float or PQU) pdf and cdf defined on interval (domainMin, domainMax),
                          if domainMax a PQU, will convert to xUnit, otherwise send in a float
        :param xLabel: (str) label for x axis (default "indep. variable")
        :param yLabel: (str) lable for y axis (default "PDF")
        :return: the instance
        """
        # Set up the x axis scale
        if isinstance(sigma, PQU):
            mean = float(sigma.inUnitsOf(xUnit).value)
        if isinstance(mu, PQU):
            stddev = float(mu.inUnitsOf(xUnit).value)
        params = {'sigma': sigma, 'mu': mu}
        self.sigma = sigma
        self.mu = mu

        # Set up domain to define pdf on if not defined
        if domainMin is None:
            domainMin = 0.0
        if domainMax is None:
            domainMax = 6.0 * self.stddev + self.mean

        # Make the pdf
        UnivariatePDF.__init__(self, xUnit=xUnit, domainMin=domainMin, domainMax=domainMax,
                               params=params, xLabel=xLabel, yLabel=yLabel, **kwds)

    @property
    def mean(self):
        return math.exp(self.mu + self.sigma * self.sigma / 2.0)

    @property
    def variance(self):
        sig2 = self.sigma * self.sigma
        return (math.exp(sig2) - 1.0) * math.exp(2.0 * self.mu + sig2)

    @staticmethod
    def getValueRawPDF(x, params):
        """
        Override this function to define the analytic form of the PDF.
        This is used to initialize self.

        :param x: dimensionless input parameter (float)
        :param params: other dimensionless parameters needed by the PDF
        :return: the value of the PDF evaluated at x (float)
        """
        sig = params['sigma']
        mu = params['mu']
        if x < 1e-8: return 0.0
        a = (numpy.log(x) - mu) / sig
        return numpy.exp(-a * a / 2.0) / x / SQRTTWOPI / sig

    @staticmethod
    def getValueRawCDF(x, params):
        """
        Override this function to define the analytic form of the CDF.
        This is used to initialize self.cdf

        :param x: dimensionless input parameter (float)
        :param params:
        :return: the value of the CDF evaluated at x (float)
        """
        import scipy.special as sp
        sig = params['sigma']
        mu = params['mu']
        if x < 1e-8: return 0.0
        a = (numpy.log(x) - mu) / sig
        return 0.5 * (1.0 + sp.erf(a / 2.0))

    def drawSample(self, size=None):
        """
        Draw a sample from the pdf.  If size is given, you'll get a numpy.ndarray of size "size" of samples.

        :param size: if not None, the size of the ndarray you'll get as a result
        :return: the sample, either a float or an ndarray of floats
        """
        return numpy.random.lognormal(mean=self.mu, sigma=self.sigma, size=size)


class Chi2Distribution(UnivariatePDF):
    """
    Chi^2 distribution (http://en.wikipedia.org/wiki/Chi-squared_distribution)

    Nuclear physicists know it as the Porter-Thomas distribution
    """

    def __init__(self, dof, xUnit='', xScale=1.0, domainMin=None, domainMax=None, xLabel="indep. variable",
                 yLabel="PDF", **kwds):
        """


        Outside of the interval (domainMin,domainMax), getValue() evaluates to None since this pdf is undefined here

        :param xUnit: units to use for the x axis
        :param xScale: (float or PQU) scale the x values by this number, if domainMin a PQU,
                       will convert to xUnit, otherwise send in a float
        :param domainMin: (float or PQU) pdf and cdf defined on interval (domainMin, domainMax),
                          if domainMin a PQU, will convert to xUnit, otherwise send in a float
        :param domainMax: (float or PQU) pdf and cdf defined on interval (domainMin, domainMax),
                          if domainMax a PQU, will convert to xUnit, otherwise send in a float
        :param xLabel: (str) label for x axis (default "indep. variable")
        :param yLabel: (str) label for y axis (default "PDF")
        :return: the instance
        """
        # Set up the x axis scale
        if isinstance(xScale, PQU):
            xScale = float(xScale.inUnitsOf(xUnit).value)
        self.xScale = xScale

        # Save the number of degrees of freedom, dof
        if dof == 0.0:
            raise ValueError("DOF->0 makes this distribution a delta function, you'd better handle it yourself")
        elif dof < 0.0:
            raise ValueError("DOF<0 is illegal")
        elif dof < 1.0:
            raise RuntimeWarning("DOF < 1.0 is very unlikely for width distribution")
        self.dof = dof

        # Load the params
        params = {'xScale': xScale, 'dof': dof}

        # Set up domain to define pdf on if not defined
        if domainMin is None:
            domainMin = -xScale / 10.0
        if domainMax is None:
            domainMax = xScale * 15.0  # 45.0, use larger domainMax when using Wikipedia's version of chi^2 dist

        # Make the pdf
        UnivariatePDF.__init__(self, xUnit=xUnit, domainMin=domainMin, domainMax=domainMax,
                               params=params, xLabel=xLabel, yLabel=yLabel, **kwds)

    @staticmethod
    def getValueRawPDF(x, params):
        """
        Override this function to define the analytic form of the PDF.
        This is used to initialize self.

        :param x: dimensionless input parameter (float)
        :param params: other dimensionless parameters needed by the PDF
        :return: the value of the PDF evaluated at x (float)
        """
        if x < 0.0: return 0.0
        xScale = params.get('xScale', 1.0)
        y = x / xScale
        nu = params['dof']
        gam = math.gamma(nu / 2.0)
        if nu < 2.0:
            ycut = pow(1e4 * gam, -1.0 / (1.0 - nu / 2.0))
            if y < ycut:
                y = ycut
        return math.exp(-y) * pow(y, nu / 2.0 - 1.0) / gam  # Froehner's eq. (277)

    @property
    def mean(self):
        return self.xScale * self.dof / 2.0

    @property
    def variance(self):
        return self.dof * self.xScale * self.xScale / 2.0

    def drawSample(self, size=None):
        """
        Draw a sample from the pdf.  If size is given, you'll get a numpy.ndarray of size "size" of samples.

        Since we use numpy's random.chisquare function, samples are draw from

            return math.exp(-y/2.0)*pow(y,nu/2.0-1.0)/math.gamma(nu/2.0)/pow(2.0,nu/2.0) #Wikipedia

        so results need to be multiplied a factor of 2 and xScale.  This, however, eliminates troubles
        sampling near y->0 which drives the interpolation in XYs1d nuts.

        :param size: if not None, the size of the ndarray you'll get as a result
        :return: the sample, either a float or an ndarray of floats
        """
        return self.xScale * numpy.random.chisquare(self.dof, size=size) / 2.0
