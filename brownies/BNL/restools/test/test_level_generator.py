import unittest
from xData import axes
from xData import XYs1d
from brownies.BNL.restools.level_generator import *

DEBUG = True
PROFILE = False


class AssertMixIn(unittest.TestCase):

    def assertWithinNSigma(self, A, B, Sigma, N=2, additional_msg=None):
        msg = "|%f-%f| = %f < %i*%f = %f not satisfied" % (A, B, abs(A - B), N, Sigma, N * Sigma)
        if additional_msg is not None:
            msg += additional_msg
        if Sigma < 1e-7:
            self.assertAlmostEqual(A, B, 6, msg=msg)
        else:
            self.assertTrue(abs(A - B) < N * Sigma, msg=msg)


# ---------------------------------------------------------------------------------
#
#   Unit tests
#
# ---------------------------------------------------------------------------------

class TestLevelGeneration(AssertMixIn):

    def setUp(self):
        self.T = 8.0e6
        self.Emin = 0.0
        self.Emax = 1.0e6
        self.nE = 1000
        self.dE = self.Emax - self.Emin

        # Make a level density rho(E) = e^(E/T)
        self.ld = XYs1d.XYs1d.createFromFunction(
            axes=axes.Axes(2, labelsUnits={0: ('rho(E)', '1/eV'), 1: ("E", "eV")}),
            Xs=[x * self.dE / self.nE for x in range(self.nE)],
            func=lambda x, par: math.exp(x / self.T),
            parameters=None,
            accuracy=1e-8,
            biSectionMax=-1)

    def test_semicircle(self):
        """
        Cummulative level distribution (CLD) of the GOE distribution should be a Wigner semicircle.
        We'll test this assertion by computing the "radius" computed from the x & y(x) values of the CLD.

        Since we've Monte-Carlo'd the sample of GOE levels, we have to approach this with some measure
        of uncertainty ;)

        Wigner's semicircle law is the following probability density:

            ..math::

                f(x)=\left{
                    \begin{array}
                        \frac{2}{\pi R^2}\sqrt{R^2-x^2} & -R\le x \le R \\
                        0                               & |x| > R
                    \end{array}\right.
        """
        from brownies.BNL.utilities.XYs import makeHistoXYs
        nSamples = 1000
        sample = sample_goe_eigenvalues(nSamples)  # picked R = 1.0
        cld = makeHistoXYs(sample, 50, normfactor=1.0 / math.sqrt(
            math.pi * nSamples))  # renormalized circle so peak value at x=0 is 2
        R = [math.sqrt(xy[0] * xy[0] + xy[1] * xy[1]) for xy in cld]
        if False:
            cld.plot()
            import fudge.core.math.pdf as pdf
            pdf.GOEDistribution(a=1.0).plot()

        # Test that "radius" of each entry in the CLD is as expected, namely 1/sqrt(pi)
        self.assertWithinNSigma(
            numpy.average(R),
            1.0 / math.sqrt(math.pi),
            numpy.std(R),
            additional_msg="; for radius of Wigner semicircle")

    def test_level_density(self):
        self.assertAlmostEqual(self.ld.evaluate(5e5), math.exp(5e5 / self.T))

        # Compute number of levels
        Ntot = self.ld.integrate()
        self.assertAlmostEqual(Ntot / 1e6, self.T * (math.exp(self.Emax / self.T) - 1.0) / 1e6, 2)

    def test_getCLDInverse(self):
        # Get the inverse cumulative level distribution function.
        # For y=CDF(E) should be y=T*ln(y*Ntot/T+1) since CDF(E)=(T/Ntot)*(exp(E/T)-1)
        invcld = getCLDInverse(self.ld)
        self.assertAlmostEqual(invcld.evaluate(0.5) / 1e5,
                               self.T * math.log(0.5 * self.ld.integrate() / self.T + 1.0) / 1e5, 5)

    def test_getGOEFakeLevelSequence(self):
        E0 = 4.0
        totalNumLevels = 400
        levelDensity = self.ld
        seq0 = getGOEFakeLevelSequence(
            E0,
            totalNumLevels,
            levelDensity,
            paddingNumLevels=100,
            keepLevelsAboveDomainMax=False,
            DOPLOTS=False)
        self.assertEqual(seq0[0], E0)
        self.assertLessEqual(len(seq0), totalNumLevels)
        self.assertGreater(len(seq0), totalNumLevels - 5)
        seq1 = getGOEFakeLevelSequence(
            E0,
            totalNumLevels,
            levelDensity,
            paddingNumLevels=100,
            keepLevelsAboveDomainMax=False,
            DOPLOTS=False)
        self.assertEqual(seq1[0], E0)
        self.assertLessEqual(len(seq1), totalNumLevels)
        self.assertGreater(len(seq1), totalNumLevels - 5)
        self.assertNotEqual(seq0, seq1)


# ---------------------------------------------------------------------------------
#
#   Main, used for running unit tests.
#   If you want to run something real, run __main__
#
# ---------------------------------------------------------------------------------

if __name__ == "__main__":

    # Actually run the tests
    if PROFILE:
        import cProfile

        cProfile.run("unittest.main()", "burr__init__.profile")
    else:
        unittest.main()
        print()
