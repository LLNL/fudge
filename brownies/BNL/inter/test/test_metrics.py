import unittest
from brownies.BNL.inter.metrics import *
from brownies.BNL.utilities.io import *


class Test_With_ApproxDiff(unittest.TestCase):

    def assertIsClose(self, a, b, abs_tol=1e-9, rel_tol=1e-8, method='average'):
        if isinstance(a, PQU.PQU) or isinstance(b, PQU.PQU):
            if not (isinstance(a, PQU.PQU) and isinstance(b, PQU.PQU)):
                raise TypeError("Cannot mix types, %s with %s" % (str(type(a)), str(type(b))))
            a = float(a.inUnitsOf(b.unit))
            b = float(b)
        delta = abs(a - b)
        absave = abs(a + b) / 2.0
        if absave != 0.0:
            reldelta = delta / absave
        else:
            reldelta = 0.0
        self.assertTrue(math.isclose(a, b, abs_tol=abs_tol, rel_tol=rel_tol),
                        "%s is not close to %s, absdiff is %s, reldiff is %s" % (
                        str(a), str(b), str(delta), str(reldelta)))


class Test_INTER_Quantities(Test_With_ApproxDiff):
    """
    Results of n-001_H_001.endf ran through legacy INTER ::

                                     PROGRAM INTER VERSION 8.08
                                     --------------------------

         Selected Integrations of ENDF File 3 and File 10 Cross Sections

         Thermal cross section : Sig(2200) = Sig(Eth)
         Thermal energy        : Eth=  2.53000E-02 (eV)

         Ezero cross section   : Sig(Ezero)
         Ezero energy (input)  : E0 =  2.53000E-02 (eV)

         Maxwellian average    : Avg-Sigma = (2/sqrt(Pi)) Intg[E1:E2] Sig(E) Phi_m(E) dE / Intg[E1:E2] Phi_m(E) dE
         Maxwellian spectrum   : Phi_m(E)  = (E/E0^2) exp(-E/E0)
         Spectrum Temperature  : E0 =  2.53000E-02 (eV)
         Integration Limits    : E1 =  1.00000E-05 (eV)  E2 =  1.00000E+01 (eV)
         Integral of Spectrum       =  1.00000E+00

         Westcott g-factor     : G-fact = 2/sqrt(Pi)  Avg-Sigma / Sig(2200)

         Resonance Integral    : Res.Integ = Intg[E3:E4] Sig(E)/E dE
         Integration Limits    : E3 =  5.00000E-01 (eV)  E4 =  1.00000E+05 (eV)
         Integral of Spectrum       =  1.22061E+01

         Fiss.spect. average   : Sig(Fiss) = Intg[E1:E2] Sig(E) Phi_f(E) dE / Intg[E1:E2] Phi_f(E) dE
         Fission spectrum      : Phi_f(E)  = sqrt(E/Ef)/Ef exp(-E/E0)
         Spectrum Temperature  : Ef =  1.35000E+06 (eV)
         Integration Limits    : E1 =  1.00000E+03 (eV)  E2 =  2.00000E+07 (eV)
         Integral of Spectrum       =  1.00000E+00

         E14 cross-section     : Sig(E14)
         Selected Energy       : E14 =  1.40000E+07 eV


           Z   A LISO  LFS  MT  Reaction    Sig(2200)   Sig(Ezero)  Avg-Sigma  G-fact   Res Integ   Sig(Fiss)    Sig(E14)   MAT
         --- ------------- ---  ---------- ----------- ----------- ---------- -------  ----------- ----------- ----------- ----
           1   1             1  Total      2.07683E+01 2.07683E+01 2.3402E+01 1.12680  2.39596E+02 3.98828E+00 6.87144E-01  125
           1   1             2  Elastic    2.04363E+01 2.04363E+01 2.3060E+01 1.12838  2.39452E+02 3.98824E+00 6.87114E-01  125
           1   1           102  n,gamma    3.32013E-01 3.32013E-01 3.3214E-01 1.00040  1.48914E-01 3.95532E-05 2.98051E-05  125
    """

    def setUp(self):
        self.eval = read_evaluation(INTERDIR + '/test/n-001_H_001.endf')

    def test_read_evaluation(self):

        listAssertEqual = getattr(self, 'assertCountEqual', None)
        if listAssertEqual is None:
            listAssertEqual = getattr(self, 'assertItemsEqual')

        listAssertEqual(list(self.eval.keys()), ['info', 'reactionSuite', 'errors', 'covarianceSuite'])

    def test_check_is_cross_section(self):
        self.assertIsNone(check_is_cross_section(self.eval['reactionSuite'].getReaction('capture').crossSection))

    def test_cross_section_evaluateWithUncertainty_function(self):
        x = self.eval['reactionSuite'].getReaction('capture').crossSection.evaluateWithUncertainty(
            PQU.PQU('14.2 MeV'))
        self.assertAlmostEqual(float(x.value), 2.9719802000000002e-05, 6)
        self.assertAlmostEqual(float(x.uncertainty), 2.2832109016510004e-06, 6)
        self.assertEqual(str(x.unit), 'b')

    def test_cross_section_integrateTwoFunctionsWithUncertainty_function(self):
        answers = {1: 263.88094687, 2: 263.726362691, 102: 0.149383616104}
        for r in self.eval['reactionSuite']:
            if r.ENDF_MT in [1, 2, 102]:
                a = computeRI(r.crossSection, useCovariance=False)
                self.assertIsClose(float(a.value), answers[r.ENDF_MT], rel_tol=0.005)

    def test_RI_vs_INTER(self):
        for r in self.eval['reactionSuite']:
            if r.ENDF_MT == 1:
                self.assertIsClose(computeRI(r.crossSection, domainMax=PQU.PQU(1.00000E+05, 'eV'), useCovariance=False),
                                   PQU.PQU(2.39596E+02, 'b'), rel_tol=1e-5)
            if r.ENDF_MT == 2:
                self.assertIsClose(computeRI(r.crossSection, domainMax=PQU.PQU(1.00000E+05, 'eV'), useCovariance=False),
                                   PQU.PQU(2.39452E+02, 'b'), rel_tol=1e-6)
            if r.ENDF_MT == 102:
                self.assertIsClose(computeRI(r.crossSection, domainMax=PQU.PQU(1.00000E+05, 'eV'), useCovariance=False),
                                   PQU.PQU(1.48914E-01, 'b'), abs_tol=5e-6)

    def test_14MeV(self):
        for r in self.eval['reactionSuite']:
            # Just test the mean values first
            if r.ENDF_MT == 1:
                self.assertIsClose(computeFourteenMeVPoint(r.crossSection, "1.40000E+07 eV", useCovariance=False),
                                   PQU.PQU(6.87144E-01, 'b'), rel_tol=5e-7)
            if r.ENDF_MT == 2:
                self.assertIsClose(computeFourteenMeVPoint(r.crossSection, "1.40000E+07 eV", useCovariance=False),
                                   PQU.PQU(6.87114E-01, 'b'), rel_tol=5e-7)
            if r.ENDF_MT == 102:
                self.assertIsClose(computeFourteenMeVPoint(r.crossSection, "1.40000E+07 eV", useCovariance=False),
                                   PQU.PQU(2.98051E-05, 'b'), rel_tol=5e-7)
            # Now turn on covariance, these tests should still all pass OK
            if r.ENDF_MT == 1:
                self.assertIsClose(computeFourteenMeVPoint(r.crossSection, "1.40000E+07 eV", useCovariance=True),
                                   PQU.PQU(6.87144E-01, 'b'), rel_tol=5e-7)
            if r.ENDF_MT == 2:
                self.assertIsClose(computeFourteenMeVPoint(r.crossSection, "1.40000E+07 eV", useCovariance=True),
                                   PQU.PQU(6.87114E-01, 'b'), rel_tol=5e-7)
            if r.ENDF_MT == 102:
                self.assertIsClose(computeFourteenMeVPoint(r.crossSection, "1.40000E+07 eV", useCovariance=True),
                                   PQU.PQU(2.98051E-05, 'b'), rel_tol=5e-7)

    def test_RoomTemp(self):
        for r in self.eval['reactionSuite']:
            if r.ENDF_MT == 1:
                self.assertIsClose(computeRoomTempCS(r.crossSection), PQU.PQU(2.07687319E+01, 'b'), rel_tol=5e-7)
            if r.ENDF_MT == 2:
                self.assertIsClose(computeRoomTempCS(r.crossSection), PQU.PQU(2.043633E+01, 'b'), rel_tol=5e-7)
            if r.ENDF_MT == 102:
                self.assertIsClose(computeRoomTempCS(r.crossSection), PQU.PQU(3.322811E-01, 'b'), rel_tol=5e-7)

    def test_Cf252Analytic(self):
        for r in self.eval['reactionSuite']:
            if r.ENDF_MT == 1:
                self.assertIsClose(computeCf252SpectrumAveAnalytic(r.crossSection, useCovariance=False),
                                   PQU.PQU(3.98828E+00, 'b'), rel_tol=0.02)
            if r.ENDF_MT == 2:
                self.assertIsClose(computeCf252SpectrumAveAnalytic(r.crossSection, useCovariance=False),
                                   PQU.PQU(3.98824E+00, 'b'), rel_tol=0.02)
            if r.ENDF_MT == 102:
                self.assertIsClose(computeCf252SpectrumAveAnalytic(r.crossSection, useCovariance=False),
                                   PQU.PQU(3.95532E-05, 'b'), rel_tol=0.01)

    def test_Westcott(self):
        """
        Test against values from INTER (recalculated to a higher precision).
        Note: INTER ignores recoil of the target, so uses a=1.0 rather than a=mTarget/(mTarget+mProjectile)
        """
        for r in self.eval['reactionSuite']:
            if r.ENDF_MT == 1:
                self.assertIsClose(computeWestcottFactor(r.crossSection, a=1.0, useCovariance=False),
                                   PQU.PQU(1.12677222244, ''))
            if r.ENDF_MT == 2:
                self.assertIsClose(computeWestcottFactor(r.crossSection, a=1.0, useCovariance=False),
                                   PQU.PQU(1.12837902942, ''))
            if r.ENDF_MT == 102:
                self.assertIsClose(computeWestcottFactor(r.crossSection, a=1.0, useCovariance=False),
                                   PQU.PQU(1.00034149453, ''))

    def test_MACS(self):
        """
        These test values come from B. Pritychenko, S.F. Mughabghab arXiv:1208.2879v3.
        """
        for r in self.eval['reactionSuite']:
            if r.ENDF_MT == 102:
                self.assertIsClose(computeMACS(r.crossSection, PQU.PQU(30., 'keV'), useCovariance=False),
                                   PQU.PQU(1.525E-4, 'b'), abs_tol=1e-7)
                self.assertIsClose(computeMACS(r.crossSection, PQU.PQU(1420., 'keV'), useCovariance=False),
                                   PQU.PQU(3.892E-5, 'b'), abs_tol=2e-8)

                # Turn on covariance, what happens then?
                self.assertIsClose(computeMACS(r.crossSection, PQU.PQU(30., 'keV'), useCovariance=True),
                                   PQU.PQU("1.524e-4 +/- 5.6e-6 b"), abs_tol=1e-7)
                break

    def test_scatteringRadius(self):
        """
        These test values come from B. Pritychenko, S.F. Mughabghab arXiv:1208.2879v3.
        """
        self.assertIsClose(computeScatteringRadius(self.eval['reactionSuite']), PQU.PQU(12.76553, 'fm'),
                           abs_tol=1e-7)

    def test_ARR(self):
        """
        These test values come from B. Pritychenko, S.F. Mughabghab arXiv:1208.2879v3,
        """
        for r in self.eval['reactionSuite']:
            if r.ENDF_MT == 102:
                self.assertIsClose(computeAstrophysicalReactionRate(r.crossSection, PQU.PQU(30., 'keV'),
                                                                    useCovariance=False),
                                                                    PQU.PQU(3.113E4, "cm**3/s/mol"), rel_tol=0.001)
                break

    def test_computeCf252SpectrumAve(self):
        """
        These results were compared to the legacy INTER fission spectrum average and we thought our results were
        reasonable so we took our calculated results as the "answer" in the comparison of this test.
        """
        for r in self.eval['reactionSuite']:
            if r.ENDF_MT == 1:
                self.assertIsClose(computeCf252SpectrumAve(r.crossSection, useCovariance=False),
                                   PQU.PQU(3.98828E+00, 'b'), rel_tol=0.09)
            if r.ENDF_MT == 2:
                self.assertIsClose(computeCf252SpectrumAve(r.crossSection, useCovariance=False),
                                   PQU.PQU(3.98824E+00, 'b'), rel_tol=0.09)
            if r.ENDF_MT == 102:
                self.assertIsClose(computeCf252SpectrumAve(r.crossSection, useCovariance=False),
                                   PQU.PQU(3.95532E-05, 'b'), rel_tol=0.09)

    def test_computeGodivaSpectrumAve(self):
        """
        These results are taken from our results.  We think their OK; they are certainly in the ball pack
        (compared to the Cf252 spectrum average above)
        """
        for r in self.eval['reactionSuite']:
            if r.ENDF_MT == 102:
                self.assertIsClose(computeGodivaSpectrumAve(r.crossSection, useCovariance=False),
                                   PQU.PQU(3.60870104207e-05, 'b'))

    def test_computeJezebelSpectrumAve(self):
        """
        These results are taken from our results.  We think their OK; they are certainly in the ball pack
        (compared to the Cf252 spectrum average above)
        """
        for r in self.eval['reactionSuite']:
            if r.ENDF_MT == 102:
                self.assertIsClose(computeJezebelSpectrumAve(r.crossSection, useCovariance=False),
                                   PQU.PQU(3.58646505181e-05, 'b'))

    def test_computeBigTenSpectrumAve(self):
        """
        These results are taken from our results.  We think their OK; they are certainly in the ball pack
        (compared to the Cf252 spectrum average above)
        """
        for r in self.eval['reactionSuite']:
            if r.ENDF_MT == 102:
                self.assertIsClose(computeBigTenSpectrumAve(r.crossSection, useCovariance=False),
                                   PQU.PQU(3.98659484215e-05, 'b'))


# -------------------------------------------------------------------------------
# Unit tests
# -------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
