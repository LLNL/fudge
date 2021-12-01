from pqu.PQU import PhysicalQuantityWithUncertainty as PQU
from fudge.core.math.pdf import WignerDistribution, Chi2Distribution
from brownies.BNL.restools import level_generator
from brownies.BNL.restools.level_analysis import *
from brownies.BNL.utilities.XYs import makeHistoXYs, makeXYs
import unittest
import math
import numpy

DEBUG = False
TESTALL = True
PROFILE = False


def Wigner_surmise_distribution(x):
    constant = (math.pi / 2)
    exponent = ((-1) * (math.pi / 4) * x * x)
    return (constant * x) * math.exp(exponent)


class AssertMixIn(unittest.TestCase):

    def assertWithinNSigma(self, A, B, Sigma, N=2, additional_msg=None):
        msg = "|%f-%f| = %f < %i*%f = %f not statisfied" % (A, B, abs(A - B), N, Sigma, N * Sigma)
        if additional_msg is not None: msg += additional_msg
        if Sigma < 1e-7:
            self.assertAlmostEqual(A, B, 6, msg=msg)
        else:
            self.assertTrue(abs(A - B) < N * Sigma, msg=msg)


# ---------------------------------------------------------------------------------
#
#   Unit tests
#
# ---------------------------------------------------------------------------------

class TestLevelAnalyzer_simple(AssertMixIn):

    def setUp(self):
        self.E0 = PQU('1 eV')
        self.D = PQU('10 eV')
        self.units = 'eV'
        self.nsamples_big = 10000
        self.nsamples_medium = 1000
        self.nsamples_small = 100
        self.nsamples_vsmall = 10
        self.debug = False
        self.expectedAverageD = 3581.25583223
        self.expectedNLevels = 442
        self.expectedLinearCorrelationCoefficient = -0.27
        self.expectedAverageU = 0.365
        self.expectedAverageQ = 0.0  # .340
        self.picketFence = level_generator.getFakeLevelSequence(
            E0=10.0, style='picket fence', aveD=self.expectedAverageD, numLevels=self.expectedNLevels)
        self.picketFenceSequenceAnalyzer = LevelSequenceAnalyzer(levels=self.picketFence)

    def test_getAverageSpacingLevelDensity(self):
        E = self.picketFenceSequenceAnalyzer.levels[4]
        self.assertWithinNSigma(
            self.picketFenceSequenceAnalyzer.getAverageSpacingLevelDensity().evaluate(E),
            1.0 / self.picketFenceSequenceAnalyzer.mean_spacing, 1e-6)

    def test_level_unfolding(self):
        for i, lev in enumerate(self.picketFenceSequenceAnalyzer.unfolded_levels):
            self.assertAlmostEqual(i, lev)

    def test_statics_for_wigner_distribution(self):
        """
        Levels should follow Wigner distribution.
        Mean is D, the average level spacing (set to 10 eV).
        Variance is D^2*(4-pi)/pow(2,3/2)
        """
        la = LevelSequenceAnalyzer(
            levels=[float(x.inUnitsOf(self.units).value) for x in
                    level_generator.getFakeLevelSequence(self.E0, self.D, self.nsamples_big, style='wigner')])
        self.assertAlmostEqual(la.mean_spacing, float(self.D.inUnitsOf(self.units).value), 0)
        self.assertAlmostEqual(la.ave_spacing, float(self.D.inUnitsOf(self.units).value), 0)
        self.assertAlmostEqual(
            float(la.var_spacing / self.D.inUnitsOf(self.units).value / self.D.inUnitsOf(self.units).value),
            (4.0 - math.pi) / 2.0 / math.sqrt(2.0), 1)
        self.assertAlmostEqual(float(la.stddev_spacing / self.D.inUnitsOf(self.units).value),
                               math.sqrt((4.0 - math.pi) / 2.0 / math.sqrt(2.0)), 1)

    def test_picket_fence_statistics_basic(self):
        """
        Test that the sequence is VERY evenly spaced
        We can't test the getNNSDist code because a picket fence NNSDist is a delta-function at aveD
        """
        self.assertAlmostEqual(numpy.average(self.picketFenceSequenceAnalyzer.spacings), self.expectedAverageD)
        self.assertAlmostEqual(numpy.std(self.picketFenceSequenceAnalyzer.spacings), 0.0)
        self.assertAlmostEqual(numpy.average(self.picketFenceSequenceAnalyzer.normalized_spacings), 1.0)
        self.assertAlmostEqual(numpy.std(self.picketFenceSequenceAnalyzer.normalized_spacings), 0.0)

    def test_picket_fence_statistics_number_levels(self):
        self.assertEqual(self.picketFenceSequenceAnalyzer.num_levels, self.expectedNLevels)

    def test_picket_fence_statistics_rho(self):
        self.assertRaises(ValueError,
                          self.picketFenceSequenceAnalyzer.get2SpacingCorrelation)  # Is undefined for picket fence

    def test_picket_fence_statistics_D(self):
        D = self.picketFenceSequenceAnalyzer.ave_spacing
        self.assertWithinNSigma(D, self.expectedAverageD, 0.00001 * D,
                                additional_msg="; for aveD")  # Should be what we put in

    @unittest.skip("FIXME")
    def test_picket_fence_statistics_U(self):
        U = self.picketFenceSequenceAnalyzer.getThermodynamicEnergyU()
        self.assertWithinNSigma(U, 0.0, 0.0000001, additional_msg="; for U")  # Should be zero for picket fence

    @unittest.skip("FIXME")
    def test_picket_fence_statistics_Q(self):
        Q = self.picketFenceSequenceAnalyzer.getDysonMehtaQ()
        self.assertWithinNSigma(Q, 0.0, 0.0001, additional_msg="; for Q")

    def test_picket_fence_statistics_Delta3(self):
        """
        Compare Delta3 statistic for a picket fence to the expected systematics

        10% agreement is lowsy, but the agreement improves as L-> big.
        If I wanted to be fancy, I could make an L dependent tolerance
        """
        testDelta3Values = self.picketFenceSequenceAnalyzer.getDysonMehtaDelta3_vs_L()
        for L in range(len(testDelta3Values[0])):  # here L is the length of the list of unfolded levels
            if L < self.picketFenceSequenceAnalyzer.Lmin:
                continue
            if L == self.picketFenceSequenceAnalyzer.Lmax:
                continue
            delta3 = testDelta3Values[0][L]
            stddev_delta3 = testDelta3Values[1][L]
            expected_delta3 = getDysonMehtaDelta3Systematics_for_PicketFence(L)
            self.assertWithinNSigma(delta3, expected_delta3, 0.1 * expected_delta3, additional_msg="; for L=%i" % L)


@unittest.skipIf(not TESTALL, "")
class TestLevelAnalyzer_90Zr(AssertMixIn):

    nsamples_big = 10000
    nsamples_medium = 1000
    nsamples_small = 100
    nsamples_vsmall = 10
    nsamples = nsamples_small

    def setUp(self):
        self.E0 = PQU('1 eV')
        self.D = PQU('10 eV')
        self.unit = 'eV'
        self.expectedAverageD = TestLevelAnalyzer_90Zr.aveD  # 3581.25583223
        self.expectedNLevels = TestLevelAnalyzer_90Zr.numLevels  # 442
        self.expectedLinearCorrelationCoefficient = -0.27
        self.expectedAverageU = 0.365
        self.expectedAverageQ = 0.340

    @classmethod
    def setUpClass(cls):
        """
        Expensive set of tests that generate stochastic realizations of resonances
        with full GOE goodness, then computes various metrics and compares them
        to expected results.

        Data for testing is from 90Zr, l=0, j=1/2 URR set

        :return: None
        """
        cls._testData = """
        <!--   energy | levelSpacing | neutronWidth | captureWidth  -->
               200000.0        8655.91      0.5280105      0.1416092
               300000.0       7811.708      0.4765142      0.1476385
               400000.0       7054.611      0.4303313      0.1538193
               500000.0       6375.102      0.3888812      0.1601527
               600000.0       5764.775      0.3516513      0.1666401
               700000.0       5216.188      0.3181875      0.1732827
               800000.0       4722.749      0.2880877      0.1800816
               900000.0       4278.602      0.2609947      0.1870381
              1000000.0       3878.548      0.2365914      0.1941535
              1100000.0       3517.969      0.2145961      0.2014288
              1200000.0       3192.751      0.1947578      0.2088653
              1300000.0       2899.241      0.1768537      0.2164641
              1400000.0        2634.18       0.160685      0.2242264
              1500000.0        2394.66      0.1460743      0.2321535
              1600000.0       2178.092      0.1328636      0.2402462
              1700000.0       1982.151      0.1209112       0.248506
              1780460.0       1838.063      0.1121218      0.2552738"""

        # Set up level density and spacing distribution
        energies = []
        cls.meanSpacing = []
        for row in cls._testData.split('\n')[2:]:
            srow = row.split()
            energies.append(float(srow[0]))
            cls.meanSpacing.append(float(srow[1]))
        meanSpacingFunc = makeXYs(
            data=[energies, cls.meanSpacing],
            xLabel="Ex", yLabel="D(Ex)", yUnit="eV", dataForm="xsandys")
        cls.levelDensity = meanSpacingFunc.applyFunction(
            lambda x, par: 1.0 / meanSpacingFunc.evaluate(x),
            None)  # applyFunc supposed to operate on y values, but really works on x
        cls.Emin = cls.levelDensity.domainMin
        cls.Emax = cls.levelDensity.domainMax
        cls.numLevels = int(cls.levelDensity.integrate().getValue())
        cls.aveD = (cls.Emax - cls.Emin) / cls.numLevels

        # Set up plotting data (if needed)
        if DEBUG:
            # The realizations
            cldOutputFile = open('CLD_calc.txt', mode='w')
            nnsdOutputFile = open('NNSD_calc.txt', mode='w')
            delta3OutputFile = open('Delta3_calc.txt', mode='w')

        # Variables for storing realizations for later analysis
        cls._testCLDValues = {E: [] for E in ["2.3e5", "7.5e5", "1.3e6"]}
        cls._testDValues = []
        cls._testUValues = []
        cls._testQValues = []
        cls._testNumberLevels = []
        cls._testNNSDValues = {D: [] for D in ["0.5", "1.0", "1.5", "2.0"]}  # really test values are D/aveD
        cls._testDelta3Values = {L: [] for L in [10, 20, 30]}
        cls._testDelta3StdValues = {L: [] for L in [10, 20, 30]}
        cls._testLinearCorrelationCoefficientValues = []

        # Make many realizations
        for i in range(cls.nsamples):

            # Set up the level scheme analyzer
            # In order to compare to systematics, we don't use the E0 offset
            levelSequence = level_generator.getFakeLevelSequence(style='goe',
                                                                 levelDensity=cls.levelDensity)  # don't need Emin since levelDensity takes care of specifying the start
            sequenceAnalyzer = LevelSequenceAnalyzer(levels=levelSequence, level_density=cls.levelDensity)

            # Cumulative level distribution testing
            # In order to compare to systematics, we don't assume a ground state level or
            # any other levels below the start of the level sequence (so c0=1, the default value)
            thisCLD = sequenceAnalyzer.getCumulativeLevelDistribution(domainMin=cls.Emin, domainMax=cls.Emax)

            for sE in cls._testCLDValues:
                E = float(sE)
                if E < cls.Emin or E > cls.Emax: 
                    raise ValueError("Test energy out of bounds of level density")
                if E < thisCLD.domainMin or E > thisCLD.domainMax: 
                    raise ValueError("Test energy out of bounds of cummulative level distribution")
                cls._testCLDValues[sE].append(thisCLD.evaluate(E))
            if DEBUG:
                cldOutputFile.write('\n'.join([str(x[0]) + ' ' + str(x[1]) for x in thisCLD]) + '\n')
                cldOutputFile.write('&\n')

            # Mean level spacing testing
            cls._testDValues.append(sequenceAnalyzer.mean_spacing)

            # Total number of level testing
            cls._testNumberLevels.append(sequenceAnalyzer.num_levels)

            # Delta3(L) testing
            dmList, dmVarList = sequenceAnalyzer.getDysonMehtaDelta3_vs_L()
            for L in cls._testDelta3Values:
                cls._testDelta3Values[L].append(dmList[L])
                cls._testDelta3StdValues[L].append(dmVarList[L])
            if DEBUG:
                delta3OutputFile.write(
                    '\n'.join(
                        [' '.join([str(x).ljust(20) for x in (L, dmList[L], math.sqrt(dmVarList[L]))]) for L in
                         range(5, 50)]))
                delta3OutputFile.write('\n&\n')

            # NNSD testing
            thisNNSD = sequenceAnalyzer.getNNSDist(
                normalizeByMeanLevelSpacing=True,
                normalizeDistribution=True,
                numBins=None)
            for sD_over_aveD in cls._testNNSDValues:
                D_over_aveD = float(sD_over_aveD)
                cls._testNNSDValues[sD_over_aveD].append(
                    makeXYs(thisNNSD[0], interpolation=standards.interpolation.flatToken).evaluate(D_over_aveD))
            if DEBUG:
                nnsdOutputFile.write('\n'.join([str(x[0]) + ' ' + str(x[1]) for x in thisNNSD]) + '\n')
                nnsdOutputFile.write('&\n')

            # Linear Correlation Coefficient testing
            cls._testLinearCorrelationCoefficientValues.append(sequenceAnalyzer.get2SpacingCorrelation())

            # Thermodynamic energy U testing
            cls._testUValues.append(0.0)  # FIXME: sequenceAnalyzer.getThermodynamicEnergyU())

            # Q statistic testing
            cls._testQValues.append(0.0)  # FIXME: sequenceAnalyzer.getDysonMehtaQ())

        # Set up plotting systematics/results (if needed)
        if DEBUG:
            # Systematics/original data (the answers)
            open('CLD_orig.txt', mode='w').write(
                '\n'.join([str(x[0]) + ' ' + str(x[1]) for x in cls.levelDensity.indefiniteIntegral()]) + '\n')
            open('NNSD_sys.txt', mode='w').write(
                '\n'.join([str(x[0]) + ' ' + str(x[1]) for x in WignerDistribution()]) + '\n')
            open('Delta3_sys_goe.txt', mode='w').write(
                '\n'.join(
                    [str(L) + ' ' + str(getDysonMehtaDelta3Systematics_for_GOE(L)) for L in range(5, 50)]) + '\n')
            open('Delta3_sys_poisson.txt', mode='w').write(
                '\n'.join([str(L) + ' ' + str(getDysonMehtaDelta3Systematics_for_Poisson(L)) for L in
                           range(5, 50)]) + '\n')
            open('Delta3_sys_sho.txt', mode='w').write(
                '\n'.join([str(L) + ' ' + str(getDysonMehtaDelta3Systematics_for_PicketFence(L)) for L in
                           range(5, 50)]) + '\n')

    # ----------------------
    # The actual unit tests
    # ----------------------

    def test_number_levels(self):
        """Test NLevels"""
        self.assertWithinNSigma(
            numpy.average(TestLevelAnalyzer_90Zr._testNumberLevels),
            self.expectedNLevels,
            numpy.std(TestLevelAnalyzer_90Zr._testNumberLevels),
            additional_msg="; for number of levels")

    def test_rho(self):
        """Test Linear Correlation Coefficient"""
        self.assertWithinNSigma(
            numpy.average(TestLevelAnalyzer_90Zr._testLinearCorrelationCoefficientValues),
            self.expectedLinearCorrelationCoefficient,
            numpy.std(TestLevelAnalyzer_90Zr._testLinearCorrelationCoefficientValues),
            additional_msg="; for rho")

    @unittest.expectedFailure
    def test_U(self):
        """Test U FIXME U messed up a little"""
        self.assertWithinNSigma(
            numpy.average(TestLevelAnalyzer_90Zr._testUValues),
            self.expectedAverageU,
            numpy.std(TestLevelAnalyzer_90Zr._testUValues),
            additional_msg="; for aveU")

    @unittest.expectedFailure
    def test_Q(self):
        """Test Q FIXME Q hosed!"""
        if DEBUG:
            print("testing Q")
            print(TestLevelAnalyzer_90Zr._testQValues)
            print(numpy.average(TestLevelAnalyzer_90Zr._testQValues), '+/-',
                  numpy.std(TestLevelAnalyzer_90Zr._testQValues))
            print(self.expectedAverageQ)
        self.assertWithinNSigma(
            numpy.average(TestLevelAnalyzer_90Zr._testQValues),
            self.expectedAverageQ,
            numpy.std(TestLevelAnalyzer_90Zr._testQValues),
            additional_msg="; for aveQ")

    def test_D(self):
        """Test D, should be nearly exactly dE/numLevels with almost zero variance"""
        self.assertWithinNSigma(
            numpy.average(TestLevelAnalyzer_90Zr._testDValues),
            self.expectedAverageD,
            5.0,  # should be within 5 eV...
            additional_msg="; for aveD")
        self.assertWithinNSigma(
            numpy.std(TestLevelAnalyzer_90Zr._testDValues),
            0.0,
            0.000001,
            additional_msg="; for aveD")

    def test_CLD(self):
        """Test CLD"""
        expectedCLD = TestLevelAnalyzer_90Zr.levelDensity.indefiniteIntegral()
        for sE in TestLevelAnalyzer_90Zr._testCLDValues:
            E = float(sE)
            self.assertWithinNSigma(
                numpy.average(TestLevelAnalyzer_90Zr._testCLDValues[sE]),
                expectedCLD.evaluate(E),
                numpy.std(TestLevelAnalyzer_90Zr._testCLDValues[sE]),
                additional_msg="; for E=%f in CLD calculation" % E)

    def test_NNSD(self):
        """Test NNSD FIXME expected answer function is broke?"""
        for sD_over_aveD in TestLevelAnalyzer_90Zr._testNNSDValues:
            D_over_aveD = float(sD_over_aveD)
            self.assertWithinNSigma(
                numpy.average(TestLevelAnalyzer_90Zr._testNNSDValues[sD_over_aveD]),
                WignerDistribution().evaluate(D_over_aveD),
                numpy.std(TestLevelAnalyzer_90Zr._testNNSDValues[sD_over_aveD]),
                additional_msg="; for D/aveD=%f in NNSD calculation" % D_over_aveD)

    def test_Delta3(self):
        """Test Delta3"""
        for L in TestLevelAnalyzer_90Zr._testDelta3Values:
            self.assertWithinNSigma(
                numpy.average(TestLevelAnalyzer_90Zr._testDelta3Values[L]),
                getDysonMehtaDelta3Systematics_for_GOE(L),
                numpy.std(TestLevelAnalyzer_90Zr._testDelta3Values[L]),
                additional_msg="; for L=%f in Delta3 calculation" % L)

    @unittest.skip("Following test really doesn't belong in this test suite")
    def test_resonance_generation(self):
        """
        This test uses the full machinery of the fake resonance region generation
        :return:
        """

        SHOWFISSION = False
        SHOWCOMPETITIVE = False
        DOPRINTRESONANCES = False

        # This first part here does the problem setup, just like it would be done in a production code
        # I'll indent it to set it off as a logical unit
        if True:
            neutronDOF = 1
            captureDOF = 0

            # Hard coded default fudge width keys (a bad choice)
            widthKeys = ('neutronWidth', 'captureWidth')  # ,'fissionWidthA','competitiveWidth')

            # Set up the variables
            energies = []
            levelSpacings = []
            widths = {}
            widthData = {widthKey: [] for widthKey in widthKeys}

            # We know the DOFs
            DOFs = {widthKey: 0 for widthKey in widthKeys}
            DOFs['neutronWidth'] = neutronDOF

            # parse the test data
            for row in self.testData.split('\n')[2:]:
                srow = row.split()
                energies.append(float(srow[0]))
                levelSpacings.append(float(srow[1]))
                widthData['neutronWidth'].append(float(srow[2]))
                widthData['captureWidth'].append(float(srow[3]))
                if SHOWFISSION: widthData['fissionWidthA'].append(0.0)
                if SHOWCOMPETITIVE: widthData['competitiveWidth'].append(0.0)

            # Turn the widths into interpolatable things
            for widthKey in widthKeys:
                widths[widthKey] = makeXYs(
                    data=[energies, widthData[widthKey]],
                    xLabel="Ex", yLabel=widthKey, yUnit="eV", dataForm="xsandys")

            # define the level density/mean level spacing
            if True:
                levelDensityData = [energies, [1.0 / x for x in levelSpacings]]  # real spacing distribution
            else:
                levelDensityData = [energies, [1.0 / levelSpacings[0] for x in levelSpacings]] # constant spacing dist.
            levelDensity = makeXYs(
                data=levelDensityData,
                xLabel="Ex", yLabel="rho(Ex)", yUnit="1/eV", dataForm="xsandys")

        # Set up plotting data (if needed)
        if self.debug:
            # The realizations
            nwidthsOutputFile = open('nwidths_calc.txt', mode='w')

        # Variables for storing realizations for later analysis
        testNeutronWidthValues = {G: [] for G in [4.0e-05, 0.0001, 0.0003]}

        # We are going to make a bunch of realizations, then do statistics on them
        for iRealization in range(self.nsamples_vsmall):
            # make the resonances
            resTable = resonance_generator.getFakeResonanceSet(
                E0=0.0,  # for GOE tests, better set this to 0.0
                style='goe',
                levelDensity=levelDensity,
                aveWidthFuncs=widths,
                L=0.0, J=0.5,
                DOFs=DOFs,
                widthKeys=widthKeys)

            # assume tests of generated energies are kosher,
            # so focus on tests of generated widths
            nwidthList = [x[3] / self.expectedNLevels for x in resTable.data]
            nhisto = makeHistoXYs(nwidthList, 15)
            for G in testNeutronWidthValues:
                testNeutronWidthValues[G].append(nhisto.evaluate(G))
            if self.debug:
                nwidthsOutputFile.write(
                    '\n'.join([str(x[0]) + ' ' + str(x[1]) for x in nhisto]) + '\n')
                nwidthsOutputFile.write('&\n')

        if DOPRINTRESONANCES:
            print('\n'.join(resTable.toXMLList()))

        eMin = widths['neutronWidth'].domainMin
        eMax = widths['neutronWidth'].domainMax
        aveWidth = widths['neutronWidth'].integrate() / abs(eMax - eMin)

        # Set up plotting systematics/results (if needed)
        if self.debug:
            # Systematics/original data (the answers)
            open('nwidths_orig.txt', mode='w').write(
                '\n'.join(
                    [
                        str(2.0 * aveWidth * x[0] / DOFs['neutronWidth']) + ' ' + str(
                            x[1] * (2.0 * aveWidth / DOFs['neutronWidth'])) for x in
                        Chi2Distribution(DOFs['neutronWidth']) if x[0] > 0.0]) + '\n')

        # ----------------------
        # The actual unit tests
        # ----------------------

        # Test neutron widths
        for G in testNeutronWidthValues:
            self.assertWithinNSigma(
                numpy.average(testNeutronWidthValues[G]),
                Chi2Distribution(DOFs['neutronWidth']).evaluate(G) * (2.0 * aveWidth / DOFs['neutronWidth']),
                numpy.std(testNeutronWidthValues[G]), additional_msg="; for Gamma=%f" % G)


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
        cProfile.run("unittest.main()", "blurr__init__.profile")
    else:
        unittest.main()
        print()
