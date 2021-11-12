import math
import numpy
import scipy.optimize
from xData import axes, XYs, standards
from brownies.BNL.utilities.fitting import bounded_nonlinear_curve_fit

# ---------------------------------------------------------------------------------
#
#  Level sequence analysis
#
# ---------------------------------------------------------------------------------


def getDysonMehtaDelta3Systematics_for_GOE(L):
    return (math.log(2.0 * math.pi * L) + (0.577215665) - 5. / 4. - math.pi * math.pi / 8.0) / math.pi / math.pi


def getDysonMehtaDelta3Systematics_for_Poisson(L):
    return L / 15.0


def getDysonMehtaDelta3Systematics_for_PicketFence(L):
    return 1.0 / 12.0


class LevelSequenceAnalyzer:

    def __init__(self,
                 levels,
                 level_density=None,
                 energyUnit='eV',
                 Lmin=5,
                 Lmax=50,
                 warn_if_no_level_density=False):
        self.levels = levels
        self.levels.sort()
        self.energyUnit = energyUnit
        if level_density is None:
            if warn_if_no_level_density:
                print("WARNING: Level density not provided, computing from mean spacing "
                      "computed from entire list of levels assuming that the mean spacing "
                      "is constant as a function of energy.")
            self.level_density = self.getAverageSpacingLevelDensity()
        else:
            self.level_density = level_density
        if self.level_density.domainMin > min(self.levels):
            raise ValueError("Level sequence has value %s which is smaller than level density's domainMin of %s" %
                             (str(min(self.levels)), str(self.level_density.domainMin)))
        if self.level_density.domainMax < max(self.levels):
            raise ValueError("Level sequence has value %s which is larger than level density's domainMax of %s" %
                             (str(max(self.levels)), str(self.level_density.domainMax)))
        self.Lmin = Lmin
        self.Lmax = Lmax

    @property
    def num_levels(self):
        return len(self.levels)

    @property
    def unfolded_levels(self):
        theCLD = self.level_density.indefiniteIntegral()
        return [theCLD.evaluate(x) for x in self.levels]

    @property
    def spacings(self):
        return [self.levels[i + 1] - self.levels[i] for i in range(self.num_levels - 1)]

    @property
    def num_spacings(self):
        return len(self.spacings)

    @property
    def ave_spacing(self):
        return numpy.average(self.spacings)

    @property
    def mean_spacing(self):
        return numpy.mean(self.spacings)

    @property
    def stddev_spacing(self):
        return numpy.std(self.spacings)

    @property
    def var_spacing(self):
        return numpy.var(self.spacings)

    @property
    def normalized_spacings(self):
        result = []
        for i in range(self.num_levels - 1):
            D = self.levels[i + 1] - self.levels[i]
            E = 0.5 * (self.levels[i + 1] + self.levels[i])
            result.append(D / self.mean_spacing_at_E(E))
        return result

    def mean_spacing_at_E(self, E):
        try:
            return 1.0 / self.level_density.evaluate(E)
        except TypeError as err:
            raise TypeError(str(err) + " for E=" + str(E))

    def getAverageSpacingLevelDensity(self):
        """
        Generate a level density corresponding to the average spacing

        :return: An XYs1d containing the level density
        """
        myAxes = axes.axes(rank=2, labelsUnits={0: ('Level density', '1/' + self.energyUnit),
                                                1: ('Excitation energy', self.energyUnit)})
        theData = []
        for thisEnergy in self.levels:
            if len(theData) > 0 and abs(thisEnergy - theData[-1][0]) < 1e-6 * thisEnergy:
                continue  # Deal with corner case where two levels are degenerate
            else:
                theData.append([thisEnergy, 1.0 / self.mean_spacing])
        return XYs.XYs1d(axes=myAxes, data=theData)

    def getCumulativeLevelDistribution(self, c0=1.0, domainMin=None, domainMax=None):
        """
        Generate a cummulative level distribution from the levels in self.levels

        :param c0: 1st bin's starting number of levels (if this one is tacked onto a lower energy sequence)
        :param domainMin:
        :param domainMax:
        :return: An XYs1d containing the histogram of the cummulative level distribution
        """
        myAxes = axes.axes(rank=2, labelsUnits={0: ('Number of levels', ''), 1: ('Excitation energy', self.energyUnit)})
        theData = []
        if domainMin is not None and domainMin < self.levels[0]:
            theData.append([domainMin, c0])
            theData.append([self.levels[0], c0 + 1])
        else:
            theData.append([self.levels[0], c0])
        for thisEnergy in self.levels[1:]:
            if abs(thisEnergy - theData[-1][0]) < 1e-6 * thisEnergy:
                theData[-1][1] = theData[-1][1] + 1  # Deal with corner case where two levels are degenerate
            else:
                theData.append([thisEnergy, theData[-1][1] + 1])
        if domainMax is not None and domainMax > self.levels[-1]:
            theData.append([domainMax, theData[-1][1]])
        return XYs.XYs1d(axes=myAxes, data=theData, interpolation=standards.interpolation.flatToken).domainSlice(
            domainMin=domainMin, domainMax=domainMax)

    def getNNSDist(self, normalizeByMeanLevelSpacing=True, normalizeDistribution=True, numBins=None):
        if normalizeByMeanLevelSpacing:
            if numBins is not None:
                raw_histogram = numpy.histogram(self.normalized_spacings, bins=numBins, density=normalizeDistribution)
            else:
                raw_histogram = numpy.histogram(self.normalized_spacings, density=normalizeDistribution)
        else:
            if numBins is not None:
                raw_histogram = numpy.histogram(self.spacings, bins=numBins, density=normalizeDistribution)
            else:
                raw_histogram = numpy.histogram(self.spacings, density=normalizeDistribution)
        table = []
        for i in range(len(raw_histogram[0])):
            table.append([0.5 * (raw_histogram[1][i] + raw_histogram[1][i + 1]), raw_histogram[0][i]])
        return table, raw_histogram

    def getDysonMehtaDelta3_vs_L(self):
        if self.num_levels < self.Lmin:
            raise ValueError(
                "Not enough levels (%i) to compute Delta_3, need at least %i" % (self.num_levels, self.Lmin))
        Lmax = min(self.num_levels, self.Lmax)
        P_L = [0.0 for L in range(Lmax + 1)]
        P2_L = [0.0 for L in range(Lmax + 1)]
        N_L = [0.0 for L in range(Lmax + 1)]
        levelSequence = self.unfolded_levels
        # Get all the L-dists for each energy
        for iE in range(self.num_levels - Lmax):
            for iL, p in enumerate(self.getDysonMehtaDelta3_vs_L_at_E(iE, levelSequence)):
                P_L[iL + self.Lmin] += p
                P2_L[iL + self.Lmin] += p * p
                N_L[iL + self.Lmin] += 1
        # Normalize the results
        for L in range(Lmax + 1):
            if N_L[L] > 0.0:
                P_L[L] /= N_L[L]
        # Compute variance
        for L in range(Lmax + 1):
            if N_L[L] > 1.0:
                P2_L[L] = (P2_L[L] - P_L[L] * P_L[L]) / N_L[L] / (N_L[L] - 1)
        return P_L, P2_L

    def getDysonMehtaDelta3_vs_L_at_E(self, iE, levelSequence):
        if self.num_levels < self.Lmin:
            raise ValueError(
                "Not enough levels (%i) to compute Delta_3, need at least %i" % (self.num_levels, self.Lmin))
        Lmax = min(self.num_levels, self.Lmax)
        return [self.getDysonMehtaDelta3(L, iE, levelSequence) for L in range(self.Lmin, Lmax + 1)]

    def getDysonMehtaDelta3(self, L, iE, levelSequence):
        """
        Compute the Dyson Mehta Delta3 statistic using Declan's algorithm (D. Mulhall, Phys. Rev. 83, 05321 (2011))

        :return:
        """
        C, W, Y = 0.0, 0.0, 0.0
        Z = (levelSequence[iE + L] - levelSequence[iE])
        X = Z * (levelSequence[iE + L] + levelSequence[iE])
        V = (pow(levelSequence[iE + L], 3.0) - pow(levelSequence[iE], 3.0)) / 3.0
        for j in range(iE, iE + L):
            jEdiff = j * (levelSequence[j + 1] - levelSequence[j])
            C += j * jEdiff
            W += -jEdiff * (levelSequence[j + 1] + levelSequence[j])
            Y += -2.0 * jEdiff
        den = (4.0 * V * Z - X * X)
        A = (X * Y - 2.0 * W * Z) / den
        B = (W * X - 2.0 * V * Y) / den
        return (C + V * A * A + W * A + X * A * B + Y * B + Z * B * B) / (L + 1)

    def get2SpacingCorrelation(self, epsD=1e-9):
        """
        Get the spacing-spacing correlation function.

            ..math::
                D_i=E_{i+1}-E_i

            ..math::

                \rho(D_i,D_{i+1}) = \frac{\sum_i(D_i-\overline{D})(D_{i+1}-\overline{D})}{\left[\sum_i
                                    (D_i-\overline{D})^2\sum_j(D_{j+1}-\overline{D})^2\right]^{1/2}}

        If you stare at this for a while, it is just the level-level level spacing correlation coefficient.

        For GOE level spacings, should see :math:`\rho=-0.27`

        :returns : value and variance in a tuple
        """
        diff_spacings = []
        for i in range(self.num_levels - 1):
            aveE = 0.5 * (self.levels[i + 1] + self.levels[i])
            diff_spacings.append(self.spacings[i] - self.mean_spacing_at_E(aveE))
        if len([x for x in diff_spacings if x > epsD]) < 2:
            raise ValueError("Level-level spacing correlation undefined, is this a picket fence?")
        correlation_matrix = numpy.corrcoef(
            numpy.array([[diff_spacings[i], diff_spacings[i + 1]] for i in range(self.num_spacings - 1)]).T)
        return correlation_matrix[0, 1]

    def getThermodynamicEnergyU(self):
        """
        Get Dyson-Mehta's estimate of the thermodynamic energy of the system, U

        Use Eqs. (19-21) from G.E. Mitchell, J.F. Shriner, "Missing Level Corrections using Neutron Spacings",
        IAEA NDS Report INDC(NDS)-0561 (2009)
        """
        L = self.num_levels / 2.0

        # normalize energies to average spacing
        # dE=self.energies[-1]-self.energies[0]
        # eps=[(2*L-1)*E/dE for E in self.levels] # normalize to average spacing (option 1)
        eps = [E / self.mean_spacing for E in self.levels]  # normalize to average spacing (option 2)

        # center energies on interval
        E0 = eps[0] + L - 0.5
        eps = [E - E0 for E in eps]

        # "picket fence" energies
        peps = [-L + (i + 1) - 0.5 for i in
                range(self.num_levels)]  # note the "1" index offset -- Mitchell and Shriner use Fortran style indexing

        # define "potential energy"
        def V(E):
            return (L - E) * (0.5 + math.log((L - E) / (2.0 * L))) + (L + E) * (0.5 + math.log((L + E) / (2.0 * L)))

        # Now compute the U statistic
        term1 = 0.0
        term2 = 0.0
        for i in range(self.num_levels):
            term2 += V(eps[i]) - V(peps[i])
            for j in range(i + 1, self.num_levels):
                term1 += math.log((eps[j] - eps[i]) / (j - i))
        return -(term1 + term2) / self.num_levels

    def getDysonMehtaQ(self):
        """
        Get Dyson-Mehta's Q statistic (billed as a lower-variance version of the thermodynamic energy of the system)

        Use Eqs. (25-32) from G.E. Mitchell, J.F. Shriner, "Missing Level Corrections using Neutron Spacings",
        IAEA NDS Report INDC(NDS)-0561 (2009)
        """
        # Mitchell and Shriner pick M=4, so ...
        M = 4
        R = M * self.mean_spacing

        # renormalize energies
        Eave = 0.5 * (self.levels[0] + self.levels[-1])
        eps = [E - Eave for E in self.levels]

        # compute scaling variable
        Emax = max(self.levels)  # FIXME: check this, not clearly defined in section D
        Emin = min(self.levels)  # FIXME: check this, not clearly defined in section D
        if Emax - self.levels[-1] < self.levels[0] - Emin:
            L = 0.5 * (self.levels[-1] + Emax)
        else:
            L = 0.5 * (Emin + self.levels[0])

        def f(x, y):
            """Eq. (29)"""
            if abs(x - y) < R and abs(x) < L and abs(y) < L:
                return 1
            return 0

        # Eq. (30)
        U0 = -R / L + R * R / 8.0 / L / L

        def F(x):
            """Eq. (30)"""
            if abs(x) < L:
                return 1
            return 0

        def U(x):
            """Eq. (30)"""
            if abs(x) < L - R:
                return -R / L
            if L > abs(x) > L - R:
                return -(R + (L - abs(x)) * (1.0 - math.log((L - abs(x)) / R))) / 2.0 / L
            raise ValueError("shouldn't get here, x=%s, L=%s, R=%s" % (str(x), str(L), str(R)))

        # Compute Q statistic
        term1 = 0.0
        term2 = 0.0
        term3 = 0.0
        term4 = 0.0
        sumF = 0.0
        intF = 2.0 * L
        for i in range(self.num_levels):
            sumF += F(eps[i])
            for j in range(self.num_levels):
                if i < j:
                    term1 -= f(eps[i], eps[j]) * math.log((eps[j] - eps[
                        i]) / R)  # suspect typo in term #1 of Eq. (25), i & j should be swapped in logrithm
                term2 += F(eps[i]) * U(eps[j])
                term3 -= 0.5 * U0 * F(eps[i]) * F(eps[j])
        term4 = 0.5 * sumF * math.log(2.0 * math.pi * R * sumF / intF)
        n = 2.0 * L / self.mean_spacing
        return (term1 + term2 + term3 + term4) / n

    def getFractionMissing(self):
        """Estimate the fraction of missing levels using an inference based on all of the metrics given below"""
        raise NotImplementedError("write me")

    def getFractionMissingFromDysonMehtaQ(self):
        """
        See G.E. Mitchell, J.F. Shriner, "Missing Level Corrections using Neutron Spacings",
        IAEA NDS Report INDC(NDS)-0561 (2009)
        """
        raise NotImplementedError("write me")

    def getFractionMissingFromDysonMehtaDelta3(self, use_b_covariance=False):
        """
        See G.E. Mitchell, J.F. Shriner, "Missing Level Corrections using Neutron Spacings",
        IAEA NDS Report INDC(NDS)-0561 (2009)
        """
        d3, dd3 = self.getDysonMehtaDelta3_vs_L()
        Ls = numpy.array(list(range(self.Lmin, self.Lmax + 1)))
        d3 = numpy.array(d3[5:])
        dd3 = numpy.array(dd3[5:])
        bmean = numpy.array([0.920, 0.978])
        bcov = numpy.array([[6.64e-3, -9.99e-4], [-9.99e-4, 4.15e-3]])

        def d3func(L, *f, **kwrds):
            f1 = 1. - f[0]
            LL = L  # /2.0  #don't need factor of 2 used in Mitchell and Shriner Eq. (39)
            return kwrds['oparms'][0] * f[0] * LL / 15.0 + kwrds['oparms'][1] * f1 * f1 * (
                        numpy.log(LL / f1) - 0.687) / math.pi / math.pi

        if use_b_covariance:
            result = bounded_nonlinear_curve_fit(func=d3func, xdata=Ls, ydata=d3,
                                                 yerr=dd3, fitparms0=[0.8], oparms=bmean, oparmscov=bcov,
                                                 bounds=([0.0], [1.0]), size=100)
        else:
            result = scipy.optimize.curve_fit(lambda L, *f: d3func(L, *f, oparms=bmean), Ls, d3, p0=[0.8], sigma=dd3,
                                              bounds=([0.0], [1.0]))
        return float(result[0]), math.sqrt(float(result[1]))

    def getFractionMissingThermodynamicEnergyU(self):
        """
        See G.E. Mitchell, J.F. Shriner, "Missing Level Corrections using Neutron Spacings",
        IAEA NDS Report INDC(NDS)-0561 (2009)
        """
        raise NotImplementedError("write me")

    def getFractionMissing2SpacingCorrelation(self):
        """
        See G.E. Mitchell, J.F. Shriner, "Missing Level Corrections using Neutron Spacings",
        IAEA NDS Report INDC(NDS)-0561 (2009)
        """
        rho = self.get2SpacingCorrelation()
        emean = numpy.array([-0.251, 0.428])
        ecov = numpy.array([[2.67e-4, -7.44e-4], [-7.44e-4, 3.22e-3]])
        dfde0 = -1.0 / emean[1]
        dfde1 = -(rho - emean[0]) / emean[1] / emean[1]
        f = (rho - emean[0]) / emean[1]
        df = math.sqrt(dfde0 * dfde0 * ecov[0][0] + dfde1 * dfde1 * ecov[1][1] + 2.0 * dfde0 * dfde1 * ecov[0][1])
        return max(min(1.0, f), 0.0), df

    def getFractionMissingNNSDist(self):
        """
        See G.E. Mitchell, J.F. Shriner, "Missing Level Corrections using Neutron Spacings",
        IAEA NDS Report INDC(NDS)-0561 (2009)
        """
        raise NotImplementedError("write me")

