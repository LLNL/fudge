from pqu.PQU import PhysicalQuantityWithUncertainty as PhysicalQuantity
from xData import axes, XYs

cSquared = PhysicalQuantity(1.00000000000000000000, 'c**2')
hbar = PhysicalQuantity("6.58211928e-16", "eV*s")
hbarc = PhysicalQuantity("197.32697", "MeV*fm")
e2 = PhysicalQuantity("1.43998e-10", "keV*cm")
muN2 = PhysicalQuantity("1.59234e-38", "keV*cm**3")
b = PhysicalQuantity('1e-24', 'cm**2')


# ---------------------------------------------------------------------------------
#
# Basic data types
#
# ---------------------------------------------------------------------------------

class LevelDensity(XYs.XYs1d):

    def __init__(self, **kw):
        """
        :param kw: Only 'data' and 'energyUnit' keywords are supported.
                   'energyUnit' is the unit of the x-axis of self.
                   'data' is the the list of (x,y) pairs to initialize the XYs1d that really is self.
        """
        energyUnit = kw.get('energyUnit', 'MeV')
        self.axes = axes.axes(
            labelsUnits={0: ('Level density', '1/' + energyUnit), 1: ('Excitation energy', energyUnit)})
        XYs.XYs1d.__init__(self, axes=self.axes, data=kw.get('data', []))

    def get_mean_level_spacing(self):
        """
        The mean level spacing D is basically the inverse of the level density

        :return: an XYs1d which is 1/(level density)
        """
        return 1.0 / self

    @property
    def total_number_levels(self):
        """
        The total number of levels over the whole domain of self

        :return: The total number of levels, integrated over the whole domain of self
        """
        return int(float(self.integrate().value))

    def get_CLD(self, normalize=False):
        """
        Cummulative level distribution (basically running integral of "self").
        The maximum value of the CLD should equal total_number_levels unless the "normalize" flag is set.

        :param normalize: normalize the CLD to the total number of levels, so the max. value of the CLD is 1.0
        :return: An XYs1d containing the CLD.
        """
        # Normalize level density so is more of cumulative probability distribution
        if normalize:
            return (1.0 / self.total_number_levels) * self.indefiniteIntegral()
        return self.indefiniteIntegral()

    def get_CLD_inverse(self, normalize=False):
        """
        Inverse of the cummulative level distribution.  In other words, if y=CLD(x), invCLD(y)=x.

        If the CLD is not normalized, invCLD's domain is [0,total_number_levels] and the range is the domain of the
        original LevelDensity instance ("self")

        If hte CLD is normalized, the invCLD's domain is [0,1.0], but the range unchanged.

        :param normalize: normalize the CLD to the total number of levels, so the max. value of the CLD is 1.0, before
                          doing anything
        :return: An XYs1d containing the inverse of the CLD.
        """
        # Invert it to make a probability -> energy converter
        return self.get_CLD(normalize=normalize).inverse()

    def evaluate_mean_level_spacing(self, Ex):
        """
        The mean level spacing D is basically the inverse of the level density

        :param Ex: The energy to evaluate the mean level spacing at.  It is assumed to have the same units as the
                   appropriate axis of self
        :return: 1/(level density), evaluated at energy Ex
        """
        return 1.0 / self.evaluate(Ex)


class NuclearLevelDensity:

    def __init__(self, **kw):
        self.axes = axes.axes(labelsUnits={0: ('Level density', '1/MeV'), 1: ('Excitation energy', 'MeV')})
        self.spin_dep_level_density = {1: {}, -1: {}}
        self.associated_quantities = {}

    def add_spin_dependent_level_density(self, data, J, Pi):
        if Pi not in [1, -1]:
            raise ValueError("Parity must be either 1 or -1")
        self.spin_dep_level_density[Pi][self.two_J(J)] = LevelDensity(data=data)

    def two_J(self, J):
        """
        Spins are used often as keys in maps.
        So, within the nuclearLevelDensity heirarchy use 2*J everywhere instead of J.
        That way, 2*J is an integer and can be used as a map key.
        We can also do checking on the value of J.
        """
        if int(2 * J) != 2 * J:
            raise ValueError("J must be integer or half-integer")
        return int(2 * J)

    def get_J_range(self):
        """
        The list of J's returned should be the actual J's, not 2*J
        """
        twoJList = []
        for twoJ in list(self.spin_dep_level_density[1].keys()) + list(self.spin_dep_level_density[-1].keys()):
            if twoJ not in twoJList:
                twoJList.append(twoJ)
        return [x / 2 for x in twoJList]

    def add_associated_quantity(self, data, key, yAxisLabel=None):
        """
        if you want to override the y axis labeling and/or units, do: yAxisLabel=( 'Level density', '1/MeV' )
        """
        if yAxisLabel is None:
            myAxes = self.axes
        else:
            myAxes = axes.axes(labelsUnits={0: yAxisLabel, 1: ('Excitation energy', 'MeV')})
        self.associated_quantities[key] = XYs.XYs1d(axes=myAxes, data=data)

    def get_total_level_density_for_parity(self, Pi):
        s = None
        for J in self.spin_dep_level_density[Pi]:
            if s is None:
                s = self.spin_dep_level_density[Pi][J]
            else:
                s += self.spin_dep_level_density[Pi][J]
        return s

    def get_total_level_density(self):
        return self.get_total_level_density_for_parity(1) + self.get_total_level_density_for_parity(-1)

    def get_CLD(self, J=None, Pi=None):
        if J is not None and Pi is not None:
            return self.spin_dep_level_density[Pi][J].get_CLD()
        if Pi is not None:
            return self.get_total_level_density_for_parity(Pi).get_CLD()
        else:
            return self.get_total_level_density().get_CLD()

    def get_mean_level_spacing(self, J, Pi):
        if J is not None and Pi is not None:
            return self.spin_dep_level_density[Pi][J].get_mean_level_spacing()
        if Pi is not None:
            return self.get_total_level_density_for_parity(Pi).get_mean_level_spacing()
        else:
            return self.get_total_level_density().get_mean_level_spacing()

    def evaluate(self, Ex, J, Pi):
        Pi, twoJ = int(Pi), self.two_J(J)
        if Pi not in self.spin_dep_level_density:
            return 0.0
        if twoJ not in self.spin_dep_level_density[Pi]:
            return 0.0
        return self.spin_dep_level_density[Pi][twoJ].evaluate(Ex) #.getValue_units(str(Ex) + ' MeV')

    def evaluate_mean_level_spacing(self, Ex, J, Pi):
        return 1.0 / self.evaluate(Ex, J, Pi)


class NuclearLevelDensityHFBM(NuclearLevelDensity):

    def __init__(self, **kw):
        NuclearLevelDensity.__init__(self, **kw)
        self.aTilde_correction = kw.get('aTildeCorrection', 0.0)
        self.pairing_shift_correction = kw.get('pairingShiftCorrection', 0.0)
        self.rescaled = False

    def apply_rescaling(self):
        """
        "Irreversible" rescaling of the spin dependent level densities (it is irreversible since I don't feel like coding the reverse operation)
        """
        import math
        # print 'p =', self.pairingShiftCorrection
        # print 'c =', self.aTildeCorrection
        for Pi in self.spin_dep_level_density:
            for J in self.spin_dep_level_density[Pi]:
                newData = []
                for i in range(len(self.spin_dep_level_density[Pi][J])):
                    E, val = self.spin_dep_level_density[Pi][J][i]
                    U = E - self.pairing_shift_correction
                    shiftedVal = self.spin_dep_level_density[Pi][J].evaluate(U)
                    if U > 0 and shiftedVal is not None:
                        newData.append([E, math.exp(self.aTilde_correction * math.sqrt(U)) * shiftedVal])
                    else:
                        newData.append([E, 0.0])
                for i in range(len(self.spin_dep_level_density[Pi][J])):
                    self.spin_dep_level_density[Pi][J][i] = newData[i]


# ---------------------------------------------------------------------------------
#
#    Level scheme readers for different formats
#
# ---------------------------------------------------------------------------------
class ParseError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def readRIPLMeanLevelSpacing(f, verbose=False):
    """
    The resonance parameters can be gotten from RIPL-3 resonances files resonances0.dat and
    resonances1.dat.  They have the following format (from the resonances.readme file):

        #                 Resonance data for s-wave resonances
        # Z El  A   Io     Bn        D0        dD      S0    dS    Gg   dG   Com.
        #                 [MeV]     [keV]     [keV]  [1E-4][1E-4][meV][meV]  &Ref
        #========================================================================
         11 Na  23  1.5   6.960  1.00E+02  2.00E+01   0.25  0.15    -     -  *07I
         12 Mg  24  0.    7.331  4.80E+02  7.00E+01   0.00  0.00    -     -  *07I

        #                  Resonance data for p-wave resonances
        # Z El  A    Io    Bn        D1        dD      S1   dS    Gg    dG    Ref
        #                 [MeV]    [keV]     [keV]  [1E-4][1E-4] [meV] [meV]
        #========================================================================
          9 F   19  0.5   6.601  6.00E+01  1.00E+01   7.50  3.40              06M
         11 Na  23  1.5   6.960  5.00E+01  1.00E+01   1.7   0.5               07I

        The recommended parameters are given separately for s- and p-wave resonances in
        files resonances0.dat and resonances1.dat respectively. The file for s-wave
        resonances contents 300 nuclides and the file for p-wave resonances includes 119
        nuclides.


        Format
        ------
        Each record of the both files contains:

           Z       : charge number of the target nucleus
           A       : mass number of the target nucleus
           El      : element symbol of the target nucleus
           Io      : spin of the ground state of the target nucleus
           Bn      : neutron binding energy for the corresponding compound
                     nucleus in MeV
           D0 or D1: average resonance spacing in keV
           dD      : uncertainty of DL in keV
           S0 or S1: s-wave or p-wave neutron strength function in 10**(-4)
           dSo     : uncertainty of So in 10**(-4)
           Gam     : average radiative width in meV
           dGam    : uncertainty of Gam in meV
           Com&Ref : reference to the work in which the analysis was performed,
                     asterisks marks nuclei for which D0 and dD were estimated from D1

        The corresponding FORTRAN format is
        (i3,1x,a2,1x,i3,2x,f3.1,2x,f6.3,2x,2(e8.2,2x),1x,2(f4.2,2x),2(f4.1,1x),2x,a4).

    """
    result = {}

    lines = f.split('\n')
    while lines:
        line = lines.pop(0)

        if line.strip() == "" or line.startswith("#"):
            continue

        # Read essential data
        sline = line.split()
        Z = int(sline[0])
        A = int(sline[2])
        name = sline[1] + sline[2]
        Io = sline[3]
        Bn = sline[4]
        Dn = sline[5]
        dDn = sline[6]

        # Read optional data; parsing rest of lines requires formatting reading, and we don't need them
        try:
            Sn = sline[7]
            dSn = sline[8]
            Gg = sline[9]
            dGg = sline[10]
            reference = sline[11]
        except:
            pass

        # Now grok the level spacing
        if Dn == '-':
            Dn = None
        else:
            Dn = PhysicalQuantity(float(Dn), "keV")
        if dDn == '-':
            dDn = None
        else:
            dDn = PhysicalQuantity(float(dDn), "keV")
        result[name] = (Dn, dDn)

    return result


def readHFBMLevelDensityTable(f, Z, A, verbose=False):
    """
    From the HFBM readme in RIPL, we have the following lines describing the format.
    These (as near as I can tell) are correct, even if the description of the deformations, etc. are incorrect.

        - 60 data lines including the energy-, spin-dependent level density for
        positive parities
             U     : excitation energy in MeV
             T     : nuclear temperature in MeV
             Ncumul: cumulative number of positive-parity levels
             Rhoobs: total level density in 1/MeV (positive parity)
             Rhotot: total state density in 1/MeV (positive parity)
             Rho   : spin-dependent level density (in 1/MeV) for the first 50 spins (positive parity)

        The corresponding fortran format is (f7.2,f7.3,1x,1p,53e9.2)

        - 60 data lines including the energy-, spin-dependent level density for
        negative parities
             U     : excitation energy in MeV
             T     : nuclear temperature in MeV
             Ncumul: cumulative number of negative-parity levels
             Rhoobs: total level density in 1/MeV (negative parity)
             Rhotot: total state density in 1/MeV (negative parity)
             Rho   : spin-dependent level density (in 1/MeV) for the first 50 spins (negative parity)

        The corresponding fortran format is (f7.2,f7.3,1x,1p,53e9.2)
    """

    negParityToken = "Z=" + str(Z).rjust(3) + " A=" + str(A).rjust(3) + ": Negative-parity"
    posParityToken = "Z=" + str(Z).rjust(3) + " A=" + str(A).rjust(3) + ": Positive-Parity"
    foundPosParity, foundNegParity=False, False

    theTable = {}

    lines = f.split('\n')

    # Read in the data tables
    while lines:
        line = lines.pop(0)

        # Positive parity entries come first
        if posParityToken in line:
            foundPosParity=True
            theTable[1] = []
            lines.pop(0)  # the "****" line
            line = lines.pop(0)  # The column headings
            for colHeading in line.split():
                theTable[1].append([colHeading])
            while line.strip() != '':
                line = lines.pop(0)
                sline = line.split()
                if len(sline) == 0:
                    continue
                nCols = len(theTable[1])
                for icol in range(nCols):
                    theTable[1][icol].append(float(sline[icol]))

        # Negative parity entries come second, we're done once we've read these
        if negParityToken in line:
            foundNegParity=True
            theTable[-1] = []
            lines.pop(0)  # the "****" line
            line = lines.pop(0)  # The column headings
            for colHeading in line.split():
                theTable[-1].append([colHeading])
            while line.strip() != '':
                line = lines.pop(0)
                sline = line.split()
                if len(sline) == 0:
                    continue
                nCols = len(theTable[1])
                for icol in range(nCols):
                    theTable[-1][icol].append(float(sline[icol]))
            break

    if verbose and not foundPosParity:
        print("WARNING: positive parity table not found for Z=%i, A=%i"%(Z,A))

    if verbose and not foundNegParity:
        print("WARNING: negative parity table not found for Z=%i, A=%i"%(Z,A))

    # Make the level density and populate it
    myLevDens = NuclearLevelDensityHFBM()
    for Pi in [1, -1]:
        UList = list(map(float, theTable[Pi][0][1:]))
        myLevDens.add_associated_quantity(data=list(zip(UList, list(map(float, theTable[Pi][1][1:])))), key=(Pi, "T"),
                                          yAxisLabel=("T", "MeV"))
        myLevDens.add_associated_quantity(data=list(zip(UList, list(map(float, theTable[Pi][2][1:])))),
                                          key=(Pi, "NCUMUL"),
                                          yAxisLabel=("Num. Levels", ""))
        myLevDens.add_associated_quantity(data=list(zip(UList, list(map(float, theTable[Pi][3][1:])))),
                                          key=(Pi, "RHOOBS"),
                                          yAxisLabel=("rho", "1/MeV"))
        myLevDens.add_associated_quantity(data=list(zip(UList, list(map(float, theTable[Pi][4][1:])))),
                                          key=(Pi, "RHOTOT"),
                                          yAxisLabel=("rho", "1/MeV"))
        JString = ''
        for col in theTable[Pi][5:]:
            try:
                JString = col[0].lstrip('J=').lstrip('0')
                if '/' in JString:
                    J = float(JString.split('/')[0]) / 2.0
                elif col[0] == "J=00":
                    J = 0
                else:
                    J = float(JString)
            except ValueError:
                raise ValueError("Problem parsing: " + col[0] + ' ' + JString)
            myLevDens.add_spin_dependent_level_density(J=J, Pi=Pi, data=list(zip(UList, list(map(float, col[1:])))))
    return myLevDens


def readHFBMLevelDensityCorrections(f, verbose=False):
    """
    File is essentially undocumented, so going by what I see in the EMPIRE lev-dens.f file
    """
    result = {}
    for line in f.split('\n'):
        if line.strip() == '':
            continue
        sline = line.split()
        Z = int(sline[0])
        A = int(sline[1])
        x = sline[2]  # ?
        y = sline[3]  # ?
        aTildeCorrection = float(sline[4])
        pairingShiftCorrection = float(sline[5])
        name = sline[6]  # wrong format, so ignore it
        result[(Z, A)] = (aTildeCorrection, pairingShiftCorrection)
    return result
