from pqu import PQU
from brownies.BNL.inter.spectra import *
from brownies.BNL.utilities.XYs import *
from brownies.BNL.utilities.units import *
from brownies.BNL.utilities.gnds import *

# -------------------------------------------------------------------------------
# spectrum "macros"
# -------------------------------------------------------------------------------

DEFAULTONEOVEREGRID = [1e-5,
                       0.99999 * float(CADMIUMCUTTOFF.value)] + \
                      list(equal_lethargy_bins(8000,
                             domainMin=float(CADMIUMCUTTOFF.value),
                             domainMax=1.5e8))


def DEFAULTONEOVEREFUNCTION(E, *args):
    if E < float(CADMIUMCUTTOFF.value):
        return 0.0
    return 1.0 / E


DEFAULTONEOVEREXYs = function_to_XYs(DEFAULTONEOVEREFUNCTION, [], Egrid=DEFAULTONEOVEREGRID)


def CF252SPECTRUMAVEFUNCTION(E, *fpars): return math.exp(-0.789 * E / 1e6) * math.sinh(math.sqrt(0.447 * E / 1e6))


CF252SPECTRUMAVE = function_to_XYs(CF252SPECTRUMAVEFUNCTION, [])


def getProjectileAndTargetMasses(_xs):
    targetID = _xs.rootAncestor.target
    if targetID in _xs.rootAncestor.PoPs.aliases:
        targetID = _xs.rootAncestor.PoPs[targetID].pid
    targ = _xs.rootAncestor.PoPs[targetID]
    proj = _xs.rootAncestor.PoPs[_xs.rootAncestor.projectile]
    try:
        m2 = targ.mass[0].float('amu')
    except IndexError:
        m2 = targ.getMass('amu')
    m1 = proj.getMass('amu')
    return m1, m2


def convolve(funcA, funcB, useCovariance=False, covariance=None, normalize=True):
    if funcA.domainUnit != funcB.domainUnit:
       raise ValueError("Units differ")

    if funcA.domainMin != funcB.domainMin:
        if funcA.domainMin < funcB.domainMin:
            raise ValueError( 'funcA has a domainMin lower than funcB, I cannot extrapolate')
        funcB.domainSlice(domainMin=funcA.domainMin)

    if funcA.domainMax != funcB.domainMax:
        if funcA.domainMax < funcB.domainMax:
            funcB.domainSlice(domainMax=funcA.domainMax)
        else:
            funcB.setValue(1.000001*funcB.domainMax, 0.0)
            funcB.setValue(funcA.domainMax, 0.0)

    return funcA.integrateTwoFunctionsWithUncertainty(funcB,
                                                      useCovariance=useCovariance,
                                                      covariance=covariance,
                                                      normalize=normalize)


# -------------------------------------------------------------------------------
# Main functions
# -------------------------------------------------------------------------------
def computeRI(xs, Ecut=None, domainMax=None, useCovariance=True, covariance=None):
    """
    Compute the resonance integral (RI):

    .. math::
        RI = \int_{Ec}^{\infty} \sigma(E) dE/E

    Where the lower cut-off is usually taken to be the Cadmium cut-off energy of
    Ecut = 0.5 eV (see S. Mughabghab, Atlas of Neutron Resonances).  If the covariance
    is provided, the uncertainty on the resonance integral will be computed.

    :param xs: the cross section
    :param PQU Ecut: lower bound of integration.
        Usually taken to be the Cadmium cut-off energy (default: '0.5 eV')
    :param PQU domainMax: an alternative upper energy integration cut-off used for cross checking against INTER mainly
    :param useCovariance: use this to override covarance usage
    :type useCovariance: bool
    :param covariance: covariance to use when computing uncertainty on the spectral average.
        If None (default: None), no uncertainty is computed.
    :type covariance: covariance instance or None
    :rtype: PQU

    """
    check_is_cross_section(xs)
    if domainMax is None:
        domainMax = PQU.PQU(min(xs.domainMax, DEFAULTONEOVEREXYs.domainMax), xs.domainUnit)

    if Ecut is None:
        return xs.integrateTwoFunctionsWithUncertainty(DEFAULTONEOVEREXYs,
                                                       domainMax=domainMax,
                                                       useCovariance=useCovariance,
                                                       covariance=covariance,
                                                       normalize=False)

    Ecut = PQU.PQU(Ecut).getValueAs(xs.domainUnit)

    # Construct an XYs instance that spans the same domain as self.
    # The y values are all 1/E and the x values are all E.
    def oneOverEFunc(E, *args):
        if E < Ecut:
            return 0.0
        return 1.0 / E

    Egrid = [1e-5, 0.99999 * Ecut] + list(equal_lethargy_bins(5000, domainMin=Ecut))
    oneOverE = function_to_XYs(oneOverEFunc, [], Egrid=Egrid)
    return xs.integrateTwoFunctionsWithUncertainty(oneOverE,
                                                   domainMax=domainMax,
                                                   useCovariance=useCovariance,
                                                   covariance=covariance,
                                                   normalize=False)


def computeMACS(xs, T, a=None, useCovariance=True, covariance=None, normalize=False):
    """
    Compute the Maxwellian average cross section (MACS):

    .. math::
        MACS(kT) = {\\frac{2}{\sqrt{\pi}}} {\\frac{a^2}{(kT)^2}} \int_0^{\infty}dE_n^{lab} \sigma(E_n^{lab})*
                   E_n^{lab}*\exp{(-a E_n^{lab}/kT)}

    Here :math:`E_n^{lab}` is the incident neutron energy in the lab frame and :math:`a=m_2/(m_1+m_2)`.
    If the covariance is provided, the uncertainty on the MACS will be computed.

    :param xs: the cross section
    :param PQU T: Temperature (as a physical quantity of the Maxwellian with which to average
    :param float a: mass ratio m2/(m1+m2) where m1 is the projectile mass and m2 is the target mass.
        If set to None, will try to determine from evaluation itself.
    :param useCovariance: use this to override covarance usage
    :type useCovariance: bool
    :param covariance: covariance to use when computing uncertainty on the spectral average.
        If None (default: None), no uncertainty is computed.
    :type covariance: covariance instance or None
    :param normalize: flag to normalize
    :rtype: PQU
    """
    check_is_cross_section(xs)

    kT = convert_units_to_energy(T)

    if a is None:
        # targetID = xs.rootAncestor.target
        # if targetID in xs.rootAncestor.PoPs.aliases:
        #     targetID = xs.rootAncestor.PoPs[targetID].pid
        # targ = xs.rootAncestor.PoPs[targetID]
        # proj = xs.rootAncestor.PoPs[xs.rootAncestor.projectile]
        # try:
        #     m2 = targ.mass[0].float('amu')
        # except IndexError:
        #     m2 = targ.getMass('amu')
        # m1 = proj.getMass('amu')
        m1, m2 = getProjectileAndTargetMasses(xs)
        a = float(m2 / (m1 + m2))
    elif isinstance(a, PQU.PQU):
        a = float(a.getValueAs(""))

    # Construct an XYs instance that spans the same domain as self.
    # The x values are the E values and the y values are a Maxwellian evaluated at x.
    a_over_kT = a / float(kT.getValueAs(xs.domainUnit))
    norm = 2.0 * a_over_kT * a_over_kT / math.sqrt(math.pi)

    def maxwellianFunc(E, *args):
        return norm * E * math.exp(-a_over_kT * E)

    maxwellian = function_to_XYs(maxwellianFunc, [])

    return convolve(xs, maxwellian, useCovariance=useCovariance, covariance=covariance, normalize=normalize)


def computeCf252SpectrumAveAnalytic(xs, useCovariance=True, covariance=None):
    """
    Compute the Cf252 spontaneous fission neutron spectrum average:

    .. math::
        ave = \exp( -0.789 * E ) * \sinh( \sqrt{ 0.447 * E } )

    (where did this come from?  Found it in the sigma QA system I think. DAB 12/1/2012)

     This is what INTER uses:
         Fiss.spect. average   : Sig(Fiss) = Intg[E1:E2] Sig(E) Phi_f(E) dE / Intg[E1:E2] Phi_f(E) dE
         Fission spectrum      : Phi_f(E)  = sqrt(E/Ef)/Ef exp(-E/E0)
         Spectrum Temperature  : Ef =  1.35000E+06 (eV)
         Integration Limits    : E1 =  1.00000E+03 (eV)  E2 =  2.00000E+07 (eV)
         Integral of Spectrum       =  1.00000E+00

    If the covariance is provided, the uncertainty on the average will be computed.

    :param xs: the cross section
    :param useCovariance: use this to override covarance usage
    :type useCovariance: bool
    :param covariance: covariance to use when computing uncertainty on the spectral average.
        If None (default: None), no uncertainty is computed.
    :type covariance: covariance instance or None
    :rtype: PQU

    """
    check_is_cross_section(xs)

    return convolve(xs, CF252SPECTRUMAVE, useCovariance=useCovariance, covariance=covariance, normalize=True)


def computeAverageOverWindow(xs, windowCenter, windowWidth, useCovariance=True, covariance=None):
    lo = max(PQU.PQU(xs.domainMin, unit=xs.domainUnit), (windowCenter - windowWidth / 2.0))
    hi = min(PQU.PQU(xs.domainMax, unit=xs.domainUnit), (windowCenter + windowWidth / 2.0))
    check_is_cross_section(xs)
    return xs.toPointwise_withLinearXYs().integrate(domainMin=lo, domainMax=hi) / windowWidth.inUnitsOf(xs.domainUnit)


def computeFourteenMeVPoint(xs, E14='14.2 MeV', useCovariance=True, covariance=None):
    """
    Compute the value of the cross section at 14.2 MeV.
    If the covariance is provided, the uncertainty on the 14.2 MeV point will be computed.

    :param xs: reference to the cross section
    :param E14: the 14 MeV point to use (in case you want to override the default of 14.2 MeV)
    :param useCovariance: use this to override covariance usage
    :type useCovariance: bool
    :param covariance: covariance to use when computing uncertainty on the spectral average.
        If None (default: None), no uncertainty is computed.
    :type covariance: covariance instance or None
    :rtype: PQU
    """
    return computeValueAtAPoint(xs, E14, useCovariance=useCovariance, covariance=covariance)


def computeValueAtAPoint(xs, point, useCovariance=True, covariance=None):
    """
    Compute the value of the cross section at a point.
    If the covariance is provided, the uncertainty will be computed.

    :param xs: reference to the cross section
    :param point: the point to evaluate this at, something like '14.2 MeV'
    :param useCovariance: use this to override covariance usage
    :type useCovariance: bool
    :param covariance: covariance to use when computing uncertainty on the spectral average.
        If None (default: None), no uncertainty is computed.
    :type covariance: covariance instance or None
    :rtype: PQU
    """
    check_is_cross_section(xs)

    return xs.evaluateWithUncertainty(PQU.PQU(point),
                                      useCovariance=useCovariance,
                                      covariance=covariance)


def computeRoomTempCS(xs, T=None, useCovariance=True, covariance=None):
    """
    The cross section exactly at room temperature E=0.0235339347844 eV or v=2200 m/s

    If the covariance is provided, the uncertainty on the 0.0235339347844 eV point will be computed.

    RoomTemp = 293.15 degK (== 20. degC)

    What temperature does INTER use?  2.53000E-02 eV.   What about Boris? dunno

    :param xs: the cross section
    :param T: the temperature
    :param useCovariance: use this to override covarance usage
    :type useCovariance: bool
    :param covariance: covariance to use when computing uncertainty on the spectral average.
        If None (default: None), will try to get from evaluation itself.
    :type covariance: covariance instance or None
    :rtype: PQU
    """
    check_is_cross_section(xs)

    if T is None:
        ERT = convert_units_to_energy(ROOMTEMPERATURE)
    else:
        ERT = convert_units_to_energy(T)
    return xs.evaluateWithUncertainty(ERT, useCovariance=useCovariance, covariance=covariance)


def computeWestcottFactor(xs, a=None, T=None, useCovariance=True, covariance=None, normalize=True):
    """
    Wescott G-factor is the ratio of Maxwellian averaged cross section
    (at room temparature) and the room temperature cross section.
    Should be pretty close to 1 if cross section goes like 1/v.

    INTER also only integrates on interval (1e-5,10) eV

    :param xs: the cross section
    :param a: mass ratio m2/(m1+m2) where m2 is target mass and m1 is projectile mass (if None, will try to get from
              evaluation itself)
    :param T: the temperature
    :param useCovariance: use this to override covariance usage
    :type useCovariance: bool
    :param covariance: covariance to use when computing uncertainty on the spectral average.
        If None (default: None), will try to get from evaluation itself.
    :type covariance: covariance instance or None
    :param normalize: FIXME
    :rtype: PQU
    """
    check_is_cross_section(xs)

    if T is None:
        ERT = convert_units_to_energy(ROOMTEMPERATURE)
    else:
        ERT = convert_units_to_energy(T)
    macs = computeMACS(xs, a=a, T=ERT, useCovariance=useCovariance, covariance=covariance, normalize=normalize)
    rt = computeRoomTempCS(xs, T=ERT, useCovariance=useCovariance, covariance=covariance)
    if macs.value == 0.0:
        return PQU.PQU(0.0, '')
    westcott = (2.0 / math.sqrt(math.pi)) * macs / rt
    return westcott


def computeAstrophysicalReactionRate(xs, T, useCovariance=True, covariance=None):
    """
    The astrophysical reaction rate R is defined as :math:`R = N_A <\sigma v>`, where
    :math:`N_A` is Avagadro's number.  It is typically reported in units of [cm^3/mole s].

    :param xs: the cross section
    :param PQU T: temperature as a physical quantity
    :param useCovariance: use this to override covarance usage
    :type useCovariance: bool
    :param covariance: covariance to use when computing uncertainty on the spectral average.
        If None (default: None), no uncertainty is computed.
    :type covariance: covariance instance or None
    :rtype: PQU

    .. warning:: Not implemented yet
    """
    check_is_cross_section(xs)

    if T is None:
        kT = convert_units_to_energy(ROOMTEMPERATURE)
    else:
        kT = convert_units_to_energy(T)

    # m2 = xs.rootAncestor.PoPs[xs.rootAncestor.target].getMass('amu')
    # m1 = xs.rootAncestor.PoPs[xs.rootAncestor.projectile].getMass('amu')
    m1, m2 = getProjectileAndTargetMasses(xs)
    mu = PQU.PQU(m1 * m2 / (m1 + m2), 'amu')

    vT = (2.0 * kT / mu).sqrt()
    return (AVAGADROSNUMBER * computeMACS(xs, kT,
                                          useCovariance=useCovariance,
                                          covariance=covariance) * vT).inUnitsOf("cm**3/mol/s")


def computeALPHA(capxs, fisxs, E=convert_units_to_energy(ROOMTEMPERATURE), useCovariance=True):
    return capxs.evaluateWithUncertainty(E, useCovariance=useCovariance) / \
           fisxs.evaluateWithUncertainty(E, useCovariance=useCovariance)


def computeETA(capxs, fisxs, nubar, E=convert_units_to_energy(ROOMTEMPERATURE), useCovariance=True):
    return nubar.evaluate(E.getValueAs('eV')) / (1.0 + computeALPHA(capxs, fisxs, E=E, useCovariance=useCovariance))


def computeScatteringRadius(rs, useCovariance=True, covariance=None, verbose=False):
    """
    Compute R', the scattering radius

    We take the E=0 point (or at least the lowest available energy point in the elastic cross section).

    with this,
        ..math::
            \sigma_{el}=4\pi(R')^2

    :param rs: the reactionSuite
    :param useCovariance: FIXME
    :param covariance: FIXME
    :param verbose: FIXME
    :return:
    """
    # get RRR parameter averages
    if hasattr(rs.resonances, 'resolved') and rs.resonances.resolved is not None:
        import fudge.processing.resonances.reconstructResonances as RRReconstruct
        resCls = RRReconstruct.getResonanceReconstructionClass(rs.resonances.resolved.evaluated)
        rrr = resCls(rs.resonances.resolved.evaluated, enableAngDists=False, verbose=verbose)
        rrr.setResonanceParametersByChannel()
        R = rrr.getScatteringLength()
        return PQU.PQU(10. * R, 'fm')
    elif hasattr(rs.resonances, 'scatteringRadius') and rs.resonances.scatteringRadius is not None:
        return PQU.PQU(rs.resonances.scatteringRadius.getValueAs('fm', energy_grid=[1e-5]), 'fm')
    else:
        return PQU.PQU(0.0, 'fm')


def computeCf252SpectrumAve(xs, useCovariance=True, covariance=None):
    check_is_cross_section(xs)
    docs = """
                   *****   MT  =  18   ****

             NEUTRON SPECTRUM OF SPONTANEOUS FISSION

        DOCUMENTATION :

        1. W. MANNHART, IAEA-TECDOC-410 (1987) 158
        2. W. MANNHART, INDC(NDS)-220/L (1989) 305

                     ****   MF =  5   ****

        THE POINTWISE SPECTRUM IS IN UNITS OF [1/EV] . THE INTERPOLATION
        IS LINEAR - LINEAR. THE ORIGINAL DATA HAVE BEEN CONDENSED AND
        MODIFIED IN THE INTERPOLATION RULES.

                     ****   MF = 35   ****

        THE COVARIANCE MATRIX IS EVALUATED BETWEEN 15 KEV AND 20 MEV.
        BELOW 15 KEV THE UNCERTAINTY IS ESTIMATED.

        THE MATRIX IS GIVEN IN ABSOLUTE UNITS.
        THE DERIVATION OF THE RELATIVE COVARIANCE MATRIX REQUIRES THE
        NUMERICAL CALCULATION OF GROUP-INTEGRALS OVER THE SPECTRUM.
        FOR CONVENIENCE THESE INTEGRALS AND THE RELATIVE STANDARD
        DEVIATIONS (RSD) OF THE SPECTRUM ARE LISTED BELOW:

        GROUP ENERGY RANGE   SPECTRUM    RSD
          NO.   ( IN MEV )    INTEGRAL     %

          1   0.000   0.015  7.71080E-4  30.0
          2   0.015   0.035  1.96699E-3  10.4
          3   0.035   0.055  2.63500E-3   4.7
          4   0.055   0.075  3.12875E-3   4.3
          5   0.075   0.095  3.53267E-3   4.8
          6   0.095   0.115  3.87271E-3   3.4
          7   0.115   0.135  4.17008E-3   3.4
          8   0.135   0.165  6.73767E-3   2.7
          9   0.165   0.195  7.23000E-3   2.2
         10   0.195   0.225  7.65712E-3   2.3
         11   0.225   0.255  8.02452E-3   1.9
         12   0.255   0.305  1.40620E-2   1.8
         13   0.305   0.355  1.47536E-2   1.7
         14   0.355   0.405  1.53208E-2   1.7
         15   0.405   0.455  1.57517E-2   1.7
         16   0.455   0.505  1.61214E-2   1.7
         17   0.505   0.555  1.63615E-2   1.8
         18   0.555   0.605  1.65634E-2   1.8
         19   0.605   0.655  1.66960E-2   1.6
         20   0.655   0.705  1.67805E-2   1.8
         21   0.705   0.755  1.68040E-2   1.9
         22   0.755   0.805  1.67857E-2   1.8
         23   0.805   0.855  1.67669E-2   1.8
         24   0.855   0.905  1.66970E-2   1.9
         25   0.905   0.955  1.65920E-2   1.7
         26   0.955   1.050  3.12024E-2   1.6
         27   1.050   1.150  3.22050E-2   1.2
         28   1.150   1.250  3.14902E-2   1.2
         29   1.250   1.350  3.05967E-2   1.2
         30   1.350   1.450  2.96733E-2   1.2
         31   1.450   1.550  2.87347E-2   1.2
         32   1.550   1.650  2.77040E-2   1.2
         33   1.650   1.750  2.66580E-2   1.3
         34   1.750   1.850  2.56120E-2   1.2
         35   1.850   1.950  2.45660E-2   1.2
         36   1.950   2.150  4.60401E-2   1.2
         37   2.150   2.350  4.20150E-2   1.1
         38   2.350   2.550  3.81019E-2   1.1
         39   2.550   2.750  3.44613E-2   1.2
         40   2.750   2.950  3.10800E-2   1.2
         41   2.950   3.250  4.07617E-2   1.2
         42   3.250   3.550  3.44418E-2   1.2
         43   3.550   3.850  2.89682E-2   1.2
         44   3.850   4.150  2.42333E-2   1.2
         45   4.150   4.450  2.02003E-2   1.3
         46   4.450   4.750  1.67786E-2   1.4
         47   4.750   5.050  1.39144E-2   1.4
         48   5.050   5.550  1.80559E-2   1.5
         49   5.550   6.050  1.31172E-2   1.5
         50   6.050   6.550  9.48581E-3   1.6
         51   6.550   7.050  6.83338E-3   1.6
         52   7.050   7.550  4.90622E-3   1.8
         53   7.550   8.050  3.51801E-3   1.9
         54   8.050   8.550  2.52038E-3   2.1
         55   8.550   9.050  1.80664E-3   2.2
         56   9.050   9.550  1.29426E-3   2.2
         57   9.550  10.050  9.26906E-4   2.5
         58  10.050  10.550  6.62406E-4   2.7
         59  10.550  11.050  4.73684E-4   2.9
         60  11.050  11.550  3.38507E-4   2.9
         61  11.550  12.050  2.41879E-4   3.4
         62  12.050  12.550  1.72563E-4   3.6
         63  12.550  13.050  1.23026E-4   5.0
         64  13.050  13.550  8.74750E-5   4.5
         65  13.550  14.050  6.21894E-5   6.3
         66  14.050  14.600  4.77954E-5   9.6
         67  14.600  15.900  6.26595E-5  12.2
         68  15.900  16.900  2.15184E-5  14.5
         69  16.900  17.900  1.08596E-5  19.0
         70  17.900  19.100  6.11592E-6  31.8
         71  19.100  20.000  2.17986E-6  75.7
    """
    energies = [1e-5]
    fluxes = []
    for line in docs.split('\n')[30:-1]:
        sline = line.split()
        energies.append(float(sline[2]) * 1e6)
        fluxes.append(float(sline[3]) / (float(sline[2]) - float(sline[1])) / 1e6)
    spectrum = grouped_values_to_XYs(energies, fluxes, domainUnit='eV', rangeUnit='1/eV')
    return convolve(xs, spectrum, useCovariance=useCovariance, covariance=covariance, normalize=True)


def computeGodivaSpectrumAve(xs, useCovariance=True, covariance=None):
    """
    Godiva spectrum average

    FIXME: THIS IS USELESS UNLESS WE DEFINE WHAT WE MEAN BY "SPECTRUM AVERAGE IN GODIVA" BETTER!! ARE WE
    IN THE MIDDLE OF THE ASSEMBLY?  ON THE OUTSIDE?  HOW LONG DO WE IRRADIATE?

    :param xs: cross section to average
    :param useCovariance: a flag
    :param covariance: a reference to the covariance of xs, if it cannot be determined otherwise
    :return:
    """
    check_is_cross_section(xs)
    #    spectrum=get_KENO_flux(os.sep.join([INTERDIR,'spectra','HMF001.001']))
    spectrum = get_spectrum("Godiva")
    return convolve(xs, spectrum, useCovariance=useCovariance, covariance=covariance, normalize=True)


def computeJezebelSpectrumAve(xs, useCovariance=True, covariance=None):
    """
    Jezebel spectrum average

    FIXME: THIS IS USELESS UNLESS WE DEFINE WHAT WE MEAN BY "SPECTRUM AVERAGE IN JEZEBEL" BETTER!! ARE WE
    IN THE MIDDLE OF THE ASSEMBLY?  ON THE OUTSIDE?  HOW LONG DO WE IRRADIATE?

    :param xs: cross section to average
    :param useCovariance: a flag
    :param covariance: a reference to the covariance of xs, if it cannot be determined otherwise
    :return:
    """
    check_is_cross_section(xs)
    #    spectrum=get_KENO_flux(os.sep.join([INTERDIR,'spectra','PMF001.001']))
    spectrum = get_spectrum("Jezebel")
    return convolve(xs, spectrum, useCovariance=useCovariance, covariance=covariance, normalize=True)


def computeBigTenSpectrumAve(xs, useCovariance=True, covariance=None):
    """
    BigTen spectrum average

    FIXME: THIS IS USELESS UNLESS WE DEFINE WHAT WE MEAN BY "SPECTRUM AVERAGE IN BIGTEN" BETTER!! ARE WE
    IN THE MIDDLE OF THE ASSEMBLY?  ON THE OUTSIDE?  HOW LONG DO WE IRRADIATE?

    :param xs: cross section to average
    :param useCovariance: a flag
    :param covariance: a reference to the covariance of xs, if it cannot be determined otherwise
    :return:
    """
    check_is_cross_section(xs)
    #    spectrum=get_KENO_flux(os.sep.join([INTERDIR,'spectra','IMF007.001']))
    spectrum = get_spectrum("BigTen")
    return convolve(xs, spectrum, useCovariance=useCovariance, covariance=covariance, normalize=True)


def computeLookupSpectrumAve(xs, key, useCovariance=True, covariance=None):
    """
    Compute pectrum average of one of the spectra stored in the spectra module

    :param xs: cross section to average
    :param key: the lookup key of the spectra
    :param useCovariance: a flag
    :param covariance: a reference to the covariance of xs, if it cannot be determined otherwise
    :return:
    """
    check_is_cross_section(xs)
    spectrum = get_spectrum(key)
    if not isinstance(spectrum, XYs.XYs1d):
        raise TypeError("Not instance of XYs.XYs1d")
    return convolve(xs, spectrum, useCovariance=useCovariance, covariance=covariance, normalize=True)


def computeUserDefinedSpectrumAve(xs, path, key, useCovariance=True, covariance=None):
    """
    User defined spectrum average

    :param xs: cross section to average
    :param path: path to the user's spectra file
    :param key: key of the spectra in the user's spectra file
    :param useCovariance: a flag
    :param covariance: a reference to the covariance of xs, if it cannot be determined otherwise
    :return:
    """
    check_is_cross_section(xs)
    spectrum = get_user_defined_flux(path, key)
    if not isinstance(spectrum, XYs.XYs1d):
        raise TypeError("Not instance of XYs.XYs1d")
    return convolve(xs, spectrum, useCovariance=useCovariance, covariance=covariance, normalize=True)
