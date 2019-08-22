from fudge.gnds.reactionData.crossSection import lowerEps
from pqu import PQU
from spectra import *

# -------------------------------------------------------------------------------
# spectrum "macros"
# -------------------------------------------------------------------------------

DEFAULTONEOVEREGRID = [1e-5,0.99999*float(CADMIUMCUTTOFF.value)]+list(equal_lethargy_bins(2000,domainMin=float(CADMIUMCUTTOFF.value)))
def DEFAULTONEOVEREFUNCTION(E,*args):
    if E < float(CADMIUMCUTTOFF.value): return 0.0
    return 1.0/E
DEFAULTONEOVEREXYs = function_to_XYs( DEFAULTONEOVEREFUNCTION, [], Egrid=DEFAULTONEOVEREGRID )


def CF252SPECTRUMAVEFUNCTION(E,*fpars): return math.exp( -0.789 * E/1e6 ) * math.sinh( math.sqrt( 0.447 * E/1e6 ) )
CF252SPECTRUMAVE = function_to_XYs( CF252SPECTRUMAVEFUNCTION, [] )



# -------------------------------------------------------------------------------
# Main functions
# -------------------------------------------------------------------------------

def computeRI( xs, Ecut=None, domainMax=None, useCovariance=True, covariance=None ):
    '''
    Compute the resonance integral (RI):
    
    .. math::
        RI = \int_{Ec}^{\infty} \sigma(E) dE/E
        
    Where the lower cut-off is usually taken to be the Cadmium cut-off energy of 
    Ecut = 0.5 eV (see S. Mughabghab, Atlas of Neutron Resonances).  If the covariance
    is provided, the uncertainty on the resonance integral will be computed.
    
    :param PQU Ecut: lower bound of integration.
        Usually taken to be the Cadmium cut-off energy (default: '0.5 eV')
    :param PQU domainMax: an alternative upper energy integration cut-off used for cross checking against INTER mainly
    :param useCovariance: use this to override covarance usage
    :type useCovariance: bool
    :param covariance: covariance to use when computing uncertainty on the spectral average.
        If None (default: None), no uncertainty is computed.
    :type covariance: covariance instance or None
    :rtype: PQU
    
    '''
    check_is_cross_section( xs )
    if Ecut is None:
        return xs.integrateTwoFunctionsWithUncertainty( DEFAULTONEOVEREXYs, domainMax=domainMax, useCovariance=useCovariance, covariance=covariance, normalize=False )

    Ecut = PQU.PQU(Ecut).getValueAs( xs.domainUnit )
    # Construct an XYs instance that spans the same domain as self.
    # The y values are all 1/E and the x values are all E.
    def oneOverEFunc(E,*args): 
        if E < Ecut: return 0.0
        return 1.0/E
    Egrid=[1e-5,0.99999*Ecut]+list(equal_lethargy_bins(5000,domainMin=Ecut))
    oneOverE = function_to_XYs( oneOverEFunc, [], Egrid=Egrid )
    return xs.integrateTwoFunctionsWithUncertainty( oneOverE, domainMax=domainMax, useCovariance=useCovariance, covariance=covariance, normalize=False )


def computeMACS( xs, T, a=None, useCovariance=True, covariance=None, normalize=False ):
    '''
    Compute the Maxwellian average cross section (MACS):

    .. math::
        MACS(kT) = {\\frac{2}{\sqrt{\pi}}} {\\frac{a^2}{(kT)^2}} \int_0^{\infty}dE_n^{lab} \sigma(E_n^{lab})*E_n^{lab}*\exp{(-a E_n^{lab}/kT)}

    Here :math:`E_n^{lab}` is the incident neutron energy in the lab frame and :math:`a=m_2/(m_1+m_2)`.
    If the covariance is provided, the uncertainty on the MACS will be computed.

    :param PQU T: Temperature (as a physical quantity of the Maxwellian with which to average
    :param float a: mass ratio m2/(m1+m2) where m1 is the projectile mass and m2 is the target mass.  
        If set to None, will try to determine from evaluation itself.
    :param useCovariance: use this to override covarance usage
    :type useCovariance: bool
    :param covariance: covariance to use when computing uncertainty on the spectral average.
        If None (default: None), no uncertainty is computed.
    :type covariance: covariance instance or None
    :rtype: PQU
    '''
    check_is_cross_section( xs )
    
    kT=convert_units_to_energy(T)

    if a is None:
        m2 = xs.getRootAncestor().PoPs[xs.getRootAncestor().target].getMass('amu')
        m1 = xs.getRootAncestor().PoPs[xs.getRootAncestor().projectile].getMass('amu')
        a = float(m2/(m1+m2))
    elif isinstance(a,PQU.PQU): a=float(a.getValueAs(""))
    
    # Construct an XYs instance that spans the same domain as self.
    # The x values are the E values and the y values are a Maxwellian evaluated at x.
    a_over_kT = a/float(kT.getValueAs( xs.domainUnit ) )
    norm = 2.0*a_over_kT*a_over_kT/math.sqrt(math.pi)
    def maxwellianFunc(E,*args): return norm*E*math.exp(-a_over_kT*E)
    maxwellian = function_to_XYs( maxwellianFunc, [] )
    return xs.integrateTwoFunctionsWithUncertainty( maxwellian, useCovariance=useCovariance, covariance=covariance, normalize=normalize )


def computeCf252SpectrumAveAnalytic( xs, useCovariance=True, covariance=None ):
    '''
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
    
    :param useCovariance: use this to override covarance usage
    :type useCovariance: bool
    :param covariance: covariance to use when computing uncertainty on the spectral average.
        If None (default: None), no uncertainty is computed.
    :type covariance: covariance instance or None
    :rtype: PQU

    '''
    check_is_cross_section( xs )
    return xs.integrateTwoFunctionsWithUncertainty( CF252SPECTRUMAVE, useCovariance=useCovariance, covariance=covariance, normalize=True )


def computeAverageOverWindow( xs, windowCenter, windowWidth, useCovariance=True, covariance=None ):
    from xData.XYs import XYs1d
    lo=max(PQU.PQU(xs.domainMin,unit=xs.domainUnit),(windowCenter-windowWidth/2.0))
    hi=min(PQU.PQU(xs.domainMax,unit=xs.domainUnit),(windowCenter+windowWidth/2.0))
    check_is_cross_section( xs )
    return  xs.toPointwise_withLinearXYs().integrate(domainMin=lo,domainMax=hi)/windowWidth.inUnitsOf(xs.domainUnit)


def computeFourteenMeVPoint( xs, E14='14.2 MeV', useCovariance=True, covariance=None ):
    '''
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
    '''
    return computeValueAtAPoint(xs, E14, useCovariance=useCovariance,covariance=covariance)


def computeValueAtAPoint( xs, point, useCovariance=True, covariance=None ):
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
    check_is_cross_section( xs )
    return xs.evaluateWithUncertainty(PQU.PQU(point),useCovariance=useCovariance,covariance=covariance)


def computeRoomTempCS( xs, T=None, useCovariance=True, covariance=None ):
    '''
    The cross section exactly at room temperature E=0.0235339347844 eV or v=2200 m/s

    If the covariance is provided, the uncertainty on the 0.0235339347844 eV point will be computed.
    
    RoomTemp = 293.15 degK (== 20. degC)

    What temperature does INTER use?  2.53000E-02 eV.   What about Boris? dunno

    :param useCovariance: use this to override covarance usage
    :type useCovariance: bool
    :param covariance: covariance to use when computing uncertainty on the spectral average.
        If None (default: None), will try to get from evaluation itself.
    :type covariance: covariance instance or None
    :rtype: PQU
    '''
    check_is_cross_section( xs )
    if T is None: ERT=convert_units_to_energy(ROOMTEMPERATURE)
    else:         ERT=convert_units_to_energy(T)
    return xs.evaluateWithUncertainty(ERT,useCovariance=useCovariance,covariance=covariance)


def computeWestcottFactor( xs, a=None, T=None, useCovariance=True, covariance=None, normalize=True ):
    '''
    Wescott G-factor is the ratio of Maxwellian averaged cross section 
    (at room temparature) and the room temperature cross section.  
    Should be pretty close to 1 if cross section goes like 1/v.
    
    INTER also only integrates on interval (1e-5,10) eV

    :param a: mass ratio m2/(m1+m2) where m2 is target mass and m1 is projectile mass (if None, will try to get from evaluation itself)
    :param useCovariance: use this to override covarance usage
    :type useCovariance: bool
    :param covariance: covariance to use when computing uncertainty on the spectral average.
        If None (default: None), will try to get from evaluation itself.
    :type covariance: covariance instance or None
    :rtype: PQU
    '''
    check_is_cross_section( xs )
    if T is None: ERT=convert_units_to_energy(ROOMTEMPERATURE)
    else:         ERT=convert_units_to_energy(T)
    macs = computeMACS(xs, a=a, T=ERT, useCovariance=useCovariance, covariance=covariance, normalize=normalize)
    rt = computeRoomTempCS(xs, T=ERT, useCovariance=useCovariance, covariance=covariance)
    if macs.value == 0.0: return PQU.PQU(0.0,'')
    westcott = (2.0/math.sqrt(math.pi))*macs/rt
    return westcott


def computeAstrophysicalReactionRate( xs, T, useCovariance=True, covariance=None ):
    '''
    The astrophysical reaction rate R is defined as :math:`R = N_A <\sigma v>`, where 
    :math:`N_A` is Avagadro's number.  It is typically reported in units of [cm^3/mole s].

    :param PQU T: temperature as a physical quantity
    :param useCovariance: use this to override covarance usage
    :type useCovariance: bool
    :param covariance: covariance to use when computing uncertainty on the spectral average.
        If None (default: None), no uncertainty is computed.
    :type covariance: covariance instance or None
    :rtype: PQU
    
    .. warning:: Not implemented yet
    '''
    check_is_cross_section( xs )

    if T is None: kT=convert_units_to_energy(ROOMTEMPERATURE)
    else:         kT=convert_units_to_energy(T)

    m2 = xs.getRootAncestor().PoPs[xs.getRootAncestor().target].getMass('amu')
    m1 = xs.getRootAncestor().PoPs[xs.getRootAncestor().projectile].getMass('amu')
    mu = PQU.PQU(m1*m2/(m1+m2),'amu')

    vT=(2.0*kT/mu).sqrt()
    return (AVAGADROSNUMBER*computeMACS(xs,kT,useCovariance=useCovariance,covariance=covariance)*vT).inUnitsOf( "cm**3/mol/s" )


def computeALPHA(capxs, fisxs, E=convert_units_to_energy(ROOMTEMPERATURE), useCovariance=True):
    return capxs.evaluateWithUncertainty(E, useCovariance=useCovariance)/\
            fisxs.evaluateWithUncertainty(E, useCovariance=useCovariance)


def computeETA(capxs, fisxs, nubar, E=convert_units_to_energy(ROOMTEMPERATURE), useCovariance=True):
    return nubar.evaluate(E.getValueAs('eV'))/(1.0+computeALPHA(capxs, fisxs, E=E, useCovariance=useCovariance))


def computeScatteringRadius( rs, useCovariance=True, covariance=None, verbose=False ):
    """
    Compute R', the scattering radius

    We take the E=0 point (or at least the lowest available energy point in the elastic cross section).

    with this,
        ..math::
            \sigma_{el}=4\pi(R')^2

    :param elxs: the elastic scattering cross section
    :param useCovariance:
    :param covariance:
    :return:
    """
    # get RRR parameter averages
    if (hasattr(rs.resonances,'resolved') and rs.resonances.resolved is not None):
        import fudge.processing.resonances.reconstructResonances as RRReconstruct
        resCls = RRReconstruct.getResonanceReconstructionClass(rs.resonances.resolved.evaluated)
        rrr = resCls( rs, None, enableAngDists=False, verbose=verbose )
        rrr.setResonanceParametersByChannel()
        R=rrr.getScatteringLength()
        return PQU.PQU(10.*R,'fm')
    elif (hasattr(rs.resonances,'scatteringRadius') and rs.resonances.scatteringRadius is not None):
        return PQU.PQU(rs.resonances.scatteringRadius.getValueAs( 'fm', energy_grid=[1e-5], L=0 ),'fm')
    else:
        return PQU.PQU(0.0,'fm')


def computeCf252SpectrumAve( xs, useCovariance=True, covariance=None ):
    check_is_cross_section( xs )
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
    fluxes=[]
    for line in docs.split('\n')[30:-1]:
        sline=line.split()
        energies.append(float(sline[2])*1e6)
        fluxes.append(float(sline[3])/(float(sline[2])-float(sline[1]))/1e6)
    spectrum=grouped_values_to_XYs( energies, fluxes, domainUnit='eV', rangeUnit='1/eV' )
    return xs.integrateTwoFunctionsWithUncertainty( spectrum, useCovariance=useCovariance, covariance=covariance, normalize=True )


def computeGodivaSpectrumAve( xs, useCovariance=True, covariance=None ):
    """
    Godiva spectrum average

    FIXME: THIS IS USELESS UNLESS WE DEFINE WHAT WE MEAN BY "SPECTRUM AVERAGE IN GODIVA" BETTER!! ARE WE
    IN THE MIDDLE OF THE ASSEMBLY?  ON THE OUTSIDE?  HOW LONG DO WE IRRADIATE?

    :param xs: cross section to average
    :param useCovariance: a flag
    :param covariance: a reference to the covariance of xs, if it cannot be determined otherwise
    :return:
    """
    check_is_cross_section( xs )
#    spectrum=get_KENO_flux(os.sep.join([INTERDIR,'spectra','HMF001.001']))
    spectrum = get_spectrum("Godiva")
    return xs.integrateTwoFunctionsWithUncertainty( spectrum, useCovariance=useCovariance, covariance=covariance, normalize=True )


def computeJezebelSpectrumAve( xs, useCovariance=True, covariance=None ):
    """
    Jezebel spectrum average

    FIXME: THIS IS USELESS UNLESS WE DEFINE WHAT WE MEAN BY "SPECTRUM AVERAGE IN JEZEBEL" BETTER!! ARE WE
    IN THE MIDDLE OF THE ASSEMBLY?  ON THE OUTSIDE?  HOW LONG DO WE IRRADIATE?

    :param xs: cross section to average
    :param useCovariance: a flag
    :param covariance: a reference to the covariance of xs, if it cannot be determined otherwise
    :return:
    """
    check_is_cross_section( xs )
#    spectrum=get_KENO_flux(os.sep.join([INTERDIR,'spectra','PMF001.001']))
    spectrum = get_spectrum("Jezebel")
    return xs.integrateTwoFunctionsWithUncertainty( spectrum, useCovariance=useCovariance, covariance=covariance, normalize=True )


def computeBigTenSpectrumAve( xs, useCovariance=True, covariance=None ):
    """
    BigTen spectrum average

    FIXME: THIS IS USELESS UNLESS WE DEFINE WHAT WE MEAN BY "SPECTRUM AVERAGE IN BIGTEN" BETTER!! ARE WE
    IN THE MIDDLE OF THE ASSEMBLY?  ON THE OUTSIDE?  HOW LONG DO WE IRRADIATE?

    :param xs: cross section to average
    :param useCovariance: a flag
    :param covariance: a reference to the covariance of xs, if it cannot be determined otherwise
    :return:
    """
    check_is_cross_section( xs )
#    spectrum=get_KENO_flux(os.sep.join([INTERDIR,'spectra','IMF007.001']))
    spectrum=get_spectrum("BigTen")
    return xs.integrateTwoFunctionsWithUncertainty( spectrum, useCovariance=useCovariance, covariance=covariance, normalize=True )


def computeLookupSpectrumAve( xs, key, verbose=True, useCovariance=True, covariance=None ):
    """
    Compute pectrum average of one of the spectra stored in the spectra module

    :param xs: cross section to average
    :param key: the lookup key of the spectra
    :param useCovariance: a flag
    :param covariance: a reference to the covariance of xs, if it cannot be determined otherwise
    :return:
    """
    check_is_cross_section( xs )
    spectrum=get_spectrum( key )
    if not isinstance( spectrum, XYs.XYs1d ): raise TypeError( "Not instance of XYs.XYs1d" )
    return xs.integrateTwoFunctionsWithUncertainty( spectrum, useCovariance=useCovariance, covariance=covariance, normalize=True )


def computeUserDefinedSpectrumAve( xs, path, key, useCovariance=True, covariance=None ):
    """
    User defined spectrum average

    :param xs: cross section to average
    :param path: path to the user's spectra file
    :param key: key of the spectra in the user's spectra file
    :param useCovariance: a flag
    :param covariance: a reference to the covariance of xs, if it cannot be determined otherwise
    :return:
    """
    check_is_cross_section( xs )
    spectrum=get_user_defined_flux( path, key )
    if not isinstance( spectrum, XYs.XYs1d ): raise TypeError( "Not instance of XYs.XYs1d" )
    return xs.integrateTwoFunctionsWithUncertainty( spectrum, useCovariance=useCovariance, covariance=covariance, normalize=True )





# -------------------------------------------------------------------------------
# Unit tests
# -------------------------------------------------------------------------------
if __name__=="__main__":
    import unittest
    from xData.isclose import isclose

    class Test_With_ApproxDiff(unittest.TestCase):

        def assertIsClose(self, a, b, abs_tol=1e-9, rel_tol=1e-8, method='average'):
            if isinstance(a,PQU.PQU) or isinstance(b,PQU.PQU):
                if not ( isinstance(a,PQU.PQU) and isinstance(b,PQU.PQU) ):
                    raise TypeError("Cannot mix types, %s with %s"%(str(type(a)),str(type(b))))
                a=float(a.inUnitsOf(b.unit))
                b=float(b)
            delta = abs(a-b)
            absave = abs(a+b)/2.0
            if absave != 0.0: reldelta = delta/absave
            else: reldelta = 0.0
            self.assertTrue(isclose( a, b, abs_tol=abs_tol, rel_tol=rel_tol, method=method),
                            "%s is not close to %s, absdiff is %s, reldiff is %s"%( str(a), str(b), str(delta), str(reldelta)) )


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
            self.eval=read_evaluation(INTERDIR+'/testData/n-001_H_001.endf')

        def test_read_evaluation(self):
            self.assertItemsEqual(self.eval.keys(),['info', 'reactionSuite', 'errors', 'covarianceSuite'])

        def test_check_is_cross_section(self):
            self.assertIsNone(check_is_cross_section(self.eval['reactionSuite'].getReaction('capture').crossSection))

        def test_cross_section_evaluateWithUncertainty_function(self):
            x = self.eval['reactionSuite'].getReaction('capture').crossSection.evaluateWithUncertainty(PQU.PQU('14.2 MeV'))
            self.assertAlmostEqual(float(x.value),2.9719802000000002e-05,6)
            self.assertAlmostEqual(float(x.uncertainty),2.2832109016510004e-06,6)
            self.assertEqual(str(x.unit),'b')

        def test_cross_section_integrateTwoFunctionsWithUncertainty_function(self):
            for r in self.eval['reactionSuite']:
                if r.ENDF_MT in [1,2,102]:
                    ptwise = r.crossSection.hasLinearForm()
                    if ptwise is None: ptwise = r.crossSection.toPointwise_withLinearXYs(lowerEps=lowerEps, upperEps=upperEps)
                    a = computeRI(r.crossSection,useCovariance=False)
                    b = PQU.PQU(ptwise.integrateWithFunction(DEFAULTONEOVEREFUNCTION,1e-8,[]),'b')
                    self.assertIsClose(a,b)

        def test_RI_vs_INTER(self):
            for r in self.eval['reactionSuite']:
                if r.ENDF_MT==1:
                    self.assertIsClose(computeRI(r.crossSection,domainMax=PQU.PQU(1.00000E+05,'eV'),useCovariance=False), PQU.PQU(2.39596E+02,'b'), rel_tol=1e-5)
                if r.ENDF_MT==2:
                    self.assertIsClose(computeRI(r.crossSection,domainMax=PQU.PQU(1.00000E+05,'eV'),useCovariance=False), PQU.PQU(2.39452E+02,'b'), rel_tol=1e-6)
                if r.ENDF_MT==102:
                    self.assertIsClose(computeRI(r.crossSection,domainMax=PQU.PQU(1.00000E+05,'eV'),useCovariance=False), PQU.PQU(1.48914E-01,'b'), abs_tol=5e-6)

        def test_14MeV(self):
            for r in self.eval['reactionSuite']:
                # Just test the mean values first
                if r.ENDF_MT==1: self.assertIsClose(computeFourteenMeVPoint(r.crossSection,"1.40000E+07 eV",useCovariance=False), PQU.PQU(6.87144E-01,'b'), rel_tol=5e-7)
                if r.ENDF_MT==2: self.assertIsClose(computeFourteenMeVPoint(r.crossSection,"1.40000E+07 eV",useCovariance=False), PQU.PQU(6.87114E-01,'b'), rel_tol=5e-7)
                if r.ENDF_MT==102: self.assertIsClose(computeFourteenMeVPoint(r.crossSection,"1.40000E+07 eV",useCovariance=False), PQU.PQU(2.98051E-05,'b'), rel_tol=5e-7)
                # Now turn on covariance, these tests should still all pass OK
                if r.ENDF_MT==1: self.assertIsClose(computeFourteenMeVPoint(r.crossSection,"1.40000E+07 eV",useCovariance=True), PQU.PQU(6.87144E-01,'b'), rel_tol=5e-7)
                if r.ENDF_MT==2: self.assertIsClose(computeFourteenMeVPoint(r.crossSection,"1.40000E+07 eV",useCovariance=True), PQU.PQU(6.87114E-01,'b'), rel_tol=5e-7)
                if r.ENDF_MT==102: self.assertIsClose(computeFourteenMeVPoint(r.crossSection,"1.40000E+07 eV",useCovariance=True), PQU.PQU(2.98051E-05,'b'), rel_tol=5e-7)

        def test_RoomTemp(self):
            for r in self.eval['reactionSuite']:
                if r.ENDF_MT==1: self.assertIsClose(computeRoomTempCS(r.crossSection), PQU.PQU(2.07687319E+01,'b'), rel_tol=5e-7)
                if r.ENDF_MT==2: self.assertIsClose(computeRoomTempCS(r.crossSection), PQU.PQU(2.043633E+01,'b'), rel_tol=5e-7)
                if r.ENDF_MT==102: self.assertIsClose(computeRoomTempCS(r.crossSection), PQU.PQU(3.322811E-01,'b'), rel_tol=5e-7)

        def test_Cf252Analytic(self):
            for r in self.eval['reactionSuite']:
                if r.ENDF_MT==1: self.assertIsClose(computeCf252SpectrumAveAnalytic(r.crossSection,useCovariance=False), PQU.PQU(3.98828E+00,'b'), rel_tol=0.02)
                if r.ENDF_MT==2: self.assertIsClose(computeCf252SpectrumAveAnalytic(r.crossSection,useCovariance=False), PQU.PQU(3.98824E+00,'b'), rel_tol=0.02)
                if r.ENDF_MT==102: self.assertIsClose(computeCf252SpectrumAveAnalytic(r.crossSection,useCovariance=False), PQU.PQU(3.95532E-05,'b'), rel_tol=0.01)

        def test_Westcott(self):
            """
            Test against values from INTER (recalculated to a higher precision).
            Note: INTER ignores recoil of the target, so uses a=1.0 rather than a=mTarget/(mTarget+mProjectile)
            """
            for r in self.eval['reactionSuite']:
                if r.ENDF_MT==1: self.assertIsClose(computeWestcottFactor(r.crossSection,a=1.0,useCovariance=False), PQU.PQU(1.12677222244,'') )
                if r.ENDF_MT==2: self.assertIsClose(computeWestcottFactor(r.crossSection,a=1.0,useCovariance=False), PQU.PQU(1.12837902942,'') )
                if r.ENDF_MT==102: self.assertIsClose(computeWestcottFactor(r.crossSection,a=1.0,useCovariance=False), PQU.PQU(1.00034149453,'') )

        def test_MACS(self):
            """
            These test values come from B. Pritychenko, S.F. Mughabghab arXiv:1208.2879v3.
            """
            for r in self.eval['reactionSuite']:
                if r.ENDF_MT==102:
                    self.assertIsClose(computeMACS(r.crossSection, PQU.PQU(30.,'keV'),useCovariance=False), PQU.PQU(1.525E-4,'b'), abs_tol=1e-7)
                    self.assertIsClose(computeMACS(r.crossSection, PQU.PQU(1420.,'keV'),useCovariance=False), PQU.PQU(3.892E-5,'b'), abs_tol=2e-8)

                    # Turn on covariance, what happens then?
                    self.assertIsClose(computeMACS(r.crossSection, PQU.PQU(30.,'keV'),useCovariance=True), PQU.PQU("1.524e-4 +/- 5.6e-6 b"), abs_tol=1e-7)
                    break

        def test_scatteringRadius(self):
            """
            These test values come from B. Pritychenko, S.F. Mughabghab arXiv:1208.2879v3.
            """
            self.assertIsClose(computeScatteringRadius(self.eval['reactionSuite']), PQU.PQU(12.76553,'fm'), abs_tol=1e-7)

        def test_ARR(self):
            """
            These test values come from B. Pritychenko, S.F. Mughabghab arXiv:1208.2879v3,
            """
            for r in self.eval['reactionSuite']:
                if r.ENDF_MT==102:
                    self.assertIsClose(computeAstrophysicalReactionRate(r.crossSection, PQU.PQU(30.,'keV'),useCovariance=False), PQU.PQU(3.113E4,"cm**3/s/mol"), rel_tol=0.001)
                    break

        def test_computeCf252SpectrumAve(self):
            """
            These results were compared to the legacy INTER fission spectrum average and we thought our results were
            reasonable so we took our calculated results as the "answer" in the comparison of this test.
            """
            for r in self.eval['reactionSuite']:
                if r.ENDF_MT==1: self.assertIsClose(computeCf252SpectrumAve(r.crossSection,useCovariance=False), PQU.PQU(3.98828E+00,'b'),rel_tol=0.09)
                if r.ENDF_MT==2: self.assertIsClose(computeCf252SpectrumAve(r.crossSection,useCovariance=False), PQU.PQU(3.98824E+00,'b'),rel_tol=0.09)
                if r.ENDF_MT==102: self.assertIsClose(computeCf252SpectrumAve(r.crossSection,useCovariance=False), PQU.PQU(3.95532E-05,'b'),rel_tol=0.09)

        def test_computeGodivaSpectrumAve(self):
            """
            These results are taken from our results.  We think their OK; they are certainly in the ball pack
            (compared to the Cf252 spectrum average above)
            """
            for r in self.eval['reactionSuite']:
                if r.ENDF_MT==102:
                    self.assertIsClose(computeGodivaSpectrumAve(r.crossSection, useCovariance=False), PQU.PQU(3.60870104207e-05,'b'))

        def test_computeJezebelSpectrumAve(self):
            """
            These results are taken from our results.  We think their OK; they are certainly in the ball pack
            (compared to the Cf252 spectrum average above)
            """
            for r in self.eval['reactionSuite']:
                if r.ENDF_MT==102:
                    self.assertIsClose(computeJezebelSpectrumAve(r.crossSection, useCovariance=False), PQU.PQU(3.58646505181e-05,'b'))

        def test_computeBigTenSpectrumAve(self):
            """
            These results are taken from our results.  We think their OK; they are certainly in the ball pack
            (compared to the Cf252 spectrum average above)
            """
            for r in self.eval['reactionSuite']:
                if r.ENDF_MT==102:
                    self.assertIsClose(computeBigTenSpectrumAve(r.crossSection, useCovariance=False), PQU.PQU(3.98659484215e-05,'b'))

    unittest.main()
