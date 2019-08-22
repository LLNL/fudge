from fudge.gnd.reactionData.crossSection import component, pointwise, lowerEps, upperEps
import math, os, cmath, datatables, collections
from xData import XYs, standards
from numpy import logspace
from pqu import PQU
from fudge.core.utilities.brb import banner, winged_banner

ROOMTEMPERATURE    = PQU.PQU(20.,"degC")
ZERODEGREESCELCIUS = PQU.PQU(0.,"degC")
CADMIUMCUTTOFF     = PQU.PQU(0.5, 'eV')
SPEEDOFLIGHT       = PQU.PQU(299792458., 'm / s')
AVAGADROSNUMBER    = PQU.PQU(6.0221415e23,'1/mol')
INTERDIR = os.path.dirname(__file__)

# switches to turn on not-yet-implemented things
TURNONNEUTRONNSOURCES=False
TURNONFUNDIPPE=False


# -------------------------------------------------------------------------------
# Utilities
# -------------------------------------------------------------------------------

def isclose(a,
            b,
            rel_tol=1e-9,
            abs_tol=0.0,
            method='weak'):
    """
    returns True if a is close in value to b. False otherwise
    :param a: one of the values to be tested
    :param b: the other value to be tested
    :param rel_tol=1e-8: The relative tolerance -- the amount of error
                         allowed, relative to the magnitude of the input
                         values.
    :param abs_tol=0.0: The minimum absolute tolerance level -- useful for
                        comparisons to zero.
    :param method: The method to use. options are:
                  "asymmetric" : the b value is used for scaling the tolerance
                  "strong" : The tolerance is scaled by the smaller of
                             the two values
                  "weak" : The tolerance is scaled by the larger of
                           the two values
                  "average" : The tolerance is scaled by the average of
                              the two values.
    NOTES:
    -inf, inf and NaN behave similar to the IEEE 754 standard. That
    -is, NaN is not close to anything, even itself. inf and -inf are
    -only close to themselves.
    Complex values are compared based on their absolute value.
    The function can be used with Decimal types, if the tolerance(s) are
    specified as Decimals::
      isclose(a, b, rel_tol=Decimal('1e-9'))
    See PEP-0485 for a detailed description

    Test implementation for an isclose() function, for possible inclusion in
    the Python standard library -- PEP0485
    This version has multiple methods in it for experimentation and testing.
    The "final" version can be found in isclose.py
    This implementation is the result of much discussion on the python-ideas list
    in January, 2015:
       https://mail.python.org/pipermail/python-ideas/2015-January/030947.html
       https://mail.python.org/pipermail/python-ideas/2015-January/031124.html
       https://mail.python.org/pipermail/python-ideas/2015-January/031313.html
    Copyright: Christopher H. Barker
    License: Apache License 2.0 http://opensource.org/licenses/apache2.0.php
    """
    if method not in ("asymmetric", "strong", "weak", "average"):
        raise ValueError('method must be one of: "asymmetric",'
                         ' "strong", "weak", "average"')

    if rel_tol < 0.0 or abs_tol < 0.0:
        raise ValueError('error tolerances must be non-negative')

    if a == b:  # short-circuit exact equality
        return True
    # use cmath so it will work with complex or float
    if cmath.isinf(a) or cmath.isinf(b):
        # This includes the case of two infinities of opposite sign, or
        # one infinity and one finite number. Two infinities of opposite sign
        # would otherwise have an infinite relative tolerance.
        return False
    diff = abs(b - a)
    if method == "asymmetric":
        return (diff <= abs(rel_tol * b)) or (diff <= abs_tol)
    elif method == "strong":
        return (((diff <= abs(rel_tol * b)) and
                 (diff <= abs(rel_tol * a))) or
                (diff <= abs_tol))
    elif method == "weak":
        return (((diff <= abs(rel_tol * b)) or
                 (diff <= abs(rel_tol * a))) or
                (diff <= abs_tol))
    elif method == "average":
        return ((diff <= abs(rel_tol * (a + b) / 2) or
                (diff <= abs_tol)))
    else:
        raise ValueError('method must be one of:'
                         ' "asymmetric", "strong", "weak", "average"')


def read_evaluation(endfFile):
    '''
    Read in the evaluation
    :param endfFile: the file name.  This routine will attempt to determine the file type itself.
    :return: the evaluation as a dict of results like {'reactionSuite':the reaction suite, 'covarianceSuite': the covariances}
    '''
    if not os.path.exists( endfFile ): raise IOError( "File named "+endfFile+" doesn't exist!" )
    if open(endfFile).readline().startswith( "<?xml" ):
        from fudge.gnd import reactionSuite
        return {'reactionSuite':reactionSuite.readXML( endfFile ),'covarianceSuite':None,'info':{},'errors':[]}
    else:
        from fudge.legacy.converting import endfFileToGND
        return endfFileToGND.endfFileToGND( endfFile, toStdOut = False, skipBadData = True )


def equal_lethargy_bins(numBins,domainMin=1e-5,domainMax=20.0e6,reverse=False):
    seq=logspace(start=math.log10(domainMin),stop=math.log10(domainMax),num=numBins)
    if reverse: seq.reverse()
    return seq


def convert_units_to_energy(T,energyUnits='eV'):
    if type(T)==str: T=PQU.PQU(T)
    if not isinstance(T,PQU.PQU): raise TypeError( "argument to convert_units_to_energy must be a PQU, got type=%s"%str(type(T)))
    if T.isEnergy(): return T.inUnitsOf(energyUnits)
    elif T.isTemperature(): return (T.inUnitsOf('K')*PQU.PQU(8.6173324e-5, 'eV/K')).inUnitsOf(energyUnits)
#    elif T.isMass(): return
#    elif T.isLength(): return
#    elif T.isTime(): return
    else: raise Exception("Cannot convert %s to energy units" %str(T))


def function_to_XYs( func, fpars,
               Egrid=equal_lethargy_bins(1000),
               domainUnit='eV', domainName='energy_in', rangeUnit='1/eV', rangeName='flux', accuracy=upperEps ):
    '''
    Helper function to convert a user-created function (OK, one of the spectra below) into an XYs instance
    that can be integrated, grouped, whatever.  We pre-defined a energy grid (in eV) that should work well
    even for pathological "spectra" like the problematic 1/E for the resonance integral.
    '''
    return XYs.XYs.createFromFunction( \
            pointwise.defaultAxes(labelsUnits={ \
                XYs.yAxisIndex : ( rangeName, rangeUnit ), \
                XYs.xAxisIndex : ( domainName, domainUnit ) }), \
            Xs=Egrid,\
            func=func, \
            parameters=fpars, \
            accuracy=accuracy, \
            biSectionMax=20, \
            checkForRoots = False, \
            infill = 1, \
            safeDivide = 1 )


def grouped_values_to_XYs( groupBdries, valueList, \
                           domainUnit='eV', domainName='energy_in', rangeUnit='1/eV', rangeName='flux', \
                           accuracy=upperEps, domainMin=1e-5, domainMax=20.0):
    if len(groupBdries) != len(valueList)+1:
        raise ValueError("Group boundries and value lists have incompatable lengths: len(bdries)=%i, len(vals)=%i"%(len(groupBdries),len(valueList)))
    return XYs.XYs(\
        data=[groupBdries, [valueList[0]]+valueList], \
        dataForm="xsandys", \
        interpolation=standards.interpolation.flatToken, \
        axes=pointwise.defaultAxes(labelsUnits={ \
                XYs.yAxisIndex : ( rangeName, rangeUnit ), \
                XYs.xAxisIndex : ( domainName, domainUnit ) }) )

def check_is_cross_section( x ):
    if not isinstance( x, component ):
        raise TypeError( "Not instance of fudge.gnd.reactionData.crossSection.component" )



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

    Ecut = PQU.PQU(Ecut).getValueAs( xs.domainUnit() )
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
        m2 = xs.getRootAncestor().target.getMass('amu')
        m1 = xs.getRootAncestor().projectile.getMass('amu')
        a = float(m2/(m1+m2))
    elif isinstance(a,PQU.PQU): a=float(a.getValueAs(""))
    
    # Construct an XYs instance that spans the same domain as self.
    # The x values are the E values and the y values are a Maxwellian evaluated at x.
    a_over_kT = a/float(kT.getValueAs( xs.domainUnit() ) )
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


def computeFourteenMeVPoint( xs, E14='14.2 MeV', useCovariance=True, covariance=None ):
    '''
    Compute the value of the cross section at 14.2 MeV.

    If the covariance is provided, the uncertainty on the 14.2 MeV point will be computed.
    
    :param useCovariance: use this to override covarance usage
    :type useCovariance: bool
    :param covariance: covariance to use when computing uncertainty on the spectral average.
        If None (default: None), no uncertainty is computed.
    :type covariance: covariance instance or None
    :rtype: PQU
    '''
    check_is_cross_section( xs )
    return xs.evaluateWithUncertainty(PQU.PQU(E14),useCovariance=useCovariance,covariance=covariance)


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

    m2 = xs.getRootAncestor().target.getMass('amu')
    m1 = xs.getRootAncestor().projectile.getMass('amu')
    mu = PQU.PQU(m1*m2/(m1+m2),'amu')

    vT=(2.0*kT/mu).sqrt()
    return (AVAGADROSNUMBER*computeMACS(xs,kT,useCovariance=useCovariance,covariance=covariance)*vT).inUnitsOf( "cm**3/mol/s" )


def computeALPHA(capxs, fisxs, E=convert_units_to_energy(ROOMTEMPERATURE), useCovariance=True):
    return capxs.evaluateWithUncertainty(E, useCovariance=useCovariance)/\
            fisxs.evaluateWithUncertainty(E, useCovariance=useCovariance)


def computeETA(capxs, fisxs, nubar, E=convert_units_to_energy(ROOMTEMPERATURE), useCovariance=True):
    return nubar.evaluate(E)/(1.0+computeALPHA(capxs, fisxs, E=E, useCovariance=useCovariance))


def computeScatteringRadius( elxs, useCovariance=True, covariance=None ):
    """
    Compute R', the scattering radius, from the elastic cross section, assumed to be passed in via the xs argument.

    We take the E=0 point (or at least the lowest available energy point in the elastic cross section).
    Here the cross section is away from resonances and dominated by the potential scattering (==shape elastic).

    with this,
        ..math::
            \sigma_{el}=4\pi(R')^2

    :param elxs: the elastic scattering cross section
    :param useCovariance:
    :param covariance:
    :return:
    """
    check_is_cross_section( elxs )
    return (elxs.evaluateWithUncertainty(PQU.PQU(1e-5,'eV'),useCovariance=useCovariance,covariance=covariance).sqrt()/math.sqrt(4.0*math.pi)).inUnitsOf('fm')


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


def get_KENO_flux( file ):
    energies=[]
    energyUnits='eV'
    fluxes=[]
    fluxUnits='1/eV'
    for line in open(file,mode='r').readlines()[6:305]:
        sline=line.split()
        energies.append(float(sline[1]))
        fluxes.append(float(sline[2]))
    energies.append(1e-5)
    energies.reverse()
    fluxes.reverse()
    return grouped_values_to_XYs( energies, fluxes, domainUnit=energyUnits, rangeUnit=fluxUnits )


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
    spectrum=get_KENO_flux(os.sep.join([INTERDIR,'spectra','HMF001.001']))
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
    spectrum=get_KENO_flux(os.sep.join([INTERDIR,'spectra','PMF001.001']))
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
    spectrum=get_KENO_flux(os.sep.join([INTERDIR,'spectra','IMF007.001']))
    return xs.integrateTwoFunctionsWithUncertainty( spectrum, useCovariance=useCovariance, covariance=covariance, normalize=True )


def computeFUNDIPPESpectrumAve( xs, useCovariance=True, covariance=None ):
    """
    FUND-IPPE spectrum average

    FIXME: THIS IS USELESS UNLESS WE DEFINE WHAT WE MEAN BY "SPECTRUM AVERAGE IN FUND-IPPE" BETTER!! ARE WE
    IN THE MIDDLE OF THE ASSEMBLY?  ON THE OUTSIDE?  HOW LONG DO WE IRRADIATE?

    :param xs: cross section to average
    :param useCovariance: a flag
    :param covariance: a reference to the covariance of xs, if it cannot be determined otherwise
    :return:
    """
    check_is_cross_section( xs )
    spectrum=get_KENO_flux(os.sep.join([INTERDIR,'spectra','FIXME']))
    return xs.integrateTwoFunctionsWithUncertainty( spectrum, useCovariance=useCovariance, covariance=covariance, normalize=True )


def computeDDSourceSpectrumAve( xs, E, useCovariance=True, covariance=None ):
    """
    d(d,n) spectrum average
    :param xs: cross section to average
    :param E: incoming energy to use to induce the source reaction, providing an outgoing spectrum to average with
    :param useCovariance: a flag
    :param covariance: a reference to the covariance of xs, if it cannot be determined otherwise
    :return:
    """
    raise NotImplementedError("%s requires converting P(mu|E) in CM frame to P(E'|E) in lab frame, have P(mu|E) as Legendre moments, so this is very hard" %'computeDDSourceSpectrumAve')
    check_is_cross_section( xs )
    EE=float(E.inUnitsOf('eV'))
    rs = read_evaluation(os.sep.join([INTERDIR,'spectra','d-001_H_002.endf']))['reactionSuite']
    print type(rs.getReaction('n + He3').outputChannel.getProductWithName('n').distributions.toPointwise_withLinearXYs().getAtEnergy(EE))
    raise Exception()
    return PQU.PQU(0.0,'b',0.0)


def computeDTSourceSpectrumAve( xs, E, useCovariance=True, covariance=None ):
    """
    d(t,n) spectrum average
    :param xs: cross section to average
    :param E: incoming energy to use to induce the source reaction, providing an outgoing spectrum to average with
    :param useCovariance: a flag
    :param covariance: a reference to the covariance of xs, if it cannot be determined otherwise
    :return:
    """
    raise NotImplementedError("%s requires converting P(mu|E) in CM frame to P(E'|E) in lab frame, have P(mu|E) as Legendre moments, so this is very hard" %'computeDTSourceSpectrumAve')
    check_is_cross_section( xs )
    EE=float(E.inUnitsOf('eV'))
    rs = read_evaluation(os.sep.join([INTERDIR,'spectra','d-001_H_003.endf']))['reactionSuite']
    print type(rs.getReaction('n + He4').outputChannel.getProductWithName('n').distributions.toPointwise_withLinearXYs().getAtEnergy(EE))
    raise Exception()
    return PQU.PQU(0.0,'b',0.0)


def computePTSourceSpectrumAve( xs, E, useCovariance=True, covariance=None ):
    """
    p(t,n) spectrum average
    :param xs: cross section to average
    :param E: incoming energy to use to induce the source reaction, providing an outgoing spectrum to average with
    :param useCovariance: a flag
    :param covariance: a reference to the covariance of xs, if it cannot be determined otherwise
    :return:
    """
    raise NotImplementedError("%s requires converting P(mu|E) in CM frame to P(E'|E) in lab frame, have P(mu|E) as Legendre moments, so this is very hard" %'computePTSourceSpectrumAve')
    check_is_cross_section( xs )
    EE=float(E.inUnitsOf('eV'))
    rs = read_evaluation(os.sep.join([INTERDIR,'spectra','p-001_H_003.endf']))['reactionSuite']
    print type(rs.getReaction('n + He3').outputChannel.getProductWithName('n').distributions.toPointwise_withLinearXYs().getAtEnergy(EE))
    return PQU.PQU(0.0,'b',0.0)


def computePLiSourceSpectrumAve( xs, E, useCovariance=True, covariance=None ):
    """
     7Li(p,n) Source spectrum average
    :param xs: cross section to average
    :param E: incoming energy to use to induce the source reaction, providing an outgoing spectrum to average with
    :param useCovariance: a flag
    :param covariance: a reference to the covariance of xs, if it cannot be determined otherwise
    :return:
    """
    raise NotImplementedError("%s requires converting P(mu|E) in CM frame to P(E'|E) in lab frame, have P(mu|E) as Legendre moments, so this is very hard" %'computePLiSourceSpectrumAve')
    check_is_cross_section( xs )
    EE=float(E.inUnitsOf('eV'))
    rs = read_evaluation(os.sep.join([INTERDIR,'spectra','p-003_Li_007.endf']))['reactionSuite']
    print type(rs.getReaction('n + Be7').outputChannel.getProductWithName('n').distributions.toPointwise_withLinearXYs().getAtEnergy(EE))
    raise Exception()
    return PQU.PQU(0.0,'b',0.0)



# -------------------------------------------------------------------------------
# Simple HTML helpers
# -------------------------------------------------------------------------------

def XML( tag, content, **kwds ):
    result = ' '.join( [ "<"+str(tag) ] + [ attributeFixer(k)+'="'+str(v)+'"' for k,v in kwds.items() ] )
    if content == None: return result+'/>'
    return result+">"+str(content)+"</"+str(tag)+">"

def H1( x, **kwds ): return XML( 'H1', x, **kwds )

def H2( x, **kwds ): return XML( 'H2', x, **kwds )

def H3( x, **kwds ): return XML( 'H3', x, **kwds )

def H4( x, **kwds ): return XML( 'H4', x, **kwds )

def SCRIPT( x, **kwds ): return XML( 'SCRIPT', x, **kwds )

def OBJECT2HTML(d, nametag="<strong>%s: </strong>", itemtag='<li>%s</li>', valuetag="%s", blocktag=('<ul>', '</ul>')):
    r=[]
    if isinstance(d, dict):
        r.append(blocktag[0])
        for k, v in d.iteritems():
            name = nametag % k
            if isinstance(v, dict) or isinstance(v, list):
                r.append(itemtag % name)
                OBJECT2HTML(r, v)
            else:
                value = valuetag % v
                r.append(itemtag % (name + value))
        r.append(blocktag[1])
    elif isinstance(d, list):
        r.append(blocktag[0])
        for i in d:
            if isinstance(i, dict) or isinstance(i, list):
                r.append(itemtag % " - ")
                OBJECT2HTML(r, i)
            else:
                r.append(itemtag % i)
        r.append(blocktag[1])
    return r



# -------------------------------------------------------------------------------
# Report generator
# -------------------------------------------------------------------------------

class EvaluationReport:

    def __init__( self, endfFile, mtList ):

        rce=read_evaluation(endfFile)

        self.reactionSuite = rce['reactionSuite']
        self.covarianceSuite = rce['covarianceSuite']
        self.info = rce['info']
        self.read_errors = rce['errors']

        self.ENDFFile = endfFile
        self.MTList = mtList

        self.reactionMetricsTable=None
        self.globalMetrics=collections.OrderedDict()
        self.globalMetadata=collections.OrderedDict()
        self.resonanceMetricsTable=None

        self.get_global_metadata()

    def get_global_metadata(self, saveDoc=False):
        if hasattr( self.reactionSuite, 'attributes' ): self.globalMetadata['attributes']=self.reactionSuite.attributes
        else: self.globalMetadata['attributes']={}
        self.globalMetadata['projectile']=self.reactionSuite.projectile
        self.globalMetadata['target']=self.reactionSuite.target
        if saveDoc: self.globalMetadata['endfDoc']=self.reactionSuite.documentation['endfDoc']
        self.globalMetadata['projectileFrame']=self.reactionSuite.getProjectileFrame()
        self.globalMetadata['temperature']=self.reactionSuite.getTemperature('eval')
        if hasattr(self.reactionSuite,'resonances'):
            if hasattr(self.reactionSuite.resonances,'resolved') and self.reactionSuite.resonances.resolved is not None:
                self.globalMetadata['RRRFormat']=str(self.reactionSuite.resonances.resolved.evaluated).split('resolved/')[-1]
                self.globalMetadata['RRRElow']=self.reactionSuite.resonances.resolved.lowerBound
                self.globalMetadata['RRREhigh']=self.reactionSuite.resonances.resolved.upperBound
            if hasattr(self.reactionSuite.resonances,'unresolved') and self.reactionSuite.resonances.unresolved is not None:
                self.globalMetadata['URRFormat']=str(self.reactionSuite.resonances.unresolved.moniker)
                self.globalMetadata['URRElow']=self.reactionSuite.resonances.unresolved.lowerBound
                self.globalMetadata['URREhigh']=self.reactionSuite.resonances.unresolved.upperBound

    def get_reaction_metrics(self, args):
        """

        :param args:
        :return:
        """

        # Final results in a datatable for later report building
        columns=[]
        units=[]
        names=[]
        data=[]

        # Helper function to manage columns and units
        def add_update_columns_units_data(__key,__unit,__point=None,idata=-1,doUncertainty=True):
            if __key not in columns:
                columns.append(__key)
                units.append(__unit)
                if doUncertainty:
                    columns.append('d('+__key+')')
                    units.append(__unit)
            if data is not None and __point is not None:
                if isinstance(__point,PQU.PQU):
                    data[idata].append(__point.value)
                    if doUncertainty:data[idata].append(__point.uncertainty.value)
                else:
                    data[idata].append(__point)
                    if doUncertainty:data[idata].append(0.0)

        import fudge.gnd.reactions.reaction, fudge.gnd.sums

        # ---------------------------------
        #  Process each requested reaction
        # ---------------------------------
        for r in self.reactionSuite:

            # Get rid of all reactions that are not plain reactions or sums of reactions (e.g. no fission components of production stuff)
            if not isinstance(r, (fudge.gnd.reactions.reaction.reaction, fudge.gnd.sums.crossSectionSum) ): continue

            # Keep only the MT's we actually want to test
            mt=r.ENDF_MT
            if args.report == 'legacy' and mt not in args.MTList:
                continue
            elif args.MTList is not None and mt not in args.MTList:
                continue

            # Initialize this row
            names.append(str(r))
            data.append([])

            # Check if we need to add the default columns
            add_update_columns_units_data('MT','',mt,doUncertainty=False)

            # Set Q
            if hasattr( r, 'getQ' ): Q=r.getQ(unit='eV')
            else: Q=0.0
            add_update_columns_units_data('Q','eV',Q,doUncertainty=False)

            # Set EThreshold
            if hasattr( r, "getThreshold" ): EThreshold=r.getThreshold(unit='eV')
            else:                            EThreshold=0.0
            add_update_columns_units_data('Threshold','eV',EThreshold,doUncertainty=False)

            # Resonance Integral
            if args.RI or args.report in ['engineering','summary','legacy']:
                if args.report == 'legacy':
                    RI = computeRI(r.crossSection, domainMax=PQU.PQU(1.00000E+05,'eV'), useCovariance=args.useCovariance)
                else:
                    RI = computeRI(r.crossSection, useCovariance=args.useCovariance)
                add_update_columns_units_data('RI',RI.unit,RI)

            # 14 MeV Point
            if args.fourteen or args.report in ['integral','summary','legacy']:
                E14='14.2 MeV'
                if args.report=='legacy': E14='14.0 MeV'
                fourteen = computeFourteenMeVPoint(r.crossSection, E14, useCovariance=args.useCovariance)
                add_update_columns_units_data(E14,fourteen.unit,fourteen)

            # Room temperature cross section
            if args.thermal or args.report in ['engineering','summary','legacy']:
                xsRT = computeRoomTempCS(r.crossSection, useCovariance=args.useCovariance)
                add_update_columns_units_data("2200 m/s",xsRT.unit,xsRT)

            # Westcott factor
            if args.Westcott or args.report in ['engineering','summary','legacy']:
                if args.report == 'legacy':
                    # Legacy INTER ignores the recoil of the target,
                    # so a = 1.0 rather than a = mTarget/(mTarget+mProjectile)
                    westcott = computeWestcottFactor(r.crossSection,a=PQU.PQU(1.,""), useCovariance=args.useCovariance)
                else:
                    westcott = computeWestcottFactor(r.crossSection, useCovariance=args.useCovariance)
                add_update_columns_units_data("Westcott factor",westcott.unit,westcott)

            # Godiva spectrum average
            if args.Godiva or args.report in ['integral','summary']:
                x=computeGodivaSpectrumAve( r.crossSection, useCovariance=args.useCovariance )
                add_update_columns_units_data("Godiva",x.unit,x)

            # Jezebel spectrum average
            if args.Jezebel or args.report in ['integral','summary']:
                x=computeJezebelSpectrumAve( r.crossSection, useCovariance=args.useCovariance )
                add_update_columns_units_data("Jezebel",x.unit,x)

            # BigTen spectrum average
            if args.BigTen or args.report in ['integral','summary']:
                x=computeBigTenSpectrumAve( r.crossSection, useCovariance=args.useCovariance )
                add_update_columns_units_data("BigTen",x.unit,x)

            # FUND-IPPE spectrum average
            if TURNONFUNDIPPE and (args.FUNDIPPE or args.report in ['integral','summary']):
                x=computeFUNDIPPESpectrumAve( r.crossSection, useCovariance=args.useCovariance )
                add_update_columns_units_data("FUND-IPPE",x.unit,x)

            # 252Cf spectrum average (using analytic approximation from INTER)
            if args.CfSpectAnalytic or args.report in ['legacy']:
                x=computeCf252SpectrumAveAnalytic(r.crossSection, useCovariance=args.useCovariance)
                add_update_columns_units_data("252Cf (analytic)",x.unit,x)

            # 252Cf spectrum average
            if args.CfSpect or args.report in ['integral','summary']:
                x=computeCf252SpectrumAve(r.crossSection, useCovariance=args.useCovariance)
                add_update_columns_units_data("252Cf",x.unit,x)

            if TURNONNEUTRONNSOURCES:

                # d(d,n) spectrum average
                if args.ddSource or args.report in ['nsources','summary']:
                    x=computeDDSourceSpectrumAve( r.crossSection, args.ddSource, useCovariance=args.useCovariance )
                    add_update_columns_units_data("d(d,n) source",x.unit,x)

                # d(t,n) spectrum average
                if args.dtSource or args.report in ['nsources','summary']:
                    x=computeDDSourceSpectrumAve( r.crossSection, args.dtSource, useCovariance=args.useCovariance )
                    add_update_columns_units_data("d(t,n) source",x.unit,x)

                # p(t,n) spectrum average
                if args.ptSource or args.report in ['nsources','summary']:
                    x=computeDDSourceSpectrumAve( r.crossSection, args.ptSource, useCovariance=args.useCovariance )
                    add_update_columns_units_data("p(t,n) source",x.unit,x)

                # 7Li(p,n) Source spectrum average
                if args.pLiSource or args.report in ['nsources','summary']:
                    x=computeDDSourceSpectrumAve( r.crossSection, args.pLiSource, useCovariance=args.useCovariance )
                    add_update_columns_units_data("7Li(p,n) source",x.unit,x)

            # MACS
            if args.MACS or args.report in ['astrophysics','summary']:
                if args.MACS: kT=PQU.PQU(args.MACS,'keV')
                else: kT=PQU.PQU(30.,'keV')
                x=computeMACS(r.crossSection,T=kT, useCovariance=args.useCovariance)
                add_update_columns_units_data("MACS(%s)"%str(kT),x.unit,x)

            # ARR
            if args.ARR or args.report in ['astrophysics','summary']:
                if args.ARR: kT=PQU.PQU(args.ARR,'keV')
                else: kT=PQU.PQU(30.,'keV')
                x=computeAstrophysicalReactionRate(r.crossSection,T=kT, useCovariance=args.useCovariance)
                add_update_columns_units_data("ARR(%s)"%str(kT),x.unit,x)

        self.reactionMetricsTable = datatables.DataTable(data=data,columns=columns,index=names,units=units)

    def get_global_metrics(self, args):

        if args.scatteringRadius or args.report in ['summary']:
            self.globalMetrics["R'"]=computeScatteringRadius(self.reactionSuite.getReaction('elastic').crossSection)

        if args.ALF or args.ETA or args.report in ['engineering','summary']:

            try:
                rfis=self.reactionSuite.getReaction('fission')
                rcap=self.reactionSuite.getReaction('capture')
            except KeyError: return None

            # ALF, ratio of thermal capture to fission
            if args.ALF or args.report in ['engineering','summary'] :
                self.globalMetrics["ALF"]=computeALPHA(rcap.crossSection, rfis.crossSection, useCovariance=args.useCovariance)

            # ETA, nubar_tot*fission cs/absorption cs, not the expanded definition used by EXFOR
            if args.ETA or args.report in ['engineering','summary']:
                totalNubar = sum( [ nubar.multiplicity.toPointwise_withLinearXYs(lowerEps=1e-8, upperEps=1e-8) for nubar in rfis.outputChannel.getProductsWithName('n') ] )
                self.globalMetrics['ETA']=computeETA(rcap.crossSection, rfis.crossSection, totalNubar, useCovariance=args.useCovariance)

    def get_resonance_metrics(self, args):
        """

        :param args:
        :return:
        """

        # Final results in a datatable for later report building
        columns=[]
        units=[]
        names=[]
        data=[]

        # ---------------------------------
        # Process the resonances
        # ---------------------------------
        if args.averageWidth or args.meanLevelSpacing or args.transmissionCoeff or args.report in ['resonance','summary']:

            import fudge.processing.resonances.reconstructResonances as RRReconstruct

            # get URR parameters
            if hasattr(self.reactionSuite.resonances,'unresolved') and self.reactionSuite.resonances.unresolved is not None:
                urr=RRReconstruct.URRcrossSection(self.reactionSuite)
                egrid,flag=urr.generateEnergyGrid()
                urrAveQuant = urr.getWidthsAndSpacings(egrid,interpolateWidths=flag)
            else:
                urr=None
                urrAveQuant=None

            # get RRR parameter averages
            if hasattr(self.reactionSuite.resonances,'resolved') and self.reactionSuite.resonances.resolved is not None:
                resCls = RRReconstruct.getResonanceReconstructionClass(self.reactionSuite.resonances.resolved.evaluated.moniker)
                rrr = resCls( self.reactionSuite, None, enableAngDists=False, verbose=args.verbose )
                rrr.setResonanceParametersByChannel()
                rrrAveQuant = rrr.getAverageQuantities(computeUncertainty=True,resonancesPerBin=50)
            else:
                rrr=None
                rrrAveQuant=None

            def get_key_from_channel(c):
                if int(c.s)==c.s:s=c.s/2.0
                else: s=c.s
                return '%s (j=%s,l=%s,s=%s)'%(str(c.reaction),str(c.J),str(c.l),str(s))

            # Helper function to manage columns and units
            def add_update_columns_units_data(__row,__col,__unit,datum=None):
                if __row not in names:
                    names.append(__row)
                    data.append([])
                if __col not in columns:
                    columns.append(__col)
                    units.append(__unit)
                irow=names.index(__row)
                if data is not None and datum is not None:
                    data[irow].append(datum)

            reactions=[]

            # Initialize the rows
            if rrr is not None:
                for c in rrr.channels:
                    if not c.reaction in reactions: reactions.append(c.reaction)
                    key = get_key_from_channel(c)
                    if key not in names:
                        names.append(key)
                        data.append([])
            if urr is not None:
                for reaction in reactions:
                    for lj in urr.averageWidthFuncs:
                        key = '%s (j=%s,l=%s,s=%s)'%(str(reaction),str(lj[1]),str(lj[0]),'0.5')
                        if key not in names:
                            names.append(key)
                            data.append([])

            if args.neutronStrengthFunction or args.report in ['resonance','summary']:
                print NotImplementedError("neutronStrengthFunction: FIXME")

            if args.gammaStrengthFunction or args.report in ['resonance','summary']:
                print NotImplementedError("gammaStrengthFunction: FIXME")

            if args.DysonMehtaDelta3 or args.report in ['resonance','summary']:
                print NotImplementedError("DysonMehtaDelta3: FIXME")

            if args.effectiveDOF or args.report in ['resonance','summary']:
                print NotImplementedError("effectiveDOF: FIXME")

            if args.averageWidth or args.report in ['resonance','summary']:
                if rrr is not None:
                    for c in rrrAveQuant:
                        add_update_columns_units_data(get_key_from_channel(c),'RRR width','eV','FIXME')
                if urr is not None:
                    for reaction in reactions:
                        for lj in urr.averageWidthFuncs:
                            key = '%s (j=%s,l=%s,s=%s)'%(str(reaction),str(lj[1]),str(lj[0]),'0.5')
                            add_update_columns_units_data(key,'URR width','eV','FIXME')

            if args.meanLevelSpacing or args.report in ['resonance','summary']:
                if rrr is not None:
                    for c in rrrAveQuant:
                        add_update_columns_units_data(get_key_from_channel(c),'RRR D','eV','FIXME')
                if urr is not None:
                    for reaction in reactions:
                        for lj in urr.averageWidthFuncs:
                            key = '%s (j=%s,l=%s,s=%s)'%(str(reaction),str(lj[1]),str(lj[0]),'0.5')
                            add_update_columns_units_data(key,'URR D','eV','FIXME')

            if args.transmissionCoeff or args.report in ['resonance','summary']:
                if rrr is not None:
                    rrrTc = rrr.getTransmissionCoefficientsFromSumRule(computeUncertainty=True)
                    for c in rrrTc:
                        add_update_columns_units_data(get_key_from_channel(c),'RRR Tc','',rrrTc[c])
                if urr is not None:
                    urrTc = urr.getTransmissionCoefficientsFromSumRule()
                    for c in urrTc:
                        add_update_columns_units_data(get_key_from_channel(c),'URR Tc','',urrTc[c])

            if args.neutronPoleStrength or args.report in ['resonance','summary']:
                if rrr is not None:
                    rrrTc = rrr.getNeutronPoleStrength(computeUncertainty=True)
                    for c in rrrTc:
                        add_update_columns_units_data(get_key_from_channel(c),'RRR sc','',rrrTc[c])
                if urr is not None:
                    urrTc = urr.getNeutronPoleStrength()
                    for c in urrTc:
                        add_update_columns_units_data(get_key_from_channel(c),'URR sc','',urrTc[c])

        self.resonanceMetricsTable=datatables.DataTable(data=data,columns=columns,index=names,units=units)

    def text_report(self, args):
        report = banner(self.ENDFFile)
        if args.report is not None: report += '\n'+winged_banner( 'global' )+'\n'
        else: report += '\n'+winged_banner( 'global information' )+'\n'
        for key in self.globalMetadata:
            if key in ['attributes','endfDoc']: continue
            report += key + ' = ' + str(self.globalMetadata[key]) + '\n'
        for key in self.globalMetrics: report += key + ' = ' + str(self.globalMetrics[key]) + '\n'
        report += '\n'

        if self.reactionMetricsTable is not None:
            if args.report is not None: report += winged_banner( args.report+' reactions report' )+'\n'
            else: report += winged_banner('reactions report')+'\n'
            report += str(self.reactionMetricsTable)+'\n'
        report += '\n'

        if self.resonanceMetricsTable is not None:
            if args.report is not None: report += '\n'+winged_banner( args.report+' resonances report' )+'\n'
            else: report += '\n'+winged_banner( 'resonances report' )+'\n'
            report += str(self.resonanceMetricsTable)+'\n'

        return report

    def json_report(self,args):
        import json
        report = []
        if self.globalMetrics: report.append('"metadata":'+json.dumps({ x:str(self.globalMetadata[x]) for x in self.globalMetadata}))
        if self.reactionMetricsTable is not None: report.append('"reaction":'+self.reactionMetricsTable.to_json())
        if self.globalMetrics: report.append('"addendum":'+json.dumps({ x:str(self.globalMetrics[x]) for x in self.globalMetrics}))
        if self.resonanceMetricsTable is not None: report.append('"resonance":'+self.resonanceMetricsTable.to_json())
        report = '{'+','.join(report)+'}'
        return report

    def csv_report(self,args):
        results = None
        if self.reactionMetricsTable is not None: resuls = self.reactionMetricsTable.to_csv()
        if self.globalMetrics: print "WARNING: csv_report cannot output global metrics"
        if self.resonanceMetricsTable is not None: print "WARNING: csv_report cannot output resonance metrics"
        return results

    def html_report(self,args):
        report = [H1(self.ENDFFile)]
        if self.globalMetrics:
            report.append(H2("Global information"))
            report.append('\n'.join(OBJECT2HTML(self.globalMetadata)))
            report.append('\n'.join(OBJECT2HTML(self.globalMetrics)))
        if self.reactionMetricsTable is not None:
            report.append(H2("Reaction report"))
            report.append(self.reactionMetricsTable.to_html())
        if self.resonanceMetricsTable is not None:
            report.append(H2("Resonance report"))
            report.append(self.resonanceMetricsTable.to_html())
        #report = '\n'.join(report)
        return report



# -------------------------------------------------------------------------------
# Unit tests
# -------------------------------------------------------------------------------
if __name__=="__main__":
    import unittest

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


    class Test_Utilities(Test_With_ApproxDiff):

        def test_equal_lethargy_bins(self):
            self.maxDiff=None
            answer = [1.0e-05, 0.00037606030930863936, 0.014142135623730951, 0.53182958969449889, 20.0]
            calc=equal_lethargy_bins(5,domainMax=20.0)
            for i in range(len(calc)):
                self.assertAlmostEqual(answer[i],calc[i],7)

        def test_function_to_XYs(self):
            def oneOverEFunc(E,*args):
               return 1.0/E
            x = function_to_XYs(oneOverEFunc,[])
            for E in [1e-4, 2.3e-3, 33e-2, 0.1, 19.0]:
                self.assertIsClose(x.evaluate(E),1.0/E)
            def myExp(E,*args):
                return math.exp(-E/1e3)
            y = function_to_XYs(myExp,[])
            for E in [33e-2, 0.1, 20.]:
                self.assertIsClose(y.evaluate(E),math.exp(-E/1e3))

        def test_grouped_values_to_XYs(self):
            self.assertIsClose(0.0,0.0)

        def test_integrateTwoFunctions(self):
            def myExp1(E,*args):
                return math.exp(-E)
            def myExp2(E,*args):
                return math.exp(-E/2.0)
            a=function_to_XYs(myExp1,[]).integrateTwoFunctions(function_to_XYs(myExp2,[])).inUnitsOf('1/eV')
            b=PQU.PQU(0.666657,'1/eV')
            self.assertIsClose(a,b,abs_tol=5e-7)

        def test_one_over_E(self):
            """
            Mathematica says the answer is \int_0.5^20e6 dE/E = 17.5044
            :return:
            """
            def oneFunc(x,*args): return 1.0
            one=function_to_XYs(oneFunc,[])
            self.assertIsClose(DEFAULTONEOVEREXYs.integrateTwoFunctions(one), PQU.PQU(17.5044,'1/eV'),rel_tol=1e-6)


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
            self.eval=read_evaluation(__file__.split(os.sep)[0]+'/testData/n-001_H_001.endf')

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
                    print computeMACS(r.crossSection, PQU.PQU(30.,'keV'),useCovariance=True)
                    break

        def test_scatteringRadius(self):
            """
            These test values come from B. Pritychenko, S.F. Mughabghab arXiv:1208.2879v3.
            """
            for r in self.eval['reactionSuite']:
                if r.ENDF_MT==2:
                    self.assertIsClose(computeScatteringRadius(r.crossSection,useCovariance=True), PQU.PQU(12.7525380409,'fm'), abs_tol=1e-7)
                    break

        def test_ARR(self):
            """
            These test values come from B. Pritychenko, S.F. Mughabghab arXiv:1208.2879v3,
            """
            for r in self.eval['reactionSuite']:
                if r.ENDF_MT==102:
                    self.assertIsClose(computeAstrophysicalReactionRate(r.crossSection, PQU.PQU(30.,'keV'),useCovariance=False), PQU.PQU(3.113E4,"cm**3/s/mol"), rel_tol=0.001)
                    break

        @unittest.skip("Nothing to test yet, the code isn't written")
        def test_SEF(self):
            """
            FIXME: "answers" here are fictitous, please get real values and put them in here
            """
            for r in self.eval['reactionSuite']:
                if r.ENDF_MT==1: self.assertIsClose(computeStellerEnhancementFactor(r.crossSection,useCovariance=False), 2.39596E+02)
                if r.ENDF_MT==2: self.assertIsClose(computeStellerEnhancementFactor(r.crossSection,useCovariance=False), 2.39452E+02)
                if r.ENDF_MT==102: self.assertIsClose(computeStellerEnhancementFactor(r.crossSection,useCovariance=False), 1.48914E-01)

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

        @unittest.skip("not written, probably won't ever be")
        def test_computeDDSourceSpectrumAve(self):
            for r in self.eval['reactionSuite']:
                if r.ENDF_MT==102:
                    self.assertIsClose(computeDDSourceSpectrumAve(r.crossSection, PQU.PQU(1.,'MeV'),useCovariance=False), PQU.PQU(0.0,'b'))

        @unittest.skip("not written, probably won't ever be")
        def test_computeDTSourceSpectrumAve(self):
            for r in self.eval['reactionSuite']:
                if r.ENDF_MT==102:
                    self.assertIsClose(computeDTSourceSpectrumAve(r.crossSection, PQU.PQU(1.,'MeV'),useCovariance=False), PQU.PQU(0.0,'b'))

        @unittest.skip("not written, probably won't ever be")
        def test_computePTSourceSpectrumAve(self):
            for r in self.eval['reactionSuite']:
                if r.ENDF_MT==102:
                    self.assertIsClose(computePTSourceSpectrumAve(r.crossSection, PQU.PQU(1.,'MeV'),useCovariance=False), PQU.PQU(0.0,'b'))

        @unittest.skip("not written, probably won't ever be")
        def test_computePLiSourceSpectrumAve(self):
            for r in self.eval['reactionSuite']:
                if r.ENDF_MT==102:
                    self.assertIsClose(computePLiSourceSpectrumAve(r.crossSection, PQU.PQU(1.,'MeV'),useCovariance=False), PQU.PQU(0.0,'b'))


    unittest.main()