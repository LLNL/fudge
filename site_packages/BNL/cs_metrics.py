from fudge.gnd.reactionData.crossSection import component, pointwise, upperEps
import math
from fudge.core.math.xData import XYs
from pqu import PQU

def funcToXYs( func, fpars,
               Egrid=[1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,1e-1,0.5,1.,5.,10.,50.,1e2,5e2,1e3,5e3,1e4,5e4,1e5,5e5,1e6,5e6,1e7,2e7],
               Eunit='eV', Yunit='1/eV', accuracy=upperEps ):
    '''
    Helper function to convert a user-created function (OK, one of the spectra below) into an XYs instance
    that can be integrated, grouped, whatever.  We pre-defined a energy grid (in eV) that should work well
    even for pathological "spectra" like the problematic 1/E for the resonance integral.
    '''
    return XYs.XYs.createFromFunction( \
            pointwise.defaultAxes(energyUnit=Eunit,crossSectionUnit=Yunit), \
            Xs=Egrid,\
            func=func, \
            parameters=fpars, \
            accuracy=accuracy, \
            biSectionMax=20, \
            checkForRoots = False, \
            infill = 1, \
            safeDivide = 1 )



def checkForCrossSectionType( x ):
    if not isinstance( x, component ):
        raise TypeError( "Not instance of fudge.gnd.reactionData.crossSection.component" )


def computeRI( xs, Ec='0.5 eV', covariance=None ):
    '''
    Compute the resonance integral (RI):
    
    .. math::
        RI = \int_{Ec}^{\infty} \sigma(E) dE/E
        
    Where the lower cut-off is usually taken to be the Cadmium cut-off energy of 
    Ec = 0.5 eV (see S. Mughabghab, Atlas of Neutron Resonances).  If the covariance 
    is provided, the uncertainty on the resonance integral will be computed.
    
    :param string Ec: lower bound of integration.  
        Usually taken to be the Cadmium cut-off energy (default: '0.5 eV')
    :param covariance: covariance to use when computing uncertainty on the spectral average.  
        If None (default: None), no uncertainty is computed.
    :type covariance: covariance instance or None
    :rtype: PQU
    
    '''
    checkForCrossSectionType( xs )
    Ecut = PQU.PQU(Ec).getValueAs(xs.domainUnit() )
    # Construct an XYs instance that spans the same domain as self.
    # The y values are all 1/E and the x values are all E.
    def oneOverEFunc(E,*args): 
        if E < Ecut: return 0.0
        return 1.0/E
    Emin=PQU.PQU(min( xs.domain(xs.domainUnit()) ), xs.domainUnit() )
    Emax=PQU.PQU(max( xs.domain(xs.domainUnit()) ), xs.domainUnit() )
    oneOverE = funcToXYs( oneOverEFunc, [] )
    return xs.computeSpectrumIntegral( oneOverE, Emin=Emin, Emax=Emax, covariance=covariance )


def computeMACS( xs, kT, a=None, covariance=None ):
    '''
    Compute the Maxwellian average cross section (MACS):

    .. math::
        MACS(kT) = {\\frac{2}{\sqrt{\pi}}} {\\frac{a^2}{(kT)^2}} \int_0^{\infty}dE_n^{lab} \sigma(E_n^{lab})*E_n^{lab}*\exp{(-a E_n^{lab}/kT)}

    Here :math:`E_n^{lab}` is the incident neutron energy in the lab frame and :math:`a=m_2/(m_1+m_2)`.
    If the covariance is provided, the uncertainty on the MACS will be computed.

    :param PQU kT: Temperature (as a physical quantity of the Maxwellian with which to average
    :param float a: mass ratio m2/(m1+m2) where m1 is the projectile mass and m2 is the target mass.  
        If set to None, will try to determine from evaluation itself.
    :param covariance: covariance to use when computing uncertainty on the spectral average.  
        If None (default: None), no uncertainty is computed.
    :type covariance: covariance instance or None
    :rtype: PQU
    '''
    if not isinstance( kT, PQU.PQU ): raise TypeError( 'kT must be a PQU, got type %s'%str(type(kT)))
    checkForCrossSectionType( xs )
    
    if a is None: 
        m2 = xs.getRootAncestor().target.getMass('amu')
        m1 = xs.getRootAncestor().projectile.getMass('amu')
        a = m2/(m1+m2)
    
    # Construct an XYs instance that spans the same domain as self.
    # The x values are the E values and the y values are a Maxwellian evaluated at x.
    a_over_kT = a/(kT.getValueAs( xs.domainUnit() ) )
    def maxwellianFunc(E,*args): return E*math.exp(-a_over_kT*E)
    Emin=PQU.PQU(min( xs.domain(xs.domainUnit()) ),xs.domainUnit() )
    Emax=PQU.PQU(max( xs.domain(xs.domainUnit()) ),xs.domainUnit() )
    maxwellian = funcToXYs( maxwellianFunc, [] )
    return xs.computeSpectrumAverage( maxwellian, Emin=Emin, Emax=Emax, covariance=covariance )


def computeCf252SpectrumAve( xs, covariance=None ):
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
    
    :param covariance: covariance to use when computing uncertainty on the spectral average.  
        If None (default: None), no uncertainty is computed.
    :type covariance: covariance instance or None
    :rtype: PQU

    '''
    checkForCrossSectionType( xs )
    # Construct an XYs instance that spans the same domain as self.
    # The x values are the E values and the y values are a Cf252 spectrum evaluated at x.
    def Cf252SpectrumAveFunc(E,*fpars): return math.exp( -0.789 * E/1e6 ) * math.sinh( math.sqrt( 0.447 * E/1e6 ) )
    Emin=PQU.PQU(min( xs.domain(xs.domainUnit()) ), xs.domainUnit())
    Emax=PQU.PQU(max( xs.domain(xs.domainUnit()) ), xs.domainUnit())
    Cf252SpectrumAve = funcToXYs( Cf252SpectrumAveFunc, [] )
    return xs.computeSpectrumAverage( Cf252SpectrumAve, Emin=Emin, Emax=Emax, covariance=covariance )


def computeFourteenMeVPoint( xs, E14='14.2 MeV', covariance=None ):
    '''
    Compute the value of the cross section at 14.2 MeV.

    If the covariance is provided, the uncertainty on the 14.2 MeV point will be computed.
    
    :param covariance: covariance to use when computing uncertainty on the spectral average.  
        If None (default: None), no uncertainty is computed.
    :type covariance: covariance instance or None
    :rtype: PQU
    '''
    checkForCrossSectionType( xs )
    return xs.getValueAtEnergy(PQU.PQU(E14),covariance)


def computeRoomTempCS( xs, ERT='0.0235339347844 eV', covariance=None ): 
    '''
    The cross section exactly at room temperature E=0.0235339347844 eV or v=2200 m/s

    If the covariance is provided, the uncertainty on the 0.0235339347844 eV point will be computed.
    
    RoomTemp = 273.15 degK = 0.0235339347844 eV since kB = 8.6173324e-5 eV/degK

    What temperature does INTER use?  2.53000E-02 eV.   What about Boris? dunno

    :param covariance: covariance to use when computing uncertainty on the spectral average.  
        If None (default: None), no uncertainty is computed.
    :type covariance: covariance instance or None
    :rtype: PQU
    '''
    checkForCrossSectionType( xs )
    return xs.getValueAtEnergy(PQU.PQU(ERT),covariance)


def computeWestcottFactor( xs, a=None, covariance=None ):
    '''
    Wescott G-factor is the ratio of Maxwellian averaged cross section 
    (at room temparature) and the room temperature cross section.  
    Should be pretty close to 1 if cross section goes like 1/v.
    
    INTER also only integrates on interval (1e-5,10) eV

    :param a: mass ratio m2/(m1+m2) where m2 is target mass and m1 is projectile mass (if None, will try to get from evaluation itself)
    :param covariance: covariance to use when computing uncertainty on the spectral average.  
        If None (default: None), no uncertainty is computed.
    :type covariance: covariance instance or None
    :rtype: PQU
    '''
    checkForCrossSectionType( xs )
    kT = PQU.PQU("273.15 K")*PQU.PQU('8.6173324e-5 eV/K')
    macs = xs.computeMACS( a=a, kT=kT, covariance=covariance )
    rt = xs.computeRoomTempCS( str(kT), covariance=covariance ) 
    westcott = (2.0/math.sqrt(math.pi))*macs/rt
    return westcott


def computeAstrophysicalReactionRate( xs, kT, covariance=None ):
    '''
    The astrophysical reaction rate R is defined as :math:`R = N_A <\sigma v>`, where 
    :math:`N_A` is Avagadro's number.  It is typically reported in units of [cm^3/mole s].

    :param PQU T: temperature as a physical quantity
    :param covariance: covariance to use when computing uncertainty on the spectral average.  
        If None (default: None), no uncertainty is computed.
    :type covariance: covariance instance or None
    :rtype: PQU
    
    .. warning:: Not implemented yet
    '''
    checkForCrossSectionType( xs )
    raise NotImplementedError()


def computeStellerEnhancementFactor( xs, covariance=None ):
    '''
    :param covariance: covariance to use when computing uncertainty on the spectral average.  
        If None (default: None), no uncertainty is computed.
    :type covariance: covariance instance or None
    :rtype: PQU
    
    .. warning:: Not implemented yet
    '''
    checkForCrossSectionType( xs )
    raise NotImplementedError()
  
def computeScatteringRadius( xs, covariance=None ):
    """
    R' the scattering radius
    
    :param covariance: covariance to use when computing uncertainty on the spectral average.  
        If None (default: None), no uncertainty is computed.
    :type covariance: covariance instance or None
    :rtype: PQU
    
    .. warning:: Not implemented yet
    """
    checkForCrossSectionType( xs )
    raise NotImplementedError()


def computeSWaveScatteringLength( xs, covariance=None ): 
    '''
    S0, the S-wave scattering length???
    
    :param covariance: covariance to use when computing uncertainty on the spectral average.  
        If None (default: None), no uncertainty is computed.
    :type covariance: covariance instance or None
    :rtype: PQU
    
    .. warning:: Not implemented yet

    '''
    checkForCrossSectionType( xs )
    raise NotImplemented()


def computePWaveScatteringLength( xs, covariance=None ): 
    '''
    S1, the P-wave scattering length???
    
    :param covariance: covariance to use when computing uncertainty on the spectral average.  
        If None (default: None), no uncertainty is computed.
    :type covariance: covariance instance or None     
    :rtype: PQU
    
    .. warning:: Not implemented yet
    '''
    checkForCrossSectionType( xs )
    raise NotImplemented()

