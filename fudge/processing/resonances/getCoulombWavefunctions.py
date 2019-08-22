#!/usr/bin/env python

# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
# <<END-copyright>>

import math
import numpy
import _getCoulombWavefunctions
from numericalFunctions import specialFunctions as nf_sf

def getCoulombWavefunctions(rho, eta, L):
    # .... Check arguments:
    NE = len(rho)
    for array in (rho,eta):
        if type(array) != numpy.ndarray: raise TypeError("rho and eta must be numpy arrays!")
        if array.dtype != numpy.float64: raise TypeError("rho and eta must be float arrays!")
        if array.shape not in ((NE,), (NE,1)):
            raise TypeError("rho and eta must be 1-d arrays, of the same length!")
    if type(L) is not int: raise TypeError("L must be an integer!")

    F,G = _getCoulombWavefunctions.getCoulombWavefunctions( rho, eta, L )
    F.shape = G.shape = rho.shape
    return (F,G)

SQRTTWOPI = math.sqrt(2.0*math.pi)

def coulombNormalizationFactor(L, eta):
    """
    Coulomb wavefunction normalization factor (see DLMF Eq. 33.2.5 or Abramowitz & Stegun Eq. 14.1.7),
    :math:`C_\\ell(\\eta)={2^\\ell e^{2\pi\\eta}|\Gamma(\\ell+1+i\\eta)|}/{(2\\ell+1)!}`

    Because we're working with numpy arrays of :math:`\\eta` (for parallelization ala' Caleb's tricks), the logic
    of numpy masked arrays is at work here to handle different regimes of eta.

    :param int L: the L to evaluate at
    :param numpy.array(type=float) eta: array of :math:`\eta` values
    :return: array of normalizations
    :rtype: numpy.array(type=float)
    """
    assert isinstance(eta, numpy.ndarray)
    cnorm = numpy.zeros_like(eta)
    ETA_SMALLCUT = 1e-16
    ETA_BIGCUT = 20.0
    small_eta = eta < ETA_SMALLCUT
    medium_eta = numpy.logical_and( eta >= ETA_SMALLCUT, eta < ETA_BIGCUT )
    big_eta  = eta >= ETA_BIGCUT

    # Small eta's are simple, basically just set them to zero in the equation
    if any(small_eta):
        cnorm[small_eta] = pow(2.0,L)*abs(nf_sf.gamma(L+1.0))/nf_sf.gamma(2.0*L+2.0)

    # Medium eta's ... gotta do the real thing (See DLMF, Eq. 33.2.6)
    if any(medium_eta):
        twopieta = 2.0*math.pi*eta[medium_eta]
        eta2=numpy.power(eta[medium_eta], 2)
        stuffinsqrt=1.0
        prefactor=numpy.power(2.0,L)/nf_sf.gamma(2.0*L+2.0)
        for k in range(1,L+1): stuffinsqrt*=eta2+k*k
        cnorm[medium_eta] = prefactor*numpy.sqrt(twopieta*stuffinsqrt/(numpy.exp(twopieta)-1.0))

    # Big eta's, hafta use Stirling's formulas for the gamma function (See Abramowitz & Stegun, Eq. 6.1.39)
    if any(big_eta):
        stirlingsGamma = SQRTTWOPI * numpy.abs( numpy.exp(-1j*eta[big_eta]) * numpy.power(1j*eta[big_eta],1j*eta[big_eta]+L+0.5) )
        cnorm[big_eta] = numpy.power(2.0,L)*numpy.exp(-math.pi*eta[big_eta]/2.0)*stirlingsGamma/nf_sf.gamma(2.0*L+2.0)
        cnorm[numpy.logical_and(big_eta,abs(cnorm)<1e-100)]=1e-100

    return cnorm

def coulombPenetrationFactor(L, rho, eta):
    """
    Compute the Coulomb penetrability, :math:`P_\\ell`,
    :math:`P_\\ell(\\rho,\\eta)={\\rho}/{(A_\\ell(\\rho,\\eta))^2}`
    where :math:`(A_\\ell(\\rho,\\eta))^2=(G_\\ell(\\rho,\\eta))^2+(F_\\ell(\\rho,\\eta))^2`

    Because we're working with numpy arrays of :math:`\eta` (for parallelization ala' Caleb's tricks), the logic
    of numpy masked arrays is at work here to handle different regimes of eta.

    Here we use an external subroutine `coulfg2` from Thompson et al, converted to c and wrapped

    :param L: the L to evaluate at
    :type L: int
    :param rho: array of :math:`\\rho` values
    :type rho: numpy.array(type=float)
    :param eta: array of :math:`\eta` values
    :type eta: numpy.array(type=float)
    :return: array of penetrabilities
    :rtype: numpy.array(type=float)
    """
    _checkRhoAndEta(rho, eta)

    RHOCUT=0.035*(L+1)
    pen = numpy.zeros_like(rho)
    bad_rho = rho  <  RHOCUT
    good_rho = rho >= RHOCUT

    # For small rho, coulfg barfs.  We'll use the asymptotic forms of F & G from DLMF Eqs. 33.5.1 & 33.5.2.
    # We can neglect contribution from F^2 to A^2, since goes like rho^2l for small rho
    if any(bad_rho):
        cnorm = coulombNormalizationFactor(L,eta[bad_rho])
        pen[bad_rho] = numpy.power(rho[bad_rho],2.0*L+1.0) * numpy.power( (2.0*L+1.0)*cnorm, 2.0 )

    # For moderate rho, we're good to go, I think
    if any(good_rho):
        F, G  = getCoulombWavefunctions( rho[good_rho], eta[good_rho], L )
        pen[good_rho] = rho[good_rho]/(F*F+G*G)

    if numpy.any(numpy.isnan(pen)):
        msg = "WARNING: Pen messed up! For L=%i, (rho,eta,pen)=%s"%(L,str(zip(rho.flatten(),eta.flatten(),pen.flatten())))
        if False: raise ValueError(msg)
        pen[numpy.isnan(pen)]=0.0

    return pen

def coulombPhi(L, rho, eta):
    """
    Compute the Coulomb phase, :math:`\\varphi_\\ell`, :math:`\\varphi_\\ell(\\rho,\\eta)=\\arg(G_\\ell(\\rho,\\eta)+iF_\\ell(\\rho,\\eta))`

    Because we're working with numpy arrays of :math:`\\eta` (for parallelization ala' Caleb's tricks), the logic
    of numpy masked arrays is at work here to handle different regimes of eta.

    Here we use an external subroutine 'coulfg2' from Thompson et al, converted to c and wrapped

    :param int L: the L to evaluate at
    :param numpy.array(type=float) rho: array of :math:`\\rho` values
    :param numpy.array(type=float) eta: array of :math:`\\eta` values
    :return: numpy.array(type=float) of phases
    """
    _checkRhoAndEta(rho, eta)

    RHOCUT=0.04*(L+1)
    BIGRHOETA=30.0
    phi = numpy.zeros_like(rho, dtype=float)
    small_rho = rho  <  RHOCUT
    good_rhoeta = numpy.logical_and(rho >= RHOCUT, rho*eta<=BIGRHOETA)
    big_rhoeta = numpy.logical_and( rho >= RHOCUT, rho*eta>BIGRHOETA) #rho*eta > BIGRHOETA

    # For small rho, coulfg barfs.  We'll use the asymptotic forms of F & G from DLMF Eqs. 33.5.1 & 33.5.2.
    if any(small_rho):
        cnorm = coulombNormalizationFactor( L, eta[small_rho] )
        F = cnorm * numpy.power( rho[small_rho], L+1 )
        G = numpy.power( cnorm * numpy.power( rho[small_rho],L)*(2.0*L+1.0), -1.0 ) # by using numpy.power, we can safely handle 1/0
        phi[small_rho] = F/G

    # For moderate rho, we're good to go, I think
    if any(good_rhoeta):
        F, G = getCoulombWavefunctions( rho[good_rhoeta], eta[good_rhoeta], L )
        phi[good_rhoeta] = numpy.arctan2(F,G) #numpy.angle( G+1j*F )

    # For super big rho, have to use asymptotic expression for phase
    # F & G asymptotic forms are from DLMF Eqs. 33.10.1 & 33.10.2
    # arg(Gamma(L+1+ieta) uses A&S Eqs. 6.1.27 and 6.3.2
    if any(big_rhoeta):
        nMax=10*(L+1)*numpy.ones_like(eta[big_rhoeta])
      #  if any(eta[big_rhoeta]>0.0): nMax[eta[big_rhoeta]>0.0] = nMax[eta[big_rhoeta]>0.0] / numpy.round( max(eta[ big_rhoeta ]) )
        nMax=int(max(nMax))

        # Coulomb asymptotic phase
        sigmaL = eta[big_rhoeta]*digamma(L+1) + sum( [
                                                      eta[big_rhoeta]/(L+1.0+n) -
                                                      numpy.arctan( eta[big_rhoeta]/(L+1.0+n) )
                                                      for n in range(0,nMax) ] )

        # and finally phi(rho), asymptotically
        phi[big_rhoeta] = rho[big_rhoeta] - eta[big_rhoeta]*numpy.log(2.0*rho[big_rhoeta]) - numpy.pi*L/2.0 + sigmaL

    # Quality control on pathological phases
    if numpy.any(numpy.isnan(phi)):
        msg = "WARNING: Phase messed up! For L=%i, (rho,eta,phi)=%s"%(L,str(zip(rho.flatten(),eta.flatten(),phi.flatten())))
        if False: raise ValueError(msg)
        phi[ numpy.isnan(phi) ] = 0.0

    return phi

def coulombShiftFactor(L, rho, eta):
    """
    Compute the Coulomb shift, :math:`S_\ell`, :math:`S_\\ell(\\rho,\\eta)=({\\rho}/{A_\\ell(\\rho,\\eta)})({\partial A_\\ell(\\rho,\\eta)}/{\partial\\rho})`
    where :math:`(A_\\ell(\\rho,\\eta))^2=(G_\\ell(\\rho,\\eta))^2+(F_\\ell(\\rho,\\eta))^2`

    Because we're working with numpy arrays of :math:`\eta` (for parallelization ala' Caleb's tricks), the logic
    of numpy masked arrays is at work here to handle different regimes of eta.

    Here we use an external subroutine 'coulfg2' from Thompson et al, converted to c and wrapped

    :param L: the L to evaluate at
    :type L: int
    :param rho: array of :math:`\\rho` values
    :type rho: numpy.array(type=float)
    :param eta: array of :math:`\eta` values
    :type eta: numpy.array(type=float)
    :return: array of shifts
    :rtype: numpy.array(type=float)
    """
    _checkRhoAndEta(rho, eta)
    RHOCUT=0.035*(L+1)
    shift = numpy.zeros_like(rho)
    bad_rho  = rho <  RHOCUT
    good_rho = rho >= RHOCUT

    # For small rho, coulfg barfs.  We'll use the asymptotic forms of F & G from DLMF Eqs. 33.5.1 & 33.5.2.
    # If you work it out, noting that A^2 is dominated by the 1/rho^L behavior of G, it gets real simple
    if any(bad_rho):
        shift[bad_rho] = ( (L+1.0) + eta[bad_rho]*rho[bad_rho]/(L+1.0) -
                           numpy.sqrt(1.0+numpy.power(eta[bad_rho]/(L+1.0),2)) * ((2*L+1.0)/(2*L+3.0)) *
                           coulombNormalizationFactor(L,eta[bad_rho]) /
                           coulombNormalizationFactor(L+1,eta[bad_rho]) )

    # For moderate rho, we're good to go, I think.  These use the derivatives of F & G in the DLMF Eq. 33.4.4
    if any(good_rho):
        F, G   = getCoulombWavefunctions( rho[good_rho], eta[good_rho], L )
        F1, G1 = getCoulombWavefunctions( rho[good_rho], eta[good_rho], L+1 )
        S1 = (L+1.0)/rho[good_rho] + (eta[good_rho])/(L+1.0)
        R1 = numpy.sqrt( 1.0 + numpy.power(eta[good_rho]/(L+1),2.0) )
        shift[good_rho] = rho[good_rho]*( S1 - R1*(G*G1+F*F1)/(F*F+G*G) )

    if numpy.any(numpy.isnan(shift)):
        msg = "WARNING: Shift messed up!  For L=%i, (rho,eta,shift)=%s"%(L,str(zip(rho.flatten(),eta.flatten(),shift.flatten())))
        if False: raise ValueError(msg)
        shift[numpy.isnan(shift)]=-(L+1+eta[numpy.isnan(shift)]*rho[numpy.isnan(shift)]/(L+1)-numpy.sqrt(1+numpy.power(eta[numpy.isnan(shift)]/(L+1),2)))

    return shift

def digamma(n):
    """
    Simple digamma function, tested against A&S tables in chapter 6, up to n=20
    :param n: an integer
    :return:
    """
    return -numpy.euler_gamma + sum( [ 1.0/k for k in range(1, n) ] )

def _checkRhoAndEta(rho, eta):
    """
    Check argument type and size
    """
    NE = len(rho)
    for array in (rho,eta):
        assert isinstance(array, numpy.ndarray) and array.dtype == numpy.float64, "rho and eta must be numpy float arrays"
        if array.shape not in ((NE,), (NE,1)):
            raise TypeError("rho and eta must be 1-d arrays with the same shape!")

#==== Tests ====

def test_getCoulombWavefunctions():
    rho = numpy.arange(1.0,6.0)
    for l in range(0,21):
        print (" L =          %2i"%l)
        for i in range(1,21):
            eta = numpy.array( [0.5*i] * len(rho) )
            F,G = getCoulombWavefunctions(rho,eta,l)
            print ("%25.17E"*len(rho)) % tuple(G)   # print F or G
        print ' '

if __name__ == '__main__':
    test_getCoulombWavefunctions()
