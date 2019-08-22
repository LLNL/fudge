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
import bisect

from xData.XYs import XYs1d

class NeedShortCollisionTime( Exception ):
    pass

def toENDL( incoherentInelastic, energyMin_MeV, energyMax_MeV, temperature_MeV ) :

    tolerance = 0.001

    gridded3d = incoherentInelastic.S_alpha_beta.gridded3d
    temperatures = gridded3d.axes[3].copy([])
    temperatures.convertToUnit('MeV/k')
    temperatures = list(temperatures.values)
    betas = numpy.array( gridded3d.axes[2].values )
    beta_interpolation = gridded3d.axes[2].interpolation
    alphas = numpy.array( gridded3d.axes[1].values )
    alpha_interpolation = gridded3d.axes[1].interpolation
    Sab_arr = gridded3d.array.constructArray()

    t_index = bisect.bisect(temperatures, temperature_MeV)
    if abs(temperatures[t_index - 1] - temperature_MeV) < abs(temperatures[t_index] - temperature_MeV):
        t_index -= 1

    # S_alpha_beta does not support temperature interpolation. Override temperature_MeV with the closest match in the file,
    # and print a warning if the match isn't very close.
    t_Kelvin = gridded3d.axes[3].values[t_index]
    if not numpy.allclose(temperatures[t_index], temperature_MeV, atol=1e-12):
        print( "WARNING: requested temperature %s MeV/k not found, using %s MeV/k instead" %
                (temperature_MeV, temperatures[t_index]) )
    temperature_MeV = temperatures[t_index]

    Sab_arr = Sab_arr[t_index,:,:]      # table at desired temperature

    # convert 2-d Sab array into a set of interpolation functions, one for each value of beta.
    # Also store flags for when we need to switch to SCT approximation
    Sab = {}
    SCTflags = {}
    truncatedSab = []
    for idx, beta in enumerate(betas):
        Sab_atBeta = Sab_arr[idx,:]
        if all(Sab_atBeta==0):
            SCTflags[beta] = 1.1 * alphas[-1]   # always use SCT for this beta
        elif any(Sab_atBeta==0):
            nonzero = numpy.nonzero(Sab_atBeta)[0][0] # first non-zero element
            SCTflags[beta] = alphas[nonzero]
            Sab_atBeta = Sab_atBeta[nonzero:]

            zeroes = Sab_atBeta == 0
            if any( zeroes ):
                truncatedSab.append( (beta, sum(zeroes)) )
                Sab_atBeta[ zeroes ] = 1e-90 # to support log-lin or log-log interpolation

            Sab[beta] = XYs1d(zip(alphas[nonzero:], Sab_atBeta), interpolation=alpha_interpolation)
        else:
            Sab[beta] = XYs1d(zip(alphas, Sab_arr[idx,:]), interpolation=alpha_interpolation)

    if truncatedSab:
        betasWithZeroes, nzeroes = zip(*truncatedSab)
        print("WARNING: %d zeroes (in %d betas) were changed to 1e-90" % (sum(nzeroes), len(betasWithZeroes)))

    # also generate parameters for extrapolating to smaller alpha (unless we already switched to SCT)
    # Per Jesse Holmes' email, for small alpha S = C * alpha^N.  Compute C and N for each beta:
    alpha_extrapolation = {}
    for beta in Sab:
        if beta in SCTflags: continue
        xys = Sab[beta]
        a1, s1 = map(math.log, xys[0])
        a2, s2 = map(math.log, xys[1])
        C = (s1 * a2 - s2 * a1) / (a2 - a1)
        N = (s1 - C) / a1
        alpha_extrapolation[beta] = (C,N)

    principalAtom = incoherentInelastic.scatteringAtoms[0]
    secondaryAtoms = list(incoherentInelastic.scatteringAtoms)[1:]

    principalAWR = principalAtom.mass.value    # FIXME GND files should store actual mass, not ratio to neutron mass
    prefactor = principalAtom.numberPerMolecule * principalAtom.boundCrossSection() / (
                2*temperature_MeV)


    def getS_alpha_beta( alpha, beta ):
        """
        Find Sab, either by interpolating on the grid or extrapolating below the alpha grid.
        Raises NeedShortCollisionTime if SCT approximation should be used.
        """
        if alpha <= 0: return 0

        if beta > betas[-1] or alpha > alphas[-1]:  # FIXME what if S is not symmetric in beta?
            raise NeedShortCollisionTime()

        beta_index = bisect.bisect(betas, beta) - 1
        beta1 = betas[beta_index]

        if beta1 in SCTflags and alpha < SCTflags[beta1]:
            raise NeedShortCollisionTime()

        if alpha < alphas[0]:
            C1, N1 = alpha_extrapolation[beta1]
            S1 = math.exp(C1) * alpha ** N1
        else:
            S1 = Sab[beta1].evaluate(alpha)

        if beta1 == beta:
            return S1
        else:
            beta2 = betas[beta_index+1]
            if beta2 in SCTflags and alpha < SCTflags[beta2]:
                raise NeedShortCollisionTime()
            if alpha < alphas[0]:
                C2, N2 = alpha_extrapolation[beta2]
                S2 = math.exp(C2) * alpha ** N2
            else:
                S2 = Sab[beta2].evaluate(alpha)

            xys = XYs1d([[beta1, S1], [beta2, S2]], interpolation='log-lin')    # FIXME hard-coded interpolation
            return xys.evaluate(beta)

    def evaluate( energy, mu, eprime ):
        """
        Get double-differential cross section
        :param energy: incident energy
        :param mu: cos(scattering angle)
        :param eprime: outgoing energy
        :return: sigma(E,mu,Eprime)
        """
        beta = (eprime - energy) / temperature_MeV
        alpha = (eprime + energy - 2 * mu * math.sqrt(energy * eprime)) / (principalAWR * temperature_MeV)

        alpha_lookup, beta_lookup = alpha, beta
        if not incoherentInelastic.asymmetric:  # S(-beta) = S(beta)
            beta_lookup = abs(beta_lookup)
        if incoherentInelastic.calculatedAtThermal:
            factor = temperature_MeV / 2.53e-8
            beta_lookup *= factor
            alpha_lookup *= factor

        try:
            Sab_term = math.exp(-beta/2) * getS_alpha_beta( alpha_lookup, beta_lookup )
        except NeedShortCollisionTime:
            Sab_term = principalAtom.shortCollisionTime(alpha, beta, t_Kelvin)

        result = prefactor * math.sqrt(eprime/energy) * Sab_term

        for atom in secondaryAtoms:
            if Sab_term is None:     # need to use short-collision approximation
                Sab_term = atom.shortCollisionTime(alpha, beta, t_Kelvin)
            result += atom.numberPerMolecule * atom.boundCrossSection() / (
                2*temperature_MeV) * math.sqrt(eprime/energy) * Sab_term
        return result

    def getSpectrum( energy, mu, accuracy = 0.001, epsilon = 1e-6 ):
        """
        Generate appropriate outgoing energy grid, and return spectrum.
        :param energy: incident energy
        :param mu: cos(scattering angle)
        :param accuracy, epsilon: parameters for refining interpolation grid
        :return: list [[eprime0, sigma0],[eprime1,sigma1]...]
        """
        eprimes = set( betas * temperature_MeV + energy )
        eprimes.update( [val for val in -betas * temperature_MeV + energy if val >= 0] )
        if min(eprimes) > 1e-11:
            logmin = math.log10(min(eprimes)) / 2
            ndecades = math.ceil(logmin + 11)
            eprimes.update( list(numpy.logspace(-11,logmin,ndecades * 5) ) )
        eprimes.update( [0] )
        eprimes = sorted(eprimes)

        initialGrid = [[ep,evaluate(energy,mu,ep)] for ep in eprimes]
        yMax = max([v[1] for v in initialGrid])
        yCutoff = 10**math.floor( math.log10( yMax ) - 6 )
        while initialGrid[-2][-1] < yCutoff: initialGrid.pop()
        while initialGrid[1][-1] < yCutoff: initialGrid.pop(0)

        def bisectionTest( res, x1, y1, x2, y2, accuracy, level = 0 ):

            if level > 3: return
            xMid = 0.5 * ( x1 + x2 )
            yMid = evaluate( energy, mu, xMid )
            yEstimate = 0.5 * ( y1 + y2 )
            if abs( yEstimate - yMid ) <= accuracy * yMid: return
            if (xMid-x1) > (2*epsilon*x1): bisectionTest( res, x1, y1, xMid, yMid, accuracy, level + 1 )
            res.append( [xMid, yMid] )
            if (x2-xMid) > (2*epsilon*xMid): bisectionTest( res, xMid, yMid, x2, y2, accuracy, level + 1 )

        res = [ initialGrid.pop(0) ]
        x1,y1 = res[0]
        for x2,y2 in initialGrid:
            bisectionTest( res, x1, y1, x2, y2, accuracy )
            res.append( [x2,y2] )
            x1,y1 = x2,y2

        return res

    elist = numpy.logspace( math.log10(energyMin_MeV), math.log10(energyMax_MeV), 60 )
    mulist = numpy.linspace(-1,1,41)

    I0, I1, I3 = [],[],[]
    for e in elist:
        dsigma_dmu = []
        doublediff = []
        for mu in mulist:
            res = getSpectrum( e, mu )

            xys = XYs1d( data=res )
            dsigma_dmu.append( [ mu, float(xys.integrate()) ] )
            # thin along Eprime
            xys = xys.thin( tolerance )
            doublediff.append( [mu, xys.copyDataToXYs()] )

        xys = XYs1d( data = dsigma_dmu )
        I0.append( [e, float(xys.integrate())] )

        # thin along mu
        xys = xys.thin( tolerance )
        dsigma_dmu = [v for v in dsigma_dmu if v[0] in xys.domainGrid]
        doublediff = [v for v in doublediff if v[0] in xys.domainGrid]

        I1.append( [e, dsigma_dmu] )
        I3.append( [e, doublediff] )

    # TODO add extra incident energies to achieve tolerance

    return ( I0, I1, I3 )
