#!/usr/bin/env python
#encoding: utf-8

# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# <<END-copyright>>

"""
reconstruct resonance region
cmattoon, 12/1/2010

See ENDF-102 (endf documentation) appendix D for equations

Basic usage:
------------

    if x is a reactionSuite instance,::

        >>xsecs = fudgeReconstructResonances.reconstructResonances(x, tolerance=0.01)

    This reconstructs all resolved/unresolved resonance sections. In each section,
    the results are accurate under linear interpolation to tolerance of 1% or better.
    All sections are summed together.::

        * xsecs = dictionary containing cross sections: {'total':XYs, 'elastic':XYs, ...}
    
    each cross section is an XYs class instance, containing data and also axes with units


Alternate use:
--------------

    If desired, reconstruct a single section::

        >> resCls = fudgeReconstructResonances.RMcrossSection(x)   # for Reich_Moore
        >> energy_grid = s.generateEnergyGrid()
        >> crossSections = s.getCrossSection( energy_grid )
        # the input to getCrossSection is the energy (or list of energies) in eV
        # crossSections are returned as a dictionary {'total':,'elastic':,'capture':,'fission':,}

        # improve grid to desired tolerance for linear interpolation:
        >> new_energy_grid, new_crossSections, messages = s.refineInterpolation(energy_grid, crossSections, tolerance=0.01)

"""
import numpy, collections
from fudge.gnd.reactionData import crossSection

__metaclass__ = type


ChannelDesignator = collections.namedtuple( 'ChannelDesignator', [ 'l', 'J', 'reaction', 'index', 's', 'gfact', 'particleA', 'particleB', 'useCoulomb', 'useRelativistic', 'eliminated' ] ) 


def getAllowedTotalSpins( L, S, useFactor2Trick=True ): 
    '''
    Returns a list of allowed J values from summing angular momenta L and S, where 
        .. math::
            \vec{J}=\vec{L}+\vec{S}
        
    which implies
        ..math::
            |L-S| \leq J \leq L+S 
    
    The useFactor2Trick flag tells the routine whether we are summing real angular momenta or momenta * 2.
    If the useFactor2Trick flag is true, then momenta are really momenta*2, meaning they can be pure integers, 
    even if the real momenta refer to 1/2-integer values (e.g. spin).  The default is to useFactor2Trick because
    most C/C++/Fortran codes that compute angular momentum-stuff use the trick so they can use integer math.  
    Also, it makes the use of the Python range() function possible. 
    '''
    if useFactor2Trick: return range( abs(L-S), L+S+2, 2 )
    else: return [ j/2.0 for j in getAllowedTotalSpins( int(L*2), int(S*2) ) ]


def reconstructResonances(reactionSuite, tolerance=None, enableAngDists=False, verbose=True):
    """ reconstruct all resonance sections (resolved/unresolved) in reactionSuite, 
    add results together for full (resonance region) pointwise cross section.
    If tolerance is specified, refine grid to required tolerance for (lin-lin) interpolation """
    
    if not reactionSuite.resonances.reconstructCrossSection:
        print ("      Resonance reconstruction was already done. Reconstruct again? y/N")
        if not raw_input().lower().startswith('y'): # bug!!!
            return xsecs
    
    egrids, xsecs, angdists = [], [], []
    if reactionSuite.resonances.resolved:
        def resolvedReconstruct( nativeData, sectionIndex = None ):
            if nativeData.moniker == 'SingleLevel_BreitWigner': resCls = SLBWcrossSection
            elif nativeData.moniker == 'MultiLevel_BreitWigner': resCls = MLBWcrossSection
            elif nativeData.moniker == 'Reich_Moore': resCls = RMcrossSection
            elif nativeData.moniker == 'R_Matrix_Limited': resCls = RMatrixLimitedcrossSection
            else:
                raise TypeError, "Don't recognize resonance type %s" % nativeData.moniker
            reconstructClass = resCls( reactionSuite, sectionIndex, enableAngDists=enableAngDists, verbose=verbose )
            egrid = reconstructClass.generateEnergyGrid()
            xsecs_now = reconstructClass.getCrossSection( egrid )
            if enableAngDists:
                estep = 10
                angdists_now = [ reconstructClass.getAngularDistribution(E) for E in egrid[::estep] ]
            else: angdists_now = None
            if tolerance:
                egrid, xsecs_now, messages = reconstructClass.refineInterpolation(numpy.array(egrid), xsecs_now, tolerance)
                if verbose:
                    for message in messages: print (message)
            return egrid, xsecs_now, angdists_now
        
        if reactionSuite.resonances.resolved.multipleRegions:
            if verbose:
                print( "      WARNING! Multiple resolved/unresolved energy regions are deprecated\n"
                        +"        and should be consolidated into single section")
            for RRidx in range(len(reactionSuite.resonances.resolved.regions)):
                nativeData = reactionSuite.resonances.resolved.regions[RRidx].nativeData
                egrid_now, xsecs_now, angdists_now = resolvedReconstruct( nativeData, RRidx )
                # merge with 'full' energy grid:
                egrids.append( egrid_now )
                xsecs.append( xsecs_now )
                if enableAngDists: angdists.append( angdists_now )
        else:
            nativeData = reactionSuite.resonances.resolved.nativeData
            egrid_now, xsecs_now, angdists_now = resolvedReconstruct( nativeData )
            # merge with 'full' energy grid:
            egrids.append( egrid_now )
            xsecs.append( xsecs_now )
            if enableAngDists: angdists.append( angdists_now )

    if reactionSuite.resonances.unresolved:
        nativeData = reactionSuite.resonances.unresolved.nativeData
        reconstructClass = URRcrossSection( reactionSuite, verbose )
        if nativeData.forSelfShieldingOnly: # don't reconstruct
            if verbose:
                print ("Skipping unresolved: for self shielding only")
        else:
            xsecs_now = reconstructClass.getCrossSection( )

            if False and tolerance:
                # Disabling this section since ENDF manual now states that we should interpolate cross sections rather than parameters for URR.
                egrid, xsecs_now, messages = reconstructClass.refineInterpolation(numpy.array(egrid), xsecs_now,
                        tolerance)
                if verbose:
                    for message in messages: print (message)

            # meld with resolved region:
            egrids.append( None )
            xsecs.append( xsecs_now )

    # 'xsecs' list now holds one or more regions. Convert to pointwise or piecewise cross section instance:
    xsecs_final = {}
    if len(xsecs)==1:   # only one region: treat as pointwise
        for key in xsecs[0]:
            if isinstance( xsecs[0][key], crossSection.pointwise ): xsecs_final[key] = xsecs[0][key]
            else:
                axes_ = crossSection.pointwise.defaultAxes()
                xsecs_final[key] = crossSection.pointwise( axes_, (egrids[0], xsecs[0][key]), dataForm="XsAndYs",
                        accuracy=(tolerance or 0.01) )
    else:               # multiple regions: treat as piecewise
        for key in xsecs[0]:
            pwxs = crossSection.piecewise( crossSection.piecewise.defaultAxes() )
            for idx in range(len(xsecs)):
                if isinstance( xsecs[idx][key], crossSection.pointwise ): xys = xsecs[idx][key]
                else:
                    axes_ = crossSection.pointwise.defaultAxes()
                    xys = crossSection.pointwise( axes_, (egrids[idx], xsecs[idx][key]), dataForm="XsAndYs",
                            accuracy=(tolerance or 0.01) )
                pwxs.append( xys )
            xsecs_final[key] = pwxs
    
    # 'angdists' list holds one or more regions.  Convert to one Legendre coefficient table
    if enableAngDists:
        angdists_final = {}

    # note that this does not add pointwise (ENDF-MF3) portion.
    # Use gnd.reactionSuite.reconstructResonances for that
    return xsecs_final


# base class common to resolved and unresolved resonance reconstruction
class resonanceReconstructionBaseClass:

    def __init__(self, reactionSuite):
        self.projectile = reactionSuite.projectile
        self.target = reactionSuite.target
        self.spin = self.target.getSpin().value
        neutronMass = reactionSuite.getParticle( 'n' ).getMass('amu')
        self.targetToNeutronMassRatio = self.target.getMass('amu') / neutronMass
    
    def k(self, energy):
        '''
        For an incident neutron, with energy in eV, the ENDF manual states

            ..math ::
                k = \frac{\sqrt{2m_n}}{\hbar}\frac{AWRI}{AWRI+1}\sqrt{|E|}

        sqrt(2*neutronMass)/hbar == 2.196807e-3 (eV*barn)**-1/2. Thus for energy in eV, k is in b**-1/2
        ''' 
        return (2.196807122623e-3 * self.targetToNeutronMassRatio / 
                (self.targetToNeutronMassRatio+1) * numpy.sqrt(energy))
    
    """ refer to SAMMY manual page 9 for penetration/shift/phase factor equations: """
    def penetrationFactor(self, L, rho):
        if L==0: return rho
        elif L==1: return rho**3 / (1+rho**2)
        elif L==2: return rho**5 / (9+3*rho**2+rho**4)
        elif L==3: return rho**7 / (225+45*rho**2+6*rho**4+rho**6)
        elif L==4: return rho**9 / (11025+1575*rho**2+135*rho**4+10*rho**6+rho**8)
        elif L>4:
            # find it recursively
            P_Lminus1 = self.penetrationFactor(L-1,rho)
            S_Lminus1 = self.shiftFactor(L-1,rho)
            return rho**2 * P_Lminus1 / ((L-S_Lminus1)**2 + P_Lminus1**2)
    
    def shiftFactor(self, L, rho):
        """ calculate shift factor used in SLBW and MLBW formalisms """
        if L==0: return 0.0*rho   # rho may be an array
        elif L==1: return -1.0 / (1+rho**2)
        elif L==2: return -(18.0+3*rho**2) / (9+3*rho**2+rho**4)
        elif L==3: return -(675.0+90*rho**2+6*rho**4) / (225+45*rho**2+6*rho**4+rho**6)
        elif L==4: return -( (44100.0+4725*rho**2+270*rho**4+10*rho**6) / 
                (11025+1575*rho**2+135*rho**4+10*rho**6+rho**8) )
        elif L>4:
            P_Lminus1 = self.penetrationFactor(L-1,rho)
            S_Lminus1 = self.shiftFactor(L-1,rho)
            return (rho**2 * (L-S_Lminus1) / ((L-S_Lminus1)**2+P_Lminus1**2) - L)
    
    def phi(self, L, rho):
        # calculate hard-sphere phase-shift
        if L==0: return rho
        elif L==1: return rho-numpy.arctan(rho)
        elif L==2: return rho-numpy.arctan(3*rho/(3-rho**2))
        elif L==3: return rho-numpy.arctan(rho*(15-rho**2)/(15-6*rho**2))
        elif L==4: return rho-numpy.arctan(rho*(105-10*rho**2)/(105-45*rho**2+rho**4))
        elif L>4:
            P_Lminus1 = self.penetrationFactor(L-1,rho)
            S_Lminus1 = self.shiftFactor(L-1,rho)
            phi_Lminus1 = self.phi(L-1,rho)
            return phi_Lminus1 - numpy.arctan(P_Lminus1 / (L-S_Lminus1))
    
    # for use after the cross section has been calculated on the initial grid:
    def refineInterpolation(self, egrid, xsecs, tolerance=0.01):
        """ generateEnergyGrid may not give a fine enough grid to linearly interpolate to desired tolerance.
        My solution to that: for all consecutive points (x0,y0), (x1,y1) and (x2,y2) do a linear interpolation between
        (x0,y0) and (x2,y2). If the interpolation doesn't agree with (x1,y1) within tolerance,
        subdivide up the region by adding two more calculated points.  Keep iterating until interpolation agrees within tolerance.
        
        This means that in the end we will have more points than required for given tolerance.
        The results can be thinned (implemented in XYs)
        """
        def checkInterpolation(x,y):
            """ do a linear interpolation of each point from its two nearest neighbors
            then check where the interpolation is insufficient
            x,y must be numpy arrays
            """
            m = (y[2:]-y[:-2])/(x[2:]-x[:-2])
            delta_x = x[1:-1] - x[:-2]
            b = y[:-2]

            # add first and last points back in, for easier comparison with
            # original y values:
            interpolated = numpy.zeros_like(x)
            interpolated[1:-1] = m*delta_x+b
            interpolated[0] = y[0]; interpolated[-1] = y[-1]

            # find where original and interpolated grids differ by more than tolerance
            # silence div/0 warnings for this step, since xsc = 0 case is explicitly handled below
            olderr = numpy.seterr( divide='ignore', invalid='ignore' )
            delt = interpolated / y
            numpy.seterr( **olderr )    # re-enable div/0 warnings

            mask = (delt>1+tolerance) + (delt<1-tolerance) # boolean array
            badindices = numpy.arange(len(mask))[ mask ]
            # at these points, we need a finer mesh

            # use absolute convergence condition at or near xsc = 0:
            zeros = (y[badindices-1]==0) + (y[badindices]==0) + (y[badindices+1]==0)
            if any(zeros):
                ignore = []
                for idx in badindices[ zeros ]:
                    if abs(y[idx]-y[idx-1])<1e-3 and abs(y[idx+1]-y[idx])<1e-3:
                        mask[ idx ] = False
                        ignore.append( idx )
                badindices = list(badindices)
                for idx in ignore[::-1]:
                    badindices.remove(idx)

            midpoints = (x[:-1]+x[1:])/2
            energies_needed = sorted( set( list( midpoints[mask[1:]]) + list( midpoints[mask[:-1]]) ) )
            return badindices, energies_needed

        def mergeArrays( newIdx, x1, x2, y1, y2 ):
            # from the list of 'bad' points in original array x1, obtain two new integer arrays, a1 and a2,
            # containing indices to insert x1 and x2 into merged array. Use these to quickly merge all data:
            newLen = len(x1)+len(x2)
            a1,a2 = [],[]
            nic = iter( newIdx )
            if len(newIdx)==1: idnow, idnext = nic.next(), -1
            else: idnow, idnext = nic.next(), nic.next()
            current, offset = idnow, 0
            extraPoint = False

            for val in xrange(newLen):
                if val==current:
                    a2.append(val)
                    offset += 1
                    if idnext==idnow+1:
                        current=idnext + offset
                        idnow = idnext
                        try: idnext = nic.next()
                        except StopIteration: idnext = -1
                    elif extraPoint:
                        extraPoint = False
                        current=idnext + offset
                        idnow = idnext
                        try: idnext = nic.next()
                        except StopIteration: idnext = -1
                    else:
                        extraPoint = True
                        current += 2
                else:
                    a1.append(val)

            # a1 and a2 now hold the indices for inserting x1 and x2 into merged array
            xnew = numpy.zeros(newLen)
            xnew[a1] = x1; xnew[a2] = x2
            ynew = {}
            for key in y1:
                ynew[key] = numpy.zeros(newLen)
                ynew[key][a1] = y1[key]; ynew[key][a2] = y2[key]
            return xnew, ynew

        messages = []
        reactionDone = dict.fromkeys(xsecs)  # values are None
        n_iter = 0
        addedPoints = 0
        while True:
            newX, newIdx = set(), set()
            for key in xsecs:
                if not any(xsecs[key]):
                    continue
                if not reactionDone[key]:
                    badindices, energies_needed = checkInterpolation(egrid,xsecs[key])
                    if len(energies_needed)==0: reactionDone[key] = True
                    newIdx.update( badindices )
                    newX.update( energies_needed )
            newX = sorted(newX)
            newIdx = sorted(newIdx)
            if len(newX)==0:    # success!
                break
            if n_iter > 20:
                messages.append("Iteration limit exceeded when refining interpolation grid!")
                break
            n_iter += 1
            addedPoints += len(newX)
            newY = self.getCrossSection( newX )
            # merge new x/y values with original list:
            egrid, xsecs = mergeArrays( newIdx, egrid, newX, xsecs, newY )

        messages.append("%i points were added (for total of %i) to achieve tolerance of %s%%" % 
            (addedPoints, len(egrid), tolerance*100))
        return egrid, xsecs, messages


#### base class for resolved resonance reconstruction ####
class RRBaseClass(resonanceReconstructionBaseClass):

    def __init__(self, reactionSuite, sectionIndex=None):
        super(RRBaseClass,self).__init__(reactionSuite)
        """ store resonance parameters in convenient structure for quick cross section calculations: """

        if sectionIndex is not None:
            # energy boundaries for this region (multiple regions are deprecated):
            energyRegion = reactionSuite.resonances.resolved.regions[sectionIndex]
            self.RR = energyRegion.nativeData
            self.lowerBound = energyRegion.lowerBound.getValueAs('eV')
            self.upperBound = energyRegion.upperBound.getValueAs('eV')
        else:
            # only one energy region:
            resolved = reactionSuite.resonances.resolved
            self.RR = resolved.nativeData
            self.lowerBound = resolved.lowerBound.getValueAs('eV')
            self.upperBound = resolved.upperBound.getValueAs('eV')
    
    def sortLandJ(self):
        """ 
        SLBW, MLBW and Reich_Moore formalisms have similar structure
        it's convenient to sort their resonances by L and J
        This method should NOT be used for R-Matrix Limited
        """
        # make a unified table of resonance parameter data, convert to eV if necessary
        nRes = len( self.RR.resonanceParameters )
        params = ('energy','L','J','channelSpin','totalWidth','neutronWidth','captureWidth',
                'fissionWidthA','fissionWidthB')
        units = ('eV','','','','eV','eV','eV','eV','eV')
        data = [self.RR.resonanceParameters.getColumn( quant,unit ) for quant,unit in zip(params,units) ]
        for i in range(len(data)):
            if data[i] is None: data[i] = [0]*nRes
        table = numpy.array( data )

        # sort resonances by L and J, store parameters in numpy arrays
        self.Ls = []
        Llist = table[ params.index('L') ]
        for L in range(int(max( Llist ))+1):
            Lres = table[ :, table[ params.index('L') ]==L ]
            Jlist = Lres[ params.index('J') ]
            Js = []
            for J in sorted(set(Jlist)):
                LJres = Lres[ :, Lres[ params.index('J') ]==J ]
                spinList = LJres[ params.index('channelSpin') ]
                spins = []
                for spin in sorted(set(spinList)):
                    spinRes = LJres[ :, LJres[ params.index('channelSpin') ]==spin ]
                    energies = spinRes[ params.index('energy') ]
                    neutronWidth = ( spinRes[ params.index('neutronWidth') ] / 
                        self.penetrationFactor(L,self.rho(abs(energies),L)) )
                    # totalWidth = numpy.array([v(res.totalWidth) for res in LJres])
                    spindict = {
                            'channelSpin': spin,
                            'energy': energies,
                            'neutronWidth': neutronWidth,
                            'captureWidth': spinRes[params.index('captureWidth')],
                            'fissionWidthA': spinRes[params.index('fissionWidthA')],
                            'fissionWidthB': spinRes[params.index('fissionWidthB')],
                            'shiftFactor': 0.5*self.shiftFactor(L,self.rho(numpy.abs(energies))),
                            }
                    spins.append( spindict )
                Jdict = {
                        'J': J,
                        'gfact': (2.0*abs(J)+1)/(2*(2*self.spin+1)),
                        'channelSpins': spins,
                        }
                Js.append( Jdict )
            
            Ldict = { 'L': L, 'Js': Js }
            self.Ls.append( Ldict )

        # for each L, gfactors should sum to 2*L+1.
        self.missingGfactor = {}
        for L in self.Ls:
            gsum = sum( [J['gfact'] * len(J['channelSpins']) for J in L['Js']] )
            self.missingGfactor[L['L']] = 2*L['L']+1 - gsum

        # for energy grid generation:
        self._energies = table[ params.index('energy') ]
        totalWidths = table[ params.index('totalWidth') ]
        if not any(totalWidths): # total widths aren't specified in Reich_Moore case
            totalWidths = ( table[params.index('neutronWidth')] + table[params.index('captureWidth')]
                    + table[params.index('fissionWidthA')] + table[params.index('fissionWidthB')] )
        self._widths = totalWidths
    
    def rho(self, E, L=None):
        '''get the channel radius, rho. If L is specified try to get L-dependent value'''
        if self.RR.calculateChannelRadius:
            a = 0.123 * self.target.getMass('amu')**(1./3.) + 0.08  # eq. D.14 in ENDF manual, a in b^-1/2
        else:
            if self.RR.scatteringRadius.isEnergyDependent():
                a = self.RR.scatteringRadius.getValueAs('10**-12*cm', E[:,0]) # a in b^-1/2
            else: a = self.RR.getScatteringRadius(L).getValueAs('10**-12*cm') # a in b^-1/2
        return self.k(E) * a # dimensionless
    
    def generateEnergyGrid(self):
        """ Create an initial energy grid by merging a rough mesh for the entire region (~10 points / decade)
        with a denser grid around each resonance. For the denser grid, multiply the total resonance width by
        the 'resonancePos' array defined below. """
        energies, widths = self._energies, self._widths
        lowBound, highBound = self.lowerBound, self.upperBound
        # ignore negative resonances
        for lidx in range(len(energies)):
            if energies[lidx] > 0: break
        energies = energies[lidx:]
        widths = widths[lidx:]
        # generate grid for a single peak, should be good to 1% using linear interpolation:
        resonancePos = numpy.array([ 
            0.000e+00, 1.000e-03, 2.000e-03, 3.000e-03, 4.000e-03, 5.000e-03, 6.000e-03, 7.000e-03, 8.000e-03, 9.000e-03, 1.000e-02, 2.000e-02,
            3.000e-02, 4.000e-02, 5.000e-02, 6.000e-02, 7.000e-02, 8.000e-02, 9.000e-02, 1.000e-01, 1.100e-01, 1.200e-01, 1.300e-01, 1.400e-01,
            1.500e-01, 1.600e-01, 1.700e-01, 1.800e-01, 1.900e-01, 2.000e-01, 2.100e-01, 2.200e-01, 2.300e-01, 2.400e-01, 2.500e-01, 2.600e-01,
            2.800e-01, 3.000e-01, 3.200e-01, 3.400e-01, 3.600e-01, 3.800e-01, 4.000e-01, 4.200e-01, 4.400e-01, 4.600e-01, 4.800e-01, 5.000e-01,
            5.500e-01, 6.000e-01, 6.500e-01, 7.000e-01, 7.500e-01, 8.000e-01, 8.500e-01, 9.000e-01, 9.500e-01, 1.000e+00, 1.050e+00, 1.100e+00,
            1.150e+00, 1.200e+00, 1.250e+00, 1.300e+00, 1.350e+00, 1.400e+00, 1.450e+00, 1.500e+00, 1.550e+00, 1.600e+00, 1.650e+00, 1.700e+00,
            1.750e+00, 1.800e+00, 1.850e+00, 1.900e+00, 1.950e+00, 2.000e+00, 2.050e+00, 2.100e+00, 2.150e+00, 2.200e+00, 2.250e+00, 2.300e+00,
            2.350e+00, 2.400e+00, 2.450e+00, 2.500e+00, 2.600e+00, 2.700e+00, 2.800e+00, 2.900e+00, 3.000e+00, 3.100e+00, 3.200e+00, 3.300e+00,
            3.400e+00, 3.600e+00, 3.800e+00, 4.000e+00, 4.200e+00, 4.400e+00, 4.600e+00, 4.800e+00, 5.000e+00, 5.200e+00, 5.400e+00, 5.600e+00,
            5.800e+00, 6.000e+00, 6.200e+00, 6.400e+00, 6.500e+00, 6.800e+00, 7.000e+00, 7.500e+00, 8.000e+00, 8.500e+00, 9.000e+00, 9.500e+00,
            1.000e+01, 1.050e+01, 1.100e+01, 1.150e+01, 1.200e+01, 1.250e+01, 1.300e+01, 1.350e+01, 1.400e+01, 1.450e+01, 1.500e+01, 1.550e+01,
            1.600e+01, 1.700e+01, 1.800e+01, 1.900e+01, 2.000e+01, 2.100e+01, 2.200e+01, 2.300e+01, 2.400e+01, 2.500e+01, 2.600e+01, 2.700e+01,
            2.800e+01, 2.900e+01, 3.000e+01, 3.100e+01, 3.200e+01, 3.300e+01, 3.400e+01, 3.600e+01, 3.800e+01, 4.000e+01, 4.200e+01, 4.400e+01,
            4.600e+01, 4.800e+01, 5.000e+01, 5.300e+01, 5.600e+01, 5.900e+01, 6.200e+01, 6.600e+01, 7.000e+01, 7.400e+01, 7.800e+01, 8.200e+01,
            8.600e+01, 9.000e+01, 9.400e+01, 9.800e+01, 1.020e+02, 1.060e+02, 1.098e+02, 1.140e+02, 1.180e+02, 1.232e+02, 1.260e+02, 1.300e+02,
            1.382e+02, 1.550e+02, 1.600e+02, 1.739e+02, 1.800e+02, 1.951e+02, 2.000e+02, 2.100e+02, 2.189e+02, 2.300e+02, 2.456e+02, 2.500e+02,
            2.600e+02, 2.756e+02, 3.092e+02, 3.200e+02, 3.469e+02, 3.600e+02, 3.892e+02, 4.000e+02, 4.200e+02, 4.367e+02, 4.600e+02, 4.800e+02,
            5.000e+02, 6.000e+02, 7.000e+02, 8.000e+02, 9.000e+02, 1.000e+03 ])
        
        grid = []
        # get the midpoints (on log10 scale) between each resonance:
        # emid = [lowBound] + list(10**( ( numpy.log10(energies[1:])+numpy.log10(energies[:-1]) ) / 2.0)) + [highBound]
        # or get midpoints on linear scale:
        emid = [lowBound] + [(e1+e2)/2.0 for e1,e2 in zip(energies[1:],energies[:-1])] + [highBound]
        for e,w,lowedge,highedge in zip(energies,widths,emid[:-1],emid[1:]):
            points = e-w*resonancePos
            grid += [lowedge] + list(points[points>lowedge])
            points = e+w*resonancePos[1:]
            grid += list(points[points<highedge])
        # also add rough grid, to cover any big gaps between resonances, should give at least 10 points per decade:
        npoints = numpy.ceil( numpy.log10(highBound)-numpy.log10(lowBound) ) * 10
        grid += list( numpy.logspace(numpy.log10(lowBound),numpy.log10(highBound), npoints) )[1:-1]
        grid += [lowBound, highBound, 0.0253]   # region boundaries + thermal
        grid = sorted(set(grid))
        # toss any points outside of energy bounds:
        grid = grid[ grid.index(lowBound) : grid.index(highBound)+1 ]
        return grid

    def setResonanceParametersByChannel( self, multipleSScheme = 'NJOY', useReichMooreApproximation = False ):
        '''
        Reorganize member data into channels (relies heavily on groundwork in sortLandJ).
        '''
        self.nResonances = len( self.RR.resonanceParameters )
        self.lMax = self.RR.LvaluesNeededForConvergence 
        self.channelConstantsBc = []
        allowedSs = getAllowedTotalSpins( 0.5, self.spin, useFactor2Trick=False )

        # Make a unified table of resonance parameter data, convert to eV if necessary
        lList = []
        params = ('energy','L','J','channelSpin','totalWidth','neutronWidth','captureWidth', 'fissionWidthA','fissionWidthB')
        units = ('eV','','','','eV','eV','eV','eV','eV')
        data = [ self.RR.resonanceParameters.getColumn( quant,unit ) for quant,unit in zip(params,units) ]
        for i in range(len(data)):
            if data[i] is None: data[i] = [0]*self.nResonances    

        def addOrUpdateDict( theDict, theKey, theValue ):
            if not theKey in theDict: theDict[ theKey ] = {}
            theDict[ theKey ][ theValue[0] ] = theValue[1]

        # Now make the dictionary of channels
        channelDict = collections.OrderedDict()
        for iR in range( self.nResonances ):
            ER = data[0][iR]
            L =  data[1][iR]
            J =  data[2][iR]
            S =  None # not given in MLBW, RM or RML, I just list it here for completeness
            GT = data[4][iR]
            GN = data[5][iR]
            GG = data[6][iR]
            GF = data[7][iR]
            GX = data[8][iR]
            JisAllowed = False
            self.lMax = max( self.lMax, L )
            if L not in lList: lList.append(L)
            allowedSByThisJl = []
            for SS in allowedSs:
                if J in getAllowedTotalSpins( L, SS, useFactor2Trick=False ): 
                    JisAllowed = True
                    allowedSByThisJl.append( SS )
            if not JisAllowed: raise ValueError( "Invalid total angular momentum and/or orbital momentum: cannot couple up to J = "+str(J)+" with I = "+str(self.spin)+", i = 0.5 and L = "+str(L) )
            gfact = (2.0*abs(J)+1)/(2*(2*self.spin+1))
            numAllowedSByThisJl = len( allowedSByThisJl )
            
            # Way recommended by NJOY:  
            #   * One of valid S's get width, other(s) get zero width
            #   * Zero width channels left in sum to get potential scattering correct
            #   * Gets good (cs from U) to (cs from Fudge) comparison; good (cs from U) to (cs from BB) 
            if multipleSScheme == 'NJOY':  
                addOrUpdateDict( channelDict, ChannelDesignator( L, J, 'elastic', 0, allowedSByThisJl[0], gfact, None, None, False, False, False ), ( iR, GN ) )
                addOrUpdateDict( channelDict, ChannelDesignator( L, J, 'capture', 0, allowedSByThisJl[0], gfact, None, None, False, False, useReichMooreApproximation ), ( iR, GG ) )
                if GF != 0.0: addOrUpdateDict( channelDict, ChannelDesignator( L, J, 'fission', 0, allowedSByThisJl[0], gfact, None, None, False, False, False ), ( iR, GF ) )
                if GX != 0.0: addOrUpdateDict( channelDict, ChannelDesignator( L, J, 'competitive', 0, allowedSByThisJl[0], gfact, None, None, False, False, False ), ( iR, GX ) )
                for SS in allowedSByThisJl[1:]:
                    addOrUpdateDict( channelDict, ChannelDesignator( L, J, 'elastic', 0, SS, gfact, None, None, False, False, False ), ( iR, 0.0 ) )
                    addOrUpdateDict( channelDict, ChannelDesignator( L, J, 'capture', 0, SS, gfact, None, None, False, False, useReichMooreApproximation ), ( iR, 0.0 ) )
                    if GF != 0.0: addOrUpdateDict( channelDict, ChannelDesignator( L, J, 'fission', 0, SS, gfact, None, None, False, False, False ), ( iR, 0.0 ) )
                    if GX != 0.0: addOrUpdateDict( channelDict, ChannelDesignator( L, J, 'competitive', 0, SS, gfact, None, None, False, False, False ), ( iR, 0.0 ) )
            
            # Way hinted at by ENDF manual: 
            #   * Width divided equally between channels with valid S's
            #   * All channels left in sum to get potential scattering correct
            #   * Gets poor (cs from U) to (cs from Fudge) comparison; poor (cs from U) to (cs from BB) 
            elif multipleSScheme == 'ENDF':   
                for SS in allowedSByThisJl:
                    addOrUpdateDict( channelDict, ChannelDesignator( L, J, 'elastic', 0, SS, gfact, None, None, False, False, False ), ( iR, GN/numAllowedSByThisJl ) )
                    addOrUpdateDict( channelDict, ChannelDesignator( L, J, 'capture', 0, SS, gfact, None, None, False, False, useReichMooreApproximation ), ( iR, GG/numAllowedSByThisJl ) )
                    if GF != 0.0: addOrUpdateDict( channelDict, ChannelDesignator( L, J, 'fission', 0, SS, gfact, None, None, False, False, False ), ( iR, GF/numAllowedSByThisJl ) )
                    if GX != 0.0: addOrUpdateDict( channelDict, ChannelDesignator( L, J, 'competitive', 0, SS, gfact, None, None, False, False, False ), ( iR, GX/numAllowedSByThisJl ) )
            
            # Ignore problem: 
            #   * Ignore spin of channels and the fact may be multiple valid spins
            #   * Gets best (cs from U) to (cs from Fudge) comparison; poor (cs from U) to (cs from BB) 
            else :  
                addOrUpdateDict( channelDict, ChannelDesignator( L, J, 'elastic', 0, S, gfact, None, None, False, False, False ), ( iR, GN ) )
                addOrUpdateDict( channelDict, ChannelDesignator( L, J, 'capture', 0, S, gfact, None, None, False, False, useReichMooreApproximation ), ( iR, GG ) )
                if GF != 0.0: addOrUpdateDict( channelDict, ChannelDesignator( L, J, 'fission', 0, S, gfact, None, None, False, False, False ), ( iR, GF ) )
                if GX != 0.0: addOrUpdateDict( channelDict, ChannelDesignator( L, J, 'competitive', 0, S, gfact, None, None, False, False, False ), ( iR, GX ) )

        # Take care of extra l's needed to get potential scattering correct
        # Make a single elastic channel with a zero width resonance at the energy of the first real resonance
        # This won't add to the resonance scattering part, but it will make a potential scattering part.
        # We add in channels for all spins and total angular momenta to make sure get the angular momentum algebra correct.
        for l in range( 0, self.lMax+1 ):
            for s in allowedSs:
                for j in getAllowedTotalSpins( l, s, useFactor2Trick=False ):
                    newChannel = ChannelDesignator( l, j, 'elastic', 0, s, (2.0*abs(j)+1)/(2*(2*self.spin+1)), None, None, False, False, False )
                    if not newChannel in channelDict: addOrUpdateDict( channelDict, newChannel, ( 0, 0 ) )

        # Set up the channel->resonance mappings
        self.allChannels = channelDict
        self.channels = collections.OrderedDict()  # just the kept ones
        self.eliminatedChannels = collections.OrderedDict()
        for c in self.allChannels:
            if c.eliminated:    self.eliminatedChannels[ c ] = self.allChannels[ c ]
            else:               self.channels[ c ] = self.allChannels[ c ]
        self.nChannels = len( self.channels )
        self.identityMatrix = numpy.identity( self.nChannels, dtype=complex )

    def getChannelConstantsBc( self ):
        raise NotImplementedError( "Override in derived classes if used" )
        
    def getRMatrix( self, Ein ): raise NotImplementedError( "Override in derived classes if used" )

    def getKMatrix( self, Ein ): raise NotImplementedError( "Override in derived classes if used" )

    def getAMatrix( self, Ein ): raise NotImplementedError( "Override in derived classes if used" )

    def getL0Matrix( self, Ein ):
        '''
        Get the L0 matrix of Froehner:
        
            ..math::
                
                {\bf L^0}_{cc'} = \delta_{cc'} (L_c-B_c) 

        where 
    
            ..math::
                
                L_c = S_c + i P_c

        '''
        if self.channelConstantsBc == []: self.channelConstantsBc = self.getChannelConstantsBc()
        L0 = []
        for ic,c in enumerate(self.channels):
            rho = self.rho( Ein, c.l )
            L0.append( complex( self.shiftFactor( c.l, rho ), self.penetrationFactor( c.l, rho ) ) - self.channelConstantsBc[ ic ] )
        return numpy.diag( L0 )
        
    def getXMatrix( self, Ein ): 
        '''
        Get the X matrix for use in computing W.  X is:

            ..math::
                
                {\bf X}_{cc'} = P^{-1/2}_c ( ( {\bf I} - {\bf R}{\bf L^0} )^{-1}{\bf R} )_{cc'} P_{c'}^{-1/2}\delta_{JJ'}
        '''
        L0 = self.getL0Matrix( Ein )
        
        # Precompute phase factor, sqrt(Penetrability) 
        penFact = []
        for ic, c in enumerate( self.channels ):
            if c.reaction == 'elastic': penFact.append( numpy.sqrt( L0[ic,ic].imag ) )
            else:                       penFact.append( 1.0 )

        # Compute the (unscaled) X matrix
        R = self.getRMatrix( Ein )
        X = numpy.linalg.inv( self.identityMatrix - R * L0 ) * R
    
        # Apply the sqrt(Penetrability) rescaling
        for ic1, c1 in enumerate( self.channels ):
            for ic2, c2 in enumerate( self.channels ):
                if c1.J != c2.J: continue
                X[ ic1, ic2 ] *= penFact[ ic1 ] * penFact[ ic2 ]
                
        return X

    def getWMatrix( self, Ein ): 
        '''
        Get the W matrix for use in computing U.  W is:
        
            ..math::
                
                {\bf W} = {\bf I} + 2i{\bf X}
        '''
        return self.identityMatrix + complex( 0.0, 2.0 ) * self.getXMatrix( Ein )
                  
    def getEiPhis( self, Ein, useTabulatedScatteringRadius = True ):
        '''
        The phase factor for the collision matrix:
        
            ..math::
                
                \Omega_c = e^{-\varphi_c}
        '''
 
        # Precompute phase factor
        eiphis = []
        for ic, c in enumerate( self.channels ):
            if c.reaction == 'elastic': 
                # For calculating phi, default is to use tabulated scattering radius:
                if useTabulatedScatteringRadius:
                    if self.RR.scatteringRadius.isEnergyDependent():
                        rhohat = self.RR.scatteringRadius.getValueAs('10*fm', Ein) * self.k(Ein)
                    else:
                        rhohat = self.RR.getScatteringRadius(L=c.l).getValueAs('10*fm') * self.k(Ein)
                else: rhohat = self.rho(Ein, L=c.l)
                eiphis.append( numpy.exp( complex( 0.0, -self.phi( c.l, rhohat ) ) ) )
            else:
                eiphis.append( complex( 1.0, 0.0 ) )
        return eiphis
                
    def getScatteringMatrixU( self, Ein, useTabulatedScatteringRadius = True ):
        '''
        Compute the scattering matrix using 
        '''
        # Initialize lists and matrices
        eiphis =  self.getEiPhis( Ein, useTabulatedScatteringRadius=useTabulatedScatteringRadius )        
        
        # This adds in the phase factors to get us U
        U = numpy.zeros( ( self.nChannels, self.nChannels ), dtype=complex )
        W = self.getWMatrix( Ein )
        for ic1 in range( self.nChannels ):
            for ic2 in range( self.nChannels ):
                U[ ic1, ic2 ] = W[ ic1, ic2 ] * eiphis[ ic1 ] * eiphis[ ic2 ]
        return U

    def getAngularDistribution(self, E, keepL0Term=True, renormalize=True, reaction1='elastic', reaction2='elastic' ):
        from numericalFunctions import angularMomentumCoupling as nf_amc
        k = self.k(E)                                           # Relative momentum, comes out in b^-1/2
        projSpin = 0.5                                          # Assuming that we have a neutron as a projectile, projSpin is 1/2
        targSpin = self.spin                                    # Target spin from baseclass
        prefactor = 1.0/k/k/(2*projSpin+1)/(2*targSpin+1)       # Prefactor for cross section, gets units, main energy dependence, and spin dep. prefactor
        B = []                                                  # List of Blatt-Biedenharn coefficients, used in Legendre series expansion of the angular distribution
        relTol = 1e-5                                           # Relative tolerance used to decide when to cut off Legendre series
        Lmax = min( max( [ xL['L'] for xL in self.Ls ] ), 10 )  # Lmax also used to decide when to cut off Legendre series
        allowedSpins = getAllowedTotalSpins( targSpin, projSpin, False )    # Possible values of the total channel spin, based on the projectile spin (projSpin) and target spin (targSpin)
        U = self.getScatteringMatrixU( E )                      # The scattering matrix
        L = 0
        while len(B)<=2 or ( B[-1]>=relTol*B[0] and L<=Lmax ):
            B.append(0.0)
            for S in allowedSpins:
                for Sp in allowedSpins:
                    spinFactor =  pow( -1, S - Sp ) / 4.0
                    for i1, c1 in enumerate( self.channels ):
                        if c1.eliminated: continue
                        if c1.reaction != reaction1: continue
                        if c1.s != None and c1.s != S: continue
                        for i2, c2 in enumerate( self.channels ):
                            if c2.eliminated: continue
                            if c2.reaction != reaction2: continue
                            if c2.s != None and c2.s != S: continue
                            Z = nf_amc.z_coefficient( int(2.*c1.l), int(2.*c1.J), int(2.*c2.l), int(2.*c2.J), int(2.*S), int(2.*L) )
                            if Z == 0.0: continue
                            for i1p, c1p in enumerate( self.channels ):
                                if c1p.eliminated: continue
                                if c1p.reaction != reaction1: continue
                                if c1.J != c1p.J: continue
                                if c1p.s != None and c1p.s != Sp: continue
                                for i2p, c2p in enumerate( self.channels ):
                                    if c2p.eliminated: continue
                                    if c2p.reaction != reaction2: continue
                                    if c2.J != c2p.J: continue
                                    if c2p.s != None and c2p.s != Sp: continue
                                    Zp = nf_amc.z_coefficient( int(2*c1p.l), int(2*c1p.J), int(2*c2p.l), int(2*c2p.J), int(2*Sp), int(2*L) )
                                    if Zp == 0.0: continue
                                    thisB = spinFactor * Z * Zp * ( ( (c1==c1p) - U[i1,i1p] ).conj() * ( (c2==c2p) - U[i2,i2p] ) ).real
                                    B[-1] += thisB
            L += 1
        return [ prefactor * x for x in B ]


        
"""
@blockwise: function decorator for improving performance in resolved region.
Each 'getCrossSection' method is wrapped by this function.
If we have lots of incident energies, this splits up the calculation using the multiprocessing module.
May still need to tweak the 'NE' and 'nprocesses' variables for best performance.
"""
def blockwise(function):
    def wrapped(self,E):
        if numpy.isscalar(E):
            E = numpy.array([[E]])
            return function(self,E)
        else:
            NE = len(E)
            # turn E into a column vector
            E = numpy.array(E).reshape(NE,1)
            if NE < 1000: # faster to run directly
                return function(self,E)
            else:
                from multiprocessing import Process, Queue
                queue = Queue()

                def enqueue_result(elist, Slice, queue):
                    # perform the calculation, put result in the queue
                    result = function(self, elist)
                    queue.put( (Slice, result) )
                    return

                jobs = []
                nprocesses = 8  # number of processes to spawn
                # how many energies does each process calculate?
                chunk = int(numpy.ceil((NE+0.) / nprocesses))
                for i in range(nprocesses): # start the calculations
                    Slice = slice(i*chunk, (i+1)*chunk)
                    p = Process(target=enqueue_result, args=(E[Slice], Slice, queue))
                    jobs.append(p)
                    p.start()

                for i in range(nprocesses): # collect results as they finish
                    Slice, vals = queue.get()
                    if i==0:
                        result = dict.fromkeys( vals )
                        for key in result: result[key] = numpy.zeros(NE)
                    for key in result:
                        result[key][Slice] = vals[key]

                # allow child processes to exit:
                for job in jobs: job.join()
                return result
    return wrapped


#### Single-level Breit-Wigner ###
class SLBWcrossSection(RRBaseClass):
    """
    given a resonance region in SLBW format,
    create a class with all data required to reconstruct
    
    Note, the resonances in the SLBW format each correspond to different "levels" and so do not
    interfere.  This is unlike all of the other RRR formats.  So, while one resonance energy 
    and one width go with one channel and all the channels associated with one resonance energy 
    go together to make one reaction, the reactions for each resonance are added incoherently (no
    interference effects at all).
    """
    def __init__(self, reactionSuite, sectionIndex=None, enableAngDists=False, verbose=True):
        super(SLBWcrossSection, self).__init__(reactionSuite, sectionIndex)
        self.verbose = verbose
        self.sortLandJ()
        if self.verbose:
            print ("From %f to %f eV, reconstructing using Single-level Breit-Wigner" % 
                    (self.lowerBound,self.upperBound))
        if enableAngDists: self.setResonanceParametersByChannel()
    
    def setResonanceParametersByChannel( self, multipleSScheme = 'NJOY' ):
        '''
        Reorganize member data into channels (relies heavily on groundwork in sortLandJ).
        
        Note, unlike the getResonanceParametersByChannel() function in MLBW, RM or RML, 
        the fact that different resonances are entirely different reactions means that the 
        channelDicts have to have an additional layer of sorting that corresponds to the SLBW "level".
        '''
        # make a unified table of resonance parameter data, convert to eV if necessary
        self.nResonances = len( self.RR.resonanceParameters )
        params = ('energy','L','J','channelSpin','totalWidth','neutronWidth','captureWidth',
                'fissionWidthA','fissionWidthB')
        units = ('eV','','','','eV','eV','eV','eV','eV')
        data = [self.RR.resonanceParameters.getColumn( quant,unit ) for quant,unit in zip(params,units) ]
        for i in range(len(data)):
            if data[i] is None: data[i] = [0]*self.nResonances

        allowedSs = getAllowedTotalSpins( 0.5, self.spin, useFactor2Trick=False )

        # Now make the dictionary of channels, grouped by independent level
        channelDicts = []
        for iR in range( self.nResonances ):
            ER = data[0][iR]
            L =  data[1][iR]
            J =  data[2][iR]
            S =  None # not given in SLBW, I just list it here for completeness
            GT = data[4][iR]
            GN = data[5][iR]
            GG = data[6][iR]
            GF = data[7][iR]
            GX = data[8][iR]
            JisAllowed = False
            allowedSByThisJl = []
            for SS in allowedSs:
                if J in getAllowedTotalSpins( L, SS, useFactor2Trick=False ): 
                    JisAllowed = True
                    allowedSByThisJl.append( SS )
            if not JisAllowed: raise ValueError( "Invalid total angular momentum and/or orbital momentum: cannot couple up to J = "+str(J)+" with I = "+str(self.spin)+", i = 0.5 and L = "+str(L) )
            gfact = (2.0*abs(J)+1)/(2*(2*self.spin+1))
            channelDict=collections.OrderedDict()
            numAllowedSByThisJl = len( allowedSByThisJl )
            
            # Way recommended by NJOY:  
            #   * One of valid S's get width, other(s) get zero width
            #   * Zero width channels left in sum to get potential scattering correct
            #   * Gets good (cs from U) to (cs from Fudge) comparison; good (cs from U) to (cs from BB) 
            if multipleSScheme == 'NJOY':  
                channelDict[ ChannelDesignator( L, J, 'elastic', iR, allowedSByThisJl[0], gfact, None, None, False, False, False )] = [ ( ER, GN ) ]
                channelDict[ ChannelDesignator( L, J, 'capture', iR, allowedSByThisJl[0], gfact, None, None, False, False, False )] = [ ( ER, GG ) ]
                if GF != 0.0: channelDict[ ChannelDesignator( L, J, 'fission', iR, allowedSByThisJl[0], gfact, None, None, False, False, False )]     = [ ( ER, GF ) ]
                if GX != 0.0: channelDict[ ChannelDesignator( L, J, 'competitive', iR, allowedSByThisJl[0], gfact, None, None, False, False, False )] = [ ( ER, GX ) ]
                for SS in allowedSByThisJl[1:]:
                    channelDict[ ChannelDesignator( L, J, 'elastic', iR, SS, gfact, None, None, False, False, False )] = [ ( ER, 0.0 ) ]
                    channelDict[ ChannelDesignator( L, J, 'capture', iR, SS, gfact, None, None, False, False, False )] = [ ( ER, 0.0 ) ]
                    if GF != 0.0: channelDict[ ChannelDesignator( L, J, 'fission', iR, SS, gfact, None, None, False, False, False )]     = [ ( ER, 0.0 ) ]
                    if GX != 0.0: channelDict[ ChannelDesignator( L, J, 'competitive', iR, SS, gfact, None, None, False, False, False )] = [ ( ER, 0.0 ) ]
            
            # Way hinted at by ENDF manual: 
            #   * Width divided equally between channels with valid S's
            #   * All channels left in sum to get potential scattering correct
            #   * Gets poor (cs from U) to (cs from Fudge) comparison; poor (cs from U) to (cs from BB) 
            elif multipleSScheme == 'ENDF':   
                for SS in allowedSByThisJl:
                    channelDict[ ChannelDesignator( L, J, 'elastic', iR, SS, gfact, None, None, False, False, False )] = [ ( ER, GN/numAllowedSByThisJl ) ]
                    channelDict[ ChannelDesignator( L, J, 'capture', iR, SS, gfact, None, None, False, False, False )] = [ ( ER, GG/numAllowedSByThisJl ) ]
                    if GF != 0.0: channelDict[ ChannelDesignator( L, J, 'fission', iR, SS, gfact, None, None, False, False, False )]     = [ ( ER, GF/numAllowedSByThisJl ) ]
                    if GX != 0.0: channelDict[ ChannelDesignator( L, J, 'competitive', iR, SS, gfact, None, None, False, False, False )] = [ ( ER, GX/numAllowedSByThisJl ) ]
            
            # Ignore problem: 
            #   * Ignore spin of channels and the fact may be multiple valid spins
            #   * Gets best (cs from U) to (cs from Fudge) comparison; poor (cs from U) to (cs from BB) 
            else :  
                channelDict[ ChannelDesignator( L, J, 'elastic', iR, S, gfact, None, None, False, False, False )] = [ ( ER, GN ) ]
                channelDict[ ChannelDesignator( L, J, 'capture', iR, S, gfact, None, None, False, False, False )] = [ ( ER, GG ) ]
                if GF != 0.0: channelDict[ ChannelDesignator( L, J, 'fission', iR, S, gfact, None, None, False, False, False )]     = [ ( ER, GF ) ]
                if GX != 0.0: channelDict[ ChannelDesignator( L, J, 'competitive', iR, S, gfact, None, None, False, False, False )] = [ ( ER, GX ) ]
            channelDicts.append( channelDict )
        
        self.channels = channelDicts
        
    def getScatteringMatrixU( self, Ein, useTabulatedScatteringRadius = True ):
        '''
        Compute the scattering matrix
        
        Note, unlike the getScatteringMatrixU() function in other resonance classes, 
        the fact that different resonances are entirely different reactions means that the 
        channelDicts have to have an additional layer of sorting that corresponds to the SLBW "level".
        '''
        U = collections.OrderedDict( \
            [ ( iL, numpy.diag(len(cs)*[complex(1.0)]) ) for iL,cs in enumerate(self.channels)])
        rho = self.rho(Ein)

        # For calculating phi, default is to use tabulated scattering radius:
        if useTabulatedScatteringRadius:
            if self.RR.scatteringRadius.isEnergyDependent():
                rhohat = self.RR.scatteringRadius.getValueAs('10*fm', Ein) * self.k(Ein)
            else:
                rhohat = self.RR.scatteringRadius.getValueAs('10*fm') * self.k(Ein)
        else: rhohat = self.rho(Ein)

        for iL,channels in enumerate(self.channels):
            ER = self._energies[iL]
            ERp = ER
            GamTot = 0.0
            eiphis = [ ]
            sqrtgams = [ ] # all are real
            for c in channels:
                if c.reaction == 'elastic': 
                    Gam = channels[c][0][1] * \
                        self.penetrationFactor( c.l, self.rho( Ein,       c.l ) ) / \
                        self.penetrationFactor( c.l, self.rho( abs( ER ), c.l ) )
                    ERp += ( \
                            self.shiftFactor( c.l, self.rho( abs( ER ), c.l ) ) - \
                            self.shiftFactor( c.l, self.rho( Ein,       c.l ) ) \
                        ) * channels[c][0][1] / (2.0 * self.penetrationFactor( c.l, self.rho( abs(ER), c.l) ) )
                    eiphis.append( numpy.exp( complex( 0.0, -self.phi( c.l, rhohat ) ) ) )
                else:
                    Gam = channels[c][0][1]
                    eiphis.append( complex( 1.0, 0.0 ) )
                sqrtgams.append( numpy.sqrt( Gam ) )
                GamTot += Gam
            for i1, c1 in enumerate(channels):
                denominator = complex( ERp - Ein, -GamTot/2.0 ) # c1 & c2 have same l
                for i2, c2 in enumerate(channels):
                    if c2.index != c1.index: continue
                    U[iL][i1][i2] += complex( 0.0, sqrtgams[i1]*sqrtgams[i2] ) / denominator
                    U[iL][i1][i2] *= eiphis[i1] * eiphis[i2]
        return U
        
    def getAngularDistribution(self, E, keepL0Term=True, renormalize=True, reaction1='elastic', reaction2='elastic' ):
        from numericalFunctions import angularMomentumCoupling as nf_amc
        k = self.k(E)                                           # Relative momentum, comes out in b^-1/2
        projSpin = 0.5                                          # Assuming that we have a neutron as a projectile, projSpin is 1/2
        targSpin = self.spin                                    # Target spin from baseclass
        prefactor = 1.0/k/k/(2*projSpin+1)/(2*targSpin+1)       # Prefactor for cross section, gets units, main energy dependence, and spin dep. prefactor
        B = []                                                  # List of Blatt-Biedenharn coefficients, used in Legendre series expansion of the angular distribution
        relTol = 1e-5                                           # Relative tolerance used to decide when to cut off Legendre series
        Lmax = min( max( [ xL['L'] for xL in self.Ls ] ), 10 )  # Lmax also used to decide when to cut off Legendre series
        allowedSpins = getAllowedTotalSpins( targSpin, projSpin, False )    # Possible values of the total channel spin, based on the projectile spin (projSpin) and target spin (targSpin)
        U = self.getScatteringMatrixU( E )                      # The scattering matrix
        L = 0
        while len(B)<=2 or ( B[-1]>=relTol*B[0] and L<=Lmax ):
            B.append(0.0)
            for iLev, channels in enumerate( self.channels ):
                for S in allowedSpins:
                    for Sp in allowedSpins:
                        spinFactor =  pow( -1, S - Sp ) / 4.0
                        for i1, c1 in enumerate(channels):
                            if c1.reaction != reaction1: continue
                            if c1.s != None and c1.s != S: continue
                            for i2, c2 in enumerate(channels):
                                if c2.reaction != reaction2: continue
                                if c2.s != None and c2.s != S: continue
                                Z = nf_amc.z_coefficient( int(2.*c1.l), int(2.*c1.J), int(2.*c2.l), int(2.*c2.J), int(2.*S), int(2.*L) )
                                if Z == 0.0: continue
                                for i1p, c1p in enumerate(channels):
                                    if c1p.reaction != reaction1: continue
                                    if c1.J != c1p.J: continue
                                    if c1p.s != None and c1p.s != Sp: continue
                                    for i2p, c2p in enumerate(channels):
                                        if c2p.reaction != reaction2: continue
                                        if c2.J != c2p.J: continue
                                        if c2p.s != None and c2p.s != Sp: continue
                                        Zp = nf_amc.z_coefficient( int(2*c1p.l), int(2*c1p.J), int(2*c2p.l), int(2*c2p.J), int(2*Sp), int(2*L) )
                                        if Zp == 0.0: continue
                                        thisB = spinFactor * Z * Zp * ( ( (c1==c1p) - U[iLev][i1][i1p] ).conj() * ( (c2==c2p) - U[iLev][i2][i2p] ) ).real
                                        B[-1] += thisB
            L += 1
        return [ prefactor * x for x in B ]

    @blockwise
    def getCrossSection(self, E):
        captureSum = 0
        elasticSum = 0
        fissionSum = 0
        rho = self.rho(E)
        # for calculating phi, always use tabulated scattering radius:
        rhohat = self.RR.scatteringRadius.getValueAs('10*fm') * self.k(E)
        for L in self.Ls:
            l = L['L']
            phi = self.phi(l,rhohat)
            P = self.penetrationFactor(l,rho)
            S = 0.5*self.shiftFactor(l,rho)
            for J in L['Js']:
                gfactor = J['gfact']
                for spin in J['channelSpins']:
                    dE = (E-(spin['energy']+spin['neutronWidth']*(spin['shiftFactor']-S)))
                    totalWidth = P*spin['neutronWidth'] + spin['captureWidth'] + spin['fissionWidthA']
                    denominator = dE**2 + totalWidth**2 / 4
                    captureSum += gfactor * numpy.sum( ( P * spin['neutronWidth'] * spin['captureWidth'] ) / denominator , axis=1)
                    fissionSum += gfactor * numpy.sum( ( P * spin['neutronWidth'] * spin['fissionWidthA'] ) / denominator , axis=1)
                    elasticSum += gfactor * numpy.sum( ( P*spin['neutronWidth'] * ( P*spin['neutronWidth'] 
                        - 2 * numpy.sin(phi)**2 * totalWidth + 4 * numpy.sin(phi)*numpy.cos(phi) * (E-spin['energy']) ) ) 
                        / denominator , axis=1)
            # numpy.sum(..., axis=1) returns row vector, so also convert first term to row vector:
            elasticSum += 4*(2*l+1)*numpy.sin(phi)[:,0]**2
        
        # get common factor 'beta' as a row vector:
        beta = numpy.pi / self.k(E)[:,0]**2
        capture = beta * captureSum
        elastic = beta * elasticSum
        fission = beta * fissionSum
        for reaction in (capture,elastic,fission):
            reaction[ reaction<=0 ] = 0
        total = elastic + capture + fission
        result = {'total':total, 'elastic':elastic, 'capture':capture, 'fission':fission}
        return result


#### Multi-level Breit-Wigner ###
class MLBWcrossSection(RRBaseClass):
    """
    Given a resonance region in MLBW format,
    create a class with all data required to reconstruct
    
    Only the elastic channel differs from SLBW
    """
    def __init__(self, reactionSuite, sectionIndex=None, enableAngDists=False, verbose=True):
        super(MLBWcrossSection, self).__init__(reactionSuite, sectionIndex)
        self.verbose = verbose
        self.sortLandJ()
        if self.verbose:
            print ("From %f to %f eV, reconstructing using Multi-level Breit-Wigner" % 
                    (self.lowerBound,self.upperBound))
        if enableAngDists: self.setResonanceParametersByChannel()
    
    def getChannelConstantsBc( self ):
        '''
        For ENDF's MLBW, should be 
        
            ..math::
                B_c = S_\ell(|E_\lambda|)
                
        where $ell$ is the channel angular momentum and $\lambda$ is the resonances index for the channel
        '''
        raise NotImplementedError( "write me" )
        
    def getScatteringMatrixU( self, Ein, useTabulatedScatteringRadius = True ):
        '''
        Compute the scattering matrix.  We could have used the generic U function in the base class, 
        but Froehner has "simplifications" that we took advantage of here (that and I don't know what the 
        R matrix is exactly in the case of MLBW).
        '''
        # Initialize lists and matrices
        sqrtGam = numpy.zeros( ( self.nChannels, self.nResonances ) )
        DEN  = self.nResonances * [ complex( 1.0 ) ]
        U = numpy.identity( self.nChannels, dtype=complex )
        Gtots = self.nResonances * [ 0.0 ]
        ERp = [ ER for ER in self._energies ]
        eiphis = self.getEiPhis(  Ein, useTabulatedScatteringRadius = useTabulatedScatteringRadius )
        
        # Precompute sqrt(Gamma) and sum up to get the total Gamma
        # Total Gamma may be different from what is given in the ENDF data file
        # Also, work out energy shifts
        for ic, c in enumerate( self.channels ):
            if c.reaction == 'elastic': 
                for iR, G in self.channels[c].items():
                    ER = self._energies[iR] 
                    Gam = G * \
                        self.penetrationFactor( c.l, self.rho( Ein,       c.l ) ) / \
                        self.penetrationFactor( c.l, self.rho( abs( ER ), c.l ) )
                    ERp[iR] += ( \
                            self.shiftFactor( c.l, self.rho( abs( ER ), c.l ) ) - \
                            self.shiftFactor( c.l, self.rho( Ein,       c.l ) ) \
                        ) * G / 2.0 / self.penetrationFactor( c.l, self.rho( abs(ER), c.l) ) 
                    sqrtGam[ ic, iR ] = numpy.sqrt(Gam)
                    Gtots[ iR ] += Gam                    
            else:
                for iR, G in self.channels[c].items(): 
                    ER = self._energies[iR] 
                    sqrtGam[ ic, iR ] = numpy.sqrt(G)
                    Gtots[ iR ] += G
                
        # If the width in the ENDF file is bigger than the sum of parts, there is a competitive width
        # otherwise, the ENDF value is wrong, use the one we get from the sum of parts
        for iR in range( self.nResonances ): Gtots[ iR ] =  max( Gtots[ iR ], self._widths[ iR ] )
        
        # Compute the denominator 
        for iR,Gtot in enumerate( Gtots ): DEN[ iR ] = complex( ERp[ iR ] - Ein, -Gtot/2.0 ) 

        # Compute U itself.  Can we accelerate this with numpy routines?
        for i1 in range( self.nChannels ):
            for i2 in range( self.nChannels ):
                for iR in range( self.nResonances ): U[i1,i2] += complex( 0.0, sqrtGam[i1,iR] * sqrtGam[i2,iR] ) / DEN[ iR ]
                U[i1,i2] *= eiphis[i1] * eiphis[i2]
        return U
        
    @blockwise
    def getCrossSection(self, E):
        captureSum = 0
        elasticSum = 0
        fissionSum = 0
        rho = self.rho(E)
        # for phi:
        if self.RR.scatteringRadius.isEnergyDependent():
            rhohat = numpy.array(self.RR.scatteringRadius.getValueAs('10*fm', E[:,0]))[:,numpy.newaxis] * self.k(E)
        else:
            rhohat = self.RR.scatteringRadius.getValueAs('10*fm') * self.k(E)
        for L in self.Ls:
            l = L['L']
            phi = self.phi(l,rhohat)
            P = self.penetrationFactor(l,rho)
            S = 0.5*self.shiftFactor(l,rho)
            for J in L['Js']:
                gfactor = J['gfact']
                for spin in J['channelSpins']:
                    dE = (E-(spin['energy']+spin['neutronWidth']*(spin['shiftFactor']-S)))
                    totalWidth = P*spin['neutronWidth'] + spin['captureWidth'] + spin['fissionWidthA']
                    denominator = dE**2 + totalWidth**2 / 4
                    commonFactor = P * spin['neutronWidth'] / denominator
                    captureSum += gfactor * numpy.sum( commonFactor * spin['captureWidth'] , axis=1)
                    fissionSum += gfactor * numpy.sum( commonFactor * spin['fissionWidthA'] , axis=1)

                    # simple elastic method, Eq D.19 - D.21 in ENDF 102. This has numerical issues, however:
                    # U_nn = numpy.exp(-2*1j*phi[:,0]) * (1 + numpy.sum( 1j*P*spin['neutronWidth'] /
                    #     (spin['energy']-E-1j*totalWidth/2) , axis=1))
                    # elasticSum += gfactor * abs( (1-U_nn)**2 )

                    # Instead of the above, use the following from RECENT:
                    elasticTerm1 = numpy.sum( totalWidth/2 * commonFactor, axis=1 )
                    elasticTerm2 = numpy.sum( dE * commonFactor, axis=1 )
                    sin2ps = 2*numpy.sin(phi[:,0])*numpy.cos(phi[:,0])
                    sinps2 = 2*numpy.sin(phi[:,0])**2
                    elasticSum += gfactor * ((sinps2-elasticTerm1)**2 + (sin2ps+elasticTerm2)**2)
            if True: #addMissingGfactor:
                elasticSum += 2*numpy.sin(phi[:,0])**2 * 2*self.missingGfactor[l]
        
        # get common factor 'beta' as a row vector:
        beta = numpy.pi / self.k(E)[:,0]**2
        capture = beta * captureSum
        elastic = beta * elasticSum
        fission = beta * fissionSum
        for reaction in (capture,elastic,fission):
            reaction[ reaction<=0 ] = 0
        total = elastic + capture + fission
        return {'total':total, 'elastic':elastic, 'capture':capture, 'fission':fission}


###### Reich_Moore and R-Matrix Limited ######

# some helper functions:

def getRI_SI(E, Eres, captureWidth, widths, penetrabilities):
    """
    Both versions of Reich_Moore formalisms (LRF=3 and 7) rely on building matrices R and S,
    which represent symmetric and anti-symmetric scattering respectively
    
    matrix elements R[i,j] and S[i,j] =
    
    (summed over resonances) partialWidth[i]*partialWidth[j] * coefficient/(dE**2+captureWidth**2),
    
    for the ith/jth open channel. For S, the coefficient is 'dE', for R 'captureWidth'
    and partialWidth[i] = widths[i] * penetrabilities[i]
    
        ( widths[i] is a row vector of resonance widths, and penetrabilities[i] is a column vector 
        of the penetrability for each incident energy )
    
    Then, invert to find RI and SI such that (I+R+jS)*(I+RI+jSI) = I
    where I is the identity matrix, and j = sqrt(-1)
    
    Incident energy dependence appears in both in E and the penetrabilities
    
    Additional documentation from RECENT::
    
      THE CROSS SECTIONS ARE DEFINED TO BE,
 
      TOTAL        =2*GJ*REAL(I - U(N,N))
      ABSORPTION   =4*GJ*(REAL(RHO(N,N)) - RHO(N,N)**2)
                   =  GJ*(I - U(N,N)**2)
      ELASTIC      =  GJ*(I - U(N,N))**2
      FISSION      =4*GJ*(SUM OVER C)(RHO(N,C)**2)
      CAPTURE      = ABSORPTION - FISSION
 
      WHICH ARE COMPLETELY DEFINED IN TERMS OF U(N,N) AND RHO(N,C),
 
      RHO(N,C)     =I - INVERSE(I - K)
      U(N,N)       =EXP(-I*2*PS)*(2*INVERSE(I - K) - I)
                   =(COS(2*PS) - I*SIN(2*PS))*(2*INVERSE(I - K) - I)
                   =(COS(2*PS) - I*SIN(2*PS))*(1 - 2*RHO(N,N))
 
                   =COS(2*PS)*(I-2*REAL(RHO)) - I*2*SIN(2*PS)*IM(RHO)
 
      ... (cmattoon: cutting some of the documentation here) ...
      
      (I - K)      = (R + I) - I*S
 
      R            = SQRT(GAM(C)/2*GAM(C*)/2)*(GAM/2)/DEN
      S            = SQRT(GAM(C)/2*GAM(C*)/2)*(ER-E)/DEN
      GAM(R)       = ELIMINATED RADIATIVE WIDTH
      GAM(C)       = PARTIAL WIDE FOR CHANNEL C FOR A RESONANCE
      DEN          = ((ER - E)**2 + (GAM/2)**2)
 
      SUMMED OVER RESONANCES FOR EACH (L,J) SEQUENCE.
      
    """
    NE = len(E)         # number of energies
    dim = len(widths)   # dimension of matrix at each energy
    try:
        # use wrapped c version if available:
        import getScatteringMatrices
        # convert any constant width/penetrability data to array w/correct dimensions:
        nRes = len(Eres)
        for i in range(dim):
            if len(widths[i]) != nRes:
                widths[i] = numpy.ones((1,nRes)) * widths[i][0,0]
            if len(penetrabilities[i]) != NE:
                penetrabilities[i] = numpy.ones((NE,1)) * penetrabilities[i][0,0]
        R, S = getScatteringMatrices.getScatteringMatrices( E, Eres.copy(), captureWidth,
                numpy.array(widths), numpy.array(penetrabilities) )
    except ImportError:
        # can't import c version, so use numpy (memory-intensive)
        # define some common factors.
        # These are all 2d arrays with len(E) rows, len(Eres) columns:
        dE = Eres - E
        DEN = 1/(dE**2 + captureWidth**2)
        captOverDEN = captureWidth*DEN
        dEoverDEN = dE*DEN
        del dE, DEN

        R = numpy.zeros((NE,dim,dim))
        S = numpy.zeros((NE,dim,dim))
        for i in range(dim):
            for j in range(i,dim):
                width_ij = widths[i]*widths[j]
                p_ij = penetrabilities[i]*penetrabilities[j]
                R[:,i,j] = p_ij.T * numpy.sum( width_ij * captOverDEN , axis=1)
                S[:,i,j] = p_ij.T * numpy.sum( width_ij * dEoverDEN , axis=1)
                # symmetrize:
                R[:,j,i] = R[:,i,j]
                S[:,j,i] = S[:,i,j]
        del captOverDEN, dEoverDEN

    # now we have complex matrix I+R+j*S
    # invert to find RI and SI such that (I+R+j*S) * (I+RI+j*SI) = I
    # see comments for subroutine FROBNS3 in recent.f:
    if dim==1:
        # only have neutron width (and 'eliminated' capture width).
        # Can obtain quicker solution by manually inverting:
        R,S = R.squeeze(),S.squeeze()   # cut extra dimensions
        # invert:
        DET = (R+1)**2 + S**2
        SI = -S/DET
        RI = -(R*(R+1)+S**2)/DET
        return RI,SI
    else:
        # have competition and/or fission, must invert matrices:
        def dotproduct(*args):
            # helper function: dotproduct(A,B,C) == numpy.dot(A, numpy.dot(B,C))
            def quickdot(a,b):
                # take dot product of two Nx(MxM) arrays (containing MxM matrices):
                result = numpy.zeros( a.shape )
                dim = a.shape[1]
                for i in range(dim):
                    for j in range(dim):
                        result[:,i,j] = numpy.sum( a[:,i,:] * b[:,:,j], axis=1 )
                return result
            return reduce( quickdot, args )
        
        def vectorInv(arr):
            # invert all MxM matrices in Nx(MxM) array
            dim = arr.shape[1]
            if dim==2:
                arrinv = numpy.zeros(arr.shape)
                det = (arr[:,0,0]*arr[:,1,1]-arr[:,0,1]*arr[:,1,0])
                arrinv[:,0,0] = arr[:,1,1]/det
                arrinv[:,0,1] = -arr[:,0,1]/det
                arrinv[:,1,0] = -arr[:,1,0]/det
                arrinv[:,1,1] = arr[:,0,0]/det
                return arrinv
            elif dim==3:
                arrinv = numpy.zeros(arr.shape)
                det = (arr[:,0,0]*(arr[:,1,1]*arr[:,2,2]-arr[:,1,2]*arr[:,2,1])
                    -arr[:,0,1]*(arr[:,1,0]*arr[:,2,2]-arr[:,1,2]*arr[:,2,0])
                    +arr[:,0,2]*(arr[:,1,0]*arr[:,2,1]-arr[:,1,1]*arr[:,2,0]) )
                arrinv[:,0,0] = (arr[:,1,1]*arr[:,2,2]-arr[:,1,2]*arr[:,2,1])/det
                arrinv[:,0,1] = (arr[:,0,2]*arr[:,2,1]-arr[:,1,0]*arr[:,2,2])/det
                arrinv[:,0,2] = (arr[:,0,1]*arr[:,1,2]-arr[:,0,2]*arr[:,1,1])/det
                arrinv[:,1,0] = (arr[:,1,2]*arr[:,2,0]-arr[:,1,0]*arr[:,2,2])/det
                arrinv[:,1,1] = (arr[:,0,0]*arr[:,2,2]-arr[:,0,2]*arr[:,2,0])/det
                arrinv[:,1,2] = (arr[:,0,2]*arr[:,1,0]-arr[:,0,0]*arr[:,1,2])/det
                arrinv[:,2,0] = (arr[:,1,0]*arr[:,2,1]-arr[:,1,1]*arr[:,2,0])/det
                arrinv[:,2,1] = (arr[:,0,1]*arr[:,2,0]-arr[:,0,0]*arr[:,2,1])/det
                arrinv[:,2,2] = (arr[:,0,0]*arr[:,1,1]-arr[:,0,1]*arr[:,1,0])/det
                return arrinv
            else:
                result = numpy.zeros( arr.shape )
                NE = arr.shape[0]
                for i in range(NE):
                    result[i] = numpy.linalg.inv(arr[i])
                return result
        
        identity = numpy.zeros( (NE, dim, dim) )
        for i in range(dim):
            identity[:,i,i] = 1
        
        invRplusI = vectorInv( R + identity )
        SRS = dotproduct( S, invRplusI, S )
        RI = dotproduct( -vectorInv( R+identity + SRS ), R+SRS )
        SI = dotproduct( -invRplusI, S, RI+identity )
        return RI, SI


#### Reich_Moore ####
class RMcrossSection(RRBaseClass):
    """
    simplified Reich_Moore (LRF=3 in ENDF)
    More complex than Breit-Wigner approximations, but we can still use same __init__
    """
    def __init__(self, reactionSuite, sectionIndex=None, enableAngDists=False, verbose=True):
        super(RMcrossSection, self).__init__(reactionSuite, sectionIndex)
        self.verbose = verbose
        self.sortLandJ()
        if self.verbose:
            print ("From %f to %f eV, reconstructing using Reich_Moore (LRF=3)" % 
                    (self.lowerBound,self.upperBound))
        if enableAngDists: self.setResonanceParametersByChannel( useReichMooreApproximation = True )

    def getChannelConstantsBc( self ):
        '''
        For ENDF's Reich-Moore, should be 
        
            ..math::
                B_c = -\ell
                
        where $ell$ is the channel angular momentum
        
        what is not clearly stated in the ENDF manual is that Bc = 0 for this case
        '''
        return [ 0.0 for c in self.channels ]
        
    def getL0Matrix( self, Ein, Bc=None ):
        '''
        Get the L0 matrix of Froehner:
        
            ..math::
                
                {\bf L^0}_{cc'} = \delta_{cc'} (L_c-B_c) 

        where 
    
            ..math::
                
                L_c = S_c + i P_c
                
        But.... ENDF's RM formulation uses a shift of 0 and a Bc of 0!!!

        '''
        L0 = []
        for ic,c in enumerate(self.channels):
            rho = self.rho( Ein, c.l )
            L0.append( complex( 0.0, self.penetrationFactor( c.l, rho ) )  )
        return numpy.diag( L0 )

    def getRMatrix( self, Ein ): 
        '''
        The R matrix in the Reich-Moore approximation:
        
            ..math::
            
                R_{cc'} = \sum_\lambda \frac{\gamma_{\lambda c}\gamma_{\lambda c'}{E_\lambda - E - i\Gamma_{\lambda\gamma}/2}
        '''
        R = numpy.zeros( ( self.nChannels, self.nChannels), dtype=complex )
        
        # Loop through all resonances
        for iR, ER in enumerate( self._energies ):
        
            # Extract the gamma width for the first gamma channel that has this resonance
            gamWidth = 0.0
            for cg in self.eliminatedChannels:
                if iR in self.eliminatedChannels[ cg ]:
                    gamWidth = self.eliminatedChannels[ cg ][ iR ]
                    if gamWidth != 0.0: break
            
            # Precompute  the reduced widths
            redWidth = []
            for c in self.channels:
                if iR in self.channels[ c ]: 
                    if c.reaction == 'elastic': pen = numpy.copysign( self.penetrationFactor( c.l, self.rho( abs( ER ), c.l ) ), ER )
                    else:                       pen = numpy.copysign( 1.0, ER )
                    redWidth.append( numpy.copysign( numpy.sqrt( self.channels[ c ][ iR ] / ( 2.0*abs( pen ) ) ), pen ) )
                else: redWidth.append( 0.0 )
                                    
            # Loop through all channels to accumulate the R Matrix elements
            for ic1, c1 in enumerate( self.channels ): 
                if not iR in self.channels[ c1 ]: continue
                for ic2, c2 in enumerate( self.channels ):
                    if not iR in self.channels[ c2 ]: continue
                    R[ ic1, ic2 ] += redWidth[ ic1 ] * redWidth[ ic2 ] / complex( ER - Ein, -gamWidth/2.0 ) 
                                
        return R
            
    @blockwise
    def getCrossSection(self, E):
        elasticSum = 0
        absorbtionSum = 0
        fissionSum = 0
        haveFission = self.RR.resonanceParameters.getColumn('fissionWidthA',units='eV') is not None
        # for calculating phi, always use tabulated scattering radius:
        for L in self.Ls:
            l = L['L']
            rho = self.rho(E,l)
            rhohat = self.RR.getScatteringRadius(l).getValueAs('10*fm') * self.k(E)
            phi = self.phi(l, rhohat)
            P = numpy.sqrt( self.penetrationFactor(l,rho) )

            for J in L['Js']:
                gfactor = J['gfact']
                for spin in J['channelSpins']:
                    captureWidth = spin['captureWidth'] * 0.5
                    widths = []
                    penetrabilities = []
                    widths.append( numpy.sqrt( spin['neutronWidth'] * 0.5 ) )
                    penetrabilities.append(P)
                    if haveFission:
                        fissWidthA = numpy.sqrt(numpy.abs(spin['fissionWidthA'] * 0.5))
                        fissWidthB = numpy.sqrt(numpy.abs(spin['fissionWidthB'] * 0.5))
                        fissWidthA[ spin['fissionWidthA']<0 ] *= -1
                        fissWidthB[ spin['fissionWidthB']<0 ] *= -1
                        widths += [fissWidthA, fissWidthB]
                        penetrabilities += [numpy.array([[1]]),numpy.array([[1]])]
                        
                        RI,SI = getRI_SI(E, spin['energy'], captureWidth, widths, penetrabilities)
                        
                        elasticSum += gfactor * ((2*numpy.sin(phi[:,0])**2+2*RI[:,0,0])**2 + 
                                (2*numpy.sin(phi[:,0])*numpy.cos(phi[:,0]) + 2*SI[:,0,0])**2)
                        absorbtionSum += -4*gfactor * (RI[:,0,0]+RI[:,0,0]**2+SI[:,0,0]**2)
                        fissionSum += 4*gfactor * (RI[:,0,1]**2+RI[:,0,2]**2+SI[:,0,1]**2+SI[:,0,2]**2)
                    
                    else: # faster method when we don't have fission:
                        RI, SI = getRI_SI(E, spin['energy'], captureWidth, widths, penetrabilities)
                        elasticSum += gfactor * ((2*numpy.sin(phi[:,0])**2+2*RI)**2 + 
                                (2*numpy.sin(phi[:,0])*numpy.cos(phi[:,0]) + 2*SI)**2)
                        absorbtionSum += -4*gfactor * (RI + RI**2 + SI**2)

            if True: #addMissingGfactor:
                elasticSum += 2*numpy.sin(phi[:,0])**2 * 2*self.missingGfactor[l]

        # get common factor 'beta' as a row vector:
        beta = numpy.pi / self.k(E)[:,0]**2
        elastic = beta * elasticSum
        absorbtion = beta * absorbtionSum
        fission = beta * fissionSum
        capture = absorbtion - fission
        total = elastic + absorbtion

        return {'total':total, 'elastic':elastic, 'capture':capture, 'fission':fission}


#### R-Matrix Limited ####
class RMatrixLimitedcrossSection(RRBaseClass):
    """
    extended Reich_Moore (LRF=7 in ENDF)
    Here, resonances are sorted primarily by J: within each 'spin group', total J is conserved
    One or more competitive channels may be used in this case.
    Also, each resonance may have contributions from multiple l-waves
    """
    def __init__(self, reactionSuite, sectionIndex=None, enableAngDists=False, verbose=True):
        super(RMatrixLimitedcrossSection, self).__init__(reactionSuite, sectionIndex)
        self.verbose = verbose
        if self.verbose:
            print ("From %f to %f eV, reconstructing using R-Matrix Limited (LRF=7)" % 
                    (self.lowerBound,self.upperBound))
        
        # for RML, need info on residual nuclei from each channel:
        for ch in self.RR.openChannels:
            reaction = reactionSuite.getReaction(ch.channel)
            # elastic, capture or competitive?
            if reaction == reactionSuite.getReaction('capture'): ch.tag = 'capture'
            elif reaction == reactionSuite.getReaction('elastic'): ch.tag = 'elastic'
            elif 'fission' in ch.channel: ch.tag = ch.channel
            else: ch.tag = 'competitive'
            
            Q = ch.Qvalue.getValueAs('eV')
            if not Q:
                Q = reaction.getQ('eV')
                # adjust Q value for residual in excited state:
                for particle in reaction.outputChannel.particles:
                    Q -= particle.getLevelAsFloat('eV')
            Q /= self.targetToNeutronMassRatio / (self.targetToNeutronMassRatio + 1.0)
            particles = reaction.outputChannel.particles
            ch.reactionInfo = { 'Q':Q, 'particles': particles }
        
        for sg in self.RR.spinGroups:
            sg.energy = numpy.array( sg.resonanceParameters.getColumn('energy','eV') )
            
            for column in sg.resonanceParameters.columns:
                if column.name=='energy': continue
                channelName = column.name.split(' width')[0]
                openChannel = [ch for ch in self.RR.openChannels if ch.channel==channelName]
                assert len(openChannel)==1, ("Could not find unique openChannel named '%s'" % channelName)
                column.tag = openChannel[0].tag
        
        # for energy grid generation:
        energies, totalWidths = [],[]
        for sg in self.RR.spinGroups:
            energies.extend( sg.resonanceParameters.getColumn('energy','eV') )
            widths = [sg.resonanceParameters.getColumn( col.name, 'eV' ) for col in
                    sg.resonanceParameters.columns if col.name != 'energy']
            totalWidths.extend( numpy.sum( widths, axis=0 ) )
        zipped = sorted(zip(energies,totalWidths))
        self._energies, self._widths = zip(*zipped)

    def rho(self, E):
        # get the channel radius rho for each channel:
        if self.RR.calculateChannelRadius:
            a = [0.123 * self.target.getMass('amu')**(1./3.) + 0.08     # eq D.14 in ENDF manual
                    for i in range(len(self.RR.channels))]
        else:
            a = [ch.scatteringRadius.getValueAs('10*fm')
                    for ch in self.RR.channels]
        return [self.k(E) * b for b in a]
    
    def k_competitive(self, E, pA, pB):
        """ calculate k for any 2-body output channel. Note that if pA and pB are target and neutron,
        this reduces to self.k(E) as defined above in the resonanceReconstructionBaseClass """
        pA_mass, pB_mass = [p.particle.getMass('amu') for p in pA,pB]
        mn = self.projectile.getMass('amu')
        # sqrt(2*neutronMass)/hbar == 2.196807 (eV*barn)**-1/2. Thus for energy in eV, k is in b**-1/2
        return 2.196807122623e-3 * numpy.sqrt( pB_mass**2 * pA_mass / (pB_mass + mn) /
                (pB_mass + pA_mass) / mn ) * numpy.sqrt(E)
    
    def eta(self, E, pA, pB):
        '''for competitive channels with 2 charged particles,
        parameter eta used to find penetrability. eq D.79 in ENDF manual:'''
        pA_z, pB_z = [p.particle.getZ_A_SuffixAndZA()[0] for p in pA, pB]
        if pA_z * pB_z == 0: return numpy.array([0])
        pA_mass, pB_mass = [p.particle.getMass('amu') for p in pA,pB]
        mn = self.projectile.getMass('amu')
        eSq_over_hbarSq = 3.4746085579272e-1    # units?
        eta = (eSq_over_hbarSq * pA_z * pB_z * (pA_mass/mn) * (pB_mass/(pB_mass+pA_mass)) 
                / self.k_competitive(E, pA, pB))
        return eta
    
    def getChannelConstantsBc( self ):
        '''
        For ENDF's Reich-Moore, should be 
        
            ..math::
                B_c = -\ell
                
        where $ell$ is the channel angular momentum
        '''
        return [ -c.l for c in self.channels ] # what ENDF manual says

    def coulombPenetrationFactor(self, L, rho, eta):
        """ for competitive channels with 2 charged particles, 
        calculate coulomb wavefunctions and penetrability.
        Here we use an external subroutine 'coulfg2' from Thompson et al, converted to c and wrapped """
        import getCoulombWavefunctions
        if any(rho<0.02) and self.verbose:
            print ("      WARNING: small rho values encountered, may need modified coulomb penetrability!")
        penetrability = numpy.zeros( rho.shape )
        F,G = getCoulombWavefunctions.getCoulombWavefunctions( rho, eta, L )
        penetrability = rho / (F**2 + G**2)
        return penetrability
    
    @blockwise
    def getCrossSection(self, E):
        """
        similar to Reich_Moore (LRF=3). Refer to comments there, and in
        subroutine SIGRM2 in recent.f
        """
        elasticSum = 0
        absorbtionSum = 0
        fissionSum = 0
        nCompetitive = len( [ch for ch in self.RR.openChannels if ch.tag=='competitive'] )
        competitiveSum = [0,] * nCompetitive
        haveFission = any(['fission' in a.channel for a in self.RR.openChannels])
        
        for sg in self.RR.spinGroups:
            j = sg.spin.value
            gfact = (2*j+1.0)/(2*(2*self.spin+1))
            
            widths = []
            penetrabilities = []
            phis = []    # phi is l-dependant
            qvalues = []
            for chan in sg.resonanceParameters.columns:
                if chan.name=='energy': continue
                chanWidths = numpy.array( sg.resonanceParameters.getColumn( chan.name, 'eV' ) )
                channelName = chan.name.split(' width')[0]
                parentChannel = [ch for ch in self.RR.openChannels if ch.channel == channelName][0]
                
                if chan.tag == 'capture':
                    captureWidth = chanWidths / 2
                    continue
                # else:
                l = chan['L']
                pA,pB = parentChannel.reactionInfo['particles']
                Qvalue = parentChannel.reactionInfo['Q']
                qvalues.append( Qvalue )
                
                # channel radius:
                def rho( E ):
                    if self.RR.calculateChannelRadius:
                        a = 0.123 * self.target.getMass('amu')**(1./3.) + 0.08
                    elif chan.attributes.get('effectiveRadius') is not None:
                        a = chan['effectiveRadius'].getValueAs('10*fm')
                    else:
                        a = self.RR.scatteringRadius.getValueAs('10*fm')
                    return self.k_competitive(E, pA, pB) * a
                
                # for calculating phi, always use tabulated (not calculated) scattering radius:
                if chan.attributes.get('effectiveRadius') is not None:
                    rhohat = chan['effectiveRadius'].getValueAs('10*fm') * self.k(E)
                else:
                    rhohat = self.RR.scatteringRadius.getValueAs('10*fm') * self.k(E)
                phinow = self.phi(l, rhohat)
                phinow[ phinow/rhohat < 1e-6 ] = 0  # see subroutine 'facphi' in RECENT
                phis.append( phinow )
                
                # penetrability:
                Ex1 = abs(sg.energy)+Qvalue     # evaluated at resonances
                Ex2 = E+Qvalue; Ex2[Ex2<0] = 0  # evaluated at each incident energy
                eta1 = self.eta(Ex1, pA, pB)
                eta2 = self.eta(Ex2, pA, pB)
                if any(eta1 > 0):
                    # channel has two charged particles, need Coulomb penetrability:
                    penetrabilityAtResonances = self.coulombPenetrationFactor(l, rho(abs(Ex1)), eta1)
                    penetrabilityAtEin = numpy.sqrt( self.coulombPenetrationFactor(l,rho(Ex2), eta2) )
                else:
                    # no coulomb contribution:
                    penetrabilityAtResonances = self.penetrationFactor(l, rho(abs(Ex1)))
                    penetrabilityAtEin = numpy.sqrt( self.penetrationFactor(l, rho(Ex2)) )
                
                reducedWidths = numpy.sqrt( abs(chanWidths) / (2 * penetrabilityAtResonances) )
                reducedWidths[ chanWidths<0 ] *= -1
                widths.append( reducedWidths )
                penetrabilities.append( penetrabilityAtEin )
            
            # are there any threshold reactions (negative Q-values)?
            # If so, we must break up the calculation above/below each threshold
            if any( [q<0 for q in qvalues] ):
                # sanity check: channels should be sorted by increasing threshold
                assert qvalues==sorted(qvalues)[::-1]
                
                RI = numpy.zeros((len(E), len(widths), len(widths)))
                SI = numpy.zeros( RI.shape )
                
                # get energy index for each threshold:
                import bisect
                thresholds = ([bisect.bisect(E+Q,0) for Q in qvalues])
                thresholds.append( len(E) )
                threshSet = sorted(set(thresholds))
                for i in range(len(threshSet)-1):
                    low = threshSet[i]
                    high = threshSet[i+1]
                    nOpen = len( [th for th in thresholds if th-low <= 0] )
                    RI_now, SI_now = getRI_SI( E[low:high], sg.energy, captureWidth, 
                            widths[:nOpen], [p[low:high] for p in penetrabilities[:nOpen]] )
                    if nOpen==1:  # numpy insists on same number of dimensions for copy:
                        RI_now.shape = SI_now.shape = (len(RI_now),1,1)
                    
                    RI[low:high, :nOpen, :nOpen] = RI_now
                    SI[low:high, :nOpen, :nOpen] = SI_now
            
            else:   # no threshold reactions
                RI,SI = getRI_SI(E, sg.energy, captureWidth, widths, penetrabilities)
            
            # reconstruct:
            chanIds = []
            for col in sg.resonanceParameters.columns:
                if col.name=='energy' or col.tag=='capture': continue
                thisChannel = [ch for ch in self.RR.openChannels if 
                        ch.channel==col.name.split(' width')[0] ][0]
                chanIds.append( self.RR.openChannels.index( thisChannel ) )
            chanIds = numpy.array( chanIds )
            elasID = [tmp.index for tmp in self.RR.openChannels if tmp.tag=='elastic'][0] # should be '1'
            competitiveIDs = [tmp.index for tmp in self.RR.openChannels if tmp.tag=='competitive']
            for chan in self.RR.openChannels:
                if chan.tag == 'capture':
                    continue    # 'eliminated' channel in Reich_Moore
                # which matrix elements correspond to this channel?
                thisChanIds = numpy.where( chanIds==chan.index )[0]
                if len(thisChanIds)==0:
                    continue
                id1, id2 = min(thisChanIds), max(thisChanIds)+1
                assert len(thisChanIds) == id2-id1
                if chan.tag == 'elastic':
                    elas = 0; absorb = 0
                    # diagonal component
                    for i in range(id1,id2):
                        phi = phis[i] [:,0] # convert to row vector
                        sinsqr = numpy.sin(phi)**2
                        sincos = numpy.sin(phi)*numpy.cos(phi)
                        elas += ( sinsqr * (sinsqr+2*RI[:,i,i]) + sincos * (sincos + 2*SI[:,i,i]) )
                        absorb += -RI[:,i,i]
                    # add cross-terms:
                    elas += numpy.sum(numpy.sum( RI[:,id1:id2,id1:id2]**2 + SI[:,id1:id2,id1:id2]**2,
                            axis=1), axis=1)
                    absorb -= numpy.sum(numpy.sum( RI[:,id1:id2,id1:id2]**2 + SI[:,id1:id2,id1:id2]**2, 
                            axis=1), axis=1)
                    elasticSum += gfact * 4 * elas
                    absorbtionSum += gfact * 4 * absorb
                
                else: # competitive or fission:
                    # sum of cross-terms between this channel and elastic:
                    if 'fission' in chan.tag:
                        print ("      WARNING: haven't dealt with R-Matrix fission yet... proceed with caution!")
                    elasticIds = numpy.where( chanIds==elasID )[0]
                    el1, el2 = min(elasticIds), max(elasticIds)+1
                    comp = numpy.sum(numpy.sum(SI[:,id1:id2,el1:el2]**2 + RI[:,id1:id2,el1:el2]**2,
                            axis=1),axis=1)
                    # may have more than one competitive: add to correct channel
                    competitiveSum[ competitiveIDs.index( chan.index ) ] += gfact * 4 * comp
            
        # get common factor 'beta' as a row vector:
        beta = numpy.pi / self.k(E)[:,0]**2
        elastic = beta * elasticSum
        absorbtion = beta * absorbtionSum
        total = elastic + absorbtion
        competitive = [ beta * comp for comp in competitiveSum ]
        fission = beta * fissionSum
        capture = absorbtion - fission - sum(competitive)
        
        retDict = {'total':total, 'elastic':elastic, 'capture':capture, 'fission':fission}
        competitiveNames = [ch.channel for ch in self.RR.openChannels if ch.tag=='competitive']
        for key,val in zip(competitiveNames,competitive):
            retDict[ key ] = val
        return retDict


##### unresolved resonance region. Only one formalism here: #####
class URRcrossSection(resonanceReconstructionBaseClass):
    def __init__(self, reactionSuite, verbose=True):
        super(URRcrossSection,self).__init__(reactionSuite)
        self.verbose = verbose
        
        # energy boundaries for this region:
        urr = reactionSuite.resonances.unresolved
        self.URR = urr.nativeData
        self.lowerBound = urr.lowerBound.getValueAs('eV')
        self.upperBound = urr.upperBound.getValueAs('eV')
        if self.verbose:
            print ("Unresolved from %f to %f eV" % (self.lowerBound,self.upperBound))
    
    def rho(self, E):
        # always calculate channel radius for unresolved:
        a = 0.123 * self.target.getMass('amu')**(1./3.) + 0.08
        return self.k(E) * a
    
    def getFluctuationIntegrals(self, widths, DOF):
        """ from subroutine GNRL3 in RECENT. If possible, this will be replaced
        with more basic approach (rather than using lookup table)... not finding 
        appropriate equations right now """
        XX = [[ 3.0013465e-03,1.3219203e-02,1.0004488e-03,1.3219203e-02,1.0e+0],
            [7.8592886e-02,7.2349624e-02,2.6197629e-02,7.2349624e-02,0.0e+0],
            [4.3282415e-01,1.9089473e-01,1.4427472e-01,1.9089473e-01,0.0e+0],
            [1.3345267e+00,3.9528842e-01,4.4484223e-01,3.9528842e-01,0.0e+0],
            [3.0481846e+00,7.4083443e-01,1.0160615e+00,7.4083443e-01,0.0e+0],
            [5.8263198e+00,1.3498293e+00,1.9421066e+00,1.3498293e+00,0.0e+0],
            [9.9452656e+00,2.5297983e+00,3.3150885e+00,2.5297983e+00,0.0e+0],
            [1.5782128e+01,5.2384894e+00,5.2607092e+00,5.2384894e+00,0.0e+0],
            [2.3996824e+01,1.3821772e+01,7.9989414e+00,1.3821772e+01,0.0e+0],
            [3.6216208e+01,7.5647525e+01,1.2072069e+01,7.5647525e+01,0.0e+0],]
        WW = [[1.1120413e-01,3.3773418e-02,3.3376214e-04,1.7623788e-03,1.0e+0],
            [2.3546798e-01,7.9932171e-02,1.8506108e-02,2.1517749e-02,0.0e+0],
            [2.8440987e-01,1.2835937e-01,1.2309946e-01,8.0979849e-02,0.0e+0],
            [2.2419127e-01,1.7652616e-01,2.9918923e-01,1.8797998e-01,0.0e+0],
            [1.0967668e-01,2.1347043e-01,3.3431475e-01,3.0156335e-01,0.0e+0],
            [3.0493789e-02,2.1154965e-01,1.7766657e-01,2.9616091e-01,0.0e+0],
            [4.2930874e-03,1.3365186e-01,4.2695894e-02,1.0775649e-01,0.0e+0],
            [2.5827047e-04,2.2630659e-02,4.0760575e-03,2.5171914e-03,0.0e+0],
            [4.9031965e-06,1.6313638e-05,1.1766115e-04,8.9630388e-10,0.0e+0],
            [1.4079206e-08,0.0000000e+00,5.0989546e-07,0.0000000e+00,0.0e+0],]
        
        if numpy.all(widths['captureWidth']<=0):
            return 0,0,0
        shape = widths['neutronWidth'].shape
        RN,RC,RF = numpy.zeros(shape), numpy.zeros(shape), numpy.zeros(shape)

        # We may store non-integer DOF, but we just end up converting to int:
        MUN = int( DOF['neutronDOF'] )
        if MUN<1 or MUN>4: MUN=5
        MUF = int( DOF['fissionDOF'] )
        if MUF<1 or MUF>5: MUF=5
        MUX=int( DOF['competitiveDOF'] )
        if MUX<1 or MUX>4: MUX=5
        for j in range(10):
            xj = XX[j][MUN-1]
            wj = WW[j][MUN-1]
            if numpy.any(widths['fissionWidthA']) and numpy.any(widths['competitiveWidth']):
                effj = widths['neutronWidth']*xj+widths['captureWidth']
                for k in range(10):
                    xk = XX[k][MUF-1]
                    wk = WW[k][MUF-1]
                    effjk = effj+widths['fissionWidthA']*xk
                    for i in range(10):
                        xi = XX[i][MUX-1]
                        wi = WW[i][MUX-1]
                        factor = wi*wk*wj*xj/(effjk+widths['competitiveWidth']*xi)
                        RN += xj*factor
                        RC += factor
                        RF += xk*factor
            elif numpy.any(widths['fissionWidthA']):
                effj = widths['neutronWidth']*xj+widths['captureWidth']
                for k in range(10):
                    xk = XX[k][MUF-1]
                    wk = WW[k][MUF-1]
                    factor = wk*wj*xj/(effj+widths['fissionWidthA']*xk)
                    RN += xj*factor
                    RC += factor
                    RF += xk*factor
            elif numpy.any(widths['competitiveWidth']):
                effj = widths['neutronWidth']*xj + widths['captureWidth']
                for k in range(10):
                    xk = XX[k][MUX-1]
                    wk = WW[k][MUX-1]
                    factor = wk*wj*xj/(effj+widths['competitiveWidth']*xk)
                    RN += xj*factor
                    RC += factor
            else:   # only neutron and capture widths:
                factor = wj*xj/(widths['neutronWidth']*xj+widths['captureWidth'])
                RN += xj*factor
                RC += factor
        for arr in (RN,RC,RF):
            arr[ widths['neutronWidth']<=0 ] = 0
        return RN,RC,RF

    def generateEnergyGrid(self):
        """Get the energy grid for reconstructing unresolved resonances. Usually this is just the energy grid
        chosen by the evaluator for storing energy-dependent widths."""
        egrid = set()
        for L in self.URR.L_values:
            for J in L.J_values:
                try: egrid.update( J.energyDependentWidths.getColumn('energy','eV') )
                except: egrid.update( [self.lowerBound, self.upperBound] )
        egrid = sorted( egrid )
        if egrid[0] > self.lowerBound or egrid[-1] < self.upperBound:
            raise Exception("Energy grid in unresolved region doesn't span stated energy range:" +
                    " (%e - %e) vs (%e - %e)" % (egrid[0], egrid[-1], self.lowerBound, self.upperBound) )
        return egrid

    def getCrossSection(self, interpolateWidths=False):
        """Get unresolved region cross section for each incident energy in E.
        By default, calculate on the energy grid supplied and then interpolate the resulting cross section.
        If interpolateWidths==True, interpolate the average widths during reconstruction instead."""
        E = numpy.array( self.generateEnergyGrid() )

        # below is a nasty hack since some ENDF files don't supply a dense enough energy grid in the URR.
        # These files should be fixed, but in the meantime add extra points and interpolate the widths:
        energyGridGaps = numpy.argwhere( E[1:] / E[:-1] > 3 ).flatten()
        if len( energyGridGaps ):
            import math
            if self.verbose:
                print "WARNING: filling large gaps in unresolved energy grid!"
            e_set = set( E )
            fixGrid = numpy.array( [1.0, 1.25, 1.5, 1.7, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.2, 8.5] )
            for idx in energyGridGaps:
                lowend, highend = E[idx], E[idx+1]
                lowPow, highPow = int( math.log10( lowend ) ), int( math.ceil( math.log10( highend ) ) )
                for i in xrange( lowPow, highPow + 1 ):
                    points_to_add = [point for point in 10**i * fixGrid if lowend < point < highend]
                    e_set.update( points_to_add )
            interpolateWidths=True
            E = numpy.array( sorted( e_set ) )
        # end of nasty hack
        
        captureSum = 0
        elasticSum = 0
        fissionSum = 0
        rho = self.rho(E)
        rhohat = self.URR.scatteringRadius.getValueAs('10*fm') * self.k(E)
        for L in self.URR.L_values:
            l = L.L
            phi = self.phi(l,rhohat)
            #S = 0.5*self.shiftFactor(l,rho)
            for J in L.J_values:
                j = J.J.value
                gfactor = (2.0*abs(j)+1)/(2*(2*self.spin+1))
                # dicts for widths and degrees of freedom.
                # widths may be constant or energy-dependent:
                # convert here from PhysicalQuantity to floats:
                widths, DOF = {}, {}
                for wid in ('neutronWidth','captureWidth','fissionWidthA',
                        'competitiveWidth','levelSpacing'):
                    widths[wid] = getattr(J.constantWidths, wid, None)
                    if type(widths[wid]) != type(None):
                        widths[wid] = widths[wid].getValueAs('eV')
                    elif (bool(J.energyDependentWidths)
                            and J.energyDependentWidths.getColumn( wid, units='eV' )):
                        if interpolateWidths:
                            from numericalFunctions import pointwiseXY_C
                            Nenergies = len(J.energyDependentWidths.getColumn('energy','eV'))
                            interp = pointwiseXY_C.pointwiseXY_C(
                                data = [J.energyDependentWidths.getColumn('energy','eV'), J.energyDependentWidths.getColumn(wid,'eV')],
                                dataForm = "XsAndYs", interpolation = self.URR.interpolation )
                            widthsNow = []
                            for en in E:
                                try:
                                    widthsNow.append( interp.getValue(en) )
                                except Exception as e:
                                    # grrr... data contains invalid log-log or lin-log interpolation
                                    print ("WARNING: unresolved resonance widths contain an invalid interpolation!")
                                    from fudge.core.math.xData.axes import linearToken
                                    interp2 = pointwiseXY_C.pointwiseXY_C( interp.copyDataToXYs(),
                                            interpolation = "%s,%s" % (linearToken,linearToken) )
                                    widthsNow.append( interp2.getValue(en) )
                            widths[wid] = numpy.array( widthsNow )
                        else:
                            widths[wid] = numpy.array( J.energyDependentWidths.getColumn(wid,'eV') )
                            if len(widths[wid]) != len(E):
                                raise Exception("Inconsistent energy arrays encountered in unresolved region!")
                for dof in ('neutronDOF','fissionDOF','competitiveDOF'):
                    DOF[dof] = getattr(J,dof,None)
                # penetration factor, different from resolved:
                VL = DOF['neutronDOF'] * self.penetrationFactor(l,rho) / rho
                widths['neutronWidth'] *= VL*numpy.sqrt(E)
                
                RN,RC,RF = self.getFluctuationIntegrals(widths,DOF)
                
                # common factor for all reactions:
                comfac = 2*numpy.pi*gfactor*widths['neutronWidth'] / widths['levelSpacing']
                
                captureSum += RC * widths['captureWidth'] * comfac
                fissionSum += RF * widths['fissionWidthA'] * comfac
                elasticSum += (RN * widths['neutronWidth']-2*numpy.sin(phi)**2) * comfac
            
            elasticSum += 4*(2*l+1)*numpy.sin(phi)**2
        
        # common factor 'beta' as row vector:
        beta = numpy.pi / self.k(E)**2
        capture = beta * captureSum
        elastic = beta * elasticSum
        fission = beta * fissionSum
        for reaction in (capture,elastic,fission):
            reaction[ reaction<=0 ] = 0
        total = elastic + capture + fission

        xscs = {'total':total, 'elastic':elastic, 'capture':capture, 'fission':fission}

        # convert to crossSection.pointwise instances, using interpolation specified in the evaluation:
        eInterp, xsInterp = self.URR.interpolation.split(',')
        for key in xscs:
            axes_ = crossSection.pointwise.defaultAxes( energyInterpolation=eInterp, crossSectionInterpolation=xsInterp )
            xscs[key] = crossSection.pointwise( axes_, zip(E,xscs[key]), 0.001 )
            xscs[key] = xscs[key].xSlice( self.lowerBound, self.upperBound )

        return xscs



