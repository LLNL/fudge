#!/usr/bin/env python
#encoding: utf-8

# <<BEGIN-copyright>>
# <<END-copyright>>

"""
reconstruct resonance region
cmattoon, 12/1/2010

See ENDF-102 (endf documentation) appendix D for equations

Basic usage:

    if x is a reactionSuite instance,

        >>xsecs = fudgeReconstructResonances.reconstructResonances(x, tolerance=0.01)

    This reconstructs all resolved/unresolved resonance sections. In each section,
    the results are accurate under linear interpolation to tolerance of 1% or better.
    All sections are summed together.

    * xsecs = dictionary containing cross sections: {'total':XYs, 'elastic':XYs, ...}
    each cross section is an XYs class instance, containing data and also axes with units


Alternate use:

    If desired, reconstruct a single section:

    >> resCls = fudgeReconstructResonances.RMcrossSection(x)   # for Reich_Moore
    >> energy_grid = s.generateEnergyGrid()
    >> crossSections = s.getCrossSection( energy_grid )
    # the input to getCrossSection is the energy (or list of energies) in eV
    # crossSections are returned as a dictionary {'total':,'elastic':,'capture':,'fission':,}

    # improve grid to desired tolerance for linear interpolation:
    >> new_energy_grid, new_crossSections, messages = s.refineInterpolation(energy_grid, crossSections, tolerance=0.01)
"""
import numpy
from fudge.gnd.reactionData import crossSection

__metaclass__ = type


def reconstructResonances(reactionSuite, tolerance=None, verbose=True):
    """ reconstruct all resonance sections (resolved/unresolved) in reactionSuite, 
    add results together for full (resonance region) pointwise cross section.
    If tolerance is specified, refine grid to required tolerance for (lin-lin) interpolation """
    
    if not reactionSuite.resonances.reconstructCrossSection:
        print ("      Resonance reconstruction was already done. Reconstruct again? y/N")
        if not raw_input().lower().startswith('y'):
            return xsecs
    
    egrids, xsecs = [], []
    if reactionSuite.resonances.resolved:
        def resolvedReconstruct( nativeData, sectionIndex = None ):
            if nativeData.moniker == 'SingleLevel_BreitWigner': resCls = SLBWcrossSection
            elif nativeData.moniker == 'MultiLevel_BreitWigner': resCls = MLBWcrossSection
            elif nativeData.moniker == 'Reich_Moore': resCls = RMcrossSection
            elif nativeData.moniker == 'R_Matrix_Limited': resCls = RMatrixLimitedcrossSection
            else:
                raise TypeError, "Don't recognize resonance type %s" % nativeData.moniker
            reconstructClass = resCls( reactionSuite, sectionIndex, verbose )
            egrid = reconstructClass.generateEnergyGrid()
            xsecs_now = reconstructClass.getCrossSection( egrid )
            if tolerance:
                egrid, xsecs_now, messages = reconstructClass.refineInterpolation(numpy.array(egrid), xsecs_now,
                        tolerance)
                if verbose:
                    for message in messages: print (message)
            return egrid, xsecs_now
        
        if reactionSuite.resonances.resolved.multipleRegions:
            if verbose:
                print( "      WARNING! Multiple resolved/unresolved energy regions are deprecated\n"
                        +"        and should be consolidated into single section")
            for RRidx in range(len(reactionSuite.resonances.resolved.regions)):
                nativeData = reactionSuite.resonances.resolved.regions[RRidx].nativeData
                egrid_now, xsecs_now = resolvedReconstruct( nativeData, RRidx )
                # merge with 'full' energy grid:
                egrids.append( egrid_now )
                xsecs.append( xsecs_now )
        else:
            nativeData = reactionSuite.resonances.resolved.nativeData
            egrid_now, xsecs_now = resolvedReconstruct( nativeData )
            # merge with 'full' energy grid:
            egrids.append( egrid_now )
            xsecs.append( xsecs_now )

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
        # sqrt(2*neutronMass)/hbar == 2.196807 (eV*barn)**-1/2. Thus for energy in eV, k is in b**-1/2
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
        # get the channel radius, rho. If L is specified try to get L-dependent value
        if self.RR.calculateChannelRadius:
            a = 0.123 * self.target.getMass('amu')**(1./3.) + 0.08
        else:
            if self.RR.scatteringRadius.isEnergyDependent():
                a = self.RR.scatteringRadius.getValueAs('10**-12*cm', E[:,0])
            else: a = self.RR.getScatteringRadius(L).getValueAs('10**-12*cm')
        return self.k(E) * a
    
    def generateEnergyGrid(self):
        """ create the energy grid by merging central resonance energies, fine mesh between resonances,
        and a rough mesh covering the entire region. use the total widths """
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
    """
    def __init__(self, reactionSuite, sectionIndex=None, verbose=True):
        super(SLBWcrossSection, self).__init__(reactionSuite, sectionIndex)
        self.verbose = verbose
        self.sortLandJ()
        if self.verbose:
            print ("From %f to %f eV, reconstructing using Single-level Breit-Wigner" % 
                    (self.lowerBound,self.upperBound))
    
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
    given a resonance region in MLBW format,
    create a class with all data required to reconstruct
    
    only the elastic channel differs from SLBW
    """
    def __init__(self, reactionSuite, sectionIndex=None, verbose=True):
        super(MLBWcrossSection, self).__init__(reactionSuite, sectionIndex)
        self.verbose = verbose
        self.sortLandJ()
        if self.verbose:
            print ("From %f to %f eV, reconstructing using Multi-level Breit-Wigner" % 
                    (self.lowerBound,self.upperBound))
    
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
    
    incident energy dependence appears in both in E and the penetrabilities
    
    additional documentation from RECENT::
    
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
    def __init__(self, reactionSuite, sectionIndex=None, verbose=True):
        super(RMcrossSection, self).__init__(reactionSuite, sectionIndex)
        self.verbose = verbose
        self.sortLandJ()
        if self.verbose:
            print ("From %f to %f eV, reconstructing using Reich_Moore (LRF=3)" % 
                    (self.lowerBound,self.upperBound))

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
    def __init__(self, reactionSuite, sectionIndex=None, verbose=True):
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
            a = [0.123 * self.target.getMass('amu')**(1./3.) + 0.08 
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
                                dataForm = "XsAndYs", interpolation = J.interpolation )
                            widthsNow = []
                            for en in E:
                                try:
                                    widthsNow.append( interp.getValue(en) )
                                except Exception as e:
                                    # grrr... data contains invalid log-log or lin-log interpolation
                                    print ("WARNING: unresolved resonance widths contain an invalid interpolation!")
                                    from fudge.core.math.xData.axes import linearToken
                                    interp2 = pointwiseXY_C.pointwiseXY_C( interp, interpolation = "%s,%s" % (linearToken,linearToken) )
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



