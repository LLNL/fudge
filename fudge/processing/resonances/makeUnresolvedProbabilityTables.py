#! /usr/bin/env python3
#encoding: utf-8

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import math
import numpy

from fudge.reactionData import crossSection as crossSectionModule
from fudge.resonances import common as commonResonancesModule, resolved as resolvedModule

from pqu import PQU as PQUModule

from xData import table as tableModule
from xData import axes as axesModule
from xData import XYs as XYsModule
from xData import multiD_XYs as multiD_XYsModule

from . import reconstructResonances

"""
Classes for generating probability tables P(crossSection | energy_in, T) from unresolved resonance parameters.
Methodology: create multiple random sets of unresolved resonance parameters, reconstruct cross sections,
optionally heat to one or more temperatures, and convert the results into probability tables.
The unresolved region is divided into sub-regions, with one pdf generated for each sub-region.

Usage example:
tableGenerator = ProbabilityTableGenerator( reactionSuite, verbose=True )
pdfs = tableGenerator.generatePDFs( resClass = reconstructResonances.SLBWcrossSection, nSamples = 10 )
"""

__metaclass__ = type

debug = False   # set True before debugging (disables multiprocessing)
if not debug:
    import matplotlib
    matplotlib.use('Agg')   # avoid trouble with matplotlib and multiprocessing

class ProbabilityTableGenerator:

    def __init__(self, reactionSuite, tolerance=None, verbose=False):

        self.reactionSuite = reactionSuite
        self.RRR = self.reactionSuite.resonances.resolved.evaluated.copy()
        self.URR = reconstructResonances.URRcrossSection( reactionSuite.resonances.unresolved.evaluated )

        if getattr(self.URR, 'averageWidths', None) is None:  # do some pre-sorting for convenience
            self.URR.getWidthsAndSpacings()

        self.tolerance = tolerance
        self.verbose = verbose

        # obtain background cross sections (will be added to each realization). These include both
        # the tabulated background and the contribution of resolved resonances in the URR
        self.backgrounds = {}
        if not self.URR.URR.useForSelfShieldingOnly:
            for reaction in self.URR.URR.resonanceReactions:
                xsc = reaction.reactionLink.link.crossSection.evaluated
                assert isinstance(xsc, crossSectionModule.resonancesWithBackground)
                key = reaction.label
                for simpleKey in ('elastic', 'capture', 'fission', 'total'):
                    if reaction.reactionLink.link is reactionSuite.getReaction(simpleKey):
                        key = simpleKey

                self.backgrounds[key] = xsc.background.unresolvedRegion.data

        # disabling for now: looks like URR parameters stand on their own, no need to add tails from RRR
        """
        if self.reactionSuite.resonances.resolved is not None:
            import copy
            self.RRR = copy.copy(self.reactionSuite.resonances.resolved.evaluated)     # FIXME should be a deep copy instead
            self.truncateResolvedRegion()   # drop resolved resonances in the unresolved region

            # get contribution of resolved resonances in resolved and unresolved regions:
            rrClass = reconstructResonances.getResonanceReconstructionClass(self.RRR)
            rrCrossSection = rrClass(self.RRR)
            egrid1 = rrCrossSection.generateEnergyGrid(highBound=self.lowerBound)   # extend up to bottom edge of URR

            egrid2 = list(self.levelDensities.values())[0].domainGrid       # extends past original URR range due to 'addPoints'
            egrid = numpy.array(sorted(set(egrid1 + egrid2)))
            xsecs = rrCrossSection.getCrossSection(egrid)
            egrid, xsecs, message = rrCrossSection.refineInterpolation(egrid, xsecs)
            for key in xsecs:
                xsNow = crossSectionModule.XYs1d( axes=crossSectionModule.defaultAxes('eV'), data=list(zip(egrid, xsecs[key])))
                if key in self.backgrounds:
                    self.backgrounds[key] += xsNow
                else:
                    self.backgrounds[key] = xsNow
        """

        # get mass ratio for heating cross sections:
        projectileMass = reactionSuite.PoPs[ reactionSuite.projectile ].getMass('amu')

        targetID = reactionSuite.target
        if targetID in reactionSuite.PoPs.aliases: targetID = reactionSuite.PoPs[targetID].pid
        target = reactionSuite.PoPs[targetID]

        self.massRatio = target.getMass('amu') / projectileMass

    def truncateResolvedRegion(self, domainMax=None):
        """
        Remove all resolved resonances with energies above domainMax.
        If domainMax is None, use self.lowerBound  (i.e. lower limit of the unresolved region).
        Modifies the local copy of the resolved parameters (self.RRR)
        """
        if domainMax is None: domainMax = self.lowerBound

        if isinstance(self.RRR, resolvedModule.BreitWigner):
            table = self.RRR.resonanceParameters.table
            energies = table.getColumn('energy')
            while energies and energies[-1] > domainMax:
                energies.pop()
                table.data.pop()
        elif isinstance(self.RRR, resolvedModule.RMatrix):
            for spinGroup in self.RRR.spinGroups:
                table = spinGroup.resonanceParameters.table
                energies = table.getColumn('energy')
                while energies and energies[-1] > domainMax:
                    energies.pop()
                    table.data.pop()

    def getLastResolvedResonanceRegion(self):
        """
        Get the highest energy resolved resonance region. Most evaluations only contain one resolved region,
        but a few (notably Pu239 in ENDF-VII.1) break the RRR into several regions.

        :return: The highest resolved resonance region, or None if the evaluation contains no resolved region
        """
        if self.RRR:
            # For multiple regions, we need to do each region separately, then add them to the unified xs table & egrid
            if isinstance(self.RRR, commonResonancesModule.energyIntervals):
                return self.RRR[-1]
            else:  # Single region, everything goes on unified grid
                return self.RRR
        else:
            return None

    def getLastResolvedResonanceEnergy(self, l=None, j=None):
        """
        Get the last resonance energy from the resolved region, that will start all the ladders in the URR

        :param l: orbital angular momentum of the resonance to find
        :param j: total angular momentum of the resonance to find
        :return: either the energy of the last resonance with requested (l,j) or None if it can't be found
        """
        resolved = self.getLastResolvedResonanceRegion()
        if isinstance(resolved, resolvedModule.BreitWigner):
            for x in reversed(self.getLastResolvedResonanceRegion().resonanceParameters.table):
                if x[0] > self.reactionSuite.resonances.resolved.domainMax:
                    # ignore resonances added above the end of resolved region (only included for their tails)
                    continue
                if x[1] == l and x[2] == j: return x[0]
        else:   # R-Matrix
            elastic, = [reac for reac in resolved.resonanceReactions if reac.reactionLink.link is
                        self.reactionSuite.getReaction('elastic')]
            elasticLabel = elastic.label

            maxRes = -1
            for spinGroup in resolved.spinGroups:
                if j is None or spinGroup.spin.value == j:
                    elasticChannels = [chan for chan in spinGroup.channels if chan.resonanceReaction == elasticLabel]
                    if l is None or any( [chan.L == l for chan in elasticChannels] ):
                        if len(spinGroup.resonanceParameters.table) > 0:
                            maxRes = max( maxRes, max( spinGroup.resonanceParameters.table.getColumn('energy') ) )
            if maxRes > 0: return maxRes

        return None

    def extrapolate_URR_parameters(self, newDomainMin, newDomainMax):
        """
        Extrapolate widths, level densities and possibly scattering radius beyond either end
        of the unresolved region, so we can draw realizations extending to higher/lower energy
        """
        def addPoints( function1d, newDomainMin, newDomainMax ):
            xys1d = function1d.toPointwise_withLinearXYs(lowerEps=1e-8)
            span = xys1d.domainMax - xys1d.domainMin

            x1,y1 = xys1d[0]
            x2,y2 = xys1d[1]
            newX = newDomainMin
            newY = y1 - (x1-newX) * (y2-y1)/(x2-x1)
            xys1d.setValue(newX, newY)

            x1,y1 = xys1d[-2]
            x2,y2 = xys1d[-1]
            newX = newDomainMax
            newY = y2 + (newX-x2) * (y2-y1)/(x2-x1)
            xys1d.setValue(newX, newY)
            return xys1d

        # make local copies of widths, densities and scattering radius rather than modifying the evaluation
        self.scatteringRadius = self.URR.URR.scatteringRadius.copy()
        if self.scatteringRadius.isEnergyDependent():
            self.scatteringRadius.form = addPoints( self.scatteringRadius.form, newDomainMin, newDomainMax )

        self.levelDensities = {}
        self.averageWidths = {}
        for lj in self.URR.levelDensities:

            self.levelDensities[lj] = addPoints( self.URR.levelDensities[lj].copy(), newDomainMin, newDomainMax )
            self.averageWidths[lj] = {}
            for key in self.URR.averageWidths[lj]:
                self.averageWidths[lj][key] = addPoints( self.URR.averageWidths[lj][key].copy(),
                        newDomainMin, newDomainMax )

                if key == 'elastic':
                    # ENDF neutron widths aren't actually neutron widths.
                    # Rather, they are what ENDF calls "reduced widths".
                    # We need to convert them back into widths
                    # FIXME: do this conversion at ENDF translation step instead?
                    for ip, p in enumerate(self.averageWidths[lj][key]):
                        rho = self.URR.rho(p[0])
                        f = math.sqrt(p[0]) * self.URR.penetrationFactor(lj[0], rho) / rho # See ENDF Eq. (D.95)
                        if lj[0] != 0:  # adjust widths since two channel spins are possible for this J,L
                            f *= 0.5
                        self.averageWidths[lj][key][ip] = [p[0], p[1]*f]

        self.lowerBound = newDomainMin
        self.upperBound = newDomainMax

    def sampleRR(self, lastResonanceEnergies, lowerBound=None, upperBound=None, style='goe', seed=None, verbose=True):
        """
        Generate a sample of a resolved resonance set using the average URR parameters

        :param lastResonanceEnergies: the energy of the last resolved resonance for each l/j
        :param lowerBound: optional. If supplied, resonances below lowerBound will be discarded
        :param upperBound: optional. If supplied, resonances above upperBound will be discarded
        :param style: method for generating resonance energies. Options include 'goe', 'wigner', 'picket fence', 'poisson', 'brody'
        :param seed: used to seed the random generator
        :param verbose: turn on verbose output
        :return: dictionary containing resonance parameters sorted by l/j
        """
        import brownies.BNL.restools.resonance_generator as rg

        fakeRRR={}

        # WARNING: be very careful that the columns defined in getFakeResonanceSet match
        #          in both name and order with those defined below
        for lj in self.levelDensities.keys():

            # First energy should be roughly D(E) above the last energy of the RRR
            # We add fake energies drawn from a Wigner distribution until the starting energy is above the lowerBound
            E0 = 0
            if style != 'goe':
                E0 = lastResonanceEnergies[lj]
                if E0 is None:      # resolved region doesn't have any resonances for this L/J combination
                    E0 = lowerBound - numpy.random.rayleigh(1.0/self.levelDensities[lj].evaluate(lowerBound))
                while E0 < lowerBound:
                    E0 += numpy.random.rayleigh(1.0/self.levelDensities[lj].evaluate(lowerBound))

            fakeRRR[lj]=rg.getFakeResonanceSet(
                E0=E0,
                style=style,
                L=lj[0], J=float(lj[1]),
                levelDensity=self.levelDensities[lj],
                aveWidthFuncs=self.averageWidths[lj],
                widthKeys = list( self.averageWidths[lj].keys() ), #('elastic','capture','fission'),
                DOFs=self.URR.DOFs[lj],
                domainMin=lowerBound,
                domainMax=upperBound,
                seed=seed,
                verbose=verbose)

            if style != 'goe' and fakeRRR[lj].data[0][0] != E0:
                # We hand in weird level densities that mess up the first energy point.
                # We'll fix that now
                if verbose: print("Adjusting %s,%s resonances by %f eV due to weird level densities" %
                                  (lj[0], lj[1], fakeRRR[lj].data[0][0] - E0))
                offset = fakeRRR[lj].data[0][0] - E0
                for ix, x in enumerate(fakeRRR[lj].data):
                    fakeRRR[lj].data[ix][0]-=offset

        return fakeRRR

    def generatePDFs(self, resClass, nSamples=10, temperatures=None, style='goe', interpolateWidths=True,
                     energyUnit='eV', temperatureUnit='K', verbose=False, plotSamples=0, plotDir='realizationPlots'):
        """
        Generate multiple resonance parameter realizations, reconstruct cross sections, optionally heat to one or
        more temperatures, and generate cross section PDFs at several incident energies.
        The PDFs from each sample are summed together, and only normalized after all samples are completed.

        :param resClass: The class instance to use for reconstructing realizations of the cross section set.
        :param nSamples: number of realizations generated and averaged to create the PDFs
        :param temperatures: list of temperatures for generating heated cross section PDFs.
        :param style: style of resonance realization generation (see sampleRR method for style options)
        :param interpolateWidths: interpolate the URR widths or not (just say "True")
        :param energyUnit: unit to use for energy grid in reconstructed cross sections
        :param temperatureUnit: unit for the list of temperatures. Default = 'K'
        :param verbose: enable verbose output
        :param plotSamples: number of cross section realizations to plot (helpful for debugging)
        :return: dictionary with probability tables (stored as XYs2d) for each temperature / reaction
        """

        if plotSamples:
            from matplotlib import pyplot as plt
            if not os.path.isdir(plotDir):
                os.mkdir(plotDir)

        if temperatures is None:
            temperatures = [0]
        temperaturesEnergy = sorted([PQUModule.PQU(T, temperatureUnit+'*k').getValueAs(energyUnit) for T in temperatures])

        # Initialize the widths, DOFs and level spacings.
        # Also setup the main egrid for cross section reconstruction.
        egrid, interpolateWidths = self.URR.generateEnergyGrid(interpolateWidths=interpolateWidths)
        LD=dict(self.URR.levelDensities)

        # Check that main grid spacing is compatible with parameter magnitude
        DeltaE=numpy.mean([egrid[i+1]-egrid[i] for i in range(len(egrid)-1)])
        for lj in LD:
            warnCount = 0
            for p in LD[lj]:
                if 10.0/p[1] > DeltaE:
                    warnCount += 1
                    if verbose:
                        print("WARNING: Energy grid too fine, DeltaE=%s < 10.0*D=%s at energy %s eV (L=%d, J=%g)"%(
                            str(DeltaE), str(10.0/p[1]), str(p[0]), lj[0], lj[1]))
            if warnCount:
                print("WARNING: Energy grid too fine for %d patches (L=%d, J=%g)" % (warnCount, lj[0], lj[1]))

            for rxn in self.URR.averageWidths[lj].keys():
                warnCount = 0
                widths = self.URR.averageWidths[lj][rxn]
                for incidentEnergy in widths.domainGrid:
                    widthAtEnergy = widths.evaluate(incidentEnergy)
                    if 10.0 * widthAtEnergy > DeltaE:
                        warnCount += 1
                        if verbose:
                            print("WARNING: Energy grid too fine, DeltaE=%s < 10.0*Gamma(%s)=%s at energy %s eV" % (
                                str(DeltaE), rxn, str(widthAtEnergy), str(incidentEnergy)))
                if warnCount:
                    print("WARNING: %s energy grid too fine for %d patches" % (rxn, warnCount))

        # Average number of resonances around a given energy in egrid determines the size of the
        # patch to use to enable "ergodicity trick" to accelerate convergence
        ljPatchE={}
        for lj in LD:
            ljPatchE[lj]=[]
            warnCount = 0
            for p in LD[lj]:
                # Figure for decent job at getting width distributions correct in a given sample,
                # need ~20 resonances in the patch
                patchSize=20.0/p[1]
                if patchSize > DeltaE:
                    warnCount += 1
                    if verbose:
                        print("WARNING: Patch size needed to average resonances around grid point %s for L=%i, J=%s > grid spacing"%(
                            str(p[0]), lj[0], str(lj[1])))
                    ljPatchE[lj].append(min(patchSize, DeltaE))
                else:
                    ljPatchE[lj].append(patchSize)
            if warnCount:
                print("WARNING: Patch size > grid spacing at %d energies for L=%i, J=%s" % (warnCount,lj[0],lj[1]))

        # Check all lj indexed patches have same number or energies in grid
        patchLengths=[len(ljPatchE[lj]) for lj in ljPatchE]
        if patchLengths.count(patchLengths[0])!=len(patchLengths):
            raise ValueError("One or more patch energy grids has the wrong number of values")

        # Take biggest patch at a given energy from all the lj's
        patchE=[0.0,] * patchLengths[0]
        for lj in ljPatchE:
            for i in range(patchLengths[0]):
                patchE[i]=max(patchE[i], ljPatchE[lj][i])

        patchLimits = []
        for iE, E in enumerate(egrid):
            # Make patch grid slightly larger than needed so that we can use the center region of width patchE
            # and the tails of the resonances outside the patch can do their job
            Elo=E - patchE[iE]/2.0
            Ehi=E + patchE[iE]/2.0
            patchLimits.append( (Elo, Ehi) )
        assert patchLimits[0][0] > 0

        self.extrapolate_URR_parameters( 0.95 * patchLimits[0][0], 1.05 * patchLimits[-1][0] )

        # Get highest-energy resolved resonance before start of realization for each L/J
        self.truncateResolvedRegion( self.lowerBound )
        lastResonances = {}
        for k in self.URR.levelSpacings:
            lastResonances[k] = self.getLastResolvedResonanceEnergy(l=k[0], j=k[1])

        # Loop through the samples to build the PDFs
        crossSectionSamples = []
        thePDFs={T:{} for T in temperatures}

        def to_XYs2ds(thePDFs, normalize=False):
            # helper method, used after all samples collected to pack final results.
            # May also be used to write intermediate results for debugging
            probabilityTableAxes = axesModule.axes(rank=3, labelsUnits={2:('energy_in',energyUnit),
                                                                        1:('crossSection','b'),
                                                                        0:('P(crossSection|energy_in)','')})
            newPDFs = {}
            for T in thePDFs:
                for rxn in thePDFs[T]:
                    probabilityTable = multiD_XYsModule.XYs2d(axes=probabilityTableAxes)
                    for Ein,PDF in zip(egrid,thePDFs[T][rxn]):
                        PDF = PDF.copy()
                        PDF.outerDomainValue = Ein
                        if normalize:
                            PDF.normalize(insitu=True)
                        PDF.axes = probabilityTableAxes
                        probabilityTable.append(PDF)
                    newPDFs.setdefault(T,{})[rxn] = probabilityTable

            return newPDFs


        for iSample in range(nSamples):
            print("Sample %d of %d" % (iSample + 1, nSamples))

            # Generate a realization of the resonances
            thisSample=self.sampleRR(lastResonances, style=style, lowerBound=self.lowerBound, verbose=verbose)

            # convert sample into GNDS resolved resonance:
            # FIXME currently only works for BreitWigner
            combinedData = []
            for key,table in thisSample.items():
                columnHeaders = table.columns
                combinedData.extend( table.data )

            # sort by energy (probably not necessary):
            combinedData.sort(key=lambda res: res[0])
            combinedTable = tableModule.table( columns=columnHeaders, data=combinedData )

            # FIXME following necessary since BreitWigner and URR use different naming conventions:
            renameColumns = {
                'elastic': 'neutronWidth',
                'capture': 'captureWidth',
                'fission': 'fissionWidth',
            }
            for column in combinedTable.columns:
                if column.name in renameColumns:
                    column.name = renameColumns[column.name]

            resolvedRealization = resolvedModule.BreitWigner(label="realization",
                approximation=resolvedModule.BreitWigner.singleLevel,
                resonanceParameters=commonResonancesModule.resonanceParameters( combinedTable ),
                scatteringRadius=self.scatteringRadius,
                PoPs=self.URR.URR.PoPs)

            # FIXME: next lines are required since the domain currently lives on <resolved> rather than the form
            resolvedContainer = resolvedModule.resolved(self.lowerBound, self.upperBound, self.URR.energyUnit)
            resolvedContainer.add(resolvedRealization)
            resolvedContainer.setAncestor(self.reactionSuite)

            # Setup the reconstruction class
            resonanceReconstructor = resClass(
                resolvedRealization,
                verbose=self.verbose,
                enableAngDists=False)

            # delay truncating negative cross sections until after adding the background contribution
            #resonanceReconstructor.removeNegatives = False

            egridNow = resonanceReconstructor.generateEnergyGrid(stride=5)
            zeroK_xscs = resonanceReconstructor.getCrossSection(egridNow)

            # convert to XYs1d
            zeroK_xscs = { key:crossSectionModule.XYs1d(
                axes=crossSectionModule.defaultAxes(energyUnit),
                data=(egridNow, zeroK_xscs[key]),
                dataForm="XsAndYs").dullEdges(lowerEps=1e-8)
                for key in zeroK_xscs }

            # add background cross section if applicable
            for key in zeroK_xscs:
                if key in self.backgrounds:
                    zeroK_xscs[key] += self.backgrounds[key]
                zeroK_xscs[key] = zeroK_xscs[key].clip(rangeMin=0)

            if iSample < plotSamples:
                crossSectionSamples.append(zeroK_xscs)

            if plotSamples and iSample == plotSamples-1:    # plot zero-K realizations
                figs = []
                for reaction in crossSectionSamples[0].keys():
                    if crossSectionSamples[0][reaction].rangeMax == 0:
                        continue
                    fig, ax = plt.subplots()
                    for ridx, sample in enumerate(crossSectionSamples):
                        ax.loglog(*sample[reaction].copyDataToXsAndYs(), label=str(ridx))
                    if self.reactionSuite.getReaction(reaction):
                        ax.loglog(*self.reactionSuite.getReaction(reaction).crossSection.toPointwise_withLinearXYs().copyDataToXsAndYs(),
                                label="evaluation")
                    ax.set_xlabel('Incident energy (%s)' % energyUnit)
                    ax.set_ylabel('Cross section (b)')
                    ax.set_xlim( 0.75 * self.lowerBound, 1.25 * self.upperBound )
                    ax.set_title("%s 0K cross section realizations" % reaction)

                    fig.savefig(os.path.join(plotDir, "%s.png" % reaction))
                    figs.append(fig)
                if debug:
                    plt.show()
                for fig in figs:
                    plt.close(fig)


            # next step: divide 0K cross sections into small regions, heat each region to desired temperature(s)
            # convert result to pdf(crossSection) and add to 'running total' pdf

            def heatAndMakePDFs(xs, ys, temperaturesEnergy, domainMin, domainMax, massRatio, **plotOptions):
                """
                Called for each cross section / incident energy range.
                Heats to all desired temperatures and computes pdf(crossSection) at each temperature.
                """
                initialXsc = crossSectionModule.XYs1d(
                        axes=crossSectionModule.defaultAxes(energyUnit),
                        data=(xs, ys), dataForm="XsAndYs")

                initialTemperature = 0
                results = []
                figs = []
                for nextTemp in temperaturesEnergy:
                    heatedXsc = initialXsc.heat(initialTemperature, nextTemp, massRatio, EMin=PQUModule.PQU(1e-5, 'eV'),
                            lowerlimit='constant', upperlimit=None, interpolationAccuracy=0.001, heatAllPoints=False, doNotThin=True,
                            heatBelowThreshold=True, heatAllEDomain=True, setThresholdToZero=False)

                    pdfNow = heatedXsc.pdfOfY(epsilon=1e-8, domainMin=domainMin, domainMax=domainMax)
                    results.append(pdfNow.copyDataToXsAndYs())

                    if plotOptions:
                        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(14,6))
                        ax[0].semilogy(*heatedXsc.copyDataToXsAndYs())
                        ax[0].set_xlabel("Incident energy (%s)" % plotOptions['energyUnit'])
                        ax[0].set_ylabel("Cross section (b)")
                        ax[0].axvline(domainMin, color='r', linewidth=1)
                        ax[0].axvline(domainMax, color='r', linewidth=1)

                        if 'elastic' in plotOptions['title'] or 'total' in plotOptions['title']:
                            c1 = ax[1].semilogy(*pdfNow.copyDataToXsAndYs(), label="pdf")
                        else:
                            c1 = ax[1].loglog(*pdfNow.copyDataToXsAndYs(), label="pdf")
                        ax[1].set_xlabel("Cross section (b)")
                        ax[1].set_ylabel("PDF(crossSection)")
                        ax2 = ax[1].twinx()
                        ax2.set_ylabel("CDF(crossSection)")
                        cdfNow = pdfNow.indefiniteIntegral()
                        if 'elastic' in plotOptions['title'] or 'total' in plotOptions['title']:
                            c2 = ax2.plot(*cdfNow.copyDataToXsAndYs(), color="g", label="cdf")
                        else:
                            c2 = ax2.semilogx(*cdfNow.copyDataToXsAndYs(), color="g", label="cdf")
                        curves = c1 + c2
                        ax[1].legend(curves, [c.get_label() for c in curves], loc=7)

                        fig.suptitle("%s - %g %s/k" % (plotOptions['title'], nextTemp, plotOptions['energyUnit']))
                        fig.tight_layout()
                        fig.subplots_adjust(top=0.93)
                        fig.savefig("%s_%g_%s_per_k.png" % (plotOptions['savePrefix'], nextTemp, plotOptions['energyUnit']))
                        figs.append(fig)

                    # heat incrementally rather than restarting at 0K each time:
                    initialTemperature = nextTemp
                    initialXsc = heatedXsc

                if debug:
                    plt.show()
                for fig in figs:
                    plt.close(fig)

                return results

            def unpackResults(index, label, results):
                """
                heatAndMakePDFs returns a list containing a pdfOfY at each temperature.
                Unpack those and add them to thePDFs
                """
                for temperature, result in zip(temperatures, results):
                    pdfNow = XYsModule.XYs1d(data=result, dataForm="XsAndYs")
                    patches = thePDFs[temperature].setdefault(label, [None] * len(patchLimits))
                    if patches[index] is None:      # 1st realization
                        patches[index] = pdfNow
                    else:                           # subsequent realizations
                        patches[index], current = patches[index].mutualify(1e-8, 1e-8, False,
                                pdfNow, 1e-8, 1e-8, False)
                        patches[index] += current
                        # thin to prevent average pdfs from consuming too much memory
                        patches[index] = patches[index].thin(1e-3)

            reactionsToHeat = set(zeroK_xscs.keys())
            for key in ('competitive', 'nonelastic'):
                reactionsToHeat.discard(key)
            if 'fission' in zeroK_xscs and zeroK_xscs['fission'].rangeMax == 0:
                reactionsToHeat.discard('fission')

            parameters = [(index, label, Elo, Ehi) for index, (Elo, Ehi) in enumerate(patchLimits) for label in reactionsToHeat]
            njobs = len(parameters)

            plotOptions = {}
            if iSample < plotSamples:   # plot realizations
                plotOptions['energyUnit'] = energyUnit

            # Add a buffer (factor of N * thermal velocity) before / after each cross section slice to account for heating edge effects.
            # The buffer is included when heating but excluded when computing pdfOfY.
            N = 10

            if debug:   # multiprocessing disabled
                for index, label, Elo, Ehi in parameters:
                    if plotOptions:
                        plotOptions['title'] = "%s patch %d" % (label, index)
                        plotOptions['savePrefix'] = os.path.join(plotDir, "r%d_patch%d_%s" % (iSample, index, label))

                    # temperaturesEnergy[-1] is the average thermal energy (in energyUnit) of highest temperature.
                    paddingLow = N * numpy.sqrt(2 * Elo * temperaturesEnergy[-1] / self.massRatio)
                    paddingHigh = N * numpy.sqrt(2 * Ehi * temperaturesEnergy[-1] / self.massRatio)

                    xs, ys = zeroK_xscs[label].domainSlice(Elo - paddingLow, Ehi + paddingHigh).copyDataToXsAndYs()
                    results = heatAndMakePDFs(xs, ys, temperaturesEnergy, Elo, Ehi, self.massRatio, **plotOptions)
                    unpackResults(index, label, results)

            else:
                from multiprocessing import Process, Queue
                import time
                from fudge.defaults import numTasks
                queue = Queue()

                # helper functions for running in parallel:
                def runAndEnqueue(queue, index, label, xs, ys, temperaturesEnergy, domainMin, domainMax, massRatio, plotOptions):
                    try:
                        results = heatAndMakePDFs(xs, ys, temperaturesEnergy, domainMin, domainMax, massRatio, **plotOptions)
                        queue.put((index, label, results))
                    except Exception as err:
                        queue.put((index, label, err))

                def startJob(job_index):
                    index, label, Elo, Ehi = parameters[job_index]
                    if plotOptions:
                        plotOptions['title'] = "%s patch %d" % (label, index)
                        plotOptions['savePrefix'] = os.path.join(plotDir, "r%d_patch%d_%s" % (iSample, index, label))

                    # temperaturesEnergy[-1] is the average thermal energy (in energyUnit) of highest temperature.
                    paddingLow = N * numpy.sqrt(2 * Elo * temperaturesEnergy[-1] / self.massRatio)
                    paddingHigh = N * numpy.sqrt(2 * Ehi * temperaturesEnergy[-1] / self.massRatio)

                    xs, ys = zeroK_xscs[label].domainSlice(Elo - paddingLow, Ehi + paddingHigh).copyDataToXsAndYs()
                    p = Process(target=runAndEnqueue, args=(queue, index, label, xs, ys, temperaturesEnergy, Elo, Ehi, self.massRatio, plotOptions))
                    p.start()
                    return p

                job_index = 0
                jobs = []
                for i1 in range(min(njobs, numTasks)):
                    jobs.append( startJob(job_index) )
                    job_index += 1

                done = 0
                while(done < njobs):
                    while queue.qsize():
                        index, label, results = queue.get()
                        if isinstance(results, Exception):
                            print("Error encountered in child process:")
                            raise results

                        unpackResults(index, label, results)
                    new_jobs = []
                    for process in jobs:
                        if process.is_alive():
                            new_jobs.append(process)
                        else:
                            done += 1
                            if process.exitcode != 0:
                                print("Job failed!")
                            if job_index < njobs:
                                new_jobs.append( startJob(job_index) )
                                job_index += 1
                    time.sleep(0.5)
                    jobs = new_jobs

            if False and iSample % 10 == 0:  # save intermediate results for debugging
                outDir = os.path.join(plotDir, "intermediateResults")
                if not os.path.exists(outDir):
                    os.makedirs(outDir)

                pdfsNow = to_XYs2ds(thePDFs, normalize=False)
                for temp in pdfsNow:
                    for reac in pdfsNow[temp]:
                        with open(os.path.join(outDir, "iteration%s_%s_%sK.xml" % (str(iSample).zfill(3), reac, temp)), "w") as fout:
                            fout.write('\n'.join(pdfsNow[temp][reac].toXMLList()))


        # all realizations complete. Convert results for each reaction / temperature into XYs2d
        return to_XYs2ds(thePDFs, normalize=True)

