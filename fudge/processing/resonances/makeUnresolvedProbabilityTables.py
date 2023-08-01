#! /usr/bin/env python3
# encoding: utf-8

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import math
import time
import bisect
import numpy

from fudge.reactionData import crossSection as crossSectionModule
from fudge.resonances import common as commonResonancesModule, resolved as resolvedModule

from pqu import PQU as PQUModule

from xData import table as tableModule
from xData import axes as axesModule
from xData import multiD_XYs as multiD_XYsModule
from xData import XYs1d as XYs1dModule

from . import reconstructResonances

"""
Classes for generating cross section probability densities P(crossSection | energy_in, T)
and probability tables (pdf for the total cross section along with conditional probabilities for each reaction)
from unresolved resonance parameters.
Methodology: create multiple random sets of unresolved resonance parameters, reconstruct cross sections,
optionally heat to one or more temperatures, and convert the results into pdfs and/or probability tables.
Results are given at various incident energies spanning the unresolved region.

Usage example:
tableGenerator = ProbabilityTableGenerator( reactionSuite, verbose=True )
pdfs = tableGenerator.generatePDFs( resClass = reconstructResonances.SLBWcrossSection, nSamples = 10 )
"""

debug = False   # set True before debugging (disables multiprocessing)
if not debug:
    import matplotlib

    matplotlib.use('Agg')  # avoid trouble with matplotlib and multiprocessing


class ProbabilityTableGenerator:

    def __init__(self, reactionSuite, tolerance=None, verbose=False):

        self.reactionSuite = reactionSuite
        self.RRR = None
        if self.reactionSuite.resonances.resolved:
            self.RRR = self.reactionSuite.resonances.resolved.evaluated.copy()
        self.URR = reconstructResonances.URRcrossSection(reactionSuite.resonances.unresolved.evaluated)

        # currently always use SLBW for URR
        self.reconstructionClass = reconstructResonances.SLBWcrossSection

        if getattr(self.URR, 'averageWidths', None) is None:  # do some pre-sorting for convenience
            self.URR.getWidthsAndSpacings()

        self.tolerance = tolerance
        self.verbose = verbose

        """ # disabled: adding background cross section to realizations biases samples
        # obtain background cross sections (will be added to each realization).
        self.backgrounds = {}
        if not self.URR.URR.useForSelfShieldingOnly:
            for reaction in self.URR.URR.resonanceReactions:
                xsc = reaction.link.link.crossSection.evaluated
                assert isinstance(xsc, crossSectionModule.ResonancesWithBackground)
                key = reaction.label
                for simpleKey in ('elastic', 'capture', 'fission', 'total'):
                    if reaction.link.link is reactionSuite.getReaction(simpleKey):
                        key = simpleKey

                self.backgrounds[key] = xsc.background.unresolvedRegion.data
        """

        # get mass ratio for heating cross sections:
        projectileMass = reactionSuite.PoPs[reactionSuite.projectile].getMass('amu')

        targetID = reactionSuite.target
        if targetID in reactionSuite.PoPs.aliases: targetID = reactionSuite.PoPs[targetID].pid
        target = reactionSuite.PoPs[targetID]

        self.massRatio = target.getMass('amu') / projectileMass

    def truncateResolvedRegion(self, domainMax=None):
        """
        Remove all resolved resonances with energies above domainMax.
        If domainMax is None, use self.lowerBound  (i.e. lower limit of the unresolved region).
        This method modifies the local copy of the resolved parameters (self.RRR).
        """
        if domainMax is None: domainMax = self.URR.lowerBound

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
        else:
            raise NotImplementedError("truncateResolvedRegion for resonances of type %s" % type(self.RRR))

    def getLastResolvedResonanceRegion(self):
        """
        Get the highest energy resolved resonance region. Most evaluations only contain one resolved region,
        but a few (notably Pu239 in ENDF-VII.1) break the RRR into several regions.

        :return: The highest resolved resonance region, or None if the evaluation contains no resolved region
        """
        if self.RRR:
            # For multiple regions, we need to do each region separately, then add them to the unified xs table & egrid
            if isinstance(self.RRR, commonResonancesModule.EnergyIntervals):
                return self.RRR[-1]
            else:  # Single region, everything goes on unified grid
                return self.RRR
        else:
            return None

    def getLastResolvedResonanceEnergy(self, l, j):
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
        elif resolved is not None:  # R-Matrix
            elastic, = [reac for reac in resolved.resonanceReactions if reac.link.link is
                        self.reactionSuite.getReaction('elastic')]
            elasticLabel = elastic.label

            maxRes = -1
            for spinGroup in resolved.spinGroups:
                if j is None or spinGroup.spin.value == j:
                    elasticChannels = [chan for chan in spinGroup.channels if chan.resonanceReaction == elasticLabel]
                    if l is None or any([chan.L == l for chan in elasticChannels]):
                        if len(spinGroup.resonanceParameters.table) > 0:
                            maxRes = max(maxRes, max(spinGroup.resonanceParameters.table.getColumn('energy')))
            if maxRes > 0: return maxRes

        return None

    def extrapolate_URR_parameters(self, newDomainMin, newDomainMax):
        """
        Extrapolate widths, level densities and possibly scattering radius beyond either end
        of the unresolved region, so we can draw realizations extending to higher/lower energy
        """

        def extrapolate(xnew, xy1, xy2, interpolation):
            """
            Use designated interpolation rule to extrapolate using 2 points xy1 and xy2 to get y-value at xnew.
            """
            #FIXME should be implemented in LUPY or numericalFunctions instead?
            import math
            from xData import enums
            if interpolation is enums.Interpolation.flat:
                if xnew > xy2[0]: return xy2[1]
                elif xnew < xy1[0]: return xy1[1]
                raise Exception("Oops!")

            xlog = interpolation in (enums.Interpolation.linlog, enums.Interpolation.loglog)
            ylog = interpolation in (enums.Interpolation.loglin, enums.Interpolation.loglog)

            x1,y1 = xy1
            x2,y2 = xy2
            if xlog:
                x1, x2, xnew = map(math.log, (x1, x2, xnew))
            if ylog:
                y1, y2 = map(math.log, (y1, y2))

            m = (y2-y1)/(x2-x1)
            ynew = y1 + m * (xnew-x1)
            if ylog:
                ynew = math.exp(ynew)
            return ynew

        def addPoints(function1d, newDomainMin, newDomainMax):
            from xData import constant, regions
            function1d = function1d.copy()  # don't modify original evaluation
            if isinstance(function1d, XYs1dModule.XYs1d):
                newY = extrapolate(newDomainMin, function1d[0], function1d[1], function1d.interpolation)
                function1d.setValue(newDomainMin, newY)

                newY = extrapolate(newDomainMax, function1d[-2], function1d[-1], function1d.interpolation)
                function1d.setValue(newDomainMax, newY)
            elif isinstance(function1d, regions.Regions1d):
                newY = extrapolate(newDomainMin, function1d[0][0], function1d[0][1], function1d[0].interpolation)
                function1d[0].setValue(newDomainMin, newY)

                newY = extrapolate(newDomainMax, function1d[-1][-2], function1d[-1][-1], function1d[-1].interpolation)
                function1d[-1].setValue(newDomainMax, newY)
            elif isinstance(function1d, constant.Constant1d):
                function1d.domainMin = newDomainMin
                function1d.domainMax = newDomainMax

            return function1d.toPointwise_withLinearXYs(lowerEps=1e-8)

        # extrapolate local copies of widths, densities and scattering radius, all as lin-lin functions
        self.scatteringRadius = self.URR.URR.getScatteringRadius().copy()
        self.scatteringRadius.form = addPoints(self.scatteringRadius.form, newDomainMin, newDomainMax)

        self.levelSpacings = {}
        self.levelDensities = {}
        self.averageWidths = {}
        for lj in self.URR.levelSpacings:

            self.levelSpacings[lj] = addPoints(self.URR.levelSpacings[lj], newDomainMin, newDomainMax)
            self.levelDensities[lj] = 1/self.levelSpacings[lj]  # FIXME doesn't add points like it should!
            self.averageWidths[lj] = {}
            for key in self.URR.averageWidths[lj]:
                self.averageWidths[lj][key] = addPoints(self.URR.averageWidths[lj][key],
                                                        newDomainMin, newDomainMax)

                if key == 'elastic':
                    # ENDF URR neutron widths aren't actually neutron widths.
                    # Rather, they are what ENDF calls "reduced widths".
                    # We need to convert them back into widths. Careful about units!
                    energyFactor = math.sqrt(PQUModule.PQU(1, self.averageWidths[lj][key].axes[-1].unit).getValueAs('eV'))
                    for ip, p in enumerate(self.averageWidths[lj][key]):
                        rho = self.URR.rho(p[0])
                        f = math.sqrt(p[0]) * energyFactor * self.URR.penetrationFactor(lj[0], rho) / rho  # See ENDF Eq. (D.95)
                        if lj[0] != 0:  # adjust widths since two channel spins are possible for this J,L
                            f *= 0.5
                        self.averageWidths[lj][key][ip] = [p[0], p[1] * f]

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

        fakeRRR = {}

        # WARNING: be very careful that the columns defined in getFakeResonanceSet match
        #          in both name and order with those defined below
        for lj in self.levelSpacings.keys():

            # First energy should be roughly D(E) above the last energy of the RRR
            # We add fake energies drawn from a Wigner distribution until the starting energy is above the lowerBound
            E0 = 0
            if style != 'goe':
                E0 = lastResonanceEnergies[lj]
                if E0 is None:  # resolved region doesn't have any resonances for this L/J combination
                    E0 = lowerBound - numpy.random.rayleigh(1.0 * self.levelSpacings[lj].evaluate(lowerBound))
                while E0 < lowerBound:
                    E0 += numpy.random.rayleigh(1.0 * self.levelSpacings[lj].evaluate(lowerBound))

            fakeResonances = rg.getFakeResonanceSet(
                E0=E0,
                style=style,
                L=lj[0], J=float(lj[1]),
                levelDensity=self.levelDensities[lj],
                aveWidthFuncs=self.averageWidths[lj],
                widthKeys=list(self.averageWidths[lj].keys()),  # ('elastic','capture','fission'),
                DOFs=self.URR.DOFs[lj],
                domainMin=lowerBound,
                domainMax=upperBound,
                seed=seed,
                verbose=verbose)

            if style != 'goe' and len(fakeResonances) > 0 and fakeResonances.data[0][0] != E0:
                # We hand in weird level densities that mess up the first energy point.
                # We'll fix that now
                if verbose: print("Adjusting %s,%s resonances by %f eV due to weird level densities" %
                                  (lj[0], lj[1], fakeResonances.data[0][0] - E0))
                offset = fakeResonances.data[0][0] - E0
                for ix, x in enumerate(fakeResonances.data):
                    fakeResonances.data[ix][0] -= offset

            fakeRRR[lj] = fakeResonances

        return fakeRRR

    def generatePDFs(self, nSamples=10, temperatures=None, style='goe', makePDFs=True, interpolateWidths=True,
                     temperatureUnit='K', verbose=False, plotSamples=0, plotDir='realizationPlots',
                     debugFile=None):
        """
        Generate multiple resonance parameter realizations, reconstruct cross sections, optionally heat to one or
        more temperatures, and generate cross section PDFs at several incident energies.
        The PDFs from each sample are summed together, and only normalized after all samples are completed.

        :param nSamples: number of realizations generated and averaged to create the PDFs
        :param temperatures: list of temperatures for generating heated cross section PDFs.
        :param style: style of resonance realization generation (see sampleRR method for style options)
        :param makePDFs: generate cross section PDFs in addition to probability tables
        :param interpolateWidths: interpolate the URR widths or not (just say "True")
        :param temperatureUnit: unit for the list of temperatures. Default = 'K'
        :param verbose: enable verbose output
        :param plotSamples: number of cross section realizations to plot (helpful for debugging)
        :param debugFile: optional file name to save realizations (will be saved as pickled dictionary)
        :return: dictionary with probability tables and optionally also pdfs
        """

        if plotSamples:
            from matplotlib import pyplot as plt
            if not os.path.isdir(plotDir):
                os.mkdir(plotDir)

        if temperatures is None:
            temperatures = [0]
        temperaturesEnergy = sorted(
            [PQUModule.PQU(T, temperatureUnit + '*k').getValueAs(self.URR.energyUnit) for T in temperatures])

        # Initialize the widths, DOFs and level spacings.
        # Also setup the main egrid for cross section reconstruction.
        egrid, interpolateWidths = self.URR.generateEnergyGrid(interpolateWidths=interpolateWidths)
        LS = dict(self.URR.levelSpacings)

        # Check that main grid spacing is compatible with parameter magnitude
        DeltaE = numpy.mean([egrid[i + 1] - egrid[i] for i in range(len(egrid) - 1)])
        for lj in LS:
            warnCount = 0
            for p in LS[lj].toPointwise_withLinearXYs():
                if 10.0 * p[1] > DeltaE:
                    warnCount += 1
                    if verbose:
                        print("WARNING: Energy grid too fine, DeltaE=%s < 10.0*D=%s at energy %s eV (L=%d, J=%g)" % (
                            str(DeltaE), str(10.0 * p[1]), str(p[0]), lj[0], lj[1]))
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
        ljPatchE = {}
        for lj in LS:
            ljPatchE[lj] = []
            warnCount = 0
            for incidentEnergy in egrid:
                levelSpacing = LS[lj].evaluate(incidentEnergy)
                # Figure for decent job at getting width distributions correct in a given sample,
                # need ~20 resonances in the patch
                patchSize = 40.0 * levelSpacing
                if patchSize > DeltaE:
                    warnCount += 1
                    if verbose:
                        print(
                            "WARNING: Patch size needed to average resonances around grid point %s for L=%i, J=%s > grid spacing" % (
                                str(incidentEnergy), lj[0], str(lj[1])))
                    ljPatchE[lj].append(min(patchSize, DeltaE))
                else:
                    ljPatchE[lj].append(patchSize)
            if warnCount:
                print("WARNING: Patch size > grid spacing at %d energies for L=%i, J=%s" % (warnCount, lj[0], lj[1]))

        # Check all lj indexed patches have same number of energies in grid
        patchLengths = [len(ljPatchE[lj]) for lj in ljPatchE]
        if patchLengths.count(patchLengths[0]) != len(patchLengths):
            raise ValueError("One or more patch energy grids has the wrong number of values")

        # Take biggest patch at a given energy from all the lj's
        patchE = [0.0] * patchLengths[0]
        for lj in ljPatchE:
            for i in range(patchLengths[0]):
                patchE[i] = max(patchE[i], ljPatchE[lj][i])

        patchLimits = []

        # Add a buffer (factor of N * thermal velocity) before / after each cross section 'patch' to account for heating edge effects.
        # The buffer is included when heating but excluded for computing probability tables / pdfs
        N = 10  # multiples to pad
        maxTemperatureEnergy = max(temperaturesEnergy)

        for iE, E in enumerate(egrid):
            # Make patch grid slightly larger than needed so that we can use the center region of width patchE
            # and the tails of the resonances outside the patch can do their job
            Elo = E - patchE[iE] / 2.0
            Ehi = E + patchE[iE] / 2.0

            if Elo <= 0:
                print("WARNING: reduced patch size for patch %d to avoid negative energies" % iE)
                Elo = 0.01 * E  # means patch is no longer symmetric around desired incident energy

            # additional padding to account for edge effects when Doppler broadening
            Elo_heating = Elo - N * numpy.sqrt(2 * Elo * maxTemperatureEnergy / self.massRatio)
            if Elo_heating <= 0:
                Elo_heating = 0.01 * Elo
            Ehi_heating = Ehi + N * numpy.sqrt(2 * Ehi * maxTemperatureEnergy / self.massRatio)
            patchLimits.append((Elo_heating, Elo, Ehi, Ehi_heating))
            assert all([l>0 for l in patchLimits[-1]])

        self.extrapolate_URR_parameters(min([p[0] for p in patchLimits]), max([p[-1] for p in patchLimits]))

        # Get highest-energy resolved resonance before start of realization for each L/J
        if self.RRR is not None:
            self.truncateResolvedRegion(self.lowerBound)
        lastResonances = {}
        for k in self.URR.levelSpacings:
            lastResonances[k] = self.getLastResolvedResonanceEnergy(l=k[0], j=k[1])

        self.patchLimits = patchLimits

        # Loop through the samples to build the PDFs
        parameterRealizations = []
        crossSectionSamples = []
        thePDFs = {T: {} for T in temperatures}

        theProbabilityTables = {
            'egrid': egrid,
            'samplePoints': [None] * len(patchLimits),
            'temperatures': {T: {} for T in temperatures}
        }

        def smoothPDF(pdf, maxPoints=200):
            """
            Ruthlessly smooth pdf to desired number of points by integrating over intervals.
            """
            numberOfPoints = len(pdf)
            pointsPerBin = numberOfPoints / maxPoints
            if numberOfPoints < maxPoints:
                return pdf

            data = []
            start = 0
            for index1 in range(maxPoints):
                end = min(int((index1 + 1) * pointsPerBin + 0.5), numberOfPoints)
                segment = pdf[start:end]
                width = segment.domainMax - segment.domainMin
                domainMid = 0.5 * (segment.domainMax + segment.domainMin)
                if index1 == 0:
                    data.append([segment.domainMin, segment.rangeMin])
                data.append([domainMid, segment.integrate() / width])
                start = end
            data.append([segment.domainMax, segment.rangeMin])
            return XYs1dModule.XYs1d(data=data)

        def to_XYs2ds(thePDFs, normalize=False, smooth=False):
            # helper method, used after all samples collected to pack final results.
            # May also be used to write intermediate results for debugging
            xsc_pdf_axes = axesModule.Axes(3, labelsUnits={2: ('energy_in', self.URR.energyUnit),
                                                                   1: ('crossSection', 'b'),
                                                                   0: ('P(crossSection|energy_in)', '')})
            newPDFs = {}
            for T in thePDFs:
                for rxn in thePDFs[T]:
                    probabilityTable = multiD_XYsModule.XYs2d(axes=xsc_pdf_axes)
                    for Ein, PDF in zip(egrid, thePDFs[T][rxn]):
                        PDF = PDF.copy()
                        if smooth:
                            maxPoints = min(200, len(PDF)//5)
                            PDF = smoothPDF(PDF, maxPoints=maxPoints).thin(0.001)    # FIXME hard-coded parameters
                        PDF.outerDomainValue = Ein
                        if normalize:
                            PDF.normalize(insitu=True)
                        PDF.axes = xsc_pdf_axes
                        probabilityTable.append(PDF)
                    newPDFs.setdefault(T, {})[rxn] = probabilityTable

            return newPDFs

        if nSamples == 0:
            return  # go through setup but don't draw any samples, useful for debugging

        for iSample in range(nSamples):
            print("Sample %d of %d" % (iSample + 1, nSamples))

            # random points where cross section will be sampled:
            samplePoints = numpy.sort(numpy.random.random(10000))

            # Generate a realization of the resonances
            start = time.time()
            thisSample = self.sampleRR(lastResonances, style=style, lowerBound=self.lowerBound, verbose=verbose)

            # convert sample into GNDS resolved resonance:
            # FIXME currently only works for BreitWigner
            nColumns = [table.nColumns for table in thisSample.values()]
            maxIndex = nColumns.index(max(nColumns))
            columnHeaders = list(thisSample.values())[maxIndex].columns

            combinedData = [[] for header in columnHeaders]
            for key, table in thisSample.items():
                for h_index, header in enumerate(columnHeaders):
                    vals = table.getColumn(header.name)
                    if vals is None:
                        vals = [0] * len(table)
                    combinedData[h_index].extend(vals)

            combinedData = numpy.array(combinedData).T.tolist()

            # sort by energy (probably not necessary):
            combinedData.sort(key=lambda res: res[0])
            combinedTable = tableModule.Table(columns=columnHeaders, data=combinedData)

            parameterRealizations.append(combinedData)

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
                approximation=resolvedModule.BreitWigner.Approximation.singleLevel,
                resonanceParameters=commonResonancesModule.ResonanceParameters( combinedTable ),
                scatteringRadius=self.scatteringRadius,
                PoPs=self.URR.URR.getLocalPoPs())
            if verbose:
                print("  drew %d resonances, elapsed time = %.3fs" % (len(combinedTable), time.time() - start))

            # FIXME: next lines are required since the domain currently lives on <resolved> rather than the form
            resolvedContainer = resolvedModule.Resolved(self.lowerBound, self.upperBound, self.URR.energyUnit)
            resolvedContainer.add(resolvedRealization)
            resolvedContainer.setAncestor(self.reactionSuite)

            # Setup the reconstruction class
            resonanceReconstructor = self.reconstructionClass(
                resolvedRealization,
                verbose=self.verbose,
                enableAngDists=False)

            # delay truncating negative cross sections to ensure total cross section is handled properly
            resonanceReconstructor.removeNegatives = False

            start = time.time()
            egridNow = resonanceReconstructor.generateEnergyGrid(stride=5)

            # down-select points to only compute inside each 'patch', including buffer for Doppler broadening:
            revisedGrid = set()
            for patchLow, tmp, tmp, patchHigh in patchLimits:
                i1 = bisect.bisect(egridNow, patchLow)
                i2 = bisect.bisect(egridNow, patchHigh)
                revisedGrid.update(egridNow[i1-1:i2+1])
                revisedGrid.update([patchLow, patchHigh])
            revisedGrid = sorted(revisedGrid)

            zeroK_xscs = resonanceReconstructor.getCrossSection(revisedGrid)
            if verbose:
                print("  reconstructed at 0 K, elapsed time = %.3fs" % (time.time() - start))

            # convert to XYs1d
            zeroK_xscs = {key: crossSectionModule.XYs1d(
                axes=crossSectionModule.defaultAxes(self.URR.energyUnit),
                data=(revisedGrid, zeroK_xscs[key]),
                dataForm="XsAndYs").dullEdges(lowerEps=1e-8)
                          for key in zeroK_xscs}

            for key in zeroK_xscs:
                # Following lines are disabled for now: adding background cross section biases prob. distribution.
                # Probably needs more discussion with other processing codes, though.
                """
                # add background cross section if applicable
                if key in self.backgrounds:
                    recon, background = zeroK_xscs[key].mutualify(1e-8, 1e-8, False,
                            self.backgrounds[key], 1e-8, 1e-8, False)
                    zeroK_xscs[key] = recon + background
                """
                zeroK_xscs[key] = zeroK_xscs[key].clip(rangeMin=0)

            # may need to recompute total after trimming zeros:
            if zeroK_xscs['total'].rangeMin == 0:
                resummed = crossSectionModule.XYs1d(axes=crossSectionModule.defaultAxes(self.URR.energyUnit))
                for key in zeroK_xscs:
                    if key not in ('total', 'nonelastic'):
                        resummed += zeroK_xscs[key]
                zeroK_xscs['total'] = resummed

            if iSample < plotSamples:
                crossSectionSamples.append(zeroK_xscs)

            if plotSamples and iSample == plotSamples - 1:  # plot zero-K realizations
                figs = []
                for reaction in crossSectionSamples[0].keys():
                    if crossSectionSamples[0][reaction].rangeMax == 0:
                        continue
                    fig, ax = plt.subplots(figsize=(10,6))
                    for ridx, sample in enumerate(crossSectionSamples):
                        ax.loglog(*sample[reaction].copyDataToXsAndYs(), label=str(ridx), linewidth=0.5)
                    if self.reactionSuite.getReaction(reaction):
                        ax.loglog(*self.reactionSuite.getReaction(
                            reaction).crossSection.toPointwise_withLinearXYs().copyDataToXsAndYs(),
                                  label="evaluation")
                    ax.set_xlabel('Incident energy (%s)' % self.URR.energyUnit)
                    ax.set_ylabel('Cross section (b)')
                    ax.set_xlim(0.75 * self.lowerBound, 1.25 * self.upperBound)
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
                    axes=crossSectionModule.defaultAxes(self.URR.energyUnit),
                    data=(xs, ys), dataForm="XsAndYs")

                initialTemperature = 0
                results = []
                figs = []
                for nextTemp in temperaturesEnergy:
                    heatedXsc = initialXsc.heat(initialTemperature, nextTemp, massRatio, EMin=PQUModule.PQU(1e-5, 'eV'),
                                                lowerlimit='constant', upperlimit=None, interpolationAccuracy=0.001,
                                                heatAllPoints=False, doNotThin=True,
                                                heatBelowThreshold=True, heatAllEDomain=True, setThresholdToZero=False)

                    xscSamplePoints = domainMin + samplePoints * (domainMax - domainMin)
                    xscVals = [heatedXsc.evaluate(x) for x in xscSamplePoints]
                    pdfNow = None
                    if makePDFs:
                        pdfNow = heatedXsc.pdfOfY(epsilon=1e-8, domainMin=domainMin, domainMax=domainMax)
                    results.append((xscVals, pdfNow))

                    if makePDFs and plotOptions:
                        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(14, 6))
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
                        fig.savefig(
                            "%s_%g_%s_per_k.png" % (plotOptions['savePrefix'], nextTemp, plotOptions['energyUnit']))
                        figs.append(fig)

                    # heat incrementally rather than restarting at 0K each time:
                    initialTemperature = nextTemp
                    initialXsc = heatedXsc

                if plotOptions and debug:
                    plt.show()
                for fig in figs:
                    plt.close(fig)

                return results

            def unpackResults(index, label, results):
                """
                heatAndMakePDFs returns a list containing a pdfOfY and probTable at each temperature.
                Unpack those and add them to thePDFs and theProbabilityTables
                """
                for temperature, (resultPT, resultPDF) in zip(temperatures, results):
                    patches = theProbabilityTables['temperatures'][temperature].setdefault(label, [None] * len(patchLimits))
                    if patches[index] is None:  # 1st realization
                        patches[index] = [resultPT]
                    else:
                        patches[index].append(resultPT)

                    if makePDFs:
                        patches = thePDFs[temperature].setdefault(label, [None] * len(patchLimits))
                        if patches[index] is None:  # 1st realization
                            patches[index] = resultPDF
                        else:  # subsequent realizations
                            patches[index], current = patches[index].mutualify(1e-8, 1e-8, False,
                                    resultPDF, 1e-8, 1e-8, False)
                            patches[index] += current
                            # thin to prevent average pdfs from consuming too much memory
                            patches[index] = patches[index].thin(1e-3)


            reactionsToHeat = set(zeroK_xscs.keys())
            #reactionsToHeat.discard('total')  # recompute total after heating other reactions
            for key in ('competitive', 'nonelastic'):   # unused
                reactionsToHeat.discard(key)
            if 'fission' in zeroK_xscs and zeroK_xscs['fission'].rangeMax == 0:
                reactionsToHeat.discard('fission')

            # save incident energy grids for each realization (for debugging only)
            for index, (Elo_heating, Elo, Ehi, Ehi_heating) in enumerate(patchLimits):
                if theProbabilityTables['samplePoints'][index] is None:
                    theProbabilityTables['samplePoints'][index] = []
                theProbabilityTables['samplePoints'][index].append(
                    Elo + samplePoints * (Ehi - Elo) 
                )

            parameters = [(index, label, Elo_heating, Elo, Ehi, Ehi_heating) for index, (Elo_heating, Elo, Ehi, Ehi_heating) in enumerate(patchLimits)
                          for label in reactionsToHeat]
            njobs = len(parameters)

            plotOptions = {}
            if iSample < plotSamples:  # plot realizations
                plotOptions['energyUnit'] = self.URR.energyUnit

            start = time.time()
            if debug:  # multiprocessing disabled
                for index, label, Elo_heating, Elo, Ehi, Ehi_heating in parameters:
                    if plotOptions:
                        plotOptions['title'] = "%s patch %d" % (label, index)
                        plotOptions['savePrefix'] = os.path.join(plotDir, "r%d_patch%d_%s" % (iSample, index, label))

                    xs, ys = zeroK_xscs[label].domainSlice(Elo_heating, Ehi_heating).copyDataToXsAndYs()
                    results = heatAndMakePDFs(xs, ys, temperaturesEnergy, Elo, Ehi, self.massRatio, **plotOptions)
                    unpackResults(index, label, results)

            else:
                from multiprocessing import Process, Queue
                from fudge.defaults import numTasks
                queue = Queue()

                # helper functions for running in parallel:
                def runAndEnqueue(queue, index, label, xs, ys, temperaturesEnergy, domainMin, domainMax, massRatio,
                                  plotOptions):
                    try:
                        results = heatAndMakePDFs(xs, ys, temperaturesEnergy, domainMin, domainMax, massRatio,
                                                  **plotOptions)
                        queue.put((index, label, results))
                    except Exception as err:
                        print(f"Error encountered for heatAndMakePDFs for {index}, {label}, {domainMin}, {domainMax}, {xs[:4]}, {xs[-4:]}")
                        queue.put((index, label, err))

                def startJob(job_index):
                    index, label, Elo_heating, Elo, Ehi, Ehi_heating = parameters[job_index]
                    if plotOptions:
                        plotOptions['title'] = "%s patch %d" % (label, index)
                        plotOptions['savePrefix'] = os.path.join(plotDir, "r%d_patch%d_%s" % (iSample, index, label))

                    xs, ys = zeroK_xscs[label].domainSlice(Elo_heating, Ehi_heating).copyDataToXsAndYs()
                    p = Process(target=runAndEnqueue, args=(queue, index, label, xs, ys, temperaturesEnergy, Elo, Ehi, self.massRatio, plotOptions))
                    p.start()
                    return p

                job_index = 0
                jobs = []
                for i1 in range(min(njobs, numTasks)):
                    jobs.append(startJob(job_index))
                    job_index += 1

                done = 0
                errors = []
                while (done < njobs):
                    while not queue.empty():
                        index, label, results = queue.get()
                        if isinstance(results, Exception):
                            print("Error encountered in child process!")
                            errors.append(results)
                            continue

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
                                new_jobs.append(startJob(job_index))
                                job_index += 1
                    time.sleep(0.5)
                    jobs = new_jobs

                if len(errors) > 0:
                    print("%d child processes failed:" % len(errors))
                    raise errors[0]

            if False and iSample % 10 == 0:  # save intermediate results for debugging
                outDir = os.path.join(plotDir, "intermediateResults")
                if not os.path.exists(outDir):
                    os.makedirs(outDir)

                pdfsNow = to_XYs2ds(thePDFs, normalize=False, smooth=False)
                for temp in pdfsNow:
                    for reac in pdfsNow[temp]:
                        with open(os.path.join(outDir, "iteration%s_%s_%sK.xml" % (str(iSample).zfill(3), reac, temp)), "w") as fout:
                            fout.write('\n'.join(pdfsNow[temp][reac].toXML_strList()))

            if verbose:
                print("  heat and accumulate probabilities, elapsed time = %.3fs\n" % (time.time() - start))

        # convert probability table samples to numpy arrays for convenience
        for key in ('egrid', 'samplePoints'):
            theProbabilityTables[key] = numpy.array(theProbabilityTables[key])
        for temperature in theProbabilityTables['temperatures'].values():
            for label in temperature:
                temperature[label] = numpy.array(temperature[label])

        if debugFile is not None:
            theProbabilityTables['parameterRealizations'] = parameterRealizations
            import pickle
            with open(debugFile, "wb") as fout:
                pickle.dump(theProbabilityTables, fout)

        def genPTs(egrid, theProbabilityTables, normalize=True):
            resultPTs = {}
            nbins = 50
            probs = [1 / nbins] * nbins
            for temperature, probTable in theProbabilityTables['temperatures'].items():
                resultPTs[temperature] = []
                for eidx, incidentEnergy in enumerate(egrid):
                    total = probTable['total'][eidx].flatten()
                    order = numpy.argsort(total)
                    # sort all reactions in ascending order by total cross section:
                    total = total[order]
                    conditionals = {}
                    for key in probTable:
                        if key == 'total': continue
                        conditionals[key] = probTable[key][eidx].flatten()[order]

                    step = len(total) // nbins
                    PTs = {'total': []}
                    for idx in range(nbins):
                        PTs.setdefault('total', []).append(total[step * idx:step * (idx + 1)].mean())
                        for key in conditionals:
                            PTs.setdefault(key, []).append(conditionals[key][step * idx:step * (idx + 1)].mean())

                    if normalize:
                        for key in PTs:
                            mean = numpy.mean(PTs[key])
                            PTs[key] /= mean
                
                    resultPTs[temperature].append([incidentEnergy, probs, PTs])
            
            return resultPTs

        resultPTs = genPTs(egrid, theProbabilityTables, normalize=True)
        if makePDFs:
            thePDFs = to_XYs2ds(thePDFs, normalize=True, smooth=True)

        # all realizations complete, return probability tables and optionally cross section PDFs
        result = {
            'probabilityTables': resultPTs,
            'pdfs': thePDFs
        }
        return result
