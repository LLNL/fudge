from fudge.core.utilities.brb import banner, winged_banner
from fudge import physicalQuantity as physicalQuantityModule
from fudge.vis.matplotlib import plot2d
from xData import table as tableModule
import xData.enums as enumsModule
from brownies.BNL.inter.datatables import *
from brownies.BNL.inter.metrics import *
from brownies.BNL.utilities.html import *
from brownies.BNL.utilities.XYs import *

"""
===========================
Everything to be catalogued
===========================

    Evaluation (dict)
    =================

        - Author (str)
        - Lab (str)
        - Date (date)
        - MAT (int)
        - Target (str)
        - Projectile (str)
        - Compound nucleus formed (str)
        - Projectile frame (str: "lab" or "cm")
        - Temperature (PQU)

    Fissile metrics (dict)
    ======================
        - ETA (PQU)
        - ALF (PQU)

    Resonances (dict)
    =================

        - Particle pair list (list)
        - R' (PQU)

        RRR (dict)
        ----------

            - Format (str)
                -- approximation (str)
            - Lower bound (PQU)
            - Upper bound (PQU)
            - Lmax (int)
            - No. resonances (int)
            - No. channels (int)

        URR (dict)
        ----------

            - Format (str)
                -- approximation (str)
            - Lower bound (PQU)
            - Upper bound (PQU)
            - No. channels (int)

        Per Particle Pair (list)
        ------------------------

           - index (int)
           - name (str)
           - Particle1 (dict)
                - name (str)
                - mass (PQU)
                - spin (str)
                - parity (str: '+' or '-')
                - charge (int)
            - Particle2 (dict)
                - name (str)
                - mass (PQU)
                - spin (str)
                - parity (str: '+' or '-')
                - charge (int)

        Per Channel (datatable)
        ------------------

            - c (J,L,S,Pi,PP) (struct)
            - MT (int)
            - Reaction string (str)
            - Q (PQU)
            - Eth = Xi (PQU)
            - nu (float)
            - Eliminated? (bool)
            - Competative? (bool)
            - Relativisitic kinematics? (bool)
            - For potential scattering only? (bool)
            - No. resonances (dict)
                -- RRR (int)
                -- estimated URR (float)
            - List of resonances with ER < 0 (list)
            - Excitation functions
                -- Staircase: Cum. no. resonances vs. E (plot)
                -- <D> vs. E (plot)
                -- <Gamma> vs. E (plot)
                -- Transmission coefficients: Tc vs. E (plot)
                -- Scattering radii: Rc vs. E (plot)
                -- Strength functions: Sc vs. E (plot)
            - Overridden defaults
                -- Penetribility: Pc vs. E (plot)
                -- Shift: Sc vs. E (plot)
                -- Phase: Phic vs. E (plot)
                -- Background cross section present? (bool)
            - GOE plots
                -- Delta3 (plot)
                -- Porter-Thomas, nu fixed, fit (Gamma/<Gamma>)_min (plot)
                -- Porter-Thomas, (Gamma/<Gamma>)_min fixed, fit nu (plot)

    Per Reaction (datatable)
    ========================

        - MT (int)
        - Reaction string (str)
        - No. products (int)
        - Residual (str)
        - Q (PQU)
        - Eth == Xi (PQU)
        - Reaction rates
            -- Godiva (PQU)
            -- Flattop (PQU)
            -- Jezebel (PQU)
            -- BigTen (PQU)
            -- FUND-IPPE (PQU)
            -- Cf252 (analytic) (PQU)
            -- Cf252 (PQU)
        - Engineering metrics
            -- RI (PQU)
            -- Westcott factor (PQU)
            -- Thermal xs (PQU)
            -- E14 (PQU)
        - Astrophysics metrics
            -- MACS(30 keV) (PQU)
            -- ARR (PQU)
        - INTER metrics
            -- RI (PQU)
            -- Westcott factor (PQU)
            -- Thermal xs (PQU)
            -- Cf252 (analytic) (PQU)
            -- Cf252 (PQU)
            -- MACS(30 keV) (PQU)
            -- E14 (PQU)

"""


# -------------------------------------------------------------------------------
# Report formatting classes
# -------------------------------------------------------------------------------

class ComplexEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, complex):
            return [obj.real, obj.imag]
        if isinstance(obj, set):
            return tuple(obj)
        if isinstance(obj, PQU.PQU):
            return str(obj)
        if isinstance(obj, physicalQuantityModule.Temperature):
            return str(obj)
        if isinstance(obj, DataTable):
            return obj.as_dict()
        if isinstance(obj, datatables.DataTable):
            return obj.as_dict()
        if isinstance(obj, enumsModule.Frame):
            return str(obj)
        # Let the base class default method raise the TypeError
        try:
            return json.JSONEncoder.default(self, obj)
        except TypeError as err:
            print("Encoding TypeError: %s" % err)
            return str(obj)

class Report(collections.OrderedDict):

    def __init__(self, title=None, indent=4, *args, **kwds):
        self.title = title
        self.indent = indent
        collections.OrderedDict.__init__(self, *args, **kwds)

    @property
    def id(self):
        return id(self)

    @property
    def cls(self):
        return 'report %s' % ''.join([x for x in self.title if x.isalnum()])

    def text_title(self):
        return banner(self.title)

    def text_report(self, *args, **kwds):
        return self.json_report(*args, indent=self.indent, separators=(',', ': '), **kwds)

    def json_report(self, *args, **kwds):
        if self.title is None:
            return json.dumps(self, *args, cls=ComplexEncoder, **kwds)
        return json.dumps({self.title: self}, *args, cls=ComplexEncoder, **kwds)

    def html_report(self, *args, **kwds):
        # FIXME: don't work yet
        report = [H1(self.title), '\n'.join(OBJECT2HTML(self))]
        if True:
            report = [line.replace('Gamma', SYMBOL('Gamma')) for line in report]
        return report


class Plot:

    def __init__(self, title=None, datasets=None, xaxis=None, yaxis=None):
        self.title = title
        if datasets is None:
            self.datasets = []
        else:
            self.datasets = datasets
        self.xaxis = xaxis
        self.yaxis = yaxis

    @property
    def id(self): return id(self)

    @property
    def cls(self): return 'plot %s' % ''.join([x for x in self.title if x.isalnum()])

    def doit(self, outFile=None, useBokeh=False):
        plot2d.makePlot2d(self.datasets, xAxisSettings=self.xaxis, yAxisSettings=self.yaxis,
                          title=self.title, legendOn=False, outFile=outFile,
                          legendXY=(0.05, 0.95), figsize=(10, 7.5), infill=False, useBokeh=useBokeh)


# -------------------------------------------------------------------------------
# Specific reports
# -------------------------------------------------------------------------------

LEGACYINTERMODE = False


# _______ reports ________

def getEvaluationReport(rs, title=None, saveDoc=False, **otherArgs):
    """
    Evaluation (dict)
    =================

    Data to record in this report ::

        - Author (str)
        - Lab (str)
        - Date (date)
        - MAT (int)
        - Target (str)
        - Projectile (str)
        - Compound nucleus formed (str)
        - Projectile frame (str: "lab" or "cm")
        - Temperature (PQU)

    :param title:
    :param rs:
    :param saveDoc:
    :return:
    """
    endfDoc = rs.styles.getEvaluatedStyle().documentation.endfCompatible.body.split('\n')
    globalMetadata = Report(title=title)
    globalMetadata['Authors'] = endfDoc[0][33:].strip()
    globalMetadata['Lab'] = endfDoc[0][11:22].strip()
    globalMetadata['Evaluation date'] = endfDoc[0][22:33].replace('EVAL-', '').strip()
    globalMetadata['MAT'] = rs.MAT
    globalMetadata['Projectile'] = str(rs.projectile)
    globalMetadata['Target'] = str(rs.target)
    if str(rs.target) != 'n':
        try:
            globalMetadata['Compound nucleus'] = str(rs.compound_nucleus)
        except KeyError:
            pass
    if hasattr(rs, 'attributes'):
        globalMetadata['Attributes'] = rs.attributes
    if saveDoc:
        globalMetadata['ENDF documentation'] = rs.styles.getEvaluatedStyle().documentation.toXMLList()
    globalMetadata['Projectile frame'] = rs.projectileFrame
    globalMetadata['Temperature'] = rs.styles.getEvaluatedStyle().temperature
    globalMetadata['GNDS version'] = rs.formatVersion
    return globalMetadata


def getFissionMetricReport(rs, metricMenu, title="Fission metrics", useCovariance=True, **otherArgs):
    """
    Fissile metrics (dict)
    ======================
        - ETA (PQU)
        - ALF (PQU)

    :param rs:
    :param metricMenu:
    :param title:
    :param useCovariance:
    :return:
    """
    if (hasattr(metricMenu, 'ALF') and metricMenu.ALF) or (hasattr(metricMenu, 'ETA') and metricMenu.ETA):
        try:
            rfis = rs.getReaction('fission')
            rcap = rs.getReaction('capture')
        except KeyError:
            return Report()
        fissionReport = Report(title=title)
        # ALF, ratio of thermal capture to fission
        if hasattr(metricMenu, 'ALF') and metricMenu.ALF:
            fissionReport["ALF"] = computeALPHA(rcap.crossSection, rfis.crossSection,
                                                useCovariance=useCovariance)
        # ETA, nubar_tot*fission cs/absorption cs, not the expanded definition used by EXFOR
        if hasattr(metricMenu, 'ETA') and metricMenu.ETA:
            totalNubar = sum([nubar.multiplicity.toPointwise_withLinearXYs(lowerEps=1e-8, upperEps=1e-8) for nubar in
                              rfis.outputChannel.getProductsWithName('n')])
            fissionReport['ETA'] = computeETA(rcap.crossSection, rfis.crossSection, totalNubar,
                                              useCovariance=useCovariance)
        return {title: fissionReport}
    else:
        return {title: Report()}


def getResonanceReport(rs, title="Resonances", **otherArgs):
    """
    Resonances (dict)
    =================

        - Particle pair list (list)
        - RRR (dict)
        - URR (dict)
        - Channels (list)

    :param rs:
    :param title:
    :return:
    """
    rrReport = Report(title=title)
    rrReport["RRR information"] = getRRRReport(rs, title=None)
    rrReport["URR information"] = getURRReport(rs, title=None)
    rrReport["Particle pairs"] = getParticlePairsReport(rs, title=None)
    return rrReport


def getRRRReport(rs, title="RRR information", **otherArgs):
    """

    RRR (dict)
    ----------

        - Format (str)
            -- approximation (str)
        - Lower bound (PQU)
        - Upper bound (PQU)
        - Scattering length R' (PQU)
        - Lmax (int)
        - No. resonances (int)
        - No. channels (int)

    :param rs:
    :param title:
    :return:
    """
    if hasattr(rs.resonances, 'resolved') and rs.resonances.resolved is not None:
        import fudge.processing.resonances.reconstructResonances as RRReconstruct

        rrrReport = Report(title=title)
        rrrReport['Format'] = rs.resonances.resolved.evaluated.moniker
        if rrrReport['Format'] == "R_Matrix_Limited":
            rrrReport['Format'] += ' (%s approximation)' % rs.resonances.resolved.evaluated.approximation
        rrrReport["Upper bound"] = rs.resonances.resolved.domainMax
        rrrReport['Lower bound'] = rs.resonances.resolved.domainMin
        try:
            resCls = RRReconstruct.getResonanceReconstructionClass(rs.resonances.resolved.evaluated)
            rrr = resCls(rs.resonances.resolved.evaluated, enableAngDists=False, verbose=False)
            rrr.setResonanceParametersByChannel()
            rrrReport["Scattering length (R')"] = PQU.PQU(rrr.getScatteringLength() * 10, "fm")
            rrrReport['Lmax'] = rrr.getLMax()
            rrrReport['No. channels'] = rrr.nChannels
            rrrReport['No. resonances'] = rrr.nResonances
        except ValueError:
            pass
        return rrrReport
    else:
        return Report()


def getURRReport(rs, title="URR information", **otherArgs):
    """
    URR (dict)
    ----------

        - Format (str)
            -- approximation (str)
        - Lower bound (PQU)
        - Upper bound (PQU)
        - No. channels (int)

    :return:
    """
    if hasattr(rs.resonances, 'unresolved') and rs.resonances.unresolved is not None:
        urrReport = Report(title=title)
        urrReport['Format'] = rs.resonances.unresolved.moniker
        urrReport["Upper bound"] = rs.resonances.unresolved.domainMax
        urrReport['Lower bound'] = rs.resonances.unresolved.domainMin
        return urrReport
    else:
        return Report()


def getParticlePairsReport(rs, title="Particle Pairs", **otherArgs):
    """
    Per Particle Pair (list)
    ------------------------
        - name (str)
        - Particle1 (dict)
            - name (str)
            - mass (PQU)
            - spin (str)
            - parity (str: '+' or '-')
            - charge (int)
        - Particle2 (dict)
            - name (str)
            - mass (PQU)
            - spin (str)
            - parity (str: '+' or '-')
            - charge (int)


    :param rs:
    :param title:
    :return:
    """
    import PoPs
    rep = Report(title=title)
    pps = {}
    # Get the particle pairs
    if hasattr(rs.resonances, 'resolved') and rs.resonances.resolved is not None:
        import fudge.processing.resonances.reconstructResonances as RRReconstruct
        resCls = RRReconstruct.getResonanceReconstructionClass(rs.resonances.resolved.evaluated)
        rrr = resCls(rs.resonances.resolved.evaluated, enableAngDists=False, verbose=False)
        try:
            rrr.setResonanceParametersByChannel(warnOnly=True)
        except ValueError as err:
            print(err)
            return rep
        for c in rrr.allChannels:
            if c.reaction not in pps:
                pps[c.reaction] = (c.particleA, c.particleB)
    for pp in pps:
        if pp in ['fission', 'competitive'] or 'fission' in pp:
            continue
        rep[pp] = {}
        spins = rrr.getParticleSpins(pp)
        for ip, p in enumerate(pps[pp]):
            if p == 'gamma':
                particle = rs.getParticle('photon')
            else:
                particle = rs.getParticle(p)
            rep[pp][p] = {}
            if isinstance(particle, PoPs.families.baryon.Particle):
                rep[pp][p]["Mass (amu)"] = particle.mass.float('amu')
                rep[pp][p]["Spin (hbar)"] = particle.spin.float('hbar')
                rep[pp][p]["Parity"] = '+' if str(particle.parity[0]) == "1 " else '-'
                rep[pp][p]["Charge (e)"] = particle.charge.float('e')
            elif isinstance(particle, PoPs.families.gaugeBoson.Particle):
                rep[pp][p]["Mass (amu)"] = particle.mass.float('amu')
                rep[pp][p]["Spin (hbar)"] = particle.spin.float('hbar')
                rep[pp][p]["Parity"] = '+' if str(particle.parity[0]) == "1 " else '-'
                rep[pp][p]["Charge (e)"] = particle.charge.float('e')
            elif isinstance(particle, PoPs.families.nuclide.Particle):
                try:
                    rep[pp][p]["Mass (amu)"] = particle.mass.float('amu')
                except (KeyError, IndexError):
                    rep[pp][p]["Mass (amu)"] = 'Unknown'
                try:
                    rep[pp][p]["Spin (hbar)"] = particle.nucleus.spin.float('hbar')
                except (KeyError, IndexError):
                    rep[pp][p]["Spin (hbar)"] = 'Unknown'
                try:
                    rep[pp][p]["Parity"] = str(particle.nucleus.parity[0])
                except (KeyError, IndexError):
                    rep[pp][p]["Parity"] = 'Unknown'
                try:
                    rep[pp][p]["Charge (e)"] = particle.nucleus.charge.float('e')
                except (KeyError, IndexError):
                    rep[pp][p]["Charge (e)"] = 'Unknown'
            else:
                raise TypeError("Unknown particle type %s" % str(particle.__class__))
    return rep


# _______ plots ________

def getChannelPlots(rs, metricMenu, title="Channel plots", channelMatchFunction=lambda c: True, verbose=False, justTable=True):
    """
        - Excitation functions
            -- Staircase: Cum. no. resonances vs. E (plot)
            -- <D> vs. E (plot)
            -- <Gamma> vs. E (plot)
            -- Transmission coefficients: Tc vs. E (plot)
            -- Scattering radii: Rc vs. E (plot)
        - Overridden defaults
            -- Penetribility: Pc vs. E (plot)
            -- (SKIP) Shift: Sc vs. E (plot)
            -- Phase: Phic vs. E (plot)
            -- Background cross section present? (bool)
        - GOE diagnostics
            -- Delta3 (plot)
            -- Porter-Thomas, nu fixed, fit (Gamma/<Gamma>)_min (plot)
            -- Porter-Thomas, (Gamma/<Gamma>)_min fixed, fit nu (plot)
            -- spacing distribution (plot)

    :param rs:
    :param metricMenu:
    :param title:
    :param channelMatchFunction:
    :param verbose:
    :param justTable:
    :return:
    """
    # ---------------------------------
    # Process the resonances
    # ---------------------------------
    import fudge.processing.resonances.reconstructResonances as RRReconstruct

    # get URR parameters
    if hasattr(rs.resonances, 'unresolved') and rs.resonances.unresolved is not None:
        urr = RRReconstruct.URRcrossSection(rs)
        egrid, flag = urr.generateEnergyGrid()
        urrAveQuant = urr.getWidthsAndSpacings()  # egrid, interpolateWidths=flag)
    else:
        urr = None
        urrAveQuant = None

    # get RRR parameter averages
    if hasattr(rs.resonances, 'resolved') and rs.resonances.resolved is not None:
        resCls = RRReconstruct.getResonanceReconstructionClass(rs.resonances.resolved.evaluated)
        rrr = resCls(rs.resonances.resolved.evaluated, enableAngDists=False, verbose=verbose)
        rrr.setResonanceParametersByChannel()
        rrrAveQuant = rrr.getAverageQuantities()  # computeUncertainty=True, nBins=3) #resonancesPerBin=50)
    else:
        rrr = None
        rrrAveQuant = None

    def get_key_from_channel(_c):
        if int(_c.s) == _c.s:
            s = _c.s / 2.0
        else:
            s = _c.s
        return '%s (j=%s,l=%s,s=%s)' % (str(_c.reaction), str(_c.J), str(_c.l), str(s))

    # Excitation functions...
    if hasattr(metricMenu, 'aveSpacingPlot') and metricMenu.aveSpacingPlot:
        if verbose:
            print(winged_banner('aveSpacingPlot'))
        if rrr is not None:
            rrrAveWidth = rrrAveQuant['spacings']
            rrrEgrid = rrrAveQuant['energyGrid']
        if urr is not None:
            raise NotImplementedError("urr data for aveSpacingPlot: FIXME")
        raise NotImplementedError("aveSpacingPlot: FIXME")

    if hasattr(metricMenu, 'aveWidthPlot') and metricMenu.aveWidthPlot:
        if verbose:
            print(winged_banner('aveWidthPlot'))
        if rrr is not None:
            rrrAveWidth = rrrAveQuant['widths']
            rrrEgrid = rrrAveQuant['energyGrid']
        if urr is not None:
            raise NotImplementedError("urr data for aveWidthPlot: FIXME")
        raise NotImplementedError("aveWidthPlot: FIXME")

    if hasattr(metricMenu, 'scatteringRadiusPlot') and metricMenu.scatteringRadiusPlot:
        if verbose:
            print(winged_banner('scatteringRadiusPlot'))
        raise NotImplementedError("scatteringRadiusPlot: FIXME")

    if hasattr(metricMenu, 'transmissionCoeff') and metricMenu.transmissionCoeff:
        if verbose:
            print(winged_banner("transmissionCoeff"))
        if rrr is not None:
            rrrTc = rrr.getTransmissionCoefficientsFromSumRule(computeUncertainty=True)
        else:
            rrrTc = None
        if urr is not None:
            urrTc = urr.getTransmissionCoefficientsFromSumRule()
        else:
            urrTc = None
        for c in rrrTc:
            ds = []
            currm = None
            for curr in urrTc:
                if curr.J == c.J and curr.l == c.l:
                    if curr.reaction == 'elastic' and 'n +' in c.reaction:
                        currm = curr
                        break
                    if curr.reaction == 'capture' and '+ photon' in c.reaction:
                        currm = curr
                        break
            if rrrTc is not None:
                ds.append(rrrTc[c])
            if urrTc is not None and currm is not None:
                ds.append(urrTc[currm])
            if urrTc is None and rrrTc is None:
                p = None
            elif justTable:
                if currm is not None and currm.reaction == "elastic":
                    print(c)
                    print('RRR')
                    for iE in range(len(ds[0]['Tc'])):
                        print(ds[0]['energyGrid'][iE], ds[0]['Tc'][iE])
                    print('URR')
                    for x in ds[1]:
                        print(' '.join(map(str, x)))
            else:
                p = Plot(title="T_c for %s channel" % get_key_from_channel(c), datasets=ds)
                p.doit()
    #        print( NotImplementedError("transmissionCoeff: FIXME") )

    if hasattr(metricMenu, 'staircasePlot') and metricMenu.staircasePlot:
        if verbose:
            print(winged_banner("staircasePlot"))
        if rrr is not None:
            for c in rrr.channels:
                if not channelMatchFunction(c):
                    continue
                ERList = [rrr._energies[iR] for iR in rrr.channels[c]]
                if len(ERList) < 2:
                    print("Channel %c doesn't have enough resonances to generate a meaningful staircase plot" %
                          get_key_from_channel(c))
                    continue
                ERList.sort()
                theData = [[ERList[0], 1]]
                for thisEnergy in ERList[1:]:
                    if abs(thisEnergy - theData[-1][0]) < 1e-6 * thisEnergy:
                        theData[-1][1] = theData[-1][1] + 1  # Deal with corner case where two levels are degenerate
                    else:
                        theData.append([thisEnergy, theData[-1][1] + 1])
                theAxes = XYs1d.defaultAxes(labelsUnits={
                    XYs.yAxisIndex: ('No. levels', ''),
                    XYs.xAxisIndex: ('Ex', 'eV')})
                p = Plot(
                    title="Cummulative level distribution for %s channel" % get_key_from_channel(c),
                    datasets=[XYs.XYs1d(data=theData, axes=theAxes, interpolation='flat')])
                p.doit()

    # GOE related diagnostics...
    if hasattr(metricMenu, 'PorterThomas') and metricMenu.PorterThomas:
        if rrr is not None:
            rrrResults = rrr.getPorterThomasFitToWidths()
            for c in rrrResults:
                if not channelMatchFunction(c):
                    continue
                if verbose:
                    print(winged_banner("PorterThomas for %s" % get_key_from_channel(c)))
                print(rrrResults[c]['bin_edges'], rrrResults[c]['hist'])
                PD = grouped_values_to_XYs(rrrResults[c]['bin_edges'], rrrResults[c]['hist'],
                                           domainUnit='', domainName='D/<D>', rangeUnit='', rangeName='P(D/<D>)',
                                           accuracy=upperEps, domainMin=0.0, domainMax=8.0)
                p = Plot(title="Porter-Thomas analysis for %s channel" % get_key_from_channel(c), datasets=[PD])
                p.doit()
        print(NotImplementedError("PorterThomas: FIXME"))

    if hasattr(metricMenu, 'spacingDistributionPlot') and metricMenu.spacingDistributionPlot:
        if verbose:
            print(winged_banner("spacingDistributionPlot"))
        print(NotImplementedError("spacingDistributionPlot: FIXME"))

    if hasattr(metricMenu, 'DysonMehtaDelta3') and metricMenu.DysonMehtaDelta3:
        if verbose:
            print(winged_banner("DysonMehtaDelta3"))
        print(NotImplementedError("DysonMehtaDelta3: FIXME"))

    # Overriding defaults...
    if hasattr(metricMenu, 'penetrabilityPlot') and metricMenu.penetrabilityPlot:
        if rrr is not None:
            for c in rrr.channels:
                if not channelMatchFunction(c):
                    continue
                if verbose:
                    print(winged_banner("penetrabilityPlot for %s" % get_key_from_channel(c)))
                Pc = function_to_XYs(func=lambda E, x: rrr.penetrationFactorByChannel(c, E), fpars={}, rangeUnit='',
                                     rangeName='Channel penetrability, P_c(E)')
                p = Plot(title="Penetrability %s channel" % get_key_from_channel(c), datasets=[Pc])
                p.doit()

    if hasattr(metricMenu, 'phasePlot') and metricMenu.phasePlot:
        if rrr is not None:
            for c in rrr.channels:
                if not channelMatchFunction(c):
                    continue
                if verbose:
                    print(winged_banner("phasePlot for %s" % get_key_from_channel(c)))
                phic = function_to_XYs(func=lambda E, x: rrr.phiByChannel(c, E), fpars={}, rangeUnit='',
                                       rangeName='Channel phase, phi_c(E)')
                p = Plot(title="Phase for %s channel" % get_key_from_channel(c), datasets=[phic])
                p.doit()

    if hasattr(metricMenu, 'shiftPlot') and metricMenu.shiftPlot:
        if rrr is not None:
            for c in rrr.channels:
                if not channelMatchFunction(c):
                    continue
                if verbose:
                    print(winged_banner("shiftPlot for %s" % get_key_from_channel(c)))
                phic = function_to_XYs(func=lambda E, x: rrr.shiftFactorByChannel(c, E), fpars={}, rangeUnit='',
                                       rangeName='Channel shift factor, phi_c(E)')
                p = Plot(title="Shift for %s channel" % get_key_from_channel(c), datasets=[phic])
                p.doit()

    if hasattr(metricMenu, 'backgroundXSPlot') and metricMenu.backgroundXSPlot:
        if verbose:
            print(winged_banner("backgroundXSPlot"))
        print(NotImplementedError("backgroundXSPlot: FIXME"))


# _______ datatables ________

def getChannelDataTable(rs, metricMenu, title="Channels", doURR=False, verbose=False, **otherArgs):
    """
    Per Channel (dict)
    ------------------

        - c (J,L,S,Pi,PP) (struct)
        - MT (int)
        - Reaction string (str)
        - Q (PQU)
        - Eth = Xi (PQU)
        - nu (float)
        - Eliminated? (bool)
        - Competative? (bool)
        - Relativisitic kinematics? (bool)
        - For potential scattering only? (bool)
        - No. resonances (dict)
            -- RRR (int)
            -- estimated URR (float)
        - List of resonances with ER < 0 (list)
        - Pole strength (PQU)

    :param rs:
    :param title:
    :return:
    """
    # Final results in a table for later report building
    columnHeaders = collections.OrderedDict()  # columnHeader instances
    data = collections.OrderedDict()  # nested list with the actual data!
    icol = 0

    def addColumn(columns, icol, key, unit=None):
        columns[key] = tableModule.ColumnHeader(icol, key, unit=unit)
        icol += 1

    def updateColumnsAndRow(columns, row, icol, key, val, unit=None):
        if key not in columns:
            addColumn(columns, icol, key, unit=unit)
        row.append(val)

    # ---------------------------------
    # Process the resonances
    # ---------------------------------
    import fudge.processing.resonances.reconstructResonances as RRReconstruct

    # get URR parameters
    if hasattr(rs.resonances, 'unresolved') and rs.resonances.unresolved is not None:
        urr = RRReconstruct.URRcrossSection(rs.resonances.unresolved.evaluated)
        egrid, flag = urr.generateEnergyGrid()
        urr.getWidthsAndSpacings()  # egrid, interpolateWidths=flag)

    # get RRR parameter averages
    if hasattr(rs.resonances, 'resolved') and rs.resonances.resolved is not None:
        resCls = RRReconstruct.getResonanceReconstructionClass(rs.resonances.resolved.evaluated)
        rrr = resCls(rs.resonances.resolved.evaluated, enableAngDists=False, verbose=verbose)
        rrr.setResonanceParametersByChannel(warnOnly=True)
        rrrAveQuant = rrr.getAverageQuantities(computeUncertainty=True)
    else:
        rrr = None
        rrrAveQuant = None

    def get_key_from_channel(c):
        if int(c.s) == c.s:
            s = c.s / 2.0
        else:
            s = c.s
        return '%s (j=%s,l=%s,s=%s)' % (str(c.reaction), str(float(c.J)), str(c.l), str(float(s)))

    reactions = []

    # Initialize the rows
    if rrr is not None:
        for c in rrr.allChannels:
            if c.reaction not in reactions:
                reactions.append(c.reaction)
            key = get_key_from_channel(c)
            if key not in data: data[key] = []
    if doURR and urr is not None:
        for reaction in reactions:
            for lj in urr.averageWidthFuncs:
                key = '%s (j=%s,l=%s,s=%s)' % (str(reaction), str(float(lj[1])), str(lj[0]), '0.5')
                if key not in data:
                    data[key] = []

    # Compute the simple stuff
    if rrr is not None:
        for c in rrr.allChannels:
            key = get_key_from_channel(c)
            numResonances = len([x for x in rrr.allChannels[c] if rrr.allChannels[c][x] != 0.0])
            numNegEnRes = sum([rrr._energies[iER] < 0 for iER in rrr.allChannels[c] if rrr.allChannels[c][iER] != 0.0])
            updateColumnsAndRow(columnHeaders, data[key], icol, 'No. resonances', numResonances, '')
            updateColumnsAndRow(columnHeaders, data[key], icol, 'No. resonances w/ ER<0', numNegEnRes, '')
            if False:
                updateColumnsAndRow(columnHeaders, data[key], icol, 'MT', -13, '')  # FIXME
            updateColumnsAndRow(columnHeaders, data[key], icol, 'gfact', c.gfact, '')
            if False:
                updateColumnsAndRow(columnHeaders, data[key], icol, 'Q', -13, 'eV')  # FIXME
            updateColumnsAndRow(columnHeaders, data[key], icol, 'Threshold E', c.Xi, 'eV')
            updateColumnsAndRow(columnHeaders, data[key], icol, 'Eliminated?', str(c.eliminated), '')
            updateColumnsAndRow(columnHeaders, data[key], icol, 'Competative?',
                                str(c.channelClass == RRReconstruct.COMPETITIVECHANNEL), '')
            updateColumnsAndRow(columnHeaders, data[key], icol, 'Relativistic?', str(c.useRelativistic), '')
            if False:
                updateColumnsAndRow(columnHeaders, data[key], icol, 'Has background xs?', str(False), '')  # FIXME
            if False:
                updateColumnsAndRow(columnHeaders, data[key], icol, 'B_c', str(False), '')  # FIXME
            if not c.eliminated:
                updateColumnsAndRow(columnHeaders, data[key], icol, 'Pot. scatt. only?',
                                    str(len([x for x in rrr.allChannels[c] if rrr.allChannels[c][x] != 0.0]) == 0), '')
            else:
                updateColumnsAndRow(columnHeaders, data[key], icol, 'Pot. scatt. only?', str(False), '')
        for c in rrrAveQuant:
            key = get_key_from_channel(c)
            #           print( rrrAveQuant[c]['spacings'] )
            if rrrAveQuant[c]['spacings']:
                updateColumnsAndRow(columnHeaders, data[key], icol, 'RRR <D>', str(rrrAveQuant[c]['spacings'][0]), 'eV')
            if rrrAveQuant[c]['widths']:
                updateColumnsAndRow(columnHeaders, data[key], icol, 'RRR <Gamma>', str(rrrAveQuant[c]['widths'][0]), 'eV')
    if doURR and urr is not None:
        raise NotImplementedError("TODO: URR average extraction")
        for reaction in reactions:
            for lj in urr.averageWidthFuncs:
                key = '%s (j=%s,l=%s,s=%s)' % (str(reaction), str(float(lj[1])), str(lj[0]), '0.5')
                try:
                    updateColumnsAndRow(columnHeaders, data[key], icol, 'URR <D>', urr.levelSpacings[lj][str(reaction)],
                                        'eV')
                except TypeError:
                    pass
                try:
                    updateColumnsAndRow(columnHeaders, data[key], icol, 'URR <Gamma>',
                                        urr.averageWidths[lj][str(reaction)], 'eV')
                except TypeError:
                    pass

    # Tough to extract quantities
    if (hasattr(metricMenu, 'effectiveDOF') and metricMenu.effectiveDOF):
        if rrr is not None:
            try:
                rrrResults = rrr.getPorterThomasFitToWidths()
            except AttributeError:
                print("Are you using an old version of SciPy?  scipy.optimize.curve_fit is missing.")
                rrrResults = []
            for c in rrrResults:
                nu = PQU.PQU(rrrResults[c]['dof'], uncertainty=rrrResults[c]['ddof'])
                key = get_key_from_channel(c)
                updateColumnsAndRow(columnHeaders, data[key], icol, 'RRR DOF', str(nu), '')
        if doURR and urr is not None:
            raise NotImplementedError("TODO: URR average extraction")
            for reaction in reactions:
                for lj in urr.DOFs:
                    key = '%s (j=%s,l=%s,s=%s)' % (str(reaction), str(float(lj[1])), str(lj[0]), '0.5')
                    try:
                        updateColumnsAndRow(columnHeaders, data[key], icol, 'URR DOF', urr.DOFs[lj][str(reaction)], '')
                    except KeyError:
                        pass
    if (hasattr(metricMenu, 'poleStrength') and metricMenu.poleStrength):
        if rrr is not None:
            sc = rrr.getPoleStrength(computeUncertainty=True)
            for c in sc:
                key = get_key_from_channel(c)
                updateColumnsAndRow(columnHeaders, data[key], icol, 'RRR sc', str(sc[c]['sc'][0]), '')

    nCols = len(list(columnHeaders.values()))
    tmpData = []
    for irow, row in enumerate(data.values()):
        if row == []:
            tmpData.append([0 for i in range(nCols)])
        elif len(row) < nCols:
            for icell in range(len(row), nCols): row.append(0)
            tmpData.append(row)
        else:
            tmpData.append(row)
    result = DataTable(data=tmpData, columns=list(columnHeaders.values()), rows=list(data.keys()))
    result.title = title
    return result


def getReactionDataTable(rs, metricMenu, title="Reactions", MTList=[], useCovariance=True, **otherArgs):
    """
    Per Reaction (datatable)
    ========================

        - MT (int)
        - Reaction string (str)
        - No. products (int)
        - Residual (str)
        - Q (PQU)
        - Eth == Xi (PQU)
        - Reaction rates
            -- Godiva (PQU)
            -- Flattop (PQU)
            -- Jezebel (PQU)
            -- BigTen (PQU)
            -- FUND-IPPE (PQU)
            -- Cf252 (analytic) (PQU)
            -- Cf252 (PQU)
        - Engineering metrics
            -- RI (PQU)
            -- Westcott factor (PQU)
            -- Thermal xs (PQU)
            -- E14 (PQU)
        - Astrophysics metrics
            -- MACS(30 keV) (PQU)
            -- ARR (PQU)
        - INTER metrics
            -- RI (PQU)
            -- Westcott factor (PQU)
            -- Thermal xs (PQU)
            -- Cf252 (analytic) (PQU)
            -- Cf252 (PQU)
            -- MACS(30 keV) (PQU)
            -- E14 (PQU)

    :param rs:
    :param metricMenu:
    :param title:
    :param useCovariance:
    :param MTList:
    :return:
    """
    # Final results in a table for later report building
    columnHeaders = collections.OrderedDict()  # columnHeader instances
    rowNames = []  # really just row names
    data = []  # nested list with the actual data!
    icol = 0

    def addColumn(columns, icol, key, unit=None):
        columns[key] = tableModule.ColumnHeader(icol, key, unit=unit)
        icol += 1

    def updateColumnsAndRow(columns, row, icol, key, val, unit=None):
        if key not in columns:
            addColumn(columns, icol, key, unit=unit)
        row.append(val)

    import fudge.sums

    # ---------------------------------
    #  Process each requested reaction
    # ---------------------------------
    for r in rs:

        # Get rid of all reactions that are not plain reactions or sums of reactions
        # (e.g. no fission components of production stuff)
        if not isinstance(r, (fudge.reactions.reaction.Reaction, fudge.sums.CrossSectionSum)):
            continue

        # Keep only the MT's we actually want to test
        mt = r.ENDF_MT
        if mt not in MTList:
            continue

        # Initialize this row
        rowNames.append(str(r))
        row = []

        # Check if we need to add the default columns
        updateColumnsAndRow(columnHeaders, row, icol, 'MT', mt)

        # Set Q
        if hasattr(r, 'getQ'):
            Q = r.getQ(unit='eV')
        else:
            Q = 0.0
        updateColumnsAndRow(columnHeaders, row, icol, 'Q', Q, 'eV')

        # Set EThreshold
        if hasattr(r, "getThreshold"):
            EThreshold = r.getThreshold(unit='eV')
        else:
            EThreshold = 0.0
        updateColumnsAndRow(columnHeaders, row, icol, 'Threshold', EThreshold, 'eV')

        # Resonance Integral
        if hasattr(metricMenu, 'RI') and metricMenu.RI:
            if EThreshold < 0.0235339347844:
                RI = computeRI(r.crossSection, useCovariance=useCovariance, Ecut=otherArgs.get('RI_Ecut', None))
                updateColumnsAndRow(columnHeaders, row, icol, 'RI', RI, RI.unit)
            else:
                row.append(tableModule.Blank())

        # 14 MeV Point
        if hasattr(metricMenu, 'fourteen') and metricMenu.fourteen:
            if LEGACYINTERMODE:
                # Legacy INTER uses 14.0 for the peak of the d-t outgoing neutron spectrum
                E14 = PQU.PQU(14., 'MeV')
            else:
                E14 = PQU.PQU(14.2, 'MeV')
            fourteen = computeFourteenMeVPoint(r.crossSection, E14, useCovariance=useCovariance)
            updateColumnsAndRow(columnHeaders, row, icol, str(E14), fourteen, fourteen.unit)

        # Window average
        if hasattr(metricMenu, 'windowCenter') and hasattr(metricMenu, 'windowWidth') and \
                metricMenu.windowCenter and metricMenu.windowWidth:
            window = computeAverageOverWindow(r.crossSection, PQU.PQU(metricMenu.windowCenter, unit='keV'),
                                              PQU.PQU(metricMenu.windowWidth, unit='keV'), useCovariance=useCovariance)
            updateColumnsAndRow(columnHeaders, row, icol, 'Ave. over %s+/-%s keV' %
                                (str(metricMenu.windowCenter), str(metricMenu.windowWidth / 2.0)), window, window.unit)

        # User defined spectrum average
        if hasattr(metricMenu, 'userDefinedSpecAve') and metricMenu.userDefinedSpecAve:
            path = '/'.join(metricMenu.userDefinedSpecAve.split('/')[:-1])
            key = metricMenu.userDefinedSpecAve.split('/')[-1]
            ave = computeUserDefinedSpectrumAve(r.crossSection, path, key, useCovariance=useCovariance)
            updateColumnsAndRow(columnHeaders, row, icol, key, ave, ave.unit)

        # Lookup spectrum average
        if hasattr(metricMenu, 'otherSpecAve') and metricMenu.otherSpecAve and \
                metricMenu.otherSpecAve not in ['Godiva', 'Jezebel', 'BigTen']:
            ave = computeLookupSpectrumAve(r.crossSection, metricMenu.otherSpecAve, useCovariance=useCovariance)
            updateColumnsAndRow(columnHeaders, row, icol, metricMenu.otherSpecAve, ave, ave.unit)

        # Room temperature cross section
        if hasattr(metricMenu, 'thermal') and metricMenu.thermal:
            if EThreshold < 0.0235339347844:
                xsRT = computeRoomTempCS(r.crossSection, useCovariance=useCovariance)
                updateColumnsAndRow(columnHeaders, row, icol, "2200 m/s", xsRT, xsRT.unit)
            else:
                row.append(tableModule.Blank())

        # Westcott factor
        if hasattr(metricMenu, 'Westcott') and metricMenu.Westcott:
            if EThreshold < 0.0235339347844 and computeRoomTempCS(r.crossSection, useCovariance=False).value != 0.0:
                if LEGACYINTERMODE:
                    # Legacy INTER ignores the recoil of the target,
                    # so a = 1.0 rather than a = mTarget/(mTarget+mProjectile)
                    westcott = computeWestcottFactor(r.crossSection, a=PQU.PQU(1., ""), useCovariance=useCovariance)
                else:
                    westcott = computeWestcottFactor(r.crossSection, useCovariance=useCovariance)
                updateColumnsAndRow(columnHeaders, row, icol, "Westcott factor", westcott, westcott.unit)
            else:
                row.append(tableModule.Blank())

        # Godiva spectrum average
        if hasattr(metricMenu, 'Godiva') and metricMenu.Godiva:
            if EThreshold < 5e6:
                x = computeGodivaSpectrumAve(r.crossSection, useCovariance=useCovariance)
                updateColumnsAndRow(columnHeaders, row, icol, "Godiva", x, x.unit)
            else:
                row.append(tableModule.Blank())

        # Jezebel spectrum average
        if hasattr(metricMenu, 'Jezebel') and metricMenu.Jezebel:
            if EThreshold < 5e6:
                x = computeJezebelSpectrumAve(r.crossSection, useCovariance=useCovariance)
                updateColumnsAndRow(columnHeaders, row, icol, "Jezebel", x, x.unit)
            else:
                row.append(tableModule.Blank())

        # BigTen spectrum average
        if hasattr(metricMenu, 'BigTen') and metricMenu.BigTen:
            if EThreshold < 5e6:
                x = computeBigTenSpectrumAve(r.crossSection, useCovariance=useCovariance)
                updateColumnsAndRow(columnHeaders, row, icol, "BigTen", x, x.unit)
            else:
                row.append(tableModule.Blank())

        # FUND-IPPE spectrum average
        if hasattr(metricMenu, 'FUNDIPPE') and metricMenu.FUNDIPPE:
            raise NotImplementedError("FIXME")
            if EThreshold < 5e6:
                x = computeFUNDIPPESpectrumAve(r.crossSection, useCovariance=useCovariance)
                updateColumnsAndRow(columnHeaders, row, icol, "FUND-IPPE", x, x.unit)
            else:
                row.append(tableModule.Blank())

        # 252Cf spectrum average (using analytic approximation from INTER)
        if hasattr(metricMenu, 'CfSpectAnalytic') and metricMenu.CfSpectAnalytic:
            if EThreshold < 5e6:
                x = computeCf252SpectrumAveAnalytic(r.crossSection, useCovariance=useCovariance)
                updateColumnsAndRow(columnHeaders, row, icol, "252Cf (analytic)", x, x.unit)
            else:
                row.append(tableModule.Blank())

        # 252Cf spectrum average
        if hasattr(metricMenu, 'CfSpectENDF') and metricMenu.CfSpectENDF:
            if EThreshold < 5e6:
                x = computeCf252SpectrumAve(r.crossSection, useCovariance=useCovariance)
                updateColumnsAndRow(columnHeaders, row, icol, "252Cf", x, x.unit)
            else:
                row.append(tableModule.Blank())

        # MACS
        if hasattr(metricMenu, 'MACS') and metricMenu.MACS:
            if metricMenu.MACS:
                if isinstance(metricMenu.MACS, PQU.PQU):
                    kT = metricMenu.MACS
                else:
                    kT = PQU.PQU(metricMenu.MACS, 'keV')
            else:
                kT = PQU.PQU(30., 'keV')
            if EThreshold / 10 < kT.getValueAs('eV'):
                x = computeMACS(r.crossSection, T=kT, useCovariance=useCovariance).inUnitsOf('mb')  # pref. unit is mb
                updateColumnsAndRow(columnHeaders, row, icol, "MACS(%s)" % str(kT), x, 'mb')
            else:
                row.append(tableModule.Blank())

        # ARR
        if hasattr(metricMenu, 'ARR') and metricMenu.ARR:
            if metricMenu.ARR:
                if isinstance(metricMenu.ARR, PQU.PQU):
                    kT = metricMenu.ARR
                else:
                    kT = PQU.PQU(metricMenu.ARR, 'keV')
            else:
                kT = PQU.PQU(30., 'keV')
            if EThreshold / 10 < kT.getValueAs('eV'):
                x = computeAstrophysicalReactionRate(r.crossSection, T=kT, useCovariance=useCovariance)
                updateColumnsAndRow(columnHeaders, row, icol, "ARR(%s)" % str(kT), x, x.unit)
            else:
                row.append(tableModule.Blank())

        data.append(row)

    result = DataTable(data=data, columns=list(columnHeaders.values()), rows=rowNames)
    result.title = title
    return result


def getReactionEngineeringDataTable(rs, title="Engineering metrics", MTList=None, useCovariance=True, **otherArgs):
    """
    Engineering metrics ::

        - RI (PQU)
        - Westcott factor (PQU)
        - Thermal xs (PQU)
        - E14 (PQU)

    :param rs:
    :param title:
    :param MTList:
    :return:
    """
    if MTList is None:
        MTList = []
    metricMenu = collections.namedtuple('metricMenu',
                                        "RI Westcott thermal fourteen")  # emulates results from ArgumentParser.parse_args() function
    results = Report()
    results[title] = getReactionDataTable(rs,
                                          metricMenu=metricMenu(RI=True, Westcott=True, thermal=True, fourteen=True),
                                          title=title, MTList=MTList, useCovariance=useCovariance)
    return results


def getReactionAstrophysicsDataTable(rs, title="Astrophysics metrics", MTList=None, useCovariance=True, **otherArgs):
    """
    Astrophysics metrics ::
        - MACS(30 keV) (PQU)
        - ARR (PQU)

    :param rs:
    :param title:
    :param MTList:
    :return:
    """
    if MTList is None:
        MTList = []

    metricMenu = collections.namedtuple('metricMenu',
                                        "MACS ARR")  # emulates results from ArgumentParser.parse_args() function
    results = Report()
    results[title] = getReactionDataTable(rs,
                                          metricMenu=metricMenu(MACS=PQU.PQU(30, 'keV'), ARR=PQU.PQU(30, 'keV')),
                                          title=title, MTList=MTList, useCovariance=useCovariance)
    return results


def getReactionLegacyINTERDataTable(rs, title="Legacy INTER metrics", MTList=(1, 2, 18, 102, 16), useCovariance=True, **otherArgs):
    """
    INTER metrics ::

        - RI (PQU)
        - Westcott factor (PQU)
        - Thermal xs (PQU)
        - Cf252 (analytic) (PQU)
        - Cf252 (PQU)
        - MACS(30 keV) (PQU)
        - E14 (PQU)

    :param rs:
    :param title:
    :param MTList:
    :return:
    """
    metricMenu = collections.namedtuple('metricMenu',
                                        "RI Westcott thermal CfSpectAnalytic CfSpectENDF MACS fourteen")  # emulates results from ArgumentParser.parse_args() function
    results = Report()
    results[title] = getReactionDataTable(rs, metricMenu=metricMenu(RI=True, Westcott=True, thermal=True,
                                                                    CfSpectAnalytic=True, CfSpectENDF=True,
                                                                    MACS=PQU.PQU(30., 'keV'), fourteen=True),
                                          title=title, MTList=MTList, useCovariance=useCovariance)
    return results


def getReactionAIntegralDataTable(rs, title="Integral metrics", MTList=None, useCovariance=True, **otherArgs):
    """
    Reaction rates ::

        - Godiva (PQU)
        - Flattop (PQU)
        - Jezebel (PQU)
        - BigTen (PQU)
        - FUND-IPPE (PQU)
        - Cf252 (analytic) (PQU)
        - Cf252 (PQU)

    :param rs:
    :param title:
    :param MTList:
    :return:
    """
    if MTList is None:
        MTList = []

    metricMenu = collections.namedtuple('metricMenu',
                                        "Godiva Flattop Jezebel BigTen FUNDIPPE CfSpectAnalytic CfSpectENDF")  # emulates results from ArgumentParser.parse_args() function
    results = Report()
    results[title] = getReactionDataTable(rs,
                                          metricMenu=metricMenu(Godiva=True, Flattop=False, Jezebel=True, BigTen=True,
                                                                FUNDIPPE=False, CfSpectAnalytic=True, CfSpectENDF=True),
                                          title=title, MTList=MTList, useCovariance=useCovariance)
    return results


# -------------------------------------------------------------------------------
# obsolete code
# -------------------------------------------------------------------------------

class OtherReport(Report):

    def html_report(self, args):
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
        if True:
            report = [line.replace('Gamma', SYMBOL('Gamma')) for line in report]
        return report
