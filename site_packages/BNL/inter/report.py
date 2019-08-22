import json
from fudge.vis.matplotlib import plot2d
from fudge.core.utilities.brb import banner, winged_banner
from xData.table import blank
from datatables import *
from fudge.gnds.physicalQuantity import temperature
from html import *
from metrics import *
from utils import *

'''
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

'''

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
         if isinstance(obj, temperature):
             return str(obj)
         if isinstance(obj, DataTable):
             return '\n'.join(obj.toStringList())
         if isinstance(obj, datatables.DataTable):
             return obj.to_string()
         # Let the base class default method raise the TypeError
         return json.JSONEncoder.default(self, obj)


class Report( collections.OrderedDict ):

    def __init__(self, title=None, indent=4, *args, **kwds):
        self.title = title
        self.indent = indent
        collections.OrderedDict.__init__(self, *args, **kwds)

    @property
    def id(self): return id(self)

    @property
    def cls(self): return 'report %s'%''.join([x for x in self.title if x.isalnum()])

    def text_title(self): return banner(self.title)

    def text_report(self, *args, **kwds):
        return self.json_report(*args, indent=self.indent, separators=(',', ': '), **kwds)

    def json_report(self, *args, **kwds):
        if self.title is None: return json.dumps(self, *args, cls=ComplexEncoder, **kwds)
        return json.dumps({self.title:self}, *args, cls=ComplexEncoder, **kwds)

    def html_report(self, *args, **kwds):
        #FIXME: don't work yet
        report = [H1(self.title)]
        report.append('\n'.join(OBJECT2HTML(self)))
        if True:
            report = [line.replace('Gamma', SYMBOL('Gamma')) for line in report]
        return report


class Plot():

    def __init__(self, title=None, datasets=[], xaxis=None, yaxis=None):
        self.title=title
        self.datasets=datasets
        self.xaxis=xaxis
        self.yaxis=yaxis

    @property
    def id(self): return id(self)

    @property
    def cls(self): return 'plot %s' % ''.join([x for x in self.title if x.isalnum()])

    def doit(self,outFile=None,useBokeh=False):
        plot2d.makePlot2d(self.datasets, xAxisSettings=self.xaxis, yAxisSettings=self.yaxis,
                          title=self.title, legendOn=False, outFile=outFile,
                          legendXY=(0.05, 0.95), figsize=(10, 7.5), infill=False, useBokeh=useBokeh)


# -------------------------------------------------------------------------------
# Specific reports
# -------------------------------------------------------------------------------

LEGACYINTERMODE=False

# _______ reports ________

def getEvaluationReport( rs, title=None, saveDoc=False ):
    '''
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

    :param endfFile:
    :param rs:
    :param saveDoc:
    :param doALF:
    :param doETA:
    :param useCovariance:
    :return:
    '''
    endfDoc=rs.documentation['endfDoc'].getLines()
    globalMetadata = Report(title=title)
    globalMetadata['Authors'] = endfDoc[0][33:].strip()
    globalMetadata['Lab'] = endfDoc[0][11:22].strip()
    globalMetadata['Evaluation date'] = endfDoc[0][22:33].replace('EVAL-','').strip()
    globalMetadata['MAT'] = rs.MAT
    globalMetadata['Projectile'] = str(rs.projectile)
    globalMetadata['Target'] = str(rs.target)
    if str(rs.target) != 'n':
        try: globalMetadata['Compound nucleus'] = str(rs.compound_nucleus)
        except KeyError: pass
    if hasattr(rs, 'attributes'):  globalMetadata['Attributes'] = rs.attributes
    if saveDoc: globalMetadata['ENDF documentation'] = rs.documentation['endfDoc']
    globalMetadata['Projectile frame'] = rs.getProjectileFrame()
    globalMetadata['Temperature'] = rs.getTemperature('eval')
    globalMetadata['GNDS version'] = rs.GNDS_version
    return globalMetadata


def getFissionMetricReport( rs, metricMenu, title="Fission metrics", useCovariance=True ):
    '''
    Fissile metrics (dict)
    ======================
        - ETA (PQU)
        - ALF (PQU)

    :param rs:
    :param metricMenu:
    :param title:
    :param useCovariance:
    :return:
    '''
    if ( hasattr(metricMenu, 'ALF') and metricMenu.ALF ) or ( hasattr(metricMenu, 'ETA') and metricMenu.ETA ):
        try:
            rfis = rs.getReaction('fission')
            rcap = rs.getReaction('capture')
        except KeyError:
            return Report()
        fissionReport=Report(title=title)
        # ALF, ratio of thermal capture to fission
        if ( hasattr(metricMenu, 'ALF') and metricMenu.ALF ):
            fissionReport["ALF"] = computeALPHA(rcap.crossSection, rfis.crossSection,
                                                     useCovariance=useCovariance)
        # ETA, nubar_tot*fission cs/absorption cs, not the expanded definition used by EXFOR
        if ( hasattr(metricMenu, 'ETA') and metricMenu.ETA ):
            totalNubar = sum([nubar.multiplicity.toPointwise_withLinearXYs(lowerEps=1e-8, upperEps=1e-8) for nubar in
                              rfis.outputChannel.getProductsWithName('n')])
            fissionReport['ETA'] = computeETA(rcap.crossSection, rfis.crossSection, totalNubar,
                                                   useCovariance=useCovariance)
        return {title:fissionReport}
    else: return {title:Report()}


def getResonanceReport( rs, title="Resonances"):
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
    rrReport=Report(title=title)
    rrReport["RRR information"] = getRRRReport( rs, None )
    rrReport["URR information"] = getURRReport( rs, None )
    rrReport["Particle pairs"] = getParticlePairsReport( rs, None )
    return rrReport


def getRRRReport( rs, title="RRR information" ):
    '''

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
    '''
    if hasattr(rs.resonances, 'resolved') and rs.resonances.resolved is not None:
        import fudge.processing.resonances.reconstructResonances as RRReconstruct

        rrrReport = Report(title=title)
        rrrReport['Format'] = rs.resonances.resolved.evaluated.moniker
        if rrrReport['Format'] == "R_Matrix_Limited":
            rrrReport['Format'] += ' (%s approximation)'%rs.resonances.resolved.evaluated.approximation
        rrrReport["Upper bound"] = rs.resonances.resolved.domainMax
        rrrReport['Lower bound'] = rs.resonances.resolved.domainMin
        try:
            resCls = RRReconstruct.getResonanceReconstructionClass(rs.resonances.resolved.evaluated)
            rrr = resCls(rs, None, enableAngDists=False, verbose=False)
            rrr.setResonanceParametersByChannel()
            rrrReport["Scattering length (R')"]=PQU.PQU(rrr.getScatteringLength()*10, "fm")
            rrrReport['Lmax']=rrr.getLMax()
            rrrReport['No. channels']=rrr.nChannels
            rrrReport['No. resonances']=rrr.nResonances
        except ValueError: pass
        return rrrReport
    else:
        return Report()


def getURRReport( rs, title="URR information" ):
    '''
    URR (dict)
    ----------

        - Format (str)
            -- approximation (str)
        - Lower bound (PQU)
        - Upper bound (PQU)
        - No. channels (int)

    :return:
    '''
    if hasattr(rs.resonances, 'unresolved') and rs.resonances.unresolved is not None:
        urrReport=Report(title=title)
        urrReport['Format'] = rs.resonances.unresolved.moniker
        urrReport["Upper bound"] = rs.resonances.unresolved.domainMax
        urrReport['Lower bound'] = rs.resonances.unresolved.domainMin
        return urrReport
    else:
        return Report()


def getParticlePairsReport( rs, title="Particle Pairs" ):
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
    rep=Report(title=title)
    pps={}
    # Get the particle pairs
    if hasattr(rs.resonances, 'resolved') and rs.resonances.resolved is not None:
        import fudge.processing.resonances.reconstructResonances as RRReconstruct
        resCls = RRReconstruct.getResonanceReconstructionClass(rs.resonances.resolved.evaluated)
        rrr = resCls(rs, None, enableAngDists=False, verbose=False)
        try:
            rrr.setResonanceParametersByChannel(warnOnly=True)
        except ValueError as err:
            print err
            return rep
        for c in rrr.allChannels:
            if c.reaction not in pps: pps[c.reaction]=(c.particleA, c.particleB)
    for pp in pps:
        if pp in ['fission','competitive']: continue
        rep[pp]={}
        spins = rrr.getParticleSpins(pp)
        for ip,p in enumerate(pps[pp]):
            if p == 'gamma': particle=rs.getParticle('photon')
            else: particle=rs.getParticle(p)
            rep[pp][p]={}

            if isinstance(particle,PoPs.families.baryon.particle):
                rep[pp][p]["Mass (amu)"] = particle.mass.float('amu')
                rep[pp][p]["Spin (hbar)"] = particle.spin.float('hbar')
                rep[pp][p]["Parity"] = '+' if str(particle.parity[0])=="1 " else '-'
                rep[pp][p]["Charge (e)"] = particle.charge.float('e')
            elif isinstance(particle,PoPs.families.gaugeBoson.particle):
                rep[pp][p]["Mass (amu)"] = particle.mass.float('amu')
                rep[pp][p]["Spin (hbar)"] = particle.spin.float('hbar')
                rep[pp][p]["Parity"] = '+' if str(particle.parity[0])=="1 " else '-'
                rep[pp][p]["Charge (e)"] = particle.charge.float('e')
            elif isinstance(particle,PoPs.groups.isotope.suite):
                try:    rep[pp][p]["Mass (amu)"] = particle.mass.float('amu')
                except: rep[pp][p]["Mass (amu)"] = 'Unknown'
                try:    rep[pp][p]["Spin (hbar)"] = particle.nucleus.spin.float('hbar')
                except: rep[pp][p]["Spin (hbar)"] = 'Unknown'
                try:    rep[pp][p]["Parity"] = str(particle.nucleus.parity[0])
                except: rep[pp][p]["Parity"] = 'Unknown'
                try:    rep[pp][p]["Charge (e)"] = particle.nucleus.charge.float('e')
                except: rep[pp][p]["Charge (e)"] = 'Unknown'
            else: raise TypeError("Unknown particle type %s"%str(particle.__class__))
    return rep


# _______ plots ________

def getChannelPlots( rs, metricMenu, title="Channels", channelMatchFunction=lambda c: True, verbose=False ):
    """
        - Excitation functions
            -- Staircase: Cum. no. resonances vs. E (plot)
            -- <D> vs. E (plot)
            -- <Gamma> vs. E (plot)
            -- Transmission coefficients: Tc vs. E (plot)
            -- Scattering radii: Rc vs. E (plot)
            -- Pole strength: Sc vs. E (plot)
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
        urrAveQuant = urr.getWidthsAndSpacings(egrid, interpolateWidths=flag)
    else:
        urr = None
        urrAveQuant = None

    # get RRR parameter averages
    if hasattr(rs.resonances, 'resolved') and rs.resonances.resolved is not None:
        resCls = RRReconstruct.getResonanceReconstructionClass(rs.resonances.resolved.evaluated)
        rrr = resCls(rs, None, enableAngDists=False, verbose=verbose)
        rrr.setResonanceParametersByChannel()
        rrrAveQuant = rrr.getAverageQuantities(computeUncertainty=True, nBins=3) #resonancesPerBin=50)
    else:
        rrr = None
        rrrAveQuant = None

    def get_key_from_channel(c):
        if int(c.s) == c.s:
            s = c.s / 2.0
        else:
            s = c.s
        return '%s (j=%s,l=%s,s=%s)' % (str(c.reaction), str(c.J), str(c.l), str(s))

    # Excitation functions...
    if (hasattr(metricMenu, 'aveSpacingPlot') and metricMenu.aveSpacingPlot):
        if verbose: print winged_banner('aveSpacingPlot')
        if rrr is not None:
            rrrAveWidth=rrrAveQuant['spacings']
            rrrEgrid=rrrAveQuant['energyGrid']
        if urr is not None:
            raise NotImplementedError("urr data for aveSpacingPlot: FIXME")
        raise NotImplementedError("aveSpacingPlot: FIXME")

    if (hasattr(metricMenu, 'aveWidthPlot') and metricMenu.aveWidthPlot):
        if verbose: print winged_banner('aveWidthPlot')
        if rrr is not None:
            rrrAveWidth=rrrAveQuant['widths']
            rrrEgrid=rrrAveQuant['energyGrid']
        if urr is not None:
            raise NotImplementedError("urr data for aveWidthPlot: FIXME")
        raise NotImplementedError("aveWidthPlot: FIXME")

    if (hasattr(metricMenu, 'scatteringRadiusPlot') and metricMenu.scatteringRadiusPlot):
        if verbose: print winged_banner('scatteringRadiusPlot')
        raise NotImplementedError("scatteringRadiusPlot: FIXME")

    if (hasattr(metricMenu, 'poleStrengthPlot') and metricMenu.poleStrengthPlot):
        if rrr is not None:
            Sc = rrr.getPoleStrength(computeUncertainty=True)
            for c in Sc:
                if not channelMatchFunction(c): continue
                if verbose: print winged_banner("poleStrengthPlot for %s" % get_key_from_channel(c))
                #print Sc[c]
                #print Sc[c]['bin_edges'], Sc[c]['hist']
                Schist=grouped_values_to_XYs(Sc[c]['bin_edges'], Sc[c]['hist'], \
                          domainUnit='eV', domainName='energy_in', rangeUnit='', rangeName='S_c)', \
                          accuracy=upperEps, domainMin=0.0, domainMax=8.0)
                exit()
                p = Plot(title="Pole strength for %s channel" % get_key_from_channel(c), datasets=[Schist])
                p.doit()
        raise NotImplementedError("poleStrengthPlot: FIXME")

    if (hasattr(metricMenu, 'transmissionCoeff') and metricMenu.transmissionCoeff):
        if verbose: print winged_banner("transmissionCoeff")
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
            if rrrTc is not None: ds.append(rrrTc[c])
            if urrTc is not None: ds.append(urrTc[c])
            if urrTc is None and rrrTc is None:
                p = None
            else:
                p = Plot(title="T_c for %s channel" % get_key_from_channel(c), datasets=ds)
                p.doit()
        print NotImplementedError("transmissionCoeff: FIXME")

    if (hasattr(metricMenu, 'staircasePlot') and metricMenu.staircasePlot):
        if verbose: print winged_banner("staircasePlot")
        if rrr is not None:
            for c in rrr.channels:
                if not channelMatchFunction(c): continue
                ERList=[ rrr._energies[iR] for iR in rrr.channels[c] ]
                if len(ERList)<2:
                    print "Channel %c doesn't have enough resonances to generate a meaningful staircase plot"%get_key_from_channel(c)
                    continue
                ERList.sort()
                theData = [[ERList[0],1]]
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
    if (hasattr(metricMenu, 'PorterThomas') and metricMenu.PorterThomas):
        if rrr is not None:
            rrrResults = rrr.getPorterThomasFitToWidths()
            for c in rrrResults:
                if not channelMatchFunction(c): continue
                if verbose: print winged_banner("PorterThomas for %s" % get_key_from_channel(c))
                print rrrResults[c]['bin_edges'], rrrResults[c]['hist']
                PD=grouped_values_to_XYs(rrrResults[c]['bin_edges'], rrrResults[c]['hist'],
                          domainUnit='', domainName='D/<D>', rangeUnit='', rangeName='P(D/<D>)',
                          accuracy=upperEps, domainMin=0.0, domainMax=8.0)
                p = Plot(title="Porter-Thomas analysis for %s channel" % get_key_from_channel(c), datasets=[PD])
                p.doit()
        print NotImplementedError("PorterThomas: FIXME")

    if (hasattr(metricMenu, 'spacingDistributionPlot') and metricMenu.spacingDistributionPlot):
        if verbose: print winged_banner("spacingDistributionPlot")
        print NotImplementedError("spacingDistributionPlot: FIXME")

    if (hasattr(metricMenu, 'DysonMehtaDelta3') and metricMenu.DysonMehtaDelta3):
        if verbose: print winged_banner("DysonMehtaDelta3")
        print NotImplementedError("DysonMehtaDelta3: FIXME")

    # Overriding defaults...
    if (hasattr(metricMenu, 'penetrabilityPlot') and metricMenu.penetrabilityPlot):
        if rrr is not None:
            for c in rrr.channels:
                if not channelMatchFunction(c): continue
                if verbose: print winged_banner("penetrabilityPlot for %s" % get_key_from_channel(c))
                Pc = function_to_XYs(func=lambda E, x: rrr.penetrationFactorByChannel(c, E), fpars={}, rangeUnit='',
                                     rangeName='Channel penetrability, P_c(E)')
                p = Plot(title="Penetrability %s channel" % get_key_from_channel(c), datasets=[Pc])
                p.doit()

    if (hasattr(metricMenu, 'phasePlot') and metricMenu.phasePlot):
        if rrr is not None:
            for c in rrr.channels:
                if not channelMatchFunction(c): continue
                if verbose: print winged_banner("phasePlot for %s" % get_key_from_channel(c))
                phic = function_to_XYs(func=lambda E, x: rrr.phiByChannel(c, E), fpars={}, rangeUnit='',
                                       rangeName='Channel phase, phi_c(E)')
                p = Plot(title="Phase for %s channel" % get_key_from_channel(c), datasets=[phic])
                p.doit()

    if (hasattr(metricMenu, 'shiftPlot') and metricMenu.shiftPlot):
        if rrr is not None:
            for c in rrr.channels:
                if not channelMatchFunction(c): continue
                if verbose: print winged_banner("shiftPlot for %s" % get_key_from_channel(c))
                phic = function_to_XYs(func=lambda E, x: rrr.shiftFactorByChannel(c, E), fpars={}, rangeUnit='',
                                       rangeName='Channel shift factor, phi_c(E)')
                p = Plot(title="Shift for %s channel" % get_key_from_channel(c), datasets=[phic])
                p.doit()

    if (hasattr(metricMenu, 'backgroundXSPlot') and metricMenu.backgroundXSPlot):
        if verbose: print winged_banner("backgroundXSPlot")
        print NotImplementedError("backgroundXSPlot: FIXME")


# _______ datatables ________

def getChannelDataTable( rs, metricMenu, title="Channels", doURR=False, verbose=False ):
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
        - Strength functions (PQU)

    :param rs:
    :param title:
    :return:
    """
    # Final results in a table for later report building
    columnHeaders = collections.OrderedDict() # columnHeader instances
    data = collections.OrderedDict() # nested list with the actual data!
    icol=0

    def addColumn(columns,icol,key,unit=None):
        columns[key]=xData.table.columnHeader(icol,key,unit=unit)
        icol+=1

    def updateColumnsAndRow(columns,row,icol,key,val,unit=None):
        if key not in columns: addColumn(columns,icol,key,unit=unit)
        row.append(val)


    # ---------------------------------
    # Process the resonances
    # ---------------------------------
    import fudge.processing.resonances.reconstructResonances as RRReconstruct

    # get URR parameters
    if hasattr(rs.resonances, 'unresolved') and rs.resonances.unresolved is not None:
        urr = RRReconstruct.URRcrossSection(rs)
        egrid, flag = urr.generateEnergyGrid()
        urr.getWidthsAndSpacings(egrid, interpolateWidths=flag)

    # get RRR parameter averages
    if hasattr(rs.resonances, 'resolved') and rs.resonances.resolved is not None:
        resCls = RRReconstruct.getResonanceReconstructionClass(rs.resonances.resolved.evaluated)
        rrr = resCls(rs, None, enableAngDists=False, verbose=verbose)
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
            if not c.reaction in reactions: reactions.append(c.reaction)
            key = get_key_from_channel(c)
            if not key in data: data[key] = []
    if doURR and urr is not None:
        for reaction in reactions:
            for lj in urr.averageWidthFuncs:
                key = '%s (j=%s,l=%s,s=%s)' % (str(reaction), str(float(lj[1])), str(lj[0]), '0.5')
                if not key in data: data[key] = []

    # Compute the simple stuff
    if rrr is not None:
        for c in rrr.allChannels:
            key=get_key_from_channel(c)
            numResonances=len([ x for x in rrr.allChannels[c] if rrr.allChannels[c][x]!=0.0])
            numNegEnRes=sum([rrr._energies[iER]<0 for iER in rrr.allChannels[c] if rrr.allChannels[c][iER]!=0.0])
            updateColumnsAndRow(columnHeaders, data[key], icol, 'No. resonances', numResonances, '')
            updateColumnsAndRow(columnHeaders, data[key], icol, 'No. resonances w/ ER<0', numNegEnRes, '')
            if False: updateColumnsAndRow(columnHeaders, data[key], icol, 'MT', -13, '') #FIXME
            updateColumnsAndRow(columnHeaders, data[key], icol, 'gfact', c.gfact, '')
            if False: updateColumnsAndRow(columnHeaders, data[key], icol, 'Q', -13, 'eV') #FIXME
            updateColumnsAndRow(columnHeaders, data[key], icol, 'Threshold E', c.Xi, 'eV')
            updateColumnsAndRow(columnHeaders, data[key], icol, 'Eliminated?', str(c.eliminated), '')
            updateColumnsAndRow(columnHeaders, data[key], icol, 'Competative?', str(c.channelClass==RRReconstruct.COMPETATIVECHANNEL), '')
            updateColumnsAndRow(columnHeaders, data[key], icol, 'Relativistic?', str(c.useRelativistic), '' )
            if False: updateColumnsAndRow(columnHeaders, data[key], icol, 'Has background xs?', str(False), '') #FIXME
            if False: updateColumnsAndRow(columnHeaders, data[key], icol, 'B_c', str(False), '') #FIXME
            if not c.eliminated:
                updateColumnsAndRow(columnHeaders, data[key], icol, 'Pot. scatt. only?', str(len([ x for x in rrr.allChannels[c] if rrr.allChannels[c][x]!=0.0])==0), '')
            else: updateColumnsAndRow(columnHeaders, data[key], icol, 'Pot. scatt. only?', str(False), '')
        for c in rrrAveQuant:
            key=get_key_from_channel(c)
 #           print rrrAveQuant[c]['spacings']
            updateColumnsAndRow(columnHeaders, data[key], icol, 'RRR <D>', str(rrrAveQuant[c]['spacings'][0]), 'eV')
            updateColumnsAndRow(columnHeaders, data[key], icol, 'RRR <Gamma>', str(rrrAveQuant[c]['widths'][0]), 'eV')
    if doURR and urr is not None:
        raise NotImplementedError("TODO: URR average extraction")
        for reaction in reactions:
            for lj in urr.averageWidthFuncs:
                key = '%s (j=%s,l=%s,s=%s)' % (str(reaction), str(float(lj[1])), str(lj[0]), '0.5')
                try:
                    updateColumnsAndRow(columnHeaders, data[key], icol, 'URR <D>',  urr.levelSpacings[lj][str(reaction)], 'eV')
                except TypeError: pass
                try:
                     updateColumnsAndRow(columnHeaders, data[key], icol, 'URR <Gamma>', urr.averageWidths[lj][str(reaction)], 'eV')
                except TypeError: pass

    # Tough to extract quantities
    if (hasattr(metricMenu, 'effectiveDOF') and metricMenu.effectiveDOF):
        if rrr is not None:
            try:
                rrrResults = rrr.getPorterThomasFitToWidths()
            except AttributeError:
                print("Are you using an old version of SciPy?  scipy.optimize.curve_fit is missing.")
                rrrResults=[]
            for c in rrrResults:
                nu = PQU.PQU(rrrResults[c]['dof'], uncertainty=rrrResults[c]['ddof'])
                key=get_key_from_channel(c)
                updateColumnsAndRow(columnHeaders, data[key], icol, 'RRR DOF', str(nu), '')
        if doURR and urr is not None:
            raise NotImplementedError("TODO: URR average extraction")
            for reaction in reactions:
                for lj in urr.DOFs:
                    key = '%s (j=%s,l=%s,s=%s)' % (str(reaction), str(float(lj[1])), str(lj[0]), '0.5')
                    try:
                        updateColumnsAndRow(columnHeaders, data[key], icol, 'URR DOF', urr.DOFs[lj][str(reaction)], '')
                    except KeyError: pass
    if (hasattr(metricMenu, 'strengthFunction') and metricMenu.strengthFunction):
        if rrr is not None:
            Sc = rrr.getStrengthFunction(computeUncertainty=True)
            for c in Sc:
                key=get_key_from_channel(c)
                updateColumnsAndRow(columnHeaders, data[key], icol, 'RRR Sc', str(Sc[c]), '')

    nCols=len(columnHeaders.values())
    tmpData=[]
    for irow,row in enumerate(data.values()):
        if row==[]: tmpData.append([0 for i in range(nCols)])
        elif len(row)< nCols:
            for icell in range(len(row),nCols): row.append(0)
            tmpData.append(row)
        else: tmpData.append(row)
    result = DataTable(data=tmpData, columns=columnHeaders.values(), rows=data.keys())
    result.title = title
    return result


def getReactionDataTable( rs, metricMenu, title="Reactions", MTList=[], useCovariance=True ):
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
    :param MTList:
    :return:
    """
    # Final results in a table for later report building
    columnHeaders = collections.OrderedDict() # columnHeader instances
    rowNames = [] # really just row names
    data = [] # nested list with the actual data!
    icol=0

    def addColumn(columns,icol,key,unit=None):
        columns[key]=xData.table.columnHeader(icol,key,unit=unit)
        icol+=1

    def updateColumnsAndRow(columns,row,icol,key,val,unit=None):
        if key not in columns: addColumn(columns,icol,key,unit=unit)
        row.append(val)

    import fudge.gnds.reactions.reaction, fudge.gnds.sums

    # ---------------------------------
    #  Process each requested reaction
    # ---------------------------------
    for r in rs:

        # Get rid of all reactions that are not plain reactions or sums of reactions (e.g. no fission components of production stuff)
        if not isinstance(r, (fudge.gnds.reactions.reaction.reaction, fudge.gnds.sums.crossSectionSum)): continue

        # Keep only the MT's we actually want to test
        mt = r.ENDF_MT
        if mt not in MTList: continue

        # Initialize this row
        rowNames.append(str(r))
        row=[]

        # Check if we need to add the default columns
        updateColumnsAndRow(columnHeaders, row, icol, 'MT', mt)

        # Set Q
        if hasattr(r, 'getQ'): Q = r.getQ(unit='eV')
        else:                  Q = 0.0
        updateColumnsAndRow(columnHeaders, row, icol, 'Q', Q, 'eV')

        # Set EThreshold
        if hasattr(r, "getThreshold"): EThreshold = r.getThreshold(unit='eV')
        else:                          EThreshold = 0.0
        updateColumnsAndRow(columnHeaders, row, icol, 'Threshold', EThreshold, 'eV')

        # Resonance Integral
        if hasattr(metricMenu, 'RI') and metricMenu.RI:
            if EThreshold<0.0235339347844:
                RI = computeRI(r.crossSection, useCovariance=useCovariance)
                updateColumnsAndRow(columnHeaders, row, icol, 'RI', RI, RI.unit)
            else: row.append(blank())

        # 14 MeV Point
        if hasattr(metricMenu, 'fourteen') and metricMenu.fourteen:
            if LEGACYINTERMODE:
                # Legacy INTER uses 14.0 for the peak of the d-t outgoing neutron spectrum
                E14 = PQU.PQU(14., 'MeV')
            else:
                E14 = PQU.PQU(14.2,'MeV')
            fourteen = computeFourteenMeVPoint(r.crossSection, E14, useCovariance=useCovariance)
            updateColumnsAndRow(columnHeaders, row, icol, str(E14), fourteen, fourteen.unit)

        # Window average
        if hasattr(metricMenu, 'windowCenter') and hasattr(metricMenu, 'windowWidth') and metricMenu.windowCenter and metricMenu.windowWidth:
            window=computeAverageOverWindow( r.crossSection, PQU.PQU(metricMenu.windowCenter,unit='keV'), PQU.PQU(metricMenu.windowWidth,unit='keV'), useCovariance=useCovariance )
            updateColumnsAndRow(columnHeaders, row, icol, 'Ave. over %s+/-%s keV'%(str(metricMenu.windowCenter), str(metricMenu.windowWidth/2.0)), window, window.unit)

        # User defined spectrum average
        if hasattr(metricMenu, 'userDefinedSpecAve') and metricMenu.userDefinedSpecAve:
            path='/'.join(metricMenu.userDefinedSpecAve.split('/')[:-1])
            key=metricMenu.userDefinedSpecAve.split('/')[-1]
            ave=computeUserDefinedSpectrumAve( r.crossSection, path, key, useCovariance=useCovariance )
            updateColumnsAndRow(columnHeaders, row, icol, key, ave, ave.unit)

        # Lookup spectrum average
        if hasattr(metricMenu, 'otherSpecAve') and metricMenu.otherSpecAve and metricMenu.otherSpecAve not in ['Godiva','Jezebel','BigTen']:
            ave=computeLookupSpectrumAve( r.crossSection, metricMenu.otherSpecAve, useCovariance=useCovariance )
            updateColumnsAndRow(columnHeaders, row, icol, metricMenu.otherSpecAve, ave, ave.unit)

        # Room temperature cross section
        if hasattr(metricMenu, 'thermal') and metricMenu.thermal:
            if EThreshold<0.0235339347844:
                xsRT = computeRoomTempCS(r.crossSection, useCovariance=useCovariance)
                updateColumnsAndRow(columnHeaders, row, icol, "2200 m/s", xsRT, xsRT.unit)
            else: row.append(blank())

        # Westcott factor
        if hasattr(metricMenu, 'Westcott') and metricMenu.Westcott:
            if EThreshold<0.0235339347844 and computeRoomTempCS(r.crossSection, useCovariance=False).value != 0.0:
                if LEGACYINTERMODE:
                    # Legacy INTER ignores the recoil of the target,
                    # so a = 1.0 rather than a = mTarget/(mTarget+mProjectile)
                    westcott = computeWestcottFactor(r.crossSection, a=PQU.PQU(1., ""), useCovariance=useCovariance)
                else:
                    westcott = computeWestcottFactor(r.crossSection, useCovariance=useCovariance)
                updateColumnsAndRow(columnHeaders, row, icol, "Westcott factor", westcott, westcott.unit)
            else: row.append(blank())

        # Godiva spectrum average
        if hasattr(metricMenu, 'Godiva') and metricMenu.Godiva:
            if EThreshold<5e6 :
                x = computeGodivaSpectrumAve(r.crossSection, useCovariance=useCovariance)
                updateColumnsAndRow(columnHeaders, row, icol, "Godiva", x, x.unit)
            else: row.append(blank())

        # Jezebel spectrum average
        if hasattr(metricMenu, 'Jezebel') and metricMenu.Jezebel:
            if EThreshold<5e6:
                x = computeJezebelSpectrumAve(r.crossSection, useCovariance=useCovariance)
                updateColumnsAndRow(columnHeaders, row, icol, "Jezebel", x, x.unit)
            else: row.append(blank())

        # BigTen spectrum average
        if hasattr(metricMenu, 'BigTen') and metricMenu.BigTen:
            if EThreshold<5e6:
                x = computeBigTenSpectrumAve(r.crossSection, useCovariance=useCovariance)
                updateColumnsAndRow(columnHeaders, row, icol, "BigTen", x, x.unit)
            else: row.append(blank())

        # FUND-IPPE spectrum average
        if hasattr(metricMenu, 'FUNDIPPE') and metricMenu.FUNDIPPE:
            raise NotImplementedError("FIXME")
            if EThreshold<5e6:
                x = computeFUNDIPPESpectrumAve(r.crossSection, useCovariance=useCovariance)
                updateColumnsAndRow(columnHeaders, row, icol, "FUND-IPPE", x, x.unit)
            else: row.append(blank())

        # 252Cf spectrum average (using analytic approximation from INTER)
        if hasattr(metricMenu, 'CfSpectAnalytic') and metricMenu.CfSpectAnalytic:
            if EThreshold<5e6:
                x = computeCf252SpectrumAveAnalytic(r.crossSection, useCovariance=useCovariance)
                updateColumnsAndRow(columnHeaders, row, icol, "252Cf (analytic)", x, x.unit)
            else: row.append(blank())

        # 252Cf spectrum average
        if hasattr(metricMenu, 'CfSpectENDF') and metricMenu.CfSpectENDF:
            if EThreshold<5e6:
                x = computeCf252SpectrumAve(r.crossSection, useCovariance=useCovariance)
                updateColumnsAndRow(columnHeaders, row, icol, "252Cf", x, x.unit)
            else: row.append(blank())

        # MACS
        if hasattr(metricMenu, 'MACS') and metricMenu.MACS:
            if metricMenu.MACS:
                if isinstance(metricMenu.MACS,PQU.PQU): kT=metricMenu.MACS
                else: kT = PQU.PQU(metricMenu.MACS, 'keV')
            else:
                kT = PQU.PQU(30., 'keV')
            if EThreshold/10<kT.getValueAs('eV'):
                x = computeMACS(r.crossSection, T=kT, useCovariance=useCovariance)
                updateColumnsAndRow(columnHeaders, row, icol, "MACS(%s)" % str(kT), x.getValueAs('mb'), 'mb') # preferred unit for MACS is mb
            else: row.append(blank())

        # ARR
        if hasattr(metricMenu, 'ARR') and metricMenu.ARR:
            if metricMenu.ARR:
                if isinstance(metricMenu.ARR,PQU.PQU): kT=metricMenu.ARR
                else: kT = PQU.PQU(metricMenu.ARR, 'keV')
            else:
                kT = PQU.PQU(30., 'keV')
            if EThreshold/10<kT.getValueAs('eV'):
                x = computeAstrophysicalReactionRate(r.crossSection, T=kT, useCovariance=useCovariance)
                updateColumnsAndRow(columnHeaders, row, icol, "ARR(%s)" % str(kT), x, x.unit)
            else: row.append(blank())

        data.append(row)

    result = DataTable(data=data, columns=columnHeaders.values(), rows=rowNames)
    result.title=title
    return result


def getReactionEngineeringDataTable( rs, title="Engineering metrics", MTList=[], useCovariance=True ):
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
    metricMenu = collections.namedtuple('metricMenu',
                                        "RI Westcott thermal fourteen")  # emulates results from ArgumentParser.parse_args() function
    results=Report()
    results[title] = getReactionDataTable( rs, metricMenu=metricMenu(RI=True, Westcott=True, thermal=True, fourteen=True), title=title, MTList=MTList, useCovariance=useCovariance )
    return results


def getReactionAstrophysicsDataTable( rs, title="Astrophysics metrics", MTList=[], useCovariance=True ):
    """
    Astrophysics metrics ::
        - MACS(30 keV) (PQU)
        - ARR (PQU)

    :param rs:
    :param title:
    :param MTList:
    :return:
    """
    metricMenu = collections.namedtuple('metricMenu',
                                        "MACS ARR")  # emulates results from ArgumentParser.parse_args() function
    results=Report()
    results[title] = getReactionDataTable(rs,
                                          metricMenu=metricMenu(MACS=PQU.PQU(30,'keV'),ARR=PQU.PQU(30,'keV')),
                                          title=title, MTList=MTList, useCovariance=useCovariance)
    return results


def getReactionLegacyINTERDataTable( rs, title="Legacy INTER metrics", MTList=[1,2,18,102,16], useCovariance=True ):
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
    results=Report()
    results[title] = getReactionDataTable(rs, metricMenu=metricMenu(RI=True, Westcott=True, thermal=True, CfSpectAnalytic=True, CfSpectENDF=True, MACS=PQU.PQU(30.,'keV'), fourteen=True), title=title, MTList=MTList, useCovariance=useCovariance)
    return results


def getReactionAIntegralDataTable( rs, title="Integral metrics", MTList=[], useCovariance=True ):
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
    metricMenu = collections.namedtuple('metricMenu',
                                        "Godiva Flattop Jezebel BigTen FUNDIPPE CfSpectAnalytic CfSpectENDF")  # emulates results from ArgumentParser.parse_args() function
    results=Report()
    results[title] = getReactionDataTable(rs, metricMenu=metricMenu( Godiva=True, Flattop=False, Jezebel=True, BigTen=True, FUNDIPPE=False, CfSpectAnalytic=True, CfSpectENDF=True ), title=title, MTList=MTList, useCovariance=useCovariance)
    return results



# -------------------------------------------------------------------------------
# obsolete code
# -------------------------------------------------------------------------------

class OtherReport(Report):

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
        if True:
            report = [line.replace('Gamma',SYMBOL('Gamma')) for line in report]
        return report



# -------------------------------------------------------------------------------
# Unit tests
# -------------------------------------------------------------------------------

if __name__=="__main__":
    import unittest, sys
    from fudge.gnds import reactionSuite
    from xData.isclose import isclose

    TEST_DATA_PATH, this_filename = os.path.split(os.path.realpath(__file__))
    TEST_DATA_PATH += '/testData/'
    VERBOSE = '-v' in sys.argv or '--verbose' in sys.argv
    DOPLOTS = '-p' in sys.argv or '--plots' in sys.argv
    if DOPLOTS:
        if '-p' in sys.argv: sys.argv.remove('-p')
        if '--plots' in sys.argv: sys.argv.remove('--plots')

    if VERBOSE: print 'reading test data...'
    H1TestData = reactionSuite.readXML(open(TEST_DATA_PATH + os.sep + 'n-001_H_001.endf.gnds.xml'))
    Fe56TestData = reactionSuite.readXML(open(TEST_DATA_PATH + os.sep + 'n-026_Fe_056.endf.gnds.xml'))
    U241TestData = reactionSuite.readXML(open(TEST_DATA_PATH + os.sep + 'n-092_U_241.endf.gnds.xml'))


    class TestCaseWithExtras(unittest.TestCase):

        def assertWithinXPercent(self, a, b, percent=1.0, absTol=1e-14):
            return isclose(a, b, rel_tol=percent / 100., abs_tol=absTol, method='weak')

        def assertAllWithinXPercent(self, aL, bL, percent=1.0, absTol=1e-14):
            if len(aL) != len(bL): return False
            return all([isclose(*x, rel_tol=percent / 100., abs_tol=absTol, method='weak') for x in zip(aL, bL)])

        def assertJSONEqual(self, a, b):
            """Compares two strings that are JSON dumps of dictionaries by converting them back to dicts"""
            self.assertDictEqual(json.loads(a),json.loads(b))

        def assertStringsEqual(self,a,b):
            aList = [x.rstrip() for x in a.split('\n')]
            bList = [x.rstrip() for x in b.split('\n')]
            return self.assertItemsEqual(aList,bList)

    class TestReport(TestCaseWithExtras):

        def setUp(self):
            self.example = Report(title='test data')
            self.example['a'] = 4
            self.example['b'] = True
            self.example['c'] = 3.2
            self.example['d'] = 'fred'
            self.example['e'] = None
            self.example['f'] = [1, 2]
            self.example['g'] = (3, 4)
            self.example['h'] = {5,6}
            self.example['i'] = collections.OrderedDict( ((2,'k'),('k','l')) )
            self.example['j'] = { 43:'k', 'k':'l' }
            self.example['k'] = Report(title='empty')
            self.example['l'] = PQU.PQU(1.2,'keV')

        def test_text_report(self):
            self.assertEqual(self.example.text_title(), '+-----------+\n| test data |\n+-----------+')
            self.assertJSONEqual(self.example.text_report(), '{\n    "test data": {\n        "a": 4,\n        "b": true,\n        "c": 3.2,\n        "d": "fred",\n        "e": null,\n        "f": [\n            1,\n            2\n        ],\n        "g": [\n            3,\n            4\n        ],\n        "h": [\n            5,\n            6\n        ],\n        "i": {\n            "2": "k",\n            "k": "l"\n        },\n        "j": {\n            "k": "l",\n            "43": "k"\n        },\n        "k": {},\n        "l": "1.2 keV"\n    }\n}' )

        def test_json_report(self):
            self.assertEqual(self.example.json_report(), '{"test data": {"a": 4, "b": true, "c": 3.2, "d": "fred", "e": null, "f": [1, 2], "g": [3, 4], "h": [5, 6], "i": {"2": "k", "k": "l"}, "j": {"k": "l", "43": "k"}, "k": {}, "l": "1.2 keV"}}')

        def test_properties(self):
            self.assertEqual( self.example.id, id(self.example) )
            self.assertEqual( self.example.cls, 'report testdata' )

        @unittest.skip("FIXME:didn't get working yet")
        def test_html_report(self):
            self.assertEqual(self.example.html_report(), '')

    class TestPlot(unittest.TestCase):

        def test_init(self):
            a=Plot()

        def test_properties(self):
            a = Plot(title='just a plot test')
            self.assertEqual(a.id, id(a))
            self.assertEqual(a.cls, 'plot justaplottest')

    class TestGetReports(TestCaseWithExtras):

        def test_getEvaluationReport(self):
            self.maxDiff = None
            self.assertJSONEqual( getEvaluationReport(Fe56TestData, title='n-026_Fe_056.endf.gnds.xml').text_report(), '''{
    "n-026_Fe_056.endf.gnds.xml": {
        "Authors": "IAEA Consortium",
        "Lab": "IAEA",
        "Evaluation date": "2016-03-22",
        "MAT": null,
        "Projectile": "n",
        "Target": "Fe56",
        "Compound nucleus": "Fe57",
        "Projectile frame": "lab",
        "Temperature": "0. K",
        "GNDS version": "GNDS 1.8"
    }
}''' )

        def test_getFissionMetricReport(self):
            metricMenu=collections.namedtuple('metricMenu',"ALF ETA") # emulates results from ArgumentParser.parse_args() function
            rpt = getFissionMetricReport(U241TestData,metricMenu(ALF=True,ETA=True))["Fission metrics"].text_report()
            a = json.loads(rpt)
            b = json.loads('''{
    "Fission metrics": {
        "ALF": "1143.149783603068",
        "ETA": "1.952768515392308e-3"
    }
}''' )
            for k in a['Fission metrics']:
                self.assertWithinXPercent(float(a['Fission metrics'][k]),float(b['Fission metrics'][k]),.01)

        def test_getResonanceReport(self):
            self.assertJSONEqual(getResonanceReport(H1TestData).text_report(), """{
    "Resonances": {
        "RRR information": {},
        "URR information": {},
        "Particle pairs": {}
    }
}""")
            self.assertJSONEqual(getResonanceReport(Fe56TestData).text_report(), '''{
    "Resonances": {
        "RRR information": {
            "Format": "RMatrix",
            "Upper bound": 850000.0,
            "Lower bound": 1e-05,
            "Scattering length (R')": "5.436999999987526 fm",
            "Lmax": 4,
            "No. channels": 5,
            "No. resonances": 313
        },
        "URR information": {},
        "Particle pairs": {
            "Fe57 + photon": {
                "Fe57": {
                    "Charge (e)": 26.0,
                    "Parity": "Unknown",
                    "Mass (amu)": "Unknown",
                    "Spin (hbar)": "Unknown"
                },
                "photon": {
                    "Charge (e)": 0.0,
                    "Parity": "+",
                    "Mass (amu)": 0.0,
                    "Spin (hbar)": 1.0
                }
            },
            "n + Fe56": {
                "Fe56": {
                    "Charge (e)": 26.0,
                    "Parity": "Unknown",
                    "Mass (amu)": 55.9349379634,
                    "Spin (hbar)": 0.0
                },
                "n": {
                    "Charge (e)": 0.0,
                    "Parity": "+",
                    "Mass (amu)": 1.00866491574,
                    "Spin (hbar)": 0.5
                }
            }
        }
    }
}
''')
            self.assertJSONEqual(getResonanceReport(U241TestData).text_report(), """{
    "Resonances": {
        "RRR information": {
            "Format": "SingleLevel_BreitWigner",
            "Upper bound": 102.5,
            "Lower bound": 1e-05
        },
        "URR information": {
            "Format": "unresolved",
            "Upper bound": 10000.0,
            "Lower bound": 102.5
        },
        "Particle pairs": {}
    }
}""")

        def test_getRRRReport(self):
            self.assertJSONEqual(getRRRReport(H1TestData).text_report(), '''{}''')
            self.assertJSONEqual(getRRRReport(U241TestData).text_report(), '''{
                "RRR information": {
                    "Format": "SingleLevel_BreitWigner",
                    "Upper bound": 102.5,
                    "Lower bound": 1.0e-5
                }
            }''')
            k="RRR information"
            a = json.loads(getRRRReport(Fe56TestData).text_report())
            b = json.loads( '''{
                "RRR information": {
                    "Format": "RMatrix",
                    "Upper bound": 8.5e5,
                    "Lower bound": 1.0e-5,
                    "Scattering length (R\')": "5.436999999987528 fm",
                    "Lmax": 4,
                    "No. channels": 5,
                    "No. resonances": 313
                }
            }''')
            for kk in a[k]:
                if kk =="Scattering length (R\')":
                    self.assertWithinXPercent(float(a[k][kk].replace('fm','')), float(b[k][kk].replace('fm','')), 0.01)
                elif 'bound' in kk:
                    self.assertWithinXPercent(a[k][kk],b[k][kk],0.01)
                else:
                    self.assertEqual(a[k][kk],b[k][kk])

        def test_getURRReport(self):
            self.assertJSONEqual( getURRReport(U241TestData).text_report(), '''{
    "URR information": {
        "Format": "unresolved",
        "Upper bound": 10000.0,
        "Lower bound": 102.5
    }
}''' )

        def test_getParticlePairsReport(self):
            self.assertJSONEqual( getParticlePairsReport(Fe56TestData).text_report(), """{
    "Particle Pairs": {
        "Fe57 + photon": {
            "Fe57": {
                "Charge (e)": 26.0,
                "Parity": "Unknown",
                "Mass (amu)": "Unknown",
                "Spin (hbar)": "Unknown"
            },
            "photon": {
                "Charge (e)": 0.0,
                "Parity": "+",
                "Mass (amu)": 0.0,
                "Spin (hbar)": 1.0
            }
        },
        "n + Fe56": {
            "Fe56": {
                "Charge (e)": 26.0,
                "Parity": "Unknown",
                "Mass (amu)": 55.9349379634,
                "Spin (hbar)": 0.0
            },
            "n": {
                "Charge (e)": 0.0,
                "Parity": "+",
                "Mass (amu)": 1.00866491574,
                "Spin (hbar)": 0.5
            }
        }
    }
}""" )

        #@unittest.skip("Need to discuss what channels should be there now that LRF=3 translated into GNDS analog of LRF=7")
        def test_getChannelDataTable(self):
            self.maxDiff=None
            metricMenu=collections.namedtuple('metricMenu','effectiveDOF strengthFunction poleStrengthPlot scatteringRadiusPlot PorterThomas staircasePlot aveSpacingPlot aveWidthPlot shiftPlot penetrabilityPlot phasePlot DysonMehtaDelta3 transmissionCoeff')
            x = getChannelDataTable(Fe56TestData, metricMenu(effectiveDOF=True,
                                                             strengthFunction=True,
                                                             poleStrengthPlot=False,
                                                             scatteringRadiusPlot=False,
                                                             PorterThomas=False,
                                                             staircasePlot=False,
                                                             aveSpacingPlot=False,
                                                             aveWidthPlot=False,
                                                             shiftPlot=False,
                                                             penetrabilityPlot=False,
                                                             phasePlot=False,
                                                             DysonMehtaDelta3=False,
                                                             transmissionCoeff=False))
#            print '\n'.join(x.toStringList())
#            return
            self.assertStringsEqual( '\n'.join(x.toStringList()), '''name                      No. resonances No. resonances w/ ER<0 gfact Threshold E Eliminated? Competative? Relativistic? Pot. scatt. only?              RRR <D>                        RRR <Gamma>        RRR DOF               RRR Sc
unit                                                                           eV                                                                            eV                                 eV
elastic (j=0.5,l=0,s=0.5)             40                      3   1.0         0.0       False        False         False             False   2.5e4 +/- 1.2e4 eV                          12330. eV  1.e-3 +/- 4.8      2.0e4 +/- 1.4e4
capture (j=0.5,l=0,s=0.5)             40                      3   1.0         0.0       False        False         False             False   2.5e4 +/- 1.2e4 eV                            0.96 eV  1.42 +/- 0.42    1.9e-5 +/- 1.3e-5
elastic (j=0.5,l=1,s=0.5)             63                      0   1.0         0.0       False        False         False             False   1.3e4 +/- 1.1e4 eV                            166. eV  1.15 +/- 0.17      0.017 +/- 0.014
capture (j=0.5,l=1,s=0.5)             63                      0   1.0         0.0       False        False         False             False   1.3e4 +/- 1.1e4 eV  0.4899999999999999 +/- 1.1e-16 eV  0.74 +/- 0.30   1.07e-5 +/- 8.4e-6
elastic (j=1.5,l=2,s=0.5)             75                      0   2.0         0.0       False        False         False             False   6.4e3 +/- 6.5e3 eV               137.9 +/- 2.8e-14 eV  0.72 +/- 0.12                   0.
capture (j=1.5,l=2,s=0.5)             75                      0   2.0         0.0       False        False         False             False   6.4e3 +/- 6.5e3 eV                            0.72 eV  0.55 +/- 0.30                   0.
elastic (j=1.5,l=1,s=0.5)             79                      0   2.0         0.0       False        False         False             False   9.9e3 +/- 6.9e3 eV                            646. eV  1.28 +/- 0.11        0.19 +/- 0.15
capture (j=1.5,l=1,s=0.5)             79                      0   2.0         0.0       False        False         False             False   9.9e3 +/- 6.9e3 eV  0.4899999999999999 +/- 1.1e-16 eV  1.29 +/- 0.38    4.0e-5 +/- 3.1e-5
elastic (j=2.5,l=2,s=0.5)             56                      0   3.0         0.0       False        False         False             False  1.21e4 +/- 9.0e3 eV   90.30000000000003 +/- 2.8e-14 eV  1.06 +/- 0.13    0.0752 +/- 5.9e-3
capture (j=2.5,l=2,s=0.5)             56                      0   3.0         0.0       False        False         False             False  1.21e4 +/- 9.0e3 eV  0.7199999999999998 +/- 2.2e-16 eV  0.98 +/- 0.36  1.145e-5 +/- 9.0e-7
elastic (j=2.5,l=3,s=0.5)              0                      0   3.0         0.0       False        False         False              True                    0                                  0             1.                   -1
elastic (j=3.5,l=3,s=0.5)              0                      0   4.0         0.0       False        False         False              True                    0                                  0             1.                   -1''')

        @unittest.skipIf(not DOPLOTS,"Not really a unit test, really more of an interactive test")
        def test_getChannelPlots(self):
            metricMenu = collections.namedtuple('metricMenu',
                                                'poleStrengthPlot scatteringRadiusPlot PorterThomas staircasePlot aveSpacingPlot aveWidthPlot backgroundXSPlot penetrabilityPlot shiftPlot phasePlot DysonMehtaDelta3 spacingDistributionPlot transmissionCoeff')
            def channelMatchFunction(c): return True #return not c.eliminated and int(2.*c.J)==1 and c.l==0 and c.channelClass == 1
            getChannelPlots(Fe56TestData, metricMenu( poleStrengthPlot=False,
                                                      scatteringRadiusPlot=False,
                                                      PorterThomas=False,
                                                      staircasePlot=False,
                                                      aveSpacingPlot=False,
                                                      aveWidthPlot=False,
                                                      backgroundXSPlot=True,
                                                      penetrabilityPlot=False,
                                                      shiftPlot=False,
                                                      phasePlot=False,
                                                      DysonMehtaDelta3=False,
                                                      spacingDistributionPlot=False,
                                                      transmissionCoeff=False ), channelMatchFunction=channelMatchFunction, verbose=True )

        def test_getReactionLegacyINTERReport(self):
            a=json.loads(getReactionLegacyINTERDataTable( H1TestData ).text_report())
            b="<!--       name |  MT |       Q | Threshold |                RI |     14.2 MeV |          2200 m/s |   Westcott factor |    252Cf (analytic) |               252Cf |     MACS(30. keV)  -->\n<!--            |     |      eV |        eV |                 b |            b |                 b |                   |                   b |                   b |                mb  -->\n         n + H1     2         0           0     263.83668978541     0.67806042            20.43633    1.12837868958532       3.9270475011023      3.66801927980566    14533.7414524736\n    H2 + photon   102   2224631           0   0.149135213422371   2.9719802e-5   0.332281145255188   0.707206052903381   3.94192612648554e-5   3.86144745081925e-5   0.152449777631191\n          total     1         0           0    263.983267844317     0.67809014    20.7687319650388    1.12196774546104      3.92708606339206      3.66805680325653    14533.8730068349"
            alines = a["Legacy INTER metrics"].split('\n')
            blines = b.split('\n')
            for i in range(2,len(alines)):
                aline=alines[i].split()[-10:]
                bline=blines[i].split()[-10:]
                for j in range(len(aline)):
                    self.assertWithinXPercent(aline[j],bline[j],0.01)

    unittest.main()
