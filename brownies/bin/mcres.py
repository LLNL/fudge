#!/usr/bin/env python
import argparse
import json
import os
import fractions
import math
from fudge.processing.resonances.reconstructResonances import getAllowedTotalSpins as coupleToJ
from brownies.BNL.plot_evaluation.plotio import readEvaluation
import brownies.BNL.restools.resonance_generator as rg
from fudge.core.utilities.brb import banner, winged_banner, uniquify
import fudge.processing.resonances.reconstructResonances as RRReconstruct
from xData import table as tableModule
from xData import XYs1d as XYs1dModule
from PoPs.chemicalElements import misc as chemicalElementMiscPoPsModule
from PoPs.quantities.spin import Suite as spinSuite

# --------------------------------------------------------------------
#   Utility functions
# --------------------------------------------------------------------
def is_half_int(x):
    if type(x) == int:
        return True
    if type(x) not in [str, float]:
        return False
    _f = fractions.Fraction(x)
    return _f.denominator == 2

# --------------------------------------------------------------------
#   Read & quality check the input configuration file
# --------------------------------------------------------------------
def read_config(configFileName):
    """"""
    # Read the file
    __config = json.loads(open(configFileName).read())

    if __config['target'] != 'ENDF':
        if not all([x in __config['target'] for x in ["name", "Z", "A", "spin"]]):
            raise ValueError('target must be the key "ENDF" or a dictionary containing ["name", "Z", "A", "spin"]')

    if not isinstance(__config.get("ENDFFile", ''), str):
        raise TypeError("ENDFFile must be a string, got ", type(__config["ENDFFile"]))

    if __config.get('spacingDist', None) not in ['wigner', 'picket fence', 'poisson', 'brody', 'goe']:
        raise ValueError('spacingDist must be in ["wigner", "picket fence", "poisson", "brody", "goe"]')

    if __config.get('approximation', None) not in ["SLBW", "MLBW", "RM", "R", "ENDF"]:
        raise ValueError('approximation must be one of "SLBW", "MLBW", "RM", or "R" or the key "ENDF"')

    if __config['spingroups'] != 'ENDF':
        if not isinstance(__config['spingroups'], list):
            raise TypeError("spingroups should be a list")
        for sg in __config['spingroups']:
            if not all([x in sg for x in ["L", "J", "S", "D", "Gn", "Gg", "nDOF"]]):
                raise ValueError('each spingroup must be the key "ENDF" or a dictionary containing '
                                 '["L", "J", "S", "D", "Gn", "Gg", "nDOF"]')

    if not (isinstance(__config['Emin'], float) or __config['Emin'] == 'ENDF'):
        raise ValueError("Emin must be a float or the keyword 'ENDF'")

    if not (isinstance(__config['Emax'], float) or __config['Emax'] == 'ENDF'):
        raise ValueError("Emax must be a float or the keyword 'ENDF'")

    if __config['projectile'] not in ['ENDF', 'n']:
        raise NotImplementedError("Write code to handle non-neutron projectiles.")

    return __config


# --------------------------------------------------------------------
#   Command line
# --------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(description="Monte Carlo Resonances (McRes)")
    parser.add_argument('config', type=str,
                        help="JSON formatted configuration file, this lets us know what kind of resonances you want, "
                             "not what you want to do with them")
    parser.add_argument('-o', dest='outFile', default=None, type=str, help='Output file')
    parser.add_argument('-f', default='endf', choices=['endf', 'gnds'], help="Output format")
    parser.add_argument('-v', dest='verbose', default=False, action='store_true', help='Run verbosely')
    parser.add_argument('-q', dest='verbose', action='store_false', help='Run quietly')
    parser.add_argument('--appendRRR', default=False, action='store_true',
                        help="Insert generated resonances into endfOut file's RRR, leaving existing RRR resonances")
    parser.add_argument('--removeURR', default=False, action='store_true', help="Remove URR region from file")
    parser.add_argument('--overwriteRRR', default=False, action='store_true',
                        help="Overwrite current resonances in endfOut file's RRR with generated resonances")
    parser.add_argument('--overwriteURR', default=False, action='store_true',
                        help="Combines --appendRRR and --removeURR, effectively replacing an URR with a "
                             "series of stochastic RRR resonances")
    parser.add_argument('-l', dest='list', default=False, action="store_true", help="List generated resonances")
    return parser.parse_args()


# --------------------------------------------------------------------
#   Helper classes
# --------------------------------------------------------------------
class Nucleus:
    """
    Lightweight class to hold nuclear properties
    """

    def __init__(self):
        self.spin = None
        self.name = None
        self.Z = None
        self.A = None

    def __str__(self):
        return "%s (Z=%s,A=%s,s=%s)" % (self.name, str(self.Z), str(self.A), str(self.spin))


class AngMom(fractions.Fraction):
    """
    Lightweight Angular Momentum class, basically a thin wrapper around fraction

    Remember angular momentum can only be integers or half-integers.
    """

    def __new__(cls, s, dnom=None):
        """
        S can be a float (but must be an integer or half-integer),
        an int, a string (that is an integer or half-integer)
        or a PoPs-style spin suite.
        """
        if isinstance(s, spinSuite):
            s = float(s.pqu().value)
        _result = super(AngMom, cls).__new__(cls, numerator=s, denominator=dnom)
        if _result.denominator not in [1, 2]:
            raise ValueError("Angular momentum can only be integers or half-integers")
        return _result


class SpinGroup:
    """
    Lightweight wrapper around quantum number part of a channel
    """

    def __init__(self, L, J, S):
        self.L = AngMom(L)
        self.S = AngMom(S)
        self.J = AngMom(J)

    def __str__(self):
        return "(L=%s, J=%s, S=%s)" % (str(self.L), str(self.J), str(self.S))

    @property
    def lj(self):
        return (self.L, fractions.Fraction(self.J))


# --------------------------------------------------------------------
#   Initializers for data needed for fake resonance generator
# --------------------------------------------------------------------
def get_target_and_projectile(__config, __endfRS):
    __targ = Nucleus()
    __proj = Nucleus()
    if __config['target'] == 'ENDF':
        __targ.name = str(__endfRS.target)
        ZA = chemicalElementMiscPoPsModule.ZA(__endfRS.PoPs[__endfRS.target])
        __targ.Z = int(ZA / 1000)
        __targ.A = ZA - __targ.Z * 1000
        __targ.spin = AngMom(__endfRS.getParticle(__endfRS.target).spin)
    else:
        __targ.name = __config['target']['name']
        __targ.spin = __config['target']['spin']
        __targ.Z = __config['target']['Z']
        __targ.A = __config['target']['A']
    if __config['projectile'] == 'ENDF':
        __proj.name = str(__endfRS.projectile)
        ZA = chemicalElementMiscPoPsModule.ZA(__endfRS.PoPs[__endfRS.projectile])
        __proj.Z = int(ZA / 1000)
        __proj.A = ZA - __proj.Z * 1000
        __proj.spin = AngMom(__endfRS.getParticle(__endfRS.projectile).spin)
    else:
        if __config['projectile'] == 'n':
            __proj.name = 'n'
            __proj.spin = 0.5
            __proj.Z = 0
            __proj.A = 1
        else:
            raise NotImplementedError("Implement other projectiles")
    return __targ, __proj


def get_spin_groups(__config, ___urr=None, targSpin=None, projSpin=None):
    """
    Gets a list of all allowed orbital angular momentum/spin combinations assuming
    we couple the target and projectile spin up to a total spin S and then
    couple S and orbital angular momentum up to J.
    """
    __spingroups = []
    Ss = coupleToJ(int(2 * targSpin), int(2 * projSpin), True)
    if __config['spingroups'] == 'ENDF':
        # get list from urr
        for lj in dict(urr.averageWidths).keys():
            for S in Ss:
                if lj[1].numerator in coupleToJ(int(2 * lj[0]), S, True):
                    __spingroups.append(SpinGroup(lj[0], float(lj[1]), S / 2.0))
    else:
        for sgDict in __config['spingroups']:
            # Get allowed open channels/spin groups
            L = sgDict['L']
            J = sgDict['J']
            S = sgDict['S']
            Js = coupleToJ(int(2 * L), int(2*S), True)
            if int(2*J) not in Js:
                print("ERROR: Skipping spingroup: Can't couple "+str(L)+" + "+str(S)+" up to "+str(J)+
                      ", allowed values are", str([j/2 for j in Js]))
            else:
                __spingroups.append(SpinGroup(L, J, S))
    return __spingroups


def find_matching_spingroup(__configSG, __spingroups):
    for sg in __spingroups:
        #print(sg.lj, (__configSG['L'], __configSG['J']))
        #exit()
        if sg.lj == (__configSG['L'], __configSG['J']):
            return sg


def get_level_densities(__spingroups, __config, __urr):
    if __config["spingroups"] == "ENDF":
        if __urr is None:
            raise ValueError("Need URR region in evaluation to use ENDF file to initialize spacings")
        return dict(__urr.levelDensities)
    return {find_matching_spingroup(cfgsg, __spingroups).lj:
                XYs1dModule.XYs1d(data=[[__config['Emin']*0.8, 1.0/cfgsg["D"]], [__config['Emax']*1.2, 1.0/cfgsg["D"]]],
                                axes=XYs1dModule.XYs1d.defaultAxes(
                                    labelsUnits={
                                        XYs1dModule.yAxisIndex: ('level_density', '1/eV'),
                                        XYs1dModule.xAxisIndex: ('excitation_energy', 'eV')}))
            for cfgsg in __config['spingroups']}


def get_widths(__spingroups, __config, __urr):
    if __config['spingroups'] == 'ENDF':
        aveWidths = dict(urr.averageWidths)
        for lj in aveWidths:
            # ENDF neutron widths aren't actually neutron widths.
            # Rather, they are what ENDF calls "reduced widths".
            # We need to convert them back into widths
            # FIXME: There may be a rogue factor of 2 if two channel spins are possible for a given J,L
            for ip, p in enumerate(aveWidths[lj]['elastic']):
                f = math.sqrt(p[0]) * urr.penetrationFactor(lj[0], urr.rho(p[0])) / urr.rho(p[0])  # See ENDF Eq. (D.99)
                aveWidths[lj]['elastic'][ip] = [p[0], p[1] * f]
    else:
        aveWidths={}
        for cfgsg in __config['spingroups']:
            sg = find_matching_spingroup(cfgsg, __spingroups)
            aveWidths[sg.lj] = {}
            # Gamma widths
            aveWidths[sg.lj]['capture'] = XYs1dModule.XYs1d(
                                            data=[[__config['Emin']*0.8, cfgsg['Gg']], [__config['Emax']*1.2, cfgsg['Gg']]],
                                            axes=XYs1dModule.XYs1d.defaultAxes(
                                                labelsUnits={
                                                    XYs1dModule.yAxisIndex: ('width', 'eV'),
                                                    XYs1dModule.xAxisIndex: ('excitation_energy', 'eV')}))
            # Neutron widths
            aveWidths[sg.lj]['elastic'] = XYs1dModule.XYs1d(
                                            data=[[__config['Emin']*0.8, cfgsg['Gn']], [__config['Emax']*1.2, cfgsg['Gn']]],
                                            axes=XYs1dModule.XYs1d.defaultAxes(
                                                labelsUnits={
                                                    XYs1dModule.yAxisIndex: ('width', 'eV'),
                                                    XYs1dModule.xAxisIndex: ('excitation_energy', 'eV')}))
            # No fission or competitive yet
    return aveWidths


def get_DOFs(__spingroups, __config, __urr):
    if __config["spingroups"] == "ENDF":
        if __urr is None:
            raise ValueError("Need URR region in evaluation to use ENDF file to initialize DOFs")
        return dict(__urr.DOFs)
    return {find_matching_spingroup(cfgsg, __spingroups).lj: {'elastic': cfgsg["nDOF"], 'capture': 0}
            for cfgsg in __config['spingroups']}


def get_energy_range(__config, __urr, __rrr, __appendRRR, __overwriteURR):
    # Set requested energy range (if 'ENDF' option used) or
    # check for compatibility of requested energy range and what the ENDF file provides
    if __config['Emin'] == 'ENDF':
        # The order of these if-then's hopefully ensures the lowest appropriate lower bound is used
        if __appendRRR or __overwriteURR:
            __config['Emin'] = __urr.lowerBound
        if args.overwriteRRR:
            __config['Emin'] = __rrr.lowerBound
    else:
        if __config['spingroups'] == 'ENDF' and __config['Emin'] < __urr.lowerBound:
            raise ValueError("Emin < lower bound for URR")
    if __config['Emax'] == 'ENDF':
        # The order of these if-then's hopefully ensures the highest appropriate upper bound is used
        if __overwriteURR:
            __config['Emax'] = __rrr.upperBound
        if __appendRRR or __overwriteURR:
            __config['Emax'] = __urr.upperBound
    else:
        if __config['spingroups'] == 'ENDF' and __config['Emax'] > __urr.upperBound:
            raise ValueError("Emax > upper bound for URR")
    return __config['Emin'], __config['Emax']


def get_last_resonance_energy(_rrr, _l, _j):
    for c in _rrr.channels:
        if c.l == _l and c.J == _j:
            return _rrr._energies[max(_rrr.channels[c].keys())]
    return None


def get_first_resonance_energy(_rrr, _l, _j):
    for c in _rrr.channels:
        if c.l == _l and c.J == _j:
            return _rrr._energies[min(_rrr.channels[c].keys())]
    return None


def add_resonances_to_suite(_res_table, _reaction_suite, L, J, _cut_index=None):
    moniker = _reaction_suite.resonances.resolved.evaluated.moniker
    if moniker == 'BreitWigner':
        if _cut_index is not None and len(
                _reaction_suite.resonances.resolved[0].resonanceParameters.table.data) > _cut_index:
            _reaction_suite.resonances.resolved[0].resonanceParameters.table.data = \
                _reaction_suite.resonances.resolved[0].resonanceParameters.table.data[:_cut_index]
        for res in _res_table:
            _reaction_suite.resonances.resolved[0].addResonance(
                [res[0], res[1], float(res[2]), res[3] + res[4], res[3], res[4]])
            _reaction_suite.resonances.resolved[0].resonanceParameters.table.data.sort(key=lambda x: x[0])
        return
    elif moniker == 'RMatrix':
        NotImplementedError("Reconstruction class %s is not supported" % moniker)
    else:
        pass
    raise ValueError("Reconstruction class %s is not supported" % moniker)


# --------------------------------------------------------------------
#   Main!!
# --------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()

    # Take care of things we can't do yet
    if args.overwriteRRR:
        raise NotImplementedError("Write code to overwrite resonances in an ENDF file")

    # Declare variables
    endfRS, endfCS = None, None
    rrr, urr = None, None

    # Read in the configuration and make sure we can use it
    if args.verbose:
        print(banner("Reading configuration from %s" % args.config))
    config = read_config(args.config)

    # Read in the ENDF file, if there is one, so we can get the URR parameters, if there are any
    if 'ENDFFile' in config:
        endfFileName = config['ENDFFile'].replace('~', os.environ['HOME'])
        print(banner("Reading ENDF file %s" % endfFileName))
        endfRS, endfCS = readEvaluation(endfFileName, reconstructResonances=False, verbose=args.verbose,
                                        skipBadData=True)

        # Get the RRR parameters
        if not hasattr(endfRS.resonances, 'resolved'):
            raise ValueError("Evaluation %s does not have a RRR" % endfFileName)
        if args.verbose:
            print(winged_banner("Extracting RRR parameters"))
        resCls = RRReconstruct.getResonanceReconstructionClass(endfRS.resonances.resolved.evaluated)
        rrr = resCls(endfRS.resonances.resolved.evaluated, verbose=args.verbose)
        rrr.setResonanceParametersByChannel()

        # Get the URR parameters
        if not hasattr(endfRS.resonances, 'unresolved'):
            raise ValueError("Evaluation %s does not have an URR" % endfFileName)
        if args.verbose:
            print(winged_banner("Extracting URR parameters"))
        urr = RRReconstruct.URRcrossSection(endfRS.resonances.unresolved.evaluated)
        egrid, flag = urr.generateEnergyGrid()
        urr.getWidthsAndSpacings()

    # -----------  Get the target and projectile parameters -----------
    targ, proj = get_target_and_projectile(config, endfRS)

    # -----------  Get spin group -----------
    spingroups = get_spin_groups(config, urr, targSpin=targ.spin, projSpin=proj.spin)

    # -----------  Get spacings and/or level densities -----------
    LD = get_level_densities(spingroups, config, urr)

    # -----------  Get widths -----------
    aveWidths = get_widths(spingroups, config, urr)

    # -----------  Get DOFs -----------
    DOFs = get_DOFs(spingroups, config, urr)

    # -----------  Get energy range to use -----------
    Emin, Emax = get_energy_range(config, urr, rrr, args.appendRRR, args.overwriteRRR)

    # ----------- Main part -----------
    fakeRR = {}
    print('\n' + banner("Generating resonance set"))
    if True:
        print("Target:", targ)
        print("Projectile:", proj)
        print("Allowed spingroups:", ", ".join([str(x) for x in spingroups]))

    # First time through, decide what to do with existing resonances
    cut_index = None
    if args.appendRRR or args.overwriteURR:  # Append new generated resonances to existing RRR
        cut_index = len(endfRS.resonances.resolved[0].resonanceParameters.table.data)
        endfRS.resonances.resolved.domainMax = config['Emax']
    if args.overwriteRRR:
        cut_index = 0  # Steam roll existing RRR with generated ones

    for sg in spingroups:

        if args.verbose:
            print(winged_banner("Doing SpinGroup %s" % str(sg)))
        else:
            print("Doing SpinGroup %s" % str(sg))

        # Set starting energy of fake resonance region
        if args.appendRRR or args.overwriteURR:
            startEnergy = get_last_resonance_energy(rrr, *sg.lj)
            if startEnergy is None:
                startEnergy = rrr.upperBound
        elif args.overwriteRRR:
            startEnergy = get_first_resonance_energy(rrr, *sg.lj)
            if startEnergy is None:
                startEnergy = rrr.lowerBound
        else:
            startEnergy = 0.0

        # Generate the resonances
        fakeRR[sg.lj] = rg.getFakeResonanceSet(
            E0=startEnergy,
            aveD=None,
            numLevels=None,
            style=config['spacingDist'],
            L=sg.L,
            J=sg.J,
            levelDensity=LD.get(sg.lj, None),
            aveWidthFuncs=aveWidths.get(sg.lj, None),
            DOFs=DOFs.get(sg.lj, None),
            widthKeys=('elastic', 'capture'),
            domainMin=Emin,
            domainMax=Emax,
            verbose=args.verbose)

        # Add resonances
        if args.appendRRR or args.overwriteURR or args.overwriteRRR:
            add_resonances_to_suite(fakeRR[sg.lj], endfRS, sg.lj[0], sg.lj[1], _cut_index=cut_index)
            cut_index = None

    # Remove the existing URR region
    if args.removeURR or args.overwriteURR:
        endfRS.resonances.unresolved = None

    # Just print out the resonances
    if args.list:
        # Build & sort list of resonances
        data = []
        for x in fakeRR.values():
            data += x.data
        data.sort(key=lambda res: res[0])

        # Assemble final table of resonances
        result = tableModule.Table(
            columns=[
                tableModule.ColumnHeader(0, name="energy", unit="eV"),
                tableModule.ColumnHeader(1, name="L", unit=""),
                tableModule.ColumnHeader(2, name="J", unit=""),
                tableModule.ColumnHeader(4, name="elastic width", unit="eV"),
                tableModule.ColumnHeader(5, name="capture width", unit="eV")],
            data=data)

        print("\n" + banner("Set of Fake Resonances"))
        print('\n'.join(result.toXML_strList()))  # toStringList()

    # Output to evaluation
    if args.outFile is not None:
        if args.outFormat == 'gnds':
            open(args.outFile, mode='w').write('\n'.join(endfRS.toXML_strList()))
        else:
            FIRSTLINE = " $Rev:: 700      $  $Date:: 2016-01-13#$                             1 0  0    0"
            outLines = endfRS.toENDF6('eval', {'verbosity': args.verbose * 10}, covarianceSuite=endfCS).split('\n')
            endfTxt = '\n'.join([FIRSTLINE, outLines[1].replace('-1', ' 8')] + outLines[2:])
            open(args.outFile.replace('.gnds.xml', '.endf'), mode='w').write(endfTxt)
