#! /usr/bin/env python

import argparse, copy, math
from pqu import PQU
from fudge.core.utilities.brb import *
import fudge.processing.resonances.reconstructResonances as resReconstructModule

NUMERICAL_TYPES = [float, int, PQU.PQU]


def parse_args():
    """Process command line options"""
    parser = argparse.ArgumentParser(description='Reconstruct the angular distributions from an ENDF file')
    parser.add_argument('a', help='first file to diff')
    parser.add_argument('b', help='second file to diff')
    parser.add_argument('-v', dest='verbose', default=False, action='store_true',
                        help='Enable verbose output (Default: False)')
    parser.add_argument('-s', dest='skipBadData', default=False, action='store_true',
                        help='Skip bad data on ENDF read (Default: False)')
    parser.add_argument('-w', dest='warnOnly', default=False, action='store_true',
                        help="Don't fail, just warn when an illegal spin combination is encountered")
    parser.add_argument('-r', dest='reconstructResonances', default=False, action='store_true',
                        help='Reconstruct resonances on read (Default: False)')
    parser.add_argument('-m', dest='multipleSScheme', default='NJOY', choices=['ENDF', 'NJOY', 'Ignore'],
                        help="What scheme to use to set the S of a channel when there are multiple possibilities")
    return parser.parse_args()


def read_evaluation(inFile, skipBadData=True, verbose=False, reconstructResonances=False):
    """
    Read in the evaluation
    :param inFile: the file name.  This routine will attempt to determine the file type itself.
    :param skipBadData: FIXME
    :param verbose: FIXME
    :param reconstructResonances: FIXME
    :return: the evaluation as a dict of results like {'reactionSuite':the reaction suite, 'covarianceSuite':
             the covariances}
    """
    if not os.path.exists(inFile):
        raise IOError("File named " + inFile + " doesn't exist!")
    if open(inFile).readline().startswith("<?xml"):
        from fudge import reactionSuite
        return {'reactionSuite': reactionSuite.readXML(inFile), 'covarianceSuite': None, 'info': {}, 'errors': []}
    else:
        from brownies.legacy.converting import endfFileToGNDS
        return endfFileToGNDS.endfFileToGNDS(inFile, toStdOut=verbose, skipBadData=skipBadData,
                                             reconstructResonances=reconstructResonances)


def print_comparison(aVal, bVal, itemName):
    # Numerical comparison
    if type(aVal) in NUMERICAL_TYPES and type(bVal) in NUMERICAL_TYPES:
        if aVal > bVal:
            print(itemName, 'values differ: a.val = %s > b.val = %s' % (str(aVal), str(bVal)))
        elif aVal < bVal:
            print(itemName, 'values differ: a.val = %s < b.val = %s' % (str(aVal), str(bVal)))
        else:
            print(itemName, 'values equal:', aVal)
    # For non-numerical comparisons, must make sure objects of same type
    elif type(aVal) != type(bVal):
        print(itemName, 'types differ: type(a)=%s, type(b)=%s' % (str(type(a)), str(type(b))))
    elif aVal != bVal:
        print(itemName, 'values differ:', aVal, 'vs.', bVal)
    else:
        print(itemName, 'values equal:', aVal)


def get_rrr_type(reactionSuite):
    t = reactionSuite.resonances.resolved.evaluated.moniker
    if hasattr(reactionSuite.resonances.resolved.evaluated, "approximation"):
        t += ', ' + reactionSuite.resonances.resolved.evaluated.approximation + ' approximation'
    return t


def get_rrr(reactionSuite, verbose=False, warnOnly=False, multipleSScheme='ENDF'):
    resCls = resReconstructModule.getResonanceReconstructionClass(reactionSuite.resonances.resolved.evaluated)
    rrr = resCls(reactionSuite.resonances.resolved.evaluated, enableAngDists=False, verbose=verbose)
    if (isinstance(rrr, resReconstructModule.RMatrixLimitedcrossSection) and
            (rrr.RR.approximation == fudge.resonances.resolved.RMatrix.ReichMooreToken)):
        rrr.setResonanceParametersByChannel(useReichMooreApproximation=True, warnOnly=warnOnly,
                                            multipleSScheme=multipleSScheme)
    else:
        rrr.setResonanceParametersByChannel(warnOnly=warnOnly)
    return rrr


def channel_details(c, indent=8, verbose=False):
    return indent * ' ' + "%2i Reaction %s, J,L,s=%s,%s,%s, Xi=%s" % \
           (c.index, c.reaction.rjust(13), str(c.J), str(c.l), str(c.s), str(c.Xi))


def resonance_details(fname, reactionSuite, rrrReconstructor, indent=4, verbose=False):
    def getAP(_rrr):
        _res = {}
        for attr in ['getScatteringLength']:
            if hasattr(_rrr, attr):
                _res[attr.replace('get', '')] = str(eval('_rrr.' + attr + "()"))
        return _res

    resultList = [
        'Name: %s' % fname,
        'Target: %s' % reactionSuite.target,
        'Projectile: %s' % reactionSuite.projectile,
        'Compound nucleus formed: %s' % reactionSuite.compound_nucleus,
        'Format: %s' % get_rrr_type(reactionSuite),
        'Num. channels: %i' % rrrReconstructor.nChannels,
        '\n'.join(
            [channel_details(x[0], indent=int(x[1] != 0) * 4 + 4) for x in
             [(z, i) for i, z in enumerate(rrrReconstructor.channels)]]),
        'Num. resonances: %i' % rrrReconstructor.nResonances,
        'LMax: %i' % rrrReconstructor.getLMax(),
        'LowerBound: %s' % str(reactionSuite.resonances.resolved.domainMin),
        'UpperBound: %s' % str(reactionSuite.resonances.resolved.domainMax)]
    APstuff = getAP(rrrReconstructor)
    for k in APstuff:
        resultList.append('%s : %s' % (k, APstuff[k]))
    result = '\n'.join([indent * ' ' + x for x in resultList])
    result.replace('PUTCHANNELSUMMARYHERE', '')
    return result


def channels_match(ca, cb):
    # if ca==cb:
    #    print(ca , hash(ca) )
    #    print(cb , hash(cb) )
    #    print()
    return ca == cb


def compare_channels(_aRRR, _bRRR, indent=4, verbose=False):
    matches = []
    ina_notb = []
    inb_nota = []

    # First pass, look for matches and things only in "a"
    for ca in _aRRR.channels:
        gotMatch = False
        for cb in _bRRR.channels:
            if channels_match(ca, cb):
                matches.append(ca)
                gotMatch = True
                break
        if not gotMatch:
            ina_notb.append(ca)

    # Second pass, look for things only in "b"
    for cb in _bRRR.channels:
        gotMatch = False
        for ca in _aRRR.channels:
            if channels_match(ca, cb):
                gotMatch = True
                break
        if not gotMatch:
            inb_nota.append(cb)

    if verbose:
        print('Matched channels:\n', '\n'.join([indent * ' ' + str(x) for x in matches]))
        if ina_notb == [] and inb_nota == []:
            print("All channels matched")
        if ina_notb:
            print('Channels in a only:\n', '\n'.join([indent * ' ' + str(x) for x in ina_notb]))
        if inb_nota:
            print('Channels in b only:\n', '\n'.join([indent * ' ' + str(x) for x in inb_nota]))

    return {'matches': matches, 'ina_notb': ina_notb, 'inb_nota': inb_nota}


def compare_resonances(_aRRR, _bRRR, _channelComparison, abs_tol_ER=20.0, rel_tol_ER=0.001, indent=4, verbose=False):
    upBd = min(_aRRR.upperBound, _bRRR.lowerBound)
    for case in _channelComparison:
        if case in ['ina_notb', 'inb_nota']:
            if case == 'ina_notb':
                print('Channels in a only:')
                check_this = _aRRR
            if case == 'inb_nota':
                print('Channels in b only:')
                check_this = _bRRR

            oops = []
            for c in _channelComparison[case]:
                if c.channelClass == resReconstructModule.GAMMACHANNEL:
                    # Check if channel really, for sure, is not eliminated
                    if c.eliminated:
                        pass  # is OK
                    else:
                        oops.append('Gamma channel not eliminated: ' + str(c))
                else:
                    # Check if channels above upperBound
                    if c.Xi < upBd:
                        oops.append('Channel threshold below upperBound of regions: ' + str(c))
            if not oops:
                print(indent * ' ' + '...everything cool')
            else:
                print('\n'.join([indent * ' ' + oop for oop in oops]))

        else:
            import collections

            def entry_line_format(entry):
                return ('%s %s %s %s %s %s %s' % tuple([str(x).rjust(15) for x in entry[:-1]] + [entry[-1]])).replace(
                    'None', '    ')

            def get_matching_gamma(iER, gChannels):
                for _c in gChannels:
                    for res in gChannels[_c].items():
                        if res[0] == iER:
                            return res[1]

            def percent_diff(A, B):
                ave = 0.5 * (A + B)
                if ave == 0.0:
                    return 0.0
                return 100.0 * (A - B) / ave

            print('...Building effective spin groups for comparing matching channels')
            spinGroups = collections.OrderedDict()
            for c in channelComparison[case]:
                JSL = (c.J, c.s, c.l)
                if JSL not in spinGroups:
                    spinGroups[JSL] = {'channels': [], 'resonances': []}
                spinGroups[JSL]['channels'].append(c)

            print('...Comparing resonances within each effective spin group')
            for JSL in spinGroups:
                for c in spinGroups[JSL]['channels']:
                    if c.channelClass == resReconstructModule.GAMMACHANNEL:
                        continue
                    aResonances = copy.copy(list(_aRRR.channels[c].items()))
                    bResonances = copy.copy(list(_bRRR.channels[c].items()))
                    aResonances.sort()
                    bResonances.sort()
                    # first pass through a first, then b
                    bMatches = []
                    while aResonances:
                        aiR, aGamn = aResonances.pop(0)
                        aGammaChannels = _aRRR.eliminatedChannels
                        if not aGammaChannels:
                            aGammaChannels = {c: _aRRR.channels[c] for c in _aRRR.channels if
                                              c.channelClass == resReconstructModule.GAMMACHANNEL}
                        aGamg = get_matching_gamma(aiR, aGammaChannels)
                        aER = _aRRR._energies[aiR]
                        aSHF = _aRRR.shiftFactorByChannel(c, aER)
                        if aER > 1.1 * upBd:
                            continue
                        gotMatch = False
                        for ib, bPair in enumerate(bResonances):
                            if ib in bMatches:
                                continue
                            biR, bGamn = bPair
                            bER = _bRRR._energies[biR]
                            if math.isclose(aER, bER, abs_tol=abs_tol_ER, rel_tol=rel_tol_ER):
                                gotMatch = True
                                bGammaChannels = _bRRR.eliminatedChannels
                                if not bGammaChannels:
                                    bGammaChannels = {c: _bRRR.channels[c] for c in _bRRR.channels if
                                                      c.channelClass == resReconstructModule.GAMMACHANNEL}
                                bGamg = get_matching_gamma(biR, bGammaChannels)
                                bMatches.append(ib)
                                bSHF = _bRRR.shiftFactorByChannel(c, bER)
                                dSHF = 100.0 * (aSHF - bSHF) * 2.0 / (aER + bER)
                                dER = percent_diff(aER, bER)
                                dGamn = percent_diff(aGamn, bGamn)
                                dGamg = percent_diff(aGamg, bGamg)
                                comment = 'dER = %4.2f %% (SHF=%4.2f eV), dGam_n = %4.2f %%, gGam_g = %4.2f %%' % (
                                    dER, bSHF, dGamn, dGamg)
                                if math.isclose(aER, bER, abs_tol=1.0):
                                    comment = ''
                                spinGroups[JSL]['resonances'].append(
                                    [aER, aER, aGamn, aGamg, bER, bGamn, bGamg, comment])
                        if not gotMatch:
                            spinGroups[JSL]['resonances'].append(
                                [aER, aER, aGamn, aGamg, None, None, None, "In a only"])

                    bMatches.reverse()
                    for ib in bMatches:
                        bResonances.pop(ib)

                    # second pass, start with b, just to catch things in b but not a
                    while bResonances:
                        biR, bGamn = bResonances.pop(0)
                        bGamg = get_matching_gamma(biR, _bRRR.eliminatedChannels)
                        bER = _bRRR._energies[biR]
                        if bER > 1.1 * upBd:
                            continue
                        spinGroups[JSL]['resonances'].append([bER, None, None, None, bER, bGamn, bGamg, "In b only"])
                        spinGroups[JSL]['resonances'].sort()

            print("...Here's what I found:")
            JSLSorted = list(spinGroups.keys())

            def jslSorter(x):
                return x[2]*10000 + x[1]*100 + x[1]

            JSLSorted.sort(key=jslSorter)
            for JSL in JSLSorted:
                print("\nEffective spin group (J,s,L) =", JSL, ':')
                print('\n'.join(['    ' + str(c) for c in spinGroups[JSL]['channels']]))
                for c in spinGroups[JSL]['channels']:
                    print(2 * indent * ' ' + 'AP for a:', _aRRR.getAPByChannel(c))
                    print(2 * indent * ' ' + 'AP for b:', _bRRR.getAPByChannel(c))
                print()
                print(2 * indent * ' ' + '--------------------- a -----------------------   '
                                         '--------------------- b --------------------- ')
                print(2 * indent * ' ' + entry_line_format(['ER (eV)', 'Gamn (eV)', 'Gamg (eV)', 'ER (eV)',
                                                            'Gamn (eV)', 'Gamg (eV)', 'Comment']))
                print(2 * indent * ' ' + 120 * '=')
                print('\n'.join([2 * indent * ' ' + entry_line_format(x[1:]) for x in spinGroups[JSL]['resonances']]))

        print()


if __name__ == "__main__":
    args = parse_args()

    print()
    print(winged_banner("Reading %s" % args.a))
    a = read_evaluation(args.a, verbose=args.verbose, skipBadData=args.skipBadData,
                        reconstructResonances=args.reconstructResonances)
    aRRR = get_rrr(a['reactionSuite'], verbose=args.verbose, warnOnly=args.warnOnly,
                   multipleSScheme=args.multipleSScheme)

    print()
    print(winged_banner("Reading %s" % args.b))
    b = read_evaluation(args.b, verbose=args.verbose, skipBadData=args.skipBadData,
                        reconstructResonances=args.reconstructResonances)
    bRRR = get_rrr(b['reactionSuite'], verbose=args.verbose, warnOnly=args.warnOnly,
                   multipleSScheme=args.multipleSScheme)

    print()
    print(winged_banner("Comparing RRR"))
    print()
    print('Left (a) file details:')
    print(resonance_details(args.a, a['reactionSuite'], aRRR))
    print()
    print('Right (b) file details:')
    print(resonance_details(args.b, b['reactionSuite'], bRRR))
    print()
    print_comparison(a['reactionSuite'].projectile, b['reactionSuite'].projectile, 'Projectile')
    print_comparison(a['reactionSuite'].target, b['reactionSuite'].target, 'Target')
    print_comparison(a['reactionSuite'].resonances.resolved.domainMin,
                     b['reactionSuite'].resonances.resolved.domainMin, 'domainMin')
    print_comparison(a['reactionSuite'].resonances.resolved.domainMax,
                     b['reactionSuite'].resonances.resolved.domainMax, 'domainMax')
    print_comparison(aRRR.nResonances, bRRR.nResonances, 'Num. resonances')
    print_comparison(aRRR.nChannels, bRRR.nChannels, 'Num. channels')

    print()
    channelComparison = compare_channels(aRRR, bRRR, verbose=args.verbose)

    print()
    compare_resonances(aRRR, bRRR, channelComparison, verbose=args.verbose)
