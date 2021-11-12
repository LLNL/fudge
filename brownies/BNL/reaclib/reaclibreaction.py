#! /usr/bin/env python

import collections
import io

import math
import numpy as np

OneKelvinInKeV = 0.00000008617328149741  # 1 K = 8.617328149741e-8 keV
T9inKeV = 1e9 * OneKelvinInKeV  # 1 GK = 86.17328149741 keV


class ReaclibLibrary:
    def __init__(self, s=None):
        self.chapters = {}
        if s is not None:
            self.parse(s)

    def parse(self, s):
        lines = io.StringIO(s).readlines()
        chunks = {}
        thischunk = ''
        chapter = 0
        while lines:
            line = lines.pop(0)
            if (line[:1]).isdigit():
                if thischunk.strip() != '':
                    chunks.setdefault(chapter, []).append(thischunk)
                chapter = int(line[:2])
                thischunk = ''
            thischunk += line
        if thischunk.strip() != '':
            chunks.setdefault(chapter, []).append(thischunk)
        for i in chunks:
            # print( ''.join(chunks[i]) )
            w = ReaclibChapter(''.join(chunks[i]))
            self.chapters[i] = w

    def __str__(self):
        result = ''
        for i in self.chapters:
            result += str(self.chapters[i])
        return result

    def find(self, targ=None, proj=None, reaction=None, reference=None, chapter=None):
        for chapter in self.chapters:
            pass
        raise NotImplementedError("FIXME: write me")

    def add_inverse_reactions_detailed_balance(self):
        raise NotImplementedError("FIXME: write me")


class ReaclibChapter:
    def __init__(self, s=None):
        self.reacts = []
        self.chapter = None
        if s is not None:
            self.parse(s)

    def parse(self, s):
        split = s.splitlines()
        nchunks = len(split) / 4
        for ichunk in range(0, len(split), 4):
            chunk = split[ichunk:ichunk + 4]
            header = chunk.pop(0)
            set1 = chunk.pop(0)
            set2 = chunk.pop(0)
            set3 = chunk.pop(0)
            r = ReaclibReact(set1 + set2 + set3)
            self.reacts.append(r)
            self.chapter = header

    def __str__(self):
        mymap2 = ''
        for r in self.reacts:
            mymap2 += '\n'.join(list(map(str, [self.chapter, r]))) + '\n'
        return mymap2


class ReaclibReact:
    def __init__(self, s=None):
        self.original = s
        self.set_label = None
        self.rate = None
        self.reverserate = None
        self.Qval = None
        self.e = collections.OrderedDict(((i, None) for i in range(1, 6)))
        if s is not None:
            self.parse(s)

    def add_particle(self, p, slot):
        pstring = p.rjust(5)
        self.e[slot] = pstring

    def add_coefficient(self, c, slot):
        if slot == 0:
            self.a_zero = ['', c]
        elif slot == 1:
            self.a_one = ['', c]
        elif slot == 2:
            self.a_two = ['', c]
        elif slot == 3:
            self.a_three = ['', c]
        elif slot == 4:
            self.a_four = ['', c]
        elif slot == 5:
            self.a_five = ['', c]
        elif slot == 6:
            self.a_six = ['', c]
        else:
            raise ValueError("Coefficient slots are 0-6, you entered %s" % str(slot))

    def parse(self, s):
        self.e[1] = s[5:10]
        self.e[2] = s[10:15]
        self.e[3] = s[15:20]
        self.e[4] = s[20:25]
        self.e[5] = s[25:30]
        self.e[6] = s[30:35]
        self.set_label = s[43:47]
        self.rate = s[47]
        self.reverserate = s[48]
        self.Qval = [str(s[52:64]), float(s[52:64])]
        #        print( '/'+s[74:87]+'/' )
        #        print( '/'+s[87:100]+'/' )
        self.a_zero = [str(s[74:87]), float(s[74:87])]
        self.a_one = [str(s[87:100]), float(s[87:100])]
        self.a_two = [str(s[100:113]), float(s[100:113])]
        self.a_three = [str(s[113:126]), float(s[113:126])]
        self.a_four = [str(s[148:161]), float(s[148:161])]
        self.a_five = [str(s[161:174]), float(s[161:174])]
        self.a_six = [str(s[174:187]), float(s[174:187])]

    def __str__(self):
        line1 = "     " + self.e[1] + self.e[2] + self.e[3] + self.e[4] + self.e[5] + self.e[
            6] + "        " + self.set_label + self.rate + self.reverserate + "   " + self.Qval[0] + "          "
        line2 = self.a_zero[0] + self.a_one[0] + self.a_two[0] + self.a_three[0] + "                      "
        line3 = self.a_four[0] + self.a_five[0] + self.a_six[0] + "                                   "
        mymap = list(map(str, [line1, line2, line3]))
        return "\n".join(mymap)

    def evaluate(self, x, unit='T9'):
        if unit == 'T9':
            T = x
        elif unit == 'keV':
            T = x / T9inKeV
        else:
            raise NotImplementedError("unit must be 'keV' or 'T9', got %s" % unit)
        logrrate = (
                self.a_zero[1] +
                (self.a_one[1] * pow(T, -1.0)) +
                (self.a_two[1] * pow(T, -0.33333333333)) +
                (self.a_three[1] * pow(T, 0.3333333333)) +
                (self.a_four[1] * pow(T, 1)) +
                (self.a_five[1] * pow(T, 1.66666666666)) +
                (self.a_six[1] * math.log(T)))
        return math.exp(logrrate)


def gnds_crossSection_to_Reaclib_ARR_coefficients(xs, useCovariance=False, verbose=False, minTemp=0.1, maxTemp=10.0,
                                                  numTemps=15):
    from brownies.BNL.inter.metrics import computeAstrophysicalReactionRate
    import pqu.PQU as pqu

    # Set up temperature grid
    Temp_inT9 = np.logspace(start=math.log(minTemp), stop=math.log(maxTemp),
                            num=numTemps)  # [.1, .2, .3, .4, .5, 1, 1.5, 2.0, 5.0, 10.0]
    Temp_inkeV = [T9inKeV * (T) for T in Temp_inT9]

    # Set up array of reaction rates
    b_matrix = []
    for T in Temp_inkeV:
        ARR = computeAstrophysicalReactionRate(xs, pqu.PQU(T, 'keV'), useCovariance=useCovariance)
        b_matrix.append(np.log(ARR.getValue()))
    b = np.array(b_matrix)

    # Set up matrix of powers of temperatures so we can fit the coefficients
    a_matrix = []
    for T in Temp_inT9:
        a_matrix_row = [1, pow(T, -1.0), pow(T, -1.0 / 3.0), pow(T, 1.0 / 3.0), pow(T, 1.0), pow(T, 5.0 / 3.0),
                        np.log(T)]
        a_matrix.append(a_matrix_row)
    a = np.array(a_matrix)

    # a*x = b, solve for vector x
    if verbose: print('b:', b)
    if verbose: print('a:', a)
    x, residual, rank, s = np.linalg.lstsq(a, b, rcond=None)
    if verbose: print('a*x', np.dot(a, x))
    if verbose: print(residual)
    return x


def gnds_reactionSuite_to_ReaclibLibrary(rs, maxEThreshold=0.5e6, skipCrossSectionSums=True, verbose=False):
    import fudge
    projectile = str(rs.projectile)
    target = str(rs.target)
    resultChapters = []
    for r in rs:
        # Get rid of all reactions that are not plain reactions or sums of reactions
        # (e.g. no fission components of production stuff)
        if not isinstance(r, (fudge.reactions.reaction.reaction, fudge.sums.crossSectionSum)): continue
        if isinstance(r, fudge.sums.crossSectionSum) and skipCrossSectionSums: continue
        if not hasattr(r, 'outputChannel'): continue

        # Compute outgoing particle names
        prods = [str(p) for p in r.outputChannel]
        reactionName = str(r)

        # Get Ethreshold
        if hasattr(r, "getThreshold"):
            EThreshold = r.getThreshold(unit='eV')
        else:
            EThreshold = 0.0

        # Skip over some reactions for which an astrophysical reaction rate is useless
        if projectile in prods and target in prods: continue  # skip elastic
        if EThreshold > maxEThreshold: continue  # skip high threshold reactions
        isInelastic = False
        for prod in prods:
            if '(' in prod:
                isInelastic = True  # skip anything that's a discrete level excitation, we're not ready for isomers
                break
        if isInelastic: continue

        print(20 * '-', str(r), 20 * '-')

        # Get Q
        if hasattr(r, 'getQ'):
            Q = r.getQ(unit='eV')
        else:
            Q = 0.0

        # Get coefficients in ARR parameterization
        a = gnds_crossSection_to_Reaclib_ARR_coefficients(r.crossSection, verbose=False)

        # Figure out relevant products and compute chapter
        if 'photon' in prods: prods.remove('photon')
        if len(prods) == 1:
            chapter = 4
        elif len(prods) == 2:
            chapter = 5
        elif len(prods) == 3:
            chapter = 6
        else:
            continue

        thisReaction = ReaclibReact()
        thisReaction.set_label = 'endf'
        thisReaction.rate = None
        thisReaction.reverserate = ''
        thisReaction.Qval = ['', Q]
        thisReaction.add_particle(target, 1)
        thisReaction.add_particle(projectile, 2)
        for ip, product in enumerate(prods): thisReaction.add_particle(product, ip + 3)
        for ic, c in enumerate(a): thisReaction.add_coefficient(c, ic)

        print(prods, reactionName, a, chapter, Q, EThreshold)
        #       print( thisReaction )
        thisChapter = ReaclibChapter()
        thisChapter.chapter = chapter
        thisChapter.reacts.append(thisReaction)
        resultChapters.append(thisChapter)
    return resultChapters
    raise NotImplementedError("FIXME: write me")


def endfFile_to_ReaclibLibrary(filename, verbose=False):
    from brownies.legacy.converting import endfFileToGNDS
    rs = endfFileToGNDS.endfFileToGNDS(filename,
                                       toStdOut=int(verbose) * 10,
                                       skipBadData=True,
                                       reconstructResonances=True,
                                       continuumSpectraFix=True,
                                       doCovariances=False,
                                       verboseWarnings=verbose,
                                       printBadNK14=False,
                                       ignoreBadDate=True,
                                       acceptBadMF10FissionZAP=True)['reactionSuite']
    return gnds_reactionSuite_to_ReaclibLibrary(rs, verbose=verbose)


def plot_rates(d, temperatureGrid=[.1, .2, .3, .4, .5, 1, 1.5, 2, 5, 10], title="Astrophysical Reaction Rate"):
    import matplotlib.pyplot as plt
    ''' d={ "ENDF/B-VIII.0":reactlibreac, 'KaDoNiS':...}'''
    for k in d:
        plt.plot(temperatureGrid, [d[k].evaluate(T) for T in temperatureGrid], label=k)
    plt.xlabel('$Temperature (GK)$')
    plt.ylabel('$ Astrophysical Reaction Rate(1e6)(cm^{3} s^{-1} mol^{-1})$')
    plt.title(title)
    plt.legend(loc='upper center', shadow=True, fontsize='x-large')
    plt.show()


if __name__ == "__main__":
    import unittest


    class Test_ReaclibReact(unittest.TestCase):
        def setUp(self):
            self.s = "         n    p                            wc12w     7.82300e-01          \n" \
                     "-6.781610e+00 0.000000e+00 0.000000e+00 0.000000e+00                      \n" \
                     " 0.000000e+00 0.000000e+00 0.000000e+00                                   \n"

        # self.s=open('reaclibv1.txt').read()
        def test_output(self):
            r = ReaclibReact(self.s)
            self.assertEqual(str(r), self.s)

        def test_evaluate(self):
            r = ReaclibReact(self.s)
            self.assertAlmostEqual(r.evaluate(1.0), 4.0)


    class Test_ReaclibChapter(unittest.TestCase):
        def setUp(self):
            self.s = \
                """1
                         n    p                            wc12w     7.82300e-01          
                -6.781610e+00 0.000000e+00 0.000000e+00 0.000000e+00                      
                 0.000000e+00 0.000000e+00 0.000000e+00                                   
                1
                         t  he3                            wc12w     1.86000e-02          
                -2.014560e+01 0.000000e+00 0.000000e+00 0.000000e+00                      
                 0.000000e+00 0.000000e+00 0.000000e+00                                   
                1
                       he3    t                              ecw    -1.90000e-02          
                -3.246200e+01-2.133800e-01-8.215810e-01 1.112410e+01                      
                -5.773380e-01 2.904710e-02-2.627050e-01                                   
                """

        def test_output(self):
            r = ReaclibChapter(self.s)
            self.assertEqual(str(r), self.s)


    class Test_ReaclibLibrary(unittest.TestCase):
        def setUp(self):
            self.s = ""

        def test_output(self):
            r = ReaclibLibrary(self.s)
            self.assertEqual(str(r), self.s)


    unittest.main()
