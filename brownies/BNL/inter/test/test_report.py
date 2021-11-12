import unittest
import sys
from fudge import reactionSuite
from brownies.BNL.inter.report import *
from xData.unittest import TestCaseWithIsClose as TestCaseWithExtras

TEST_DATA_PATH, this_filename = os.path.split(os.path.realpath(__file__))
VERBOSE = '-v' in sys.argv or '--verbose' in sys.argv
DOPLOTS = '-p' in sys.argv or '--plots' in sys.argv
if DOPLOTS:
    if '-p' in sys.argv:
        sys.argv.remove('-p')
    if '--plots' in sys.argv:
        sys.argv.remove('--plots')

if VERBOSE:
    print('reading test data...')

with open(TEST_DATA_PATH + os.sep + 'n-001_H_001.endf.gnds.xml', mode='r') as H1file:
    H1TestData = reactionSuite.readXML(H1file)
with open(TEST_DATA_PATH + os.sep + 'n-026_Fe_056.endf.gnds.xml', mode='r') as Fe56file:
    Fe56TestData = reactionSuite.readXML(Fe56file)
with open(TEST_DATA_PATH + os.sep + 'n-092_U_241.endf.gnds.xml', mode='r') as U241file:
    U241TestData = reactionSuite.readXML(U241file)


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
        self.example['h'] = {5, 6}
        self.example['i'] = collections.OrderedDict(((2, 'k'), ('k', 'l')))
        self.example['j'] = {43: 'k', 'k': 'l'}
        self.example['k'] = Report(title='empty')
        self.example['l'] = PQU.PQU(1.2, 'keV')

    def test_text_report(self):
        self.assertEqual(self.example.text_title(), '+-----------+\n| test data |\n+-----------+')
        self.assertJSONEqual(self.example.text_report(),
                             '{\n'
                             '    "test data": {\n'
                             '        "a": 4,\n'
                             '        "b": true,\n'
                             '        "c": 3.2,\n'
                             '        "d": "fred",\n'
                             '        "e": null,\n'
                             '        "f": [\n'
                             '            1,\n'
                             '            2\n'
                             '        ],\n'
                             '        "g": [\n'
                             '            3,\n'
                             '            4\n'
                             '        ],\n'
                             '        "h": [\n'
                             '            5,\n'
                             '            6\n'
                             '        ],\n'
                             '        "i": {\n'
                             '            "2": "k",\n'
                             '            "k": "l"\n'
                             '        },\n'
                             '        "j": {\n'
                             '            "k": "l",\n'
                             '            "43": "k"\n'
                             '        },\n'
                             '        "k": {},\n'
                             '        "l": "1.2 keV"\n'
                             '    }\n}')

    def test_json_report(self):
        jsonA = self.example.json_report()
        jsonB = '{"test data": {"a": 4, "b": true, "c": 3.2, "d": "fred", "e": null, "f": [1, 2], '\
                '"g": [3, 4], "h": [5, 6], "i": {"2": "k", "k": "l"}, "j": {"k": "l", "43": "k"}, '\
                '"k": {}, "l": "1.2 keV"}}'
        self.assertJSONEqual(jsonA, jsonB)

    def test_properties(self):
        self.assertEqual(self.example.id, id(self.example))
        self.assertEqual(self.example.cls, 'report testdata')

    @unittest.skip("FIXME:didn't get working yet")
    def test_html_report(self):
        self.assertEqual(self.example.html_report(), '')


class TestPlot(unittest.TestCase):

    def test_init(self):
        Plot()

    def test_properties(self):
        a = Plot(title='just a plot test')
        self.assertEqual(a.id, id(a))
        self.assertEqual(a.cls, 'plot justaplottest')


class TestGetReports(TestCaseWithExtras):

    def test_getEvaluationReport(self):
        self.assertJSONEqual(getEvaluationReport(Fe56TestData, title='n-026_Fe_056.endf.gnds.xml').text_report(), """{
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
    "GNDS version": "1.9"
}
}""")

    def test_getFissionMetricReport(self):
        metricMenu = collections.namedtuple('metricMenu',
                                            "ALF ETA")  # emulates results from ArgumentParser.parse_args() function
        rpt = getFissionMetricReport(U241TestData, metricMenu(ALF=True, ETA=True))["Fission metrics"].text_report()
        a = json.loads(rpt)
        b = json.loads("""{
"Fission metrics": {
    "ALF": "1143.149783603068",
    "ETA": "1.952768515392308e-3"
}
}""")
        for k in a['Fission metrics']:
            self.assertWithinXPercent(float(a['Fission metrics'][k]), float(b['Fission metrics'][k]), .01)

    def test_getResonanceReport(self):
        self.maxDiff = None
        self.assertJSONAlmostEqual(getResonanceReport(H1TestData).text_report(), """{
"Resonances": {
    "RRR information": {},
    "URR information": {},
    "Particle pairs": {}
}
}""")
        self.assertJSONAlmostEqual(getResonanceReport(Fe56TestData).text_report(), """{
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
""")
        self.assertJSONAlmostEqual(getResonanceReport(U241TestData).text_report(), """{
"Resonances": {
    "RRR information": {
        "Format": "BreitWigner",
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
        self.assertJSONEqual(getRRRReport(H1TestData).text_report(), """{}""")
        self.assertJSONEqual(getRRRReport(U241TestData).text_report(), """{
            "RRR information": {
                "Format": "BreitWigner",
                "Upper bound": 102.5,
                "Lower bound": 1.0e-5
            }
        }""")
        k = "RRR information"
        a = json.loads(getRRRReport(Fe56TestData).text_report())
        b = json.loads("""{
            "RRR information": {
                "Format": "RMatrix",
                "Upper bound": 8.5e5,
                "Lower bound": 1.0e-5,
                "Scattering length (R\')": "5.436999999987528 fm",
                "Lmax": 4,
                "No. channels": 5,
                "No. resonances": 313
            }
        }""")
        for kk in a[k]:
            if kk == "Scattering length (R\')":
                self.assertWithinXPercent(float(a[k][kk].replace('fm', '')), float(b[k][kk].replace('fm', '')),
                                          0.01)
            elif 'bound' in kk:
                self.assertWithinXPercent(a[k][kk], b[k][kk], 0.01)
            else:
                self.assertEqual(a[k][kk], b[k][kk])

    def test_getURRReport(self):
        self.assertJSONEqual(getURRReport(U241TestData).text_report(), """{
"URR information": {
    "Format": "unresolved",
    "Upper bound": 10000.0,
    "Lower bound": 102.5
}
}""")

    def test_getParticlePairsReport(self):
        self.assertJSONAlmostEqual(getParticlePairsReport(Fe56TestData).text_report(), """{
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
}""")

    # @unittest.skip("Need to discuss what channels should be there since LRF=3 translated into GNDS analog of LRF=7")
    def test_getChannelDataTable(self):
        metricMenu = collections.namedtuple('metricMenu',
                                            'effectiveDOF strengthFunction poleStrengthPlot scatteringRadiusPlot '
                                            'PorterThomas staircasePlot aveSpacingPlot aveWidthPlot shiftPlot '
                                            'penetrabilityPlot phasePlot DysonMehtaDelta3 transmissionCoeff')
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
        self.assertStringsEqual('\n'.join(x.toStringList()),
                                '<!--                           name | No. resonances | No. resonances w/ ER<0 | gfact |       Threshold E | Eliminated? | Competative? | Relativistic? | Pot. scatt. only? |             RRR <D> |        RRR <Gamma> |         RRR DOF  -->\n'
                                '<!--                                |                |                        |       |                eV |             |              |               |                   |                  eV |                 eV |                  -->\n'
                                '    Fe57 + photon (j=0.5,l=0,s=0.5)               40                        3       1   -7784318.76576425          True          False           False               False    2.5e4 +/- 1.2e4 eV     1.04 +/- 0.47 eV     0.84 +/- 0.41\n'
                                '         n + Fe56 (j=0.5,l=0,s=0.5)               40                        3       1                  -0         False          False           False               False    2.5e4 +/- 1.2e4 eV   3.8e3 +/- 5.9e3 eV   1.408 +/- 0.088\n'
                                '    Fe57 + photon (j=0.5,l=1,s=0.5)               63                        0       1   -7784318.76576425          True          False           False               False    1.3e4 +/- 1.3e4 eV     0.54 +/- 0.15 eV     0.62 +/- 0.28\n'
                                '         n + Fe56 (j=0.5,l=1,s=0.5)               63                        0       1                  -0         False          False           False               False    1.3e4 +/- 1.3e4 eV    109. +/- 2.5e2 eV   0.826 +/- 0.095\n'
                                '    Fe57 + photon (j=1.5,l=1,s=0.5)               79                        0       1   -7784318.76576425          True          False           False               False   1.01e4 +/- 8.0e3 eV     0.51 +/- 0.21 eV     0.80 +/- 0.30\n'
                                '         n + Fe56 (j=1.5,l=1,s=0.5)               79                        0       2                  -0         False          False           False               False   1.01e4 +/- 8.0e3 eV    152. +/- 2.5e2 eV     0.41 +/- 0.15\n'
                                '    Fe57 + photon (j=1.5,l=2,s=0.5)               75                        0       1   -7784318.76576425          True          False           False               False    1.1e4 +/- 1.1e4 eV     0.73 +/- 0.16 eV     0.52 +/- 0.30\n'
                                '         n + Fe56 (j=1.5,l=2,s=0.5)               75                        0       2                  -0         False          False           False               False    1.1e4 +/- 1.1e4 eV     73. +/- 1.2e2 eV   0.423 +/- 0.081\n'
                                '    Fe57 + photon (j=2.5,l=2,s=0.5)               56                        0       1   -7784318.76576425          True          False           False               False    1.5e4 +/- 1.3e4 eV     0.80 +/- 0.28 eV     0.69 +/- 0.29\n'
                                '         n + Fe56 (j=2.5,l=2,s=0.5)               56                        0       3                  -0         False          False           False               False    1.5e4 +/- 1.3e4 eV    171. +/- 2.6e2 eV   0.735 +/- 0.075')

    @unittest.skipIf(not DOPLOTS, "Not really a unit test, really more of an interactive test")
    def test_getChannelPlots(self):
        metricMenu = collections.namedtuple('metricMenu',
                                            'poleStrengthPlot scatteringRadiusPlot PorterThomas staircasePlot '
                                            'aveSpacingPlot aveWidthPlot backgroundXSPlot penetrabilityPlot '
                                            'shiftPlot phasePlot DysonMehtaDelta3 spacingDistributionPlot '
                                            'transmissionCoeff')

        def channelMatchFunction(c):
            return True  # return not c.eliminated and int(2.*c.J)==1 and c.l==0 and c.channelClass == 1

        getChannelPlots(Fe56TestData, metricMenu(poleStrengthPlot=False,
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
                                                 transmissionCoeff=False),
                        channelMatchFunction=channelMatchFunction, verbose=True)

    def test_getReactionLegacyINTERReport(self):
        a = json.loads(getReactionLegacyINTERDataTable(H1TestData).text_report())
        b = "<!--       name |  MT |       Q | Threshold |                RI |     14.2 MeV |          " \
            "2200 m/s |   Westcott factor |    252Cf (analytic) |               252Cf |     MACS(30. keV)  -->" \
            "\n<!--            |     |      eV |        eV |                 b |            b |                 " \
            "b |                   |                   b |                   b |                mb  -->\n" \
            "         n + H1     2         0           0     263.83668978541     0.67806042            20.43633    " \
            "1.12837868958532       3.9270475011023      3.66801927980566    14533.7414524736\n    H2 + photon   " \
            "102   2224631           0   0.149135213422371   2.9719802e-5   0.332281145255188   0.707206052903381   " \
            "3.94192612648554e-5   3.86144745081925e-5   0.152449777631191\n          total     1         0 " \
            "          0    263.983267844317     0.67809014    20.7687319650388    1.12196774546104      " \
            "3.92708606339206      3.66805680325653    14533.8730068349"
        alines = a["Legacy INTER metrics"].split('\n')
        blines = b.split('\n')
        for i in range(2, len(alines)):
            aline = alines[i].split()[-10:]
            bline = blines[i].split()[-10:]
            for j in range(len(aline)):
                self.assertWithinXPercent(float(aline[j]), float(bline[j]), 0.01)


# -------------------------------------------------------------------------------
# Run unit tests
# -------------------------------------------------------------------------------

if __name__ == "__main__":

    unittest.main()
