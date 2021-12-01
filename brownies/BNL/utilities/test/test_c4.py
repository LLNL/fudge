import unittest
from brownies.BNL.utilities.c4 import *


# ------------------------------------------------
# Unit tests
# ------------------------------------------------

class TestFunkyFloats(unittest.TestCase):
    def setUp(self):
        self.nums = {0.314159e14: "3.141590+13", -3.14159e13: "-3.14159+13", 51450.00: '51450      ',
                     -51450.00: '-51450     ', 0.00: '0          ', 504400.: "504400     ", 4961.493: "4961.493   ",
                     8.146: "8.146      ", 0.4042: "0.4042     "}

    def test_write_read_consistency(self):
        for x in self.nums:
            a = writeFunkyFloat(x)
            b = readFunkyFloat(a)
            self.assertEqual(x, b)

    def test_write(self):
        for k, v in self.nums.items():
            self.assertEqual(v, writeFunkyFloat(k))


class TestReadC4Point(unittest.TestCase):
    def setUp(self):
        self.a = C4Point(projectile=1, target=40092, targetMetastableState=None, MF=3, MT=1,
                         productMetastableState=None, status='A', cmFlag=None, energy=504400.0, dEnergy=4961.493,
                         data=8.146000, dData=0.404200, cosMuOrLegendreOrder=None, dCosMuOrLegendreOrder=None,
                         eLevelOrHalflife=None, dELevelOrHalflife=None, idOf78=None, reference='L.GREEN,ET.AL. (73)',
                         exforEntry='10225', exforSubEntry=20, multiDimFlag=None)

    def test_a(self):
        txt = '    1 40092   3   1 A 504400   4961.493 8.146    0.4042                                          ' \
              'L.GREEN,ET.AL. (73)      10225 20 '
        b = readC4Point(txt)
        self.assertEqual(self.a, b)
        c = writeC4Point(b)
        self.assertEqual(txt, c)
        self.assertEqual(self.a, readC4Point(c))


class TestReadC4Entry(unittest.TestCase):
    def setUp(self):
        self.a = C4Entry(
            entry='40617', author1='M.V.Pasechnik+', year=1980, institute='(4CCPIJI)',
            title="TOTAL NEUTRON CROSS-SECTIONS FOR MOLYBDENUM AND ZYRCONIUM AT LOW ENERGIES",
            authors="M.V.Pasechnik, M.B.Fedorov, V.D.Ovdienko, G.A.Smetanin, T.I.Jakovenko",
            refCode="(C,80KIEV,1,304,8009)",
            reference="Conf. 5.All Union Conf.on Neutron Phys.,Kiev,15-19 Sep 1980 Vol.1, p.304, 1980",
            numDataSets=7,
            dataSets=[
                readC4DataSet("""#DATASET    40617007
#DATE       19850305
#REACTION   40-ZR-92(N,TOT),,SIG
#PROJ       1
#TARG       40092
#MF         3
#MT         1
#C4BEGIN    [    1 40092   3   1 A ]
#DATA       4
# Prj Targ M MF MT PXC  Energy  dEnergy  Data      dData   Cos/LO   dCos/LO   ELV/HL  dELV/HL I78
#---><---->o<-><-->ooo<-------><-------><-------><-------><-------><-------><-------><-------><->
    1 40092   3   1 A  442000.0          12.74000 1.700000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  507000.0          8.790000 0.570000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  572000.0          9.520000 0.200000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  637000.0          8.480000 0.110000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
#/DATA      4
#/DATASET""".split('\n'))])

    def test_a(self):
        b = readC4Entry("""#ENTRY      40617
#AUTHOR1    M.V.Pasechnik+
#YEAR       1980
#INSTITUTE  (4CCPIJI)
#TITLE      TOTAL NEUTRON CROSS-SECTIONS FOR MOLYBDENUM
#+          AND ZYRCONIUM AT LOW ENERGIES
#AUTHOR(S)  M.V.Pasechnik, M.B.Fedorov, V.D.Ovdienko,
#+          G.A.Smetanin, T.I.Jakovenko
#REF-CODE   (C,80KIEV,1,304,8009)
#REFERENCE  Conf. 5.All Union Conf.on Neutron Phys.,Kiev,15-19 Sep 1980
#+          Vol.1, p.304, 1980
#DATASETS   7
#
#DATASET    40617007
#DATE       19850305
#REACTION   40-ZR-92(N,TOT),,SIG
#PROJ       1
#TARG       40092
#MF         3
#MT         1
#C4BEGIN    [    1 40092   3   1 A ]
#DATA       4
# Prj Targ M MF MT PXC  Energy  dEnergy  Data      dData   Cos/LO   dCos/LO   ELV/HL  dELV/HL I78
#---><---->o<-><-->ooo<-------><-------><-------><-------><-------><-------><-------><-------><->
    1 40092   3   1 A  442000.0          12.74000 1.700000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  507000.0          8.790000 0.570000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  572000.0          9.520000 0.200000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  637000.0          8.480000 0.110000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
#/DATA      4
#/DATASET
#/ENTRY
""".split('\n'))
        self.assertEqual(self.a, b)


class TestReadC4DataSet(unittest.TestCase):
    def setUp(self):
        self.a = C4DataSet(
            dataSet='40617007', date='19850305',
            reaction='40-ZR-92(N,TOT),,SIG', projectile=1, target=40092, MF=3, MT=1,
            c4Begin="[    1 40092   3   1 A ]", numData=4,
            data=[
                C4Point(projectile=1, target=40092, targetMetastableState=None, MF=3, MT=1, productMetastableState=None,
                        status='A', cmFlag=None, energy=442000.0, dEnergy=None, data=12.74000, dData=1.700000,
                        cosMuOrLegendreOrder=None, dCosMuOrLegendreOrder=None, eLevelOrHalflife=None,
                        dELevelOrHalflife=None, idOf78=None, reference='M.V.PASECHNIK,ET.AL. (80)', exforEntry='40617',
                        exforSubEntry=7, multiDimFlag=None),
                C4Point(projectile=1, target=40092, targetMetastableState=None, MF=3, MT=1, productMetastableState=None,
                        status='A', cmFlag=None, energy=507000.0, dEnergy=None, data=8.790000, dData=0.570000,
                        cosMuOrLegendreOrder=None, dCosMuOrLegendreOrder=None, eLevelOrHalflife=None,
                        dELevelOrHalflife=None, idOf78=None, reference='M.V.PASECHNIK,ET.AL. (80)', exforEntry='40617',
                        exforSubEntry=7, multiDimFlag=None),
                C4Point(projectile=1, target=40092, targetMetastableState=None, MF=3, MT=1, productMetastableState=None,
                        status='A', cmFlag=None, energy=572000.0, dEnergy=None, data=9.520000, dData=0.200000,
                        cosMuOrLegendreOrder=None, dCosMuOrLegendreOrder=None, eLevelOrHalflife=None,
                        dELevelOrHalflife=None, idOf78=None, reference='M.V.PASECHNIK,ET.AL. (80)', exforEntry='40617',
                        exforSubEntry=7, multiDimFlag=None),
                C4Point(projectile=1, target=40092, targetMetastableState=None, MF=3, MT=1, productMetastableState=None,
                        status='A', cmFlag=None, energy=637000.0, dEnergy=None, data=8.480000, dData=0.110000,
                        cosMuOrLegendreOrder=None, dCosMuOrLegendreOrder=None, eLevelOrHalflife=None,
                        dELevelOrHalflife=None, idOf78=None, reference='M.V.PASECHNIK,ET.AL. (80)', exforEntry='40617',
                        exforSubEntry=7, multiDimFlag=None)])

    def test_a(self):
        readC4DataSet("""#DATASET    40617007
#DATE       19850305
#REACTION   40-ZR-92(N,TOT),,SIG
#PROJ       1
#TARG       40092
#MF         3
#MT         1
#C4BEGIN    [    1 40092   3   1 A ]
#DATA       4
# Prj Targ M MF MT PXC  Energy  dEnergy  Data      dData   Cos/LO   dCos/LO   ELV/HL  dELV/HL I78
#---><---->o<-><-->ooo<-------><-------><-------><-------><-------><-------><-------><-------><->
    1 40092   3   1 A  442000.0          12.74000 1.700000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  507000.0          8.790000 0.570000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  572000.0          9.520000 0.200000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  637000.0          8.480000 0.110000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
#/DATA      4
""".split('\n'))


class TestReadC4File(unittest.TestCase):
    def setUp(self):
        self.testData = """
#
#ENTRY      40617
#AUTHOR1    M.V.Pasechnik+
#YEAR       1980
#INSTITUTE  (4CCPIJI)
#TITLE      TOTAL NEUTRON CROSS-SECTIONS FOR MOLYBDENUM
#+          AND ZYRCONIUM AT LOW ENERGIES
#AUTHOR(S)  M.V.Pasechnik, M.B.Fedorov, V.D.Ovdienko,
#+          G.A.Smetanin, T.I.Jakovenko
#REF-CODE   (C,80KIEV,1,304,8009)
#REFERENCE  Conf. 5.All Union Conf.on Neutron Phys.,Kiev,15-19 Sep 1980
#+          Vol.1, p.304, 1980
#DATASETS   7
#
#DATASET    40617007
#DATE       19850305
#REACTION   40-ZR-92(N,TOT),,SIG
#PROJ       1
#TARG       40092
#MF         3
#MT         1
#C4BEGIN    [    1 40092   3   1 A ]
#DATA       4
# Prj Targ M MF MT PXC  Energy  dEnergy  Data      dData   Cos/LO   dCos/LO   ELV/HL  dELV/HL I78
#---><---->o<-><-->ooo<-------><-------><-------><-------><-------><-------><-------><-------><->
    1 40092   3   1 A  442000.0          12.74000 1.700000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  507000.0          8.790000 0.570000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  572000.0          9.520000 0.200000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  637000.0          8.480000 0.110000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
#/DATA      4
#/DATASET
#/ENTRY
#
#
#ENTRY      10225
#AUTHOR1    L.Green+
#YEAR       1973
#INSTITUTE  (1USABET)
#TITLE      Total cross section measurements with a 252Cf time-of-
#+          flight spectrometer.
#AUTHOR(S)  L.Green, J.A.Mitchell
#REF-CODE   (R,WAPD-TM-1073,197304)
#REFERENCE  Rept. Westinghouse Atomic Power Div.(Bettis) Reports
#+          No.1073, 1973
#DATASETS   27
#
#DATASET    10225020
#DATE       20010305
#REACTION   40-ZR-92(N,TOT),,SIG
#PROJ       1
#TARG       40092
#MF         3
#MT         1
#C4BEGIN    [    1 40092   3   1 A ]
#DATA       3
# Prj Targ M MF MT PXC  Energy  dEnergy  Data      dData   Cos/LO   dCos/LO   ELV/HL  dELV/HL I78
#---><---->o<-><-->ooo<-------><-------><-------><-------><-------><-------><-------><-------><->
    1 40092   3   1 A  504400.0 4961.493 8.146000 0.404200                                       L.GREEN,ET.AL. (73)      10225 20
    1 40092   3   1 A  506500.0 4992.510 8.027000 0.557100                                       L.GREEN,ET.AL. (73)      10225 20
    1 40092   3   1 A  508600.0 5023.591 7.656000 0.276500                                       L.GREEN,ET.AL. (73)      10225 20
#/DATA      3
#/DATASET
#/ENTRY
#
"""
        self.a = [
            readC4Entry("""
#ENTRY      40617
#AUTHOR1    M.V.Pasechnik+
#YEAR       1980
#INSTITUTE  (4CCPIJI)
#TITLE      TOTAL NEUTRON CROSS-SECTIONS FOR MOLYBDENUM
#+          AND ZYRCONIUM AT LOW ENERGIES
#AUTHOR(S)  M.V.Pasechnik, M.B.Fedorov, V.D.Ovdienko,
#+          G.A.Smetanin, T.I.Jakovenko
#REF-CODE   (C,80KIEV,1,304,8009)
#REFERENCE  Conf. 5.All Union Conf.on Neutron Phys.,Kiev,15-19 Sep 1980
#+          Vol.1, p.304, 1980
#DATASETS   7
#
#DATASET    40617007
#DATE       19850305
#REACTION   40-ZR-92(N,TOT),,SIG
#PROJ       1
#TARG       40092
#MF         3
#MT         1
#C4BEGIN    [    1 40092   3   1 A ]
#DATA       4
# Prj Targ M MF MT PXC  Energy  dEnergy  Data      dData   Cos/LO   dCos/LO   ELV/HL  dELV/HL I78
#---><---->o<-><-->ooo<-------><-------><-------><-------><-------><-------><-------><-------><->
    1 40092   3   1 A  442000.0          12.74000 1.700000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  507000.0          8.790000 0.570000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  572000.0          9.520000 0.200000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  637000.0          8.480000 0.110000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
#/DATA      4
#/DATASET
#/ENTRY
            """.split('\n')),
            readC4Entry("""
#ENTRY      10225
#AUTHOR1    L.Green+
#YEAR       1973
#INSTITUTE  (1USABET)
#TITLE      Total cross section measurements with a 252Cf time-of-
#+          flight spectrometer.
#AUTHOR(S)  L.Green, J.A.Mitchell
#REF-CODE   (R,WAPD-TM-1073,197304)
#REFERENCE  Rept. Westinghouse Atomic Power Div.(Bettis) Reports
#+          No.1073, 1973
#DATASETS   27
#
#DATASET    10225020
#DATE       20010305
#REACTION   40-ZR-92(N,TOT),,SIG
#PROJ       1
#TARG       40092
#MF         3
#MT         1
#C4BEGIN    [    1 40092   3   1 A ]
#DATA       3
# Prj Targ M MF MT PXC  Energy  dEnergy  Data      dData   Cos/LO   dCos/LO   ELV/HL  dELV/HL I78
#---><---->o<-><-->ooo<-------><-------><-------><-------><-------><-------><-------><-------><->
    1 40092   3   1 A  504400.0 4961.493 8.146000 0.404200                                       L.GREEN,ET.AL. (73)      10225 20
    1 40092   3   1 A  506500.0 4992.510 8.027000 0.557100                                       L.GREEN,ET.AL. (73)      10225 20
    1 40092   3   1 A  508600.0 5023.591 7.656000 0.276500                                       L.GREEN,ET.AL. (73)      10225 20
#/DATA      4
#/DATASET
#/ENTRY
            """.split('\n'))]
        self.b = list(map(readC4Point, """    1 40092   3   1 A  442000.0          12.74000 1.700000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  507000.0          8.790000 0.570000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  572000.0          9.520000 0.200000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  637000.0          8.480000 0.110000                                       M.V.PASECHNIK,ET.AL. (80)40617  7
    1 40092   3   1 A  504400.0 4961.493 8.146000 0.404200                                       L.GREEN,ET.AL. (73)      10225 20
    1 40092   3   1 A  506500.0 4992.510 8.027000 0.557100                                       L.GREEN,ET.AL. (73)      10225 20
    1 40092   3   1 A  508600.0 5023.591 7.656000 0.276500                                       L.GREEN,ET.AL. (73)      10225 20""".split(
            '\n')))

    def test_a(self):
        self.assertEqual(readC4File(self.testData.split('\n'), asPointList=False), self.a)

    def test_b(self):
        self.assertEqual(readC4File(self.testData.split('\n'), asPointList=True), self.b)

# ------------------------------------------------
# Main !!
# ------------------------------------------------

if __name__ == "__main__":

#    try:
#        import xmlrunner

#        unittest.main(testRunner=xmlrunner.XMLTestRunner(output='test-results'))
#    except ImportError:
        unittest.main()
        print()
