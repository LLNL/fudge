# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import unittest
import os
import sys
import numpy

from pqu import PQU
from fudge import reactionSuite
from fudge.processing.resonances import reconstructResonances, makeUnresolvedProbabilityTables


# ----------------------------------------------------------------------------------
#
#   Load the data files
#
# ----------------------------------------------------------------------------------

TEST_DATA_PATH, this_filename = os.path.split(os.path.realpath(__file__))
VERBOSE = False
DOFETESTS = True

if not os.path.exists(os.path.join(TEST_DATA_PATH, 'MLBWExampleZr90_testFile.gnds.xml')):
    print('WARNING: you need to run rebuild_test_data.py before running the test.')
    sys.exit(0)

try:
    import brownies.BNL.restools
    HAVEBRESTOOLS = True
except ImportError:
    HAVEBRESTOOLS = False


def floatsAlmostEqual(float1, float2, epsilon=1e-15):
    diff = float1 - float2
    if diff == 0.0:
        return True
    bigger = max(abs(float1), abs(float2))
    return diff <= bigger * epsilon


if VERBOSE:
    print('reading MLBWExampleZr90...')
MLBWExampleZr90 = reactionSuite.readXML(TEST_DATA_PATH+os.sep+'MLBWExampleZr90_testFile.gnds.xml')


class TestWithIsClose(unittest.TestCase):

    def assertWithinXPercent(self, a, b, percent=1.0, absTol=1e-14):
        self.assertTrue(numpy.isclose(a, b, rtol=percent / 100., atol=absTol), 'a and b differ by %s' % str(a - b))

    def assertAllWithinXPercent(self, aL, bL, percent=1.0, absTol=1e-14):
        if len(aL) != len(bL):
            return False
        self.assertTrue(all([numpy.isclose(*x, rtol=percent / 100., atol=absTol) for x in zip(aL, bL)]),
                        'Items in list not equal')

    def assertPQUEquivalent(self, a, b):
        if not isinstance(a, PQU.PQU):
            raise TypeError("First argument not a PhysicalQuantity")
        if not isinstance(b, PQU.PQU):
            raise TypeError("Second argument not a PhysicalQuantity")
        self.assertTrue(a.equivalent(b), 'Arguments not within 1-sigma uncertainties of one another')

    def test_asserts(self):
        self.assertPQUEquivalent(PQU.PQU(value=1.234e3, unit='eV', uncertainty=0.1e3),
                                 PQU.PQU(value=1.324, unit='keV', uncertainty=0.1))
        self.assertWithinXPercent(1.234567, 1.234568, percent=1.0, absTol=1e-14)


# ----------------------------------------------------------------------------------
#
#   The test cases
#
# ----------------------------------------------------------------------------------





class TestURRClassAndBaseClasses( TestWithIsClose ):

    def setUp( self ):
        self.Zr90URR = reconstructResonances.URRcrossSection( MLBWExampleZr90.resonances.unresolved.evaluated, verbose=False )
        self.Zr90_table_generator = makeUnresolvedProbabilityTables.ProbabilityTableGenerator( MLBWExampleZr90 )

    @unittest.skip("getLastResolvedResonanceRegion returns a copy. With no __eq__ defined, self.assertEqual returns False")
    def test_getLastResolvedResonanceRegion(self):
        """
        TODO: Should have multiple region test too
        """
        self.assertEqual( self.Zr90_table_generator.getLastResolvedResonanceRegion(), self.Zr90_table_generator.reactionSuite.resonances.resolved.evaluated )

    def test_getLastResolvedResonanceEnergy(self):
        answers={(0,0.5):198400.0, (1, 1.5):193400.0, (1, 0.5):189100.0, (2, 1.5):188400.0, (2, 2.5):160000.0}
        for lj in answers:
            self.assertEqual( self.Zr90_table_generator.getLastResolvedResonanceEnergy(*lj), answers[lj] )

    # def test_URRPDFTable(self):
    #     u = URRPDFTable()
    #     u.eBins = numpy.linspace(start=0.,stop=4.,num=2+1)
    #     u.xsBins = numpy.linspace(start=1.,stop=1.5,num=3+1)
    #     u['capture'] = numpy.zeros(shape=(2, 3))
    #     u['capture'][0, 0] = 1.
    #     u['capture'][1, 2] = 2.
    #     u.save('__junk__.txt')
    #     v = URRPDFTable()
    #     v.load('__junk__.txt')
    #     os.remove('__junk__.txt')
    #     self.assertEqual(str(u),str(v))

    # Currently broken due to refactoring in makeUnresolvedProbabilityTables
    @unittest.skipIf(not HAVEBRESTOOLS, "blurr not imported")
    def test_sampleRR(self):
        self.Zr90URR.getWidthsAndSpacings()
        self.Zr90_table_generator.extrapolate_URR_parameters(self.Zr90URR.lowerBound, self.Zr90URR.upperBound)
        lastResEnergies={(0,0.5):198400.0, (1, 1.5):193400.0, (1, 0.5):189100.0, (2, 1.5):188400.0, (2, 2.5):160000.0}
        fakeRR=self.Zr90_table_generator.sampleRR(lastResEnergies, lowerBound=None, upperBound=None, style='goe', verbose=False)

        # Make sure all spin groups present
        for lj in lastResEnergies:
            self.assertTrue(lj in fakeRR)

        # Make sure starting energies are kosher,
        # Just above the lowerBound of the URR and the upper end of the RRR, by 1-2 units of mean level spacing
        for lj in lastResEnergies:
            firstEnergy = fakeRR[lj].data[0][0]
            self.assertLess(
                (firstEnergy-max(lastResEnergies[lj], self.Zr90URR.lowerBound))*self.Zr90URR.levelDensities[lj].evaluate(self.Zr90URR.lowerBound), 3)

        # For next two tests, just work with one channel
        lj0 = list( fakeRR.keys() )[0]
        num_levels=len(fakeRR[lj0].data)

        # Check the mean level spacing, this is very non-trivial test, we check that the spacing-spacing
        # correlation is consistent with GOE expectations
        diff_spacings=[]
        spacings=[fakeRR[lj0].data[i+1][0]-fakeRR[lj0].data[i][0] for i in range(num_levels-1)]
        for i in range(num_levels-1):
            aveE=0.5*(fakeRR[lj0].data[i+1][0]+fakeRR[lj0].data[i][0])
            diff_spacings.append(spacings[i]-1.0/self.Zr90URR.levelDensities[lj0].evaluate(aveE))
        correlation_matrix = numpy.corrcoef(numpy.array([[diff_spacings[i], diff_spacings[i+1]] for i in range(num_levels-2)]).T)
# BRB
        print( 'Test 3', correlation_matrix[0,1], -0.27, 30.0 )
        self.assertWithinXPercent(correlation_matrix[0,1], -0.27, 30.0) # GOE says this should be -0.27, this is a Monte-Carlo test, so have big potential variation

        # Check the mean width against the average width I get by averaging the average width function over the entire URR
        # Only check 'capture' because 'elastic' in URR is "reduced" and I'm too lazy to code up the correction to full width
        widths={'capture':[], 'elastic':[]}
        colNames=[c.name for c in fakeRR[lj0].columns]
        for w in widths.keys():
            i=colNames.index(w)
            for r in fakeRR[lj0].data:
                widths[w].append(r[i])
            aveWidthFunc=self.Zr90URR.averageWidths[lj0][w]
            aveWidth=(aveWidthFunc.integrate()/(aveWidthFunc.domainMax-aveWidthFunc.domainMin)).getValue()
            if w == 'capture':
                self.assertWithinXPercent(numpy.mean(widths[w]), aveWidth, 20.0)
            #else:
                #print(w, numpy.mean(widths[w]))

    @unittest.skipIf(True, "Still in development") #not HAVEBRESTOOLS, "blurr not imported")
    def test_getURRPDF(self):
        thePDFs = self.Zr90_table_generator.generatePDFs(reconstructResonances.SLBWcrossSection, nSamples=2, verbose=False)
        #print(thePDFs)
        #egrid=thePDFs.eBins
        #for k in thePDFs:
            #print(k, thePDFs[k][1].toList())
            #thePDFs.pdf_at_energy(k, egrid[2]).plot()


# ----------------------------------------------------------------------------------
#
#   main!!!
#
# ----------------------------------------------------------------------------------

if __name__=="__main__":

    unittest.main()
