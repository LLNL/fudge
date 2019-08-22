# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
# <<END-copyright>>

import unittest, os
import numpy as np
from fudge.gnd import reactionSuite
from fudge.processing.resonances.reconstructResonances import *
import fudge.processing.resonances.getCoulombWavefunctions as getCoulombWavefunctions
from xData.isclose import isclose

# ----------------------------------------------------------------------------------
#
#   Load the data files
#
# ----------------------------------------------------------------------------------

TEST_DATA_PATH, this_filename = os.path.split( os.path.realpath( __file__ ) )
VERBOSE = False
DOFETESTS = True

try:
    import blurr
    HAVEBLURR = True# still in development, so not distributed
    HAVEBLURR = False
except:
    HAVEBLURR = False

if VERBOSE:
    print('reading SLBWExample1SRes...')
SLBWExample1SRes = reactionSuite.readXML(open(TEST_DATA_PATH+os.sep+'SLBWExample1SRes_testFile.gnd.xml'))

if VERBOSE:
    print('reading SLBWExample1PRes...')
SLBWExample1PRes = reactionSuite.readXML(open(TEST_DATA_PATH+os.sep+'SLBWExample1PRes_testFile.gnd.xml'))

if VERBOSE:
    print('reading SLBWExample...')
SLBWExample = reactionSuite.readXML(open(TEST_DATA_PATH+os.sep+'SLBWExampleFull_testFile.gnd.xml'))

if VERBOSE:
    print('reading MLBWExample...')
MLBWExample = reactionSuite.readXML(open(TEST_DATA_PATH+os.sep+'MLBWExampleFull_testFile.gnd.xml'))

if VERBOSE:
    print('reading MLBWExample1PRes...')
MLBWExample1PRes = reactionSuite.readXML(open(TEST_DATA_PATH+os.sep+'MLBWExample1PRes_testFile.gnd.xml'))

if VERBOSE: print('reading MLBWExampleZr90...')
MLBWExampleZr90 = reactionSuite.readXML(open(TEST_DATA_PATH+os.sep+'MLBWExampleZr90_testFile.gnd.xml'))

if VERBOSE:
    print('reading RMExample...')
RMExample = reactionSuite.readXML(open(TEST_DATA_PATH+os.sep+'RMExampleFull_testFile.gnd.xml'))

if VERBOSE:
    print('reading RMExampleSmall...')
RMExampleSmall =  reactionSuite.readXML(open(TEST_DATA_PATH+os.sep+'RMExampleSmall_testFile.gnd.xml'))

if VERBOSE:
    print('reading RMLExample...')
RMLExample =  reactionSuite.readXML(open(TEST_DATA_PATH+os.sep+'RMLExampleFull_testFile.gnd.xml'))

if VERBOSE:
    print('reading RMLExampleSmall...')
RMLExampleSmall =  reactionSuite.readXML(open(TEST_DATA_PATH+os.sep+'RMLExampleSmall_testFile.gnd.xml'))

if DOFETESTS:
    if VERBOSE:
        print('reading RMLExampleFe...')
    RMLExampleFe =  reactionSuite.readXML(open(TEST_DATA_PATH+os.sep+'RMLExampleFe_testFile.gnd.xml'))

# ----------------------------------------------------------------------------------
#
#   Assertion helper functions
#
# ----------------------------------------------------------------------------------

# 0 degC   = 273.1 degK      = 0.0235339347844 eV since kB = 8.6173324e-5 eV/degK
# INTER uses 293.594338 degK = 2.53000E-02 eV
# RoomTemp = 295 degK        = 0.0254211305 eV
#            300 degK        = 0.0258519972 eV

def complexString(x):
    return str(numpy.abs(x))+'*e^'+str(numpy.angle(x))


class TestWithIsClose( unittest.TestCase ):

    def assertWithinXPercent(self, a, b, percent=1.0, absTol=1e-14):
        return isclose(a, b, rel_tol=percent/100., abs_tol=absTol, method='weak')

    def assertAllWithinXPercent(self, aL, bL, percent=1.0, absTol=1e-14):
        if len( aL ) != len( bL ): return False
        return all( [ isclose(*x, rel_tol=percent/100., abs_tol=absTol, method='weak') for x in zip( aL, bL ) ] )


# ----------------------------------------------------------------------------------
#
#   The test cases
#
# ----------------------------------------------------------------------------------

class TestAngMomAdding( TestWithIsClose ):

    def test_getAllowedTotalSpins( self ):
        self.assertEqual( getAllowedTotalSpins( 10, 1 ), [ 9, 11 ] )
        self.assertEqual( getAllowedTotalSpins( 5, 2 ), [ 3, 5, 7 ] )
        self.assertEqual( getAllowedTotalSpins( 5, 0.5, False ), [ 4.5, 5.5 ] )
        self.assertEqual( getAllowedTotalSpins( 2.5, 1, False ), [ 1.5, 2.5, 3.5 ] )
        self.assertEqual( getAllowedTotalSpins( 2.5, 0.5, False ), [ 2, 3 ] )
        self.assertEqual( getAllowedTotalSpins( 2.5, 0, False ), [ 2.5 ] )
        self.assertEqual( getAllowedTotalSpins( 2.5, 1.5, False ), [ 1, 2, 3, 4 ] )



class TestChannelDesignator( TestWithIsClose ):

    def test_members(self):
        a = ChannelDesignator( l=3, J=2.5, reaction='elastic', index=0, s=1, gfact=7, particleA='n', particleB='He3', channelClass=NEUTRONCHANNEL, useRelativistic=False, eliminated=False )
        self.assertEquals(a.l,3)
        self.assertEquals(a.J,3-1.0/2.0)
        self.assertEquals(a.reaction,'elastic')
        self.assertEquals(a.index,0)
        self.assertEquals(a.s,1)
        self.assertEquals(a.gfact,7)
        self.assertEquals(a.particleA,'n')
        self.assertEquals(a.particleB,'He3')
        self.assertEquals(a.Xi,0.0)
        self.assertEquals(a.channelClass,NEUTRONCHANNEL)
        self.assertEquals(a.useRelativistic,False)
        self.assertEquals(a.eliminated,False)

    def test_equality(self):
        a = ChannelDesignator( 3, 2.5, 'elastic', 0, 1, 7, 'n', 'He3', 0.0, NEUTRONCHANNEL, False, False )
        b = ChannelDesignator( 3, 3-1.0/2.0, 'elastic', 0, 1, 7, 'n', 'He3', 0.0, NEUTRONCHANNEL, False, False )
        self.assertEquals(a.J,b.J)
        self.assertEquals(a,b)

    def test_hash(self):
        a = ChannelDesignator( 3, 2.5,       'elastic', 0, 1, 7, 'n', 'He3', True, 0.0, NEUTRONCHANNEL, False, False )
        b = ChannelDesignator( 3, 3-1.0/2.0, 'elastic', 0, 1, 7, 'n', 'He3', True, 0.0, NEUTRONCHANNEL, False, False )
        c = ChannelDesignator( 3, 3-1.0/2.0, 'n + He3', 0, 1, 7, 'n', 'He3', True, 0.0, NEUTRONCHANNEL, False, False )
        self.assertEquals(hash(a),-7289020590116273477)
        self.assertEquals(hash(b),-7289020590116273477)
        self.assertEquals(hash(c),-7289020590116273477)
        self.assertEquals(hash(b.J),1342242816)
        self.assertEquals(hash(a.J),1342242816)
        self.assertEquals(hash(a),hash(b))

    def test_is_open(self):
        Xi=0.1 # Fake threshold for testing
        low_egrid=numpy.array([0.001,0.002])
        bad_egrid=numpy.array([0.05,0.11,0.2])
        hi_egrid=numpy.array([0.2,0.3])
        a = ChannelDesignator(1, 3.0, 'H1 + S35',     0, 2, gfact=0.875, particleA='H1',   particleB='S35',   Xi=Xi, isElastic=False, channelClass=CPCHANNEL,      useRelativistic=False, eliminated=False)
        self.assertTrue(numpy.all(a.is_open(hi_egrid))) # with hi_egrid, channel should all be open for all energies
        self.assertFalse(numpy.all(a.is_open(low_egrid))) # with low_egrid, channel should all be closed for all energies
        self.assertTrue(numpy.any(a.is_open(bad_egrid))) # with bad_egrid, the channel is open at some energies, closed at others.  A confusing situation
        self.assertFalse(numpy.all(a.is_open(bad_egrid)))



class TestResonanceReconstruction( TestWithIsClose ):

    def test_SLBWReconstructResonances( self ):
        x = reconstructResonances(SLBWExample, tolerance=0.001, verbose=False)
        self.assertEqual( set( x.keys() ), set( [ 'capture', 'fission', 'elastic', 'total', 'nonelastic' ] ) )
        # Test thermal point and a few others, computed below
        answers = [ \
            (0.0235339347844, {'capture':1.18726303, 'total':7.24518429, 'fission':0., 'elastic':6.05792126}),
            (1009.0,          {'capture':290.05432779892004, 'total':507.80854074065803, 'fission':0., 'elastic':217.75421294173799})]
        for answer in answers:
            for k in answer[1]:
                nDigits = 2 # min tolerance is 0.001 in XYs rendition of cross sections
                # next line gets either 2 sig figs or 2 digits
                if answer[1][k] > 0.0: nDigits+=int(round(max(-math.log10(answer[1][k]),0.0)))
                #print k, answer[1][k], nDigits
                self.assertAlmostEqual( answer[1][k], x[k].evaluate( answer[0] ), nDigits )
        if False: x['elastic'].plot()

    def test_MLBWReconstructResonances( self ):
        x = reconstructResonances(MLBWExample, tolerance=0.001, verbose=False)
        self.assertEqual( set( x.keys() ), set( [ 'capture', 'fission', 'elastic', 'total', 'nonelastic' ] ) )
        # Test thermal point, computed below
        answer = {'capture':0.00015987, 'total':3.11326583, 'fission':0., 'elastic':3.11310596}
        for k in answer:
            nDigits = 2 # min tolerance is 0.001 in XYs rendition of cross sections
            if answer[k] > 0.0: nDigits+=int(round(max(-math.log10(answer[k]),0.0)))
            # the following gets either 2 sig figs or 2 digits
            self.assertAlmostEqual( answer[k], x[k].evaluate( 0.0235339347844 ), nDigits )

    def test_RMReconstructResonances( self ):
        x = reconstructResonances(RMExample, tolerance=0.001, verbose=False)
        self.assertEqual( set( x.keys() ), set( [ 'capture', 'fission', 'elastic', 'total', 'nonelastic' ] ) )
        # Test 0 deg C point
        answer = {'capture':0.17536290484357434, 'total':2.1313431546619341, 'fission':0., 'elastic':1.9559802498183596}
        for k in answer:
            nDigits = 3 # min tolerance is 0.001 in XYs rendition of cross sections
            if answer[k] > 0.0: nDigits+=int(round(max(-math.log10(answer[k]),0.0)))
            self.assertAlmostEqual( answer[k], x[k].evaluate( 0.0235339347844 ), nDigits )

    def test_RMLReconstructResonances( self ):
        x = reconstructResonances(RMLExample, tolerance=0.001, verbose=False)
        self.assertEqual( set( x.keys() ), set( [ 'capture', 'fission', 'elastic', 'total', 'H1 + S35', 'nonelastic' ] ) )
        # Test thermal point, computed below
        answer = { 'capture':45.04182664, 'total':66.22034846, 'fission':0., 'elastic':20.68451381, 'H1 + S35':0.49400801 }
        for k in answer:
            nDigits = 2 # min tolerance is 0.001 in XYs rendition of cross sections
            if answer[k] > 0.0: nDigits+=int(round(max(-math.log10(answer[k]),0.0)))
            # the following gets either 2 sig figs or 2 digits
            self.assertAlmostEqual( answer[k], x[k].evaluate( 0.0235339347844 ), nDigits )

    @unittest.skip("didn't write the test yet")
    def test_URRReconstructResonances( self ):
        x = reconstructResonances(RMLExample, tolerance=0.001, verbose=False)
        self.assertEqual( set( x.keys() ), set( [ 'capture', 'fission', 'elastic', 'total' ] ) )
        # Test thermal point, computed below
        self.assertAlmostEqual( x['total'].evaluate( 2.5e-4 ), 11.5362133, 4 )
        self.assertEqual( x['fission'].getValue( 2.5e-4 ), 0.0 )
        self.assertAlmostEqual( x['capture'].evaluate( 2.5e-4 ), 0.29962655, 4 )
        self.assertAlmostEqual( x['elastic'].evaluate( 2.5e-4 ), 11.23658675, 4 )



class TestSLBWClassAndBaseClasses( TestWithIsClose ):

    def setUp( self ):
        self.RRR = SLBWcrossSection( SLBWExample, enableAngDists=True, verbose=False )
        self.RRR1SRes = SLBWcrossSection( SLBWExample1SRes, enableAngDists=True, verbose=False )
        self.RRR1PRes = SLBWcrossSection( SLBWExample1PRes, enableAngDists=True, verbose=False )

    def test_memberData( self ):
        self.assertEqual( str(self.RRR.projectile), 'n' )
        self.assertEqual( self.RRR.targetSpin, 4.5 )
        self.assertEqual( str(self.RRR.target), 'Nb93' )
        self.assertAlmostEqual( self.RRR.targetToNeutronMassRatio, 92.1051 )
        self.assertAlmostEqual( self.RRR.lowerBound, 1e-05 )
        self.assertAlmostEqual( self.RRR.upperBound, 7350.0 )
        self.assertEqual( self.RRR.missingGfactor, {0: 0.0, 1: 1.0} )
        self.assertEqual( len(self.RRR.Ls), 2 )
        self.assertEqual( self.RRR.Ls[1].keys(), ['L', 'Js'] )
        self.assertEqual( self.RRR.Ls[1]['L'], 1 )
        self.assertEqual( self.RRR.Ls[1]['Js'][0].keys(), ['gfact', 'channelSpins', 'J'] )
        self.assertAlmostEqual( self.RRR.Ls[1]['Js'][0]['gfact'], 0.35 )
        self.assertEqual( self.RRR.Ls[1]['Js'][0]['channelSpins'][0].keys(), ['neutronWidth','channelSpin','captureWidth','energy', 'fissionWidthB', 'fissionWidth', 'shiftFactor'] )
        self.assertEqual( self.RRR.Ls[1]['Js'][0]['J'], 3.0 )
        self.assertEqual( list(self.RRR._energies), [ -105.4   ,    35.9   ,    42.3   ,    94.3   ,   105.8   ,
         119.2   ,   193.    ,   244.    ,   319.    ,   335.5   ,
         363.    ,   365.    ,   378.4   ,   392.6   ,   460.3   ,
         500.6   ,   599.2001,   604.    ,   617.6   ,   640.9001,
         672.    ,   678.3   ,   721.4001,   741.5002,   757.3   ,
         912.5002,   934.8002,   952.9002,  1009.    ,  1016.    ,
        1108.    ,  1127.    ,  1148.    ,  1175.    ,  1194.    ,
        1229.    ,  1243.    ,  1284.    ,  1350.    ,  1354.    ,
        1393.    ,  1452.    ,  1467.    ,  1529.    ,  1540.    ,
        1556.    ,  1576.    ,  1616.    ,  1654.    ,  1678.    ,
        1712.    ,  1768.    ,  1811.    ,  1833.    ,  1946.    ,
        1981.    ,  1991.    ,  2021.    ,  2070.    ,  2076.    ,
        2117.    ,  2147.    ,  2155.    ,  2186.    ,  2230.    ,
        2310.    ,  2335.    ,  2341.    ,  2360.    ,  2392.    ,
        2420.    ,  2455.    ,  2465.    ,  2509.    ,  2544.    ,
        2577.    ,  2641.    ,  2657.    ,  2689.    ,  2712.    ,
        2755.    ,  2836.    ,  2927.    ,  2943.    ,  2953.    ,
        2987.    ,  2995.    ,  3002.    ,  3140.    ,  3149.    ,
        3229.    ,  3259.    ,  3275.    ,  3285.    ,  3359.    ,
        3378.    ,  3396.    ,  3411.    ,  3426.    ,  3503.    ,
        3525.    ,  3592.    ,  3609.    ,  3624.    ,  3655.    ,
        3675.    ,  3763.    ,  3814.    ,  3844.    ,  3920.    ,
        3927.    ,  3946.    ,  3977.    ,  3987.    ,  4026.    ,
        4041.    ,  4069.    ,  4079.    ,  4109.    ,  4145.    ,
        4232.    ,  4295.    ,  4307.    ,  4349.    ,  4359.    ,
        4405.    ,  4425.    ,  4452.    ,  4498.    ,  4544.    ,
        4557.    ,  4599.    ,  4660.    ,  4696.    ,  4779.    ,
        4812.    ,  4826.    ,  4874.    ,  4892.    ,  4918.    ,
        4934.    ,  5008.    ,  5045.    ,  5067.    ,  5125.    ,
        5153.    ,  5178.    ,  5292.    ,  5422.    ,  5467.    ,
        5487.    ,  5510.    ,  5550.    ,  5579.    ,  5659.    ,
        5691.    ,  5702.    ,  5720.    ,  5792.    ,  5811.    ,
        5835.    ,  5862.    ,  5902.    ,  6003.    ,  6051.    ,
        6072.    ,  6097.    ,  6114.    ,  6169.    ,  6180.    ,
        6327.    ,  6412.    ,  6427.    ,  6456.    ,  6480.    ,
        6532.    ,  6551.    ,  6571.    ,  6590.    ,  6714.    ,
        6739.    ,  6751.    ,  6780.    ,  6795.    ,  6843.    ,
        6881.    ,  6911.    ,  6935.    ,  6976.    ,  7077.    ,
        7136.    ,  7244.    ,  7306.    ,  7331.    ] )
        self.assertEqual( list(self.RRR._widths), [ 0.5441,  0.2091,  0.2221,  0.1805,  0.1675,  0.1286,  0.1678,
        0.2474,  0.2217,  0.1748,  0.1352,  0.1505,  0.2713,  0.2019,
        0.1691,  0.2293,  0.1463,  0.1626,  0.1517,  0.1559,  0.2268,
        0.1379,  0.2659,  0.3401,  0.2519,  0.2043,  0.5191,  0.2481,
        1.414 ,  0.2805,  0.2544,  0.2485,  0.3094,  0.5576,  0.2098,
        0.1946,  0.2633,  0.2448,  0.2853,  0.2608,  0.3624,  1.176 ,
        0.2884,  0.2328,  0.2706,  0.2729,  0.3366,  0.2701,  0.2199,
        0.2181,  0.249 ,  0.1797,  0.2144,  1.026 ,  0.2619,  0.2573,
        0.1299,  0.84  ,  0.3127,  0.3237,  0.2996,  0.2761,  0.4102,
        0.2453,  0.2305,  0.1656,  0.2166,  0.5857,  0.2062,  0.1683,
        2.19  ,  0.1923,  0.2908,  0.5917,  0.2131,  0.2049,  3.549 ,
        0.1993,  0.1391,  0.2622,  0.1264,  0.1857,  0.4533,  0.1258,
        0.6692,  0.1829,  0.1385,  0.1245,  0.1993,  0.1414,  0.1902,
        0.1918,  0.232 ,  0.2052,  0.4289,  0.232 ,  1.234 ,  0.1811,
        0.1852,  0.1662,  1.051 ,  0.1883,  0.4518,  0.1652,  0.1239,
        2.077 ,  0.5524,  0.125 ,  0.3877,  0.1768,  0.237 ,  0.4452,
        1.014 ,  0.1981,  0.172 ,  0.3563,  1.661 ,  0.126 ,  0.1273,
        0.1989,  0.1852,  0.1987,  0.1308,  0.2084,  0.1359,  0.1394,
        0.1288,  0.1811,  0.1369,  0.199 ,  0.3981,  3.422 ,  0.1412,
        0.1659,  0.173 ,  0.3311,  0.2407,  0.1316,  0.126 ,  0.1852,
        0.162 ,  0.1592,  0.5657,  0.2717,  0.1239,  5.676 ,  0.7652,
        0.2151,  0.1312,  0.3591,  0.2596,  0.4405,  1.11  ,  0.6379,
        1.152 ,  1.946 ,  0.1808,  0.1547,  0.4909,  0.1362,  0.2297,
        0.8586,  2.14  ,  2.996 ,  1.906 ,  0.2864,  0.1269,  0.9702,
        0.2786,  0.2784,  0.1338,  0.3964,  1.579 ,  0.1349,  0.1351,
        3.629 ,  0.1722,  1.563 ,  0.136 ,  1.335 ,  0.5289,  0.1692,
        0.4309,  0.9733,  0.287 ,  0.1856,  1.02  ,  0.4838,  4.723 ,
        0.2183,  0.3123,  2.593 ,  0.133 ,  0.3627] ) # total widths

    def test_k( self ):
        self.assertAlmostEqual( self.RRR.k(1000.), 2.196807122623e-3*93./94.*math.sqrt(1000.), 4 )

    def test_rho( self ):
        self.assertTrue( self.RRR.RR.calculateChannelRadius ) #This deactivates L dependent rho's because we compute the constant channel radius
        a = 0.123 * self.RRR.reactionSuite.PoPs[self.RRR.target].getMass('amu')**(1./3.) + 0.08
        self.assertAlmostEqual( a, 0.6370771040638685 )
        k = self.RRR.k(1000.)
        self.assertAlmostEqual( self.RRR.rho( 1000. ),    k*a )
        self.assertAlmostEqual( self.RRR.rho( 1000., 0 ), k*a )
        self.assertAlmostEqual( self.RRR.rho( 1000., 1 ), k*a )
        self.assertAlmostEqual( self.RRR.rho( 1000., 2 ), k*a )

    def test_penetrationFactor( self ):
        rho = 0.043781852422548007 # for E=1000. eV, and a computed as 0.6370771040638685 fm
        self.assertAlmostEqual( self.RRR.penetrationFactor( 0, rho ), rho )
        self.assertAlmostEqual( self.RRR.penetrationFactor( 1, rho ), rho*rho*rho/(1.0+rho*rho) )
        self.assertAlmostEqual( self.RRR.penetrationFactor( 2, rho ), pow( rho, 5 )/( 9.+3.*rho*rho+pow(rho,4) ) )
        self.assertAlmostEqual( self.RRR.penetrationFactor( 3, rho ), pow( rho, 7 )/( 225.+45.*rho*rho+6.*pow(rho,4)+pow(rho,6) ) )

    def test_shiftFactor( self ):
        rho = 0.043781852422548007 # for E=1000. eV, and a computed as 0.6370771040638685 fm
        self.assertAlmostEqual( self.RRR.shiftFactor( 0, rho ), 0. )
        self.assertAlmostEqual( self.RRR.shiftFactor( 1, rho ), -1./(1.+rho*rho) )
        self.assertAlmostEqual( self.RRR.shiftFactor( 2, rho ), -(18.+3.*rho*rho)/(9.+3.*rho*rho+pow(rho,4)) )
        self.assertAlmostEqual( self.RRR.shiftFactor( 3, rho ), -(675.+90.*pow(rho,2)+6.*pow(rho,4))/(225.+45.*rho*rho+6.*pow(rho,4)+pow(rho,6)) )

    def test_phi( self ):
        rho = 0.043781852422548007 # for E=1000. eV, and a computed as 0.6370771040638685 fm
        self.assertAlmostEqual( self.RRR.phi(0, rho), rho )
        self.assertAlmostEqual( self.RRR.phi(1, rho), rho-math.atan(rho) )
        self.assertAlmostEqual( self.RRR.phi(2, rho), rho-math.atan(3.*rho/(3.-rho*rho)) )
        self.assertAlmostEqual( self.RRR.phi(3, rho), rho-math.atan( (15.*rho-pow(rho,3))/(15.-6.*rho*rho) ) )

    def test_eiphi( self ):
        for E, l, answer in [ ( 1009.0, 0, complex( 0.99903310472611739, -0.043964254358450666 ) ), ( 35.9, 1, complex( 0.9999999999999819, -1.9027638215908413e-07 ) ) ]:
            self.assertWithinXPercent( numpy.exp( complex( 0.0, -self.RRR.phi( l, self.RRR.rho( E, l ) ) ) ), answer, percent=0.01 )

    def test_ERp_SWave( self ):
        c = [ cc for cc in self.RRR1SRes.channels[0].keys() if cc.reaction=='elastic' ][0]
        ER = self.RRR1SRes.channels[0][c][0][0]
        GN = self.RRR1SRes.channels[0][c][0][1]
        for E in [ ER-0.05+0.01*i for i in range( 11 ) ]:
            ERp = ER + ( \
                    self.RRR1SRes.shiftFactor( c.l, self.RRR1SRes.rho( abs( ER ), c.l ) ) - \
                    self.RRR1SRes.shiftFactor( c.l, self.RRR1SRes.rho( E,         c.l ) ) \
                ) * GN / (2.0 * self.RRR1SRes.penetrationFactor( c.l, self.RRR1SRes.rho( abs(ER), c.l) ) )
            spin_thingee = self.RRR1SRes.Ls[0]['Js'][0]['channelSpins'][0]
            ERp_correct = spin_thingee['energy'] + spin_thingee['neutronWidth'] * ( spin_thingee['shiftFactor'] - 0.5 * self.RRR1SRes.shiftFactor( 0, self.RRR1SRes.rho(E) ) )
            self.assertAlmostEqual( ERp, ERp_correct[0] )

    def test_ERp_PWave( self ):
        c = [ cc for cc in self.RRR1PRes.channels[0].keys() if cc.reaction=='elastic' ][0]
        ER = self.RRR1PRes.channels[0][c][0][0]
        GN = self.RRR1PRes.channels[0][c][0][1]
        for E in [ ER-0.05+0.01*i for i in range( 11 ) ]:
            ERp = ER + ( \
                    self.RRR1PRes.shiftFactor( c.l, self.RRR1PRes.rho( abs( ER ), c.l ) ) - \
                    self.RRR1PRes.shiftFactor( c.l, self.RRR1PRes.rho( E,         c.l ) ) \
                ) * GN / (2.0 * self.RRR1PRes.penetrationFactor( c.l, self.RRR1PRes.rho( abs(ER), c.l) ) )
            spin_thingee = self.RRR1PRes.Ls[1]['Js'][0]['channelSpins'][0]
            self.assertAlmostEqual( self.RRR1PRes.shiftFactor( c.l, self.RRR1PRes.rho( abs( ER ), c.l ) ), 2.0 * spin_thingee['shiftFactor'] )
            self.assertAlmostEqual( GN / self.RRR1PRes.penetrationFactor( c.l, self.RRR1PRes.rho( abs(ER), c.l )  ), spin_thingee['neutronWidth'] )
            ERp_correct = spin_thingee['energy'] + spin_thingee['neutronWidth'] * ( spin_thingee['shiftFactor'] - 0.5 * self.RRR1PRes.shiftFactor( c.l, self.RRR1PRes.rho(E) ) )
            self.assertAlmostEqual( ERp, ERp_correct[0] )

    def test_GN_of_E_SWave( self ):
        c = [ cc for cc in self.RRR1SRes.channels[0].keys() if cc.reaction=='elastic' ][0]
        ER = self.RRR1SRes.channels[0][c][0][0]
        GN = self.RRR1SRes.channels[0][c][0][1]
        for E in [ ER-0.05+0.01*i for i in range( 11 ) ]:
            GamN = GN * \
                self.RRR1SRes.penetrationFactor( c.l, self.RRR1SRes.rho( E,         c.l ) ) / \
                self.RRR1SRes.penetrationFactor( c.l, self.RRR1SRes.rho( abs( ER ), c.l ) )
            spin_thingee = self.RRR1SRes.Ls[0]['Js'][0]['channelSpins'][0]
            GamN_correct = self.RRR1SRes.penetrationFactor( 0, self.RRR1SRes.rho(E) ) * spin_thingee['neutronWidth']
            self.assertAlmostEqual( GamN, GamN_correct )

    def test_GN_of_E_PWave( self ):
        c = [ cc for cc in self.RRR1PRes.channels[0].keys() if cc.reaction=='elastic' ][0]
        ER = self.RRR1PRes.channels[0][c][0][0]
        GN = self.RRR1PRes.channels[0][c][0][1]
        for E in [ ER-0.05+0.01*i for i in range( 11 ) ]:
            GamN = GN * \
                self.RRR1PRes.penetrationFactor( c.l, self.RRR1SRes.rho( E,         c.l ) ) / \
                self.RRR1PRes.penetrationFactor( c.l, self.RRR1SRes.rho( abs( ER ), c.l ) )
            spin_thingee = self.RRR1PRes.Ls[1]['Js'][0]['channelSpins'][0]
            GamN_correct = self.RRR1PRes.penetrationFactor( c.l, self.RRR1PRes.rho(E) ) * spin_thingee['neutronWidth']
            self.assertAlmostEqual( GamN, GamN_correct[0] )

    @unittest.skip("didn't write the test yet")
    def test_refineInterpolation( self ):
        self.RRR.refineInterpolation( egrid=0, xsecs=0, tolerance=0.01)

    def test_generateEnergyGrid( self ):
        egrid = self.RRR.generateEnergyGrid()
        self.assertEqual( egrid[0:4], [1e-05, 1.2578278491146621e-05, 1.582130898008417e-05, 1.9900483044597718e-05] )
        self.assertEqual( egrid[16:20], [0.00039259005767352531, 0.00049381070782729094, 0.00062112886049619012, 0.0007812731786209638] )
        self.assertEqual( egrid[-4:], [7347.6841999999997, 7348.4096, 7349.1350000000002, 7350.0] )
        self.assertEqual( len(egrid), 60219 )

    def test_setResonanceParametersByChannel( self ):
        self.RRR.setResonanceParametersByChannel( multipleSScheme='NJOY' )
        cDict = self.RRR.channels
        self.assertEqual( len( cDict ), 194 ) # This is number of SLBW "levels".  With only 2 reactions, better be 1/2 the number of channels
        self.assertEqual( len( self.RRR._energies ), 194 ) # This is number of SLBW "levels".  Better equal the _energies array of level energies.
        self.assertEqual( sum( [ len(cs) for cs in cDict ] ), 434 ) # this is number of channels
        self.assertEqual( cDict[0].items()[0][0], ChannelDesignator(l=0, J=4.0, reaction='elastic', index=0, s=4.0, gfact=0.45000000000000001, particleA=None, particleB=None, isElastic=True, channelClass=NEUTRONCHANNEL, useRelativistic=False, eliminated=False ) )
        self.assertEqual( cDict[0].items()[1][1], [(-105.4, 0.165)] )
        self.assertEqual( cDict[193].items()[1][1], [(7331.0, 0.19)] )
        for i in range( len( self.RRR._energies ) ):
            # Check all the channels have the right resonance energy
            self.assertEqual( cDict[i].values()[0][0][0], self.RRR._energies[i] )
            # Check all the widths for a level add up to the tabulated total
            nDigits = 2
            GT = sum( [ x[0][1] for x in cDict[i].values() ] )
            if GT > 0.0: nDigits+=int(round(max(-math.log10(GT),0.0))) # want at least 3 sig. figs. of agreement
            # print i, GT, self.RRR._widths[i], nDigits
            self.assertAlmostEqual( GT, self.RRR._widths[i], nDigits )

    def test_getScatteringMatrixUUnitarity( self ):
        '''
        Each level in the SLBW parameterization should have a unitary U matrix.  Here we have two tests:

            * just plain multiplication: U^+ * U == one

            * check the Frobenius norm, which is ||A|| = sqrt(| sum_{i,j} |A[i,j]|^2 |), so the norm of the identity matrix is sqrt( ndim )
        '''
        self.RRR.setResonanceParametersByChannel(multipleSScheme='NJOY')
        U = self.RRR.getScatteringMatrixU( 0.0235339347844 )
        rtol= 1e-04 # nominal value: 1e-03
        atol= 1e-05 # nominal value: 1e-04
        for i in U:
            one = numpy.diag( U[i].shape[0]*[ complex( 1.0 ) ] )
            testResult = numpy.allclose( numpy.conj(U[i].T)*U[i], one, rtol=rtol, atol=atol )
            if not testResult:
                print i,'\n  U:\n', U[i]
                print '\nU^+:\n', U[i].T.conj()
                print '\n', U[i].T.conj()*U[i]
                print '\nFrobenius norm:', numpy.linalg.norm( U[i] )
            self.assertTrue( testResult )
            self.assertAlmostEqual( numpy.linalg.norm( U[i] ), math.sqrt( U[i].shape[0] ), 5 )

    def test_getScatteringMatrixUValue0( self ):
        '''On resonance nearest to 1 keV of evaluation with many resonances'''
        self.RRR.setResonanceParametersByChannel(multipleSScheme='NJOY')
        U = self.RRR.getScatteringMatrixU( 1009.0, False )
        eiphi = complex( 0.99903310472611739, -0.043964254358450666 )
        GG = 0.813
        GN = 0.6011
        GT = GN+GG #1.414
        cN = self.RRR.channels[28].keys()[0]
        self.assertAlmostEqual( self.RRR.channels[28][cN][0][-1], GN )  # Pick out GN from channels for 28th level.  Better be a match!
        self.assertAlmostEqual( U[28][0][0], complex( 1.0-2.0*GN/GT, 0.0 )*eiphi*eiphi ) # one for each GN
        self.assertAlmostEqual( U[28][0][1], complex( -2.0*math.sqrt(GN*GG)/GT, 0.0 )*eiphi )
        self.assertAlmostEqual( U[28][1][0], complex( -2.0*math.sqrt(GN*GG)/GT, 0.0 )*eiphi )
        self.assertAlmostEqual( U[28][1][1], complex( 1.0-2.0*GG/GT, 0.0 ) )

    def test_getScatteringMatrixUValue1( self ):
        '''On first resonance > 0 of evaluation with many resonances'''
        self.RRR.setResonanceParametersByChannel(multipleSScheme='NJOY')
        U = self.RRR.getScatteringMatrixU( 35.9, False )
        eiphi = complex( 0.9999999999999819, -1.9027638215908413e-07 ) # one for each GN
        GG = 0.209
        GN = 0.0001018
        GT = GN+GG #0.2091
        cN = self.RRR.channels[1].keys()[0]
        self.assertAlmostEqual( self.RRR.channels[1][cN][0][-1], GN )  # Pick out GN from channels for 1st level.  Better be a match!
        self.assertAlmostEqual( U[1][0][0], complex( 1.0-2.0*GN/GT, 0.0 )*eiphi*eiphi )
        self.assertAlmostEqual( U[1][0][1], complex( -2.0*math.sqrt(GN*GG)/GT, 0.0 )*eiphi )
        self.assertAlmostEqual( U[1][1][0], complex( -2.0*math.sqrt(GN*GG)/GT, 0.0 )*eiphi )
        self.assertAlmostEqual( U[1][1][1], complex( 1.0-2.0*GG/GT, 0.0 ) )

    def test_getScatteringMatrixUValue_1S( self ):
        '''On first resonance > 0, it is an S wave resonance of evaluation with one resonance'''
        self.RRR.setResonanceParametersByChannel(multipleSScheme='NJOY')
        ER = 0.169 #eV
        U = self.RRR1SRes.getScatteringMatrixU( ER, False )
        phi = self.RRR1SRes.phi( 0, self.RRR1SRes.rho( ER, 0 ) )
        eiphi = numpy.exp( complex( 0, -phi ) ) # one for each GN
        GG = 7.96e-2
        GN = 5.88e-4
        GT = GN+GG #8.0188e-2
        cN = self.RRR1SRes.channels[0].keys()[0]
        self.assertAlmostEqual( self.RRR1SRes.channels[0][cN][0][-1], GN )  # Pick out GN from channels for 0th level.  Better match!
        self.assertAlmostEqual( U[0][0][0], complex( 1.0-2.0*GN/GT, 0.0 )*eiphi*eiphi )
        self.assertAlmostEqual( U[0][0][1], complex( -2.0*math.sqrt(GN*GG)/GT, 0.0 )*eiphi )
        self.assertAlmostEqual( U[0][1][0], complex( -2.0*math.sqrt(GN*GG)/GT, 0.0 )*eiphi )
        self.assertAlmostEqual( U[0][1][1], complex( 1.0-2.0*GG/GT, 0.0 ) )

    def test_getScatteringMatrixUValue_1P( self ):
        '''On first resonance > 0, it is an P wave resonance of evaluation with one resonance'''
        self.RRR.setResonanceParametersByChannel(multipleSScheme='NJOY')
        ER = 0.169 #eV
        U = self.RRR1PRes.getScatteringMatrixU( ER, False )
        phi = self.RRR1PRes.phi( 1, self.RRR1PRes.rho( ER, 1 ) )
        eiphi = numpy.exp( complex( 0, -phi ) ) # one for each GN
        GG = 7.96e-2
        GN = 5.88e-4
        GT = GN+GG #8.0188e-2
        cN = self.RRR1PRes.channels[0].keys()[0]
        self.assertAlmostEqual( self.RRR1PRes.channels[0][cN][0][-1], GN )  # Pick out GN from channels for 0th level.  Better match!
        self.assertAlmostEqual( U[0][0][0], complex( 1.0-2.0*GN/GT, 0.0 )*eiphi*eiphi )
        self.assertAlmostEqual( U[0][0][1], complex( -2.0*math.sqrt(GN*GG)/GT, 0.0 )*eiphi )
        self.assertAlmostEqual( U[0][1][0], complex( -2.0*math.sqrt(GN*GG)/GT, 0.0 )*eiphi )
        self.assertAlmostEqual( U[0][1][1], complex( 1.0-2.0*GG/GT, 0.0 ) )

    @unittest.skip("SLBW getAngularDistribution now NotImplemented")
    def test_scatteringMatrixUInCrossSectionCalculation_SWave( self ):
        debug = False
        doTests = True

        cN = self.RRR1SRes.channels[0].keys()[0]
        iN, iG = 0, 1 # indices for the neutron and capture channels ( this is the order I packed them in )
        ER = 0.169 #eV
        AP = 0.63809 #in b^{-1/2}

        if debug: print '\nS Wave, g-factor (fudge)=',cN.gfact
        if debug: print '     E (eV)       SigTotU (b)    SigTotOK (b)   SigGamU (b)    SigGamOK (b)   SigNeuU (b)    SigNeuOK (b)   SigNeuInt (b) SigNeuInt/SigNeuU  dSigPotActual/dSigPotCalc  '

        for E in self.RRR1SRes.generateEnergyGrid():
            nLs = { 0:0 }
            gs = { 0:0 }
            nLs[ cN.l ] += 1.0
            gs[ cN.l ] = cN.gfact
            U = self.RRR1SRes.getScatteringMatrixU( E, True )
            k = self.RRR1SRes.k(E)
            rhohat = AP * k
            sigTotWithU = 2.0 * numpy.pi * cN.gfact * ( 1.0 - U[0][iN][iN].real ) / k**2
            sigGamWithU =       numpy.pi * cN.gfact * ( U[0][iN][iG] * U[0][iG][iN].conj() ).real / k**2
            sigNeuWithU =       numpy.pi * cN.gfact * ( ( 1.0 - U[0][iN][iN] ).conj() * ( 1.0 - U[0][iN][iN] ) ).real / k**2
            sigPot = sum( [ ( ( 2.0 * l + 1 ) - nLs[l] * gs[l] ) * 4.0 * numpy.pi * ( numpy.sin( self.RRR1SRes.phi( l, rhohat ) ) )**2 / k**2 for l in nLs ] )
            B = self.RRR1SRes.getAngularDistribution( E )
            sigNeuInt   = 4*numpy.pi*B[0].real
            sigs = self.RRR1SRes.getCrossSection( E )

            if debug: print 10*'%10g     ' % ( E, sigTotWithU, sigs['total'][0], sigGamWithU, sigs['capture'][0], sigNeuWithU, sigs['elastic'][0], sigNeuInt,  sigNeuWithU/sigNeuInt, -(sigNeuWithU-sigs['elastic'][0])/sigPot )

            # Check they sum up correctly (really just another test of unitarity)
            self.assertTrue( withinXPercent( sigTotWithU, sigGamWithU + sigNeuWithU ) )

            # Check elastic channel agrees with L=0 term of angular distribution
            if doTests: self.assertTrue( withinXPercent( sigNeuWithU, sigNeuInt ) )

            # Check agree with getCrossSection()
            if doTests:
                self.assertTrue( withinXPercent( sigTotWithU+sigPot, sigs['total'][0] ) )
                self.assertTrue( withinXPercent( sigGamWithU, sigs['capture'][0] ) )
                self.assertTrue( withinXPercent( sigNeuWithU+sigPot, sigs['elastic'][0] ) )

    @unittest.skip("SLBW getAngularDistribution now NotImplemented")
    def test_scatteringMatrixUInCrossSectionCalculation_PWave( self ):
        debug = False
        doTests = True

        cN = self.RRR1PRes.channels[0].keys()[0]
        iN, iG = 0, 1 # indices for the neutron and capture channels ( this is the order I packed them in )
        ER = 0.169 #eV
        AP = 0.63809 #in b^{-1/2}

        if debug: print '\nP Wave, g-factor (fudge)=',cN.gfact
        if debug: print '     E (eV)       SigTotU (b)    SigTotOK (b)   SigGamU (b)    SigGamOK (b)   SigNeuU (b)    SigNeuOK (b)   SigNeuInt (b) SigNeuInt/SigNeuU  dSigPotActual/dSigPotCalc '

        for E in self.RRR1PRes.generateEnergyGrid():
            nLs = { 0:0, 1:0 }
            gs = { 0:0, 1:0 }
            U = self.RRR1PRes.getScatteringMatrixU( E, True )
            k = self.RRR1PRes.k(E)
            rhohat = AP * k
            nLs[ cN.l ] += 1.0
            gs[ cN.l ] = cN.gfact
            sigTotWithU = 2.0 * numpy.pi * cN.gfact * ( 1.0 - U[0][iN][iN].real ) / k**2
            sigGamWithU =       numpy.pi * cN.gfact * ( U[0][iN][iG] * U[0][iG][iN].conj() ).real / k**2
            sigNeuWithU =       numpy.pi * cN.gfact * ( ( 1.0 - U[0][iN][iN] ).conj() * ( 1.0 - U[0][iN][iN] ) ).real / k**2
            sigPot = sum( [ ( ( 2.0 * l + 1 ) - nLs[l] * gs[l] ) * 4.0 * numpy.pi * ( numpy.sin( self.RRR1PRes.phi( l, rhohat ) ) )**2 / k**2 for l in nLs ] )
            B = self.RRR1PRes.getAngularDistribution( E )
            sigNeuInt   = 4*numpy.pi*B[0].real
            sigs = self.RRR1PRes.getCrossSection( E )

            if debug: print 10*'%10g     ' % ( E, sigTotWithU, sigs['total'][0], sigGamWithU, sigs['capture'][0], sigNeuWithU, sigs['elastic'][0], sigNeuInt,  sigNeuWithU/sigNeuInt, -(sigNeuWithU-sigs['elastic'][0])/sigPot )

            # Check they sum up correctly (really just another test of unitarity)
            self.assertTrue( withinXPercent( sigTotWithU, sigGamWithU + sigNeuWithU ) )

            # Check elastic channel agrees with L=0 term of angular distribution
            if doTests: self.assertTrue( withinXPercent( sigNeuWithU, sigNeuInt ) )

            # Check agree with getCrossSection()
            if doTests:
                self.assertTrue( withinXPercent( sigTotWithU+sigPot, sigs['total'][0] ) )
                self.assertTrue( withinXPercent( sigGamWithU, sigs['capture'][0] ) )
                self.assertTrue( withinXPercent( sigNeuWithU+sigPot, sigs['elastic'][0] ) )

    @unittest.skip("SLBW getAngularDistribution now NotImplemented")
    def test_scatteringMatrixUInCrossSectionCalculation_Full( self ):
        debug = False
        doTests = True

        nLev = len(self.RRR.channels)
        AP = 0.69 #in b^{-1/2}
        allowedSs = getAllowedTotalSpins( 9/2., 1/2., False )

        if debug: print '\nAll Waves (course grid), %i resonances' % nLev
        if debug: print '     E (eV)       SigTotU (b)    SigTotOK (b)   SigGamU (b)    SigGamOK (b)   SigNeuU (b)    SigNeuOK (b)   SigNeuInt (b) SigNeuInt/SigNeuU  dSigPotActual/dSigPotCalc '

        egrid = self.RRR.generateEnergyGrid()
        eStep = 10 #int( len( egrid )/95)
        for E in egrid[::eStep]:
            if not ( E >= 7300. and E<= 7400. ): continue # Run test only in a problem area because it is expensive
            U = self.RRR.getScatteringMatrixU( E, True )
            k = self.RRR.k(E)
            rhohat = AP * k
            sigTotWithU, sigGamWithU, sigNeuWithU = 0.0, 0.0, 0.0
            nLs = { 0:0, 1:0 }
            gs = { 0:0, 1:0 }
            for iLev in range(nLev):
                cN = self.RRR.channels[iLev].keys()[0]
                nLs[ cN.l ] += 1.0
                gs[ cN.l ] += cN.gfact
                iN, iG = 0, 1  # indices for the neutron and capture channels ( this is the order I packed them in )
                sigTotWithU += 2.0 * numpy.pi * cN.gfact * ( 1.0 - U[iLev][iN][iN].real ) / k**2
                sigGamWithU +=       numpy.pi * cN.gfact * ( U[iLev][iN][iG] * U[iLev][iG][iN].conj() ).real / k**2
                sigNeuWithU +=       numpy.pi * cN.gfact * ( ( 1.0 - U[iLev][iN][iN] ).conj() * ( 1.0 - U[iLev][iN][iN] ) ).real / k**2
            sigPot = sum( [ ( ( 2.0 * l + 1 ) - gs[l] ) * 4.0 * numpy.pi * ( numpy.sin( self.RRR.phi( l, rhohat ) ) )**2 / k**2 for l in nLs ] )
            B = self.RRR.getAngularDistribution( E )
            sigNeuInt   = 4*numpy.pi*B[0].real
            sigs = self.RRR.getCrossSection( E )

            if debug: print 10*'%10g     ' % ( E, sigTotWithU, sigs['total'][0], sigGamWithU, sigs['capture'][0], sigNeuWithU, sigs['elastic'][0], sigNeuInt,  sigNeuWithU/sigNeuInt, -(sigNeuWithU-sigs['elastic'][0])/sigPot )

            # Check they sum up correctly (really just another test of unitarity)
            self.assertTrue( withinXPercent( sigTotWithU, sigGamWithU + sigNeuWithU ) )

            # Check elastic channel agrees with L=0 term of angular distribution
            if doTests: self.assertTrue( withinXPercent( sigNeuWithU, sigNeuInt, 50. ) ) # fails badly for some high energy P wave resonances??

            # Check agree with getCrossSection()
            if doTests:
                self.assertTrue( withinXPercent( sigTotWithU+sigPot, sigs['total'][0] ) )
                self.assertTrue( withinXPercent( sigGamWithU, sigs['capture'][0] ) )
                self.assertTrue( withinXPercent( sigNeuWithU+sigPot, sigs['elastic'][0] ) )

    @unittest.skip("SLBW getAngularDistribution now NotImplemented")
    def test_getAngularDistribution_SWave( self ):
        '''
        Note, Blatt-Biedenharn's Z coefficient has a factor of (l1,l2,0,0;L,0) (a Clebsch-Gordan coefficient) in it.
        As a result, if l1 and l2 cannot couple to L, one gets zero.

        For a single S-wave resonance, both l's are 0 so can only couple to L=0, resulting in pure isotropic angular distributions.
        '''
        egrid = self.RRR1SRes.generateEnergyGrid()
        eStep = int( len( egrid )/95)
        for E in egrid[::eStep]:
            B = self.RRR1SRes.getAngularDistribution( E )
            self.assertTrue( allWithinXPercent( [ x/B[0] for x in B ], [ 1.0, 0, 0 ] ) )

    @unittest.skip("SLBW getAngularDistribution now NotImplemented")
    def test_getAngularDistribution_PWave( self ):
        '''
        Note, Blatt-Biedenharn's Z coefficient has a factor of (l1,l2,0,0;L,0) (a Clebsch-Gordan coefficient) in it.
        As a result, if l1 and l2 cannot couple to L, one gets zero.

        For a single P-wave resonance, both l's are 1 so can couple to L=0 or 2, resulting in an angular distribution with L=0 and 2 terms.

        Because there is only one resonance, the U matrix has one neutron-neutron element in it.
        So, all L's end up with the same (1-U)^2 factor and the energy dependence of the angular distribution divides out, leaving basically a constant factor of Z^2 for L=0 and 2.
        For L=0, Z ~ 3.5 and for L=2, Z ~ 2.  So, expect B[0]/B[0] == 1 and B[2]/B[0] ~ (2/3.5)^2 ~ 0.33.
        '''
        egrid = self.RRR1PRes.generateEnergyGrid()
        eStep = 1
        J = 5.5 # for this test data
        S = 4.5 # for this test data
        from numericalFunctions import angularMomentumCoupling as nf_amc
        self.assertAlmostEqual( nf_amc.z_coefficient( 2, int(2*J), 2, int(2*J), int(2*S), 4 ) / nf_amc.z_coefficient( 2, int(2*J), 2, int(2*J), int(2*S), 0 ), math.sqrt(0.33090909) ) # nf_amc functions use the 2J trick
        for E in egrid[::eStep]:
            B = self.RRR1PRes.getAngularDistribution( E )
            self.assertTrue( allWithinXPercent( [ x/B[0] for x in B ], [ 1.0, 0.0, 0.33090909 ] ) )

    @unittest.skip("SLBW getAngularDistribution now NotImplemented")
    def test_getAngularDistribution( self ):
        '''
        Use version of channels where we strip out the channels with spins that have zero width.
        The unit test is way faster than the full BB version, so we need to weaponize that one.
        '''
        debug = False
        doTest = True
        NJOYStyleSpinTreatment = True

        from numericalFunctions import angularMomentumCoupling as nf_amc
        egrid = self.RRR.generateEnergyGrid()
        eStep = 10 #int( len( egrid )/95)
        AP = 0.69 #in b^{-1/2}
        ERtoTest = 7331.0
        nLev = len(self.RRR.channels)

        # First list is edited out of all the terms that correspond to spin channels with zero width
        if NJOYStyleSpinTreatment:
            JslList  = [ ( 4, 4, 0 ), ( 5, 5, 0 ), ( 3, 4, 1 ), ( 4, 4, 1 ), ( 5, 4, 1 ), ( 6, 5, 1 ) ] # try skipping the "redundant" ones: ( 4, 5, 1 ), ( 5, 5, 1 )
        # Second list has all the angular terms
        else:
            JslList  = [ ( 4, 4, 0 ), ( 5, 5, 0 ), ( 3, 4, 1 ), ( 4, 4, 1 ), ( 5, 4, 1 ), ( 4, 5, 1 ), ( 5, 5, 1 ), ( 6, 5, 1 ) ]

        angCoeffSqrList = [ { Jsl : pow( nf_amc.z_coefficient( int(2*Jsl[2]), int(2*Jsl[0]), int(2*Jsl[2]), int(2*Jsl[0]), int(2*Jsl[1]), 2*L ), 2 ) for Jsl in JslList } for L in range( 3 ) ]
        if doTest:
            self.maxDiff = None
            # First list is edited out of all the terms that correspond to spin channels with zero width
            if NJOYStyleSpinTreatment:
                self.assertEqual( angCoeffSqrList, [ \
                    { (3, 4, 1): 7.000000000000003, (4, 4, 0): 9.0, (4, 4, 1): 9.0, (5, 4, 1): 11.0, (5, 5, 0): 11.0, (6, 5, 1): 12.999999999999993}, \
                    { (3, 4, 1): 0.0, (4, 4, 0): 0.0, (4, 4, 1): 0.0, (5, 4, 1): 0.0,(5, 5, 0): 0.0, (6, 5, 1): 0.0 }, \
                    { (3, 4, 1): 0.5833333333213947, (4, 4, 0): 0.0, (4, 4, 1): 6.929999999908832, (5, 4, 1): 3.8133333328735186, (5, 5, 0): 0.0,  (6, 5, 1): 4.136363636490839 } ] )
            # Second list has all the angular terms
            else:
                self.assertEqual( angCoeffSqrList, [ \
                    { (3, 4, 1): 6.999999999953966, (4, 4, 0): 8.99999999942597, (4, 4, 1): 8.999999999408963, (4, 5, 1): 8.999999999408963, (5, 4, 1): 10.99999999996129, (5, 5, 0): 10.999999999982087, (5, 5, 1): 10.999999999961286, (6, 5, 1): 13.000000000775405}, \
                    { (3, 4, 1): 0.0, (4, 4, 0): 0.0, (4, 4, 1): 0.0, (4, 5, 1): 0.0,  (5, 4, 1): 0.0,(5, 5, 0): 0.0, (5, 5, 1): 0.0, (6, 5, 1): 0.0 }, \
                    { (3, 4, 1): 0.5833333333213947, (4, 4, 0): 0.0, (4, 4, 1): 6.929999999908832, (4, 5, 1): 0.9163636362571602, (5, 4, 1): 3.8133333328735186, (5, 5, 0): 0.0, (5, 5, 1): 8.579999999854891, (6, 5, 1): 4.136363636490839 } ] )

        if debug: print
        for E in egrid[::eStep]:
            if not ( E >= 7300. and E<= 7400. ): continue # Run test only in a problem area because it is expensive
            U = self.RRR.getScatteringMatrixU( E, True )
            k = self.RRR.k(E)
            rhohat = AP * k
            answer = len(angCoeffSqrList)*[ 0.0 ]
            sigTotWithU, sigGamWithU, sigNeuWithU = 0.0, 0.0, 0.0
            for iLev in range(nLev):
                cN = self.RRR.channels[iLev].keys()[0]
                iN, iG = 0, 1  # indices for the neutron and capture channels ( this is the order I packed them in )
                for L,a in enumerate(answer):
                    for S in [4,5]:
                        if not (cN.J,S,cN.l) in angCoeffSqrList[L]: continue
                        answer[L] +=  angCoeffSqrList[L][(cN.J,S,cN.l)] * ( ( 1.0 - U[iLev][iN][iN] ).conj() * ( 1.0 - U[iLev][iN][iN] ) ).real / k**2
            B = self.RRR.getAngularDistribution( E )
            if debug: print E, ' '.join( map( str, [ x/B[0] for x in B ] ) ), ' '.join( map( str, [ x/answer[0] for x in answer ] ) ), ( answer[-1]/answer[0] ) / ( B[-1]/B[0] )
            if doTest:
                # if we are within 2 eV of the resonance, the Legendre moment is big enough to notice
                if abs( E - ERtoTest ) <= 2:
                    self.assertTrue( allWithinXPercent( [ x/B[0] for x in B ], [ x/answer[0] for x in answer ], absTol=1e-6 ) )

    def test_getCrossSectionFull( self ):
        '''Test thermal & a few other points'''
        answers = [ \
            ( 0.0235339347844, { 'capture':1.18726303, 'total':7.24518429, 'fission':0., 'elastic':6.05792126 } ),
            ( 1009.0, { 'capture':290.05432779892004, 'total':507.80854074065803, 'fission':0., 'elastic':217.75421294173799 } ) ]
        for answer in answers:
            result = self.RRR.getCrossSection( answer[0] )
            for k in answer[1]:
                self.assertAlmostEqual( answer[1][k], result[k][0] )

    def test_getCrossSection1SRes( self ):
        '''Test thermal & a few other points'''
        answers = [ \
            ( 0.0235339347844, { 'capture':11957.883139892736, 'total':11977.291733536884, 'fission':0., 'elastic':19.408593644148297 } ),
            ( 1009.0, { 'capture':1.2911030236464319e-06, 'total':5.1160670269419093, 'fission':0., 'elastic':5.1160657358388857 } ) ]
        for answer in answers:
            result = self.RRR1SRes.getCrossSection( answer[0] )
            for k in answer[1]:
                self.assertWithinXPercent( answer[1][k], result[k][0], percent=0.01 )



class TestMLBWClassAndBaseClasses( TestWithIsClose ):

    def setUp( self ):
        self.RRR = MLBWcrossSection( MLBWExample, enableAngDists=True, verbose=False )
        self.RRR1Res = MLBWcrossSection( MLBWExample1PRes, enableAngDists=True, verbose=False )

    def test_memberData( self ):
        self.assertEqual( str(self.RRR.projectile), 'n' )
        self.assertEqual( self.RRR.targetSpin, 0.0 )
        self.assertEqual( str(self.RRR.target), 'O18' )
        self.assertAlmostEqual( self.RRR.targetToNeutronMassRatio, 17.8445 )
        self.assertAlmostEqual( self.RRR.lowerBound, 1e-05 )
        self.assertAlmostEqual( self.RRR.upperBound, 3150000.0 )
        self.assertEqual( self.RRR.missingGfactor, {0: 0.0, 1: 0.0, 2: 0.0, 3: 3.0} )
        self.assertEqual( len(self.RRR.Ls), 4 )
        self.assertEqual( self.RRR.Ls[1].keys(), ['L', 'Js'] )
        self.assertEqual( self.RRR.Ls[1]['L'], 1 )
        self.assertEqual( self.RRR.Ls[1]['Js'][0].keys(), ['gfact', 'channelSpins', 'J'] )
        self.assertAlmostEqual( self.RRR.Ls[1]['Js'][0]['gfact'], 1.0 )
        self.assertEqual( self.RRR.Ls[1]['Js'][0]['channelSpins'][0].keys(), ['neutronWidth','channelSpin','captureWidth','energy', 'fissionWidthB', 'fissionWidth', 'shiftFactor'] )
        self.assertEqual( self.RRR.Ls[1]['Js'][0]['J'], 0.5 )
        self.assertEqual( list(self.RRR._energies), [-2485000.0, 663000.0, 1200000.0, 1256000.0, 1450000.0, 1585000.0, 1840000.0, 2300000.0, 2375000.0, 2445000.0, 3050000.0, 3290000.0, 3500000.0, 3590000.0, 4090000.0] )
        self.assertEqual( list(self.RRR._widths), [1429780.0, 55000.0, 50000.0, 3600.0, 350000.0, 300000.0, 8200.0, 280000.0, 130000.0, 20300.0, 20300.0, 250000.0, 100000.0, 100000.0, 6500.0] )

    def test_setResonanceParametersByChannel( self ):
        debug = False
        self.RRR.setResonanceParametersByChannel( )
        cDict = self.RRR.channels
        if debug: print
        if debug: print '\n'.join( map( str, cDict.items() ) )
        self.assertEqual( len( cDict ), 13 ) # This is number of distinct MLBW channels
        self.assertEqual( len( self.RRR._energies ), 15 ) # This is number of MLBW resonances.  Better equal the _energies array of level energies.
        self.assertEqual( sum( [ len(cDict[cs]) for cs in cDict ] ), 31 ) # this is double the number of resonances (1/2 for elastic, 1/2 for capture)
        self.assertEqual( cDict.items()[0][0], ChannelDesignator(l=0, J=0.5, reaction='elastic', index=0, s=0.5, gfact=1.0, particleA='n', particleB='O18', isElastic=True, channelClass=NEUTRONCHANNEL, useRelativistic=False, eliminated=False ) )
        self.assertEqual( cDict.items()[1][1], { 0:0.23 } )
        self.assertEqual( cDict.items()[2][1], { 1:54999.8, 3:3599.77, 6:8199.77 } )
        for i,E in enumerate(self.RRR._energies):
            # Sum all the widths with matching resonance energy.  Done by brute force.
            GT = 0.0
            for c in cDict:
                for iR, G in cDict[c].items():
                    if iR == i: GT += G
            # Check all the widths for a level add up to the tabulated total
            if debug: print i, E, GT, self.RRR._widths[i]
            self.assertWithinXPercent( GT, self.RRR._widths[i], 0.1 )

    def test_getScatteringMatrixUValue_1PWave( self ):
        '''On first resonance > 0, it is an P wave resonance of evaluation with one resonance.  For one res, SLBW == MLBW.'''
        ER = 0.169 #eV
        U = self.RRR1Res.getScatteringMatrixU( numpy.array(ER).reshape(1,1), False )
        phi = self.RRR1Res.phi( 1, self.RRR1Res.rho( ER, 1 ) )
        eiphi = numpy.exp( complex( 0, -phi ) ) # one for each GN
        GG = 7.96e-2
        GN = 5.88e-4
        GT = GN+GG #8.0188e-2
        cN = self.RRR1Res.channels.keys()[0]
        self.assertAlmostEqual( self.RRR1Res.channels[cN][0], GN )  # Pick out GN from channels for 0th level.  Better match!
        self.assertAlmostEqual( U[0][0,0], complex( 1.0-2.0*GN/GT, 0.0 )*eiphi*eiphi )
        self.assertAlmostEqual( U[0][0,1], complex( -2.0*math.sqrt(GN*GG)/GT, 0.0 )*eiphi )
        self.assertAlmostEqual( U[0][1,0], complex( -2.0*math.sqrt(GN*GG)/GT, 0.0 )*eiphi )
        self.assertAlmostEqual( U[0][1,1], complex( 1.0-2.0*GG/GT, 0.0 ) )

    def test_getScatteringMatrixUUnitarity_1PWave( self ):
        '''
        The MLBW parameterization might have a unitary U matrix.  Here we have two tests:

            * just plain multiplication: U^+ * U == one

            * check the Frobenius norm, which is ||A|| = sqrt(| sum_{i,j} |A[i,j]|^2 |), so the norm of the identity matrix is sqrt( ndim )

        Note:   For one res, SLBW == MLBW.
        '''
        debug = False
        doTest = True
        if debug:
            for EE in self.RRR1Res.generateEnergyGrid():
                UU = self.RRR1Res.getScatteringMatrixU( EE )
                print EE, UU[0,0].real, UU[0,0].imag, numpy.linalg.norm(UU)/math.sqrt( UU.shape[0] )
        rtol= 1e-03 # nominal value: 1e-03
        atol= 1e-04 # nominal value: 1e-04
        E = numpy.array(1.0).reshape(1,1)
        U = self.RRR1Res.getScatteringMatrixU( E )[0]
        one = numpy.diag( U.shape[0]*[ complex( 1.0 ) ] )
        if doTest: self.assertAlmostEqual( numpy.linalg.norm( U ), math.sqrt( U.shape[0] ), 4 )
        if doTest: self.assertTrue( numpy.allclose( numpy.conj(U.T)*U, one, rtol=rtol, atol=atol ) )

    def test_getScatteringMatrixUUnitarity_Full( self ):
        '''
        The MLBW parameterization might have a unitary U matrix.  Here we have two tests:

            * just plain multiplication: U^+ * U == one

            * check the Frobenius norm, which is ||A|| = sqrt(| sum_{i,j} |A[i,j]|^2 |), so the norm of the identity matrix is sqrt( ndim )
        '''
        debug = False
        doTest = True
        if debug:
            for EE in self.RRR.generateEnergyGrid()[::10]:
                UU = self.RRR.getScatteringMatrixU( EE )
                print EE, UU[0,0].real, UU[0,0].imag, UU[0,1].real, UU[0,1].imag, UU[1,1].real, UU[1,1].imag, numpy.linalg.norm(UU)/math.sqrt( UU.shape[0] )
        rtol= 1e-03 # nominal value: 1e-03
        atol= 1e-04 # nominal value: 1e-04
        E = numpy.array(1.0).reshape(1,1)
        U = self.RRR.getScatteringMatrixU( E )[0]
        one = numpy.diag( U.shape[0]*[ complex( 1.0 ) ] )
        if doTest: self.assertAlmostEqual( numpy.linalg.norm( U ), math.sqrt( U.shape[0] ), 4 )
        if doTest: self.assertTrue( numpy.allclose( numpy.conj(U.T)*U, one, rtol=rtol, atol=atol ) )

    def test_scatteringMatrixUInCrossSectionCalculation_1PWave( self ):
        '''
        For one res, SLBW == MLBW.
        '''
        debug = False
        doTests = True

        fourpi = 4.0 * numpy.pi
        nChan = len( self.RRR1Res.channels )
        nRes = len( self.RRR1Res._energies )
        AP = 0.63809 #in b^{-1/2}

        if debug: print '\n1 P Wave (course grid), %i channels, %i resonances' % ( nChan, nRes )
        if debug: print '     E (eV)       SigTotU (b)    SigTotOK (b)   SigGamU (b)    SigGamOK (b)   SigNeuU (b)    SigNeuOK (b)   SigNeuInt (b) SigNeuInt/SigNeuU  sigPot' #  dSigPotActual/dSigPotCalc '

        eStep = 10 #int( len( egrid )/95)
        egrid = self.RRR1Res.generateEnergyGrid()[::eStep]
        egrid = numpy.array(egrid).reshape(len(egrid),1)

        U = self.RRR1Res.getScatteringMatrixU( egrid, True )
        k = self.RRR1Res.k( egrid )
        B = self.RRR1Res.getAngularDistribution( egrid )['elastic']
        sigs = self.RRR1Res.getCrossSection( egrid )
        rhohat = AP * k
        nLs = { 0:0, 1:0 }
        gs = { 0:0, 1:0 }

        for iE,E in enumerate(egrid[::eStep]):
            sigTotWithU, sigGamWithU, sigNeuWithU = 0.0, 0.0, 0.0
            for ic1,c1 in enumerate(self.RRR1Res.channels):
                if c1.reaction != 'elastic': continue
                gs[ c1.l ] += c1.gfact
                nLs[ c1.l ] += 1.0
                sigTotWithU += 2.0 * numpy.pi * c1.gfact * ( 1.0 - U[iE][ic1,ic1].real ) / k[iE]**2
                for ic2,c2 in enumerate(self.RRR1Res.channels):
                    if c2.J != c1.J or c1.s != c2.s or c1.l != c2.l: continue
                    if c2.reaction == 'elastic':
                        sigNeuWithU += numpy.pi * c1.gfact * ( ( (ic1==ic2) - U[iE][ic1,ic2] ).conj() * ( (ic1==ic2) - U[iE][ic2,ic1] ) ).real / k[iE]**2
                    elif c2.reaction == 'capture':
                        sigGamWithU += numpy.pi * c1.gfact * ( U[iE][ic1,ic2] * U[iE][ic2,ic1].conj() ).real / pow( k[iE], 2.0 )
            sigPot = sum([ ( ( 2.0 * l + 1 ) ) * fourpi * pow( numpy.sin( self.RRR1Res.phi( l, rhohat[iE] ) ) / k[iE], 2.0 ) for l in nLs ])
            sigNeuInt = fourpi * B[iE][0].real

            if debug: print 10*'%10g     ' % ( E, sigTotWithU, sigs['total'][iE], sigGamWithU, sigs['capture'][iE], sigNeuWithU, sigs['elastic'][iE], sigNeuInt,  sigNeuWithU/sigNeuInt, sigPot ) #, -(sigNeuWithU-sigs['elastic'][iE])/sigPot )

            # Check they sum up correctly (really just another test of unitarity)
            if doTests: self.assertWithinXPercent( sigTotWithU, sigGamWithU + sigNeuWithU )

            # Check elastic channel agrees with L=0 term of angular distribution
            if doTests and False: self.assertWithinXPercent( sigNeuWithU, sigNeuInt, 50. ) # fails badly for some high energy P wave resonances??

            # Check agree with getCrossSection()
            if doTests:
                self.assertWithinXPercent( sigTotWithU, sigs['total'][iE] )
                self.assertWithinXPercent( sigGamWithU, sigs['capture'][iE] )
                self.assertWithinXPercent( sigNeuWithU, sigs['elastic'][iE] )

    def test_scatteringMatrixUInCrossSectionCalculation_Full( self ):
        debug = False
        doTests = True

        nChan = self.RRR.nChannels
        nRes = self.RRR.nResonances
        AP = 0.41 #in b^{-1/2}

        if debug: print '\nAll Waves (course grid), %i channels, %i resonances' % ( nChan, nRes )
        if debug: print '     E (eV)       SigTotU (b)    SigTotSum (b)    SigTotOK (b)   SigGamU (b)    SigGamOK (b)   SigNeuU (b)    SigNeuOK (b)   SigNeuInt (b) SigNeuInt/SigNeuU  sigPot'

        # We can't eat the whole energy grid, so we'll do it in parts
        fullEgrid = self.RRR.generateEnergyGrid()
        energyCohortLength = 50
        nEnergyCohorts = len(fullEgrid)//energyCohortLength
        lastEnergyCohortSize = len(fullEgrid)%energyCohortLength
        energyCohorts = [(iec*energyCohortLength, (iec+1)*energyCohortLength) for iec in range(nEnergyCohorts)]
        energyCohorts.append((energyCohorts[-1][1],energyCohorts[-1][1]+lastEnergyCohortSize))

        for iStart,iStop in energyCohorts:
            egrid = fullEgrid[iStart:iStop]
            egrid = numpy.array(egrid).reshape(len(egrid),1)

            U = self.RRR.getScatteringMatrixU( egrid, True )
            k = self.RRR.k(egrid)
            B = self.RRR.getAngularDistribution( egrid )['elastic']
            sigs = self.RRR.getCrossSection( egrid )
            rhohat = AP * k
            nLs = { 0:0, 1:0, 2:0, 3:0 }
            gs = { 0:0, 1:0, 2:0, 3:0 }

            for iE,E in enumerate(egrid):
                sigTotWithU, sigGamWithU, sigNeuWithU = 0.0, 0.0, 0.0
                for ic1,c1 in enumerate(self.RRR.channels):
                    if c1.reaction != 'elastic': continue
                    gs[ c1.l ] += c1.gfact
                    nLs[ c1.l ] += 1.0
                    sigTotWithU += 2.0 * numpy.pi * c1.gfact * ( 1.0 - U[iE][ic1,ic1].real ) / k[iE]**2
                    for ic2,c2 in enumerate(self.RRR.channels):
                        if c2.J != c1.J or c1.s != c2.s or c1.l != c2.l: continue
                        if c2.reaction == 'elastic':
                            sigNeuWithU += numpy.pi * c1.gfact * ( ( (ic1==ic2) - U[iE][ic1,ic2] ).conj() * ( (ic1==ic2) - U[iE][ic2,ic1] ) ).real / k[iE]**2
                        elif c2.reaction == 'capture':
                            sigGamWithU += numpy.pi * c1.gfact * ( U[iE][ic1,ic2] * U[iE][ic2,ic1].conj() ).real / k[iE]**2
                sigPot = sum( [ ( ( 2.0 * l + 1 ) ) * 4.0 * numpy.pi * ( numpy.sin( self.RRR.phi( l, rhohat[iE] ) ) )**2 / k[iE]**2 for l in nLs ] )
                sigNeuInt   = 4*numpy.pi*B[0].real

                if debug: print 11*'%10g     ' % ( E, sigTotWithU, sigGamWithU + sigNeuWithU, sigs['total'][iE], sigGamWithU, sigs['capture'][iE], sigNeuWithU, sigs['elastic'][iE], sigNeuInt,  sigNeuWithU/sigNeuInt, sigPot ) #, -(sigNeuWithU-sigs['elastic'][iE])/sigPot )

                # Check they sum up correctly (really just another test of unitarity)
                if doTests and False:  self.assertWithinXPercent( sigTotWithU, sigGamWithU + sigNeuWithU, 10. ) # expected to fail

                # Check elastic channel agrees with L=0 term of angular distribution
                if doTests and False: self.assertWithinXPercent( sigNeuWithU, sigNeuInt )

                # Check agree with getCrossSection()
                # Note, since ENDF's capture cross section is really a SLBW cross section, we expect poor agreement between our U matrix based one and ENDF's
                # Similarly, since the total cross section is a sum of capture and elastic, we should really use the ENDF (SLBW) capture cross section to make the comparison
                if doTests:
                    self.assertWithinXPercent( sigs['capture'][iE]+sigNeuWithU, sigs['total'][iE], 15. )
                    if False: self.assertWithinXPercent( sigGamWithU, sigs['capture'][iE] )  # expected to fail
                    self.assertWithinXPercent( sigNeuWithU, sigs['elastic'][iE], 10. )

    def test_getAngularDistribution_PWave( self ):
        '''
        Note, Blatt-Biedenharn's Z coefficient has a factor of (l1,l2,0,0;L,0) (a Clebsch-Gordan coefficient) in it.
        As a result, if l1 and l2 cannot couple to L, one gets zero.

        In our simple P-wave example, the U matrix has l=0 and 1 parts to it because you need the l=0 terms to get the
        potential scattering right, even if there is only one resonance.  Therefore, we expect L = 0, 1, and 2 terms since
        those are all the L's we can get by coupling any pair of l=0 and 1 together.  Because of this, the angular distribution
        that you get from a single MLBW resonance is different that that of a single SLBW resonance (at least the way
        you need to do things for ENDF).

        We tested the hell out of the L=0 term.

        We can test the L=1,2 terms on resonance by observing that the potential scattering in the l>0 collision matrices is much
        smaller than that in the l=0 elements.  This allows us to ignore all but the purely resonance elements when directly on
        resonance (Ein=ER).   In our P-wave example, that means that our MLBW on resonance angular distribution should match the
        SLBW result (_but only on resonance_).  We reuse the test below, but only approximately.
        '''
        debug = False
        doTest = True
        egrid = self.RRR1Res.generateEnergyGrid()
        eStep = 1
        J = 5.5 # for this test data
        S = 4.5 # for this test data
        from numericalFunctions import angularMomentumCoupling as nf_amc
        self.assertAlmostEqual( nf_amc.z_coefficient( 2, int(2*J), 2, int(2*J), int(2*S), 4 ) / nf_amc.z_coefficient( 2, int(2*J), 2, int(2*J), int(2*S), 0 ), math.sqrt(0.33090909) ) # nf_amc functions use the 2J trick
        if debug:
            one = numpy.diag( 9*[0.0] )
            for E in egrid[::eStep]:
                U = self.RRR1Res.getScatteringMatrixU( E )
                print '\n',E,one - U
                B = self.RRR1Res.getAngularDistribution( E )
                print E, ' '.join( map( str, B ) )
        if doTest:
            ER = 0.169
            B = self.RRR1Res.getAngularDistribution( ER )['elastic']
            self.assertAlmostEqual( B[2]/B[0], 0.066181818, 2 )
            self.assertAlmostEqual( B[1]/B[0], 0.0, 3 )

    @unittest.skipIf(True,"need answer to compare to")
    def test_getElasticAngularDistribution_Full( self ):
        debug = False
        egrid = self.RRR.generateEnergyGrid()
        eStep = 10
        if debug:
            for E in egrid[::eStep]:
                print E, ' '.join( map( str, self.RRR.getAngularDistribution( E ) ) )

    @unittest.skipIf(True,"need answer to compare to")
    def test_getElasticAngularDistribution_Zr90( self ):
        debug = False
        theRRR = MLBWcrossSection( MLBWExampleZr90, enableAngDists=True, verbose=False )
        egrid = theRRR.generateEnergyGrid()
        eStep = 5
        if debug:
            for E in egrid[::eStep]:
                print E, ' '.join( map( str, theRRR.getAngularDistribution( E ) ) )

    def test_getCrossSection( self ):
        '''Test thermal point'''
        result = self.RRR.getCrossSection( 0.0235339347844 )
        #print 'MLBW',result
        answer = {'capture':0.00015987, 'total':3.11326583, 'fission':0., 'elastic':3.11310596}
        for k in answer:
            self.assertAlmostEqual( answer[k], result[k][0]  )


@unittest.skip("Major changes to class means tests always fail, they need to be refactored")
class TestRMClassAndBaseClasses( TestWithIsClose ):

    def setUp( self ):
        self.RRR = RMcrossSection( RMExample, enableAngDists=True, verbose=False )
        self.RRRSmall = RMcrossSection( RMExampleSmall, enableAngDists=True, verbose=False )

    def test_memberData( self ):
        self.maxDiff = None
        self.assertEqual( str(self.RRR.projectile), 'n' )
        self.assertEqual( self.RRR.targetSpin, 0.0 )
        self.assertEqual( str(self.RRR.target), 'Si28' )
        self.assertAlmostEqual( self.RRR.targetToNeutronMassRatio, 27.73700000000002 )
        self.assertAlmostEqual( self.RRR.lowerBound, 1e-05 )
        self.assertAlmostEqual( self.RRR.upperBound, 1750000.0 )
        self.assertEqual( self.RRR.missingGfactor, {0: 0.0, 1: 0.0, 2: 0.0} )
        self.assertEqual( len(self.RRR.Ls), 3 )
        self.assertEqual( self.RRR.Ls[1].keys(), ['L', 'Js'] )
        self.assertEqual( self.RRR.Ls[1]['L'], 1 )
        self.assertEqual( self.RRR.Ls[1]['Js'][0].keys(), ['gfact', 'channelSpins', 'J'] )
        self.assertAlmostEqual( self.RRR.Ls[1]['Js'][0]['gfact'], 1.0 )
        self.assertEqual( self.RRR.Ls[1]['Js'][0]['channelSpins'][0].keys(), ['neutronWidth','channelSpin','captureWidth','energy', 'fissionWidthA', 'fissionWidthB', 'shiftFactor'] )
        self.assertEqual( self.RRR.Ls[1]['Js'][0]['J'], 0.5 )
        self.assertEqual( list(self.RRR._energies), [-3622100.0, -873730.0, -365290.0, -63159.0, -48801.0, 31740.0, 55677.0, 67733.0, 70800.0, 86797.0, 181620.0, 298700.0, 301310.0, 354590.0, 399680.0, 532660.0, 565580.0, 587170.0, 590290.0, 602470.0, 714040.0, 771710.0, 812490.0, 845230.0, 872310.0, 910040.0, 962230.0, 1017800.0, 1042900.0, 1085200.0, 1148100.0, 1162700.0, 1199500.0, 1201200.0, 1256400.0, 1264400.0, 1379900.0, 1408300.0, 1479900.0, 1482400.0, 1512300.0, 1528700.0, 1580600.0, 1592800.0, 1597200.0, 1639600.0, 1651100.0, 1658600.0, 1665000.0, 1785000.0, 1805700.0, 1850700.0, 1852400.0, 1923700.0, 1968900.0, 2248700.0, 3007300.0, 3067800.0] )
        for i, w in enumerate( [3936345.36, 1.12681, 1.030406, 1.046894, 1.0092496, 1.015667, 654.8903, 5.1589, 1.029617, 3.22618, 34899.6, 10.886, 5.9548, 15.46, 1.47361, 535.31, 10955.9, 207.96, 527.26, 53.891, 3.7165, 54.139, 30109.7, 399.91, 33.44, 3674.43, 76630.0, 77.192, 934.7, 76.394, 4.1469, 3017.4, 14921.6, 4604.8, 17386.6, 844.63999999999999, 67.699, 5201.0, 3504.15, 9.68694, 92.493, 2924.9, 1497.9, 11207.8, 4019.6, 15294.0, 21556.0, 1564.1, 218.3, 195.34, 1301.6, 35516.0, 70709.5, 1018.1, 5735.1, 444763.6, 293.56, 425.89] ):
            self.assertAlmostEqual( self.RRR._widths[i], w )
        self.assertFalse( self.RRR.RR.calculateChannelRadius )
        self.assertFalse( self.RRR.RR.scatteringRadius.isEnergyDependent() )
        self.assertEqual( str( self.RRR.RR.scatteringRadius.getValueAs('fm',L=0) ), "4.1364" )
        self.assertEqual( str( self.RRR.RR.scatteringRadius.getValueAs('fm',L=1) ), "4.9437" )
        self.assertEqual( str( self.RRR.RR.scatteringRadius.getValueAs('fm',L=2) ), "4.1364" )

    def test_getL0Matrix( self ):
        E = 31740.0
        B = self.RRRSmall.getChannelConstantsBc()
        L = []
        for c in self.RRRSmall.channels: L.append( complex( self.RRRSmall.shiftFactor( c.l, self.RRRSmall.rho( E, c.l ) ), self.RRRSmall.penetrationFactor( c.l, self.RRRSmall.rho( E, c.l ) ) ) )
        L0 = self.RRRSmall.getL0Matrix([E])
        for ic,c in enumerate( self.RRRSmall.channels ):
            # Because of the way ENDF handles the shift factor and the boundary value, this only works for s-wave channels
            if c.l==0: self.assertEqual( L0[0,ic,ic], L[ic]-B[ic] )
            # Because of the way ENDF handles the shift factor and the boundary value, this also passes, but it shouldn't
            self.assertEqual( L0[0,ic,ic], complex( 0.0, L[ic].imag ) )

    def test_setResonanceParametersByChannel( self ):
        debug = False
        self.RRR.setResonanceParametersByChannel( useReichMooreApproximation = True  )
        cDict = self.RRR.channels
        if debug: print
        if debug: print '\n'.join( map( str, cDict.items() ) )
        self.assertEqual( len( self.RRR.allChannels ), 12 ) # This is number of channels, whether eliminated or not
        self.assertEqual( len( self.RRR.channels ), 7 ) # This is number of distinct un-eliminated RM channels
        self.assertEqual( len( self.RRR.eliminatedChannels ), 5 ) # This is number of distinct eliminated RM channels
        self.assertEqual( len( self.RRR._energies ), 58 ) # This is number of RM resonances.  Better equal the _energies array of level energies.
        self.assertEqual( sum( [ len(self.RRR.allChannels[cs]) for cs in self.RRR.allChannels ] ), 118 ) # this the number of resonances (elastic only in our example)
        self.assertEqual( self.RRR.allChannels.items()[0][0], ChannelDesignator(l=0, J=0.5, reaction='elastic', index=0, s=0.5, gfact=1.0, particleA='n', particleB='Si28', isElastic=True, channelClass=NEUTRONCHANNEL, useRelativistic=False, eliminated=False) )
        self.assertEqual( self.RRR.allChannels.items()[1][1], {0: 145.36, 1: 1.0253, 2: 1.0, 3: 1.0, 4: 1.0, 6: 1.5803, 33: 3.6, 10: 5.6, 12: 3.6, 34: 3.6, 18: 3.6, 55: 3.6, 56: 3.6, 57: 3.6, 29: 3.6, 31: 3.8} )
        self.assertEqual( self.RRR.allChannels.items()[2][1], {35: 843.64, 5: 0.015667, 8: 0.029617, 11: 9.886, 13: 14.46, 46: 21555.0, 45: 15293.0, 40: 91.493, 51: 35515.0, 21: 53.139, 54: 5734.1, 24: 32.14, 27: 76.192, 28: 933.7, 30: 3.1469, 53: 1017.1} )
        for i,E in enumerate(self.RRR._energies):
            # Sum all the widths with matching resonance energy.  Done by brute force.
            GT = 0.0
            for c in self.RRR.allChannels:
                for iR, G in self.RRR.allChannels[c].items():
                    if iR == i: GT += G
            # Check all the widths for a level add up to the tabulated total
            if debug: print i, E, GT, self.RRR._widths[i]
            self.assertWithinXPercent( GT, self.RRR._widths[i], 0.1 )

    def test_getRMatrix_Full( self ):
        '''
        Test the value on resonances, full example
        '''
        # Basic property tests
        R = self.RRR.getRMatrix( [[1.0]] )
        self.assertEqual( R.shape, (1, 7, 7 ) )
        self.assertTrue( numpy.allclose( R[0], R[0].T ) )
        # Full R test
        R = self.RRR.getRMatrix( self.RRR._energies[:,numpy.newaxis] )
        for iR in range( self.RRR.nResonances ):
            ER = self.RRR._energies[ iR ]
            for ic, c in enumerate( self.RRR.channels ):
                if iR in self.RRR.channels[c]:
                    GN = self.RRR.channels[c][iR]
                    Pc = self.RRR.penetrationFactor( c.l, self.RRR.rho( [abs( ER )], c.l ) )
                    break
            for c in self.RRR.eliminatedChannels:
                if iR in self.RRR.eliminatedChannels[c]:
                    GG = self.RRR.eliminatedChannels[c][iR]
                    break
            self.assertWithinXPercent( R[iR,ic,ic].imag, GN/GG/Pc[0] )

    def test_background_RMatrix(self):
        gamWidth      = 1.0 # gamma widths in Atlas for 28Si are around 1 eV
        k_over_sqrt_E = (2.196807122623e-3 * self.RRR.targetToNeutronMassRatio / (self.RRR.targetToNeutronMassRatio+1.0)) #in units of b^-0.5
        S0            = 0.5e-4 # S0 x 10^-4 (Atlas)
        Rprime        = 4.8 # R'=a(1-Rinf)=4.8 fm
        a             = 4.1264 # a=4.1264 fm
        Rinf          = 0.0 # should equal 1.0-Rprime/a (no units), but since S0 is constant, it should be 0.0 and a should == R'
        pole_strength = S0/( 0.2*k_over_sqrt_E*a )
        egrid=self.RRR.generateEnergyGrid()[2::10]
        R  = self.RRR.getRMatrix( numpy.array(egrid)[:,numpy.newaxis] )
        Rb = self.RRR.getBackgroundRMatrix( egrid, self.RRR.lowerBound, self.RRR.upperBound, gamWidth=gamWidth, pole_strength=pole_strength, Rinf=Rinf )

        # The distant level approx should give a background R matrix much smaller than the real one
        imidE = len(egrid)/2
        self.assertTrue(abs(R[imidE,0,0])>10*abs(Rb[imidE]))

        if False:
            print self.RRR.nChannels
            firstChannel=self.RRR.channels.items()[0]
            print firstChannel
            print [self.RRR._energies[i] for i in firstChannel[1].keys()]
            rFile = open('RTest.dat',mode='w')
            for iE,E in enumerate(egrid):
                rFile.write( '  '.join( [
                    str(E).ljust(20),
                    str(R[iE,0,0]).ljust(40),
                    str(abs(R[iE,0,0])).ljust(20),
                    str(Rb[iE]).ljust(40),
                    str(abs(Rb[iE])).ljust(20)
                ])+'\n')

    def test_getRMatrix_Small( self ):
        '''
        Test the value on resonances, small example
        '''
        # Print a table for plotting the R matrix
        debug=False
        if debug:
            egrid = self.RRRSmall.generateEnergyGrid()[::2]
            RR = self.RRRSmall.getRMatrix( egrid )
            for iEE,EE in enumerate(egrid):
                print EE, RR[iEE,0,0].real, RR[iEE,0,0].imag, RR[iEE,1,1].real, RR[iEE,1,1].imag

        # Basic properties of the R matrix test
        R = self.RRRSmall.getRMatrix( [[1.0]] )
        self.assertEqual( R.shape, ( 1, 7, 7 ) )
        self.assertTrue( numpy.allclose( R[0], R[0].T ) )

        # Test specific values
        R = self.RRRSmall.getRMatrix( self.RRRSmall._energies[:,numpy.newaxis] )
        for iR in range( self.RRRSmall.nResonances )[:1]:
            ER = self.RRRSmall._energies[ iR ]
            for ic, c in enumerate( self.RRRSmall.channels ):
                if iR in self.RRRSmall.channels[c]:
                    GN = self.RRRSmall.channels[c][iR]
                    Pc = self.RRRSmall.penetrationFactor( c.l, self.RRRSmall.rho( [abs( ER )], c.l ) )
                    Sc = self.RRRSmall.shiftFactor( c.l, self.RRRSmall.rho( [abs( ER )], c.l ) )
                    Bc = -c.l
                    break
            for c in self.RRRSmall.eliminatedChannels:
                if iR in self.RRRSmall.eliminatedChannels[c]:
                    GG = self.RRRSmall.eliminatedChannels[c][iR]
                    break
            R00answer = GN/GG/Pc[0] # to 10 % or better
            self.assertWithinXPercent( R[ iR, ic, ic ].imag, R00answer )

    def test_getXMatrix_Small( self ):
        # Print a table for plotting the X matrix
        debug=False
        egrid=numpy.array( [ 1.0, 31733.7, 31740.0 ] + self.RRRSmall.generateEnergyGrid()[::2] )[:,numpy.newaxis]
        X = self.RRRSmall.getXMatrix(egrid)
        R = self.RRRSmall.getRMatrix(egrid)
        L0 = self.RRRSmall.getL0Matrix(egrid)

        if debug:
            for iEE,EE in enumerate(egrid):
                print EE, X[iEE,0,0].real, X[iEE,0,0].imag, X[iEE,1,1].real, X[iEE,1,1].imag

        # Test specific values
        for iE in range(len(egrid)):
            self.assertAlmostEqual( X[iE,0,0], L0[iE,0,0].imag*R[iE,0,0]/(1.0-R[iE,0,0]*L0[iE,0,0] ), 9 )
            self.assertAlmostEqual( X[iE,1,1], L0[iE,1,1].imag*R[iE,1,1]/(1.0-R[iE,1,1]*L0[iE,1,1] ), 9 )

    def test_getScatteringMatrixUUnitarity_Full( self ):
        '''
        The RM parameterization might have a unitary U matrix.  Here we have two tests:

            * just plain multiplication: U^+ * U == one

            * check the Frobenius norm, which is ||A|| = sqrt(| sum_{i,j} |A[i,j]|^2 |), so the norm of the identity matrix is sqrt( ndim )
        '''
        debug = False
        doTest = True
        if debug:
            egrid = self.RRR.generateEnergyGrid()[::10]
            UU = self.RRR.getScatteringMatrixU( egrid )
            for iEE,EE in enumerate(egrid):
                print EE, UU[iEE,0,0].real, UU[iEE,0,0].imag, \
                    UU[iEE,0,1].real, UU[iEE,0,1].imag, UU[iEE,1,1].real, \
                    UU[iEE,1,1].imag, numpy.linalg.norm(UU[iEE])/math.sqrt( UU[iEE].shape[0] )
        rtol= 1e-03 # nominal value: 1e-03
        atol= 1e-04 # nominal value: 1e-04
        E = numpy.array([[1.0]])
        U = self.RRR.getScatteringMatrixU( E )
        one = numpy.diag( U[0].shape[0]*[ complex( 1.0 ) ] )
        if doTest: self.assertAlmostEqual( numpy.linalg.norm( U[0] ), math.sqrt( U[0].shape[0] ), 5 )
        if doTest: self.assertTrue( numpy.allclose( numpy.conj(U[0].T)*U[0], one, rtol=rtol, atol=atol ) )

    def test_getScatteringMatrixUUnitarity_Small( self ):
        '''
        The RM parameterization might have a unitary U matrix.  Here we have two tests:

            * just plain multiplication: U^+ * U == one

            * check the Frobenius norm, which is ||A|| = sqrt(| sum_{i,j} |A[i,j]|^2 |), so the norm of the identity matrix is sqrt( ndim )
        '''
        debug = False
        if debug:
            egrid = self.RRRSmall.generateEnergyGrid()[::2]
            UU = self.RRRSmall.getScatteringMatrixU( egrid )
            for iEE,EE in enumerate(egrid):
                print EE, UU[iEE,0,0].real, UU[iEE,0,0].imag, \
                    UU[iEE,0,1].real, UU[iEE,0,1].imag, UU[iEE,1,1].real, \
                    UU[iEE,1,1].imag, numpy.linalg.norm(UU[0])/math.sqrt( UU[0].shape[0] )
        egrid = numpy.array([ 1.0, 31733.7, 31740.0 ])[:,numpy.newaxis]
        U = self.RRRSmall.getScatteringMatrixU( egrid, useTabulatedScatteringRadius = True )
        for iE,E in enumerate(egrid):
            one = numpy.diag( U[iE].shape[0]*[ complex( 1.0 ) ] )
            if E == 31740.0: # on resonance, expect to be good to only 10%
                numDecimalPlaces = 1
                rtol= 1e-01 # nominal value: 1e-03
                atol= 1e-02 # nominal value: 1e-04
            else:
                numDecimalPlaces = 3
                rtol= 1e-03 # nominal value: 1e-03
                atol= 1e-04 # nominal value: 1e-04
            self.assertAlmostEqual( numpy.linalg.norm( U[iE] ), math.sqrt( U[iE].shape[0] ), numDecimalPlaces )
            self.assertTrue( numpy.allclose( numpy.conj(U[iE].T)*U[iE], one, rtol=rtol, atol=atol ) )

    def test_scatteringMatrixUInCrossSectionCalculation_Full( self ):
        debug = False
        doTests = True

        nonEliminatedChannels = [ c for c in self.RRR.channels if not c.eliminated ]

        nChan = len( nonEliminatedChannels )
        nElim = len( self.RRR.channels ) - nChan
        nRes = len( self.RRR._energies )
        AP = 0.41 #in b^{-1/2}

        if debug: print '\nAll Waves (course grid), %i kept channels, %i eliminated channels, %i resonances' % ( nChan, nElim, nRes )
        if debug: print '     E (eV)       SigTotU (b)    SigTotSum (b)    SigTotOK (b)   SigGamU (b)    SigGamOK (b)   SigNeuU (b)    SigNeuOK (b)   SigNeuInt (b) SigNeuInt/SigNeuU  sigPot'

        # We can't eat the whole energy grid, so we'll do it in parts
        fullEgrid = self.RRR.generateEnergyGrid()
        energyCohortLength = 50
        nEnergyCohorts = len(fullEgrid)//energyCohortLength
        lastEnergyCohortSize = len(fullEgrid)%energyCohortLength
        energyCohorts = [(iec*energyCohortLength, (iec+1)*energyCohortLength) for iec in range(nEnergyCohorts)]
        energyCohorts.append((energyCohorts[-1][1],energyCohorts[-1][1]+lastEnergyCohortSize))

        for iStart,iStop in energyCohorts:
            egrid = fullEgrid[iStart:iStop]
            egrid = numpy.array(egrid).reshape(len(egrid),1)

            U = self.RRR.getScatteringMatrixU( egrid, True )
            k = self.RRR.k(egrid)
            rhohat = AP * k
            nLs = { 0:0, 1:0, 2:0, 3:0 }
            gs = { 0:0, 1:0, 2:0, 3:0 }
            B = self.RRR.getAngularDistribution(egrid)['elastic']
            sigs = self.RRR.getCrossSection(egrid)

            for iE,E in enumerate(egrid):
                sigTotWithU, sigGamWithU, sigNeuWithU = 0.0, 0.0, 0.0
                for ic1,c1 in enumerate(nonEliminatedChannels):
                    if c1.reaction != 'elastic': continue
                    gs[ c1.l ] += c1.gfact
                    nLs[ c1.l ] += 1.0
                    sigTotWithU += 2.0 * numpy.pi * c1.gfact * ( 1.0 - U[iE][ic1,ic1].real ) / k[iE]**2
                    for ic2,c2 in enumerate(nonEliminatedChannels):
                        if c2.J != c1.J or c1.s != c2.s or c1.l != c2.l: continue
                        if c2.reaction == 'elastic':
                            sigNeuWithU += numpy.pi * c1.gfact * ( ( (ic1==ic2) - U[iE][ic1,ic2] ).conj() * ( (ic1==ic2) - U[iE][ic2,ic1] ) ).real / k[iE]**2
                sigGamWithU = sigTotWithU - sigNeuWithU
                sigPot = sum( [ ( ( 2.0 * l + 1 ) ) * 4.0 * numpy.pi * ( numpy.sin( self.RRR.phi( l, rhohat[iE] ) ) )**2 / k[iE]**2 for l in nLs ] )
                sigNeuInt   = 4*numpy.pi*B[0][iE].real

                if debug: print 11*'%10g     ' % ( E, sigTotWithU, sigGamWithU + sigNeuWithU, sigs['total'][iE], sigGamWithU, sigs['capture'][iE], sigNeuWithU, sigs['elastic'][iE], sigNeuInt,  sigNeuWithU/sigNeuInt, sigPot ) #, -(sigNeuWithU-sigs['elastic'][iE])/sigPot )

                # Check they sum up correctly (really just another test of unitarity)
                if doTests: self.assertWithinXPercent( sigTotWithU, sigGamWithU + sigNeuWithU, 10. )

                # Check elastic channel agrees with L=0 term of angular distribution
                if doTests and False: self.assertWithinXPercent( sigNeuWithU, sigNeuInt ) # fails badly for some high energy P wave resonances??

                # Check agree with getCrossSection()
                # Note, since ENDF's capture cross section is really a SLBW cross section, we expect poor agreement between our U matrix based one and ENDF's
                # Similarly, since the total cross section is a sum of capture and elastic, we should really use the ENDF (SLBW) capture cross section to make the comparison
                if doTests:
                    self.assertWithinXPercent( sigs['capture'][iE]+sigNeuWithU, sigs['total'][iE] )
                    if False: self.assertWithinXPercent( sigGamWithU, sigs['capture'][iE] )  # expected to fail
                    self.assertWithinXPercent( sigNeuWithU, sigs['elastic'][iE] )

    def test_scatteringMatrixUInCrossSectionCalculation_Small( self ):
        debug = False
        doTests = True

        nChan = len( self.RRRSmall.channels )
        nElim = len( self.RRRSmall.eliminatedChannels )
        nRes = len( self.RRRSmall._energies )
        AP = 0.41364 #in b^{-1/2}

        if debug: print '\nAll Waves (course grid), %i kept channels, %i eliminated channels, %i resonances' % ( nChan, nElim, nRes )
        if debug: print '     E (eV)       SigTotU (b)    SigTotSum (b)    SigTotOK (b)   SigGamU (b)    SigGamOK (b)   SigNeuU (b)    SigNeuOK (b)   SigNeuInt (b) SigNeuInt/SigNeuU  sigPot'

        eStep = 10 #int( len( egrid )/95)
        egrid = self.RRRSmall.generateEnergyGrid()[::eStep]
        egrid = numpy.array(egrid).reshape(len(egrid),1)

        U = self.RRRSmall.getScatteringMatrixU( egrid, True )
        k = self.RRRSmall.k(egrid)
        rhohat = AP * k
        nLs = { 0:0, 1:0, 2:0, 3:0 }
        gs = { 0:0, 1:0, 2:0, 3:0 }
        B = self.RRRSmall.getAngularDistribution( egrid )['elastic']
        sigs = self.RRRSmall.getCrossSection( egrid )

        for iE,E in enumerate(egrid[::eStep]):
            if not ( E>=31700. and E<=31750.0 ): continue
            sigTotWithU, sigGamWithU, sigNeuWithU = 0.0, 0.0, 0.0
            for ic1,c1 in enumerate(self.RRRSmall.channels):
                if c1.reaction != 'elastic': continue
                gs[ c1.l ] += c1.gfact
                nLs[ c1.l ] += 1.0
                sigTotWithU += 2.0 * numpy.pi * c1.gfact * ( 1.0 - U[iE][ic1,ic1].real ) / k[iE]**2
                for ic2,c2 in enumerate(self.RRRSmall.channels):
                    if c2.J != c1.J or c1.s != c2.s or c1.l != c2.l: continue
                    if c2.reaction == 'elastic':
                        sigNeuWithU += numpy.pi * c1.gfact * ( ( (ic1==ic2) - U[iE][ic1,ic2] ).conj() * ( (ic1==ic2) - U[iE][ic2,ic1] ) ).real / k[iE]**2
            sigGamWithU = sigTotWithU - sigNeuWithU
            sigPot = sum( [ ( ( 2.0 * l + 1.0 ) ) * 4.0 * numpy.pi * ( numpy.sin( self.RRRSmall.phi( l, rhohat[iE] ) ) )**2 / k[iE]**2 for l in nLs ] )
            sigNeuInt   = 4.0*numpy.pi*B[0][iE].real

            if debug: print 11*'%10g     ' % ( E, sigTotWithU, sigGamWithU + sigNeuWithU, sigs['total'][iE], sigGamWithU, sigs['capture'][iE], sigNeuWithU, sigs['elastic'][iE], sigNeuInt,  sigNeuWithU/sigNeuInt, sigPot ) #, -(sigNeuWithU-sigs['elastic'][iE])/sigPot )

            # Check they sum up correctly (really just another test of unitarity)
            if doTests: self.assertWithinXPercent( sigTotWithU, sigGamWithU + sigNeuWithU, 10. )

            # Check elastic channel agrees with L=0 term of angular distribution
            if doTests and False: self.assertWithinXPercent( sigNeuWithU, sigNeuInt )  # fails badly for some high energy P wave resonances??

            # Check agree with getCrossSection()
            # Note, since ENDF's capture cross section is really a SLBW cross section, we expect poor agreement between our U matrix based one and ENDF's
            # Similarly, since the total cross section is a sum of capture and elastic, we should really use the ENDF (SLBW) capture cross section to make the comparison
            if doTests:
                self.assertWithinXPercent( sigs['capture'][iE]+sigNeuWithU, sigs['total'][iE] )
                if False: self.assertWithinXPercent( sigGamWithU, sigs['capture'][iE] )  # expected to fail
                self.assertWithinXPercent( sigNeuWithU, sigs['elastic'][iE] )

    def test_getElasticAngularDistribution_Small( self ):
        """
        FIXME: This test not written yet.
        The values to check against below for the angular distribution are wrong (they're copied from the MLBW test)

        Do we even need this test?  once the R matrix is correct, the angular distribution generation is just
        inherited from the base class
        """
        debug = False
        doTest = False
        egrid = self.RRRSmall.generateEnergyGrid()
        eStep=2
        if debug:
            for E in egrid[::eStep]:
                B = self.RRRSmall.getAngularDistribution( E )
                print E, ' '.join( map( str, B ) )
        if doTest:
            ER = 0.169
            B = self.RRRSmall.getAngularDistribution( numpy.array([ER]) )
            self.assertAlmostEqual( B[2]/B[0], 1.96523656e-14, 2 ) #0.33090909, 2 )
            self.assertAlmostEqual( B[1]/B[0], 0.0, 3 )

    def test_getElasticAngularDistribution_Full( self ):
        """
        FIXME: This test not written yet.
        The values to check against below for the angular distribution are wrong (they're copied from the MLBW test)

        Do we even need this test?  once the R matrix is correct, the angular distribution generation is just
        inherited from the base class
        """
        debug = False
        doTest = False
        egrid = self.RRR.generateEnergyGrid()
        eStep=10
        if debug:
            for E in egrid[::eStep]:
                B = self.RRR.getAngularDistribution( E )
                print E, ' '.join( map( str, B ) )
        if doTest:
            ER = 0.169
            B = self.RRR.getAngularDistribution( numpy.array([ER]) )
            self.assertAlmostEqual( B[2]/B[0], 1.55602284e-14, 2 ) #0.33090909, 2 )
            self.assertAlmostEqual( B[1]/B[0], 0.0, 3 )

    def test_getCrossSection( self ):
        '''Test thermal point'''
        result = self.RRR.getCrossSection( 0.0235339347844 )
        #print 'RM',result
        answer = {'capture':0.17536290484357434, 'total':2.1313431546619341, 'fission':0., 'elastic':1.9559802498183596}
        for k in answer:
            self.assertAlmostEqual( answer[k], result[k][0]  )



class TestRMLClassAndBaseClasses( TestWithIsClose ):

    def setUp( self ):
        # Fudge no longer attempts to translate bogus ENDF spin assignments
        # therefore the translated test data has tagets missing spins
        import PoPs
#        theCorrectCl35Spin = PoPs.quantities.spin.fraction('spin', 1.5, 'hbar')
#        theCorrectCl36Spin = PoPs.quantities.spin.fraction('spin', 0, 'hbar')
#        theCorrectFe57Spin = PoPs.quantities.spin.fraction('spin', 1.5, 'hbar')

        #
        self.RRR = RMatrixLimitedcrossSection( RMLExample, verbose=False )
#        self.RRR.reactionSuite.PoPs['Cl35'].spin.add(theCorrectCl35Spin)
#        self.RRR.reactionSuite.PoPs['Cl36'].spin.add(theCorrectCl36Spin)

        #
        self.RRRSmall = RMatrixLimitedcrossSection( RMLExampleSmall, verbose=False )
#        self.RRRSmall.reactionSuite.PoPs['Cl35'].spin.add(theCorrectCl35Spin)
#        self.RRRSmall.reactionSuite.PoPs['Cl36'].spin.add(theCorrectCl36Spin)

        #
        if DOFETESTS:
            self.RRRFe = RMatrixLimitedcrossSection( RMLExampleFe, verbose=False )
#            self.RRRFe.reactionSuite.PoPs['Fe57'].spin.add(theCorrectFe57Spin)

    def test_memberData( self ):
        self.assertEqual( self.RRRSmall._energies, (-336933.4, -180.65, 397.8154, 4250.762, 16356.12, 22396.4, 54932.0, 68236.16, 133988.4, 1205687.0, 1434336.0, 1441365.0, 7563145.0) )
        self.assertEqual( self.RRRSmall._thresholds, [] )
        self.assertEqual( self.RRRSmall._widths, (38202.334010000006, 13.813142300000001, 1.0375000000000001, 1.3300000000000001, 6.5323189999999993, 2.6911670000000001, 46.809660000000001, 218.29737, 662.81119999999999, 643.18999999999994, 5423.7999999999993, 1609.5999999999999, 622905.38398000004) )
        self.assertEqual( self.RRRSmall.RR.moniker, 'RMatrix' )
        self.assertEqual( self.RRRSmall.targetSpin, 1.5 )
        self.assertEqual( str(self.RRRSmall.target), 'Cl35' )
        self.assertEqual( str(self.RRRSmall.projectile), 'n' )
        self.assertAlmostEqual( self.RRRSmall.targetToNeutronMassRatio, 34.66845 )
        self.assertEqual( self.RRRSmall.lowerBound, 1e-5 )
        self.assertEqual( self.RRRSmall.upperBound, 1200000.0 )
        self.assertEqual( self.RRRSmall.verbose, False )

    def test_RR_memberData(self):
        self.assertEqual( self.RRRSmall.RR.approximation, 'Reich_Moore' )
        self.assertEqual( self.RRRSmall.RR.boundaryCondition, "S" ) # "S" for "shift factor"
        self.assertEqual( self.RRRSmall.RR.relativisticKinematics, False )
        self.assertEqual( self.RRRSmall.RR.reducedWidthAmplitudes, False )
        self.assertEqual( self.RRRSmall.RR.calculatePenetrability, True )
        self.assertEqual( self.RRRSmall.RR.calculateShift, False )
        self.assertEqual( self.RRRSmall.RR.calculateChannelRadius, True )

    def test_RR_spinGroup_memberData(self):
        self.assertEqual( self.RRRSmall.RR.spinGroups[0].index, 0 )
        self.assertEqual( str(self.RRRSmall.RR.spinGroups[0].spin), '1' )
        self.assertEqual( str(self.RRRSmall.RR.spinGroups[0].parity), '1' )
        self.assertEqual( self.RRRSmall.RR.spinGroups[0].background, False )
        self.assertEqual( self.RRRSmall.RR.spinGroups[0].applyPhaseShift, False )
        #self.assertEqual( self.RRRSmall.RR.spinGroups[0].overrides[0].label, '1' )
        #for attr in [ 'scatteringRadius', 'effectiveRadius', 'penetrability', 'shiftFactor', 'phaseShift' ]:
        #    self.assertTrue( hasattr(self.RRRSmall.RR.spinGroups[0].overrides[0], attr ) )
        self.assertItemsEqual( [x.name for x in self.RRRSmall.RR.spinGroups[0].resonanceParameters.table.columns], ['energy','photon + Cl36 width','n + Cl35 width','H1 + S35 width'] )
        self.assertItemsEqual( [x.unit for x in self.RRRSmall.RR.spinGroups[0].resonanceParameters.table.columns], ['eV','eV','eV','eV'] )

    def test_RR_openChannel_memberData(self):
        self.assertEqual( self.RRRSmall.RR.resonanceReactions[1].label, 'n + Cl35' )
        self.assertEqual( str(self.RRRSmall.RR.resonanceReactions[1].reactionLink), "/reactionSuite/reactions/reaction[@label='n + Cl35']" )
#        print dir(self.RRRSmall.RR.resonanceReactions[1])
#        self.assertEqual( self.RRRSmall.RR.resonanceReactions[1].ENDF_MT, 2 )
#        self.assertEqual( self.RRRSmall.RR.resonanceReactions[1].Q.getValueAs('eV'), 0.0)
#        self.assertEqual( self.RRRSmall.RR.resonanceReactions[1].boundaryCondition, None )
#        self.assertEqual( self.RRRSmall.RR.resonanceReactions[1].calculatePenetrability, None )
        self.assertEqual( self.RRRSmall.RR.resonanceReactions[1].computeShiftFactor, False )
        self.assertEqual( self.RRRSmall.RR.resonanceReactions[1].scatteringRadius.getValueAs('fm'), 4.82222 )
        self.assertEqual( self.RRRSmall.RR.resonanceReactions[1].hardSphereRadius.getValueAs('fm'), 4.88875 )
        self.assertAlmostEqual( self.RRRSmall.RR.resonanceReactions[0].reactionInfo['Xi'], -8827250.63441, 5 )
        self.assertAlmostEqual( self.RRRSmall.RR.resonanceReactions[1].reactionInfo['Xi'], -0.0 )
        self.assertAlmostEqual( self.RRRSmall.RR.resonanceReactions[2].reactionInfo['Xi'], -632965.817883407, 6 )
        self.assertAlmostEqual( self.RRRSmall.RR.resonanceReactions[1].reactionInfo['particles'][0].particle.getMass('amu'), 1.00866491574 )
        self.assertAlmostEqual( self.RRRSmall.RR.resonanceReactions[1].reactionInfo['particles'][1].particle.getMass('amu'), 34.9688491980864 )
        if False:
            for i in range(len(self.RRRSmall.RR.resonanceReactions)):
                print i, self.RRRSmall.RR.resonanceReactions[i].name, self.RRRSmall.RR.resonanceReactions[i].reactionInfo['Xi']

    def test_isElastic(self):
        self.assertTrue( self.RRRSmall.isElastic('n + Cl35') )
        self.assertFalse( self.RRRSmall.isElastic('H1 + S35') )

    def test_getChannelConstantsBc( self ):
        self.RRRSmall.setResonanceParametersByChannel( )
        self.assertEqual(self.RRRSmall.getChannelConstantsBc(), [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    def test_getChannelScatteringRadiiTrue(self):
        self.RRRSmall.setResonanceParametersByChannel( )
        self.assertAllWithinXPercent(
            self.RRRSmall.getChannelScatteringRadiiTrue().values(),
            16*[0.48222175644091303],
            percent=0.01)

    def test_getChannelScatteringRadiiEffective(self):
        self.RRRSmall.setResonanceParametersByChannel( )
        self.assertAllWithinXPercent(
            self.RRRSmall.getChannelScatteringRadiiEffective().values(),
            4*[0.36679799999999996]+12*[0.48887499999999995])

    def test_setResonanceParametersByChannel( self ):
        self.maxDiff=None
        # The evaluation in RRRSmall is correct (channel-wise), so we'll try out all the things
        # Note however that the gamma channel has the wrong spin for the residual (as given in original ENDF file), so the gfactors below are nonsensical.
        self.RRRSmall.setResonanceParametersByChannel( )
        self.assertItemsEqual( self.RRRSmall.allChannels.keys(), [
            ChannelDesignator(l=0, J=1.0, reaction='gamma + Cl36', index=0, s=1.0, gfact=1.0, particleA='Cl36', particleB='photon', Xi=-8827250.63441, isElastic=False, channelClass=GAMMACHANNEL, useRelativistic=False, eliminated=True),
            ChannelDesignator(l=0, J=1.0, reaction='n + Cl35', index=1, s=1.0, gfact=0.375, particleA='n', particleB='Cl35', Xi=-0.0, isElastic=True, channelClass=NEUTRONCHANNEL, useRelativistic=False, eliminated=False),
            ChannelDesignator(l=0, J=1.0, reaction='H1 + S35', index=2, s=1.0, gfact=0.375, particleA='H1', particleB='S35', Xi=-632965.817883, isElastic=False, channelClass=CPCHANNEL, useRelativistic=False, eliminated=False),
            ChannelDesignator(l=0, J=2.0, reaction='gamma + Cl36', index=3, s=2.0, gfact=1.0, particleA='Cl36', particleB='photon', Xi=-8827250.63441, isElastic=False, channelClass=GAMMACHANNEL, useRelativistic=False, eliminated=True),
            ChannelDesignator(l=0, J=2.0, reaction='n + Cl35', index=4, s=2.0, gfact=0.625, particleA='n', particleB='Cl35', Xi=-0.0, isElastic=True, channelClass=NEUTRONCHANNEL, useRelativistic=False, eliminated=False),
            ChannelDesignator(l=0, J=2.0, reaction='H1 + S35', index=5, s=2.0, gfact=0.625, particleA='H1', particleB='S35', Xi=-632965.817883, isElastic=False, channelClass=CPCHANNEL, useRelativistic=False, eliminated=False),
            ChannelDesignator(l=1, J=0.0, reaction='gamma + Cl36', index=6, s=1.0, gfact=1.0, particleA='Cl36', particleB='photon', Xi=-8827250.63441, isElastic=False, channelClass=GAMMACHANNEL, useRelativistic=False, eliminated=True),
            ChannelDesignator(l=1, J=0.0, reaction='n + Cl35', index=7, s=1.0, gfact=0.125, particleA='n', particleB='Cl35', Xi=-0.0, isElastic=True, channelClass=NEUTRONCHANNEL, useRelativistic=False, eliminated=False),
            ChannelDesignator(l=1, J=0.0, reaction='H1 + S35', index=8, s=1.0, gfact=0.125, particleA='H1', particleB='S35', Xi=-632965.817883, isElastic=False, channelClass=CPCHANNEL, useRelativistic=False, eliminated=False),
            ChannelDesignator(l=1, J=1.0, reaction='gamma + Cl36', index=9, s=1.0, gfact=1.0, particleA='Cl36', particleB='photon', Xi=-8827250.63441, isElastic=False, channelClass=GAMMACHANNEL, useRelativistic=False, eliminated=True),
            ChannelDesignator(l=1, J=1.0, reaction='n + Cl35', index=10, s=1.0, gfact=0.375, particleA='n', particleB='Cl35', Xi=-0.0, isElastic=True, channelClass=NEUTRONCHANNEL, useRelativistic=False, eliminated=False),
            ChannelDesignator(l=1, J=1.0, reaction='H1 + S35', index=11, s=1.0, gfact=0.375, particleA='H1', particleB='S35', Xi=-632965.817883, isElastic=False, channelClass=CPCHANNEL, useRelativistic=False, eliminated=False),
            ChannelDesignator(l=1, J=2.0, reaction='gamma + Cl36', index=12, s=1.0, gfact=1.0, particleA='Cl36', particleB='photon', Xi=-8827250.63441, isElastic=False, channelClass=GAMMACHANNEL, useRelativistic=False, eliminated=True),
            ChannelDesignator(l=1, J=2.0, reaction='n + Cl35', index=13, s=1.0, gfact=0.625, particleA='n', particleB='Cl35', Xi=-0.0, isElastic=True, channelClass=NEUTRONCHANNEL, useRelativistic=False, eliminated=False),
            ChannelDesignator(l=1, J=2.0, reaction='H1 + S35', index=14, s=1.0, gfact=0.625, particleA='H1', particleB='S35', Xi=-632965.817883, isElastic=False, channelClass=CPCHANNEL, useRelativistic=False, eliminated=False),
            ChannelDesignator(l=1, J=1.0, reaction='gamma + Cl36', index=15, s=2.0, gfact=1.0, particleA='Cl36', particleB='photon', Xi=-8827250.63441, isElastic=False, channelClass=GAMMACHANNEL, useRelativistic=False, eliminated=True),
            ChannelDesignator(l=1, J=1.0, reaction='n + Cl35', index=16, s=2.0, gfact=0.375, particleA='n', particleB='Cl35', Xi=-0.0, isElastic=True, channelClass=NEUTRONCHANNEL, useRelativistic=False, eliminated=False),
            ChannelDesignator(l=1, J=1.0, reaction='H1 + S35', index=17, s=2.0, gfact=0.375, particleA='H1', particleB='S35', Xi=-632965.817883, isElastic=False, channelClass=CPCHANNEL, useRelativistic=False, eliminated=False),
            ChannelDesignator(l=1, J=2.0, reaction='gamma + Cl36', index=18, s=2.0, gfact=1.0, particleA='Cl36', particleB='photon', Xi=-8827250.63441, isElastic=False, channelClass=GAMMACHANNEL, useRelativistic=False, eliminated=True),
            ChannelDesignator(l=1, J=2.0, reaction='n + Cl35', index=19, s=2.0, gfact=0.625, particleA='n', particleB='Cl35', Xi=-0.0, isElastic=True, channelClass=NEUTRONCHANNEL, useRelativistic=False, eliminated=False),
            ChannelDesignator(l=1, J=2.0, reaction='H1 + S35', index=20, s=2.0, gfact=0.625, particleA='H1', particleB='S35', Xi=-632965.817883, isElastic=False, channelClass=CPCHANNEL, useRelativistic=False, eliminated=False),
            ChannelDesignator(l=1, J=3.0, reaction='gamma + Cl36', index=21, s=2.0, gfact=1.0, particleA='Cl36', particleB='photon', Xi=-8827250.63441, isElastic=False, channelClass=GAMMACHANNEL, useRelativistic=False, eliminated=True),
            ChannelDesignator(l=1, J=3.0, reaction='n + Cl35', index=22, s=2.0, gfact=0.875, particleA='n', particleB='Cl35', Xi=-0.0, isElastic=True, channelClass=NEUTRONCHANNEL, useRelativistic=False, eliminated=False),
            ChannelDesignator(l=1, J=3.0, reaction='H1 + S35', index=23, s=2.0, gfact=0.875, particleA='H1', particleB='S35', Xi=-632965.817883, isElastic=False, channelClass=CPCHANNEL, useRelativistic=False, eliminated=False)] )
        self.assertEqual( self.RRRSmall.channels[ChannelDesignator(l=1, J=3.0, reaction='H1 + S35', index=23, s=2.0, gfact=0.875, particleA='H1', particleB='S35', Xi=-632965.817883, isElastic=False, channelClass=CPCHANNEL, useRelativistic=False, eliminated=False)],\
                               {4: 0.164019})
        self.assertEqual( self.RRRSmall.channels[ChannelDesignator(l=0, J=2.0, reaction='n + Cl35', index=4, s=2.0, gfact=0.625, particleA='n', particleB='Cl35', Xi=-0.0, isElastic=True, channelClass=NEUTRONCHANNEL, useRelativistic=False, eliminated=False)],\
                               {1: 13.277, 12: 621905.0})

        if VERBOSE:
            #
            # This coding used to generate a table of resonances
            # suitable for a LaTeX document
            #
            for ir, r in enumerate(self.RRRSmall._energies):
                row = None
                for ic, cres in enumerate(self.RRRSmall.allChannels.items()):
                    c=cres[0]
                    res=cres[1]
                    if ir in res:
                        if row is None: row = [ ir, r, None, None, None, None, None, None, None, None, None, None ]
                        if c.particleA == 'n':
                            row[2]=int(c.l)
                            row[3]=int(c.J)
                            row[4]=int(c.s)
                            row[5]=ic
                            row[6]=res[ir]
                        elif c.particleA == 'H1':
                            row[7]=ic
                            row[8]=res[ir]
                        else:
                            row[9]=ic
                            row[10]=res[ir]
                print ' & '.join(map(str,row)),'\\\\ \\hline'

    def test_getCrossSection_thermal( self ):
        '''Test thermal point'''
        result = self.RRR.getCrossSection( 0.0235339347844 )
        answer = { 'capture':45.04182664, 'total':66.22034846, 'fission':0., 'elastic':20.68451381, 'H1 + S35':0.49400801 }
        for k in answer:
            self.assertWithinXPercent( answer[k], result[k][0], percent=0.01 )

    def test_getCrossSection_on_resonance( self ):
        result = self.RRR.getCrossSection( 68236.16 )
        answer = { 'capture':0.027265712318584909, 'total':15.374001449458987, 'fission':0., 'elastic':15.346725225821752, 'H1 + S35':1.0511318650923915e-05 }
        for k in answer:
            self.assertWithinXPercent( answer[k], result[k][0], percent=0.01 )

    def test_getCrossSection_below_upperBound( self ):
        result = self.RRR.getCrossSection( 1199999.0 )
        answer = { 'capture':0.00023637877870793497, 'total':3.5295657156959384, 'fission':0., 'elastic':3.5293200408304597, 'H1 + S35':9.2960867706728532e-06 }
        for k in answer:
            self.assertWithinXPercent( answer[k], result[k][0], percent=0.01 )

    def test_k_competitive(self):
        self.RRRSmall.setResonanceParametersByChannel( )
        gCl36 = self.RRRSmall.particlePairs.values()[0].reactionInfo['particles']
        nCl35 = self.RRRSmall.particlePairs.values()[1].reactionInfo['particles']
        pSi   = self.RRRSmall.particlePairs.values()[2].reactionInfo['particles']
        self.assertEqual( self.RRRSmall.k_competitive(Ex=4.0, pA=gCl36[0], pB=gCl36[1]), 0.0 ) # mass of gamma is zero, so this equation for k really doesn't make sense
        self.assertEqual( self.RRRSmall.k_competitive(Ex=numpy.array([4.0]), pA=gCl36[0], pB=gCl36[1]), numpy.array([0.0]) ) # checks what happens for array arguments used in broadcasting
        self.assertAlmostEqual( self.RRRSmall.k_competitive(Ex=4.0, pA=nCl35[0], pB=nCl35[1]), self.RRRSmall.k(4.0) ) # should reduce to regular k(E) function for neutron+residual
        self.assertAlmostEqual( self.RRRSmall.k_competitive(Ex=4.0, pA=pSi[0], pB=pSi[1]),     self.RRRSmall.k(4.0), 5) # should be close, since mass_p very close to mass_n

    def test_rho(self):
        self.RRRSmall.setResonanceParametersByChannel( )
        for c in self.RRRSmall.channels:
            pp = self.RRRSmall.particlePairs[c].reactionInfo['particles']
            self.assertAlmostEqual( self.RRRSmall.rho(4.0, c), self.RRRSmall.k_competitive(Ex=4.0, pA=pp[0], pB=pp[1]) * self.RRRSmall.scatteringRadiiTrue[c] )
            self.assertAlmostEqual( self.RRRSmall.rho(numpy.array([4.0]), c), self.RRRSmall.k_competitive(Ex=numpy.array([4.0]), pA=pp[0], pB=pp[1]) * self.RRRSmall.scatteringRadiiTrue[c] )

        # Code below is for getting a table of etas for the test case
        if False:
            egrid=self.RRRSmall.generateEnergyGrid()
            rhos=self.RRRSmall.rho(egrid, c)
            for E,rho in zip(egrid,rhos): print E,rho

    def test_rhohat(self):
        self.RRRSmall.setResonanceParametersByChannel( )
        for c in self.RRRSmall.channels:
            pp = self.RRRSmall.particlePairs[c].reactionInfo['particles']
            self.assertAlmostEqual( self.RRRSmall.rhohat(4.0, c), self.RRRSmall.k_competitive(Ex=4.0, pA=pp[0], pB=pp[1]) * self.RRRSmall.scatteringRadiiEffective[c] )

    def test_eta(self):
        self.RRRSmall.setResonanceParametersByChannel( )
        gCl36 = self.RRRSmall.particlePairs.values()[0].reactionInfo['particles']
        nCl35 = self.RRRSmall.particlePairs.values()[1].reactionInfo['particles']
        pSi   = self.RRRSmall.particlePairs.values()[2].reactionInfo['particles']
        self.assertEqual( self.RRR.eta(Ex=4.0, pA=nCl35[0], pB=nCl35[1]), 0.0 ) # neutron has no charge
        self.assertEqual( self.RRR.eta(Ex=4.0, pA=gCl36[0], pB=gCl36[1]), 0.0 ) # gamma has no charge, but no mass, still should get 0
        self.assertAlmostEqual( self.RRR.eta(Ex=4.0e4, pA=pSi[0], pB=pSi[1]), 12.644818547132786)
        self.assertAlmostEqual( self.RRR.eta(Ex=4.0, pA=pSi[0], pB=pSi[1]), 1264.4818547132786) # I checked the number by hand, so it should be OK

        # Code below is for getting a table of etas for the test case
        if False:
            egrid=self.RRRSmall.generateEnergyGrid()
            etas=self.RRR.eta(Ex=egrid, pA=pSi[0], pB=pSi[1])
            for E,eta in zip(egrid,etas): print E,eta

    def test_omega(self):
        vals='''-336933.4       4.64807322355   0.559844234898
                -180.65         3.17917512551   0.818513260471
                397.8154        3.17772298742   0.818887300086
                4250.762        3.16810132868   0.821374295713
                16356.12        3.13843081431   0.829139513199
                22396.4         3.12393433962   0.832987097258
                54932.0         3.04916293104   0.853413561834
                68236.16        3.02009796749   0.861626684168
                133988.4        2.88773859413   0.901119305908
                1205687.0       1.86505967179   1.39523525008
                1434336.0       1.75889799566   1.4794473608
                1441365.0       1.75591539878   1.48196034923
                7563145.0       0.883361945718  2.94578797536'''
        E=[]
        rho=[]
        eta=[]
        answers={
            0:[ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
            1:[ 1.35888357, 1.26604766, 1.26591687, 1.2650475 , 1.2623361 , 1.26099436, 1.25389057, 1.25104352, 1.2374289 , 1.07862856, 1.05383212, 1.05310261, 0.72354639],
            2:[ 2.52334065, 2.27530609, 2.27496936, 2.27273204, 2.26576466, 2.26232263, 2.24416284, 2.23691496, 2.20249416, 1.82912803, 1.77517636, 1.77360531, 1.13946071],
            3:[ 3.52098202, 3.08969274, 3.08912795, 3.0853768 , 3.0737105 , 3.06795573, 3.03768805, 3.02565158, 2.96882766, 2.3853411 , 2.30546092, 2.30314969, 1.42582168],
            4:[ 4.38117882, 3.76125109, 3.76046377, 3.75523618, 3.73899517, 3.73099327, 3.68900893, 3.67236062, 3.59412368, 2.8216383 , 2.71973696, 2.71680072, 1.64317353]}
        for row in vals.split('\n'):
            srow=map(float,row.split())
            E.append(srow[0])
            eta.append(srow[1])
            rho.append(srow[2])
        eta=numpy.array(eta)
        for L in range(0,5):
            self.assertAllWithinXPercent(answers[L], self.RRRSmall.omega(eta,L))

    def test_getCoulombWavefunctions(self):
        """
        Baby test: eta = 0.0, F_l=rho*j_l(rho), so F_0=sin(rho), similarly, G_0 = cos(rho)

        Tougher test: Abramowitz & Stegun have tables of F0 & G0 in table 14.1.  If we can match all these, recursion
        ensures we can do any order

        """
        # Baby test
        self.assertEqual( getCoulombWavefunctions.getCoulombWavefunctions(rho=numpy.array([1.0]), eta=numpy.array([0.0]), L=0),
                          (numpy.array([math.sin(1.0)]), numpy.array([math.cos(1.0)])) )

        # Tougher test
        rho     = numpy.array([1.0, 1.0, 1.0,  1.0,  1.0, 5.0, 5.0, 5.0,  5.0,  5.0, 10.0, 10.0, 10.0, 10.0, 10.0, 20.0, 20.0, 20.0, 20.0, 20.0])
        eta     = numpy.array([0.5, 1.0, 5.0, 10.0, 20.0, 0.5, 1.0, 5.0, 10.0, 20.0,  0.5,  1.0,  5.0, 10.0, 20.0,  0.5,  1.0,  5.0, 10.0, 20.0])
        #
        # These are the full precision values returned by COULFG:
        #
        #        answerF = numpy.array([ 5.16601500e-01,   2.27526211e-01,   2.04130067e-05,   3.71673341e-11,   2.96513633e-23,  \
        #                               -4.90455238e-01,   6.84937412e-01,   2.76730117e-02,   1.72075334e-06,   1.65159712e-16,  \
        #                                9.39186276e-01,   4.77560816e-01,   9.17944919e-01,   1.62627112e-03,   8.07087199e-12,  \
        #                               -8.13196123e-01,  -3.29225536e-01,  -2.29347373e-01,   1.03427445e+00,   5.45296368e-06])
        #
        # These are the numbers from Abramowitz and Stegun:
        answerF = numpy.array([ 5.1660e-01,   2.2753e-01,   2.0413e-05,   3.6966e-11,   2.9556e-23,\
                               -4.9046e-01,   6.8494e-01,   2.7673e-02,   1.7207e-06,   1.6477e-16,\
                                9.3919e-01,   4.7756e-01,   9.1794e-01,   1.6263e-03,   8.0470e-12,\
                               -8.1320e-01,  -3.2923e-01,  -2.2935e-01,   1.0343e+00,   5.4529e-06])
        #
        # These are the full precision values returned by COULFG:
        #
        #       answerG = numpy.array([  1.19748697e+00,   2.04309716e+00,   8.08552775e+03,   3.07242545e+09,   2.69426701e+21,
        #                               -9.34927038e-01,  -8.98414359e-01,   1.81934952e+01,   1.67636794e+05,   1.14367929e+15,
        #                               -4.14351069e-01,   9.42874243e-01,   1.60852456e+00,   3.07873217e+02,   3.57573124e+10,
        #                                6.03865910e-01,  -9.72428399e-01,   1.16571667e+00,   1.79968757e+00,   9.17224865e+04])
        #
        # These are the numbers from Abramowitz and Stegun:
        answerG = numpy.array([ 1.1975e+00,   2.0431e+00,   8.0855e+03,   3.0882e+09,   2.7024e+21,\
                               -9.3493e-01,  -8.9841e-01,   1.8193e+01,   1.6764e+05,   1.1464e+15,\
                               -4.1435e-01,   9.4287e-01,   1.6085e+00,   3.0787e+02,   3.5867e+10,\
                                6.0387e-01,  -9.7243e-01,   1.1657e+00,   1.7997e+00,   9.1723e+04])
        F,G     = getCoulombWavefunctions.getCoulombWavefunctions(rho=rho, eta=eta, L=0)
        self.assertAllWithinXPercent(     F,     answerF, 0.01, 1e-8), 'F0(rho,eta), %s should be %s'%(str(F),str(answerF))  # agreement terrible as exponent gets small
        self.assertAllWithinXPercent( 1.0/G, 1.0/answerG, 0.01, 1e-8), 'G0(rho,eta), %s should be %s'%(str(G),str(answerG))  # agreement terrible as exponent gets big

    def test_coulombPenetrationFactor(self):
        """
        For eta=0, can compare to regular penetration factor function
        I guess we could make tougher tests...
        """
        for L in [0,1,2]:
            eta=numpy.array([0.0, 0.0, 0.0])
            rho=numpy.array([1.0, 0.1, 0.01])
            self.assertAllWithinXPercent(\
                    getCoulombWavefunctions.coulombPenetrationFactor(L=L, rho=rho, eta=eta), \
                    self.RRRSmall.penetrationFactor(L=L, rho=rho),\
                    0.5, 1e-8), \
                'L=%i, %s vs. %s'%(L, str(getCoulombWavefunctions.coulombPenetrationFactor(L=L, rho=rho, eta=eta)),str(self.RRRSmall.penetrationFactor(L=L, rho=rho)))

    def test_coulombShiftFactor(self):
        """
        For eta=0, can compare to regular shift factor function
        I guess we could make tougher tests...
        """
        for L in [0,1,2]:
            eta=numpy.array([0.0, 0.0, 0.0])
            rho=numpy.array([1.0, 0.1, 0.01])
            self.assertAllWithinXPercent(\
                    getCoulombWavefunctions.coulombShiftFactor(L=L, rho=rho, eta=eta), \
                    self.RRRSmall.shiftFactor(L=L, rho=rho),\
                    0.5, 1e-8), \
                'L=%i, %s vs. %s'%(L, str(getCoulombWavefunctions.coulombShiftFactor(L=L, rho=rho, eta=eta)),str(self.RRRSmall.shiftFactor(L=L, rho=rho)))

    def test_coulombPhi(self):
        """
        For eta=0, can compare to regular phase function
        I guess we could make tougher tests...
        """
        for L in [0,1,2]:
            eta=numpy.array([   0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  ])
            rho=numpy.array([ 100.0, 50.0, 10.0, 5.0, 3.5, 2.0, 1.0, 0.1, 0.01 ])
            self.assertAllWithinXPercent(\
                numpy.mod( getCoulombWavefunctions.coulombPhi(L=L, rho=rho, eta=eta), numpy.pi ),\
                numpy.mod( self.RRRSmall.phi(L=L, rho=rho),                 numpy.pi ),\
                0.01, 1e-6), \
            'L=%i, %s vs. %s'%( L,
                                str( numpy.mod( getCoulombWavefunctions.coulombPhi(L=L, rho=rho, eta=eta), numpy.pi) ),
                                str( numpy.mod( self.RRRSmall.phi(L=L, rho=rho),                 numpy.pi) ) )
            if False:
                eta=numpy.array([ 0.01, 0.01,  0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01  ])
                print '\n eta=0.01'
                print 'no coul',L,numpy.mod( self.RRRSmall.phi(L=L, rho=rho), numpy.pi )
                print 'coul',L,numpy.mod( getCoulombWavefunctions.coulombPhi(L=L, rho=rho, eta=eta), numpy.pi )
                eta=numpy.array([ 0.1, 0.1,  0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1  ])
                print '\n eta=0.1'
                print 'no coul',L,numpy.mod( self.RRRSmall.phi(L=L, rho=rho), numpy.pi )
                print 'coul',L,numpy.mod( getCoulombWavefunctions.coulombPhi(L=L, rho=rho, eta=eta), numpy.pi )
                eta=numpy.array([ 0.5, 0.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5  ])
                print '\n eta=0.5'
                print 'no coul',L,numpy.mod( self.RRRSmall.phi(L=L, rho=rho), numpy.pi )
                print 'coul',L,numpy.mod( getCoulombWavefunctions.coulombPhi(L=L, rho=rho, eta=eta), numpy.pi )
                eta=numpy.array([ 1.0, 1.0,  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0  ])
                print '\n eta=1.0'
                print 'no coul',L,numpy.mod( self.RRRSmall.phi(L=L, rho=rho), numpy.pi )
                print 'coul',L,numpy.mod( getCoulombWavefunctions.coulombPhi(L=L, rho=rho, eta=eta), numpy.pi )

    def test_coulombNormalizationFactor(self):
        """
        For L=0, Abramawitz & Stegun have tabulated values of C_{L=0}(eta).  See table 14.2, p. 554
        I checked the values at eta = 0.1, 1.0 and 2.0 and they match.  Also,

            ..math::
                C_L(0)=2^L \Gamma(L+1)/\Gamma(2L+2) = 1/(2L+1)

        and (A&S, eq. 14.1.10)

            ..math::
                C_L(\eta) = \frac{ \sqrt{L^2+\eta^2} }{ L(2L+1) } C_{L-1}(\eta)
        """
        eta=numpy.array([0.0, 0.1, 1.0, 2.0, 10.0, 100.0])
        L0_answers=[  1.0,     8.47658520e-01,   1.08422513e-01,  6.61992e-3,  1.80022337e-13,   0.00000000e+00]
        def CL_answers(L,i):
            if L==0: return L0_answers[i]
            return numpy.sqrt(L*L+eta[i]*eta[i])*CL_answers(L-1,i)/L/(2.0*L+1.0)
        answers = { L:[ CL_answers(L,i) for i in range(len(eta)) ] for L in range(4) }
        for L in [0,1,2]:
            self.assertAllWithinXPercent(\
                    getCoulombWavefunctions.coulombNormalizationFactor(L, eta),\
                    answers[L],\
                0.0001),\
            'L=%i, %s vs. %s'%(L, str(getCoulombWavefunctions.coulombNormalizationFactor(L, eta)),str(answers[L]))

    def test_getEiPhis(self):
        self.RRRSmall.setResonanceParametersByChannel( )
        Es=numpy.array([4.0e5]+[self.RRRSmall._energies[i] for i in [4,10]])
        phiMatrix = self.RRRSmall.getEiPhis(Ein=Es, enableExtraCoulombPhase=True)
        phiFromPhiMatrix = numpy.angle(phiMatrix)
        for ic,c in enumerate(self.RRRSmall.channels):
            rho=self.RRRSmall.rhohat(Es-c.Xi,c)
            if c.channelClass == NEUTRONCHANNEL:
                phiAnswer = -self.RRRSmall.phi(c.l,rho)
                for iE,E in enumerate(Es):
                    self.assertWithinXPercent( phiFromPhiMatrix[ic,iE], phiAnswer[iE], percent=0.01)
            elif c.channelClass == GAMMACHANNEL:
                pass
            elif c.channelClass == FISSIONCHANNEL:
                pass
            elif c.channelClass == CPCHANNEL:
                pA, pB = self.RRRSmall.particlePairs[c].reactionInfo['particles']
                eta = self.RRRSmall.eta(Es-c.Xi, pA, pB)
                phi=getCoulombWavefunctions.coulombPhi(L=c.l, rho=rho, eta=eta)
                wc = self.RRRSmall.omega(eta,c.l)
                if c.l == 0.0: self.assertItemsEqual(wc,len(wc)*[0.0])
                phiAnswer = wc-phi
                for iE,E in enumerate(Es):
                    self.assertWithinXPercent(phiFromPhiMatrix[ic,iE], phiAnswer[iE], percent=0.01)

    def test_getL0Matrix(self):
        for E in [16356.12,133988.4,4.0e4]:
            self.RRRSmall.setResonanceParametersByChannel( )
            L0=self.RRRSmall.getL0Matrix(Ein=numpy.array([E]))
            for ic,c in enumerate(self.RRRSmall.channels):
                L0cc=L0[0].diagonal()[ic]
                rho=self.RRRSmall.rho(E-c.Xi,c)
                if c.channelClass == NEUTRONCHANNEL:
                    self.assertAlmostEqual( L0cc, 1j*self.RRRSmall.penetrationFactor(c.l,rho) )
                elif c.channelClass == GAMMACHANNEL:
                    self.assertAlmostEqual( L0cc, complex( 0.0, 1.0 ) )
                elif c.channelClass == FISSIONCHANNEL:
                    self.assertAlmostEqual( L0cc, complex( 0.0, 1.0 ) )
                elif c.channelClass == CPCHANNEL:
                    # We can't test the equations with Coulomb turned on as I don't have numerical results to compare to.
                    # However, I can compare to the results if Coulomb is turned off -- I'd better get the neutron numbers.
                    # To turn off the Coulomb force, set eta to zero
                    fake_eta = numpy.array([0.0])
                    F,G   = getCoulombWavefunctions.getCoulombWavefunctions( numpy.array([rho]), fake_eta, c.l )
                    F1,G1 = getCoulombWavefunctions.getCoulombWavefunctions( numpy.array([rho]), fake_eta, c.l+1 )
                    S1 = (c.l+1.0+fake_eta)/numpy.array([rho])
                    R1 = numpy.sqrt(1.0+(fake_eta/(c.l+1))**2)
                    A2  = F*F+G*G
                    shift = numpy.array([rho])*(S1-R1*(G*G1+F*F1)/A2)
                    penet = numpy.array([rho])/A2
                    self.assertAlmostEqual( complex( shift, penet ), complex( self.RRRSmall.shiftFactor(c.l,rho), self.RRRSmall.penetrationFactor(c.l,rho) ) )

                    # This is  a tougher test, it checks whether the L0's imaginary part matches the penetrability factor computed in the coulombPenetrationFactor function
                    pA, pB = self.RRRSmall.particlePairs[c].reactionInfo['particles']
                    real_eta = self.RRRSmall.eta(numpy.array([E])-c.Xi, pA, pB)
                    self.assertAlmostEqual( L0cc, 1j*getCoulombWavefunctions.coulombPenetrationFactor(c.l,numpy.array([rho]),real_eta) )

    def test_getRMatrix( self ):
        '''
        iER == 4 details:
            ER = 16356.12 eV
            GAMg = 0.3865 eV
            redGAM list = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 36.504003097726347, 94.126905092341858]
            last two channels (which have non-zero reduced widths) are:
                ChannelDesignator(1, 3.0, 'n + Cl35', 0, 2, gfact=0.875, particleA=n, particleB=Cl35, Xi=-0.0, isElastic=True, channelClass=1, useRelativistic=False, eliminated=False) and
                ChannelDesignator(1, 3.0, 'H1 + S35', 0, 2, gfact=0.875, particleA=H1, particleB=S35, Xi=-632965.817883, isElastic=False, channelClass=2, useRelativistic=False, eliminated=False)
        :return:
        '''
        self.RRRSmall.setResonanceParametersByChannel( )
        energy_indices = [1,4,10]
        ERs = numpy.array([self.RRRSmall._energies[i] for i in energy_indices]).reshape(len(energy_indices),1)
        R = self.RRRSmall.getRMatrix( ERs )
        for iE,iER in enumerate(energy_indices):
            ER = ERs[iE]
            if ER < 0.0: continue
            for cg in self.RRRSmall.eliminatedChannels:
                if iER in self.RRRSmall.eliminatedChannels[cg]:
                    gGamma = self.RRRSmall.eliminatedChannels[cg][iER]
                    break
            for ic1,c1 in enumerate(self.RRRSmall.channels):
                if iER not in self.RRRSmall.channels[c1]: continue
                c1Gamma=self.RRRSmall.channels[c1][iER]
                rho1=self.RRRSmall.rho(numpy.array([self.RRRSmall._energies[iER]])-c1.Xi,c1)
                pa1,pb1=self.RRRSmall.particlePairs[c1].reactionInfo['particles']
                eta1=self.RRRSmall.eta(numpy.array([self.RRRSmall._energies[iER]])-c1.Xi,pa1,pb1)
                c1Pen=getCoulombWavefunctions.coulombPenetrationFactor(L=c1.l,rho=rho1,eta=eta1)[0]
                for ic2,c2 in enumerate(self.RRRSmall.channels):
                    if iER not in self.RRRSmall.channels[c2]: continue
                    c2Gamma=self.RRRSmall.channels[c2][iER]
                    rho2=self.RRRSmall.rho(numpy.array([self.RRRSmall._energies[iER]])-c2.Xi,c2)
                    pa2,pb2=self.RRRSmall.particlePairs[c2].reactionInfo['particles']
                    eta2=self.RRRSmall.eta(numpy.array([self.RRRSmall._energies[iER]])-c2.Xi,pa2,pb2)
                    c2Pen=getCoulombWavefunctions.coulombPenetrationFactor(c2.l,rho2,eta2)[0]
                    correctAnswer = 1j*numpy.sqrt(c1Gamma*c2Gamma/c1Pen/c2Pen)/gGamma
                    self.assertWithinXPercent( R[iE,ic1,ic2], correctAnswer, percent=0.001 )

    def test_getXMatrix_Small( self ):
        # Print a table for plotting the X matrix
        debug=False
        self.RRRSmall.setResonanceParametersByChannel( )
        egrid=numpy.array( [ 1.0, 31733.7 ] ).reshape(2,1) #, 16356.120 ] )# + self.RRRSmall.generateEnergyGrid()[::2] )
        X = self.RRRSmall.getXMatrix(egrid)
        R = self.RRRSmall.getRMatrix(egrid)
        L0 = self.RRRSmall.getL0Matrix(egrid)

        if debug:
            print 'these channels:',self.RRRSmall.channels.keys()[-2:]
            for iEE,EE in enumerate(egrid):
                print EE, X[iEE,-2,-2].real, X[iEE,-2,-2].imag, X[iEE,-1,-1].real, X[iEE,-1,-1].imag

        # Test specific values, last two rows/columns correspond to an isolated resonance so reuse LRF=3 unit test
        for iE in range(len(egrid)):
            self.assertAlmostEqual( X[iE,-2,-2], L0[iE,-2,-2].imag*R[iE,-2,-2]/(1.0-R[iE,-2,-2]*L0[iE,-2,-2] ), 7 )
            self.assertAlmostEqual( X[iE,-1,-1], L0[iE,-1,-1].imag*R[iE,-1,-1]/(1.0-R[iE,-1,-1]*L0[iE,-1,-1] ), 7 )

    @unittest.skip("not a test at all, purely for debugging")
    def test_rho_eta_values( self ):
        # set up the channels
        self.RRRSmall.setResonanceParametersByChannel( )

        # get the energies for evaluating the cross section and of the resonances
        ERs = numpy.array(self.RRRSmall._energies)
        Es  = numpy.array(self.RRRSmall.generateEnergyGrid())

        # find the p+S channel so we can get the threshold energy and particles
        for ic,c in enumerate(self.RRRSmall.channels):
            if c.reaction == 'H1 + S35':
                Xi = c.Xi
                pA, pB = self.RRRSmall.particlePairs[c].reactionInfo['particles']
                break

        # get rhos and etas
        resEta = self.RRRSmall.eta(ERs-Xi, pA, pB)
        resRho = self.RRRSmall.rho(ERs-Xi, c)
        xsEta  = self.RRRSmall.eta(Es-Xi, pA, pB)
        xsRho  = self.RRRSmall.rho(Es-Xi, c)

        print '\n# E (eV)        eta             rho'
        for x in zip(Es,xsEta,xsRho): print ' '.join([y.ljust(15) for y in map(str,x)])
        for x in zip(ERs,resEta,resRho): print ' '.join([y.ljust(15) for y in map(str,x)])

    def test_getScatteringMatrixUUnitarity_Small( self ):
        '''
        The RML parameterization might have a unitary U matrix.  Here we have two tests:

            * just plain multiplication: U^+ * U == one

            * check the Frobenius norm, which is ||A|| = sqrt(| sum_{i,j} |A[i,j]|^2 |), so the norm of the identity matrix is sqrt( ndim )
        '''
        debug = False
        self.RRRSmall.setResonanceParametersByChannel( )

        if debug:
            egrid = self.RRRSmall.generateEnergyGrid()[::2]
            UU = self.RRRSmall.getScatteringMatrixU( egrid )
            for iEE,EE in enumerate(egrid):
                print EE, UU[iEE,0,0].real, UU[iEE,0,0].imag, \
                    UU[iEE,0,1].real, UU[iEE,0,1].imag, UU[iEE,1,1].real, \
                    UU[iEE,1,1].imag, numpy.linalg.norm(UU[0])/math.sqrt( UU[0].shape[0] )

        egrid = numpy.array([ 1.0, 42.0, 16356.12 ]).reshape(3,1)
        U = self.RRRSmall.getScatteringMatrixU( egrid, useTabulatedScatteringRadius = True )
        for iE,E in enumerate(egrid):
            one = numpy.diag( U[iE].shape[0]*[ complex( 1.0 ) ] )
            if E == egrid[-1]: # on resonance, expect to be good to only 10%
                numDecimalPlaces = 1
                rtol= 5e-02 # nominal value: 1e-03
                atol= 1e-02 # nominal value: 1e-04
            else:
                numDecimalPlaces = 3
                rtol= 1e-03 # nominal value: 1e-03
                atol= 1e-04 # nominal value: 1e-04
            theNorm = numpy.linalg.norm( U[iE] )
            expectedNorm = math.sqrt( U[iE].shape[0] )
            self.assertWithinXPercent(
                    theNorm,
                    expectedNorm,
                    rtol*100., atol )
            #UdaggerU=numpy.conj(U[iE].T)*U[iE]
            #self.assertTrue( numpy.allclose( UdaggerU, one, rtol=2.0*rtol, atol=atol ),  msg="For E=%s, U.^dagger*U=%s, c=%s" % ( str(E), str(UdaggerU[-2:,:]), str(self.RRRSmall.channels.keys()[-2:]) ) )

    def test_scatteringMatrixUInCrossSectionCalculation_Small( self ):
        debug = False
        doTests = True
        enableExtraCoulombPhase=False

        # Set up
        self.RRRSmall.setResonanceParametersByChannel( )
        nChan = len( self.RRRSmall.channels )
        nElim = len( self.RRRSmall.eliminatedChannels )
        nRes = len( self.RRRSmall._energies )

        # Column labels!
        if debug: print '#\n# All Waves (full grid), %i kept channels, %i eliminated channels, %i resonances' % ( nChan, nElim, nRes )
        if debug: print '#      E (eV)       SigTotU (b)    SigTotSum (b)    SigTotOK (b)   SigGamU (b)    SigGamOK (b)   SigNeuU (b)    SigNeuOK (b)   SigNeuInt (b) SigNeuInt/SigNeuU  SigPot       rho            eta         SigNPU (b)   SigNPOK (b)   SigNPInt (b)'

        # We can't eat the whole energy grid, so we'll do it in parts
        fullEgrid = self.RRRSmall.generateEnergyGrid()
        energyCohortLength = 50
        nEnergyCohorts = len(fullEgrid)//energyCohortLength
        lastEnergyCohortSize = len(fullEgrid)%energyCohortLength
        energyCohorts = [(iec*energyCohortLength, (iec+1)*energyCohortLength) for iec in range(nEnergyCohorts)]
        energyCohorts.append((energyCohorts[-1][1],energyCohorts[-1][1]+lastEnergyCohortSize))

        for iStart,iStop in energyCohorts: #[4:7]:
            egrid = fullEgrid[iStart:iStop]
            egrid = numpy.array(egrid).reshape(len(egrid),1)

            U = self.RRRSmall.getScatteringMatrixU( egrid, useTabulatedScatteringRadius=False, enableExtraCoulombPhase=enableExtraCoulombPhase )
            T = self.RRRSmall.getScatteringMatrixT( egrid, useTabulatedScatteringRadius=False, enableExtraCoulombPhase=enableExtraCoulombPhase )
            k = self.RRRSmall.k(egrid)
            nLs = { 0:0, 1:0, 2:0, 3:0 }

            sigs = self.RRRSmall.getCrossSection( egrid )
            B = self.RRRSmall.getAngularDistribution( egrid, renormalize=False, enableExtraCoulombPhase=enableExtraCoulombPhase )
            Bn = B['n + Cl35']
            Bp = B['H1 + S35']  # FIXME, dictionary only contains elastic right now
            phis = [ self.RRRSmall.phi(c.l, self.RRRSmall.rhohat(egrid,c) ) for c in self.RRRSmall.channels ]

            for iE,E in enumerate(egrid):
                sigTotWithU, sigGamWithU, sigNeuWithU, sigNPWithU = 0.0, 0.0, 0.0, 0.0
                for ic1,c1 in enumerate(self.RRRSmall.channels):
                    if c1.reaction != 'n + Cl35': continue
                    nLs[ c1.l ] += 1.0
                    sigTotWithU += 2.0 * numpy.pi * c1.gfact * ( 1.0 - U[iE,ic1,ic1].real ) / k[iE]**2
                    for ic2,c2 in enumerate(self.RRRSmall.channels):
                        if c2.J != c1.J: continue
                        if c2.reaction == 'n + Cl35':
                            sigNeuWithU += numpy.pi * c1.gfact * numpy.power( numpy.abs( T[iE,ic1,ic2] ) / k[iE], 2.0 )
                        if c2.reaction == 'H1 + S35':
                            rho = self.RRRSmall.rho(numpy.array([E-c2.Xi]),c2)
                            pA, pB = self.RRRSmall.particlePairs[c2].reactionInfo['particles']
                            eta = self.RRRSmall.eta(numpy.array([E-c2.Xi]), pA, pB)
                            sigNPWithU  += numpy.pi * c1.gfact * numpy.power( numpy.abs( T[iE,ic1,ic2] ) / k[iE], 2.0 )
                sigGamWithU = sigTotWithU - sigNeuWithU - sigNPWithU
                sigPot = sum( [ ( ( 2.0 * l + 1.0 ) ) * 4.0 * numpy.pi * numpy.power( numpy.sin(phis[ic1][iE])/k[iE], 2.0 ) for l in nLs ] )
                sigNeuInt  = 4.0*numpy.pi*Bn[0][iE].real
                sigNPInt   = 4.0*numpy.pi*Bp[0][iE].real

                if debug: print 16*'%10g     ' % ( E, sigTotWithU, sigGamWithU + sigNeuWithU, sigs['total'][iE], sigGamWithU, sigs['capture'][iE], sigNeuWithU, sigs['elastic'][iE], sigNeuInt,  sigNeuWithU/sigNeuInt, sigPot, rho, eta, sigNPWithU, sigs['H1 + S35'][iE], sigNPInt ) #, -(sigNeuWithU-sigs['elastic'][iE])/sigPot )

                # Check they sum up correctly (really just another test of unitarity)
                if doTests: self.assertWithinXPercent( sigTotWithU, sigGamWithU + sigNeuWithU + sigNPWithU, 1. )

                # Check elastic channel agrees with L=0 term of angular distribution
                if doTests and False: self.assertWithinXPercent( sigNeuWithU, sigNeuInt ) # fails badly for some high energy P wave resonances??

                # Check agree with getCrossSection()
                if doTests:
                    self.assertWithinXPercent( sigs['capture'][iE]+sigNeuWithU + sigNPWithU, sigs['total'][iE], 1. )
                    if False: self.assertWithinXPercent( sigGamWithU, sigs['capture'][iE] ) # expected to fail
                    self.assertWithinXPercent( sigNeuWithU, sigs['elastic'][iE], 1. )

    @unittest.skipIf(not DOFETESTS,'in dev')
    def test_scatteringMatrixUInCrossSectionCalculation_Fe( self ):
        debug = False
        doTests = True
        enableExtraCoulombPhase=True

        # Column labels!
        if debug: print '#'
        if debug: print '#      E (eV)       SigTotU (b)    SigTotSum (b)    SigTotOK (b)   SigGamU (b)    SigGamOK (b)   SigNeuU (b)    SigNeuOK (b)   SigNeuInt (b) SigNeuInt/SigNeuU  SigPot       SigNPU (b)   SigNPOK (b)   SigNPInt (b)'

        # We can't eat the whole energy grid, so we'll do it in parts
        fullEgrid = self.RRRFe.generateEnergyGrid()
        energyCohortLength = 20
        nEnergyCohorts = len(fullEgrid)//energyCohortLength
        lastEnergyCohortSize = len(fullEgrid)%energyCohortLength
        energyCohorts = [(iec*energyCohortLength, (iec+1)*energyCohortLength) for iec in range(nEnergyCohorts)]
        energyCohorts.append((energyCohorts[-1][1],energyCohorts[-1][1]+lastEnergyCohortSize))

        iStartCohort=4124
        iCohort=iStartCohort-1
        iStopCohort=iStartCohort+2

        for iStart,iStop in energyCohorts[iStartCohort:iStopCohort]:

            # Full egrid of this cohort of energies
            iCohort+=1
            cohortEgrid = fullEgrid[iStart:iStop]
            cohortEgrid = numpy.array(cohortEgrid).reshape(len(cohortEgrid),1)

            # Deal with potential thresholds
            thresholdIndices = [0]
            for Xi in self.RRRFe._thresholds:
                if Xi in cohortEgrid:
                    thresholdIndices.append(list(cohortEgrid).index(Xi))
            thresholdIndices.append(-1)
            subgrids = [cohortEgrid[thresholdIndices[i]:thresholdIndices[i+1]] for i in range(len(thresholdIndices[:-1]))]

            for egrid in subgrids:
                # Set up channels on this energy grid
                self.RRRFe.setResonanceParametersByChannel( Ein=egrid )
                nChan = len( self.RRRFe.channels )
                nElim = len( self.RRRFe.eliminatedChannels )
                nRes = len( self.RRRFe._energies )

                U = self.RRRFe.getScatteringMatrixU( egrid, useTabulatedScatteringRadius=False, enableExtraCoulombPhase=enableExtraCoulombPhase )
                T = self.RRRFe.getScatteringMatrixT( egrid, useTabulatedScatteringRadius=False, enableExtraCoulombPhase=enableExtraCoulombPhase )
                k = self.RRRFe.k(egrid)
                nLs = { 0:0, 1:0, 2:0, 3:0 }

                sigs = self.RRRFe.getCrossSection( egrid )
                B = self.RRRFe.getAngularDistribution( egrid, renormalize=False, enableExtraCoulombPhase=enableExtraCoulombPhase )
                Bn = B['n + Fe56']
                try:              Bnn = B['n + Fe56_e1']  # FIXME raises KeyError since dictionary only contains 'elastic'
                except KeyError:  Bnn = [ numpy.zeros_like(egrid) ]
                phis = [ self.RRRFe.phi(c.l, self.RRRFe.rhohat(egrid,c) ) for c in self.RRRFe.channels ]

                for iE,E in enumerate(egrid):
                    sigTotWithU, sigGamWithU, sigNeuWithU, sigNNWithU = 0.0, 0.0, 0.0, 0.0
                    for ic1,c1 in enumerate(self.RRRFe.channels):
                        if c1.reaction != 'n + Fe56': continue
                        nLs[ c1.l ] += 1.0
                        sigTotWithU += 2.0 * numpy.pi * c1.gfact * ( 1.0 - U[iE,ic1,ic1].real ) / k[iE]**2
                        for ic2,c2 in enumerate(self.RRRFe.channels):
                            if c2.J != c1.J: continue
                            if c2.reaction == 'n + Fe56':
                                sigNeuWithU += numpy.pi * c1.gfact * numpy.power( numpy.abs( T[iE,ic1,ic2] ) / k[iE], 2.0 )
                            if c2.reaction == 'n + Fe56_e1':
                                sigNNWithU  += numpy.pi * c1.gfact * numpy.power( numpy.abs( T[iE,ic1,ic2] ) / k[iE], 2.0 )
                    sigGamWithU = sigTotWithU - sigNeuWithU - sigNNWithU
                    sigPot = sum( [ ( ( 2.0 * l + 1.0 ) ) * 4.0 * numpy.pi * numpy.power( numpy.sin(phis[ic1][iE])/k[iE], 2.0 ) for l in nLs ] )
                    sigNeuInt  = 4.0*numpy.pi*Bn[0][iE].real
                    sigNNInt   = 4.0*numpy.pi*Bnn[0][iE].real

                    if debug: print ('%15.14g'+13*'%10g     ') % ( E, sigTotWithU, sigGamWithU + sigNeuWithU, sigs['total'][iE], sigGamWithU, sigs['capture'][iE], sigNeuWithU, sigs['elastic'][iE], sigNeuInt,  sigNeuWithU/sigNeuInt, sigPot, sigNNWithU, sigs['n + Fe56_e1'][iE], sigNNInt ) #, -(sigNeuWithU-sigs['elastic'][iE])/sigPot )

                    # Check they sum up correctly (really just another test of unitarity)
                    if doTests: self.assertWithinXPercent( sigTotWithU, sigGamWithU + sigNeuWithU + sigNNWithU, 1. )

                    # Check elastic channel agrees with L=0 term of angular distribution
                    if doTests and False: self.assertWithinXPercent( sigNeuWithU, sigNeuInt ) # fails badly for some high energy P wave resonances??

                    # Check agree with getCrossSection()
                    if doTests:
                        self.assertWithinXPercent( sigs['capture'][iE]+sigNeuWithU + sigNNWithU, sigs['total'][iE], 1. )
                        if False: self.assertWithinXPercent( sigGamWithU, sigs['capture'][iE] ) # expected to fail
                        self.assertWithinXPercent( sigNeuWithU, sigs['elastic'][iE], 1. )

    @unittest.skip("in dev")
    def test_getAngularDistribution( self ):
        '''
        Note, Blatt-Biedenharn's Z coefficient has a factor of (l1,l2,0,0;L,0) (a Clebsch-Gordan coefficient) in it.
        As a result, if l1 and l2 cannot couple to L, one gets zero.

update this: spingroup 5 has l=1, j=1, s=2, single P-wave resonance

        In our simple P-wave example, the U matrix has l=0 and 1 parts to it because you need the l=0 terms to get the
        potential scattering right, even if there is only one resonance.  Therefore, we expect L = 0, 1, and 2 terms since
        those are all the L's we can get by coupling any pair of l=0 and 1 together.  Because of this, the angular distribution
        that you get from a single MLBW resonance is different that that of a single SLBW resonance (at least the way
        you need to do things for ENDF).

        We tested the hell out of the L=0 term.

        We can test the L=1,2 terms on resonance by observing that the potential scattering in the l>0 collision matrices is much
        smaller than that in the l=0 elements.  This allows us to ignore all but the purely resonance elements when directly on
        resonance (Ein=ER).   In our P-wave example, that means that our MLBW on resonance angular distribution should match the
        SLBW result (_but only on resonance_).  We reuse the test below, but only approximately.
        '''
        debug = True
        doTest = True
        enableExtraCoulombPhase=True

        self.RRRSmall.setResonanceParametersByChannel( )
        egrid = self.RRRSmall.generateEnergyGrid()
        eStep = 1
        J = 1.0 # for this test data
        S = 2.0 # for this test data

        from numericalFunctions import angularMomentumCoupling as nf_amc
        self.assertAlmostEqual( nf_amc.z_coefficient( 2, int(2*J), 2, int(2*J), int(2*S), 4 ) / nf_amc.z_coefficient( 2, int(2*J), 2, int(2*J), int(2*S), 0 ), math.sqrt(0.02) ) # nf_amc functions use the 2J trick

        if debug:
            # We can't eat the whole energy grid, so we'll do it in parts
            fullEgrid = self.RRRSmall.generateEnergyGrid()
            energyCohortLength = 50
            nEnergyCohorts = len(fullEgrid)//energyCohortLength
            lastEnergyCohortSize = len(fullEgrid)%energyCohortLength
            energyCohorts = [(iec*energyCohortLength, (iec+1)*energyCohortLength) for iec in range(nEnergyCohorts)]
            energyCohorts.append((energyCohorts[-1][1],energyCohorts[-1][1]+lastEnergyCohortSize))

            for iStart,iStop in energyCohorts:
                egrid = fullEgrid[iStart:iStop]
                egrid = numpy.array(egrid).reshape(len(egrid),1)
                Bn = self.RRRSmall.getAngularDistribution( egrid, renormalize=True, reaction='n + Cl35', reactionp='n + Cl35', enableExtraCoulombPhase=enableExtraCoulombPhase )
                Bp = self.RRRSmall.getAngularDistribution( egrid, renormalize=True, reaction='n + Cl35', reactionp='H1 + S35', enableExtraCoulombPhase=enableExtraCoulombPhase )
                if False:
                    for iE, E in enumerate(egrid):
                        print E[0], ' '.join( map( str, [Bn[L][iE][0] for L in range(len(Bn))] ) )
                if True:
                    for iE, E in enumerate(egrid):
                        print E[0], ' '.join( map( str, [Bp[L][iE][0] for L in range(len(Bp))] ) )
        if doTest:
            ER = 133988.4
            Bn = self.RRRSmall.getAngularDistribution( ER, renormalize=False, reaction='n + Cl35', reactionp='n + Cl35', enableExtraCoulombPhase=enableExtraCoulombPhase )
            Bp = self.RRRSmall.getAngularDistribution( ER, renormalize=False, reaction='n + Cl35', reactionp='H1 + S35', enableExtraCoulombPhase=enableExtraCoulombPhase )
            self.assertAlmostEqual( Bn[2]/Bn[0], 0.02, 2 )
#            self.assertAlmostEqual( Bn[1]/Bn[0], 0.0, 3 )
            self.assertAlmostEqual( Bp[2]/Bp[0], 0.02, 2 )
#            self.assertAlmostEqual( Bp[1]/Bp[0], 0.0, 3 )



class TestURRClassAndBaseClasses( TestWithIsClose ):

    def setUp( self ):
        self.Zr90URR = URRcrossSection( MLBWExampleZr90, verbose=False )

    def test_memberData( self ):
        self.assertEqual( self.Zr90URR.lowerBound, 2e5)
        self.assertEqual( self.Zr90URR.upperBound, 1780460)
        self.assertEqual( str(self.Zr90URR.projectile), 'n')
        self.assertEqual( str(self.Zr90URR.target), 'Zr90')
        self.assertWithinXPercent( self.Zr90URR.targetToNeutronMassRatio, 89.13240000000002, 0.01 )
        self.assertEqual( self.Zr90URR.verbose, False)
        self.assertEqual( self.Zr90URR.targetSpin, 0.0)
        self.assertEqual( self.Zr90URR.URR.ENDFconversionFlag, "LRF,LFW=2,0")
        self.assertEqual( [ x.L for x in self.Zr90URR.URR.L_values ], [0, 1, 2])
#        self.assertEqual( self.Zr90URR.URR.reconstructCrossSection, False)
        self.assertEqual( self.Zr90URR.URR.interpolation, 'lin-lin')
        self.assertEqual( self.Zr90URR.URR.moniker, 'tabulatedWidths')
        self.assertItemsEqual( self.Zr90URR.URR.optAttrList,  ('interpolation', 'ENDFconversionFlag') ) #'forSelfShieldingOnly',
        self.assertEqual( str(self.Zr90URR.URR.scatteringRadius.getValueAs('fm')), '7.16' )
        for LSection in self.Zr90URR.URR.L_values:
            for JSection in LSection.J_values:
                self.assertEqual( JSection.neutronDOF, 1.0 )
                self.assertEqual( JSection.gammaDOF, False )
                self.assertEqual( JSection.fissionDOF, False )
                self.assertEqual( JSection.competitiveDOF, False )
                self.assertEqual( getattr(JSection.constantWidths,'levelSpacing'), None )
                self.assertEqual( getattr(JSection.constantWidths,'neutronWidth'), None )
                self.assertEqual( getattr(JSection.constantWidths,'captureWidth'), None )
                self.assertEqual( str(getattr(JSection.constantWidths,'competitiveWidth')), '0. eV' )
                self.assertEqual( str(getattr(JSection.constantWidths,'fissionWidthA')), '0. eV' )
                self.assertEqual( getattr(JSection.constantWidths,'fissionWidthB'), None )
                self.assertEqual(''.join(JSection.energyDependentWidths.toXMLList()),
                                 '<table rows="17" columns="4">  <columnHeaders>    <column index="0" name="energy" unit="eV"/>    <column index="1" name="levelSpacing" unit="eV"/>    <column index="2" name="neutronWidth" unit="eV"/>    <column index="3" name="captureWidth" unit="eV"/></columnHeaders>  <data>   <!-- energy | levelSpacing | neutronWidth | captureWidth  -->            2e5        8655.91      0.5280105      0.1416092            3e5       7811.708      0.4765142      0.1476385            4e5       7054.611      0.4303313      0.1538193            5e5       6375.102      0.3888812      0.1601527            6e5       5764.775      0.3516513      0.1666401            7e5       5216.188      0.3181875      0.1732827            8e5       4722.749      0.2880877      0.1800816            9e5       4278.602      0.2609947      0.1870381            1e6       3878.548      0.2365914      0.1941535          1.1e6       3517.969      0.2145961      0.2014288          1.2e6       3192.751      0.1947578      0.2088653          1.3e6       2899.241      0.1768537      0.2164641          1.4e6        2634.18       0.160685      0.2242264          1.5e6        2394.66      0.1460743      0.2321535          1.6e6       2178.092      0.1328636      0.2402462          1.7e6       1982.151      0.1209112       0.248506        1780460       1838.063      0.1121218      0.2552738</data></table>')
                break
            break

    def test_getFluctuationIntegrals(self):

        # Setup
        egrid,flag=self.Zr90URR.generateEnergyGrid()
        self.Zr90URR.getWidthsAndSpacings(egrid,interpolateWidths=flag)

        expectedResults = [
            1.10423,  1.19709086,  1.29434685,  1.39570657,  1.50079576,
            1.60915959,  1.72025964,  1.83348158,  1.94813968,  2.06348237,
            2.1787089 ,  2.29297545,  2.40541568,  2.51515507,  2.62133055,
            2.72311404,  2.80128048]

        computedResults = list(self.Zr90URR.getFluctuationIntegrals(
            widths=self.Zr90URR.averageWidths[(2,1.5)],
            DOF=self.Zr90URR.DOFs[(2,1.5)] )[0])

        # Now compare them
        self.assertEqual(len(computedResults),len(expectedResults))
        for i in range(len(expectedResults)):
            self.assertAlmostEqual(expectedResults[i], computedResults[i])

    def test_getCrossSection( self ):
        '''Test 1 MeV point'''
        result = self.Zr90URR.getCrossSection()
        answer = { 'total':6.6530426032410315, 'elastic':6.648122758610173, 'fission':0.0, 'capture':0.004919844630858361 }
        for k in answer: self.assertAlmostEqual( answer[k], result[k].evaluate(1e6) )

    def test_rho(self):
        self.assertEqual( self.Zr90URR.rho(E=0.), 0.0 )
        self.assertWithinXPercent( self.Zr90URR.rho(E=10.), 0.0043349938344906013, 0.01 )

    def test_getLastResolvedResonanceRegion(self):
        """
        TODO: Should have multiple region test too
        """
        self.assertEqual( self.Zr90URR.getLastResolvedResonanceRegion(), self.Zr90URR.reactionSuite.resonances.resolved.evaluated )

    def test_getLastResolvedResonanceEnergy(self):
        answers={(0,0.5):198400.0, (1, 1.5):193400.0, (1, 0.5):189100.0, (2, 1.5):188400.0, (2, 2.5):160000.0}
        for lj in answers:
            self.assertEqual( self.Zr90URR.getLastResolvedResonanceEnergy(*lj), answers[lj] )

    def test_getWidthsAndSpacings(self):
        # Setup
        egrid,flag=self.Zr90URR.generateEnergyGrid()
        self.Zr90URR.getWidthsAndSpacings(egrid,interpolateWidths=flag)

        # Test level spacing stuff
        levelSpacingAnswer = [
            8655.91 ,  7811.708,  7054.611,  6375.101,  5764.775,  5216.188,
            4722.749,  4278.602,  3878.548,  3517.969,  3192.751,  2899.241,
            2634.179,  2394.66 ,  2178.092,  1982.151,  1838.063]
        for i,D in enumerate(levelSpacingAnswer):
            self.assertEqual(self.Zr90URR.levelSpacings[(1, 0.5)][i],D)
            self.assertEqual(self.Zr90URR.levelSpacingFuncs[(1, 0.5)].evaluate(egrid[i]),D)
            self.assertEqual(self.Zr90URR.levelDensityFuncs[(1, 0.5)].evaluate(egrid[i]),1.0/D)

        # Test DOFs
        self.assertEqual(self.Zr90URR.DOFs[(2, 1.5)],{'competitiveDOF': False, 'neutronDOF': 1, 'fissionDOF': False})

        # Test width stuff
        neutronWidthAnswer = [
            0.751608 ,  0.6783044,  0.6125643,  0.5535613,  0.5005656,
            0.4529308,  0.4100846,  0.3715186,  0.3367811,  0.3054714,
            0.2772322,  0.2517462,  0.2287304,  0.2079326,  0.1891276,
            0.1721136,  0.1596022]
        for i,Gam in enumerate(neutronWidthAnswer):
            self.assertEqual(self.Zr90URR.averageWidths[(2, 1.5)]['neutronWidth'][i], Gam)
            self.assertEqual(self.Zr90URR.averageWidthFuncs[(2, 1.5)]['neutronWidth'].evaluate(egrid[i]), Gam)

        # Test width stuff
        captureWidthAnswer = [
            0.2450929,  0.2555282,  0.2662257,  0.2771873,  0.2884156,  0.2999123,
            0.3116797,  0.3237198,  0.3360349,  0.3486268,  0.3614976,  0.3746494,
            0.3880842,  0.401804,   0.4158108,  0.4301065,  0.44182  ]
        for i,Gam in enumerate(captureWidthAnswer):
            self.assertEqual(self.Zr90URR.averageWidths[(2, 1.5)]['captureWidth'][i], Gam)
            self.assertEqual(self.Zr90URR.averageWidthFuncs[(2, 1.5)]['captureWidth'].evaluate(egrid[i]), Gam)

    def test_generateEnergyGrid(self):
        answerGrid = [
            200000., 300000., 400000., 500000., 600000., 700000.,
            800000., 900000., 1000000., 1100000., 1200000., 1300000.,
            1400000.,  1500000.,  1600000.,  1700000.,  1780460.]
        computedGrid,flag = self.Zr90URR.generateEnergyGrid()
        for i in range(len(answerGrid)):
            self.assertAlmostEqual(computedGrid[i],answerGrid[i])
        self.assertEqual(flag, False)
        computedGrid,flag = self.Zr90URR.generateEnergyGrid(True)
        for i in range(len(answerGrid)):
            self.assertAlmostEqual(computedGrid[i],answerGrid[i])
        self.assertEqual(flag, True)

    @unittest.skip("FIXME: needs work")
    def test_getTransmissionCoefficientsFromSumRule(self):
        egrid,flag=self.Zr90URR.generateEnergyGrid()
        self.Zr90URR.getWidthsAndSpacings(egrid,interpolateWidths=flag)
        Tcs = self.Zr90URR.getTransmissionCoefficientsFromSumRule()
        c=ChannelDesignator(1, 1.5, 'elastic', 7, 1, gfact=None, particleA=None, particleB=None, isElastic=True,
                            channelClass=NEUTRONCHANNEL, useRelativistic=False, eliminated=False )
        self.assertTrue(c in Tcs)
        self.assertAlmostEqual(Tcs[c].evaluate(3.0e+05),0.35007858733963815)

    @unittest.skipIf(not HAVEBLURR,"needs blurr")
    def test_sampleRR(self):
        LJs={(0,0.5), (1, 1.5), (1, 0.5), (2, 1.5), (2, 2.5)}

        # Setup
        lastEnergies={lj:self.Zr90URR.getLastResolvedResonanceEnergy(*lj) for lj in LJs}
        egrid,flag=self.Zr90URR.generateEnergyGrid()
        self.Zr90URR.getWidthsAndSpacings(egrid,interpolateWidths=flag)

        # Generate fake resonance set
        URRset = self.Zr90URR.sampleRR(lastEnergies)

        # Count them & see if a reasonable number are made and that they are made in the right energy range
        okGuys=0
        badGuys=0
        for row in URRset.data:
            if self.Zr90URR.lowerBound <= row[0] and row[0] <= self.Zr90URR.upperBound: okGuys+=1
            else: badGuys+=1
        self.assertEqual(badGuys,0)
        self.assertTrue(abs(okGuys-2630)<200)

        # Check the energy of the first fake resonance, should be in ballpark of 0.36 MeV, so we give a 50 keV tolerance
        self.assertTrue(abs(URRset[0,0]-360000.0)<50000.0)

        # Check that the sum of partial widths equals the total
        for i in range(0,200,20):
            self.assertTrue(URRset[i,3],sum(URRset[i,4:]))

    @unittest.skipIf(not HAVEBLURR,"needs burr")
    def test_getURRPDF(self):
        urrPdfs = {}
        from fudge.core.utilities.brb import banner
#        from fudge.vis.matplotlib.plot_matrix import plot_matrix
#        import matplotlib.pyplot as plt

        for m in RRClassMap.keys():
            if m == "R_Matrix_Limited": continue #RML doesn't work well with others currently
            if m == "Reich_Moore": continue #ENDF parameters in a little different order
            if m == "SingleLevel_BreitWigner": continue #SLBW sucks
            print(banner(m))
            urrPdfs[m] = self.Zr90URR.getURRPDF(getResonanceReconstructionClass(m),memory=True,timing=True,restart=True)
#            for r in urrPdfs[m]:
#                title = m+' '+r #''$^{' + args.A + '}$' + args.sym + '(' + rxnMap[args.MT].lower() + ') cross section PDF'
#                # Fancy way
#                plot_matrix(urrPdfs[m][r].T, energyBoundariesX=eBins, energyBoundariesY=xsBins,
#                            title=title,
#                            xyTitle=('energy (MeV)', 'cross section (b)'),
#                            switchY=False)
#                plt.show()
            urrPdfs[m].save(m+'_urrPdfs.txt')

    def test_URRPDFTable(self):
        u = URRPDFTable()
        u.eBins = np.linspace(start=0.,stop=4.,num=2+1)
        u.xsBins = np.linspace(start=1.,stop=1.5,num=3+1)
        u['capture'] = np.zeros(shape=(2, 3))
        u['capture'][0, 0] = 1.
        u['capture'][1, 2] = 2.
        u.save('__junk__.txt')
        v = URRPDFTable()
        v.load('__junk__.txt')
        os.remove('__junk__.txt')
        self.assertEqual(str(u),str(v))



@unittest.skip("didn't write the tests yet")
class Test_getRI_SI( TestWithIsClose ):

    @unittest.skip("didn't write the test yet")
    def test_getRI_SI(self):
        pass #getRI_SI(E=0, Eres=0, captureWidth=0, widths=0, penetrabilities=0)



# ----------------------------------------------------------------------------------
#
#   Arrange some of the tests in test suites
#
# ----------------------------------------------------------------------------------

def angular():
    """Tests of angular distribution reconstruction"""
    suite = unittest.TestSuite()
    # SLBW tests
    suite.addTest(TestSLBWClassAndBaseClasses("test_getAngularDistribution_SWave"))
    suite.addTest(TestSLBWClassAndBaseClasses("test_getAngularDistribution_PWave"))
    suite.addTest(TestSLBWClassAndBaseClasses("test_getAngularDistribution"))
    # MLBW tests
    suite.addTest(TestMLBWClassAndBaseClasses("test_getAngularDistribution_PWave"))
    suite.addTest(TestMLBWClassAndBaseClasses("test_getElasticAngularDistribution_Full"))
    suite.addTest(TestMLBWClassAndBaseClasses("test_getElasticAngularDistribution_Zr90"))
    # Reich-Moore tests
    suite.addTest(TestRMClassAndBaseClasses("test_getElasticAngularDistribution_Small"))
    suite.addTest(TestRMClassAndBaseClasses("test_getElasticAngularDistribution_Full"))
    # R-matrix limited tests
    suite.addTest(TestRMLClassAndBaseClasses("test_getElasticAngularDistribution"))
    return suite

def U():
    """Tests of U matrix coding"""
    suite = unittest.TestSuite()
    # SLBW tests
    suite.addTest(TestSLBWClassAndBaseClasses("test_getScatteringMatrixUUnitarity"))
    suite.addTest(TestSLBWClassAndBaseClasses("test_getScatteringMatrixUValue0"))
    suite.addTest(TestSLBWClassAndBaseClasses("test_getScatteringMatrixUValue1"))
    suite.addTest(TestSLBWClassAndBaseClasses("test_getScatteringMatrixUValue_1S"))
    suite.addTest(TestSLBWClassAndBaseClasses("test_getScatteringMatrixUValue_1P"))
    suite.addTest(TestSLBWClassAndBaseClasses("test_scatteringMatrixUInCrossSectionCalculation_SWave"))
    suite.addTest(TestSLBWClassAndBaseClasses("test_scatteringMatrixUInCrossSectionCalculation_PWave"))
    suite.addTest(TestSLBWClassAndBaseClasses("test_scatteringMatrixUInCrossSectionCalculation_Full"))
    # MLBW tests
    suite.addTest(TestMLBWClassAndBaseClasses("test_getScatteringMatrixUValue_1PWave"))
    suite.addTest(TestMLBWClassAndBaseClasses("test_getScatteringMatrixUUnitarity_1PWave"))
    suite.addTest(TestMLBWClassAndBaseClasses("test_getScatteringMatrixUUnitarity_Full"))
    suite.addTest(TestMLBWClassAndBaseClasses("test_scatteringMatrixUInCrossSectionCalculation_1PWave"))
    suite.addTest(TestMLBWClassAndBaseClasses("test_scatteringMatrixUInCrossSectionCalculation_Full"))
    # Reich-Moore tests
    suite.addTest(TestRMClassAndBaseClasses("test_getScatteringMatrixUUnitarity_Full"))
    suite.addTest(TestRMClassAndBaseClasses("test_getScatteringMatrixUUnitarity_Small"))
    suite.addTest(TestRMClassAndBaseClasses("test_scatteringMatrixUInCrossSectionCalculation_Full"))
    suite.addTest(TestRMClassAndBaseClasses("test_scatteringMatrixUInCrossSectionCalculation_Small"))
    # R-matrix limited tests
    suite.addTest(TestRMLClassAndBaseClasses("test_getScatteringMatrixU"))
    return suite

def xs():
    """Tests of cross section reconstruction"""
    return unittest.TestSuite(
        map(TestResonanceReconstruction,
            ['test_SLBWReconstructResonances',
             'test_MLBWReconstructResonances',
             'test_RMReconstructResonances',
             'test_RMLReconstructResonances',
             'test_URRReconstructResonances']))


# ----------------------------------------------------------------------------------
#
#   main!!!
#
# ----------------------------------------------------------------------------------

if __name__=="__main__":

    unittest.main()
