# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from fudge.core.math.pdf import *

# ---------------------------------------------------------------------------------
#
#   Main, used for running unit tests.
#   If you want to run something real, run __main__
#
# ---------------------------------------------------------------------------------

if __name__=="__main__":
    import unittest

    DEBUG=False
    PROFILE=False

    # ---------------------------------------------------------------------------------
    #
    #   Unit tests
    #
    # ---------------------------------------------------------------------------------

    class TestNumpy(unittest.TestCase):

        def test_average(self):
            self.assertEqual(numpy.average(range(1,5)), 2.5)

        def test_mean(self):
            self.assertEqual(numpy.mean(range(1,5)), 2.5)


    @unittest.skipIf(DEBUG,'debugging')
    class TestPDF(unittest.TestCase):

        @unittest.expectedFailure
        def test_init(self):
            UnivariatePDF()


    @unittest.skipIf(DEBUG,'debugging')
    class TestWignerDistribution(unittest.TestCase):

        def test_simple_init(self):
            x=WignerDistribution()
            self.assertIsNone(x.evaluate(-1.0))
            self.assertAlmostEqual(x.evaluate(0.0),  0.0)
            self.assertAlmostEqual(x.evaluate(0.5),  0.6453812728576525)
            self.assertAlmostEqual(x.evaluate(1.0),  0.7161859363405692)
            self.assertAlmostEqual(x.evaluate(2.0),  0.1357605281502967)
            self.assertAlmostEqual(x.evaluate(3.0),  0.004012308664112029)
            self.assertIsNone(x.evaluate(11.0))

        def test_less_simple_init(self):
            y=WignerDistribution(xScale=2.0)
            self.assertIsNone(y.evaluate(-1.0))
            self.assertAlmostEqual(y.evaluate(0.0),  0.0)
            self.assertAlmostEqual(y.evaluate(0.5),  0.18694399076760995)
            self.assertAlmostEqual(y.evaluate(1.0),  0.3226906364288262)
            self.assertAlmostEqual(y.evaluate(2.0),  0.3580929681702846)
            self.assertAlmostEqual(y.evaluate(3.0),  0.2012423783795464)
            self.assertIsNone(y.evaluate(21.0))

        def test_fancy_init(self):
            x=WignerDistribution(xUnit='MeV',xScale=PQU("24.3 keV"),xLabel="D", yLabel="P(D)")
            self.assertIsNone(x.evaluate(-1.0))
            self.assertAlmostEqual(x.evaluate(0.0),  0.0, 6)
            self.assertAlmostEqual(x.evaluate(0.01), 23.28856576595401)
            self.assertAlmostEqual(x.evaluate(0.02), 31.252084478962146)
            self.assertAlmostEqual(x.evaluate(0.04), 12.668756479803134)
            self.assertAlmostEqual(x.evaluate(0.06), 1.3290263130187128)
            self.assertIsNone(x.evaluate(0.3))

        def test_getValueRawPDF(self):
            self.assertAlmostEqual(WignerDistribution.getValueRawPDF(-1.0,{'xScale':1.0}), 0.0)
            self.assertAlmostEqual(WignerDistribution.getValueRawPDF(0.0,{'xScale':1.0}),  0.0)
            self.assertAlmostEqual(WignerDistribution.getValueRawPDF(0.5,{'xScale':1.0}),  0.6453812728576525)
            self.assertAlmostEqual(WignerDistribution.getValueRawPDF(1.0,{'xScale':1.0}),  0.7161859363405692)
            self.assertAlmostEqual(WignerDistribution.getValueRawPDF(2.0,{'xScale':1.0}),  0.1357605281502967)
            self.assertAlmostEqual(WignerDistribution.getValueRawPDF(3.0,{'xScale':1.0}),  0.004012308664112029)

        def test_getValueRawPDF_fancy(self):
            self.assertAlmostEqual(WignerDistribution.getValueRawPDF(-1.0,{'xScale':2.0}), 0.0)
            self.assertAlmostEqual(WignerDistribution.getValueRawPDF(0.0,{'xScale':2.0}),  0.0)
            self.assertAlmostEqual(WignerDistribution.getValueRawPDF(0.5,{'xScale':2.0}),  0.18694399076760995)
            self.assertAlmostEqual(WignerDistribution.getValueRawPDF(1.0,{'xScale':2.0}),  0.3226906364288262)
            self.assertAlmostEqual(WignerDistribution.getValueRawPDF(2.0,{'xScale':2.0}),  0.3580929681702846)
            self.assertAlmostEqual(WignerDistribution.getValueRawPDF(3.0,{'xScale':2.0}),  0.2012423783795464)

        def test_invCDF(self):
            x=WignerDistribution(xUnit='MeV',xScale=PQU("24.3 keV"),xLabel="D", yLabel="P(D)")
            x.setCDFInverse()
#            x.invcdf.plot()
            sample=x.drawSample(20000)
            for i,xx in enumerate(sample):
                if xx is None: print( "bad sample:", i, xx )
                self.assertIsNotNone( xx )
            self.assertAlmostEqual( numpy.mean( sample ), x.mean, 3 )
            self.assertAlmostEqual( numpy.std( sample ), x.stddev, 3 )


    class TestGOEDistribution(unittest.TestCase):

        def test_simple_init(self):
            x=GOEDistribution()
            self.assertAlmostEqual(x.evaluate(-1.0), 0.0)
            self.assertAlmostEqual(x.evaluate(0.0),  2.0/math.pi)
            self.assertAlmostEqual(x.evaluate(0.8),  (2.0/math.pi)*math.sqrt(1.0-0.8*0.8) )
            self.assertAlmostEqual(x.evaluate(1.0),  0.0)
            self.assertAlmostEqual(x.evaluate(2.0),  0.0)

        def test_less_simple_init(self):
            x=GOEDistribution(a=2.0)
            self.assertAlmostEqual(x.evaluate(-2.0), 0.0)
            self.assertAlmostEqual(x.evaluate(0.0),  2.0/math.pi/2.0)
            self.assertAlmostEqual(x.evaluate(1.6),  (2.0/math.pi/4.0)*math.sqrt(4.0-1.6*1.6) )
            self.assertAlmostEqual(x.evaluate(2.0),  0.0)
            self.assertAlmostEqual(x.evaluate(4.0),  0.0)

        def test_fancy_init(self):
            x=GOEDistribution(xUnit='MeV',a=PQU("24.3 keV"))
            self.assertAlmostEqual(x.evaluate(0.0),  26.198344541875773)
            self.assertAlmostEqual(x.evaluate(0.01), 23.877159126363704)
            self.assertAlmostEqual(x.evaluate(0.02), 14.88002515374847)
            self.assertAlmostEqual(x.evaluate(0.04), 0.0)
            self.assertAlmostEqual(x.evaluate(0.06), 0.0)

        def test_getValueRawPDF(self):
            self.assertAlmostEqual(GOEDistribution.getValueRawPDF(-1.0,{'a':1.0}), 0.0)
            self.assertAlmostEqual(GOEDistribution.getValueRawPDF(0.0,{'a':1.0}),  2.0/math.pi)
            self.assertAlmostEqual(GOEDistribution.getValueRawPDF(0.8,{'a':1.0}),  (2.0/math.pi)*math.sqrt(1.0-0.8*0.8))
            self.assertAlmostEqual(GOEDistribution.getValueRawPDF(1.0,{'a':1.0}),  0.0)
            self.assertAlmostEqual(GOEDistribution.getValueRawPDF(2.0,{'a':1.0}),  0.0)
            self.assertAlmostEqual(GOEDistribution.getValueRawPDF(3.0,{'a':1.0}),  0.0)

        def test_getValueRawPDF_fancy(self):
            self.assertAlmostEqual(GOEDistribution.getValueRawPDF(-2.0,{'a':2.0}), 0.0)
            self.assertAlmostEqual(GOEDistribution.getValueRawPDF(0.0,{'a':2.0}),  2.0/math.pi/2.0)
            self.assertAlmostEqual(GOEDistribution.getValueRawPDF(1.6,{'a':2.0}),  (2.0/math.pi/4.0)*math.sqrt(4.0-1.6*1.6))
            self.assertAlmostEqual(GOEDistribution.getValueRawPDF(2.0,{'a':2.0}),  0.0)
            self.assertAlmostEqual(GOEDistribution.getValueRawPDF(3.0,{'a':2.0}),  0.0)

        def test_invCDF(self):
            x=GOEDistribution(xUnit='MeV',a=PQU("24.3 keV"), xLabel="E", yLabel="P(E)")
#            x.cdf.plot()
            x.setCDFInverse()
#            x.invcdf.plot()
            sample=x.drawSample(2000)
            print( sample )
            for i,xx in enumerate(sample):
                if xx is None: print( "bad sample:",i,xx )
                self.assertIsNotNone( xx )
            self.assertAlmostEqual( numpy.std( sample, dtype=float ), x.stddev, 3 )
            self.assertAlmostEqual( numpy.mean( sample, dtype=float ), x.mean, 3 )


    @unittest.skipIf(DEBUG,'debugging')
    class TestBrodyDistribution(unittest.TestCase):

        def test_simple_init(self):
            x=BrodyDistribution()
            self.assertIsNone(x.evaluate(-1.0))
            self.assertAlmostEqual(x.evaluate(0.001),  0.04068434485353072)
            self.assertAlmostEqual(x.evaluate(0.5),  0.6717747679027428)
            self.assertAlmostEqual(x.evaluate(1.0),  0.5456750040739453)
            self.assertAlmostEqual(x.evaluate(2.0),  0.1608239521290362)
            self.assertAlmostEqual(x.evaluate(3.0),   0.025846797494801642)
            self.assertIsNone(x.evaluate(11.0))
            #x.plot()

        def test_fancy_init_Wigner(self):
            """Should recover Wigner Distribution results"""
            x=BrodyDistribution(xUnit='MeV', xScale=PQU("24.3 keV"), w=1.0, xLabel="D", yLabel="P(D)")
            self.assertIsNone(x.evaluate(-1.0))
            self.assertAlmostEqual(x.evaluate(0.0),  0.0, 6)
            self.assertAlmostEqual(x.evaluate(0.01), 23.28856576595401)
            self.assertAlmostEqual(x.evaluate(0.02), 31.252084478962146)
            self.assertAlmostEqual(x.evaluate(0.04), 12.668756479803134)
            self.assertAlmostEqual(x.evaluate(0.06), 1.3290263130187128)
            self.assertIsNone(x.evaluate(0.3))

        def test_fancy_init_Poisson(self):
            x=BrodyDistribution(xUnit='MeV', xScale=PQU("24.3 keV"), w=0.0, xLabel="D", yLabel="P(D)")
#            x.plot()
            self.assertIsNone(x.evaluate(-1.0))
            self.assertAlmostEqual(x.evaluate(0.0001),  40.98326054045971)
            self.assertAlmostEqual(x.evaluate(0.01), 27.269157473159645)
            self.assertAlmostEqual(x.evaluate(0.02), 18.069648752690604)
            self.assertAlmostEqual(x.evaluate(0.04), 7.934246553488046)
            self.assertAlmostEqual(x.evaluate(0.06), 3.483867849515246)
            self.assertIsNone(x.evaluate(0.3))

        def test_getValueRawPDF_Wigner(self):
            """Should recover Wigner Distribution results for w=1"""
            # Wigner values
            self.assertAlmostEqual(BrodyDistribution.getValueRawPDF(-1.0,{'xScale':1.0, 'w':1.0}), 0.0)
            self.assertAlmostEqual(BrodyDistribution.getValueRawPDF(0.0,{'xScale':1.0, 'w':1.0}),  0.0, 6)
            self.assertAlmostEqual(BrodyDistribution.getValueRawPDF(0.5,{'xScale':1.0, 'w':1.0}),  0.6453812728576525)
            self.assertAlmostEqual(BrodyDistribution.getValueRawPDF(1.0,{'xScale':1.0, 'w':1.0}),  0.7161859363405692)
            self.assertAlmostEqual(BrodyDistribution.getValueRawPDF(2.0,{'xScale':1.0, 'w':1.0}),  0.1357605281502967)
            self.assertAlmostEqual(BrodyDistribution.getValueRawPDF(3.0,{'xScale':1.0, 'w':1.0}),  0.004012308664112029)

        def test_getValueRawPDF_Poisson(self):
            """Should recover Poisson Distribution results for w=0"""
            # Poisson values
            self.assertAlmostEqual(BrodyDistribution.getValueRawPDF(-1.0,{'xScale':1.0, 'w':0.0}),    0.0)
            self.assertAlmostEqual(BrodyDistribution.getValueRawPDF(0.0,{'xScale':1.0, 'w':0.0}),     1.0)
            self.assertAlmostEqual(BrodyDistribution.getValueRawPDF(0.0001,{'xScale':1.0, 'w':0.0}),  0.9999000049998333)
            self.assertAlmostEqual(BrodyDistribution.getValueRawPDF(0.5,{'xScale':1.0, 'w':0.0}),     0.6065306597126334)
            self.assertAlmostEqual(BrodyDistribution.getValueRawPDF(1.0,{'xScale':1.0, 'w':0.0}),     0.36787944117144233)
            self.assertAlmostEqual(BrodyDistribution.getValueRawPDF(2.0,{'xScale':1.0, 'w':0.0}),     0.1353352832366127)
            self.assertAlmostEqual(BrodyDistribution.getValueRawPDF(3.0,{'xScale':1.0, 'w':0.0}),     0.049787068367863944)


    @unittest.skipIf(DEBUG,'debugging')
    class TestPoissonDistribution(unittest.TestCase):

        def test_fancy_init(self):
            x=PoissonDistribution(xUnit='MeV', xScale=PQU("24.3 keV"), xLabel="D", yLabel="P(D)")
#            x.plot()
            self.assertIsNone(x.evaluate(-1.0))
            self.assertAlmostEqual(x.evaluate(0.0001),  40.98326054045971)
            self.assertAlmostEqual(x.evaluate(0.01), 27.269157473159645)
            self.assertAlmostEqual(x.evaluate(0.02), 18.069648752690604)
            self.assertAlmostEqual(x.evaluate(0.04), 7.934246553488046)
            self.assertAlmostEqual(x.evaluate(0.06), 3.483867849515246)
            self.assertIsNone(x.evaluate(0.3))

        def test_getValueRawPDF(self):
            self.assertAlmostEqual(PoissonDistribution.getValueRawPDF(-1.0,{'xScale':1.0}),    0.0)
            self.assertAlmostEqual(PoissonDistribution.getValueRawPDF(0.0,{'xScale':1.0}),     1.0)
            self.assertAlmostEqual(PoissonDistribution.getValueRawPDF(0.0001,{'xScale':1.0}),  0.9999000049998333 )
            self.assertAlmostEqual(PoissonDistribution.getValueRawPDF(0.5,{'xScale':1.0}),     0.6065306597126334)
            self.assertAlmostEqual(PoissonDistribution.getValueRawPDF(1.0,{'xScale':1.0}),     0.36787944117144233)
            self.assertAlmostEqual(PoissonDistribution.getValueRawPDF(2.0,{'xScale':1.0}),     0.1353352832366127)
            self.assertAlmostEqual(PoissonDistribution.getValueRawPDF(3.0,{'xScale':1.0}),     0.049787068367863944)


    @unittest.skipIf(DEBUG,'debugging')
    class TestUniformDistribution(unittest.TestCase):

        def test_plain_init(self):
            x=UniformDistribution()
            self.assertIsNone(x.evaluate(-11.0))
            self.assertAlmostEqual(x.evaluate(-0.1), 0.0)
            self.assertAlmostEqual(x.evaluate(0.0),  1.0)
            self.assertAlmostEqual(x.evaluate(0.5),  1.0)
            self.assertAlmostEqual(x.evaluate(1.0),  1.0)
            self.assertAlmostEqual(x.evaluate(1.1),  0.0)
            self.assertIsNone(x.evaluate(10.1))

        def test_fancy_init(self):
            xMin = 0.0
            xMax = 20.0
            p = 1.0/(xMax-xMin)
            x=UniformDistribution(xUnit='MeV', supportMin=xMin, supportMax=xMax, domainMin=-1.0, domainMax=21.0)
            self.assertIsNone(x.evaluate(-10.0))
            self.assertAlmostEqual(x.evaluate(-0.1),  0.0)
            self.assertAlmostEqual(x.evaluate(0.0001),  p)
            self.assertAlmostEqual(x.evaluate(0.5),  p)
            self.assertAlmostEqual(x.evaluate(1.0),  p)
            self.assertAlmostEqual(x.evaluate(19.0), p)
            self.assertAlmostEqual(x.evaluate(21.0), 0.0)
            self.assertIsNone(x.evaluate(31.1))

        @unittest.expectedFailure
        def test_failed_init(self):
            UniformDistribution(supportMax=-1.0,supportMin=0.0)

        def test_getValueRawPDF(self):
            self.assertAlmostEqual(UniformDistribution.getValueRawPDF(-1.0,{}), 0.0)
            self.assertAlmostEqual(UniformDistribution.getValueRawPDF(0.0,{}),  1.0)
            self.assertAlmostEqual(UniformDistribution.getValueRawPDF(0.5,{}),  1.0)
            self.assertAlmostEqual(UniformDistribution.getValueRawPDF(1.0,{}),  1.0)
            self.assertAlmostEqual(UniformDistribution.getValueRawPDF(1.1,{}),  0.0)

        def test_mean(self):
            xMin = 0.0
            xMax = 20.0
            p = 1.0/(xMax-xMin)
            x=UniformDistribution(xUnit='MeV', supportMin=xMin, supportMax=xMax, domainMin=-1.0, domainMax=21.0)
            self.assertAlmostEqual(x.mean,1.0/(2.0*p))

        def test_kurtosis(self):
            xMin = 0.0
            xMax = 20.0
            p = 1.0/(xMax-xMin)
            x=UniformDistribution(xUnit='MeV', supportMin=xMin, supportMax=xMax, domainMin=-1.0, domainMax=21.0)
            self.assertAlmostEqual(x.kurtosis,3.0-6.0/5.0)

        def test_skew(self):
            xMin = 0.0
            xMax = 20.0
            x=UniformDistribution(xUnit='MeV', supportMin=xMin, supportMax=xMax, domainMin=-1.0, domainMax=21.0)
            self.assertAlmostEqual(x.skew,0.0)

        def test_variance(self):
            xMin = 0.0
            xMax = 20.0
            p = 1.0/(xMax-xMin)
            x=UniformDistribution(xUnit='MeV', supportMin=xMin, supportMax=xMax, domainMin=-1.0, domainMax=21.0)
            self.assertAlmostEqual(x.variance,pow(xMax-xMin,2.0)/12.0)

        def test_stddev(self):
            xMin = 0.0
            xMax = 20.0
            p = 1.0/(xMax-xMin)
            x=UniformDistribution(xUnit='MeV', supportMin=xMin, supportMax=xMax, domainMin=-1.0, domainMax=21.0)
            self.assertAlmostEqual(x.stddev,math.sqrt(pow(xMax-xMin,2.0)/12.0))

        def test_getCDF(self):
            """cdf(x)=x since pdf(x)=1.0"""
            x=UniformDistribution()
            self.assertIsNone(x.cdf.evaluate(-11.0))
            self.assertAlmostEqual(x.cdf.evaluate(-10.0),0.0)
            self.assertAlmostEqual(x.cdf.evaluate(0.0),  0.0, 6)
            self.assertAlmostEqual(x.cdf.evaluate(0.5),  0.5, 6)
            self.assertAlmostEqual(x.cdf.evaluate(1.0),  1.0, 6)
            self.assertAlmostEqual(x.cdf.evaluate(1.1),  1.0, 6)
            self.assertIsNone(x.cdf.evaluate(11.1))

        def test_invCDF(self):
            x=UniformDistribution()
            x.setCDFInverse()
#            x.invcdf.plot()
            self.assertIsNone(x.invcdf.evaluate(-11.0))
            self.assertIsNone(x.invcdf.evaluate(-10.0))
            self.assertAlmostEqual(x.invcdf.evaluate(0.0),  0.0)
            self.assertAlmostEqual(x.invcdf.evaluate(0.5),  0.5)
            self.assertAlmostEqual(x.invcdf.evaluate(1.0),  1.0)
            self.assertIsNone(x.invcdf.evaluate(1.1))
            self.assertIsNone(x.invcdf.evaluate(11.1))

        def test_drawSample(self):pass

        def test_isSampleConsistentWithPDF(self):pass


    class TestNormalDistribution(unittest.TestCase):

        def test_simple_init(self):
            x=NormalDistribution(mean=1.8,stddev=4.5)
            #x.plot()
            #x.cdf.plot()

        def test_invCDF(self):
            x=NormalDistribution(mean=1.0,stddev=2.0)
            for y in [0.0001, 0.25, 0.5, 0.75, 1.0]: # can't test 0
                self.assertAlmostEqual( x.cdf.evaluate(x.invcdf.evaluate(y)), y, 3 )

        def test_sampling(self):
            x=NormalDistribution(mean=1.8,stddev=1.5)
            #x.plot()
            #x.cdf.plot()
            #x.invcdf.plot()
            sample=x.drawSample(2000000)
            for i,xx in enumerate(sample):
                if xx is None: print( "bad sample:",i,xx )
                self.assertIsNotNone( xx )
            self.assertAlmostEqual(numpy.mean(sample), x.mean, 2)
            self.assertAlmostEqual(numpy.std(sample), x.stddev, 2)
            self.assertAlmostEqual(numpy.var(sample), x.variance, 2)

        def test_getValueRawPDF(self):
            self.assertAlmostEqual(NormalDistribution.getValueRawPDF(-1.0,{'mean':1.0, 'stddev':2.0}), 0.12098536225957168)
            self.assertAlmostEqual(NormalDistribution.getValueRawPDF( 0.0,{'mean':1.0, 'stddev':2.0}), 0.17603266338214976)


    class TestLogNormalDistribution(unittest.TestCase):

        def test_simple_init(self):
            x=LogNormalDistribution(mu=1.8,sigma=4.5)

        def test_invCDF(self):
            x=LogNormalDistribution(mu=1.0,sigma=2.0)
            for y in [0.0001, 0.25, 0.5, 0.75, 1.0]: # can't test 0
                self.assertAlmostEqual( x.cdf.evaluate(x.invcdf.evaluate(y)), y, 3 )

        def test_sampling(self):
            x=LogNormalDistribution(mu=1.8,sigma=1.5)
            #x.plot()
            #x.cdf.plot()
            #x.invcdf.plot()
            sample=x.drawSample(19000000)
            for i,xx in enumerate(sample):
                if xx is None: print( "bad sample:",i,xx )
                self.assertIsNotNone( xx )
            self.assertAlmostEqual(numpy.mean(sample), x.mean, 1)
            self.assertAlmostEqual(numpy.std(sample)/100., x.stddev/100., 1)
            self.assertAlmostEqual(numpy.var(sample)/10000., x.variance/10000., 2)

        def test_getValueRawPDF(self):
            self.assertAlmostEqual(LogNormalDistribution.getValueRawPDF(-1.0, {'mu':1.0, 'sigma':2.0}), 0.0)
            self.assertAlmostEqual(LogNormalDistribution.getValueRawPDF( 0.0, {'mu':1.0, 'sigma':2.0}), 0.0)
            self.assertAlmostEqual(LogNormalDistribution.getValueRawPDF( 4.0, {'mu':1.0, 'sigma':2.0}), 0.048946227003151078)


    @unittest.skipIf(DEBUG,'debugging')
    class TestChi2Distribution(unittest.TestCase):

        def test_simple_init(self):
            x=Chi2Distribution(1.8)
            # x.plot()
            # x.cdf.plot()

        def test_invCDF(self):
            x=Chi2Distribution(1.8)
            x.setCDFInverse()
            # x.invcdf.plot()

        def test_sampling_nu_is_1(self):
            x=Chi2Distribution(dof=1.0)
            #x.plot()
            #x.cdf.plot()
            #x.invcdf.plot()
            sample=x.drawSample(2000000)
            for i,xx in enumerate(sample):
                if xx is None: print( "bad sample:",i,xx )
                self.assertIsNotNone( xx )
            self.assertAlmostEqual( numpy.mean( sample ), x.mean, 2 )
            self.assertAlmostEqual( numpy.std( sample ), x.stddev, 2 )

        def test_sampling_nu_is_2(self):
            x=Chi2Distribution(2.0)
            sample=x.drawSample(2000000)
            for i,xx in enumerate(sample):
                if xx is None: print( "bad sample:",i,xx )
                self.assertIsNotNone( xx )
            self.assertAlmostEqual( numpy.mean( sample ), x.mean, 2 )
            self.assertAlmostEqual( numpy.std( sample ), x.stddev, 2 )

        def test_sampling_nu_is_3(self):
            x=Chi2Distribution(3.0)
            sample=x.drawSample(2000000)
            for i,xx in enumerate(sample):
                if xx is None: print( "bad sample:",i,xx )
                self.assertIsNotNone( xx )
            self.assertAlmostEqual( numpy.mean( sample ), x.mean, 2 )
            self.assertAlmostEqual( numpy.std( sample ), x.stddev, 2 )

        def test_sampling_nu_is_4(self):
            x=Chi2Distribution(4.0)
            sample=x.drawSample(200000)
            for i,xx in enumerate(sample):
                if xx is None: print( "bad sample:",i,xx )
                self.assertIsNotNone( xx )
            self.assertAlmostEqual( numpy.mean( sample ), x.mean, 2 )
            self.assertAlmostEqual( numpy.std( sample ), x.stddev, 2 )

        def test_getValueRawPDF(self):
            # For nu=2.0, our chi^2 distribution == our Poisson/exponential distribution
            self.assertAlmostEqual(Chi2Distribution.getValueRawPDF(-1.0,{'xScale':1.0, 'dof':2.0}),    0.0)
            self.assertAlmostEqual(Chi2Distribution.getValueRawPDF(0.0,{'xScale':1.0, 'dof':2.0}),     1.0)
            self.assertAlmostEqual(Chi2Distribution.getValueRawPDF(0.0001,{'xScale':1.0, 'dof':2.0}),  0.9999000049998333 )
            self.assertAlmostEqual(Chi2Distribution.getValueRawPDF(0.5,{'xScale':1.0, 'dof':2.0}),     0.6065306597126334)
            self.assertAlmostEqual(Chi2Distribution.getValueRawPDF(1.0,{'xScale':1.0, 'dof':2.0}),     0.36787944117144233)
            self.assertAlmostEqual(Chi2Distribution.getValueRawPDF(2.0,{'xScale':1.0, 'dof':2.0}),     0.1353352832366127)
            self.assertAlmostEqual(Chi2Distribution.getValueRawPDF(3.0,{'xScale':1.0, 'dof':2.0}),     0.049787068367863944)




    # Actually run the tests
    if PROFILE:
        import cProfile
        cProfile.run("unittest.main()", "burr__init__.profile")
    else:
        #try:
        #    import xmlrunner
        #    unittest.main(testRunner=xmlrunner.XMLTestRunner(output='test-results'))
        #except ImportError:
            unittest.main()
            print()
