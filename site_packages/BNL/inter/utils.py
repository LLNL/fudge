from fudge.gnds.reactionData.crossSection import component, XYs1d, upperEps
import math, os
from xData import XYs, standards
from numpy import logspace
from pqu import PQU


# -------------------------------------------------------------------------------
# Utilities
# -------------------------------------------------------------------------------

def read_evaluation(endfFile):
    '''
    Read in the evaluation
    :param endfFile: the file name.  This routine will attempt to determine the file type itself.
    :return: the evaluation as a dict of results like {'reactionSuite':the reaction suite, 'covarianceSuite': the covariances}
    '''
    if not os.path.exists( endfFile ): raise IOError( "File named "+endfFile+" doesn't exist!" )
    if open(endfFile).readline().startswith( "<?xml" ):
        from fudge.gnds import reactionSuite
        return {'reactionSuite':reactionSuite.readXML( endfFile ),'covarianceSuite':None,'info':{},'errors':[]}
    else:
        from fudge.legacy.converting import endfFileToGNDS
        return endfFileToGNDS.endfFileToGNDS( endfFile, toStdOut=False, skipBadData=True, continuumSpectraFix=True )


def equal_lethargy_bins(numBins,domainMin=1e-5,domainMax=20.0e6,reverse=False):
    seq=logspace(start=math.log10(domainMin),stop=math.log10(domainMax),num=numBins)
    if reverse: seq.reverse()
    return seq


def convert_units_to_energy(T,energyUnits='eV'):
    if type(T)==str: T=PQU.PQU(T)
    if not isinstance(T,PQU.PQU): raise TypeError( "argument to convert_units_to_energy must be a PQU, got type=%s"%str(type(T)))
    if T.isEnergy(): return T.inUnitsOf(energyUnits)
    elif T.isTemperature(): return (T.inUnitsOf('K')*PQU.PQU(8.6173324e-5, 'eV/K')).inUnitsOf(energyUnits)
#    elif T.isMass(): return
#    elif T.isLength(): return
#    elif T.isTime(): return
    else: raise Exception("Cannot convert %s to energy units" %str(T))


def function_to_XYs( func, fpars,
                     Egrid=equal_lethargy_bins(1000),
                     domainUnit='eV', domainName='energy_in', rangeUnit='1/eV', rangeName='flux',
                     accuracy=upperEps ):
    '''
    Helper function to convert a user-created function (OK, one of the spectra below) into an XYs instance
    that can be integrated, grouped, whatever.  We pre-defined a energy grid (in eV) that should work well
    even for pathological "spectra" like the problematic 1/E for the resonance integral.
    '''
    return XYs.XYs1d.createFromFunction(
            XYs1d.defaultAxes( labelsUnits={
                XYs.yAxisIndex : ( rangeName, rangeUnit ),
                XYs.xAxisIndex : ( domainName, domainUnit ) } ),
            Xs=Egrid,
            func=func,
            parameters=fpars,
            accuracy=accuracy,
            biSectionMax=20,
            checkForRoots = False,
            infill = 1,
            safeDivide = 1 )


def grouped_values_to_XYs( groupBdries, valueList,
                           domainUnit='eV', domainName='energy_in', rangeUnit='1/eV', rangeName='flux',
                           accuracy=upperEps, domainMin=1e-5, domainMax=20.0):
    if len(groupBdries) != len(valueList)+1:
        raise ValueError("Group boundries and value lists have incompatable lengths: len(bdries)=%i, len(vals)=%i"%(len(groupBdries),len(valueList)))
    curve = XYs.XYs1d(
        data=[groupBdries, [valueList[0]]+valueList],
        dataForm="xsandys",
        interpolation=standards.interpolation.flatToken,
        axes=XYs1d.defaultAxes( labelsUnits={
                XYs.yAxisIndex : ( rangeName, rangeUnit ),
                XYs.xAxisIndex : ( domainName, domainUnit ) } ) )
    return curve

def check_is_cross_section( x ):
    if not isinstance( x, component ):
        raise TypeError( "Not instance of fudge.gnds.reactionData.crossSection.component" )


# -------------------------------------------------------------------------------
# Unit tests
# -------------------------------------------------------------------------------
if __name__=="__main__":
    import unittest
    from xData.isclose import isclose
    from .metrics import DEFAULTONEOVEREXYs

    class Test_With_ApproxDiff(unittest.TestCase):

        def assertIsClose(self, a, b, abs_tol=1e-9, rel_tol=1e-8, method='average'):
            if isinstance(a,PQU.PQU) or isinstance(b,PQU.PQU):
                if not ( isinstance(a,PQU.PQU) and isinstance(b,PQU.PQU) ):
                    raise TypeError("Cannot mix types, %s with %s"%(str(type(a)),str(type(b))))
                a=float(a.inUnitsOf(b.unit))
                b=float(b)
            delta = abs(a-b)
            absave = abs(a+b)/2.0
            if absave != 0.0: reldelta = delta/absave
            else: reldelta = 0.0
            self.assertTrue(isclose( a, b, abs_tol=abs_tol, rel_tol=rel_tol, method=method),
                            "%s is not close to %s, absdiff is %s, reldiff is %s"%( str(a), str(b), str(delta), str(reldelta)) )


    class Test_Utilities(Test_With_ApproxDiff):

        def test_equal_lethargy_bins(self):
            self.maxDiff=None
            answer = [1.0e-05, 0.00037606030930863936, 0.014142135623730951, 0.53182958969449889, 20.0]
            calc=equal_lethargy_bins(5,domainMax=20.0)
            for i in range(len(calc)):
                self.assertAlmostEqual(answer[i],calc[i],7)

        def test_function_to_XYs(self):
            def oneOverEFunc(E,*args):
               return 1.0/E
            x = function_to_XYs(oneOverEFunc,[])
            for E in [1e-4, 2.3e-3, 33e-2, 0.1, 19.0]:
                self.assertIsClose(x.evaluate(E),1.0/E)
            def myExp(E,*args):
                return math.exp(-E/1e3)
            y = function_to_XYs(myExp,[])
            for E in [33e-2, 0.1, 20.]:
                self.assertIsClose(y.evaluate(E),math.exp(-E/1e3))

        def test_grouped_values_to_XYs(self):
            self.assertIsClose(0.0,0.0)

        def test_integrateTwoFunctions(self):
            def myExp1(E,*args):
                return math.exp(-E)
            def myExp2(E,*args):
                return math.exp(-E/2.0)
            a=function_to_XYs(myExp1,[]).integrateTwoFunctions(function_to_XYs(myExp2,[])).inUnitsOf('1/eV')
            b=PQU.PQU(0.666657,'1/eV')
            self.assertIsClose(a,b,abs_tol=5e-7)

        def test_one_over_E(self):
            """
            Mathematica says the answer is \int_0.5^20e6 dE/E = 17.5044
            :return:
            """
            def oneFunc(x,*args): return 1.0
            one=function_to_XYs(oneFunc,[])
            self.assertIsClose(DEFAULTONEOVEREXYs.integrateTwoFunctions(one), PQU.PQU(17.5044,'1/eV'),rel_tol=1e-6)

    unittest.main()
