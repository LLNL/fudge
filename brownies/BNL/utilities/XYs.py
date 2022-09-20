from xData import enums as xDataEnumsModule
from xData import XYs1d as XYs1dModule
from xData import axes as axesModule

from fudge.reactionData.crossSection import XYs1d, upperEps
from brownies.BNL.utilities.bins import *

# -------------------------------------------------------------------------------
# Utilities
# -------------------------------------------------------------------------------


def makeXYs(data, xUnit='eV', xLabel='energy_in', yUnit='b', yLabel="crossSection", dataForm="xys", interpolation=xDataEnumsModule.Interpolation.linlin):
    return XYs1dModule.XYs1d(
        axes=axesModule.Axes(2, labelsUnits={0: (yLabel, yUnit), 1: (xLabel, xUnit)}),
        data=data,
        dataForm=dataForm,
        interpolation=interpolation
    )


def makeHistoXYs(listOfPoints, nBins, normfactor=1.0, xLabel="Ex", xUnit='eV', yLabel="rho(Ex)", yUnit="1/eV"):
    import numpy
    histo = numpy.histogram(listOfPoints, nBins)
    return makeXYs(data=[histo[1][0:-1], histo[0] * normfactor], xLabel=xLabel, yLabel=yLabel, xUnit=xUnit, yUnit=yUnit,
                   dataForm="xsandys")


def function_to_XYs(func, fpars,
                    Egrid=equal_lethargy_bins(1000),
                    domainUnit='eV', domainName='energy_in', rangeUnit='1/eV', rangeName='flux',
                    accuracy=upperEps):
    """
    Helper function to convert a user-created function (OK, one of the spectra below) into an XYs instance
    that can be integrated, grouped, whatever.  We pre-defined a energy grid (in eV) that should work well
    even for pathological "spectra" like the problematic 1/E for the resonance integral.
    """
    return XYs1dModule.XYs1d.createFromFunction(
        XYs1d.defaultAxes(labelsUnits={
            XYs1dModule.yAxisIndex: (rangeName, rangeUnit),
            XYs1dModule.xAxisIndex: (domainName, domainUnit)}),
        Xs=Egrid,
        func=func,
        parameters=fpars,
        accuracy=accuracy,
        biSectionMax=20,
        checkForRoots=False,
        infill=1,
        safeDivide=1)


def grouped_values_to_XYs(groupBdries, valueList,
                          domainUnit='eV', domainName='energy_in', rangeUnit='1/eV', rangeName='flux'):
    if len(groupBdries) != len(valueList) + 1:
        raise ValueError("Group boundries and value lists have incompatable lengths: "
                         "len(bdries)=%i, len(vals)=%i" % (len(groupBdries), len(valueList)))
    curve = XYs1dModule.XYs1d(
        data=[groupBdries, [valueList[0]] + valueList],
        dataForm="xsandys",
        interpolation=xDataEnumsModule.Interpolation.flat,
        axes=XYs1d.defaultAxes(labelsUnits={XYs1dModule.yAxisIndex: (rangeName, rangeUnit), XYs1dModule.xAxisIndex: (domainName, domainUnit)}))

    return curve
