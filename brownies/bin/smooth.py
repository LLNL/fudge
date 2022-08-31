#! /usr/bin/env python
import argparse
import copy
import os
import math
import pqu.PQU as PQU

from xData import enums as xDataEnumsModule
from xData import XYs1d as XYs1d
from xData import axes as axesModule
from xData import values as valuesModule
from xData import link as linkModule
from xData import xDataArray as arrayModule
from xData import gridded as griddedModule

from fudge.reactionData.crossSection import XYs1d as pointwise

# -------------------------------------------
# Variables ---------------------------------
# -------------------------------------------
defaultDeltaE = '200. keV'
defaultTarget = "Fe56"
defaultProjectile = 'n'
accuracy = 1e-8
ENDF_NeutronSublibrary_EMin = 1e-5  # in eV
ENDF_NeutronSublibrary_EMax = 2e7  # in eV


# -------------------------------------------
# Smoothing routines ------------------------
# -------------------------------------------
def defaultAxes(xName, xUnit, yName, yUnit):
    """
    Prepare generic axes, but with the right name and units
    :param xName: the x axis name
    :param xUnit: the PQU unit of the x axis
    :param yName: the y axis name
    :param yUnit: the PQU unit of the y axis
    :return: the axes instance, ready to be attached to an XYs1d
    """
    axes = axesModule.Axes(2)
    axes[0] = axesModule.Axis(yName, 0, yUnit)
    axes[1] = axesModule.Axis(xName, 1, xUnit)
    return axes


def make_egrid(E0, dE):
    ddE = dE/100
    _egrid = [1e-5] + \
             [E0 - i * 100.1 * ddE for i in range(3, 10)] + \
             [E0 - i * 10.1 * ddE for i in range(2, 20)] + \
             [E0 - i * 2.1 * ddE for i in range(1, 50)] + \
             [E0 - i * ddE for i in range(1, 100)] + \
             [E0 + i * ddE for i in range(101)] + \
             [E0 + i * 2.1 * ddE for i in range(1, 50)] + \
             [E0 + i * 10.1 * ddE for i in range(2, 20)] + \
             [E0 + i * 100.1 * ddE for i in range(3, 10)] + \
             [2e7]
    __egrid = []
    for x in _egrid:
        if x >= ENDF_NeutronSublibrary_EMin and x not in __egrid:
            __egrid.append(x)
    __egrid.sort()
    return __egrid


def InelasticThreshold(Ecut, Ethreshold, xsAtEcut):
    """
    Ecut in eV, this is point where we switch from averaging to the sqrt thingee
    Ethreshold in eV 
    xsAtEcut in b
    """
    print(Ecut, Ethreshold, xsAtEcut)
    __the_const = xsAtEcut / math.sqrt(Ecut - Ethreshold)

    def __sqrtThingee(E, *args):
        if E < Ethreshold:
            return 0.0
        else:
            return __the_const * math.sqrt(E - Ethreshold)
    __this_grid = copy.copy(Egrid) + [Ecut, Ethreshold]
    __this_grid.sort()
    return XYs1d.XYs1d.createFromFunction(
        defaultAxes(xName='E', xUnit='eV', yName='Sigma(E)', yUnit='b'),
        Xs=__this_grid,
        func=__sqrtThingee,
        parameters=[],
        accuracy=accuracy,
        biSectionMax=20,
        checkForRoots=False,
        infill=1,
        safeDivide=1)


def Lorentzian(E0, dE):
    def __lorentzian(E, *args):
        return 1.0 / ((E - E0) ** 2 + dE * dE)
    return XYs1d.XYs1d.createFromFunction(
        defaultAxes(xName='E', xUnit='eV', yName="Lorentzian(E)", yUnit='1/eV'),
        Xs=make_egrid(E0, dE),
        func=__lorentzian,
        parameters=[],
        accuracy=accuracy,
        biSectionMax=20,
        checkForRoots=False,
        infill=1,
        safeDivide=1)


def Gaussian(E0, dE):
    def __gaussian(E, *args):
        return math.exp(-(E - E0) ** 2 / 2.0 / dE / dE)
    return XYs1d.XYs1d.createFromFunction(
        defaultAxes(xName='E', xUnit='eV', yName="Gaussian(E)", yUnit='1/eV'),
        Xs=make_egrid(E0, dE),
        func=__gaussian,
        parameters=[],
        accuracy=accuracy,
        biSectionMax=20,
        checkForRoots=False,
        infill=1,
        safeDivide=1)


def WindowFunction(E0, dE):
    def __window(E, *args):
        if abs(E-E0) <= dE:
            return 0.5 / dE
        else:
            return 0.0
    return XYs1d.XYs1d.createFromFunction(
        defaultAxes(xName='E', xUnit='eV', yName="Window(E)", yUnit='1/eV'),
        Xs=make_egrid(E0, dE),
        func=__window,
        parameters=[],
        accuracy=accuracy,
        biSectionMax=20,
        checkForRoots=False,
        infill=1,
        safeDivide=1)


profileOptions = {'lorentzian': Lorentzian, 'window': WindowFunction, 'gaussian': Gaussian}


def smooth(xs, _cov, profile=None, _dE=defaultDeltaE, _EThreshold=0.0, _nE=201, _Emax=None, _dataType='crossSection',
           verbose=True):

    _Emin = ENDF_NeutronSublibrary_EMin
    if _Emax is None:
        _Emax = ENDF_NeutronSublibrary_EMax
    _Emin_asPQU = PQU.PQU(str(_Emin) + ' eV')
    _Emax_asPQU = PQU.PQU(str(_Emax) + ' eV')
    result = []
    eStep = (_Emax - _Emin) / (_nE - 1)
    norm = None

    if verbose:
        print('#', profile.capitalize(), _dE, 'eV')

    # Smooth over entire energy range
    for i in range(_nE):
        E0 = max(_Emin, eStep * i)
        s = profileOptions[profile](E0, _dE)

        # Get the norm of the averaging function since none of them are normalized properly
        if norm is None:
            norm = s.integrate(_Emin, _Emax)
            print('# Smoothing function norm is ', str(norm))

        if _cov is not None:
            aveXS = xs.integrateTwoFunctionsWithUncertainty(s,
                                                            domainMin=_Emin_asPQU, domainMax=_Emax_asPQU,
                                                            useCovariance=_cov is not None, covariance=_cov,
                                                            normalize=True)
            meanValue = aveXS.getValueAs('b')
            uncertainty = aveXS.getUncertaintyValueAs('b')

        else:
            # I don't understand why I have to divide by 2 to get this right?
            meanValue = (xs.integrateTwoFunctions(s, domainMin=_Emin, domainMax=_Emax)/norm).value/2.0
            uncertainty = 0.0

        # it worked, save result
        result.append([E0, meanValue, uncertainty])

    # Compute the analytic shape of the inelastic threshold
    if _EThreshold > _Emin:
        iCut = None
        Ecut = None
        xsAtEcut = None
        dxsAtEcut = None
        for iCut, triplet in enumerate(result):
            print(iCut, triplet)
            if triplet[0] >= _EThreshold + 2.5 * _dE and triplet[1] > 0.0:
                Ecut = triplet[0]
                xsAtEcut = triplet[1]
                dxsAtEcut = triplet[2]
                break
        if _dataType == 'crossSection':
            threshXS = InelasticThreshold(Ecut, _EThreshold, xsAtEcut).domainSlice(domainMax=Ecut)
        else:
            threshXS = [[_EThreshold, 0.0, 0.0]]
        return [[xy[0], xy[1], xy[1] * dxsAtEcut / xsAtEcut] for xy in threshXS] + result[iCut:]

    return result


# -------------------------------------------
# Reading routines --------------------------
# -------------------------------------------
def readXYdYData(filename):
    if not os.path.exists(filename):
        raise IOError("XYdY data file %s not found" % filename)
    results = []
    for line in open(filename).readlines():
        try:
            results.append(map(float, line.split()[0:3]))
        except ValueError:
            pass
    return results


def readXYdXdYData(filename):
    if not os.path.exists(filename):
        raise IOError("XYdXdY data file %s not found" % filename)
    results = []
    for line in open(filename).readlines():
        try:
            results.append(map(float, line.split()[0:4]))
        except ValueError:
            pass
    return results


# -------------------------------------------
# Conversion routine ------------------------
# -------------------------------------------
def convertXYdYToGNDCrossSectionAndCovariance(xydyList, doCov=True):
    # Load the cross section into a pointwise cross section object
    import fudge.reactionData.crossSection
    _theXS = fudge.reactionData.crossSection.Component()
    _theXS.add(
        pointwise(
            axes=defaultAxes(xName='energy', xUnit='eV', yName='crossSection', yUnit='b'),
            data=[x[0:2] for x in xydyList]))

    # Load the uncertainty into a covariance object
    # This is tricky because xydy tables usually mean points with associated error bars
    # while ENDF covariances are grouped values.
    if doCov:
        import fudge.covariances.covarianceMatrix

        # set up covariance group boundaries
        dE = (xydyList[1][0] - xydyList[0][0]) / 2.0
        energyBounds = [max(xydyList[0][0] - dE, 0.0)]
        for i, triplet in enumerate(xydyContents[:-1]):
            energyBounds.append(triplet[0] + dE)
            dE = xydyContents[i + 1][0] - energyBounds[-1]
        energyBounds.append(xydyList[-1][0] + dE)

        # axes are standard covariance ones
        axes = axesModule.Axes(labelsUnits={0: ('matrix_elements', 'b**2'),
                                            1: ('column_energy_bounds', 'eV'),
                                            2: ('row_energy_bounds', 'eV')})
        axes[2] = axesModule.Grid(axes[2].label, axes[2].index, axes[2].unit,
                                  style=xDataEnumsModule.GridStyle.boundaries,
                                  values=valuesModule.Values(energyBounds))
        axes[1] = axesModule.Grid(axes[1].label, axes[1].index, axes[1].unit,
                                  style=xDataEnumsModule.GridStyle.boundaries,
                                  values=linkModule.Link(link=axes[2].values, relative=False))

        # matrix is diagonal, constructed from uncertainties
        matrix = arrayModule.Diagonal([triplet[2] ** 2 for triplet in xydyList])

        # the covariance itself
        _cov = fudge.covariances.covarianceMatrix.CovarianceMatrix(
            label='converted', matrix=griddedModule.Gridded2d(axes, matrix))

    else:
        _cov = None

    return _theXS, _cov


# -------------------------------------------
# Command line controls ---------------------
# -------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(description='Smooth the cross section in an evaluation')

    # Required command line options
    parser.add_argument('MT', type=int, help='MT of the cross section to smooth.')
    parser.add_argument('ENDF', type=str, help='ENDF file(s) whose cross section you want to smooth.')

    # Set output
    parser.add_argument('-o', dest='outFile', default=None, type=str,
                        help='Output to a file called OUTFILE, instead of printing to stdout.')

    # Verbosity
    parser.add_argument('-v', default=False, action='store_true', help="Enable verbose output.")

    # Tell the code about your file, if it can't tell by itself
    parser.add_argument('--notAnENDFFile', default=False, action="store_true",
                        help="Input ENDF file is not really an ENDF file, but rather a table of x,y,dy values to "
                             "be interpreted as the cross section corresponding to this MT.")
    parser.add_argument('--EThreshold', type=float, default=None,
                        help="Threshold energy (in MeV) of this MT, if we can't determine it from the ENDF file.")
    parser.add_argument('--noCov', dest='useCov', default=True, action='store_false',
                        help="Don't use the ENDF covariance data so don't do uncertainty calculation")

    # Special controls for angular distributions
    parser.add_argument('--angDist', default=False, action='store_true',
                        help="Work with the outgoing angular distribution (assumed for neutrons only, in MF=4)")
    parser.add_argument('--L', type=int, default=None, help="L of the angular moment to work with")

    # Misc. controls
    parser.add_argument('--printThresholds', default=False, action="store_true",
                        help="Print out the thresholds & Q values for all reactions in ENDF (if it is an ENDF "
                             "file that is)")
    parser.add_argument('--plot', default=False, action='store_true',
                        help="Plot cross section and its smoothed version (if applicable)")
    parser.add_argument('--plotProfile', default=False, action='store_true',
                        help='Plot the profile function used in smoothing (if not "group")')

    # Smoothing controls
    parser.add_argument('--smoothing', default=None, choices=['group', 'lorentzian', 'gaussian', 'window'],
                        help='Profile function to use while smoothing.  Specify the width of the smoothing function '
                             'with the "--dE" option.  Note if the "group" option is used, this code will group, '
                             'not smooth, using constant group boundaries given by the "--dE" option.  (default None')
    parser.add_argument('--dE', type=float, default=1e5,
                        help="Give width (dE) of smoothing profile in eV (default 1e5 eV)")
    return parser.parse_args()


# ---------------------------------------------------------
#                   Main Routine!
#  --------------------------------------------------------
if __name__ == "__main__":

    # -------------------------------------------
    # Parse the command line --------------------
    # -------------------------------------------
    args = parse_args()

    # -------------------------------------------
    # Clear everything before we get started ----
    # -------------------------------------------
    theRxnSuite, theCovSuite, theXS, theXSPtwise, EThreshold, Emax, dataType = None, None, None, None, None, None, None

    # -------------------------------------------
    # Read in the stuff to average --------------
    # -------------------------------------------
    if not os.path.exists(args.ENDF):
        raise IOError("File named " + args.ENDF + " doesn't exist!")

    if args.notAnENDFFile:

        # Bail if we're not doing cross sections
        if args.angDist:
            raise NotImplementedError("Don't know how to read angular distributions from a non-evaluation file")

        # Read in the file
        xydyContents = readXYdYData(args.ENDF)

        # Convert the contents into GND objects
        theXS, cov = convertXYdYToGNDCrossSectionAndCovariance(xydyContents)

        # Set EThreshold
        if args.EThreshold is not None:
            EThreshold = args.EThreshold * 1e6  # to get in eV

    else:

        # Read in the evaluation
        if open(args.ENDF).readline().startswith("<?xml"):
            from fudge import reactionSuite
            from fudge.covariances import covarianceSuite

            theRxnSuite = reactionSuite.ReactionSuite.readXML_file(args.ENDF)
            if os.path.exists(args.ENDF.replace('gnds.xml', 'gndsCov.xml')):
                theCovSuite = covarianceSuite.CovarianceSuite.readXML_file(args.ENDF.replace('gnds.xml', 'gndsCov.xml'))
            else:
                theCovSuite = None
        else:
            from brownies.legacy.converting import endfFileToGNDS
            rce = endfFileToGNDS.endfFileToGNDS(args.ENDF, toStdOut=False, skipBadData=True)
            theRxnSuite, theCovSuite = rce['reactionSuite'], rce['covarianceSuite']

        # Get the correct reaction from the reaction suite
        rxn = theRxnSuite.getReaction(args.MT)

        # Set EThreshold & dataType
        if hasattr(rxn, "getThreshold"):
            EThreshold = rxn.getThreshold(unit='eV')
        else:
            EThreshold = 1e-5
        dataType = 'crossSection'  # by default
        Emax = None  # let the data set the upper limit

        # Get the data if we're doing angular distributions
        if args.angDist:
            raise NotImplementedError("Angular distributions are tricky, please rewrite this routine")

        # Get the data for cross sections    
        else:

            # Get the correct covariance from the covariance suite
            cov = None
            if args.useCov:
                if hasattr(theCovSuite, 'sections'):
                    covList = [c for c in theCovSuite.sections if
                               hasattr(c, 'rowData') and c.rowData.attributes['ENDF_MFMT'] == '33,%i' % args.MT]
                    if not covList:
                        print("# MT = %i covariance not present in the file" % args.MT)
                    else:
                        try:
                            cov = covList[0].evaluated.toCovarianceMatrix()
                        except ValueError:
                            pass
                else:
                    print("# No sections in covariance file?")
                if cov is None:
                    print("# No covariance data used")
                else:
                    print("# Using covariance data")

            # Get the pointwise version of the cross section
            theXS = rxn.crossSection
            theXSPtwise = theXS.toPointwise_withLinearXYs()

    # -------------------------------------------
    # Output Threshold Information --------------
    # -------------------------------------------
    if args.printThresholds:
        if not args.notAnENDFFile:
            print("   MT      Q (MeV)   Ethreshold (MeV)    Reaction")
            for r in theRxnSuite.reactions:
                print(str(r.attributes['ENDF_MT']).rjust(5), str(r.getQ('MeV')).rjust(12),
                      str(r.getThreshold('MeV')).rjust(18), '  ', r.toString().strip())
        else:
            raise ValueError("Cannot extract threshold information without a valid ENDF file as input")

    # -------------------------------------------
    # Do the Smoothing---------------------------
    # -------------------------------------------
    outputList = []
    if args.smoothing in ['lorentzian', 'window', 'gaussian']:
        outputList = smooth(theXSPtwise, cov, profile=args.smoothing, _dE=args.dE, _EThreshold=EThreshold, _Emax=Emax,
                            _dataType=dataType)
    elif args.smoothing == 'group':
        raise NotImplementedError("Grouping: Please write me!")
    elif args.smoothing is None:
        if cov is None:
            for ESig in theXSPtwise:
                outputList.append([ESig[0], ESig[1]])
        else:
            unc = cov.getUncertaintyVector(theXSPtwise, relative=False)
            for ESig in theXSPtwise:
                dSig = unc.getValue(ESig[0])
                outputList.append([ESig[0], ESig[1], dSig])
    else:
        raise ValueError("No valid smoothing approach given")

    # -------------------------------------------
    # Output ------------------------------------
    # -------------------------------------------
    if args.outFile is None:
        print('\n'.join(['    '.join(map(str, x)) for x in outputList]))
    else:
        with open(args.outFile, mode='w') as outFile:
            outFile.write('\n'.join(['    '.join(map(str, x)) for x in outputList]))

    # -------------------------------------------
    # Plotting ----------------------------------
    # -------------------------------------------
    if args.plot:
        theXS.plot()

    if args.plotProfile and args.smoothing in ['lorentzian', 'window', 'gaussian']:
        profileOptions[args.smoothing](1.0e6, args.dE).plot()
