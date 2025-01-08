#! /usr/bin/env python

import argparse
from pqu import PQU
from brownies.legacy.converting import endfFileToGNDS
from fudge.core.utilities.brb import winged_banner
from fudge.productData.distributions import angular
from xData import XYs1d as XYs1dModule
from xData import axes as axesModule
import numpy.polynomial.legendre as legendreModule

ENABLEANGDISTS = True
TOOMANYENERGIES = 5000
FIRSTLINE = " $Rev:: 700      $  $Date:: 2016-01-13#$                             1 0  0    0"


def get_sad(theRxnSuite, _styleLabel, _MT, outgoingParticle='n'):
    # Get the correct reaction from the reaction suite
    rxn = None
    for r in theRxnSuite.reactions:
        if r.ENDF_MT == _MT:
            rxn = r
            break
    if rxn is None:
        raise KeyError("MT %i not present in the file" % _MT)

    return rxn \
        .outputChannel \
        .getProductWithName(outgoingParticle) \
        .distribution[_styleLabel].subforms[0]


def Legendre_fit_to_xy_data(_xys, Lmax):
    if _xys.isIsotropic():
        return [1.0]
    xys = (_xys - 0.5).copyDataToXsAndYs()  # subtract L=0 term, that better give 1/2
    Lmax = min(max(1, len(xys[1]) - 3), Lmax)
    minResidual = 100.0e12
    legs = []
    # iterate re-fitting until we converge, each step adding another L
    # if we overfit (with too many L's), the residual will get worse
    for l in range(Lmax+1):
        fit_result = legendreModule.Legendre.fit(xys[0], xys[1], l, full=True)
        legs = fit_result[0].coef
        residuals = min(fit_result[1][0])
        if residuals > 0.0:
            if minResidual < residuals:
                break
            minResidual = min(minResidual, residuals)
    legs[0] = 0.5  # set the L=0 term since it had better be 1/2 (before fixing normalization)
    return [coeff * 2. / (2. * L + 1.) for L, coeff in enumerate(legs)]  # numpy & ENDF normalization a little different


def pack_datatable(_sad, Lmax=None, Emax=None):
    from brownies.BNL.inter.datatables import DataTable
    from xData import table as tableModule

    if Lmax is None:
        Lmax = 100
    if Emax is None:
        Emax = 1e9

    def __do_a_region(__sadRegion):
        sadRegion = get_LegendreMoments(__sadRegion, verbose=False)

        # Get all the L's
        __Llist = [_L for _L in sadRegion.keys() if _L <= Lmax]
        __Llist.sort()

        # Get all the E's
        __egrid = [x[0] for x in sadRegion[__Llist[0]].copyDataToXYs() if x[0] <= Emax]

        # get data table for all L's and E's
        __data = []
        for i, E in enumerate(__egrid):
            row = [E]
            for _L in __Llist:
                row.append(sadRegion[_L].evaluate(E))
            __data.append(row)
        return __Llist, __egrid, __data

    data = []
    egrid = []
    Llist = []
    if isinstance(_sad, angular.Regions2d):
        for iregion, region in enumerate(_sad):
            _Llist, _egrid, _data = __do_a_region(region)
            egrid += _egrid
            data += _data
            for L in _Llist:
                if L not in Llist and L <= Lmax:
                    Llist.append(L)
    elif isinstance(_sad, angular.XYs2d):
        Llist, egrid, data = __do_a_region(_sad)
    cols = [tableModule.ColumnHeader(0, "E", 'eV')] + [tableModule.ColumnHeader(L+1, 'L=%i' % L, "") for L in range(0, Lmax+1)]
    return DataTable(data=data, columns=cols)


def get_LegendreMoments(theAngDist, defaultLmax=30, verbose=False):
    """

    :param theAngDist:
    :param defaultLmax:
    :param verbose:
    :return:
    """
    # Reorder the distributions so we have a table of energies for each L
    if isinstance(theAngDist, angular.Legendre):
        Lmax = len(theAngDist[-1])
    else:
        Lmax = defaultLmax
    if verbose:
        print(Lmax)

    theLegMoments = {L: [] for L in range(Lmax + 1)}
    for ePair in theAngDist:
        E = ePair.outerDomainValue

        # For Legendre moment data, we just copy it over
        if isinstance(ePair, angular.Legendre):
            for L, coeff in enumerate(ePair):
                theLegMoments[L].append((E, coeff))

        # For angle tables, do a Legendre fit to get the coefficients
        elif isinstance(ePair, angular.XYs1d):
            for L, coeff in enumerate(Legendre_fit_to_xy_data(ePair, Lmax)):
                theLegMoments[L].append((E, coeff))

        # What the heck?
        else:
            raise TypeError("must be an angular distribution, got %s" % str(type(theAngDist)))

    if verbose:
        print("L=1 moment:", theLegMoments[1])

    # Construct the axes for all the XY's we'll need below
    coeffAxes = axesModule.Axes(2)
    coeffAxes[0] = theAngDist.axes[0]
    coeffAxes[1] = theAngDist.axes[-1]

    # Make sure all moments have the same domain before turning them into XYs1d objects
    try:
        domainMin = min([min([x[0] for x in theLegMoments[L]]) for L in theLegMoments])
        domainMax = max([max([x[0] for x in theLegMoments[L]]) for L in theLegMoments])
        for L in theLegMoments:
            defaultVal = 0.0
            if L == 0:
                defaultVal = 1.0
            if theLegMoments[L][0][0] > domainMin:
                theLegMoments[L].insert(0, (domainMin, defaultVal))
            if theLegMoments[L][-1][0] < domainMax:
                theLegMoments[L].append((domainMax, defaultVal))
    except ValueError:
        pass

    # Convert them to XYs
    return {L: XYs1dModule.XYs1d(data=theLegMoments[L], axes=coeffAxes) for L in theLegMoments}


def smooth_pointwise(theAngDist, weight_by_xs=None, egrid=None, NE=TOOMANYENERGIES):
    """
    Smooth angular distribution by grouping.  Optionally weight theLegendre moments by a cross
    section to emulate grouping the full differential cross section.

    :param theAngDist: a Python dict of XYs1d containing the Legendre moments keyed by L
    :param weight_by_xs: the XYs1d of the cross section to use as a weighting function, if specified
    :param egrid: the energy grid to group with, if specified.  If not specified, use TOOMANYENERGIES equally
                  spaced bins over the domain of the angular distribution's XYs1d
    :param NE:
    :return:
    """
    Lmax = len(theAngDist[-1].coefficients) - 1

    # Set up common E grid for smoothed distributions
    emin = theAngDist.domainMin
    emax = theAngDist.domainMax
    if egrid is None:
        dE = (emax - emin) / NE
        egrid = [ie * dE + emin for ie in range(0, NE + 1)]

    theLegMomentsAsXYs = get_LegendreMoments(theAngDist)

    # Weight by relevant cross section
    if weight_by_xs is not None:
        for L in theLegMomentsAsXYs:
            theLegMomentsAsXYs[L].setSafeDivide(True)
        theLegMomentsAsXYs = {L: theLegMomentsAsXYs[L] * weight_by_xs.domainSlice(theLegMomentsAsXYs[L].domainMin,
                                                                                  theLegMomentsAsXYs[L].domainMax)
                              for L in range(Lmax + 1)}

    # Smooth them
    theLegMomentsAsXYs = {L: theLegMomentsAsXYs[L].group(egrid, norm='dx', asXYs=True) for L in range(Lmax + 1)}

    # Convert the XYs back into pointwise
    newAngDist = theAngDist.copy()
    newAngDist.functionals = []
    for iE, E in enumerate(egrid):
        # Undo xs weighting
        if weight_by_xs is not None:
            newAngDist.append(
                angular.Legendre(
                    outerDomainValue=E,
                    coefficients=[theLegMomentsAsXYs[L][iE][1] / weight_by_xs.domainSlice(theLegMomentsAsXYs[L].domainMin,
                                                                                          theLegMomentsAsXYs[L].domainMax)
                    .group(egrid, norm='dx', asXYs=True).evaluate(E) for L in range(Lmax + 1)]))
        else:
            newAngDist.append(angular.Legendre(outerDomainValue=E,
                                       coefficients=[theLegMomentsAsXYs[L][iE][1] for L in range(Lmax + 1)]))

    return newAngDist


def smooth_angular_distribution(theRxnSuite, styleLabel, MT, outgoingParticle='n', weight_by_xs=False,
                                NE=TOOMANYENERGIES):
    """

    :param theRxnSuite:
    :param styleLabel:
    :param MT:
    :param outgoingParticle:
    :param weight_by_xs:
    :param NE:
    :return:
    """
    # Get the correct reaction from the reaction suite
    rxn = None
    for r in theRxnSuite.reactions:
        if r.ENDF_MT == MT:
            rxn = r
            break
    if rxn is None:
        raise KeyError("MT %i not present in the file" % args.MT)

    # Get the cross section for weighting, if requested.  Note, this changes the meaning of the variable
    # `weight_by_xs` within this subroutine.  Usually a very bad programming practice even though it adds clarity here.
    if weight_by_xs:
        weight_by_xs = rxn.crossSection.toPointwise_withLinearXYs(lowerEps=1e-8, upperEps=1e-8)
    else:
        weight_by_xs = None

    # Compute the angular distribution 
    theAngDist = rxn \
        .outputChannel \
        .getProductWithName(outgoingParticle) \
        .distribution[styleLabel].subforms[0]
    if isinstance(theAngDist, angular.Regions2d):
        for iregion, region in enumerate(theAngDist):
            if len(region) > NE:
                print('Smooth region', iregion, 'it has', len(region), 'points')
                theAngDist[iregion] = smooth_pointwise(region, weight_by_xs=weight_by_xs, NE=NE)
    elif isinstance(theAngDist, angular.XYs2d):
        if len(theAngDist) > NE:
            print('Smooth angular distribution, it has', len(theAngDist), 'points')
            theAngDist = smooth_pointwise(theAngDist, weight_by_xs=weight_by_xs, NE=NE)
        else:
            print('No need to smooth angular distribution, the number of points ', len(theAngDist), '<= NE =', NE)
            return
    else:
        raise TypeError('distribution must be of type pointwise or piecewise to smooth, got %s' %
                        str(theAngDist.__class__))

    # Overwrite 
    rxn.outputChannel \
        .getProductWithName(outgoingParticle) \
        .distribution[styleLabel].subforms[0] = theAngDist


def parse_args():
    """Process command line options"""
    parser = argparse.ArgumentParser(description='Reconstruct the angular distributions from an ENDF file')
    parser.add_argument('-i', dest='inFile', type=str,
                        help='Input ENDF file')
    parser.add_argument('-v', dest='verbose', default=False, action='store_true',
                        help='Enable verbose output (Default: False)')
    parser.add_argument('--makeEvalStyle', default=False, action='store_true',
                        help="Make the reconstructed angular distribution the 'eval' style")
    parser.add_argument('--reconstructAngDist', default=False, action='store_true',
                        help="Reconstruct the angular distribution from RRR parameters")
    parser.add_argument('--smooth', default=False, action='store_true', help="Smooth the angular distributions")
    parser.add_argument('--MTs', default=[2], nargs='+', type=int,
                        help="List of MT's whose angular distributions you want to smooth (default: [2])")
    parser.add_argument('--weight_by_xs', default=False, action='store_true',
                        help="Weight the smoothed angular distributions by the appropriate cross section")
    parser.add_argument('--Lmax', default=None, type=int,
                        help="For txt formatted output, the maximum Legendre moment to output")
    parser.add_argument('--Emax', default=None, type=float,
                        help="For txt formatted output, the maximum incident energy to output in eV")
    parser.add_argument('--NE', default=TOOMANYENERGIES, type=int,
                        help="When smoothing, use this number of energy bins (Default: %i)" % TOOMANYENERGIES)
    parser.add_argument('-s', dest='skipBadData', default=False, action='store_true',
                        help='Do not halt if Fudge encounters bad data (Default: False)')
    parser.add_argument('-o', dest='outPrefix', default='junk', type=str,
                        help='Output file prefix for files (Default: "junk").')
    parser.add_argument('-f', dest='outFormat', default='gnds', choices=['gnds', 'endf', 'ace', 'txt'],
                        help='Output file format (Default: "gnds").  Note .gnds.xml and .gndsCov.xml extensions '
                             'are used for GNDS files, .ace for ACE files and .endf for ENDF files.')
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    # Read file
    print(winged_banner("Reading %s" % args.inFile))
    endfEval = endfFileToGNDS.endfFileToGNDS(args.inFile, toStdOut=args.verbose, skipBadData=args.skipBadData)

    # Process angular distributions
    if args.reconstructAngDist:
        if args.makeEvalStyle:
            styleLabel = 'eval'
        else:
            styleLabel = 'angularDistributionReconstructed'
        print(winged_banner("Reconstructing angular distributions for style %s" % styleLabel))
        endfEval['reactionSuite'].reconstructResonancesAngularDistributions(styleLabel,
                                                                            overwrite=(styleLabel == 'eval'))
    else:
        styleLabel = 'eval'

    if args.smooth:
        for MT in args.MTs:
            print(winged_banner("Smoothing MT=%i angular distribution" % MT))
            smooth_angular_distribution(endfEval['reactionSuite'], styleLabel, MT, weight_by_xs=args.weight_by_xs,
                                        NE=args.NE)

    # Save files
    if args.outFormat == 'gnds':
        print(winged_banner("Saving GNDS file to %s" % args.outPrefix + '.gnds.xml'))
        endfEval['reactionSuite'].saveToFile(args.outPrefix + '.gnds.xml')
        if endfEval['covarianceSuite'] is not None:
            print(winged_banner("Saving GNDS covariance file to %s" % args.outPrefix + '.gndsCov.xml'))
            endfEval['covarianceSuite'].saveToFile(args.outPrefix + '.gndsCov.xml')

    elif args.outFormat == 'ace':
        import brownies.LANL.toACE.gndsToACE as toACEModule

        print(winged_banner("Saving ACE file to %s" % args.outPrefix + '.ace'))
        toACEModule.toACE(endfEval['reactionSuite'], fileName=args.outPrefix + '.ace', evaluationId="",
                          temperature=PQU.PQU(0.0, 'MeV'), productData={}, addAnnotation={})

    elif args.outFormat == 'endf':

        print(winged_banner("Saving ENDF file to %s" % args.outPrefix + '.endf'))
        outLines = endfEval['reactionSuite'].toENDF6('eval', {'verbosity': args.verbose * 10},
                                                     covarianceSuite=endfEval['covarianceSuite']).split('\n')
        endfTxt = '\n'.join([FIRSTLINE, outLines[1].replace('-1', ' 8')] + outLines[2:])
        open(args.outPrefix + '.endf', mode='w').write(endfTxt)

    elif args.outFormat == 'txt':
        outStream = open(args.outPrefix + '.txt', mode='w')
        for MT in args.MTs:
            outStream.write(('# MT=%i\n' % MT) + '\n'.join(pack_datatable(get_sad(endfEval['reactionSuite'],
                                                                          styleLabel, MT), Lmax=args.Lmax,
                                                                          Emax=args.Emax).toStringList(na_rep='0.0')))

    else:
        raise ValueError("Unknown output format: %s" % args.outFormat)
