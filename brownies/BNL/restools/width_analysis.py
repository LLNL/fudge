import math
import numpy


def Chi2Distribution(y, nu, norm, verbose=False):
    """
    Porter-Thomas distribution
    Froehner's eq. (277).
    Really just a chi^2 distribution with nu degrees of freedom

    :param y:
    :param nu:
    :param norm:
    :param verbose:
    :return:
    """
    if verbose:
        print('        y', y)
        print('        norm', norm)
        print('        nu', nu)
    gam = math.gamma(nu / 2.0)
    if nu < 2.0:
        ycut = pow(1e4 * gam, -1.0 / (1.0 - nu / 2.0))
        y[y < ycut] = ycut
    return norm * numpy.exp(-y) * pow(y, nu / 2.0 - 1.0) / gam


def makeWidthHistogram(widthList, verbose=False):
    aveWidth = numpy.mean(widthList)
    if verbose:
        print('    scaled widths', [w / aveWidth for w in widthList])
    hist, bin_edges = numpy.histogram([w / aveWidth for w in widthList])
    if verbose:
        print('    hist', hist)
        print('    bins', bin_edges)
    return hist, bin_edges


def getPorterThomasFitToWidths(reducedWidthList, verbose=False):
    """
    Perform a channel-by-channel fit of the histogram of widths to a Porter-Thomas distribution.

    :param reducedWidthList:
    :param verbose:
    :return:
    """
    try:
        import scipy.optimize
    except ImportError:
        print("WARNING: scipy.optimize not imported, Porter-Thomas analysis not done")
        return {}

    if verbose:
        print('    widths', reducedWidthList)
    if all([w == reducedWidthList[0] for w in reducedWidthList]):
        if verbose:
            print('    widths identical')
        results = {'dof': 1.0, 'ddof': 0.0}
    else:
        norm = len(reducedWidthList)
        hist, bin_edges = makeWidthHistogram(reducedWidthList, verbose=verbose)

        if False:
            from xData import XYs1d
            xys = [[], []]
            for i in range(len(hist)):
                xys[0].append((bin_edges[i] + bin_edges[i + 1]) / 2)
                xys[1].append(hist[i])
            XYs1d.XYs1d(data=xys, dataForm='xsandys').plot()

        if verbose:
            print('    fitting:')
        popt, pcov = scipy.optimize.curve_fit(
            lambda y, nu: Chi2Distribution(y=y, nu=nu, norm=norm, verbose=verbose),
            [(bin_edges[i + 1] + bin_edges[i + 2]) / 2 for i in range(len(bin_edges) - 2)],
            # dump 1st bin, missing resonances there
            hist[1:],
            bounds=[0.0, numpy.inf])

        if verbose:
            print('    popt', popt)
            print('    pcov', pcov)
        results = {'dof': popt[0], 'ddof': numpy.sqrt(pcov[0, 0])}
    return results


def getFractionMissingPorterThomasFitToWidths():
    """
    Estimate the fraction of missing widths (resonances) from the
    quality of the fit of the width distribution
    """
    raise NotImplementedError("write me")
