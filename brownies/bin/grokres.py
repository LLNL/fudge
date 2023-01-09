#!/usr/bin/env python

import argparse
import collections
import xData
import numpy

from fudge.core.utilities import banner
from fudge.processing.resonances.reconstructResonances import spins_equal
from brownies.BNL.inter import *
from brownies.BNL.restools.level_analysis import *
#from brownies.BNL.utilities.html import TABLE, IMG
from brownies.BNL.utilities.fitting import linear_regression
from brownies.BNL.utilities.io import read_evaluation

# Default plot energy grid
dE = (1e6 - 1e-5) / 100
defaultEgrid = [dE * i for i in range(101)]
defaultLs = (0, 1, 2, 3, 4, 5, 6)


# ---------------------------------------------
# Famous distributions, good for comparisons
# ---------------------------------------------
def PorterThomas_distribution(y, nu):
    """
    Porter-Thomas distribution
    Froehner's eq. (277).
    Really just a chi^2 distribution with nu degrees of freedom

    :param y:
    :param nu:
    :return:
    """
    gam = math.gamma(nu / 2.0)
    if nu < 2.0:
        ycut = pow(1e4 * gam, -1.0 / (1.0 - nu / 2.0))
        y[y < ycut] = ycut
    return numpy.exp(-y) * pow(y, nu / 2.0 - 1.0) / gam


def Wigner_surmise_distribution(x):
    constant = (math.pi / 2)
    exponent = ((-1) * (math.pi / 4) * x * x)
    return (constant * x) * math.exp(exponent)


def Poisson(x, _k=1):
    return numpy.power(x, _k) * numpy.exp(-x) / math.factorial(_k)


# ---------------------------------------------
# Plot helpers
# ---------------------------------------------
# PlotDataSet = namedtuple( "PlotDataSet", ['data', 'type', 'legend', 'color', 'kwrds']  )

class PlotDataSet:
    def __init__(self, data, _type, legend=None, color=None, **kwrds):
        """
        :param data: holds actual data.  See under 'type' for how data should be given
        :param _type: can be 'hist', 'scatter_xy', 'scatter_xydy', 'scatter_xydxdy', 'line_xy' or 'step_xy'

                     * For _type=='hist', 'data' is a list of values that will be converted into a histogram

                     * For _type in ['scatter_xy', 'line_xy', 'step_xy'], 'data' is a list such that data=[xdata, ydata]

                     * For _type=='scatter_xydy', data=[xdata, ydata, yerr]

                     * For _type=='scatter_xydxdy', data=[xdata, ydata, xerr, yerr]

        :param legend: a string with the legend entry or None if no legend entry is to be provided
        :param color:  a matplotlib color or None if you want matplotlib to pick
        :param kwrds:  a dictionary of keywords passed into the various pyplot plot functions specific to the current
                       dataset.
        """
        self.data = data
        self.type = _type
        self.legend = legend
        self.color = color
        self._kwrds = kwrds

    @property
    def bins(self):
        return self._kwrds.get('bins', numpy.logspace(math.log10(1e-4), math.log10(10.0), self.numbins + 1))

    @property
    def log(self): return self._kwrds.get('log', False)

    @property
    def numbins(self): return self._kwrds.get('numbins', 20)


def plot_functions(datasets, title=None, xaxislabel='', yaxislabel='',
                   legend_loc='upper left', show=True, subplot=None,
                   grid=False, xscale=None, yscale=None, outfile=None,
                   clear=False, xlim=None, ylim=None, **kwrds):
    """
    :param datasets: A list of PlotDataSet instances holding the things to plot.
    :param title: Title of the plot.  Set to None to disable.
    :param xaxislabel: Label of the xaxis.
    :param yaxislabel: Label of the yaxis.
    :param legend_loc: Sets the location of the plot legend.
                       See https://matplotlib.org/api/_as_gen/matplotlib.pyplot.legend.html.
    :param show: Flag to tell matplotlib to show the current plot.  Set to False if building subplots.
    :param subplot: Argument to matplotlib's pyplot subplot() function.
                    See https://matplotlib.org/api/_as_gen/matplotlib.pyplot.subplot.html.
    :param grid: Argument to matplotlib's pyplot grid() function.  Toggles grids on the plots.
    :param xscale: Argument to matplotlib's pyplot xscale() function.
                   See
                   https://matplotlib.org/gallery/pyplots/pyplot_scales.html#sphx-glr-gallery-pyplots-pyplot-scales-py.
    :param yscale: Argument to matplotlib's pyplot yscale() function.
    :param outfile: Filename to save the current figure to, in lieu of calling show().
                    If outfile is None, then we'll use show() (this is the default behavior).
    :param clear:
    :param xlim:
    :param ylim:
    **kwrds
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for ds in datasets:
        set_kwrds = {}
        if ds.legend is not None:
            set_kwrds['label'] = ds.legend
        if ds.color is not None:
            set_kwrds['color'] = ds.color
        if ds.type == 'hist':
            plt.hist(ds.data, density=True, bins=ds.bins, **set_kwrds)
        elif ds.type == 'scatter_xy':
            plt.plot(ds.data[0], ds.data[1], marker='o', linewidth=0, **set_kwrds)
        elif ds.type == 'scatter_xydy':
            if len(ds.data[0]) != len(ds.data[1]):
                continue
            plt.errorbar(ds.data[0], ds.data[1], yerr=ds.data[2], marker='o', elinewidth=1, capsize=5, linewidth=0,
                         **set_kwrds)
        elif ds.type == 'scatter_xydxdy':
            plt.errorbar(ds.data[0], ds.data[1], xerr=ds.data[2], yerr=ds.data[3], marker='o', elinewidth=1, capsize=5,
                         linewidth=0, **set_kwrds)
        elif ds.type == 'line_xy':
            plt.plot(ds.data[0], ds.data[1], **set_kwrds)
        elif ds.type == 'step_xy':
            plt.step(ds.data[0], ds.data[1], **set_kwrds)
        else:
            raise ValueError("Unknown plot type '%s'" % ds.type)
    plt.xlabel(xaxislabel, fontsize='large')
    plt.ylabel(yaxislabel, fontsize='large')
    if title is not None:
        plt.title(title)
    plt.legend(loc=legend_loc, shadow=True, fontsize='large')
    plt.grid(grid)
    if xscale is not None:
        plt.xscale(xscale)
    if yscale is not None:
        plt.yscale(yscale)
    if subplot is not None:
        plt.subplot(subplot)
    if xlim is not None:
        plt.xlim(*xlim)
    if ylim is not None:
        plt.ylim(*ylim)
    if show:
        if outfile is None:
            plt.show()
        else:
            plt.savefig(outfile)
    if clear:
        plt.clf()


# ---------------------------------------------
# Actual things to plot
# ---------------------------------------------
def plot_penetrability(_rrr, _Ls=defaultLs, channel=None, _Egrid=defaultEgrid, title_string=None, legend_loc='best',
                       **kwrds):
    plot_functions(
        [
            PlotDataSet(
                data=[_Egrid, [_rrr.penetrationFactor(L, _rrr.rho(E, channel)) for E in _Egrid]],
                legend="L=%i" % L, _type='line_xy', color=None, kwrds={}) for L in _Ls],
        title='$P_{L}(E)$ for ' + title_string,
        xaxislabel="E (eV)",
        yaxislabel='penetrability (no-dim)',
        legend_loc=legend_loc,
        **kwrds)


def plot_phase(_rrr, _Ls=defaultLs, channel=None, _Egrid=defaultEgrid, title_string=None, legend_loc='best', **kwrds):
    plot_functions(
        [
            PlotDataSet(
                data=[_Egrid, [_rrr.phi(L, _rrr.rho(E, channel)) for E in _Egrid]],
                legend="L=%i" % L, _type='line_xy', color=None, kwrds={}) for L in _Ls],
        title='$\phi_{L}(E)$ for ' + title_string,
        xaxislabel="E (eV)",
        yaxislabel='phase (rad)',
        legend_loc=legend_loc,
        **kwrds)


def plot_shift(_rrr, _Ls=defaultLs, channel=None, _Egrid=defaultEgrid, title_string=None, legend_loc='best', **kwrds):
    plot_functions(
        [
            PlotDataSet(
                data=[_Egrid, [_rrr.shiftFactor(L, _rrr.rho(E, channel)) for E in _Egrid]],
                legend="L=%i" % L, _type='line_xy', color=None, kwrds={}) for L in _Ls],
        title='$S_{L}(E)$ for ' + title_string,
        xaxislabel="E (eV)",
        yaxislabel='shift (no-dim)',
        legend_loc=legend_loc,
        **kwrds)


def plot_widthddist(_rrr, aveWidth, expectedDOF, L, J, reaction, scaled=False, title_string=None, legend_loc='best',
                    **kwrds):
    # Get the widths
    widths = get_resonance_widths(_rrr, c)

    # Do Porter-Thomas fits
    theFit = _rrr.getPorterThomasFitToWidths()
    bestFitDOF = theFit[c]['dof']

    # Construct title
    title = "%s width distribution for L=%i" % (title_string, L)
    if J is not None:
        title += ', J=' + str(J)
    widthSubscript = {'capture': '\\gamma', 'elastic': 'n', 'fission': 'f', 'total': "tot"}

    # plot limits
    ylim = {'fission': (1e-3, 100.), 'elastic': (1e-3, 100.), 'capture': (1e-5, 1.0)}

    # Now do the plots
    if not scaled:
        yaxislabel = 'Frequency (1/eV)'
        xaxislabel = '$\Gamma_{%s}$ (eV)' % widthSubscript[reaction]
        norm = len(widths)
        plot_functions(
            [PlotDataSet(data=widths, legend='ENDF', _type='hist')],
            legend_loc=legend_loc,
            title=title,
            xaxislabel=xaxislabel,
            yaxislabel=yaxislabel,
            xscale='log',
            yscale='log',
            **kwrds)
    else:
        # dx=.1
        xs = numpy.logspace(math.log10(1e-4), math.log10(10.0), 201)  # numpy.array( [i*dx for i in range(1,101)] )
        yaxislabel = 'Frequency (no-dim)'
        xaxislabel = '$\Gamma_{%s}/\overline{\Gamma_{%s}}$' % (
            widthSubscript.get(reaction, reaction.split(' + ')[0]),
            widthSubscript.get(reaction, reaction.split(' + ')[0]))
        plot_functions(
            [
                PlotDataSet(
                    data=[xs, PorterThomas_distribution(xs, expectedDOF)],
                    legend='Porter-Thomas ($\\nu=% 6.2f$)' % expectedDOF,
                    _type='line_xy'),
                PlotDataSet(
                    data=[xs, PorterThomas_distribution(xs, bestFitDOF)],
                    legend='Fitted Porter-Thomas ($\\nu=% 6.2f$)' % bestFitDOF,
                    _type='line_xy'),
                PlotDataSet(
                    data=[G / aveWidth for G in widths if aveWidth > 0.0],
                    legend='ENDF', _type='hist')],
            legend_loc=legend_loc,
            title=title,
            xaxislabel=xaxislabel,
            yaxislabel=yaxislabel,
            xscale='log',
            yscale='log',
            ylim=ylim.get(reaction, (1e-3, 100.)),
            **kwrds)


def plot_nnsd(_lsa, L, J, title_string=None, legend_loc='best', **kwrds):
    numpts = 101
    xmax = 6.0
    dx = xmax / numpts
    xs = numpy.array([i * dx for i in range(1, numpts + 1)])
    title = "%s nearest neighbor spacing distribution for L=%i" % (title_string, L)
    if J is not None:
        title += ', J=' + str(J)
    plot_functions(
        [
            PlotDataSet(data=[xs, [Wigner_surmise_distribution(x) for x in xs]],
                        legend="Wigner's Surmise", _type='line_xy'),
            PlotDataSet(data=[xs, [Poisson(x, 0) for x in xs]],
                        legend="Poisson (k=0)", _type='line_xy'),
            PlotDataSet(data=_lsa.normalized_spacings, _type='hist', numbins=30, legend='ENDF')],
        title=title,
        xaxislabel="$D/\\overline{D}$",
        yaxislabel='Frequency',
        xlim=(0.0, xmax),
        legend_loc=legend_loc, **kwrds)


def plot_cumlev(_lsa, title_string=None, include_fits=True, legend_loc='best', yscale='linear', **kwrds):
    cld = _lsa.getCumulativeLevelDistribution()
    xs = numpy.array([p[0] for p in cld])
    ys = numpy.array([p[1] for p in cld])
    ds = [PlotDataSet(data=[xs, ys], legend="ENDF", _type='step_xy')]
    if include_fits:
        half = len(xs) // 2

        def flin(E, a, b): return a + b * E

        reslin = linear_regression(xs[:half], ys[:half])
        ds.append(
            PlotDataSet(data=[xs, flin(xs, reslin[0], reslin[1])],
                        legend="Fit to $c+E/\overline{D}$,\n$\overline{D}=%6.2f\pm%6.2f$\nFit range (0.0, %6.2f) eV" % (
                        1.0 / reslin[1], reslin[2] / reslin[1] / reslin[1], xs[half]),
                        _type='line_xy'))
        # Residuals (not plotted yet)
        res_ds = [
            PlotDataSet(data=[xs, 0.0 * ys], _type='line_xy'),
            PlotDataSet(data=[xs, ys - flin(xs, reslin[0], reslin[1])], _type='step_xy')]
    plot_functions(ds,
                   title='Cumulative level distribution for ' + title_string,
                   xaxislabel="E (eV)",
                   yaxislabel='Num. levels with $E_i < E$',
                   yscale=yscale,
                   legend_loc=legend_loc,
                   **kwrds)


def plot_Delta3(_lsa, L, J, title_string=None, legend_loc='best', **kwrds):
    Ls = list(range(_lsa.Lmin, _lsa.Lmax + 1))
    d3 = _lsa.getDysonMehtaDelta3_vs_L()
    #    print Ls, len(Ls)
    #    print d3[0], len(d3[0])
    #    print d3[1], len(d3[1])
    #   need to check lenghts of these guys, sometimes don't have enough L's see 90Zr (original)
    ds = []
    ds.append(PlotDataSet(
        data=[Ls, d3[0][_lsa.Lmin:], d3[1][_lsa.Lmin:]],
        legend="ENDF", _type='scatter_xydy'))
    ds.append(PlotDataSet(
        data=[Ls, [getDysonMehtaDelta3Systematics_for_GOE(LL) for LL in Ls]],
        legend="GOE systematics", _type='line_xy'))
    ds.append(PlotDataSet(
        data=[Ls, [getDysonMehtaDelta3Systematics_for_Poisson(LL) for LL in Ls]],
        legend="Poisson systematics", _type='line_xy'))
    ds.append(PlotDataSet(
        data=[Ls, [getDysonMehtaDelta3Systematics_for_PicketFence(LL) for LL in Ls]],
        legend="Picket fence (SHO) systematics", _type='line_xy'))
    plot_functions(ds,
                   title='%s  spectral stiffness for L=%i, J=%s' % (title_string, L, J),
                   xaxislabel='L', yaxislabel='$\\Delta_3(L)$',
                   legend_loc=legend_loc, **kwrds)


def plot_levlevcorr(_lsa, title_string=None, legend_loc='best', **kwrds):
    Dbar = _lsa.mean_spacing
    DD = [(D - Dbar) / Dbar for D in _lsa.spacings]
    plot_functions(
        [PlotDataSet(data=[DD[0:-2], DD[1:-1]], _type='scatter_xy')],
        xaxislabel="$(D_i-\\overline{D})/\\overline{D}$",
        yaxislabel="$(D_{i+1}-\\overline{D})/\\overline{D}$",
        legend_loc=legend_loc)


def plot_Delta3_by_E(_lsas, num_chunks, chunk_size, _Es, _D3L, L=None, J=None, title_string=None, legend_loc='best',
                     **kwrds):
    xs = []
    xerrs = []
    ys = []
    yerrs = []
    ygoe = []
    ypoiss = []
    ysho = []
    _Es = [_EE for _EE in _Es if _EE > 0]
    for ich in range(num_chunks):
        d3 = _lsas[ich].getDysonMehtaDelta3_vs_L()
        try:
            xs.append(0.5 * (_Es[ich * chunk_size] + _Es[min((ich + 1) * chunk_size, len(_Es))]))
            xerrs.append(0.5 * (_Es[(ich + 1) * chunk_size] - _Es[min((ich) * chunk_size, len(_Es))]))
        except:
            continue
        ys.append(d3[0][_D3L])
        yerrs.append(d3[1][_D3L])
        ygoe.append(getDysonMehtaDelta3Systematics_for_GOE(_D3L))
        ypoiss.append(getDysonMehtaDelta3Systematics_for_Poisson(_D3L))
        ysho.append(getDysonMehtaDelta3Systematics_for_PicketFence(_D3L))
    if len(xs) == 2:
        x2s = [xs[0] - xerrs[0], xs[-1] + xerrs[-1]]
    elif len(xs) < 2:
        return
    else:
        x2s = [xs[0] - xerrs[0]] + xs[1:-1] + [xs[-1] + xerrs[-1]]
    plot_functions(
        [
            PlotDataSet(
                data=[xs, ys, xerrs, yerrs],
                legend="ENDF", _type='scatter_xydxdy'),
            PlotDataSet(
                data=[x2s, ygoe],
                legend="GOE systematics", _type='line_xy'),
            PlotDataSet(
                data=[x2s, ypoiss],
                legend="Poisson systematics", _type='line_xy'),
            PlotDataSet(
                data=[x2s, ysho],
                legend="Picket fence (SHO) systematics", _type='line_xy')
        ],
        title='%s  spectral stiffness for segments of length %i for L=%i, J=%s' % (title_string, _D3L, L, J),
        xaxislabel='E (ev)', yaxislabel='$\\Delta_3$',
        legend_loc=legend_loc, **kwrds)


def plot_rho(_lsas, num_chunks, chunk_size, _Es, L=None, J=None, title_string=None, legend_loc='best', **kwrds):
    xs = []
    ys = []
    fail = False
    _Es = [_EE for _EE in _Es if _EE > 0]
    for ich in range(num_chunks):
        imin = ich * chunk_size
        imax = min((ich + 1) * chunk_size, len(_Es) - 1)
        try:
            rho = _lsas[ich].get2SpacingCorrelation()
            xs.append(_Es[imin])
            xs.append(_Es[imax])
            ys.append(rho)
            ys.append(rho)
        except ValueError:
            pass
    if len(xs) < 2 or len(xs) != len(ys):
        return
    plot_functions(
        [
            PlotDataSet(data=[xs, ys], _type="line_xy", legend="ENDF"),
            PlotDataSet(data=[xs, [-0.27 for x in xs]], _type="line_xy", legend="GOE prediction"),
            PlotDataSet(data=[xs, [0.0 for x in xs]], _type="line_xy", legend="Uncorrelated"),
            PlotDataSet(data=[xs, [1.0 for x in xs]], _type="line_xy", legend="Picket fence (100% correlated)")],
        title='%s  level-level spacing correlation for L=%i, J=%s' % (title_string, L, J),
        xaxislabel='E (eV)', yaxislabel='$\\rho$',
        legend_loc=legend_loc, **kwrds)


def plot_fraction_missing(_lsas, num_chunks, chunk_size, _Es, L=None, J=None, skipDelta3=False, title_string=None,
                          legend_loc='best', **kwrds):
    xs = []
    xerrs = []
    rhos = []
    rhoerrs = []
    d3s = []
    d3errs = []
    for ich in range(num_chunks):
        imin = ich * chunk_size
        imax = min((ich + 1) * chunk_size, len(Es) - 1)
        xs.append(0.5 * (_Es[imin] + _Es[imax]))
        xerrs.append(0.5 * (_Es[imin] - _Es[imax]))
        try:
            rho = _lsas[ich].getFractionMissing2SpacingCorrelation()
            rhos.append(rho[0])
            rhoerrs.append(rho[1])
        except ValueError:  # need better solution for bad data
            rhos.append(0.0)
            rhoerrs.append(0.0)
        if not skipDelta3:
            delta3 = _lsas[ich].getFractionMissingFromDysonMehtaDelta3()
            d3s.append(delta3[0])
            d3errs.append(delta3[1])
    ds = []
    ds.append(PlotDataSet(
        data=[xs, rhos, numpy.abs(xerrs), numpy.abs(rhoerrs)],
        _type='scatter_xydxdy', legend='Using $\\rho$'))
    if not skipDelta3:
        ds.append(PlotDataSet(
            data=[xs, d3s, numpy.abs(xerrs), numpy.abs(d3errs)],
            _type='scatter_xydxdy', legend='Using $\\Delta_3$'))
    plot_functions(ds,
                   title='%s missing resonance fraction for L=%i, J=%s' % (title_string, L, J),
                   xaxislabel='E (eV)', yaxislabel='f',
                   legend_loc=legend_loc, **kwrds)


def plot_avespacing(_rrrAveQuant, _urr, _c, title_string=None, legend_loc='best', **kwrds):
    datasets = []
    if _rrrAveQuant is not None:
        rrrxs = []
        rrrys = []
        rrrxerrs = []
        rrryerrs = []
        for i in range(len(_rrrAveQuant[_c]['spacings'])):
            rrrys.append(valueAsFloatOrZero(_rrrAveQuant[_c]['spacings'][i], 'eV'))
            rrrxs.append(0.5 * (_rrrAveQuant[_c]['energyGrid'][i + 1] + _rrrAveQuant[_c]['energyGrid'][i]))
            rrrxerrs.append(0.5 * (_rrrAveQuant[_c]['energyGrid'][i + 1] - _rrrAveQuant[_c]['energyGrid'][i]))
            rrryerrs.append(uncAsFloatOrZero(_rrrAveQuant[_c]['spacings'][i], 'eV'))
        datasets.append(
            PlotDataSet(data=[rrrxs, rrrys, rrrxerrs, rrryerrs],
                        legend='ENDF (RRR)', _type='scatter_xydxdy'))
    if _urr is not None:
        urrxs = []
        urrys = []
        LJ = (_c.l, _c.J)
        if LJ in _urr.levelSpacings:
            for p in _urr.levelSpacings[LJ]:
                urrxs.append(p[0])
                urrys.append(p[1])
            datasets.append(
                PlotDataSet(data=[urrxs, urrys],
                            legend='ENDF (URR)', _type='line_xy'))
    plot_functions(
        datasets,
        title='%s average level spacing for L=%i, J=%s' % (title_string, _c.l, _c.J),
        xaxislabel='E (eV)', yaxislabel='$\\overline{D}(E)$ (eV)',
        legend_loc=legend_loc, xscale='log', **kwrds)


def plot_avewidth(_rrrAveQuant, _urr, _c, _reaction, title_string=None, legend_loc='best', **kwrds):
    datasets = []
    if _rrrAveQuant is not None:
        rrrxs = []
        rrrys = []
        rrrxerrs = []
        rrryerrs = []
        for i in range(len(_rrrAveQuant[_c]['widths'])):
            rrrys.append(valueAsFloatOrZero(_rrrAveQuant[_c]['widths'][i], 'eV'))
            rrrxs.append(0.5 * (_rrrAveQuant[_c]['energyGrid'][i + 1] + _rrrAveQuant[_c]['energyGrid'][i]))
            rrrxerrs.append(0.5 * (_rrrAveQuant[_c]['energyGrid'][i + 1] - _rrrAveQuant[_c]['energyGrid'][i]))
            rrryerrs.append(uncAsFloatOrZero(_rrrAveQuant[_c]['widths'][i], 'eV'))
        datasets.append(
            PlotDataSet(data=[rrrxs, rrrys, rrrxerrs, rrryerrs],
                        legend='ENDF (RRR)', _type='scatter_xydxdy'))
    if _urr is not None:
        urrxs = []
        urrys = []
        if (_c.l, _c.J) in _urr.averageWidths and _reaction in _urr.averageWidths[(_c.l, _c.J)]:
            for p in _urr.averageWidths[(_c.l, _c.J)][_reaction]:
                penfact = 1.0
                if _c.isElastic:
                    penfact = math.sqrt(p[0]) * urr.penetrationFactor(_c.l, urr.rho(p[0])) / urr.rho(
                        p[0])  # See ENDF Eq. (D.99)
                urrxs.append(p[0])
                urrys.append(p[1] * penfact)
            datasets.append(PlotDataSet(data=[urrxs, urrys],
                                        legend='ENDF (URR)', _type='line_xy'))
    plot_functions(
        datasets,
        title='%s average width for L=%i, J=%s' % (title_string, _c.l, _c.J),
        xaxislabel='E (eV)', yaxislabel='$\\overline{\\Gamma}(E)$ (eV)',
        legend_loc=legend_loc, **kwrds)


def plot_transmission_coefficient(title_string=None, legend_loc='best', **kwrds):
    raise NotImplementedError("FIXME: finish coding this!")


# ---------------------------------------------
# Helper functions to extract information from resonances
# ---------------------------------------------
def get_channels(rrr, urr, L=None, S=None, J=None, rxn=None):
    """
    Find all channels in the RRR and URR that match the reaction and requested quantum numbers
    """
    channels = []
    search_list = []
    if rrr is not None:
        for c in rrr.allChannels:
            if L is not None and not spins_equal(c.l, L):
                continue
            if S is not None and not spins_equal(c.s, S):
                continue
            if J is not None and not spins_equal(c.J, J):
                continue
            if rxn is not None:
                if rxn == 'elastic' and not c.isElastic:
                    continue
                elif rxn == 'capture' and 'photon' not in [c.particleA, c.particleB]:
                    continue
                elif rxn == 'fission' and 'fission' not in c.reaction:
                    continue
            if c not in channels:
                channels.append(c)
    if urr is not None:
        print("FIXME: in get_channels Not doing URR channels")
        pass
    return channels


def get_channel(rrr, urr, L=None, S=None, J=None, rxn=None):
    """
    Find the exact match to the channel in the RRR and URR for the reaction and requested quantum numbers
    """
    channels = get_channels(rrr, urr, L=L, S=S, J=J, rxn=rxn)
    if len(channels) == 1:
        return channels[0]
    if len(channels) == 0:
        print("No matching channels")
        return None
    raise ValueError("Multiple channels match because the (L, S, J, reaction) specification is ambiguous.  "
                     "These are the matching channels:\n\t%s" % '\n\t'.join([str(x) for x in channels]))


def get_key_from_channel(c):
    reaction = c.reaction
    if 'fission' in reaction:
        reaction = 'Fission'
    if int(c.s) == c.s:
        s = c.s / 2.0
    else:
        s = c.s
    return '%s (j=%s,l=%s,s=%s)' % (str(reaction), str(float(c.J)), str(c.l), str(float(s)))


def get_resonance_widths(_rrr, _c, keepBoundLevels=False):
    widths = []
    for iR in rrr.allChannels[_c]:
        ER = rrr._energies[iR]
        if ER < 0.0 and not keepBoundLevels: continue
        widths.append(rrr.allChannels[_c][iR])
    return widths


def get_resonance_energies(_rrr, _c, keepBoundLevels=False):
    Es = [_rrr._energies[iE] for iE in _rrr.allChannels[_c].keys() if _rrr._energies[iE] > 0.0 or keepBoundLevels]
    Es.sort()
    return Es


def get_channel_report(_theEvaluation):
    _rep = collections.OrderedDict()
    metricMenu = collections.namedtuple('metricMenu', '')  # emulates results from ArgumentParser.parse_args() function
    _rep["Channels"] = report.getChannelDataTable(_theEvaluation['reactionSuite'], metricMenu())
    return _rep


def get_resonance_report(_theEvaluation):
    _rep = collections.OrderedDict()
    _rep['Main'] = report.getEvaluationReport(_theEvaluation['reactionSuite'])
    _rep.update(report.getResonanceReport(_theEvaluation['reactionSuite']))
    _rep.update(get_channel_report(_theEvaluation))
    return _rep


def get_simple_reaction_name(rxn):
    if 'fission' in rxn:
        return 'fission'
    if ' + photon' in rxn:
        return 'capture'
    if 'n + ' in rxn and '_e' not in rxn:
        return 'elastic'
    return rxn


def valueAsFloatOrZero(x, unit):
    if x is None:
        return 0.0
    return float(x.getValueAs(unit))


def uncAsFloatOrZero(x, unit):
    if x is None:
        return 0.0
    return float(x.getUncertaintyValueAs(unit))


# ---------------------------------------------
# Enums that make our life easy
# ---------------------------------------------
legend_loc_options = [
    'right',
    'center left',
    'upper right',
    'lower right',
    'best',
    'center',
    'lower left',
    'center right',
    'upper left',
    'upper center',
    'lower center']

PLOTCHOICES = [
    'nnsd', 'cumlev', 'levlevcorr',  # measures of level statistics, part 1
    'Delta3', 'Delta3_by_E', 'rho',  # measures of level statisitics, part 2
    'widthdist', 'widthdist_scaled',  # measures of width statistics
    'avewidth', 'avespacing',  # average parameter energy dependence
    'tc',  # combine the average parameters into transmission coefficients
    'hardsphere', 'penetrability', 'phase', 'shift',  # hard-sphere things
    'fraction_missing'  # the fraction of missing levels, done several ways
]

REACTIONNAMES = ['capture', 'elastic', 'fission', 'total']


# ---------------------------------------------
# Command line parsing
# ---------------------------------------------
def parse_arguments():
    """Parse the command line"""
    parser = argparse.ArgumentParser(
        description="""
            Let's grok some resonances!

            Note: the -l, -t, -p and -f options are mutually exclusive.
            Also, an output file (with -o option) must be specified if you are using the -f option.
        """)
    parser.add_argument('ENDF', type=str, help='ENDF file(s) whose cross section you want to study.')
    parser.add_argument('-l', dest='channelList', default=False, action='store_true',
                        help="Print out the list of channels")
    parser.add_argument('-t', dest='report', default=False, action='store_true',
                        help='Return the simplified resonance report')
    parser.add_argument('-f', dest='fullReport', default=False, action='store_true',
                        help='Generate the full HTML report suitable for ADVANCE.')
    parser.add_argument('-o', dest='outFile', default=None, type=str,
                        help='Output to a file called OUTFILE, instead of printing to stdout.')
    parser.add_argument('-v', dest='verbose', default=False, action='store_true', help="Enable verbose output.")
    parser.add_argument('-p', dest="plot", default=None, choices=PLOTCHOICES, help="Plot the requested quantity")
    parser.add_argument('-r', dest='reaction', default=None, choices=REACTIONNAMES,
                        help="Investigate channel with specific name")
    parser.add_argument('-J', default=None, type=float, help="Investigate channel with total angular momentum J")
    parser.add_argument('-L', default=None, type=int, help="Investigate channel with orbital angular momentum L")
    parser.add_argument('-S', default=None, type=float, help="Investigate channel with total spin S (optional for L=0)")
    parser.add_argument('--D3L', type=int, default=20,
                        help="Value of segment length to use in Delta_3 by E calculation (Default: 20)")
    parser.add_argument('-n', dest='nbins', default=2, type=int,
                        help="Number of bins to use when computing average widths and spacings (Default: 2)")
    parser.add_argument('-u', dest='computeUncertainty', default=False, action="store_true",
                        help="Compute the uncertainty on average widths and spacings")
    parser.add_argument('-c', dest='chunkSize', type=int, default=70,
                        help="Number of resonances to include in each analysis chunk when computing rho, Delta3, Q, U, "
                             "avewidth and avespacing (Default:50)")
    parser.add_argument('--legendLoc', default="best", choices=legend_loc_options,
                        help="Override default legend location and put it HERE (Default: best)")
    parser.set_defaults(reportFormat='txt')
    return parser.parse_args()


# ---------------------------------------------
# Main routine!
# ---------------------------------------------
if __name__ == "__main__":
    # Parse command line
    theArgs = parse_arguments()

    # Read in evaluation
    theEvaluation = read_evaluation(theArgs.ENDF, reconstructResonances=True)

    # Get the requested reaction suite
    rs = theEvaluation['reactionSuite']
    import fudge.processing.resonances.reconstructResonances as RRReconstruct

    # Extract target & projectile names
    target = rs.target
    projectile = rs.projectile

    # Get URR parameters
    if hasattr(rs.resonances, 'unresolved') and rs.resonances.unresolved is not None:
        rs.resonances.unresolved.evaluated.setAncestor(rs.resonances.unresolved)  # FIXME: I shouldn't need to do this!
        urr = RRReconstruct.URRcrossSection(rs.resonances.unresolved.evaluated)
        urrEgrid, flag = urr.generateEnergyGrid()
        urr.getWidthsAndSpacings()
    else:
        urr = None
        urrEgrid = None

    # Get RRR parameter averages
    if hasattr(rs.resonances, 'resolved') and rs.resonances.resolved is not None:
        resCls = RRReconstruct.getResonanceReconstructionClass(rs.resonances.resolved.evaluated)
        rrr = resCls(rs.resonances.resolved.evaluated, enableAngDists=False, verbose=theArgs.verbose)
        rrr.setResonanceParametersByChannel(warnOnly=True)
        rrrAveQuant = rrr.getAverageQuantities(nBins=theArgs.nbins, computeUncertainty=theArgs.computeUncertainty)
    else:
        rrr = None
        rrrAveQuant = None

    # Just list the channels
    if theArgs.channelList:
        ch_rep = get_channel_report(theEvaluation)["Channels"]
        if len(ch_rep) != 0:
            ch_rep.removeColumn("Threshold E")
            ch_rep.removeColumn("Eliminated?")
            ch_rep.removeColumn("Competative?")
            ch_rep.removeColumn("Pot. scatt. only?")
            ch_rep.removeColumn("RRR <Gamma>")
            ch_rep.removeColumn("RRR <D>")
            ch_rep.removeColumn("gfact")
            ch_rep.removeColumn("Relativistic?")
            ch_rep.removeColumn("No. resonances w/ ER<0")
        print('\n'.join(ch_rep.toStringList()))

    # Make the simple report
    elif theArgs.report:
        rep = get_resonance_report(theEvaluation)

        # Setup the output filestream
        if theArgs.outFile is not None:
            if theArgs.outFile.endswith(theArgs.reportFormat):
                outFile = open(theArgs.outFile, mode='w')
            else:
                outFile = open(theArgs.outFile + '.' + theArgs.reportFormat, mode='w')
        else:
            import io

            outFile = io.StringIO()

        # Output report
        for k in rep:
            if isinstance(rep[k], report.Report):
                repOut = rep[k].text_report()
            elif isinstance(rep[k], xData.table.Table):
                repOut = '\n'.join(rep[k].toStringList())
            else:
                repOut = str(rep[k])
            outFile.write(banner(k) + '\n')
            outFile.write(repOut + '\n')
        if hasattr(outFile, 'getvalue'):
            print(outFile.getvalue())

    # Make plots
    elif theArgs.plot is not None:

        # Many of the plots need to know the EXACT channel
        c = get_channel(rrr, urr, L=theArgs.L, S=theArgs.S, J=theArgs.J, rxn=theArgs.reaction)
        if c is None:
            print("Cannot make plots")
            exit()

        # Title string construction widgets
        titleString_incomming = '%s + %s' % (projectile, target)
        titleString_reaction = '%s(%s, %s)' % (target, projectile,
                                               {'elastic': 'el', 'capture': '$\\gamma$', 'fission': 'f',
                                                "total": 'tot'}[theArgs.reaction])

        # Plot hard-sphere things
        if theArgs.plot in ['hardsphere', 'penetrability', 'phase', 'shift']:
            Egrid = defaultEgrid
            if theArgs.reaction != 'elastic':
                raise ValueError('Plot "%s" only available for elastic reactions' % theArgs.plot)
            else:
                if theArgs.L is None:
                    Ls = defaultLs
                else:
                    Ls = [theArgs.L]
            if theArgs.plot == 'penetrability':
                plot_penetrability(rrr, title_string=titleString_incomming, Ls=Ls, channel=c, Egrid=defaultEgrid,
                                   legend_loc=theArgs.legendLoc, outfile=theArgs.outFile)
            elif theArgs.plot == 'phase':
                plot_phase(rrr, title_string=titleString_incomming, Ls=Ls, channel=c, Egrid=defaultEgrid,
                           legend_loc=theArgs.legendLoc, outfile=theArgs.outFile)
            elif theArgs.plot == 'shift':
                plot_shift(rrr, title_string=titleString_incomming, Ls=Ls, channel=c, Egrid=defaultEgrid,
                           legend_loc=theArgs.legendLoc, outfile=theArgs.outFile)
            else:
                plot_penetrability(rrr, title_string=titleString_incomming, Ls=Ls, channel=c, Egrid=defaultEgrid,
                                   legend_loc=theArgs.legendLoc, show=False, subplot=131)
                plot_phase(rrr, title_string=titleString_incomming, Ls=Ls, channel=c, Egrid=defaultEgrid,
                           legend_loc=theArgs.legendLoc, show=False, subplot=132)
                plot_shift(rrr, title_string=titleString_incomming, Ls=Ls, channel=c, Egrid=defaultEgrid,
                           legend_loc=theArgs.legendLoc, show=True, subplot=133, outfile=theArgs.outFile)
            exit()

        # Width distributions
        if theArgs.plot in ['widthdist', 'widthdist_scaled']:
            aveWidth = rrrAveQuant[c]['widths'][0].getValue()  # in eV, like ENDF widths
            if theArgs.reaction == 'capture':
                expectedDOF = 100
            else:
                expectedDOF = 1.0
            plot_widthddist(rrr, aveWidth, expectedDOF, theArgs.L, theArgs.J, theArgs.reaction,
                            title_string=theArgs.reaction, scaled=(theArgs.plot == 'widthdist_scaled'),
                            legend_loc=theArgs.legendLoc, outfile=theArgs.outFile)
            exit()

        # Measures of level statisitics
        if theArgs.plot in ['Delta3', 'nnsd', 'cumlev', 'levlevcorr']:
            Es = get_resonance_energies(rrr, c)
            full_lsa = LevelSequenceAnalyzer(Es)
            if theArgs.plot == 'nnsd':
                plot_nnsd(full_lsa, theArgs.L, theArgs.J, title_string=titleString_incomming,
                          legend_loc=theArgs.legendLoc, outfile=theArgs.outFile)
            elif theArgs.plot == 'cumlev':
                plot_cumlev(full_lsa)
            elif theArgs.plot == 'Delta3':
                plot_Delta3(full_lsa, theArgs.L, theArgs.J, title_string=titleString_incomming,
                            legend_loc=theArgs.legendLoc, outfile=theArgs.outFile)
            elif theArgs.plot == 'levlevcorr':
                plot_levlevcorr(full_lsa, title_string=titleString_incomming, legend_loc=theArgs.legendLoc,
                                outfile=theArgs.outFile)
            exit()

        if theArgs.plot in ['Delta3_by_E', 'rho', 'fraction_missing']:
            Es = get_resonance_energies(rrr, c)
            lsa = []
            chunkSize = min(theArgs.chunkSize, len(Es))
            nChunks = int(len(Es) / chunkSize)
            for ich in range(nChunks):
                lsa.append(LevelSequenceAnalyzer(Es[ich * chunkSize:min((ich + 1) * chunkSize, len(Es))], Lmax=20))
            if theArgs.plot == 'Delta3_by_E':
                plot_Delta3_by_E(lsa, nChunks, chunkSize, Es, theArgs.D3L, L=c.l, J=c.J,
                                 title_string=titleString_reaction, legend_loc=theArgs.legendLoc,
                                 outfile=theArgs.outFile)
            elif theArgs.plot == 'rho':
                plot_rho(lsa, nChunks, chunkSize, Es, L=c.l, J=c.J, title_string=titleString_reaction,
                         legend_loc=theArgs.legendLoc, outfile=theArgs.outFile)
            elif theArgs.plot == 'fraction_missing':
                plot_fraction_missing(lsa, nChunks, chunkSize, Es, L=c.l, J=c.J, title_string=titleString_reaction,
                                      legend_loc=theArgs.legendLoc, outfile=theArgs.outFile)
            exit()

        # Average parameter energy dependence
        if theArgs.plot in ['avewidth', 'avespacing']:
            if theArgs.plot == 'avespacing':
                plot_avespacing(rrrAveQuant, urr, c, title_string=titleString_reaction, legend_loc=theArgs.legendLoc,
                                outfile=theArgs.outFile)
            if theArgs.plot == 'avewidth':
                plot_avewidth(rrrAveQuant, urr, c, theArgs.reaction, title_string=titleString_reaction,
                              legend_loc=theArgs.legendLoc, outfile=theArgs.outFile)
            exit()

        # Combine the average parameters into transmission coefficients
        if theArgs.plot in ['tc']:
            plot_transmission_coefficient(title_string=titleString_reaction, legend_loc=theArgs.legendLoc,
                                          outfile=theArgs.outFile)

    # Make the FULL report
    elif theArgs.fullReport is not None:

        # Sort channels into spin groups and set up space for plots
        spinGroups = collections.OrderedDict()
        if rrr is not None:
            for c in rrr.allChannels:
                LJ = (c.l, c.J)
                if c.eliminated:
                    continue  # don't want eliminated channel with fake spins
                if LJ not in spinGroups:
                    spinGroups[LJ] = collections.OrderedDict()
                    spinGroups[LJ]['level_report'] = collections.OrderedDict()
                    spinGroups[LJ]['errors'] = []
                    spinGroups[LJ]['warnings'] = []
                    spinGroups[LJ]['channels'] = collections.OrderedDict()
                if c.l == 0 and c.J == 1.5:  # FIXME: Do this test right!!!
                    spinGroups[LJ]['errors'].append( "Channel with bogus total angular momentum got by: " + str(c) )
                    if False:
                        raise ValueError("Channel with bogus total angular momentum got by:", str(c))
                spinGroups[LJ]['channels'][c] = collections.OrderedDict()

        # Generate all the plots for the reports
        plot_filename_template = theArgs.outFile + '_%s' + '.png'

        def get_plot_filename(kind, L, J, rxn=None):
            if rxn is None:
                return plot_filename_template % ('_'.join([kind, str(L), str(float(J))]))
            return plot_filename_template % ('_'.join([kind, str(L), str(float(J)), rxn]))

        titleString_incomming = '%s + %s' % (projectile, target)

        for LJ in spinGroups:

            c0 = list(spinGroups[LJ]['channels'].keys())[0]
            titleString_incoming = '%s + %s' % (target, projectile)
            Es = get_resonance_energies(rrr, c0)
            if len(Es) < 5:
                if c.isElastic:
                    spinGroups[LJ]['errors'].append("Potential scattering only channel: " + get_key_from_channel(c))
                if c.eliminated:
                    spinGroups[LJ]['errors'].append("Eliminated channel: " + get_key_from_channel(c))
                spinGroups[LJ]['errors'].append(
                    "Not enough levels (%i) to analyze for spin group L=%i, J=%s" % (len(Es), c.l, str(c.J)))
                continue

            full_lsa = LevelSequenceAnalyzer(Es)
            lsa = []
            chunkSize = min(theArgs.chunkSize, len(Es))
            nChunks = int(len(Es) / chunkSize)
            for ich in range(nChunks):
                lsa.append(LevelSequenceAnalyzer(Es[ich * chunkSize:min((ich + 1) * chunkSize, len(Es))], Lmax=20))

            # Average spacing
            spinGroups[LJ]['level_report']['filename_avespacing'] = get_plot_filename('avespacing', LJ[0], LJ[1])
            plot_avespacing(rrrAveQuant, urr, c0, title_string=titleString_incoming,
                            outfile=spinGroups[LJ]['level_report']['filename_avespacing'], clear=True)

            # Cummulative level distribution
            spinGroups[LJ]['level_report']['filename_cumlev'] = get_plot_filename('cumlev', LJ[0], LJ[1])
            plot_cumlev(full_lsa, title_string=titleString_incoming,
                        outfile=spinGroups[LJ]['level_report']['filename_cumlev'], clear=True)

            # Nearest neighbor spacing distribution
            spinGroups[LJ]['level_report']['filename_nnsd'] = get_plot_filename('nnsd', LJ[0], LJ[1])
            plot_nnsd(full_lsa, LJ[0], LJ[1], title_string=titleString_incoming,
                      outfile=spinGroups[LJ]['level_report']['filename_nnsd'], clear=True)

            # Delta3 using whole energy range
            if 2 * theArgs.D3L > len(Es):
                if spinGroups[LJ]['warnings'] is None:
                    spinGroups[LJ]['warnings'] = []
                spinGroups[LJ]['warnings'].append(
                    "Not enough levels to compute Delta3, need at least %i, got %i" % (2 * theArgs.D3L, len(Es)))
            else:
                spinGroups[LJ]['level_report']['filename_Delta3'] = get_plot_filename('Delta3', LJ[0], LJ[1])
                plot_Delta3(full_lsa, LJ[0], LJ[1], title_string=titleString_incoming,
                            outfile=spinGroups[LJ]['level_report']['filename_Delta3'], clear=True)

            # Spacing-spacing correlation, rho
            spinGroups[LJ]['level_report']['filename_rho'] = get_plot_filename('rho', LJ[0], LJ[1])
            plot_rho(lsa, nChunks, chunkSize, Es, L=LJ[0], J=LJ[1], title_string=titleString_incoming,
                     outfile=spinGroups[LJ]['level_report']['filename_rho'], clear=True)

            # Delta3 as function of energy
            if False:  # FIXME: Currently broken
                if 2 * theArgs.D3L > len(Es):
                    if spinGroups[LJ]['warnings'] is None:
                        spinGroups[LJ]['warnings'] = []
                    spinGroups[LJ]['warnings'].append(
                        "Not enough levels to compute Delta3(E), need at least %i, got %i" % (2 * theArgs.D3L, len(Es)))
                else:
                    spinGroups[LJ]['level_report']['filename_Delta3_by_E'] = get_plot_filename('Delta3_by_E', LJ[0], LJ[1])
                    plot_Delta3_by_E(lsa, nChunks, chunkSize, Es, theArgs.D3L, L=LJ[0], J=LJ[1],
                                     title_string=titleString_incoming,
                                     outfile=spinGroups[LJ]['level_report']['filename_Delta3_by_E'], clear=True)

            # Missing fraction
            spinGroups[LJ]['level_report']['filename_fraction_missing'] = get_plot_filename('fraction_missing', LJ[0],
                                                                                            LJ[1])
            plot_fraction_missing(lsa, nChunks, chunkSize, Es, L=LJ[0], J=LJ[1], skipDelta3=2 * theArgs.D3L > len(Es),
                                  title_string=titleString_incoming,
                                  outfile=spinGroups[LJ]['level_report']['filename_fraction_missing'], clear=True)

            for c in spinGroups[LJ]['channels']:
                reaction = get_simple_reaction_name(c.reaction)
                expectedDOF = 1.0
                if reaction == 'capture':
                    expectedDOF = 10.0
                std_reactions = {'elastic': 'el', 'capture': '$\\gamma$', 'fission': 'f', "total": 'tot'}
                if reaction in std_reactions:
                    titleString_reaction = '%s(%s, %s)' % (target, projectile, std_reactions[reaction])
                else:
                    prods = reaction.split(' + ')
                    titleString_reaction = '%s(%s, %s)%s' % (target, projectile, prods[0], prods[1])

                # Average widths
                spinGroups[LJ]['channels'][c]['filename_avewidth'] = get_plot_filename('avewidth', LJ[0], LJ[1],
                                                                                       reaction)
                spinGroups[LJ]['channels'][c]['filename_widthdist'] = get_plot_filename('widthdist', LJ[0], LJ[1],
                                                                                        reaction)
                plot_avewidth(rrrAveQuant, urr, c, reaction, title_string=titleString_reaction,
                              legend_loc='best', outfile=spinGroups[LJ]['channels'][c]['filename_avewidth'],
                              clear=True)

                # Width distribution
                if rrrAveQuant[c]['widths'] and rrrAveQuant[c]['widths'][0] is not None:
                    plot_widthddist(rrr, rrrAveQuant[c]['widths'][0].getValue(), expectedDOF, c.l, c.J, reaction,
                                    title_string=titleString_reaction, scaled=True, legend_loc='lower left',
                                    outfile=spinGroups[LJ]['channels'][c]['filename_widthdist'], clear=True)

        # Now do the main part of the report
        showMeta = False
        rep = get_resonance_report(theEvaluation)
        for k in rep:
            if 'URR' in k and len(rep[k]) == 0:
                rep[k]['warnings']=["No URR in this evaluation"]
            elif 'RRR' in k and len(rep[k]) == 0:
                rep[k]['warnings']=["No RRR in this evaluation"]
            else:
                continue

        # Spingroup part of the report
        rep['spin_groups'] = collections.OrderedDict()
        for LJ in spinGroups:
            rep['spin_groups'][str(LJ)] = collections.OrderedDict()
            rep['spin_groups'][str(LJ)]['title'] = "Spin Group L=%s, J=%s" % (str(LJ[0]), str(LJ[1]))

            rep['spin_groups'][str(LJ)]['secular_variation'] = {
                'title':"Secular variation of levels",
                'plots':[
                    spinGroups[LJ]['level_report'].get('filename_avespacing', None),
                    spinGroups[LJ]['level_report'].get('filename_cumlev', None)
                ]
            }
            rep['spin_groups'][str(LJ)]['local_fluctuations'] = {
                'title':"Local fluctuations in levels",
                'plots':[
                    spinGroups[LJ]['level_report'].get('filename_nnsd', None),
                    spinGroups[LJ]['level_report'].get('filename_rho', None),
                    spinGroups[LJ]['level_report'].get('filename_Delta3', None),
                    spinGroups[LJ]['level_report'].get('filename_Delta3_by_E', None)
                ]
            }
            rep['spin_groups'][str(LJ)]['missing_levels'] = {
                'title':"Assessment of fraction of missing levels",
                'plots':[
                    spinGroups[LJ]['level_report'].get('filename_fraction_missing', None)
                ]
            }
            rep['spin_groups'][str(LJ)]['channels'] = collections.OrderedDict()
            for c in spinGroups[LJ]['channels']:
                rep['spin_groups'][str(LJ)]['channels'][str(c)] = collections.OrderedDict()
                rep['spin_groups'][str(LJ)]['channels'][str(c)]['title']= "Details for channel %s" % get_key_from_channel(c)
                rep['spin_groups'][str(LJ)]['channels'][str(c)]['plots']:[
                    spinGroups[LJ]['channels'][c].get('filename_avewidth', None),
                    spinGroups[LJ]['channels'][c].get('filename_widthdist', None)
                ]
                if c.eliminated:
                    rep['spin_groups'][str(LJ)]['channels'][c]['note'] = \
                        "Note: capture widths likely won't show spread since we probably only know the " + \
                        "average width to any precision and we know there are many many open capture channels."

        if theArgs.outFile is not None:
            if theArgs.outFile.endswith('.json'):
                outFileName = theArgs.outFile
            else:
                outFileName = theArgs.outFile + '.json'
            with open(outFileName, mode='w') as outFile:
                import json
                json.dump(rep, outFile, cls=report.ComplexEncoder)

    else:
        exit()
