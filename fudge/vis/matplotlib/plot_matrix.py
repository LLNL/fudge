# <<BEGIN-copyright>>
# <<END-copyright>>

import numpy
from matplotlib import pyplot
from matplotlib.colors import LogNorm

def plot_matrix( matrix, energyBoundariesX=None, energyBoundariesY=None, title="Matrix", xyTitle='Energy (MeV)',
        xylog=False, zlog=False, zRange=(), switchY=True, subplot=None ):
    """Plot a matrix. Arguments are:
        @matrix (required): 2-dimensional numpy array (or equivalent)
        @energyBoundariesX: list of group boundary edges. If omitted, use matrix indices for axis
        @energyBoundariesY: "". If omitted, use same boundaries as x-axis
        @title
        @xyTitle: specify x (and optionally y) titles. May be string or tuple of two strings
        @xylog: if True, energyBoundariesX (and Y if different) must also be specified
        @zlog
        @zRange: for plotting a sub-range of z
        @switchY: if True, plot matrix[0,0] at upper left rather than lower left
        @subplot: matplotlib.axes.AxesSubplot instance (for drawing matrix inside larger figure)"""

    matrix = numpy.array( matrix )
    if matrix.ndim != 2:
        raise Exception("plot_matrix can only handle 2-d arrays!")

    if not subplot:
        fig = pyplot.figure()
        subplot = fig.add_subplot(111)

    if type(xyTitle) is str: xtitle = ytitle = xyTitle
    elif len(xyTitle)==2: xtitle, ytitle = xyTitle
    else:
        raise Exception("Can't understand xyTitle argument: %s" % xyTitle)
    pyplot.xlabel(xtitle)
    pyplot.ylabel(ytitle)

    if energyBoundariesY is None: energyBoundariesY = energyBoundariesX

    if energyBoundariesX is not None:
        lowX, highX = energyBoundariesX[0], energyBoundariesX[-1]
        lowY, highY = energyBoundariesY[0], energyBoundariesY[-1]
    else:
        lowX = lowY = 0
        highX, highY = matrix.shape

    pyplot.xlim( lowX, highX )
    if switchY: # switch y-axis: put low energy at top
        pyplot.ylim( highY, lowY )
    else:
        pyplot.ylim( lowY, highY )

    if zRange: vmin, vmax = zRange
    else:
        vmin = min( matrix[matrix>0] ); vmax = matrix.max()

    # z-axis log scale?
    if zlog: zopts = {'norm': LogNorm( vmin=vmin, vmax=vmax ) }
    else: zopts = {'vmin':vmin, 'vmax':vmax }

    if xylog:
        if energyBoundariesX is None: raise Exception("xylog option requires energy group boundary info")
        X,Y = numpy.lib.function_base.meshgrid(energyBoundariesX,energyBoundariesY)
        pyplot.pcolor( X,Y, matrix, **zopts )
        pyplot.setp( subplot, xscale='log', yscale='log' )
    elif energyBoundariesX is not None:
        pyplot.pcolor( energyBoundariesX, energyBoundariesY, matrix, **zopts )
    else:
        pyplot.pcolor( matrix, **zopts )

    pyplot.title(title)
    pyplot.colorbar()
    pyplot.matplotlib.cm.jet.set_under('w')


