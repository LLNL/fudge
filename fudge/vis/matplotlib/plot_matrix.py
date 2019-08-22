# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# <<END-copyright>>

import numpy
from matplotlib import pyplot
from matplotlib.colors import LogNorm

def plot_matrix( matrix, energyBoundariesX=None, energyBoundariesY=None, title="Matrix", xyTitle='Energy (MeV)',
        xylog=False, zlog=False, zRange=(), switchY=True, colorMap='jet', subplot=None ):
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
        @colorMap: select color range for z-axis. Options include 'jet', 'jet_r', 'correlation', etc.
        @subplot: matplotlib.axes.AxesSubplot instance (for drawing matrix inside larger figure)
        
        This prepares a 2-d matrix plot. Use pyplot.show() to view the result."""

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
        vmin, vmax = matrix.min(), matrix.max()
        if zlog:
            try: vmin = min( matrix[matrix>0] )
            except ValueError:
                print "WARNING: switching to linear scale for z-axis"
                zlog=False

    # z-axis log scale?
    if zlog: zopts = {'norm': LogNorm( vmin=vmin, vmax=vmax ) }
    else: zopts = {'vmin':vmin, 'vmax':vmax }

    # color map:
    zopts['cmap'] = colorMap
    if colorMap=='correlation': zopts['cmap'] = generateColorMap()

    if xylog:
        if energyBoundariesX is None: raise Exception("xylog option requires energy group boundary info")
        X,Y = numpy.lib.function_base.meshgrid(energyBoundariesX,energyBoundariesY)
        pyplot.pcolor( X,Y, matrix, **zopts )
        pyplot.setp( subplot, xscale='log', yscale='log' )
    elif energyBoundariesX is not None:
        pyplot.pcolor( energyBoundariesX, energyBoundariesY, matrix, **zopts )
    else:
        pyplot.imshow( matrix, **zopts )

    pyplot.title(title)
    pyplot.colorbar()
    pyplot.matplotlib.cm.jet.set_under('w')

def generateColorMap():
    # based on jet_r, but the center is white (for 0 correlation)
    import matplotlib
    cdict = {'blue': [(0.0, 0, 0),
            (0.34999999999999998, 0, 0),
            (0.499, 0.4828709677, 1),
            (0.501, 1, 0.4828709677),
            (0.65999999999999992, 1, 1),
            (0.89000000000000001, 1, 1),
            (1.0, 0.5, 0.5)],
        'green': [(0.0, 0, 0),
            (0.089999999999999969, 0, 0),
            (0.35999999999999999, 1, 1),
            (0.499, 1, 1),
            (0.501, 1, 1),
            (0.625, 1, 1),
            (0.875, 0, 0),
            (1.0, 0, 0)],
        'red': [(0.0, 0.5, 0.5),
            (0.10999999999999999, 1, 1),
            (0.33999999999999997, 1, 1),
            (0.499, 0.516129, 1),
            (0.501, 1, 0.516129),
            (0.65000000000000002, 0, 0),
            (1.0, 0, 0)]} 
    return matplotlib.colors.LinearSegmentedColormap('correlation',cdict,512)

