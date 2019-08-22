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
                print "WARNING: no matrix elements > 0. Switching to linear scale for z-axis"
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

