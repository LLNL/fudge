# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

description = '''
    This file creates an interactive gnuplot session that plots two column (2d) data from multiple files on
    the same window.  Each file must contain at least 2 columns of numbers. One column for the x-data and
    another for the y-data. For each file, the x- and y-data can be selected using the xColumn and yColumn
    keywords (see below) with the defaults being xColumn = 1 and yColumn = 2. This file is called as,

    python multiplotFromAscii.py [Parameters] files <datafile1> [keyword/value pairs] <datafile2> [keyword/value pairs]
'''

import re
import argparse

from xData.interactivePlot import multiplot as interactivePlotModule

parser = argparse.ArgumentParser(description=description)
parser.add_argument('file', nargs='+', type=str, help='List of files with the two-column plotting data')
parser.add_argument('--xLabel', metavar='xLabel', type=str, default='x', help='Plot x-label')
parser.add_argument('--yLabel', metavar='yLabel', type=str, default='y', help='Plot y-label')
parser.add_argument('--title', metavar='title', type=str, default='', help='Plot title')

args = parser.parse_args()


def readASCII(_filename):
    _twoColumnDataRegex = re.compile(r'\s*(\S+)\s+(\S+)\s*$')
    _twoColumnData = [[], []]
    _plotLabel = _filename

    with open(_filename) as _fileObject:
        # read header
        _fileLine = _fileObject.readline()
        while not _twoColumnDataRegex.match(_fileLine):
            _fileLine = _fileObject.readline()

        # read 2-column data
        while _twoColumnDataRegex.match(_fileLine):
            _x, _y = _twoColumnDataRegex.findall(_fileLine)[0]
            _twoColumnData[0].append(float(_x))
            _twoColumnData[1].append(float(_y))
            _fileLine = _fileObject.readline()

    return {_plotLabel: _twoColumnData}


# read data from ASCII file
plotData = {}
for filename in args.file:
    plotData.update(readASCII(filename))

# axes min/max values
minMax = {'xMin': None, 'xMax': None, 'yMin': None, 'yMax': None}
for plotLabel in plotData.keys():
    minMax['xMin'] = min(plotData[plotLabel][0]) if minMax['xMin'] is None else \
        min(minMax['xMin'], min(plotData[plotLabel][0]))
    minMax['xMax'] = max(plotData[plotLabel][0]) if minMax['xMax'] is None else \
        max(minMax['xMax'], max(plotData[plotLabel][0]))
    minMax['yMin'] = min(plotData[plotLabel][1]) if minMax['yMin'] is None else \
        min(minMax['yMin'], min(plotData[plotLabel][0]))
    minMax['yMax'] = max(plotData[plotLabel][0]) if minMax['yMax'] is None else \
        max(minMax['yMax'], max(plotData[plotLabel][0]))

# plot attributes
plotAttributes = {'title': args.title, 'xLabel': args.xLabel, 'yLabel': args.yLabel}
for key in minMax.keys():
    plotAttributes[key] = str(minMax[key])

# generate plot
interactivePlotModule.MultiPlotWithPyQt5(plotAttributes, plotData)
