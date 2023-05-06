#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse
import pathlib
import shutil

from xData import XYs1d as XYs1dModule

summaryDocStringxData = '''Reads 1d data from a file and convolutes its data with data from another file or a Gaussian.'''

description = '''
Reads 1d data from a file, and convolutes with data from another file if present or a Gaussian function. That is, performs
the convolution "f * g" where "f"  is the first file specified, and "g" is either the second file, if specified,
or a Guassian function with parameters specified via the options "--stdDev" and "--offset".
All curves are thinned to lin-lin accuracy as specified by the "--accuracy" option.

The results are printed to standard output.
'''

descriptionGaussian = '''
The parameters for the Guassian function of the form "amplitude * exp( -(x - offset)**2 / stdDev**2 / 2)".
This script picks the "amplitude" so that the Guassian function has normalization of 1.0.
If a second file is specified in the arguments, these parameters are ignored.
'''

parser = argparse.ArgumentParser(description=description)

parser.add_argument('data1', type=pathlib.Path,                                         help='The file containing data to convolute.')
parser.add_argument('data2', type=pathlib.Path, nargs='?',                              help='The file containing data to convolute.')
parser.add_argument('--accuracy', type=float, default=10.e-3,                           help='The results are thinned to this accuracy. Default is 1e-3.')
group = parser.add_argument_group(title='Gaussian function parameters', description=descriptionGaussian)
group.add_argument('-s', '--stdDev', type=float, default=1.0,                           help='The standard deviation of the Gaussian function. Default is 1.0.')
group.add_argument('-o', '--offset', type=float, default=0.0,                           help='The offset of the Gaussian function. Default is 0.0.')

def read(fileName):

    file = pathlib.Path(fileName)

    with file.open() as fIn:
        lines = fIn.readlines()

    data = []
    for line in lines:
        line = line.split('#')[0]
        if len(line) != 0:
            data.append(list(map(float,line.split())))

    return XYs1dModule.XYs1d(data=data)

if __name__ == '__main__':
    args = parser.parse_args()

    function1 = read(args.data1)
    function1 = function1.thin(args.accuracy)
    if args.data2 is not None:
        function2 = read(args.data2)
        function2 = function2.thin(args.accuracy)
    else:
        domainMin = -5 * args.stdDev
        domainMax =  5 * args.stdDev
        function2 = XYs1dModule.XYs1d(data=XYs1dModule.pointwiseXY_C.gaussian(args.accuracy, domainMin, domainMax, args.offset, args.stdDev, 1.0)).normalize()

    print(len(function1), len(function2))
    convolution = function1.convolute(function2).thin(accuracy=args.accuracy)
    print(convolution.toString())
