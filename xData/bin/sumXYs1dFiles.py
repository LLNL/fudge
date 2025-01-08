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

summaryDocString__xData = '''Reads 1d data from each file listed into an XYs1d instance, sums them and prints the sum.'''

description = '''
Reads 1d data from each file listed into an XYs1d instance, sums them and prints the sum. If the files do not 
have mutual domains, use the '--mutualify' option.
The data in each file must be in columns and currently only columns 1 and 2 of each file are summed.
The '#' character is treated as a comment character and, for each line, it and all remaining character
on the line are ignored. Also, blank lines, after the comment is removed, are ignored.
'''

epsDefault = 1e-6
args = None

parser = argparse.ArgumentParser(description=description)

parser.add_argument('files', type=pathlib.Path, nargs='*',              help='The list of files whose data are plotted.')
parser.add_argument('-m', '--mutualify', action='store_true',           help='If present, domains will be mutualified if needed.')
parser.add_argument('-l', '--lowerEpsilon', action='store', type=float, default=epsDefault,
                                                                        help='If curve domains are not mutual, this is the lowerEps used to mutualify domains. Default is %s' % epsDefault)
parser.add_argument('-u', '--upperEpsilon', action='store', type=float, default=epsDefault,
                                                                        help='If curve domains are not mutual, this is the upperEps used to mutualify domains. Default is %s' % epsDefault)

def readCurve(file):
    '''Reads in data from *file* and returns and XYs1d instance of it.'''

    with open(file) as fIn:
        lines = fIn.readlines()

    data = []
    for line in lines:
        line = line.split('#')[0]
        values = list(map(float, line.split()))
        if len(values) > 0:
            data.append([values[0], values[1]])

    return XYs1dModule.XYs1d(data=data)

def add(sum, curve):
    '''Adds *curve* to *sum* and mutualifies if requested and needed.'''

    if args.mutualify:
        sum, curve = sum.mutualify(args.lowerEpsilon, args.upperEpsilon, True, curve, args.lowerEpsilon, args.upperEpsilon, True)

    return sum + curve

if __name__ == '__main__':
    args = parser.parse_args()

    sum = XYs1dModule.XYs1d()
    for file in args.files:
        try:
            sum = add(sum, readCurve(file))
        except:
            raise
            print('Reading %s failed.' % file)
            continue

    print(sum.toString(), end='')
