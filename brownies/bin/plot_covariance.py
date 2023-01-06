#! /usr/bin/env python

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse
from brownies.BNL.plot_evaluation import plotio as io

# Process command line options
parser = argparse.ArgumentParser(description='Plot the covariance data in an ENDF file')
parser.add_argument('inFile', type=str, help='The ENDF file you want to check.')
parser.add_argument('MT', type=int, help='The MT you want to plot.')
parser.add_argument('--crossMT', type=int, default=None, help='The MT of a cross correlation you want to plot')
parser.add_argument('--MF', type=int, default=33, help='MF to plot (Default is 33 for cross sections)')
parser.add_argument('--abs', dest='abs', default=None, action='store_true', help='Change covariance to absolute')
parser.add_argument('--rel', dest='rel', default=None, action='store_true', help='Change covariance to relative')
parser.add_argument('--corr', dest='corr', default=None, action='store_true',
                    help='Change covariance to a correlation matrix')
parser.add_argument("--skipBadData", action="store_true", default=False,
                    help="skip bad data, rather than throw an exception, when reading an ENDF file")
parser.add_argument('-v', dest='verbose', default=False, action='store_true', help='Enable verbose output')
parser.add_argument('-l', dest='list', default=False, action='store_true', help='List covariance data in file')
args = parser.parse_args()

# if True:
try:
    results = io.readEvaluation(args.inFile, verbose=args.verbose, skipBadData=args.skipBadData)
    evaluation = results[0]
    covariances = results[1]
except Exception as err:
    print('WARNING: ENDF READ HALTED BECAUSE ' + str(err))
    exit()

if args.list:
    print("\n\nList of covariance data")
    print("=======================")
    for c in covariances.covarianceSections:
        row = c.rowData.ENDF_MFMT
        if c.columnData is not None:
            col = c.columnData.ENDF_MFMT
        else:
            col = row
        print(c, row, col)
    exit()

made_a_plot = False

for c in covariances.parameterCovariances:
    if args.MF == 32:
        c.evaluated.plot(title='(%i,%i) x (%i,%i)' % (args.MF, args.MT, args.MF, args.MT))
        made_a_plot = True

    elif args.MF == 30:
        raise NotImplementedError('(%i,%i) x (%i,%i)' % (args.MF, args.MT, args.MF, args.MT))

for c in covariances.covarianceSections:
    if hasattr(c, 'rowData') and c.rowData.ENDF_MFMT == '%i,%i' % (args.MF, args.MT):
        if args.crossMT is None or \
                (c.columnData is not None and c.columnData.ENDF_MFMT == '%i,%i' % (
                args.MF, args.crossMT)) or \
                (c.columnData is None and (args.MF, args.MT) == (args.MF, args.crossMT)):
            otherMT = None
            if args.crossMT is not None:
                otherMT = args.crossMT
            elif c.columnData is None:
                otherMT = args.MT
            else:
                otherMT = int(c.columnData.ENDF_MFMT.split(',')[-1])

            c2 = c.evaluated.toCovarianceMatrix()
            if args.abs:
                c2.toAbsolute().plot(title='(%i,%i) x (%i,%i)' % (args.MF, args.MT, args.MF, otherMT))
            elif args.rel:
                c2.toRelative().plot(title='(%i,%i) x (%i,%i)' % (args.MF, args.MT, args.MF, otherMT))
            elif args.corr:
                c2.toAbsolute().toCorrelationMatrix().plot(
                    title='(%i,%i) x (%i,%i)' % (args.MF, args.MT, args.MF, otherMT))
            else:
                c2.plot(title='(%i,%i) x (%i,%i)' % (args.MF, args.MT, args.MF, otherMT))
            made_a_plot = True

if not made_a_plot:
    if args.crossMT is not None:
        otherMT = args.crossMT
    elif c.columnData is None:
        otherMT = args.MT
    raise KeyError("(MF0,MT0)x(MF1,MT1) = (%s,%s)x(%s,%s) not found" % (args.MF, args.MT, args.MF, otherMT))
