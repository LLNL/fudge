#! /usr/bin/env python

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
import os
import argparse

binDir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(binDir))


def process_args():
    parser = argparse.ArgumentParser(description='Translate an ENDF file to the new GNDS format')
    parser.set_defaults(verbose=True)
    parser.add_argument("-v", action="store_true", dest='verbose', help="enable verbose output")
    parser.add_argument("-q", action="store_false", dest='verbose', help="disable verbose output")
    parser.add_argument("-ifpy", action="store_true", dest='ifpy', help="show independent (prompt) fission yields")
    parser.add_argument("-cfpy", action="store_true", dest='cfpy', help="show cummulative (delayed) fission yields")
    parser.add_argument("--za", default=None, type=int, dest='za', help='the za to attempt to learn about (optional)')
    parser.add_argument("--state", default=0, type=int, dest='state',
                        help='the excitation of the za to look for (default=0 for g.s.)')
    parser.add_argument('--fudgepath', default='', dest='fudgepath', type=str,
                        help="set the path to fudge (use if not already in your PYTHONPATH)")
    parser.add_argument("inFile", type=str, help="the input endf file (should have '.endf' extension)")
    return parser.parse_args()


def grok_yields(d):
    isoList = sorted(d.keys())
    print(isoList)


def print_sum_by_A(d, iE=0):
    result = 200 * [0]
    dresult = 200 * [0]
    for key in d:
        A = int(str(key[0])[-3:])
        result[A] += d[key][iE][1]
        dresult[A] += d[key][iE][2]
    # return result, dresult
    for i in range(200)[75:175]:
        print(str(i).rjust(4), int(result[i] * 200) * '*')


def print_ZA_by_E(d, ZA, state=0):
    for E, Y, DY in d[(ZA, state)]:
        print("{:<11} {:>11} +/- {:<5}%".format(E / 1e6, Y * 100.0, (DY / Y) * 100.0))


if __name__ == "__main__":
    args = process_args()
    if args.fudgepath != '' and args.fudgepath not in sys.path: sys.path.append(args.fudgepath)
    from brownies.legacy.converting.endfFileToGNDSMisc import parseENDFByMT_MF
    from brownies.legacy.converting.endfFileToGNDS import endfFileToGNDS

    header, MAT, theData = parseENDFByMT_MF(args.inFile)

    x, c, i = endfFileToGNDS(args.inFile)

    # ifpy = ( 454, 8 ), independent (prompt) fission yields
    # cfpy = ( 459, 8 ), cumulative (delayed) fission yields

    if args.ifpy:
        theData = i.independentFissionYields
    if args.cfpy:
        theData = i.cumulativeFissionYields

    if args.za is not None:
        print_ZA_by_E(theData, args.za, state=args.state)
    else:
        print_sum_by_A(theData, iE=0)
#    else: grok_yields( theData )
