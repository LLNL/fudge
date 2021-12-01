#!/usr/bin/env python

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import collections
import os.path
import sys
import xData
import pqu
import argparse

sys.path.append(os.path.split(__file__)[0] + os.sep + '..')

from fudge.core.utilities import banner
from brownies.BNL.inter import *
from brownies.BNL.inter.report import ComplexEncoder
from brownies.BNL.utilities.io import *

def getPlottableReactions(reactionSuite, observable='crossSection'):
    """
    Retrieve references to all reactions with the observable specified.
    For now, the observable should be an attribute of the reaction.
    :param reactionSuite:
    :param observable:
    :return:
    """
    result = []
    for reactionList in reactionSuite.reactions, reactionSuite.sums.crossSectionSums, \
                        reactionSuite.sums.multiplicitySums, \
                        reactionSuite.fissionComponents, reactionSuite.productions:
        for r in reactionList:
            if hasattr(r, observable):
                result.append(r)
    return result


def getEvaluationMTs(reactionSuite, mtFilter=None):
    reactionList = getPlottableReactions(reactionSuite)
    if mtFilter is None:
        return [int(r.ENDF_MT) for r in reactionList]
    else:
        return [int(r.ENDF_MT) for r in reactionList if int(r.ENDF_MT) in mtFilter]


# -------------------------------------------------------------------------------
# Comamnd line parsing
# -------------------------------------------------------------------------------

def parse_arguments():
    """Parse the command line"""
    parser = argparse.ArgumentParser(
        description="""
            INTERrogate an evaluation!

            Super fancy replacement for legacy INTER code.

            This does everything that legacy INTER does (use the "--report legacy" option), but so much more.
        """)
    # Required command line options
    parser.add_argument('ENDF', type=str, help='ENDF file(s) whose cross section you want to study.')

    # Set output
    parser.add_argument('-o', dest='outFile', default=None, type=str,
                        help='Output to a file called OUTFILE, instead of printing to stdout.')

    # Verbosity
    parser.add_argument('-v', dest='verbose', default=False, action='store_true', help="Enable verbose output.")

    # Control covariance usage
    parser.add_argument('--noCovariance', dest='useCovariance', action='store_false',
                        help="Do not use covariance when computing values")
    parser.add_argument('--useCovariance', dest='useCovariance', action='store_true', default=True,
                        help="If an evaluation has them, use covariances when computing values (the default)")

    # Limit output to certain things, such as MT or angular momenta J, L, S, and therefore Pi
    parser.add_argument('--MT', type=int, default=None, help='If given, only work with this MT.')
    parser.add_argument('--MTList', default='major', choices=['legacy', 'major', 'all'],
                        help='List of MTs to do in the table (overrides the --MT option)')

    # Astrophysics metrics of cross section data
    parser.add_argument('--MACS', type=float, default=None, help="Compute MACS, give kT in keV as argument")
    parser.add_argument('--ARR', type=float, default=None,
                        help="Compute astrophysical reaction rate, give kT in keV as argument")

    # Nuclear enginering metrics of cross section data
    parser.add_argument('--RI', default=False, action='store_true',
                        help="Compute resonance integral, cut off at 0.5 eV")
    parser.add_argument('--RI_Ecut', default=0.5, type=float,
                        help="Cut off energy for the resonance integral in eV, (Default is 0.5 eV)")
    parser.add_argument('--thermal', default=False, action='store_true', help="Compute thermal cross section")
    parser.add_argument('--Westcott', default=False, action='store_true', help="Compute Westcott factor")
    parser.add_argument('--ALF', default=False, action='store_true', help="(capture cs)/(fission cs) at thermal")
    parser.add_argument('--ETA', default=False, action='store_true',
                        help="nubar*(fission cs)/(absorption cs) at thermal")

    # Integral metrics of cross section data
    parser.add_argument('--CfSpectAnalytic', default=False, action='store_true',
                        help="Compute 252Cf spontaneous fission spectrum average, "
                             "using analytic approximation of [FIXME]")
    parser.add_argument('--CfSpectENDF', default=False, action='store_true',
                        help="Compute 252Cf spontaneous fission spectrum average")
    parser.add_argument('--14MeV', dest='fourteen', default=False, action='store_true', help="Get 14 MeV point")
    parser.add_argument('--Godiva', default=False, action='store_true',
                        help="Compute Godiva (HMF001) assembly spectrum average")
    parser.add_argument('--Jezebel', default=False, action='store_true',
                        help="Compute Jezebel (PMF001) assembly spectrum average")
    parser.add_argument('--BigTen', default=False, action='store_true',
                        help="Compute BigTen (IMF007) assembly spectrum average")

    # Resonance metrics
    parser.add_argument('--scatteringRadius', default=False, action='store_true',
                        help="R', the potential scattering radius")
    parser.add_argument('--poleStrength', default=False, action='store_true',
                        help="The neutron/photon/whatever pole strength (related to strength functions, "
                             "but grounded in the R-matrix)")
    parser.add_argument('--averageWidth', default=False, action='store_true',
                        help="Average resonance width per channel in both the RRR and URR")
    parser.add_argument('--meanLevelSpacing', default=False, action='store_true',
                        help="Mean level spacing per channel in both the RRR and URR")
    parser.add_argument('--effectiveDOF', default=False, action='store_true',
                        help="Effective DOF from Porter Thomas distribution fit to the resonance width "
                             "distributions within each channel in both the RRR and URR")

    # Plots of resonance-related things
    parser.add_argument('--transmissionCoeff', default=False, action='store_true',
                        help="Transmision coefficient per channel in both the RRR and URR, "
                             "computed by Moldauer's sum rule")
    parser.add_argument('--DysonMehtaDelta3', default=False, action='store_true', help="FIXME")

    # Average cross section over a window (esentially grouping it)
    parser.add_argument('--windowCenter', default=None, type=float,
                        help="Window's Center: Energy in keV that, if given along with windowWidth, "
                             "compute the average of the specified reactions over that window")
    parser.add_argument('--windowWidth', default=None, type=float,
                        help="Window's Total width: Energy in keV that, if given along with windowCenter, "
                             "compute the average of the specified reactions over that window")

    # Average over a user-defined spectrum
    parser.add_argument('--userDefinedSpecAve', default=None, type=str,
                        help="Compute spectrum average using the spectrum defined in this file/key")
    parser.add_argument('--otherSpecAve', default=None, type=str,
                        help="Compute spectrum average using the spectrum corresponding to this key")
    parser.add_argument('--listDefinedSpectra', default=False, action='store_true',
                        help="List the pre-defined spectra defined, including ones not in the default reports")

    # Full report controls
    reportChoices = ['astrophysics', 'engineering', 'integral', 'resonance', 'legacy', 'fission', 'summary', 'theworks']
    parser.add_argument('--report', default=None, choices=reportChoices,
                        help='Style of report to generate (Default: None)')
    parser.add_argument('--reportFormat', default='txt', choices=['html', 'txt', 'csv', 'json'],
                        help="Format of output report (Default:txt)")

    return parser.parse_args()


# -------------------------------------------------------------------------------
# Main!!
# -------------------------------------------------------------------------------

if __name__ == "__main__":

    # Parse command line
    theArgs = parse_arguments()
    if theArgs.listDefinedSpectra:
        print('Predefined spectra are (in addition to MACS, RI and "window") are')
        print("=================================================================")
        for k in spectra.spectra.keys():
            print(k, ':', spectra.spectra[k]['reference'])
            d = spectra.spectra[k]['description']
            if type(d) in [str]:
                print('   ', d)
            else:
                for l in d:
                    print('   ', l.rstrip())
        exit()

    # Read in evaluation
    theEvaluation = read_evaluation(theArgs.ENDF)

    # Set the MTList
    if theArgs.MT is not None:
        theArgs.MTList = [theArgs.MT]
    elif theArgs.report == 'legacy' or theArgs.MTList == 'legacy':
        theArgs.MTList = [1, 2, 102, 18]
    elif theArgs.MTList == 'major':
        theArgs.MTList = [1, 2, 4, 16, 17, 22, 51, 102, 103, 107, 18]
    elif theArgs.MTList == 'all':
        theArgs.MTList = getEvaluationMTs(theEvaluation['reactionSuite'])
    else:
        raise ValueError("Unknown MTList option %s" % theArgs.MTList)

    # set up the report
    rep = collections.OrderedDict()
    rep['Main'] = report.getEvaluationReport(theEvaluation['reactionSuite'])

    # Big reports
    if theArgs.report is not None:
        if theArgs.report in ['astrophysics', 'summary', 'theworks']:
            rep.update(report.getReactionAstrophysicsDataTable(theEvaluation['reactionSuite'], MTList=theArgs.MTList,
                                                               useCovariance=theArgs.useCovariance))
        if theArgs.report in ['engineering', 'summary', 'theworks']:
            rep.update(report.getReactionEngineeringDataTable(theEvaluation['reactionSuite'], MTList=theArgs.MTList,
                                                              useCovariance=theArgs.useCovariance))
        if theArgs.report in ['integral', 'summary', 'theworks']:
            rep.update(report.getReactionAIntegralDataTable(theEvaluation['reactionSuite'], MTList=theArgs.MTList,
                                                            useCovariance=theArgs.useCovariance))
        if theArgs.report in ['legacy', 'summary', 'theworks']:
            rep.update(report.getReactionLegacyINTERDataTable(theEvaluation['reactionSuite'], MTList=theArgs.MTList,
                                                              useCovariance=theArgs.useCovariance))
        if theArgs.report in ['resonance', 'theworks']:
            rep.update(report.getResonanceReport(theEvaluation['reactionSuite'], title="Resonances"))
            metricMenu = collections.namedtuple('metricMenu',
                                                'effectiveDOF poleStrength')
            rep.update({"Channels": report.getChannelDataTable(theEvaluation['reactionSuite'],
                                                               metricMenu(effectiveDOF=True, poleStrength=True))})
        if theArgs.report in ['fission', 'summary', 'theworks']:
            metricMenu = collections.namedtuple('metricMenu',
                                                "ALF ETA")  # emulates results from ArgumentParser.parse_args() function
            rep.update(report.getFissionMetricReport(theEvaluation['reactionSuite'], metricMenu(ALF=True, ETA=True),
                                                     title="Fission metrics", useCovariance=theArgs.useCovariance))

    # Plots
    elif theArgs.transmissionCoeff or theArgs.DysonMehtaDelta3:
        report.getChannelPlots(theEvaluation['reactionSuite'], theArgs)

    # Reports, ala' carte
    else:
        rep["Results"] = report.getReactionDataTable(theEvaluation['reactionSuite'],
                                                     theArgs,
                                                     title="Reactions",
                                                     MTList=theArgs.MTList,
                                                     useCovariance=theArgs.useCovariance,
                                                     RI_Ecut=pqu.PQU.PQU(theArgs.RI_Ecut, 'eV'))

    # Output report
    if theArgs.outFile is not None:
        outFile = open(theArgs.outFile, mode='w')
    else:
        import io
        outFile = io.StringIO()

    if theArgs.reportFormat == 'json':
        import json
        repOut = json.dumps({k:json.loads(rep[k].json_report()) for k in rep}, cls=ComplexEncoder)
        outFile.write(repOut)
    else:
        for k in rep:
            if theArgs.reportFormat == 'txt':
                if isinstance(rep[k], report.Report):
                    repOut = rep[k].text_report()
                elif isinstance(rep[k], xData.table.table):
                    repOut = '\n'.join(rep[k].toStringList())
                else:
                    repOut = str(rep[k])
            elif theArgs.reportFormat == 'csv':
                repOut = rep[k].csv_report()
            elif theArgs.reportFormat == 'html':
                repOut = '\n'.join(rep[k].html_report())
            else:
                raise ValueError("Output format %s unknown" % theArgs.reportFormat)

        outFile.write(banner(k) + '\n')
        outFile.write(repOut)
    if hasattr(outFile, 'getvalue'):
        print(outFile.getvalue())
