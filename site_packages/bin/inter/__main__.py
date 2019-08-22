#!/usr/bin/env python
import os.path, sys, argparse
sys.path.append(os.path.split(__file__)[0]+os.sep+'..')
from inter import *


# -------------------------------------------------------------------------------
# Comamnd line parsing
# -------------------------------------------------------------------------------

def parse_arguments():
    '''Parse the command line'''
    parser = argparse.ArgumentParser(\
        description='''
            Super fancy replacement for legacy INTER format.

            This does everything that legacy INTER does (use the "--legacyINTERReport" option), and so much more.
        ''')
    # Required command line options
    parser.add_argument('ENDF', type=str, help='ENDF file(s) whose cross section you want to study.' )

    # Set output
    parser.add_argument('-o',   dest='outFile', default=None, type=str, help='Output to a file called OUTFILE, instead of printing to stdout.' )

    # Verbosity
    parser.add_argument('-v', dest='verbose', default=False,  action='store_true', help="Enable verbose output." )

    # Control covariance usage
    parser.add_argument('--noCovariance', dest='useCovariance', action='store_false', help="Do not use covariance when computing values")
    parser.add_argument('--useCovariance', dest='useCovariance', action='store_true', default=True, help="If an evaluation has them, use covariances when computing values (the default)")

    # Limit output to certain things, such as MT or angular momenta J, L, S, and therefore Pi
    parser.add_argument('--MT', type=int, default=None, help='If given, only work with this MT.' )
    parser.add_argument('--MTList', default=None, choices=['legacy', 'major', 'all'], help='')
    # parser.add_argument('--L', type=int, default=None, help="If given, only work with this L, orbital angular momentum of channel" )
    # parser.add_argument('--J', type=int, default=None, help="If given, only work with this J, total angular momentum of channel" )
    # parser.add_argument('--S', type=int, default=None, help="If given, only work with this S, spin of compound nucleus formed" )

    # Astrophysics metrics of cross section data
    parser.add_argument('--MACS', type=float, default=None, help="Compute MACS, give kT in keV as argument" )
    parser.add_argument('--ARR', type=float, default=None, help="Compute astrophysical reaction rate, give kT in keV as argument")

    # Nuclear enginering metrics of cross section data
    parser.add_argument('--RI', default=False, action='store_true', help="Compute resonance integral, cut off at 0.5 eV")
    parser.add_argument('--thermal', default=False, action='store_true', help="Compute thermal cross section")
    parser.add_argument('--Westcott',default=False, action='store_true', help="Compute Westcott factor")
    parser.add_argument('--ALF', default=False, action='store_true', help="(capture cs)/(fission cs) at thermal")
    parser.add_argument('--ETA', default=False, action='store_true', help="nubar*(fission cs)/(absorption cs) at thermal")

    # Integral metrics of cross section data
    parser.add_argument('--CfSpectAnalytic',default=False, action='store_true', help="Compute 252Cf spontaneous fission spectrum average, using analytic approximation of [FIXME]")
    parser.add_argument('--CfSpect',default=False, action='store_true', help="Compute 252Cf spontaneous fission spectrum average")
    parser.add_argument('--14MeV', dest='fourteen', default=False, action='store_true', help="Get 14 MeV point")
    parser.add_argument('--Godiva', default=False, action='store_true', help="Compute Godiva (HMF001) assembly spectrum average")
    parser.add_argument('--Jezebel', default=False, action='store_true', help="Compute Jezebel (PMF001) assembly spectrum average")
    parser.add_argument('--BigTen', default=False, action='store_true', help="Compute BigTen (IMF007) assembly spectrum average")
    parser.add_argument('--FUNDIPPE', default=False, action='store_true', help="Compute FUND-IPPE (FIXME) assembly spectrum average")

    # Popular neutron sources
    if TURNONNEUTRONNSOURCES: parser.add_argument('--ddSource', default=None, type=float, help="Compute spectrum average using d(d,n)3He reaction as a neutron source, argument is deuteron energy in MeV.")
    if TURNONNEUTRONNSOURCES: parser.add_argument('--dtSource', default=None, type=float, help="Compute spectrum average using d(t,n)4He reaction as a neutron source, argument is deuteron energy in MeV.")
    if TURNONNEUTRONNSOURCES: parser.add_argument('--ptSource', default=None, type=float, help="Compute spectrum average using p(t,n)3He reaction as a neutron source, argument is proton energy in MeV.")
    if TURNONNEUTRONNSOURCES: parser.add_argument('--pLiSource', default=None, type=float, help="Compute spectrum average using p(7Li,n)7Be reaction as a neutron source, argument is proton energy in MeV.")

    # Resonance metrics
    parser.add_argument('--scatteringRadius', default=False, action='store_true', help="FIXME")
    parser.add_argument('--neutronStrengthFunction', default=False, action='store_true', help="FIXME")
    parser.add_argument('--neutronPoleStrength', default=False, action='store_true', help="FIXME")
    parser.add_argument('--gammaStrengthFunction', default=False, action='store_true', help="FIXME")
    parser.add_argument('--averageWidth', default=False, action='store_true', help="Average resonance width per channel in both the RRR and URR")
    parser.add_argument('--meanLevelSpacing', default=False, action='store_true', help="Mean level spacing per channel in both the RRR and URR")
    parser.add_argument('--DysonMehtaDelta3', default=False, action='store_true', help="FIXME")
    parser.add_argument('--effectiveDOF', default=False, action='store_true', help="FIXME")
    parser.add_argument('--transmissionCoeff', default=False, action='store_true', help="Transmision coefficient per channel in both the RRR and URR, computed by Moldauer's sum rule")

    # Full report controls
    reportChoices=['astrophysics','engineering','integral','resonance','legacy','summary']
    if TURNONNEUTRONNSOURCES: reportChoices.append('nsources')
    parser.add_argument('--report', default=None, choices=reportChoices, help='Style of report to generate (Default: None)')
    parser.add_argument('--reportFormat', default='txt', choices=['html','txt','csv','json'], help="Format of output report (Default:txt)")

    return parser.parse_args()


# -------------------------------------------------------------------------------
# Main!!
# -------------------------------------------------------------------------------

if __name__ == "__main__":

    # Parse command line
    theArgs = parse_arguments()

    # Set the MTList
    if theArgs.report == 'legacy' or theArgs.MTList == 'legacy':
        theArgs.MTList = [1,2,102,18]
    elif theArgs.MT != None:
        theArgs.MTList = [theArgs.MT]
    elif theArgs.MTList == 'major':
        theArgs.MTList = [1,2,4,16,17,22,102,103,107,18]
    else:
        theArgs.MTList=None

    # Read in evaluation
    rep=EvaluationReport(theArgs.ENDF,theArgs.MTList)

    # Reaction table
    if theArgs.report is None or theArgs.report in ['astrophysics','engineering','integral','legacy','summary']:
        rep.get_reaction_metrics(theArgs)

    # ALF, ETA addendum
    if theArgs.report is None or theArgs.report in ['engineering','summary']:
        rep.get_global_metrics(theArgs)

    # Resonance table
    if theArgs.report in ['resonance','summary']:
        rep.get_resonance_metrics(theArgs)

    if theArgs.scatteringRadius:
        print computeScatteringRadius(rep.reactionSuite.getReaction('elastic').crossSection)

    # ---------------------------------
    # Generate report
    # ---------------------------------

    if theArgs.reportFormat=='txt': report = rep.text_report(theArgs)
    elif theArgs.reportFormat=='json': report = rep.json_report(theArgs)
    elif theArgs.reportFormat=='csv': report = rep.csv_report(theArgs)
    elif theArgs.reportFormat=='html': report = '\n'.join(rep.html_report(theArgs))
    else: raise ValueError("Output format %s unknown"%theArgs.reportFormat)

    # Output report
    if theArgs.outFile is None: print report
    else: open(theArgs.outFile,mode='w').write(report)
