#! /usr/bin/env python
import argparse
from site_packages.BNL.data_io import readEvaluation
from site_packages.BNL.cs_metrics import *

"""
2 Jan 2015
David Brown


The legacy NNDC data testing code INTER is a valuable part of the ENDF quality assurance process.
It is used to compute various integral and other metrics of ENDF cross section data.
This code is a massive update of the legacy INTER code, for use as a standalone code or as part of a data testing framework such as ADVANCE.

As of Jan 2015, this code is mostly a framework of a code that is scheduled to be written in FY16 as part of the BNL NCSP work plan.
"""


# Process command line options
parser = argparse.ArgumentParser(description='Compute cross section metrics of an ENDF evaluation')
parser.add_argument('inFile', type=str, help='The input file to work on' )
parser.add_argument('-v', dest='verbose', default=False, action='store_true', help='Enable verbose output' )
parser.add_argument('--reportStyle',                        choices=['inter','html','csv'], type=str, default=None, help='write me!')
parser.add_argument('--computeResonanceMetrics',            default=False, action='store_true', help='write me!')
parser.add_argument('--computeThermalScatteringParameters', default=False, action='store_true', help='write me!')
parser.add_argument('--computeNuclearEngineeringMetrics',   default=False, action='store_true', help='write me!')
parser.add_argument('--computeAstrophysicalMetrics',        default=False, action='store_true', help='write me!')
parser.add_argument('--computeIntegralValidationMetrics',   default=False, action='store_true', help='write me!')

args = parser.parse_args()

inEval = readEvaluation( args.inFile, verbose=args.verbose, skipBadData=True )
print '\n\n'


# Resonance physics metrics
if args.computeResonanceMetrics:
    # For comparison to Said Mughabghab systematics.  Works on (n,el) only
    computeR0ScatteringRadius( xs, covariance=None )
    computeR1ScatteringRadius( xs, covariance=None )
    computeS0NeutronStrengthFunction() # write me!
    computeS1NeutronStrengthFunction() # write me!
    computeS2NeutronStrengthFunction() # write me!

    # Other RR metrics
    # these two work on (n,el)
    computeSWaveScatteringLength( xs, covariance=None )
    computePWaveScatteringLength( xs, covariance=None )
    # this one works on (n,g)
    computeGammaRayStrengthFunction() # write me!

    # Stuff that aids in setting up URR
    computeResonanceAverageWidth() # write me!
    computeMeanLevelSpacing() # write me!
    computeDelta3Statistic() # write me!
    computeFStatistic() # write me!
    computeDOFestimate() # write me!


# Thermal scattering parameters:
if args.computeThermalScatteringParameters:
    computeFreeCoherentScatteringLength() # write me!
    computeFreeIncoherentScatteringLength() # write me!
    computeBoundCoherentScatteringLength() # write me!
    computeCoherentScatteringCrossSection() # write me!
    computeIncoherentScatteringCrossSection() # write me!


# Nuclear engineering metrics, many of these are in EXFOR
if args.computeNuclearEngineeringMetrics:
    for rxn in nonThresholdReactions: # write me!
        computeRI( xs, Ec='0.5 eV', covariance=None )
        computeRoomTempCS( xs, ERT='0.0235339347844 eV', covariance=None )
        computeWestcottFactor( xs, a=None, covariance=None )
    computeETA() # write me!
    computeALF() # write me!


# Astrophysical metrics, Boris Pritychenko and KaDonis have tables of these
if args.computeAstrophysicalMetrics:
    for rxn in nonThresholdReactions: # write me!
        computeMACS( xs, kT, a=None, covariance=None )
        computeAstrophysicalReactionRate( xs, kT, covariance=None )
        computeStellerEnhancementFactor( xs, covariance=None )


# Integral validation metrics
if args.computeIntegralValidationMetrics:
    for rxn in allReactions: # write me!
        computeFourteenMeVPoint( xs, E14='14.2 MeV', covariance=None )
        computeFourteenMeVPoint( xs, E14='14 MeV', covariance=None ) # for comparison with legacy INTER
        computeCf252SpectrumAve( xs, covariance=None ) # uses spectrum parameterization in legacy INTER
        computeEvaluatedCf252SpectrumAve( xs, covariance=None ) # uses neutron standards version of spectrum (write me!)
        computeGodivaSpectrumAve( xs, covariance=None ) # write me!
        computeJezebelSpectrumAve( xs, covariance=None ) # write me!
        computeBigTenSpectrumAve( xs, covariance=None ) # write me!
        computeIPPESpectrumAve( xs, covariance=None ) # write me!


# Inter style-report
# Here is an example:
'''

                             PROGRAM INTER VERSION 8.07
                             --------------------------

 Selected Integrations of ENDF File 3 and File 10 Cross Sections

 Thermal cross section : Sig(2200) = Sig(Eth)
 Thermal energy        : Eth=  2.53000E-02 (eV)

 Ezero cross section   : Sig(Ezero)
 Ezero energy (input)  : E0 =  2.53000E-02 (eV)

 Maxwellian average    : Avg-Sigma = (2/sqrt(Pi)) Intg[E1:E2] Sig(E) Phi_m(E) dE / Intg[E1:E2] Phi_m(E) dE
 Maxwellian spectrum   : Phi_m(E)  = (E/E0^2) exp(-E/E0)
 Spectrum Temperature  : E0 =  2.53000E-02 (eV)
 Integration Limits    : E1 =  1.00000E-05 (eV)  E2 =  1.00000E+01 (eV)
 Integral of Spectrum       =  1.00000E+00

 Westcott g-factor     : G-fact = 2/sqrt(Pi)  Avg-Sigma / Sig(2200)

 Resonance Integral    : Res.Integ = Intg[E3:E4] Sig(E)/E dE
 Integration Limits    : E3 =  5.00000E-01 (eV)  E4 =  1.00000E+05 (eV)
 Integral of Spectrum       =  1.22061E+01

 Fiss.spect. average   : Sig(Fiss) = Intg[E1:E2] Sig(E) Phi_f(E) dE / Intg[E1:E2] Phi_f(E) dE
 Fission spectrum      : Phi_f(E)  = sqrt(E/Ef)/Ef exp(-E/E0)
 Spectrum Temperature  : Ef =  1.35000E+06 (eV)
 Integration Limits    : E1 =  1.00000E+03 (eV)  E2 =  2.00000E+07 (eV)
 Integral of Spectrum       =  1.00000E+00

 E14 cross-section     : Sig(E14)
 Selected Energy       : E14 =  1.40000E+07 eV


   Z   A LISO  LFS  MT  Reaction    Sig(2200)   Sig(Ezero)  Avg-Sigma  G-fact   Res Integ   Sig(Fiss)    Sig(E14)   MAT
 --- ------------- ---  ---------- ----------- ----------- ---------- -------  ----------- ----------- ----------- ----
  23  50             1  Total      5.22940E+01 5.22940E+01 5.3287E+01 1.01900  3.93038E+02 3.84716E+00 2.36499E+00 2325
  23  50             2  Elastic    7.62005E+00 7.62005E+00 8.5976E+00 1.12829  3.33773E+02 2.89589E+00 1.07771E+00 2325
  23  50             4  Inelas     0.00000E+00 0.00000E+00 0.0000E+00 0.00000  0.00000E+00 9.40349E-01 4.38366E-01 2325
  23  50            16  n,2n       0.00000E+00 0.00000E+00 0.0000E+00 0.00000  0.00000E+00 4.67801E-04 6.55028E-01 2325
  23  50           102  n,gamma    4.46732E+01 4.46732E+01 4.4689E+01 1.00036  5.92615E+01 1.34565E-03 9.35782E-04 2325
  23  50           103  n,p        7.10000E-04 7.10000E-04 7.1028E-04 1.00040  3.61196E-03 8.72171E-03 7.86943E-02 2325
  23  50           104  n,d        0.00000E+00 0.00000E+00 0.0000E+00 0.00000  0.00000E+00 3.89319E-06 6.01921E-03 2325
  23  50           105  n,t        0.00000E+00 0.00000E+00 0.0000E+00 0.00000  0.00000E+00 3.61629E-08 6.02009E-05 2325
  23  50           106  n,He3      0.00000E+00 0.00000E+00 0.0000E+00 0.00000  0.00000E+00 8.17049E-11 5.42660E-11 2325
  23  50           107  n,alpha    0.00000E+00 0.00000E+00 0.0000E+00 0.00000  3.89144E-14 2.30208E-04 3.72632E-02 2325'''
if args.reportStyle=='inter': pass

# Full HTML report, useful for inclusion in ADVANCE.
# Would be better if had comparison tables from EXFOR and Boris Pritychenko.
# This maybe should have graphical and tabular variants
#
if args.reportStyle=='html': pass
if args.reportStyle=='csv': pass