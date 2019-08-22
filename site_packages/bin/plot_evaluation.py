#! /usr/bin/env python

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

from __future__ import print_function
import collections, sys, os.path, traceback
try:                import argparse
except ImportError: from fudge.core.utilities import argparse

sys.path.append(os.path.split(__file__)[0]+os.sep+'..')

from site_packages.BNL.plot_evaluation import *
from site_packages.BNL.plot_evaluation import plotio
from site_packages.BNL.plot_evaluation import plotstyles


#---------------------------------------------------
# Set up the command line parser
#---------------------------------------------------
def parseArgs():
    parser = argparse.ArgumentParser(description='Plot nuclear data from an ENDF or GNDS file')

    # Required things so we know what to plot
    parser.add_argument('mt',       metavar='mt', type=int, help='MT of the cross section to plot.  If set to 0, will try to make all plots for all open channels (requires -o option too)' )
    parser.add_argument('endf',     metavar='endf', type=str, nargs='+', help='ENDF file(s) whose cross section you want to plot.  Use "None" for no input file.  AMPX BOF files and tables of X,Y(,dY) pairs may work too.' )

    # Misc FUDGE controls
    parser.add_argument('--skipBadData', default=False, action='store_true', help='Skip bad data in an ENDF file that might otherwise cause Fudge to halt (Default: False)')
    parser.add_argument('--continuumSpectraFix', default=False, action='store_true', help='Apply a continuum spectra fix for evalautions with messed up gamma spectra, otherwise causes Fudge to halt (Default: False)')
    parser.add_argument("--skipCovariances", action="store_true", default=False, help="skip the covariance, if present")
    parser.add_argument("--verboseWarnings", action="store_true", default=False, help="print verbose warnings")
    parser.add_argument("--printBadNK14", action="store_true", default=False, help="print bad NK's if found")
    parser.add_argument("--ignoreBadDate", action="store_true", default=False, help="ignore malformed ENDF dates")
    parser.add_argument("--acceptBadMF10FissionZAP", action="store_true", default=False, help="allow MF=10 MT=18 IZAP=0")
    parser.add_argument("--traceback", action="store_true", default=False, help="print traceback on exception")

    # Plot output file controls
    parser.add_argument('-o',       dest='outFile', default=None, type=str, help='Output file for plot (disables interactive plotting)' )

    # Override the target/projectile/product of interest
    parser.add_argument('--target',     default=None, type=str, help="The target nucleus, given in GNDS notation, e.g. 'Pu239' (Default is None which means to take it from the first ENDF file)" )
    parser.add_argument('--projectile', default=None, type=str, help="The projectile particles, given in GNDS notation, e.g. 'n' (Default is None which means to take it from the first ENDF file)" )
    parser.add_argument('--product',    default=None, type=str, help="The product particle of interest, given in GNDS notation, e.g. 'n' (Default is None which means to take the first emitted particle for this observable)" )

    # Resolved resonance region reconstruction controls
    parser.add_argument('--enableRRAngDist', default=False, action='store_true', help='Reconstruct the angular distribution from the resolved resonance parameters (default: False)' )
    parser.add_argument('--noReconstruct',   dest='doResonanceReconstruction', default=True, action='store_false', help="Don't reconstruct resonances (default: False)'" )
    parser.add_argument('--showURRCloud',    default=False, action='store_true', help='Overlay a contour plot of the PDF for the cross section in the URR (default: False)' )
    parser.add_argument('--evaluationStyle', default='eval', type=str, help="Style in GNDS file to show")

    # Uncertainty related options
    parser.add_argument('--uncRatio',       default=False, action='store_true', help='Plot a ratio of the uncertainty over the data (default: False)' )
    parser.add_argument('--noUnc',          default=False, action='store_true', help='Do not plot uncertainties' )

    # Experimental data sources & controls for cross comparison
    parser.add_argument('--c4File',         default=None,  type=str, help="Optional C4 file to pull data from instead of (or in addition to) EXFOR" )
    parser.add_argument('--noX4',           default=False, action='store_true', help='Do not plot EXFOR data' )
    parser.add_argument('--showX4Evals',    default=True, action='store_false', help='Plot evaluations found in EXFOR library (with "V" SUBENT).  The default behavior is not to plot them.' )

    # Overall plot style controls.  Overrides defaults in DEFAULT_STYLE_DICT
    parser.add_argument('--style',          default=None, type=str, help='JSON file with plot style overrides, see "plot_defaults.json" in source distribution for examples')
    parser.add_argument('--figSizeX',       default=20.0, type=float, help="Width of generated figure (in cm), default=20 cm")
    parser.add_argument('--figSizeY',       default=10.0, type=float, help="Height of generated figure (in cm), default=10 cm")
    parser.add_argument('--logX',           dest='logX', default=None, action='store_true', help="Enable log scaling on x axis")
    parser.add_argument('--noLogX',         dest='logX', action='store_false', help="No log scaling on x axis")
    parser.add_argument('--logY',           dest='logY', default=None, action='store_true', help="Enable log scaling on y axis")
    parser.add_argument('--noLogY',         dest='logY', action='store_false', help="No log scaling on y axis")
    parser.add_argument('--useBokeh',       default=False, action='store_true', help="Use experimental hooks to bokeh for plotting (must specify outfile as well)")

    # Plot legend controls
    parser.add_argument('--noX4Legend',     default=False, action='store_true', help='Do not put legend for EXFOR data on plot' )
    parser.add_argument('--showX4Legend',   default=True,  action='store_true', help='Put legend for EXFOR data on plot, even if it would be automatically supressed' )
    parser.add_argument('--simpleLegend',   default=False, action='store_true', help='Do not put entries for the components of requested MT in the plot legend' )

    # Detailed control of the observable to plot
    parser.add_argument('--observable',  choices=plotstyles.DEFAULT_STYLE_DICT["plotStyles"], type=str, default="crossSection", help='The observable to plot (default is "crossSection")' )
    parser.add_argument('--showParts',   default=False, action='store_true', help='Show all of the channels that comprise the one requested (if applicable)' )
    parser.add_argument('--L',           type=int, default=0, help='Legendre moment to show (if applicable)' )
    parser.add_argument('--referenceFrame', choices=['centerOfMass','lab'], type=str, default=None, help='The frame of the observable to plot.' )

    return parser.parse_args()


#---------------------------------------------------
# main routine!
#---------------------------------------------------
if __name__ == "__main__":

    args = parseArgs()

    # Declare the various particles involved in the plot
    target = args.target
    projectile = args.projectile
    product = args.product

 #   print(product)
 #   exit()

    # Read in the plot style information
    userDefs = plotstyles.readUserStyles( args.style )
    if False:  # for debugging
        import json
        print( json.dumps(userDefs,indent=2) )
        exit()

    # Read the ENDF (& other user provided) data
    xyData = {}
    xydyData = {}
    xdxydyData = {}
    gndsMap = collections.OrderedDict()
    mtMap = {}
    print( fudge.core.utilities.brb.banner( "Reading evaluation files" ) )
    for endf in args.endf:
        if endf == "None": continue
        print( fudge.core.utilities.brb.winged_banner( "Reading "+endf ) )

        # Try special json file
        if endf.endswith( '.json' ):
            import json
            mixture=json.loads(open(endf).read())
            target=mixture['target']
            projectile=mixture['projectile']
            print( "   ",projectile,'+',target,'using mixture of evaluations:' )
            for iso in mixture['isotopes']:
                print( "       ", mixture['isotopes'][iso]['atomicFraction']*100,'% of', mixture['isotopes'][iso]['pathToFile'] )
            print()
            print('    Disabled plotting evaluation uncertainty')
            args.noUnc=True
            print()
            gndsMap['mixture']=mixture
            for iso in mixture['isotopes']:
                isoEndf=mixture['isotopes'][iso]['pathToFile']
                print( fudge.core.utilities.brb.winged_banner("Reading " + isoEndf) )
                gndsMap[isoEndf] = plotio.readEvaluation( isoEndf,
                                                         skipBadData=args.skipBadData,
                                                         reconstructResonances=args.doResonanceReconstruction,
                                                         continuumSpectraFix=args.continuumSpectraFix,
                                                         verbose=False,
                                                         skipCovariances=args.skipCovariances,
                                                         verboseWarnings=args.verboseWarnings,
                                                         printBadNK14=args.printBadNK14,
                                                         ignoreBadDate=args.ignoreBadDate,
                                                         acceptBadMF10FissionZAP=args.acceptBadMF10FissionZAP)

        # Try plain XY data
        elif endf.endswith( '.xy.dat'):       xyData[endf]   = plotio.readXYData( endf )

        # Try plain XYdY data
        elif endf.endswith( '.xydy.dat'):   xydyData[endf] = plotio.readXYdYData( endf )

        # Try plain XdXYdY data
        elif endf.endswith( '.xdxydy.dat'):   xdxydyData[endf] = plotio.readXdXYdYData( endf )

        # Must be an evaluation
        else: gndsMap[endf]   = plotio.readEvaluation( endf,
                                                      skipBadData=args.skipBadData,
                                                      reconstructResonances=args.doResonanceReconstruction,
                                                      continuumSpectraFix=args.continuumSpectraFix,
                                                      verbose=False,
                                                      skipCovariances=args.skipCovariances,
                                                      verboseWarnings=args.verboseWarnings,
                                                      printBadNK14=args.printBadNK14,
                                                      ignoreBadDate=args.ignoreBadDate,
                                                      acceptBadMF10FissionZAP=args.acceptBadMF10FissionZAP)

        if not endf in gndsMap: continue

        mtMap[endf] = getEvaluationMTs( gndsMap[endf][0] )

    # Get the target & projectile
    if target is None:     target = str(gndsMap.items()[0][1][0].target)
    if projectile is None: projectile = str(gndsMap.items()[0][1][0].projectile)

    # MT is specified by the user.  Just plot that.
    if args.mt != 0: mtList = [ args.mt ]

    # MT not specified, so assemble a list of mt's.
    # For discrete level excitations, we need to add the summed up version.
    # Also, put them in order so the first is always the sum.
    else:
        mtList = []
        if projectile=='n':
            # Neutrons... MT Map has all the MT's in every evaluation.  We just want highlights.
            for mt in mtMap[ args.endf[0] ]: 
                if   mt >= 50  and mt <= 91 :
                    if 4 not in mtList:   mtList.append( 4 )   # (*,n)
                elif mt >= 600 and mt <= 649:
                    if 103 not in mtList: mtList.append( 103 ) # (*,p)
                elif mt >= 650 and mt <= 699:
                    if 104 not in mtList: mtList.append( 104 ) # (*,d)
                elif mt >= 700 and mt <= 749:
                    if 105 not in mtList: mtList.append( 105 ) # (*,t)
                elif mt >= 750 and mt <= 799:
                    if 106 not in mtList: mtList.append( 106 ) # (*,3He)
                elif mt >= 800 and mt <= 849:
                    if 107 not in mtList: mtList.append( 107 ) # (*,a)
                elif mt >= 875 and mt <= 891:
                    if 16 not in mtList:  mtList.append( 16 )  # (*,2n)
                elif mt in [ 19, 20, 21, 38 ]:
                    if 18 not in mtList:  mtList.append( 18 )  # (*,f)
                elif mt not in mtList: mtList.append( mt )
                try: reaction = getEXFORRxn( mt )
                except KeyError:
                    if mt in range( 850, 871, 1 ): print( "Got lumped covariance for MT = "+ str( mt ) )
                    else: raise KeyError( "Unknown MT: "+str( mt ) )
        else: mtList = mtMap[args.endf[0]]
        mtList.sort()


    # Determine the names of the particles involved
    print( fudge.core.utilities.brb.banner( "Plot details" ) )
    if product is None and args.observable not in ['crossSection','crossSectionIntegrals',"energyDeposit",'energyBalance',"momentumDeposit",'momentumBalance','fissionEnergyRelease']:
        if args.observable in ["formFactor",'anomolousScatteringFactor']:
            print( "Assuming product is a 'gamma' for observable %s"%args.observable )
            product = 'gamma'
        elif args.observable in ["energyTransfer"]:
            print( "Assuming product is an 'e-' for observable %s"%args.observable )
            product = 'e-'
        elif args.observable in ["nubar"]:
            print( "Assuming product is an 'n' for observable %s"%args.observable )
            product = 'n'
        else:
            raise ValueError( "For observable='"+str(args.observable)+"', need to declare a product using --product command line argument" )
    print("Projectile is "+projectile)
    print("Target is "+target)
    print("List of MT's to plot:",mtList)


    # Determine the observable type
    if product is None: print("Observable is "+args.observable)
    else: print()

    # Loop over reactions
    for mt in mtList:

        try:

            # What is the proper name of the reaction?
            reaction = projectile.capitalize() +','+ getEXFORRxn( mt )

            # Compute the output filename
            if args.outFile is not None and not args.outFile.endswith('.png'):
                outFile = args.outFile+'_mt'+str(mt)+'_'+str(args.observable)+'.png'
            else: outFile = None

            # Determine the kind of plot and set up the axes
            print(fudge.core.utilities.brb.banner( "Generating plots for "+target+"("+reaction.lower()+')' ))
            if args.observable == "mubar":
                makeAngDistMubarPlot(
                    gndsMap, xyData, xydyData, xdxydyData,
                    mt=mt,
                    projectile=projectile, target=target, product=product,
                    referenceFrame=args.referenceFrame,
                    outFile=outFile,
                    plotStyle=userDefs,
                    figsize=(args.figSizeX,args.figSizeY),
                    useBokeh=args.useBokeh)

            elif args.observable == "formFactor":
                if mt not in [502,504]:
                    print("in ENDF, formFactor data only exists for (in)coherent photon-atom scattering (MT=502, 504), so only showing MT=502 and 504")
                makeFormFactorPlot(
                    gndsMap, xyData, xydyData, xdxydyData,
                    projectile=projectile, target=target, product=product,
                    outFile=outFile,
                    plotStyle=userDefs,
                    figsize=(args.figSizeX,args.figSizeY),
                    useBokeh=args.useBokeh )
                break

            elif args.observable == 'anomolousScatteringFactor':
                if mt != 502:
                    print("in ENDF, anomolousScatteringFactor data only exists for coherent photon-atom scattering (MT=502), so only showing MT=502")
                makeAnomolousScatteringFactorPlot(
                    gndsMap, xyData, xydyData, xdxydyData,
                    projectile=projectile, target=target, product=product,
                    outFile=outFile,
                    plotStyle=userDefs,
                    figsize=(args.figSizeX,args.figSizeY),
                    useBokeh=args.useBokeh )
                break

            elif args.observable == 'energyTransfer':
                if mt not in [527, 528]:
                    print("in ENDF, energyTransfer data only exists for bremstrahlung reactions (MT=527, 528)")
                    continue
                makeEnergyTransferPlot(
                    gndsMap, xyData, xydyData, xdxydyData,
                    mt=mt,
                    projectile=projectile, target=target, product=product,
                    outFile=outFile,
                    plotStyle=userDefs,
                    figsize=(args.figSizeX,args.figSizeY),
                    useBokeh=args.useBokeh )

            elif args.observable in ["energyDeposit",'energyBalance'] :
                makeEnergyDepositPlot(
                    gndsMap, xyData, xydyData, xdxydyData,
                    mt=mt,
                    projectile=projectile, target=target, product=product,
                    outFile=outFile,
                    plotStyle=userDefs,
                    observable =args.observable,
                    figsize=(args.figSizeX,args.figSizeY),
                    useBokeh=args.useBokeh )

            elif args.observable == "fissionEnergyRelease":
                makeFissionEnergyReleasePlot(
                    gndsMap, xyData, xydyData, xdxydyData,
                    mt=mt,
                    projectile=projectile, target=target, product=product,
                    outFile=outFile,
                    plotStyle=userDefs,
                    figsize=(args.figSizeX,args.figSizeY),
                    useBokeh=args.useBokeh )

            elif args.observable == "nubar":
                makeMultiplicityPlot(
                    gndsMap, xyData, xydyData, xdxydyData,
                    mt=mt,
                    projectile=projectile, target=target, product=product,
                    c4File=args.c4File,
                    showparts=args.showParts,
                    nox4evals=not args.showX4Evals,
                    nox4legend=args.noX4Legend,
                    evaluationStyle='eval',
                    logX=args.logX,
                    logY=args.logY,
                    figsize=(args.figSizeX, args.figSizeY),
                    outFile=outFile,
                    plotStyle=userDefs )

            elif args.observable == "crossSectionIntegrals":
                if mt==2 and projectile not in ['n','g']: 
                    print("Total elastic cross section doesn't make sense for charged particles, skipping plot\n")
                    continue
                makeCrossSectionIntegralsPlot(
                    gndsMap, xyData, xydyData, xdxydyData,
                    mt=mt,
                    projectile=projectile, target=target, product=product,
                    outFile=outFile,
                    plotStyle=userDefs,
                    figsize=(args.figSizeX,args.figSizeY),
                    useBokeh=args.useBokeh )

            elif args.observable == "crossSection":
                if mt==2 and projectile not in ['n','g']: 
                    print("Total elastic cross section doesn't make sense for charged particles, skipping plot\n")
                    continue
                if  args.uncRatio:
                    makeCrossSectionUncertaintyPlot(
                        gndsMap, xyData, xydyData, xdxydyData, c4File=args.c4File,
                        mt=0,
                        projectile='n', target='1H',
                        outFile=None,
                        plotStyle={},
                        figsize=(args.figSizeX,args.figSizeY),
                        useBokeh=args.useBokeh )

                else:
                    try:
                        makeCrossSectionPlot(
                            gndsMap, xyData, xydyData, xdxydyData,
                            mt=mt,
                            c4File=args.c4File,
                            projectile=projectile, target=target,
                            nounc=args.noUnc,
                            nox4=args.noX4,
                            showparts=args.showParts,
                            nox4evals=not args.showX4Evals,
                            nox4legend=args.noX4Legend,
                            evaluationStyle='eval',
                            logX=args.logX,
                            logY=args.logY,
                            outFile=outFile,
                            plotStyle=userDefs,
                            figsize=(args.figSizeX,args.figSizeY),
                            useBokeh=args.useBokeh )

                    except ValueError as err:
                        if 'absolute vs. relative' in err.message:
                            makeCrossSectionPlot(
                                gndsMap, xyData, xydyData, xdxydyData,
                                mt=mt,
                                projectile=projectile, target=target,
                                nounc=True,
                                nox4=args.noX4,
                                showparts=args.showParts,
                                nox4evals=not args.showX4Evals,
                                nox4legend=args.noX4Legend,
                                evaluationStyle='eval',
                                outFile=outFile,
                                plotStyle=userDefs,
                                figsize=(args.figSizeX,args.figSizeY),
                                useBokeh=args.useBokeh )
                            print("WARNING: problem in uncertainty calculation, encountered exception %s with message \"%s\", regenerating plot w/o uncertainty" % (str(type(err)).split('.')[-1].split("'")[0],err.message))
                        else:
                            traceback.print_exc(file=sys.stdout)
                            raise err

                    except AttributeError as err:
                        if "summedCovariance' object has no attribute 'axes" in err.message:
                            makeCrossSectionPlot(
                                gndsMap, xyData, xydyData, xdxydyData,
                                mt=mt,
                                projectile=projectile, target=target,
                                nounc=True,
                                nox4=args.noX4,
                                showparts=args.showParts,
                                nox4evals=not args.showX4Evals,
                                nox4legend=args.noX4Legend,
                                style='eval',
                                outFile=outFile,
                                plotStyle=userDefs,
                                figsize=(args.figSizeX,args.figSizeY),
                                useBokeh=args.useBokeh )
                            print("WARNING: problem in uncertainty calculation, encountered exception %s with message \"%s\", regenerating plot w/o uncertainty" % (str(type(err)).split('.')[-1].split("'")[0],err.message))
                        else:
                            if args.traceback: traceback.print_exc(file=sys.stdout)
                            raise err

            elif args.observable in ["momentumDeposit", 'momentumBalance']:
                makeMomentumDepositPlot(
                    gndsMap, xyData, xydyData, xdxydyData,
                    mt=mt,
                    projectile=projectile, target=target, product=product,
                    outFile=outFile,
                    plotStyle=userDefs,
                    observable =args.observable,
                    figsize=(args.figSizeX,args.figSizeY),
                    useBokeh=args.useBokeh )

            elif args.observable == "LegendreMoment":
                makeAngDistLegendreMomentPlot(
                    gndsMap, xyData, xydyData, xdxydyData,
                    mt=mt,
                    L=args.L,
                    projectile=projectile, target=target, product=product,
                    outFile=outFile,
                    plotStyle=userDefs,
                    figsize=(args.figSizeX,args.figSizeY),
                    useBokeh=args.useBokeh )
            else: raise ValueError( "Unsupported plot observable: %s" % args.observable )

        except NotImplementedError as err:
            if args.traceback: traceback.print_exc(file=sys.stdout)
            print("\nWARNING: plot not generated, a NotImplementedError was raised with message '%s'\n" % err.message)

#        except Exception, err:
#            if args.mt==0:
#                # Batch job, we must keep going at all costs!
#                print "\nWARNING: plot not generated, an %s was raised with message '%s'\n" % (str(type(err)).split('.')[-1].split("'")[0], err.message)
#            else:
#                # OK, regular run, raise an exception
#                raise err
