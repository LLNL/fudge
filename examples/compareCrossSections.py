#! /usr/bin/env python
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
compareCrossSections.py: compare the cross section for given MT number from two different evaluated files.
"""
import sys, traceback

def compare_plot( xsc1, xsc2, title="comparison plot", legend1="first file", legend2="second file",
        saveFile=None, legendXY = (0.05, 0.95) ):
    """ starting with XYs data for xsc1 and xsc2, draw a comparison plot """
    from fudge.vis.matplotlib import plot2d
    import matplotlib.pyplot as plt

    if xsc1.domain() != xsc2.domain():
        xsc1, xsc2 = xsc1.mutualify( 1e-8, 1e-8, 0, xsc2, 1e-8, 1e-8, 0 )
    diff = xsc1 - xsc2
    mean = (xsc1 + xsc2) / 2.0

    import numpy
    x1,y1 = list( map( numpy.array, diff.copyDataToXsAndYs() ) )
    x2,y2 = list( map( numpy.array, mean.copyDataToXsAndYs() ) )
    y2[ (y2==0)*(y1==0) ] = 1.0 # silence zero/zero division warnings
    relative_diff = list( zip( x1, y1 / y2 * 100 ) )

    """ # XYs division can take a long time, unnecessary in this case
    mean.nf_pointwiseXY.setSafeDivide( True )  # control divide/0 errors
    relative_diff = (xsc1 - xsc2) / mean * 100
    """

    plot1 = plot2d.DataSet2d( xsc1, legend=legend1, symbol="+" )
    plot2 = plot2d.DataSet2d( xsc2, legend=legend2, lineStyle="--", symbol="+", color="red" )
    reldiff_plot = plot2d.DataSet2d( relative_diff, legend="percent difference" )

    xAxisSettings = plot2d.AxisSettings( label="", isLog=True )
    yUnit = args.yUnit or 'barn'
    yAxisSettings = plot2d.AxisSettings( label="Cross Section (%s)" % yUnit, isLog=True )

    fig = plt.figure( figsize=(10,8) )
    fig.subplots_adjust( top=0.88, bottom=0.12, wspace=0.4 )

    ax1 = plt.subplot2grid((4,1), (0,0), rowspan=3)
    mplot = plot2d.__makePlot2d( [plot1, plot2], xAxisSettings, yAxisSettings,
            legendOn=True, legendXY=legendXY, thePlot = ax1, minY=0 )
    plt.setp( ax1.get_xticklabels(), visible=False )
    plt.setp( ax1.get_label(), visible=False )

    # also plot the relative difference (needs different y-axis):
    xUnit = args.xUnit or 'eV'
    xAxisSettings = plot2d.AxisSettings( label="$E_n$ (%s)" % xUnit, isLog=True )
    yAxisSettings = plot2d.AxisSettings( label="% diff" )

    ax2 = plt.subplot2grid((4,1), (3,0), sharex=ax1)
    plot2d.__makePlot2d( [reldiff_plot], xAxisSettings, yAxisSettings,
            legendOn=False, thePlot = ax2, minY=0)
    # tick marks may be too dense on this y-axis:
    #ax2.get_yaxis().set_ticks( [-0.2,0,0.2] )

    plt.suptitle( title, fontsize=24, fontweight='bold' )
    if saveFile: plt.savefig( saveFile )
    else: plt.show()

def process_args():
    from argparse import ArgumentParser
    parser = ArgumentParser(
            description = """Compare the same cross section in two different evaluations,
            or compare a summed cross section from one evaluation with the sum of its parts
            (using option --summed).""",
            epilog = """Input files can be in GNDS or ENDF format.
            Resonances will be reconstructed ONLY for input files in GNDS format.
            For ENDF-6 files, reconstruction should be done before-hand using RECENT or another tool.""",
            )
    parser.add_argument( "mt", type=int, help="ENDF MT of reaction to compare" )
    parser.add_argument( "file1", type=str, help="First file" )
    parser.add_argument( "file2", type=str, nargs="?",
            help="Second file (required unless option --summed supplied)" )
    parser.add_argument( "-t", "--tolerance", type=float, default=0.001,
            help="specify tolerance for reconstruction. 0.001 => 0.1%%" )
    parser.add_argument( "-l", "--legend", default=None, nargs=2,
            help="legend for each evaluation. Requires two legends" )
    parser.add_argument( "-L", "--legendLocation", default='ul',
            help="legend location: 'ul', 'ur', 'll' or 'lr'" )
    parser.add_argument( "-T", "--title", default=None,
            help="specify plot title" )
    parser.add_argument( "-o", "--outfile", default=None, help="Output file name")
    parser.add_argument( "--xUnit", type=str, help="Convert x-axes to this unit (e.g. MeV)" )
    parser.add_argument( "--yUnit", type=str, help="Convert y-axes to this unit (e.g. mb)" )
    parser.add_argument( "-S", "--summed", action='store_true', default=False,
            help="For a single evaluation, compare a summed cross section (e.g. total, inelastic) with the sum of its parts. Only one input file needed" )
    return parser.parse_args()

if __name__ == '__main__':
    from fudge import reactionSuite as reactionSuiteModule, styles as stylesModule
    from fudge.reactionData import crossSection
    from brownies.legacy.converting import endfFileToGNDS

    args = process_args()

    reconstructedStyleName = 'tmp_reconstructed'

    def getReactionSuite( filename, singleMTOnly = None ):
        try:
            RS = reactionSuiteModule.ReactionSuite.readXML_file( filename )
        except:
            try:
                rce = endfFileToGNDS.endfFileToGNDS( filename, singleMTOnly = singleMTOnly,
                                                     skipBadData = True, continuumSpectraFix = True,
                                                     parseCrossSectionOnly = True)
                RS, c = rce['reactionSuite'], rce['covarianceSuite']
            except Exception as excep:
                print("Exception raised:", excep)
                print("File %s doesn't seem to be a legal ENDF or GNDS file!" % filename)
                traceback.print_exc(file=sys.stdout)
                sys.exit()
        RS.originalFile = filename
        return RS

    def getXS( reactionSuite, MT, sumsOnly = False ):
        allReacs = list( reactionSuite.sums.crossSectionSums )
        if not sumsOnly:
            allReacs += list(reactionSuite.reactions)
        reac = [r for r in allReacs if r.ENDF_MT == MT]
        if len(reac) != 1:
            print("Couldn't find unique reaction for MT%d in %s" % (MT, reactionSuite.originalFile))
        xsc = reac[0].crossSection
        if isinstance( xsc.evaluated, crossSection.ResonancesWithBackground ):
            evalStyle = reactionSuite.styles.getEvaluatedStyle()
            reconstructedStyle = stylesModule.CrossSectionReconstructed( reconstructedStyleName, derivedFrom=evalStyle.label )
            reactionSuite.reconstructResonances( reconstructedStyle, accuracy=args.tolerance )
            pwxs = xsc[ reconstructedStyleName ]
        else:
            pwxs = xsc.toPointwise_withLinearXYs( accuracy = 1e-3, lowerEps = 1e-8 )
        return pwxs.convertAxisToUnit(1,'eV').convertAxisToUnit(0,'b')

    if args.summed:
        RS = getReactionSuite( args.file1 )
        xs1 = getXS(RS, args.mt, sumsOnly = True)
        summedReac = [r for r in (RS.sums.crossSectionSums) if int( r.ENDF_MT ) == args.mt]
        if len(summedReac) != 1:
            print("Couldn't find unique summed reaction for MT%d in %s" % (args.mt, RS.originalFile))
            sys.exit(1)
        summedReac = summedReac[0]
        if reconstructedStyleName in summedReac.summands[0].link:
            summedXsc = summedReac.summands[0].link[ reconstructedStyleName ]
        else:
            summedXsc = summedReac.summands[0].link.toPointwise_withLinearXYs( accuracy = 1e-3, lowerEps = 1e-8 )
        for summand in summedReac.summands[1:]:
            if reconstructedStyleName in summand.link:
                newXsc = summand.link[ reconstructedStyleName ]
            else:
                newXsc = summand.link.toPointwise_withLinearXYs( accuracy = 1e-3, lowerEps = 1e-8 )
            summedXsc, newXsc = summedXsc.mutualify( 1e-8,1e-8,0, newXsc, 1e-8,1e-8,0 )
            summedXsc += newXsc
        xs2 = summedXsc.convertAxisToUnit(1,'eV').convertAxisToUnit(0,'b')
        l1,l2 = ('tabulated sum','calculated sum')
    else:
        rs1 = getReactionSuite(args.file1, singleMTOnly=args.mt)
        xs1 = getXS( rs1, args.mt )
        rs2 = getReactionSuite(args.file2, singleMTOnly=args.mt)
        xs2 = getXS( rs2, args.mt )
        l1,l2 = args.file1, args.file2

    if args.xUnit:
        for xs in (xs1,xs2):
            xs.convertUnits( {xs.axes[1].unit: args.xUnit } )
    if args.yUnit:
        for xs in (xs1,xs2):
            xs.convertUnits( {xs.axes[0].unit: args.yUnit } )

    if args.legend: l1,l2 = args.legend
    if args.title: title = args.title
    else: title="MT%i xsc comparison" % args.mt

    legendXY = {'ul': (0.05, 0.95), 'ur': (0.75, 0.95),
            'll': (0.05, 0.2), 'lr': (0.75, 0.2)}.get( args.legendLocation )

    compare_plot( xs1, xs2, title=title, legend1=l1, legend2=l2, legendXY=legendXY, saveFile=args.outfile )
