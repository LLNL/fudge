# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
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
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# <<END-copyright>>

"""
compareCrossSections.py: compare the cross section for given MT number from two different evaluated files.
"""
import os
import sys

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
    x1,y1 = map(numpy.array, diff.copyDataToXsAndYs())
    x2,y2 = map(numpy.array, mean.copyDataToXsAndYs())
    y2[ (y2==0)*(y1==0) ] = 1.0 # silence zero/zero division warnings
    relative_diff = zip(x1, y1/y2 * 100)
    
    """ # XYs division can take a long time, unnecessary in this case
    mean.setSafeDivide( True )  # control divide/0 errors
    relative_diff = (xsc1 - xsc2) / mean * 100
    """

    plot1 = plot2d.DataSet2d( xsc1, legend=legend1, symbol="+" )
    plot2 = plot2d.DataSet2d( xsc2, legend=legend2, lineStyle="--", symbol="+", color="red" )
    reldiff_plot = plot2d.DataSet2d( relative_diff, legend="percent difference" )

    xAxisSettings = plot2d.AxisSettings( label="", isLog=True )
    yAxisSettings = plot2d.AxisSettings( label="Cross Section (barn)", isLog=True )

    fig = plt.figure( figsize=(10,8) )
    fig.subplots_adjust( top=0.88, bottom=0.12, wspace=0.4 )

    ax1 = subplot2grid((4,1), (0,0), rowspan=3)
    mplot = plot2d.__makePlot2d( [plot1, plot2], xAxisSettings, yAxisSettings,
            legendOn=True, legendXY=legendXY, thePlot = ax1, minY=0 )
    plt.setp( ax1.get_xticklabels(), visible=False )
    plt.setp( ax1.get_label(), visible=False )

    # also plot the relative difference (needs different y-axis):
    xAxisSettings = plot2d.AxisSettings( label="$E_n$ (eV)", isLog=True )
    yAxisSettings = plot2d.AxisSettings( label="% diff" )

    ax2 = subplot2grid((4,1), (3,0), sharex=ax1)
    plot2d.__makePlot2d( [reldiff_plot], xAxisSettings, yAxisSettings,
            legendOn=False, thePlot = ax2, minY=0 )
    # tick marks may be too dense on this y-axis:
    #ax2.get_yaxis().set_ticks( [-0.2,0,0.2] )

    plt.suptitle( title, fontsize=24, fontweight='bold' )
    if saveFile: plt.savefig( saveFile )
    else: plt.show()

# Useful function, may not be available on older matplotlib installations:
def subplot2grid(shape, loc, rowspan=1, colspan=1, **kwargs):
    from matplotlib.pyplot import GridSpec, gcf
    fig = gcf()
    s1,s2 = shape
    subplotspec = GridSpec(s1, s2).new_subplotspec(loc, rowspan=rowspan, colspan=colspan)
    a = fig.add_subplot(subplotspec, **kwargs)
    bbox = a.bbox
    byebye = []
    for other in fig.axes:
        if other==a: continue
        if bbox.fully_overlaps(other.bbox):
            byebye.append(other)
    for ax in byebye: delaxes(ax)
    return a

def process_args():
    from argparse import ArgumentParser
    parser = ArgumentParser(
            description = """Compare the same cross section in two different evaluations,
            or compare a summed cross section from one evaluation with the sum of its parts
            (using option --summed).""",
            epilog = """Input files can be in GND or ENDF format.
            Resonances will be reconstructed ONLY for input files in GND format.
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
    parser.add_argument( "-S", "--summed", action='store_true', default=False,
            help="For a single evaluation, compare a summed cross section (e.g. total, inelastic) with the sum of its parts. Only one input file needed" )
    return parser.parse_args()

if __name__ == '__main__':
    from fudge.gnd import reactionSuite as reactionSuiteModule, sums as sumsModule
    from fudge.legacy.converting import endfFileToGND

    args = process_args()

    def getReactionSuite( filename, singleMTOnly = None ):
        try:
            RS = reactionSuiteModule.readXML( filename )
        except:
            try:
                rce = endfFileToGND.endfFileToGND( filename, singleMTOnly = singleMTOnly, parseCrossSectionOnly = True, skipBadData = True )
                RS, c = rce['reactionSuite'], rce['covarianceSuite']
            except:
                print "File %s doesn't seem to be a legal ENDF or GND file!" % filename
                sys.exit()
        RS.attributes['originalFile'] = filename
        return RS

    def getXS( reactionSuite, MT, sumsOnly = False ):
        allReacs = [reac for reac in reactionSuite.sums if isinstance(reac, sumsModule.crossSectionSum)]
        if not sumsOnly:
            allReacs += list(reactionSuite.reactions)
        reac = [r for r in allReacs if r.ENDF_MT == MT]
        if len(reac) != 1:
            print "Couldn't find unique reaction for MT%d in %s" % (MT,reactionSuite.attributes['originalFile'])
        xsc = reac[0].crossSection
        if xsc.evaluated.moniker == 'resonancesWithBackground':
            reactionSuite.reconstructResonances( styleName='reconstructed', accuracy=args.tolerance )
            pwxs = xsc['reconstructed']
        else:
            try:
                pwxs = xsc.toPointwise_withLinearXYs( 1e-08 )
            except:
                pwxs = xsc.toPointwise_withLinearXYs( 1e-10, 1e-10 )
        return pwxs

    if args.summed:
        RS = getReactionSuite( args.file1 )
        xs1 = getXS(RS, args.mt, sumsOnly = True)
        summedReac = [r for r in (RS.sums) if isinstance(r, sumsModule.crossSectionSum) and int( r.ENDF_MT ) == args.mt]
        if len(summedReac) != 1:
            print "Couldn't find unique summed reaction for MT%d in %s" % (args.mt,RS.attributes['originalFile'])
            sys.exit(1)
        summedReac = summedReac[0]
        if 'reconstructed' in summedReac.summands[0].link:
            summedXsc = summedReac.summands[0].link['reconstructed']
        else:
            summedXsc = summedReac.summands[0].link.toPointwise_withLinearXYs( 1e-08 )
        for summand in summedReac.summands[1:]:
            if 'reconstructed' in summand.link:
                newXsc = summand.link['reconstructed']
            else:
                newXsc = summand.link.toPointwise_withLinearXYs( 1e-08 )
            summedXsc, newXsc = summedXsc.mutualify( 1e-8,1e-8,0, newXsc, 1e-8,1e-8,0 )
            summedXsc += newXsc
        xs2 = summedXsc
        l1,l2 = ('tabulated sum','calculated sum')
    else:
        xs1 = getXS( getReactionSuite(args.file1, singleMTOnly=args.mt), args.mt )
        xs2 = getXS( getReactionSuite(args.file2, singleMTOnly=args.mt), args.mt )
        l1,l2 = args.file1, args.file2

    if args.legend: l1,l2 = args.legend
    if args.title: title = args.title
    else: title="MT%i xsc comparison" % args.mt

    legendXY = {'ul': (0.05, 0.95), 'ur': (0.75, 0.95),
            'll': (0.05, 0.2), 'lr': (0.75, 0.2)}.get( args.legendLocation )

    compare_plot( xs1, xs2, title=title, legend1=l1, legend2=l2, legendXY=legendXY )
