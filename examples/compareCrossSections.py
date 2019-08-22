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
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# <<END-copyright>>

"""
compareCrossSections.py: compare the cross section for given MT number from two different evaluated files.
"""
import os
import sys

def compare_plot( xsc1, xsc2, title="comparison plot", legend1="first file", legend2="second file",
        saveFile=None ):
    """ starting with XYs data for xsc1 and xsc2, draw a comparison plot """
    from fudge.vis.matplotlib import plot2d
    import matplotlib.pyplot as plt

    if xsc1.getDomain() != xsc2.getDomain():
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
            legendOn=True, thePlot = ax1, minY=0 )
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
    from optparse import OptionParser
    usage = """Usage: compareCrossSections.py MT evaluation1 evaluation2

The input files (evaluation1 and evaluation2) can be in GND or ENDF format.
The cross section for the given MT is extracted from each and drawn on one plot along with the relative difference.

If the given cross section has a resonance region contribution, for GND it will be reconstructed using Fudge.
For ENDF, however, reconstruction should be done before-hand using RECENT or another tool.
"""
    parser = OptionParser(usage)
    parser.add_option( "-t", dest="tolerance", type=float, default=0.001,
            help="specify tolerance for reconstruction. 0.001 => 0.1%" )
    parser.add_option( "-l", dest="legend", default=None, nargs=2,
            help="legends to attach to each cross section. Requires two legends" )
    parser.add_option( "-T", dest="title", default=None, 
            help="specify plot title" )
    return parser.parse_args()

if __name__ == '__main__':
    from fudge.gnd import reactionSuite
    from fudge.legacy.converting import endfFileToGND

    opts, args = process_args()
    mt = int(args[0])

    def getXS( filename, MT ):
        try:
            RS = reactionSuite.readXML( filename )
        except:
            try:
                rce = endfFileToGND.endfFileToGND( filename, singleMTOnly=MT, parseCrossSectionOnly=True, skipBadData=True )
                RS, c = rce['reactionSuite'], rce['covarianceSuite']
            except:
                print "File %s doesn't seem to be a legal ENDF or GND file!" % filename
                sys.exit()
        reac = [r for r in (RS.reactions + RS.summedReactions) if int( r.attributes['ENDF_MT'] )==MT]
        if len(reac) != 1:
            print "Couldn't find unique reaction for given MT in %s" % filename
        xsc = reac[0].crossSection
        if xsc.forms.keys() == ['resonancesWithBackground']:
            RS.reconstructResonances( opts.tolerance )
        try: pwxs = xsc.toPointwise_withLinearXYs( 1e-08 )
        except:
            pwxs = xsc.toPointwise_withLinearXYs( 1e-9 )
        return pwxs

    xs1, xs2 = getXS(args[1], mt), getXS(args[2], mt)

    if opts.legend: l1,l2 = opts.legend
    else: l1,l2 = args[1:3]
    if opts.title: title = opts.title
    else: title="MT%i xsc comparison" % mt

    compare_plot( xs1, xs2, title=title, legend1=l1, legend2=l2 )
