# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

""" basic plotting example: reconstructs resonances (if necessary), and plots cross sections for requested MTs.
 The original file can be either in GNDS or ENDF format. This example requires either matplotlib or gnuplot.

 Sample use:
>python plotCrossSection.py <filename> 2 18 102   #plot elastic, fission and capture """

import os
import sys
import argparse

from fudge import physicalQuantity as physicalQuantityModule
from fudge import styles as stylesModule

exampleDir = os.path.dirname( os.path.abspath( __file__ ) )
sys.path.insert(0, exampleDir)

parser = argparse.ArgumentParser(description = """
This example plots the cross section for the given MT number(s). The 'filename' argument can be either
a GNDS or ENDF-formatted file """ )
parser.add_argument('file', metavar='file', type=str, nargs=1, help='ENDF or GNDS file with cross section data')
parser.add_argument('mt', metavar='mt', type=int, nargs='+', help='MT number(s) to plot')
parser.add_argument('--title', help="Plot title")
parser.add_argument('--temp', dest='temp', type=str, default=False, help="""Optional temperature for heating.
Temperature should be given as a string with units. Possible values are '1200 K' or '0.1 eV/k' (quotes are required)""")
parser.add_argument('--legendLabel', action='store_true', help="Use reaction label (instead of MT#) in legend")
args = parser.parse_args()

filename = args.file[0]
if open(filename).readline().startswith( "<?xml" ):
    from fudge import reactionSuite
    RS = reactionSuite.ReactionSuite.readXML_file( filename )
else:
    from brownies.legacy.converting import endfFileToGNDS
    rce = endfFileToGNDS.endfFileToGNDS(filename)
    RS, CS = rce['reactionSuite'], rce['covarianceSuite']

data = {}
for MT in args.mt:
    reac = [ r for r in RS if r.ENDF_MT == MT and hasattr(r,'crossSection') ]
    if not reac:
        print( "MT %i not present in the file" % MT )
        continue
    reac = reac[0]

    if args.temp:
        from pqu import PQU
        temp = PQU.PQU( args.temp )
        heatedStyle = stylesModule.Heated( 'heated', derivedFrom=RS.styles.getEvaluatedStyle().label,
                temperature=physicalQuantityModule.Temperature( temp.value, temp.unit ) )
        data[(MT, reac.label)] = reac.crossSection.heat( heatedStyle, EMin=PQU.PQU(1e-5,'eV') )
    else:
        xsc = [form for form in reac.crossSection if hasattr(form, 'toPointwise_withLinearXYs')]
        if len(xsc)>0:
            data[(MT, reac.label)] = xsc[0].toPointwise_withLinearXYs( accuracy = 1e-3, lowerEps = 1e-8 )
        else:
            print( "Don't know how to plot cross section form(s): %s" % list( reac.crossSection.forms.keys() ) )


# plotting:
lows, highs = list( zip( *[ d.domain() for d in data.values() ] ) )
low,high = min(lows), max(highs)
title="%s MTs=%s" % (os.path.split(filename)[-1], ','.join([str(mt) for mt in args.mt]))
if args.temp:
    title += ' heated to %s' % args.temp
if args.title:
    title = args.title
try:
    from fudge.vis.matplotlib import plot2d
    for key in sorted(data.keys()):
        if args.legendLabel:
            legend = key[1]
        else:
            legend = 'MT%d' % key[0]
        data[key] = plot2d.DataSet2d( data[key], legend=legend )
    plot2d.makePlot2d( list( data.values() ),
            xAxisSettings = plot2d.AxisSettings(label="$E_{%s}$ (eV)" % RS.projectile, isLog=True,
                axisMin = low, axisMax = high ),
            yAxisSettings = plot2d.AxisSettings(label="Cross Section (barn)", isLog=True, autoscale=True),
            title=title, legendOn=True )
except ImportError: # likely means matplotlib is missing. Try with Gnuplot instead:
    from fudge.vis.gnuplot import fudgeMultiPlots
    from brownies.legacy.endl import endl2dmath
    data = [ endl2dmath( d.copyDataToXYs() ) for d in data.values() ]
    fudgeMultiPlots.multiPlot( data, xylog=3, xMin=low, xMax=high, xLabel="E_n (eV)", title=title )
