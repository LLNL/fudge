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

""" basic plotting example: reconstructs resonances (if necessary), and plots cross sections for requested MTs.
 The original file can be either in GND or ENDF format. This example requires either matplotlib or gnuplot.

 Sample use:
>python plotCrossSection.py <filename> 2 18 102   #plot elastic, fission and capture """

import os
import sys
from fudge.core.utilities import argparse

exampleDir = os.path.dirname( os.path.abspath( __file__ ) )
sys.path.insert(0, exampleDir)

parser = argparse.ArgumentParser(description = """
This example plots the cross section for the given MT number(s). The 'filename' argument can be either
a GND or ENDF-formatted file """ )
parser.add_argument('file', metavar='file', type=str, nargs=1, help='ENDF or GND file with cross section data')
parser.add_argument('mt', metavar='mt', type=int, nargs='+', help='MT number(s) to plot')
parser.add_argument('--temp', dest='temp', type=str, default=False, help="""Optional temperature for heating.
Temperature should be given as a string with units. Possible values are '1200 K' or '0.1 eV' (quotes are required)""")
args = parser.parse_args()

filename = args.file[0]
if open(filename).readline().startswith( "<?xml" ):
    from fudge.gnd import reactionSuite
    RS = reactionSuite.readXML( filename )
else:
    from fudge.legacy.converting import endfFileToGND
    rce = endfFileToGND.endfFileToGND( filename )
    RS, CS = rce['reactionSuite'], rce['covarianceSuite']

reconstructed = False

data = {}
for MT in args.mt:
    reac = [ r for r in (RS.reactions + RS.summedReactions) if r.ENDF_MT == MT ]
    if not reac: print ("MT %i not present in the file" % MT); continue
    reac = reac[0]

    if reac.crossSection.nativeData == 'resonancesWithBackground' and not reconstructed:
        RS.reconstructResonances( styleName='reconstructed', accuracy=0.001 )
        reconstructed = True

    if args.temp:
        import pqu
        data[MT] = reac.crossSection.heat( pqu.PQU(args.temp), pqu.PQU('1e-5 eV') )
    elif 'pointwise' in reac.crossSection.forms:
        data[MT] = reac.crossSection['pointwise']
    else:
        xsc = [form for form in reac.crossSection.forms.values() if hasattr(form, 'toPointwise_withLinearXYs')]
        if len(xsc)>0:
            data[MT] = xsc[0].toPointwise_withLinearXYs(1e-8,0)
        else: print("Don't know how to plot cross section form(s): %s" % reac.crossSection.forms.keys())

    

# plotting:
lows,highs = zip(*[d.domain() for d in data.values()])
low,high = min(lows), max(highs)
title="%s MTs=%s" % (os.path.split(filename)[-1], ','.join([str(mt) for mt in args.mt]))
if args.temp:
    title += ' heated to %s' % args.temp
try:
    from fudge.vis.matplotlib import plot2d
    for key in data:
        data[key] = plot2d.DataSet2d( data[key], legend="MT%s" % key )
    plot2d.makePlot2d( data.values(),
            xAxisSettings = plot2d.AxisSettings(label="$E_n$ (eV)", isLog=True,
                axisMin = low, axisMax = high ),
            yAxisSettings = plot2d.AxisSettings(label="Cross Section (barn)", isLog=True, autoscale=True),
            title=title, legendOn=True )
except ImportError: # likely means matplotlib is missing. Try with Gnuplot instead:
    from fudge.vis.gnuplot import fudgeMultiPlots
    from fudge.legacy.endl.endl2dmathClasses import endl2dmath
    data = [ endl2dmath( d.copyDataToXYs() ) for d in data.values() ]
    fudgeMultiPlots.multiPlot( data, xylog=3, xMin=low, xMax=high, xLabel="E_n (eV)", title=title )

