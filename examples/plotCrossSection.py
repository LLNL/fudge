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

""" basic plotting example: reconstructs resonances (if necessary), and plots cross sections for requested MTs.
 The original file can be either in GND or ENDF format. This example requires either matplotlib or gnuplot.

 Sample use:
>python plotCrossSection.py <filename> 2 18 102   #plot elastic, fission and capture """

import os
import sys
from fudge.core.utilities import argparse

from fudge.gnd import styles as stylesModule
from fudge.gnd.reactionData import crossSection as crossSectionModule

exampleDir = os.path.dirname( os.path.abspath( __file__ ) )
sys.path.insert(0, exampleDir)

parser = argparse.ArgumentParser(description = """
This example plots the cross section for the given MT number(s). The 'filename' argument can be either
a GND or ENDF-formatted file """ )
parser.add_argument('file', metavar='file', type=str, nargs=1, help='ENDF or GND file with cross section data')
parser.add_argument('mt', metavar='mt', type=int, nargs='+', help='MT number(s) to plot')
parser.add_argument('--temp', dest='temp', type=str, default=False, help="""Optional temperature for heating.
Temperature should be given as a string with units. Possible values are '1200 K' or '0.1 eV/k' (quotes are required)""")
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
    reac = [ r for r in RS if r.ENDF_MT == MT and hasattr(r,'crossSection') ]
    if not reac: print ("MT %i not present in the file" % MT); continue
    reac = reac[0]

    if args.temp:
        from fudge.gnd import physicalQuantity
        from pqu import PQU
        temp = PQU.PQU( args.temp )
        heatedStyle = stylesModule.heated( 'heated', derivedFrom=RS.styles.getEvaluatedStyle().label,
                temperature=physicalQuantity.temperature( temp.value, temp.unit ) )
        data[MT] = reac.crossSection.heat( heatedStyle, EMin=PQU.PQU(1e-5,'eV') )
    else:
        xsc = [form for form in reac.crossSection if hasattr(form, 'toPointwise_withLinearXYs')]
        if len(xsc)>0:
            data[MT] = xsc[0].toPointwise_withLinearXYs( accuracy = 1e-3, lowerEps = 1e-8 )
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

