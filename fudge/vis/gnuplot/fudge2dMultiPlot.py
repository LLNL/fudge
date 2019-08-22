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
This file creates an interactive gnuplot session that plots two column (2d) data from multiple files on 
the same window.  Each file must contain at least 2 columns of numbers. One column for the x-data and 
another for the y-data. For each file, the x- and y-data can be selected using the xColumn and yColumn 
keywords (see below) with the defaults being xColumn = 1 and yColumn = 2. This file is called as, 

python fudge2dMultiPlot.py [Parameters] files <datafile1> [keyword/value pairs] <datafile2> [keyword/value pairs] ...

The parameters can be any optional combination (with the except of "files") of the following:
Parameter   | value type | Comment
------------+------------+-------------------------------------------------------------------------------------------
   xMin     | float      | The argument after xMin is the initial minimum value of the x-axis.
   xMax     | float      | The argument after xMax is the initial maximum value of the x-axis.
   yMin     | float      | The argument after yMin is the initial minimum value of the y-axis.
   yMax     | float      | The argument after yMax is the initial maximum value of the y-axis.
   xLabel   | string     | The argument after xLabel is the x-axis label (e.g., xLabel "Energy (MeV)").
            |            | If value is a single word the "'s are not required.
   yLabel   | string     | The argument after yLabel is the y-axis label.
   title    | string     | The argument after title is the title for the plot.
   xylog    | integer    | If the argument after xylog is 0 x- and y-axis are linear, if 1 x-axis is log and y-axis is 
            |            | linear, if 2 x-axis is linear and y-axis is log and if 3 x- and y-axis are log.
   delete   | none       | If present the input files are deleted when this python script exits.
windowTitle | string     | The title of the main Tix windown.
   files    | string     | All arguments after the "files" parameter are filenames - and options keyword/value pairs - 
            |            | for each dataset to plot. "files" must be the last parameter.
       The following describes the keyword/value pairs that can follow each file name:
       Keyword | Comment
       --------+----------------------------------------------------------------------------------------------------
       title   | The label to use for this curve in the legend.
       xColumn | The column to use for the x axis. The next value must be an integer from 1 to n where n 
               | is the number of columns in the file.
       yColumn | The column to use for the y axis. See xColumn for value specification.


Examples::
    python fudge2dMultiPlot.py a
    python fudge2dMultiPlot.py xLabel "x-axis" yLabel "y-axis" xylog 2 files a b title "file b" c title "Hi" yColumn 3

In the last example the user has requested that the legend for file "b" be "file b" and for file "c" its is "Hi". Also,
column 1 is used as the x-data for all files. The y-data is column 2 for "a" and "b" and column 3 for file "c".
"""

import sys
import Tkinter
root = Tkinter.Tk()
import fudge2dMultiPlotClass

options = {}
for option in fudge2dMultiPlotClass.Options : options[option] = fudge2dMultiPlotClass.Options[option]
doDelete = 0
xlog = 0
ylog = 0

i = 1
n = len( sys.argv ) - 1
while ( i < n ) :
    inc = 2
    if( sys.argv[i] == "delete" ) :
        inc = 1
        doDelete = 1
    elif( sys.argv[i] == "files" ) :
        i += 1
        break
    elif( sys.argv[i] == 'windowTitle' ) :
        root.title( sys.argv[i+1] )
    elif( sys.argv[i] == 'xylog' ) :
        ylog = int( sys.argv[i+1] )
        options['xlog'] = ylog % 2
        options['ylog'] = ( ylog / 2 ) % 2
    elif( sys.argv[i] in options ) :
        options[sys.argv[i]] = sys.argv[i+1]
    else :
        print
        print "Invalid option = <%s>" % sys.argv[i]
        print
        for o in sys.argv : print o,
        sys.exit( )
    i = i + inc

f2dmultiPlot = fudge2dMultiPlotClass.fudge2dMultiPlotClass( root, options = options )
fileNameList = []
while ( i <= n ) :
    inputFileName = sys.argv[i]
    fileNameList.append( inputFileName )
    options = { 'title' : inputFileName, 'xColumn' : 1, 'yColumn' : 2 }
    i += 1
    while( i <= n ) :
        option = sys.argv[i]
        if( option in options ) :
            i += 1
            if( i > n ) : raise Exception( "\nError from fudge2dMultiPlot.py: %s option for file %s does not have an option value" % \
                    ( option, inputFileName ) )
            options[option] = sys.argv[i]
            i += 1
        else :
            break
    for option in [ 'xColumn', 'yColumn' ] :
        try :
            options[option] = int( options[option] )
        except :
            raise Exception( "\nError from fudge2dMultiPlot.py: could not convert %s option's value = %s to an integer" % ( option, options[option] ) )
    f2dmultiPlot.add2dFile( inputFileName, title = options['title'], xColumn = options['xColumn'], yColumn = options['yColumn'], replot = False )

f2dmultiPlot.replot( )
root.mainloop( )
if ( doDelete ) : 
    import os
    for inputFileName in fileNameList : os.remove( inputFileName )
