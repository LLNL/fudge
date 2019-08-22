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

"""
This file creates an interactive gnuplot session that plots 3d data from multiple files on the same window. 
Each file must contain at least 3 columns of numbers. One column is for the x-data, one for the y-data, one
for the z-data and all others are ignored. For each file, the x-, y- and z-data can be selected using the 
xColumn, yColumn and zColumn keywords (see below) with the defaults being xColumn = 1, yColumn = 2 and
zColumn = 3. This file is called as,

python fudge3dMultiPlot.py [Parameters] files <datafile1> [keyword/value pairs] <datafile2> [keyword/value pairs] ...

The parameters can be any optional combination (with the except of "files") of the following:
Parameters | value type | Comment
-----------+------------+-------------------------------------------------------------------------------------------
   xMin    | float      | The argument after xMin is the initial minimum value of the x-axis.
   xMax    | float      | The argument after xMax is the initial maximum value of the x-axis.
   yMin    | float      | The argument after yMin is the initial minimum value of the y-axis.
   yMax    | float      | The argument after yMax is the initial maximum value of the y-axis.
   zMin    | float      | The argument after zMin is the initial minimum value of the z-axis.
   zMax    | float      | The argument after zMax is the initial maximum value of the z-axis.
   xLabel  | string     | The argument after xLabel is the x-axis label (e.g., xLabel "Energy (MeV)").
           |            | If value is a single word the "'s are not required.
   yLabel  | string     | The argument after yLabel is the y-axis label.
   zLabel  | string     | The argument after zLabel is the z-axis label.
   title   | string     | The argument after title is the title for the plot.
   xyzlog  | integer    | If the argument after xyzlog is 0 x-, y- and z-axis are linear, if xyzlog % 2 is 1 
           |            | x-axis is log otherwise linear, if ( xyzlog / 2 ) % 2 is 1 y-axis is log otherwise 
           |            | linear, if ( xyzlog / 4 ) % 2 is 1 z-axis is log otherwise linear.
   delete  | none       | If present the input files are deleted when this python script exits
   files   | string     | All arguments after the "files" parameter are filenames - and options keyword/value pairs -
           |            | for each dataset to plot. "files" must be the last parameter.
       The following describes the keyword/value pairs that can follow each file name:
       Keyword | Comment
       --------+----------------------------------------------------------------------------------------------------
       title   | The label to use for this curve in the legend
       xColumn | The column to use for the x axis. The next value must be an integer from 1 to n where n is the number of columns in the file.
       yColumn | The column to use for the y axis. See xColumn for value specification.
       zColumn | The column to use for the z axis. See xColumn for value specification.


Examples::
   python fudge3dMultiPlot.py a
   python fudge3dMultiPlot.py xLabel "x-axis" yLabel "y-axis" xyzlog 2 files a b title "file b" zColumn 2 yColumn 3

In the last example the user has requested that the legend for file "b" be "file b". Also,
the x-, y- and z-data are (1,2,3) for file "a" and (1,3,2) for file "b".
"""

import os, sys, shutil
from Tkinter import *
import tkFileDialog
import Gnuplot
from fudge.core.utilities import subprocessing

__metaclass__ = type

gnuplot_old = 0
if __name__ == '__main__' :
    returnCode, stdout, stderr = subprocessing.executeCommand( [ "gnuplot", "--version" ] )
    gnuplot_old = returnCode != 0

Options = { "xMin" : "*", "xMax" : "*", "yMin" : "*", "yMax" : "*", "zMin" : "*", "zMax" : "*", "xLabel" : "x", "yLabel" : "y", "zLabel" : "z", \
    "xrot" : 60, "zrot" : 30, "title" : "", "xyzlog" : 0 }

DoDelete = 0
i = 1
n = len( sys.argv ) - 1
while ( i < n ) :
    inc = 2
    if( sys.argv[i] == "delete" ) :
        inc = 1
        DoDelete = 1
    elif( sys.argv[i] == "files" ) :
        i += 1
        break
    elif( sys.argv[i] in Options ) :
        Options[sys.argv[i]] = sys.argv[i+1]
    else :
        print
        print "Invalid option = <%s>" % sys.argv[i]
        print
        for o in sys.argv : print o,
        sys.exit( )
    i = i + inc

xrot = float( Options["xrot"] )
zrot = float( Options["zrot"] )
zlog = int( Options["xyzlog"] )
xlog = zlog % 2
zlog /= 2
ylog = zlog % 2
zlog /= 2

fileNameArgs = ()
fileNameList = []
tmpFileNames = []
while ( i <= n ) :
    InputFileName = sys.argv[i]
    fileNameList.append( InputFileName )
    title = InputFileName
    options = { 'title' : InputFileName, 'xColumn' : 1, 'yColumn' : 2, 'zColumn' : 3 }
    i += 1
    while( i <= n ) :
        option = sys.argv[i]
        if( option in options ) :
            i += 1
            if( i > n ) : raise Exception( "\nError from fudge3dMultiPlot.py: %s option for file %s does not have an option value" % \
                                ( option, InputFileName ) )
            options[option] = sys.argv[i]
            i += 1
        else :
            break
    for option in [ 'xColumn', 'yColumn', 'zColumn' ] :
        try :
            options[option] = int( options[option] )
        except :
            raise Exception( "\nError from fudge3dMultiPlot.py: could not convert %s option's value = %s to an integer" % ( option, options[option] ) )
    f = Gnuplot.File( InputFileName, title = options['title'] )
    fileNameArgs = fileNameArgs + ( f, )

g = Gnuplot.Gnuplot( )
g( 'set ticslevel 0.' )

def Replot( extraGnuplotCommand = None ) :
    """This function is called whenever the plot needs to be redrawn (e.g., axis is changed)."""

    if( len( fileNameArgs ) == 0 ) : return
    g( "set nokey" )
    if( legendOnVar.get( ) ) : g( "set key" )

    if( not gnuplot_old ) :
        fs_ = fontSizeVar.get( )
        fs = max( 8, min( fs_, 72 ) )
        g( 'set ylabel "%s" font "times,%d"' % ( yAxisVar.get( ), fs ) )
        g( 'set terminal X11 font "times,%d"' % fs )

    ps_ = pointSizeVar.get( )
    ps = max( 1, min( ps_, 20 ) )
    if( ps != ps_ ) : pointSizeVar.set( ps )

    lw_ = lineWidthVar.get( )
    lw = max( 1, min( lw_, 20 ) )
    if( lw != lw_ ) : lineWidthVar.set( lw )

    lp = lineTypeVar.get( )
    if( lp == 'lines' ) :
        w = 'lines linewidth %d' % lw
    elif( lp == 'points' ) :
        w = 'points pointsize %d' % ps
    else :
        w = 'linespoints linewidth %d pointsize %d' % ( lw, ps )
    for d in fileNameArgs :
        d._options['with'] = ( w, 'with ' + w )
#        d.set_option( with=w )

    g( 'set xlabel "%s"' % xAxisVar.get( ) )
    g( 'set ylabel "%s"' % yAxisVar.get( ) )
    g( 'set zlabel "%s"' % zAxisVar.get( ) )
    g( 'set title "%s"' % TitleVar.get( ) )
    g( "set nologscale xyz" )
    xlog = xlogVar.get( )
    ylog = ylogVar.get( )
    zlog = zlogVar.get( )
    if ( xlog or ylog or zlog ) :
        s = "set logscale "
        if ( xlog ) : s = s + "x"
        if ( ylog ) : s = s + "y"
        if ( zlog ) : s = s + "z"
        g( s )

    g( "set view %e, %e" % ( xrot, zrot ) )
    if( extraGnuplotCommand is not None ) : g( extraGnuplotCommand )
    if( len( fileNameArgs ) ) : apply( g.splot, fileNameArgs )

SaveAsCounter = 0

def PlotSaveAsEps( ) :
    """Called when "File -> SaveAs eps" menu is selected. Uses a FileDialog to get the output 
file name and produces an eps file of the current plot."""

    global SaveAsCounter
    fn = "Fudge%d.eps" % SaveAsCounter
    FileName = tkFileDialog.asksaveasfilename( initialfile = fn )
    if( ( type( FileName ) == type( "" ) ) and ( FileName != "" ) ) :
        g.hardcopy( filename = FileName, color = 1 )
        SaveAsCounter += 1

def PlotSaveAsASCII( ) :
    """Currently not implemented."""

    global SaveAsCounter
    fn = "Fudge%d.asc" % SaveAsCounter
    FileName = tkFileDialog.asksaveasfilename( initialfile = fn )
    if( ( type( FileName ) == type( "" ) ) and ( FileName != "" ) ) :
        shutil.copyfile( InputFileName, FileName )
        SaveAsCounter += 1

def PlotPrint( ) :
    """Called when "File -> Print" menu is selected. Prints the current plot."""

    g.hardcopy( )

class logButton :
    """A class that contains a Checkbutton widget that toggles the x-, y- or z-axis between linear and log scaling."""

    def __init__( self, master, label, row, column, var ) :
        """The contructor for the logButton class."""

        self.label = label
        self.button = Checkbutton( frame, text=label, command = self.logcallback, variable = var )
        self.button.grid( row = row, column = column )

    def logcallback( self ) :
        """Called whenever the x-, y- or z-axis log checkbutton state changes."""

        Replot( )

class xzRot :
    """A class that contains a Scale widget for rotating the plot about an axis."""

    def __init__( self, master, label, to_, value, row, column, span ) :
        """The contructor for the xzRot class."""

        self.label = label[0]
        self.scale = Scale( master, label = label, length = to_, from_ = 0, to = to_, orient = "horizontal", command = self.xzcallback )
        self.scale.set( value )
        self.scale.grid( row = row, column = column, columnspan = span )

    def xzcallback( self, xz ) :
        """Called whenever the Scale widget is changed."""

        global xrot, zrot
        if ( self.label == "x" ) :
            xrot = int( xz )
        else :
            zrot = int( xz )
        Replot( )

root = Tk( )
if( Options["title"] != "" ) : root.title( Options["title"] )

menubar = Menu( root )
filemenu = Menu( menubar, tearoff = 0 )
filemenu.add_command( label = "SaveAs eps", command = PlotSaveAsEps )
filemenu.add_command( label = "Print", command = PlotPrint )
filemenu.add_separator( )
filemenu.add_command( label = "Exit", command = root.quit )
menubar.add_cascade( label = "File", menu = filemenu )
optionsmenu = Menu( menubar, tearoff = 0 )

def legendCallback( ) :
    """Called whenever the "Options -> Legend" menu state is changed."""

    Replot( )

legendOnVar = IntVar( )
optionsmenu.add_checkbutton( label = "Legend", variable = legendOnVar, command = legendCallback )
menubar.add_cascade( label = "Options", menu = optionsmenu )
root.config(menu=menubar)

frame = LabelFrame( root, text = "Plot Options" )
frame.pack( )

def Entry_Callback( e ) :
    """Called whenever the <Return> key is pressed in certain entry window."""

    Replot( )

Row = 0
lineTypes = [ ( "Lines", "lines" ), ( "Points", "points" ), ( "Lines & Points", "linespoints" ) ]
lineTypeVar = StringVar( )
lineTypeVar.set( lineTypes[0][1] )
c = 0
for t, m in lineTypes :
    b = Radiobutton( frame, text=t, variable = lineTypeVar, value = m, command = Replot )
    b.grid( row = Row, column = c )
    c += 1

Row += 1
lineWidthLabel = Label( frame, text = "point, line width = " )
lineWidthLabel.grid( row = Row, column = 0, sticky = E )

pointSizeVar = IntVar( )
pointSizeVar.set( 3 )
pointSizeEntry = Entry( frame, textvariable = pointSizeVar )
pointSizeEntry.grid( row = Row, column = 1, sticky = W + E )
pointSizeEntry.bind( "<Return>", Entry_Callback )

lineWidthVar = IntVar( )
lineWidthVar.set( 1 )
lineWidthEntry = Entry( frame, textvariable = lineWidthVar )
lineWidthEntry.grid( row = Row, column = 2, sticky = W + E )
lineWidthEntry.bind( "<Return>", Entry_Callback )

if( not gnuplot_old ) :
    Row += 1
    fontWidthLabel = Label( frame, text = "font size = " )
    fontWidthLabel.grid( row = Row, column = 0, sticky = E )

fontSizeVar = IntVar( )
fontSizeVar.set( 12 )
fontSizeEntry = Entry( frame, textvariable = fontSizeVar )
fontSizeEntry.grid( row = Row, column = 1, sticky = W + E )
fontSizeEntry.bind( "<Return>", Entry_Callback )

Row += 1
xAxisLabel = Label( frame, text = "x label = " )
xAxisLabel.grid( row = Row, column = 0, sticky = E )
xAxisVar = StringVar( )
xAxisVar.set( Options["xLabel"] )
xAxisEntry = Entry( frame, textvariable = xAxisVar )
xAxisEntry.grid( row = Row, column = 1, columnspan = 3, sticky = W + E )
xAxisEntry.bind( "<Return>", Entry_Callback )

Row += 1
yAxisLabel = Label( frame, text = "y label = " )
yAxisLabel.grid( row = Row, column = 0, sticky = E )
yAxisVar = StringVar( )
yAxisVar.set( Options["yLabel"] )
yAxisEntry = Entry( frame, textvariable = yAxisVar )
yAxisEntry.grid( row = Row, column = 1, columnspan = 3, sticky = W + E )
yAxisEntry.bind( "<Return>", Entry_Callback )

Row += 1
zAxisLabel = Label( frame, text = "z label = " )
zAxisLabel.grid( row = Row, column = 0, sticky = E )
zAxisVar = StringVar( )
zAxisVar.set( Options["zLabel"] )
zAxisEntry = Entry( frame, textvariable = zAxisVar )
zAxisEntry.grid( row = Row, column = 1, columnspan = 3, sticky = W + E )
zAxisEntry.bind( "<Return>", Entry_Callback )

Row += 1
TitleLabel = Label( frame, text = "title = " )
TitleLabel.grid( row = Row, column = 0, sticky = E )
TitleVar = StringVar( )
TitleVar.set( Options["title"] )
TitleEntry = Entry( frame, textvariable = TitleVar )
TitleEntry.grid( row = Row, column = 1, columnspan = 3, sticky = W + E )
TitleEntry.bind( "<Return>", Entry_Callback )

Row += 1
s_xrot = xzRot( frame, "x rot", 180, xrot, Row, 0, 2 )

Row += 1
s_zrot = xzRot( frame, "z rot", 360, zrot, Row, 0, 3 )

Row += 1
xlabel = Label( frame, text = "xrange" )
xlabel.grid( row = Row, column = 0, columnspan = 1, sticky = E, pady = 4 )

def xrange_Callback( e ) :
    """Called when the <Return> key is pressed in the xrange entry window."""

    s = "set xrange [ %s ]" % xrange.get( )
    Replot( s )

xrange = Entry( frame )
RangeStr = "%s : %s" % ( Options["xMin"], Options["xMax"] )
g( "set xrange [ %s ]" % RangeStr )
xrange.insert( 0, RangeStr )
xrange.grid( row = Row, column = 1, columnspan = 2, sticky = W + E )
xrange.bind( "<Return>", xrange_Callback )
xlogVar = IntVar( )
xlogVar.set( xlog )
xLogButton = logButton( frame, "xlog", Row, 3, xlogVar )

Row += 1
ylabel = Label( frame, text = "yrange" )
ylabel.grid( row = Row, column = 0, columnspan = 1, sticky = E, pady = 4 )
def yrange_Callback( e ) :
    """Called when the <Return> key is pressed in the yrange entry window."""

    s = "set yrange [ %s ]" % yrange.get( )
    Replot( s )

yrange = Entry( frame )
RangeStr = "%s : %s" % ( Options["yMin"], Options["yMax"] )
g( "set yrange [ %s ]" % RangeStr )
yrange.insert( 0, RangeStr )
yrange.grid( row = Row, column = 1, columnspan = 2, sticky = W + E )
yrange.bind( "<Return>", yrange_Callback )
ylogVar = IntVar( )
ylogVar.set( ylog )
yLogButton = logButton( frame, "ylog", Row, 3, ylogVar )

Row += 1
zlabel = Label( frame, text = "zrange" )
zlabel.grid( row = Row, column = 0, columnspan = 1, sticky = E, pady = 4 )
def zrange_Callback( e ) :
    """Called when the <Return> key is pressed in the zrange entry window."""

    s = "set zrange [ %s ]" % zrange.get( )
    Replot( s )

zrange = Entry( frame )
RangeStr = "%s : %s" % ( Options["zMin"], Options["zMax"] )
g( "set zrange [ %s ]" % RangeStr )
zrange.insert( 0, RangeStr )
zrange.grid( row = Row, column = 1, columnspan = 2, sticky = W + E )
zrange.bind( "<Return>", zrange_Callback )
zlogVar = IntVar( )
zlogVar.set( zlog )
zLogButton = logButton( frame, "zlog", Row, 3, zlogVar )

Row += 1
gnuplotCommandLabel = Label( frame, text = "plot command" )
gnuplotCommandLabel.grid( row = Row, column = 0, sticky = W, pady = 4 )

def gnuplotCommand_Callback( e ) :
    """Called when the <Return> key is pressed in the gnuplotCommand entry window."""

    Replot( gnuplotCommand.get( ) )

gnuplotCommand = Entry( frame )
gnuplotCommand.grid( row = Row, column = 1, columnspan = 3, sticky = W + E )
gnuplotCommand.bind( "<Return>", gnuplotCommand_Callback )

if __name__ == '__main__' :
    Replot( )
    root.mainloop( )
    if ( DoDelete ) : 
        for InputFileName in fileNameList : os.remove( InputFileName )
    for tmpFileName in tmpFileNames : os.remove( tmpFileName )
