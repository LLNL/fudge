# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
"""
This file creates an interactive gnuplot session that plots 2d data from a single file.
The data file must contain 2 columns (called 2d data) of numbers. The first column is 
plotted as the x data and the second column is plotted as the y data.
This file is called as, 

python endl2dplot.py [xLabel string] [yLabel string] [title string] [xylog int] [delete] <datafile>

where the parameters enclosed by [] are optional.

Parameters::
   xMin     The next argument is the initial minimum value of the x-axis
   xMax     The next argument is the initial maximum value of the x-axis
   yMin     The next argument is the initial minimum value of the y-axis
   yMax     The next argument is the initial maximum value of the y-axis
   xLabel   The next argument is the x-axis label enclosed in "'s (e.g., xLabel "Energy (MeV)")
   yLabel   The next argument is the y-axis label enclosed in "'s
   title    The next argument is the title enclosed in "'s
   xylog    If the next argument is 0 x- and y-axis are linear, if 1 x-axis is log and y-axis is linear,
                if 2 x-axis is linear and y-axis is log and if 3 x- and y-axis are log
   delete   If present the input file is delete with this python script exits

Examples::
    python endl2dplot.py a
    python endl2dplot.py xLabel "x-axis" yLabel "y-axis" xylog 2 a
"""

import sys, os, string, shutil
import math

if( sys.version_info[0] == 2 ) :
    from Tkinter import *
    import tkFileDialog
else :
    from tkinter import *
    import tkinter.filedialog as tkFileDialog

import Gnuplot
from LUPY import subprocessing

__metaclass__ = type

gnuplot_old = 0
if __name__ == '__main__' :
    returnCode, stdout, stderr = subprocessing.executeCommand( [ "gnuplot", "--version" ] )
    gnuplot_old = returnCode != 0

Options = { "xMin" : "*", "xMax" : "*", "yMin" : "*", "yMax" : "*", "xLabel" : "x", "yLabel" : "y", "title" : "", "xylog" : 0 }

def adjustRange( xy1, xy2 ) :

    if( xy1 == xy2 ) :
        dxy = 0.1 * abs( xy1 )
        if( xy1 == 0. ) : dxy = 1
        xy1 -= dxy
        xy2 += dxy
    return xy1, xy2

filesToPlot = 0
if( len( sys.argv ) > 1 ) :
    filesToPlot = 1
    InputFileName = sys.argv[-1]
    fi = open( InputFileName, "r" )
    ls = fi.readlines( )
    fi.close( )

    if( len( ls ) == 0 ) :
        xMin = yMin = -1.
        xMax = yMax =  1.
    else :
        s = ls[0].split()
        xMin = float( s[0] )
        xMax = xMin
        yMin = float( s[1] )
        yMax = yMin
    for l in ls :
        s = l.split( "#" )[0]
        s = s.split( )
        if( len( s ) > 0 ) :
            x = float( s[0] )
            xMin = min( xMin, x )
            xMax = max( xMax, x )
            y = float( s[1] )
            yMin = min( yMin, y )
            yMax = max( yMax, y )
    xMin, xMax = adjustRange( xMin, xMax )
    yMin, yMax = adjustRange( yMin, yMax )
    Options['xMin'] = xMin
    Options['xMax'] = xMax
    Options['yMin'] = yMin
    Options['yMax'] = yMax

g = Gnuplot.Gnuplot( )
g( 'set ticslevel 0.' )
g( "set nokey" )
if( filesToPlot ) : d = Gnuplot.File( InputFileName )

DoDelete = 0
i = 1
n = len( sys.argv ) - 1
while ( i < n ) :
    inc = 2
    if( sys.argv[i] == "delete" ) :
        inc = 1
        DoDelete = 1
    elif( sys.argv[i] in Options ) :
        Options[sys.argv[i]] = sys.argv[i+1].strip( '"' )
    else :
        print()
        print("Invalid option = <%s>" % sys.argv[i])
        print()
        for o in sys.argv : print(o, end="")
        sys.exit( )
    i = i + inc

ylog = int( Options["xylog"] )
xlog = ylog % 2
ylog /= 2
ylog = ylog % 2

def Replot( extraGnuplotCommand = None ) :
    """This function is called whenever the plot needs to be redrawn (e.g., axis is changed)."""

    if( not filesToPlot ) : return
    if( extraGnuplotCommand is not None ) : g( extraGnuplotCommand )

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
    d._options['with'] = ( w, 'with ' + w )
#    d.set_option( with=w )

    g( 'set xlabel "%s"' % xAxisVar.get( ) )
    g( 'set ylabel "%s"' % yAxisVar.get( ) )
    g( 'set title "%s"' % TitleVar.get( ) )
    g( "set nologscale xy" )
    xlog = xlogVar.get( )
    ylog = ylogVar.get( )
    if ( xlog or ylog ) :
        s = "set logscale "
        if ( xlog ) : s = s + "x"
        if ( ylog ) : s = s + "y"
        g( s )
    g.plot( d )

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
    """Called when "File -> SaveAs ASCII" menu is selected. Uses a FileDialog to get the output
file name and produces an ascii file of the plot data."""


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
    """A class that contains a Checkbutton widget that toggles the x- or y-axis between linear and log scaling."""

    def __init__( self, master, label, row, column, var ) :
        """The contructor for the logButton class."""

        self.label = label
        self.button = Checkbutton( frame, text=label, command = self.logcallback, variable = var )
        self.button.grid( row = row, column = column )

    def logcallback( self ) :
        """Called whenever the x- or y-axis log checkbutton state changes."""

        Replot( )

root = Tk( )
if( Options["title"] != "" ) :
    root.title( Options["title"] )
else :
    root.title( "endl2dplot.py" )

menubar = Menu( root )
filemenu = Menu( menubar, tearoff = 0 )
filemenu.add_command( label = "SaveAs eps", command = PlotSaveAsEps )
filemenu.add_command( label = "SaveAs ASCII", command = PlotSaveAsASCII )
filemenu.add_command( label = "Print", command = PlotPrint )
filemenu.add_separator( )
filemenu.add_command( label = "Exit", command = root.quit )
menubar.add_cascade( label = "File", menu = filemenu )
root.config(menu=menubar)

frame = Frame( root )
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
TitleLabel = Label( frame, text = "title = " )
TitleLabel.grid( row = Row, column = 0, sticky = E )
TitleVar = StringVar( )
TitleVar.set( Options["title"] )
TitleEntry = Entry( frame, textvariable = TitleVar )
TitleEntry.grid( row = Row, column = 1, columnspan = 3, sticky = W + E )
TitleEntry.bind( "<Return>", Entry_Callback )

Row += 1
s = "xrange [ %s : %s ]" % ( Options['xMin'], Options['xMax'] )
xlabel = Label( frame, text = s )
xlabel.grid( row = Row, column = 0, columnspan = 2, sticky = W, pady = 4 )

def xrange_Callback( e ) :
    """Called when the <Return> key is pressed in the xrange entry window."""

    s = "set xrange [ %s ]" % xrange.get( )
    Replot( s )

xrange = Entry( frame )
RangeStr = "%s : %s" % ( Options["xMin"], Options["xMax"] )
g( "set xrange [ %s ]" % RangeStr )
xrange.insert( 0, RangeStr )
xrange.grid( row = Row, column = 2 )
xrange.bind( "<Return>", xrange_Callback )
xlogVar = IntVar( )
xlogVar.set( xlog )
xLogButton = logButton( frame, "xlog", Row, 3, xlogVar )

Row += 1
s = "yrange [ %s : %s ]" % ( Options['yMin'], Options['yMax'] )
ylabel = Label( frame, text = s )
ylabel.grid( row = Row, column = 0, columnspan = 2, sticky = W, pady = 4 )

def yrange_Callback( e ) :
    """Called when the <Return> key is pressed in the yrange entry window."""

    s = "set yrange [ %s ]" % yrange.get( )
    Replot( s )

yrange = Entry( frame )
RangeStr = "%s : %s" % ( Options["yMin"], Options["yMax"] )
g( "set yrange [ %s ]" % RangeStr )
yrange.insert( 0, RangeStr )
yrange.grid( row = Row, column = 2 )
yrange.bind( "<Return>", yrange_Callback )
ylogVar = IntVar( )
ylogVar.set( ylog )
yLogButton = logButton( frame, "ylog", Row, 3, ylogVar )

if __name__ == '__main__' :
    Replot( )
    root.mainloop( )
    if( DoDelete ) : os.remove( InputFileName )
