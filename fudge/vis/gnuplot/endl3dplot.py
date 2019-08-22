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
This file creates an interactive gnuplot session that plots 3d data from a single file.
The data file must contain 3 columns (called 3d data) of numbers. The first column is 
plotted as the x data, the second column is plotted as the y data and the third is
plotted as the z data.  This file is called as, 

>>> python endl3dplot.py [xMin number] [xMax number] [yMin number] [yMax number] [zMin number] [zMax number] 
...        [xLabel string] [yLabel string] [zLabel string] [title string] [xyzlog int] [delete] <datafile>

where the parameters enclosed by [] are optional.

Parameters::
   xMin     The next argument is the initial minimum value of the x-axis
   xMax     The next argument is the initial maximum value of the x-axis
   yMin     The next argument is the initial minimum value of the y-axis
   yMax     The next argument is the initial maximum value of the y-axis
   zMin     The next argument is the initial minimum value of the z-axis
   zMax     The next argument is the initial maximum value of the z-axis
   xLabel   The next argument is the x-axis label enclosed in "'s (e.g., xLabel "Energy (MeV)")
   yLabel   The next argument is the y-axis label enclosed in "'s
   zLabel   The next argument is the z-axis label enclosed in "'s
   title    The next argument is the title enclosed in "'s
   xrot     The next argument is the xrot float value
   zrot     The next argument is the xrot float value
   xyzlog   If the next argument is even (odd) x-axis is linear (log), 
              this next argument divided by two is even (odd) y-axis is linear (log),
              this next argument divided by four is even (odd) z-axis is linear (log)
   delete   If present the input file is delete with this python script exits

Examples::
    python endl3dplot.py a
    python endl3dplot.py xLabel "x-axis" yLabel "y-axis" xyzlog 6 a    # x-axis is linear and y- and z-axis are log.
"""

import sys
import os
import string
import math
import shutil
from Tkinter import *
import tkFileDialog
import Gnuplot
from fudge.core.utilities import fudgeFileMisc

__metaclass__ = type

DoDelete = 0

filesToPlot = 0
xMin = 0
xMax = 1
yMin = 0
yMax = 1
zMin = 0
zMax = 1
if( len( sys.argv ) > 1 ) :
    filesToPlot = 1
    InputFileName = sys.argv[-1]
    fi = open( InputFileName, "r" )
    ls = fi.readlines( )
    fi.close( )
    f = fudgeFileMisc.fudgeTempFile( )

    s = string.split( ls[0] )
    lx = None
    xMin = float( s[0] )
    xMax = xMin
    yMin = float( s[1] )
    yMax = yMin
    zMin = float( s[2] )
    zMax = zMin
    for l in ls :
        s = string.split( l )
        x = float( s[0] )
        if ( x != lx ) : 
            if ( lx is not None ) : f.write( "\n" )
            lx = x
            xMin = min( xMin, x )
            xMax = max( xMax, x )
        y = float( s[1] )
        yMin = min( yMin, y )
        yMax = max( yMax, y )
        z = float( s[2] )
        zMin = min( zMin, z )
        zMax = max( zMax, z )
        s = "%15.7e %15.7e %15.7e\n" % ( x, y, z )
        f.write( s )

    f.close( )

g = Gnuplot.Gnuplot( )
g( 'set style data lines' )
g( 'set ticslevel 0.' )
g( "set nokey" )
if( filesToPlot ) : d = Gnuplot.File( f.getName( ) )

xlog = 0
ylog = 0
zlog = 0
xMinStr = "*"
xMaxStr = "*"
yMinStr = "*"
yMaxStr = "*"
zMinStr = "*"
zMaxStr = "*"
xLabel = "x"
yLabel = "y"
zLabel = "z"
Title = ""
xrot = 60
zrot = 30

i = 1
n = len( sys.argv ) - 1
while ( i < n ) :
    inc = 2
    if   ( sys.argv[i] == 'xMin' ) : xMinStr = sys.argv[i+1]
    elif ( sys.argv[i] == 'xMax' ) : xMaxStr = sys.argv[i+1]
    elif ( sys.argv[i] == 'yMin' ) : yMinStr = sys.argv[i+1]
    elif ( sys.argv[i] == 'yMax' ) : yMaxStr = sys.argv[i+1]
    elif ( sys.argv[i] == 'zMin' ) : zMinStr = sys.argv[i+1]
    elif ( sys.argv[i] == 'zMax' ) : zMaxStr = sys.argv[i+1]
    elif ( sys.argv[i] == 'xLabel' ) : xLabel = sys.argv[i+1]
    elif ( sys.argv[i] == 'yLabel' ) : yLabel = sys.argv[i+1]
    elif ( sys.argv[i] == 'zLabel' ) : zLabel = sys.argv[i+1]
    elif ( sys.argv[i] == 'title' ) : Title = sys.argv[i+1]
    elif ( sys.argv[i] == 'xrot' ) : xrot = float( sys.argv[i+1] )
    elif ( sys.argv[i] == 'zrot' ) : zrot = float( sys.argv[i+1] )
    elif ( sys.argv[i] == 'xyzlog' ) :
        zlog = int( sys.argv[i+1] )
        xlog = zlog % 2
        zlog /= 2
        ylog = zlog % 2
        zlog /= 2
        zlog = zlog % 2
    elif ( sys.argv[i] == 'delete' ) :
        inc = 1
        DoDelete = 1
    else :
        print
        print "Invalid option = <%s>" % sys.argv[i]
        print
        for o in sys.argv : print o,
        sys.exit( )
    i = i + inc

def Replot( c = None ) :
    """This function is called whenever the plot needs to be redrawn (e.g., axis is changed)."""

    if( c is not None ) : g( c )
    g( 'set xlabel "%s"' % xAxisVar.get( ) )
    g( 'set ylabel "%s"' % yAxisVar.get( ) )
    g( 'set zlabel "%s" rotate left' % zAxisVar.get( ) )
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

    s = "set view %e, %e" % ( xrot, zrot )

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

    g( s )
    g.splot( d )

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
        """Called whenever a xzRot Scale widget is changed."""

        global xrot, zrot
        if ( self.label == "x" ) :
            xrot = int( xz )
        else :
            zrot = int( xz )
        Replot( )

root = Tk( )
if( Title != "" ) :
    root.title( Title )
else :
    root.title( "endl3dplot.py" )

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
xAxisVar.set( xLabel )
xAxisEntry = Entry( frame, textvariable = xAxisVar )
xAxisEntry.grid( row = Row, column = 1, columnspan = 3, sticky = W + E )
xAxisEntry.bind( "<Return>", Entry_Callback )

Row += 1
yAxisLabel = Label( frame, text = "y label = " )
yAxisLabel.grid( row = Row, column = 0, sticky = E )
yAxisVar = StringVar( )
yAxisVar.set( yLabel )
yAxisEntry = Entry( frame, textvariable = yAxisVar )
yAxisEntry.grid( row = Row, column = 1, columnspan = 3, sticky = W + E )
yAxisEntry.bind( "<Return>", Entry_Callback )

Row += 1
zAxisLabel = Label( frame, text = "z label = " )
zAxisLabel.grid( row = Row, column = 0, sticky = E )
zAxisVar = StringVar( )
zAxisVar.set( zLabel )
zAxisEntry = Entry( frame, textvariable = zAxisVar )
zAxisEntry.grid( row = Row, column = 1, columnspan = 3, sticky = W + E )
zAxisEntry.bind( "<Return>", Entry_Callback )

Row += 1
TitleLabel = Label( frame, text = "title = " )
TitleLabel.grid( row = Row, column = 0, sticky = E )
TitleVar = StringVar( )
TitleVar.set( Title )
TitleEntry = Entry( frame, textvariable = TitleVar )
TitleEntry.grid( row = Row, column = 1, columnspan = 3, sticky = W + E )
TitleEntry.bind( "<Return>", Entry_Callback )

Row += 1
s_xrot = xzRot( frame, "x rot", 180, xrot, Row, 0, 2 )

Row += 1
s_zrot = xzRot( frame, "z rot", 360, zrot, Row, 0, 3 )

Row += 1
s = "xrange [ %e : %e ]" % ( xMin, xMax )
xlabel = Label( frame, text = s )
xlabel.grid( row = Row, column = 0, columnspan = 2, sticky = W, pady = 4 )

def xrange_Callback( e ) :
    """Called when the <Return> key is pressed in the xrange entry window."""

    s = "set xrange [ %s ]" % xrange.get( )
    Replot( s )

xrange = Entry( frame )
RangeStr = "%s : %s" % ( xMinStr, xMaxStr )
g( "set xrange [ %s ]" % RangeStr )
xrange.insert( 0, RangeStr )
xrange.grid( row = Row, column = 2 )
xrange.bind( "<Return>", xrange_Callback )
xlogVar = IntVar( )
xlogVar.set( xlog )
xLogButton = logButton( frame, "xlog", Row, 3, xlogVar )

Row += 1
s = "yrange [ %e : %e ]" % ( yMin, yMax )
ylabel = Label( frame, text = s )
ylabel.grid( row = Row, column = 0, columnspan = 2, sticky = W, pady = 4 )

def yrange_Callback( e ) :
    """Called when the <Return> key is pressed in the yrange entry window."""

    s = "set yrange [ %s ] " % yrange.get( )
    Replot( s )

yrange = Entry( frame )
RangeStr = "%s : %s" % ( yMinStr, yMaxStr )
g( "set yrange [ %s ]" % RangeStr )
yrange.insert( 0, RangeStr )
yrange.grid( row = Row, column = 2 )
yrange.bind( "<Return>", yrange_Callback )
ylogVar = IntVar( )
ylogVar.set( ylog )
yLogButton = logButton( frame, "ylog", Row, 3, ylogVar )

Row += 1
s = "zrange [ %e : %e ]" % ( zMin, zMax )
zlabel = Label( frame, text = s )
zlabel.grid( row = Row, column = 0, columnspan = 2, sticky = W, pady = 4 )

def zrange_Callback( e ) :
    """Called when the <Return> key is pressed in the zrange entry window."""

    s = "set zrange [ %s ]" % zrange.get( )
    Replot( s )

zrange = Entry( frame )
RangeStr = "%s : %s" % ( zMinStr, zMaxStr )
g( "set zrange [ %s ]" % RangeStr )
zrange.insert( 0, RangeStr )
zrange.grid( row = Row, column = 2 )
zrange.bind( "<Return>", zrange_Callback )
zlogVar = IntVar( )
zlogVar.set( zlog )
zLogButton = logButton( frame, "zlog", Row, 3, zlogVar )

if __name__ == '__main__' :
    root.mainloop( )
    os.remove( f.getName( ) )
    if( DoDelete ) : os.remove( InputFileName )
