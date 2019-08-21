# <<BEGIN-copyright>>
# <<END-copyright>>

"""
This file creates an interactive gnuplot session that plots 4d data from a single file.
The data file must contain 4 columns (called 4d data) of numbers. The first column is plotted 
as the t data, the second column is plotted as the x data, the third is plotted as the y 
data and the fourth column is plotted as the z data.  The data are broken up into a "pseudo
time" series of nTime 3d data sets where nTime is the number of distinct values in the first
column (assuming all like values are contiguous). A 3d data set is the x, y and z data for a 
given contiguous t datum.  Only one of the 3d data sets is plotted at one time. Mechanisms 
exist for stepping and animating through the nTime "time" series. This file is called as,

>>> python endl4dplot.py [xMin number] [xMax number] [yMin number] [yMax number] [zMin number] [zMax number]
...        [xLabel string] [yLabel string] [zLabel string] [title string] [xyzlog int] [delete] <datafile>

where the parameters enclosed by [] are optional.

Parameters::
   xMin        The next argument after xMin is the initial minimum value of the x-axis
   xMax        The next argument aftet xMax is the initial maximum value of the x-axis
   yMin        The next argument after yMin is the initial minimum value of the y-axis
   yMax        The next argument after yMax is the initial maximum value of the y-axis
   zMin        The next argument after zMin is the initial minimum value of the z-axis
   zMax        The next argument after zMax is the initial maximum value of the z-axis
   tLabel      The next argument after tLabel is the string used to contruct a title (i.e., "set title \"%s = %e\"" % ( tLabel, ... ) )
   xLabel      The next argument after xLabel is the x-axis label enclosed in "'s (e.g., xLabel "Energy (MeV)")
   yLabel      The next argument after yLabel is the y-axis label enclosed in "'s
   zLabel      The next argument after zLabel is the z-axis label enclosed in "'s
   title       The next argument after title is the title enclosed in "'s
   xrot        The next argument after xrot is the xrot float value
   zrot        The next argument after zrot is the xrot float value
   xyzlog      If the next argument after xyzlog is even (odd) x-axis is linear (log),
                 this next argument divided by two is even (odd) y-axis is linear (log),
                 this next argument divided by four is even (odd) z-axis is linear (log)
   tScaleLabel The next argument after tScaleLabel is the string used to contruct a the label for the tScale
   delete      If present the input file is delete with this python script exits

Examples::
    python endl4dplot.py a
    python endl4dplot.py xLabel "x-axis" yLabel "y-axis" xyzlog 6 a    # x-axis is linear and y- and z-axis are log.
"""

import sys
import os
import string
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
nE = 2
if( len( sys.argv ) > 1 ) :
    filesToPlot = 1
    InputFileName = sys.argv[-1]
    fi = open( InputFileName, "r" )
    ls = fi.readlines( )
    fi.close( )
    Enfs = []

    s = string.split( ls[0] )
    lE = None
    nE = 0
    lx = None
    xMin = float( s[1] )
    xMax = xMin
    yMin = float( s[2] )
    yMax = yMin
    zMin = float( s[3] )
    zMax = zMin
    for l in ls :
        s = string.split( l )
        Ei = float( s[0] )
        if ( Ei != lE ) :
            nE += 1
            if ( lE != None ) :
                f.close( )
            f = fudgeFileMisc.fudgeTempFile( )
            Enfs.append( ( Ei, f, Gnuplot.File( f.getName( ) ) ) )
            lE = Ei
            lx = None
        x = float( s[1] )
        if ( x != lx ) : 
            if ( lx != None ) : f.write( "\n" )
            lx = x
            xMin = min( xMin, x )
            xMax = max( xMax, x )
        y = float( s[2] )
        yMin = min( yMin, y )
        yMax = max( yMax, y )
        z = float( s[3] )
        zMin = min( zMin, z )
        zMax = max( zMax, z )
        s = "%15.7e %15.7e %15.7e\n" % ( x, y, z )
        f.write( s )

    f.close( )

g = Gnuplot.Gnuplot( )
g( 'set style data lines' )
g( 'set ticslevel 0.' )
g( "set nokey" )

xlog = 0
ylog = 0
zlog = 0
xMinStr = "*"
xMaxStr = "*"
yMinStr = "*"
yMaxStr = "*"
zMinStr = "*"
zMaxStr = "*"
tLabel = "t"
xLabel = "x"
yLabel = "y"
zLabel = "z"
Title = ""
tScaleLabel = "t scale"
xrot = 60
zrot = 30
iE = 0

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
    elif ( sys.argv[i] == 'tLabel' ) : tLabel = sys.argv[i+1]
    elif ( sys.argv[i] == 'xLabel' ) : xLabel = sys.argv[i+1]
    elif ( sys.argv[i] == 'yLabel' ) : yLabel = sys.argv[i+1]
    elif ( sys.argv[i] == 'zLabel' ) : zLabel = sys.argv[i+1]
    elif ( sys.argv[i] == 'title' ) : Title = sys.argv[i+1]
    elif ( sys.argv[i] == 'tScaleLabel' ) : tScaleLabel = sys.argv[i+1]
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

def Replot( extraGnuplotCommand = None ) :
    """This function is called whenever the plot needs to be redrawn (e.g., axis is changed)."""

    if( not filesToPlot ) : return
    if( extraGnuplotCommand != None ) : g( extraGnuplotCommand )
    g( 'set xlabel "%s"' % xAxisVar.get( ) )
    g( 'set ylabel "%s"' % yAxisVar.get( ) )
    g( 'set zlabel "%s"' % zAxisVar.get( ) )
    Title = TitleVar.get( )
    if ( Title == "" ) :
        s = "set title \"%s = %e\"" % ( tLabel, Enfs[iE][0] )
        g( s )
    else :
        g( 'set title "%s"' % Title )
    g( "set nologscale xyz" )
    hidden = hiddenVar.get( )
    g( 'set nohidden3d' )
    if( hidden ) : g( 'set hidden3d' )
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
    g( s )
    g.splot( Enfs[iE][2] )

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

class AnimateButton :
    """A class that contains a Checkbutton for toggling animations on and off."""

    def __init__( self, master, r, c, var ) :
        """The contructor for the AnimateButton class."""

        self.button = Checkbutton( frame, text="Animate", command = self.acallback, variable = var )
        self.button.grid( row = r, column = c )

    def acallback( self ) :
        """Called whenever the animation button is toggled."""

        if ( AnimateVar.get( ) ) : root.after( 0, self.DoAnimateCallback, 0 )

    def DoAnimateCallback( self, i ) :
        """Called periodically (see timeScale class for delay) by the scheduler when the animation button is depressed."""

        if ( AnimateVar.get( ) ) :
            if ( i >= nE ) : i = 0
            s_tScale.scale.set( i )
            Replot( )
            root.after( AnimateTimeVar.get( ) * 100, self.DoAnimateCallback, i + 1 )

class timeScale :
    """A class that contains a Scale widget for controlling the speed of animations."""

    def __init__( self, master, row, column, var ) :
        """The contructor for the timeScale class."""

        self.scale = Scale( master, label = "time (sec * 10)", length = 100, from_ = 1, to = 40, 
            orient = "horizontal", variable = var )
        self.scale.grid( row = row, column = column )

class tScale :
    """A class that contains a Scale widget for stepping through the t-axis data."""

    def __init__( self, master, to_, row, column ) :
        """The contructor for the tScale class."""

        l = 100
        if ( to_ > 100 ) : l = to_
        self.scale = Scale( master, label = tScaleLabel, length = l, from_ = 0, to = to_, 
            orient = "horizontal", command = self.ecallback )
        self.scale.set( 0 )
        self.scale.grid( row = row, column = column )

    def ecallback( self, s ) :
        """Called whenever a tScale Scale widget is changed."""

        global iE
        iE = int( s )
        Replot( )

class xyRot :
    """A class that contains a Scale widget for rotating the plot about an axis."""

    def __init__( self, master, label, to_, value, row, column, span ) :
        """The contructor for the xyRot class."""

        self.label = label[0]
        self.scale = Scale( master, label = label, length = to_, from_ = 0, to = to_, bigincrement = 5, 
            orient = "horizontal", command = self.xzcallback )
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
if( Title != "" ) : root.title( Title )

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
s_tScale = tScale( frame, nE - 1, Row, 0 )
AnimateVar = IntVar( )
b_Animate = AnimateButton( frame, Row, 1, AnimateVar )
AnimateTimeVar = IntVar( )
AnimateTimeVar.set( 10 )
s_timeScale = timeScale( frame, Row, 2, AnimateTimeVar )

Row += 1
s_xrot = xyRot( frame, "x rot", 180, xrot, Row, 0, 2 )
hiddenVar = IntVar( )
hiddenVar.set( 0 )
hiddenButton = logButton( frame, "Hidden", Row, 3, hiddenVar )

Row += 1
s_zrot = xyRot( frame, "z rot", 360, zrot, Row, 0, 3 )

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

    s = "set yrange [ %s ]" % yrange.get( )
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
    for f in Enfs : os.remove( f[1].getName( ) )
    if ( DoDelete ) : os.remove( InputFileName )
