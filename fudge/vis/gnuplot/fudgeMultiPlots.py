# <<BEGIN-copyright>>
# <<END-copyright>>

"""
This module contains routines for plotting multiple data sets on the same plot.
"""

import os, types
from fudge.core import fudgemisc
from fudge.core.utilities import fudgeFileMisc, subprocessing
from fudge.core.math import endl2dmathmisc, endl3dmathmisc, endl3dmathClasses
from fudge import __path__
import plotbase


def multiPlot( dataList, xylog = 0, xMin = None, xMax = None, yMin = None, yMax = None, xLabel = "Energy (MeV)",
        yLabel = "Cross section (barn)", title = None, lineWidth = 1, fontSize = None ) :
    """This routine sends all data objects of the first argument, dataList, to the module fudge2dMultiPlot for interactive plotting. 
Here, dataList must be a list or tuple, and each element must be an endl2dmath object or subclass of it. 
The legend for each data element of dataList is taken from that object's label member."

See module fudge2dMultiPlot for a description of the following parameters::
    xylog    xMin     xMax     yMin     
    yMax     xLabel   yLabel   title    

Currently, lineWitdh and fontSize are not implemented.

Examples of usages where d1, d2 and d3 are endl2dmath objects.

>>> d2.label = "Bad data from boss"                     # Sets the legend for d2 to "Bad data from boss".
>>> multiPlot( ( d1, d2, d3 ), xylog = 3 )
>>> I0Data = za.findDatas( I = 0, C = ( 12, 13, 14 ) )  # za is an endlZA object.
>>> multiPlot( I0Data )

Also see the routine qmultiPlot.
"""
    from fudge.core.math.xData import XYs

    if( len( dataList ) == 0 ) : return
    if( ( type( dataList ) != type( [] ) ) and ( type( dataList ) != type( () ) ) ) : dataList = [ dataList ]
    fs = []
    for a in dataList :
        isInstance = ( type( a ) == types.InstanceType )
        if( hasattr( a, '__class__' ) ) : isInstance |= ( type( a ) == a.__class__ )    # New style class instance.
        if( isInstance ) :
            Msg = "cannot plot object of type %s" % a.__class__.__name__
            if( isinstance( a, XYs.XYs ) ) :
                f = fudgeFileMisc.fudgeTempFile( )
                fs.append( f.getName( ) )
                for x, y in a : f.write( "%.10e %14.8e\n" % ( x, y ) )
                f.close( )
            elif( endl2dmathmisc.valid2dClassType( a, "multiPlot", Msg ) ) :
                f = fudgeFileMisc.fudgeTempFile( )
                fs += [ f.getName( ), 'title', a.label ]
                for x, y in a.data : f.write( "%.10e %14.8e\n" % ( x, y ) )
                f.close( )
            else :
                print "Warning in multiPlot: cannot plot object named %s" % a.__class__.__name__
        else :
            raise Exception( "\nError in multiPlot: cannot plot none instance object" )
    o = plotbase.parsePlotOptions( xMin, xMax, yMin, yMax, xLabel, yLabel, title )
    p = fudgemisc.findPythonFile( os.sep.join( __path__ + [ "vis", "gnuplot", "fudge2dMultiPlot.py" ] ) )
    s = [ 'python', p, 'xylog', str( xylog ) ] + o + [ 'files' ] + fs
    subprocessing.spawn( s )

def qmultiPlot( dataList, xylog = 0, xMin = None, xMax = None, yMin = None, yMax = None, xLabel = "Energy (MeV)", 
        yLabel = "Cross section (barn)", title = None, legends = [], lineWidth = 1, fontSize = None ) :
    """This routine plots the data objects of the first argument, dataList, to a non-interactive plot.
The argument lineWidth sets the width of the lines and fontSize sets the font size for the plot.
The argument legends is a deprecated argument. To specify a legend for a data set, use that data's
label member.

Examples of usages where d1, d2 and d3 are endl2dmath objects.

>>> d2.label = "Bad data from boss"                     # Sets the legend for d2 to "Bad data from boss".
>>> multiPlot( ( d1, d2, d3 ), xylog = 3 )
>>> I0Data = za.findDatas( I = 0, C = ( 12, 13, 14 ) )  # za is an endlZA object.
>>> multiPlot( I0Data )

See the routine multiPlot for additional information."""

    import Gnuplot
    xylog = int( xylog )            # Allow argument to be a string
    if( xMin != None ) : xMin = float( xMin )
    if( xMax != None ) : xMax = float( xMax )
    if( yMin != None ) : yMin = float( yMin )
    if( yMax != None ) : yMax = float( yMax )
    lineWidth = float( lineWidth )

    def qmultiPlotAddPlot( d, g, t, withLineWidth ) :

        f = fudgeFileMisc.fudgeTempFile( )
        for i in d.data :
            d = "%.12e %.12e\n" % ( i[0], i[1] )
            f.write( d )
        f.close( )
        gf = Gnuplot.File( f.getName( ), title = t ) # , with = withLineWidth )
        gf._options['with'] = ( withLineWidth, 'with ' + withLineWidth )
        g._add_to_queue( [ gf ] )

    if( len( dataList ) == 0 ) : return
    g = Gnuplot.Gnuplot( )
    withLineWidth = 'line linewidth %f' % lineWidth
    nl = len( legends )
    il = 0
    if( ( type( dataList ) != type( [] ) ) and ( type( dataList ) != type( () ) ) ) : dataList = [ dataList ]
    for a in dataList :
        t = a.label
        if( il < nl ) : t = legends[il]
        if( type( a ) == types.InstanceType ) :
            Msg = "cannot plot object of type %s" % a.__class__.__name__
            if( endl2dmathmisc.valid2dClassType( a, "qmultiPlot", Msg ) ) :
                qmultiPlotAddPlot( a, g, t, withLineWidth )
            elif( a.__class__.__name__ == "endlZA" ) :
                for F in a.filelist :
                    if( F.levels != [] ) and ( F.levels.columns( ) == 2 ) :
                        for i in F.levels : qmultiPlotAddPlot( i, g, t, withLineWidth )
            else :
                print "Warning in qmultiPlot: cannot plot object named %s" % a.__class__.__name__
            il = il + 1
    if( title == None ) : title = "Untitled"
    g( 'set title "' + title + '"' )
    g( 'set style data linespoints' )
    if  ( xylog == 1 ) : g( 'set logscale x' )
    elif( xylog == 2 ) : g( 'set logscale y' )
    elif( xylog == 3 ) : g( 'set logscale xy' )
    if( xMin != None ) or ( xMax != None ) :
        if( xMin == None ) :
            s = 'set xrange [ * to %e ]' % xMax
        elif( xMax == None ) :
            s = 'set xrange [ %e to * ]' % xMin
        else :
            s = 'set xrange [ %e to %e ]' % ( xMin, xMax )
        g( s )
    if( yMin != None ) or ( yMax != None ) :
        if( yMin == None ) :
            s = 'set yrange [ * to %e ]' %  yMax
        elif( yMax == None ) :
            s = 'set yrange [ %e to * ]' % yMin
        else :
            s = 'set yrange [ %e to %e ]' % ( yMin, yMax )
        g( s )
    if( fontSize != None ) : g( 'set terminal x11 font "times,%d"' % fontSize )
    g( "set xlabel %s" % `xLabel` )
    g( "set ylabel %s" % `yLabel` )
    g.replot( )
    return g

def multiPlot3d( dataList, xyzlog = 0, xMin = None, xMax = None, yMin = None, yMax = None, zMin = None, zMax = None, xLabel = None,
        yLabel = None, zLabel = None, title = None, lineWidth = 1, fontSize = None ) :
    """This routine sends all data objects of the first argument, dataList, to the module fudge3dMultiPlot for interactive plotting.
Here, dataList must be a list or tuple, and each element must be an endl3dmath object or subclass of it.
The legend for each data element of dataList is taken from that object's label member."

See module fudge3dMultiPlot for a description of the following parameters::
    xyzlog   xMin     xMax     yMin     yMax
    zMin     zMax     xLabel   yLabel   zLabel
    title

Currently, lineWitdh and fontSize are not implemented.

Examples of usages where d1, d2 and d3 are endl3dmath objects.

>>> d2.label = "Bad data from boss"                     # Sets the legend for d2 to "Bad data from boss".
>>> multiPlot3d( ( d1, d2, d3 ), xyzlog = 3 )           # x log, y log and z linear.
>>> I1Data = za.findDatas( I = 1, C = ( 12, 13, 14 ) )  # za is an endlZA object.
>>> multiPlot( I1Data )
"""

    if( len( dataList ) == 0 ) : return
    if( ( type( dataList ) != type( [] ) ) and ( type( dataList ) != type( () ) ) ) : dataList = [ dataList ]
    fs = []
    for a in dataList :
        if( isinstance( a, endl3dmathClasses.endl3dmath ) ) :
            f = fudgeFileMisc.fudgeTempFile( )
            fs += [ f.getName( ) ]
            if( hasattr( a, 'label' ) ) : fs += [ 'title', a.label ]
            for x, yz in a.data :
                xs = "%.10e" % x
                for y, z in yz : f.write( "%s %.10e %14.8e\n" % ( xs, y, z ) )
                f.write( "\n" )
            f.close( )
            if( a == dataList[0] ) :
                if( xLabel == None ) : xLabel = a.xLabel
                if( yLabel == None ) : yLabel = a.yLabel
                if( zLabel == None ) : zLabel = a.zLabel
        else :
            raise Exception( "\nError in multiPlot: cannot plot instance object" )
    o = plotbase.parsePlotOptions( xMin, xMax, yMin, yMax, xLabel, yLabel, title, zMin, zMax, zLabel )
    p = fudgemisc.findPythonFile( os.sep.join( __path__ + [ "vis", "gnuplot", "fudge3dMultiPlot.py" ] ) )
    s = [ 'python', p, 'xyzlog', str( xyzlog ) ] + o + [ 'files' ] + fs
    subprocessing.spawn( s )
