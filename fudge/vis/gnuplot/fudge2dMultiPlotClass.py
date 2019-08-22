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
This file contains the class fudge2dMultiPlotClass. This class constructs a new window for displaying
fudgeNdMultiPlotData and fudgeNdMultiPlotFile objects in a single plot window.
"""

import os, shutil
from _tkinter import TclError
import Tix
import Tkinter
import tkFileDialog
import Gnuplot
import fudgeNdMultiPlotMisc
from fudge.core.utilities import subprocessing

__metaclass__ = type

returnCode, stdout, stderr = subprocessing.executeCommand( [ "gnuplot", "--version" ] )
gnuplot_old = returnCode != 0

Options = { "xMin" : "*", "xMax" : "*", "yMin" : "*", "yMax" : "*", "xLabel" : "x", "yLabel" : "y", "title" : "", "xlog" : 0, "ylog" : 0 }
saveAsCounter = 0

class logButton :
    """A class that contains a Checkbutton widget that toggles the x- or y-axis between linear and log scaling."""

    def __init__( self, plotter, label, row, column, var ) :
        """The contructor for the logButton class."""

        self.plotter = plotter
        self.label = label
        self.button = Tkinter.Checkbutton( plotter.frame, text = label, command = self.logcallback, variable = var )
        self.button.grid( row = row, column = column )

    def logcallback( self ) :
        """Called whenever the x- or y-axis log checkbutton state changes."""

        self.plotter.replot( )

class fudge2dMultiPlotClass( Tkinter.Frame ) :

    def __init__( self, rootWindow = None, options = {} ) :

        self.lastTitle = None
        if( rootWindow is None ) : rootWindow = Tkinter.Toplevel( )
        self.rootWindow = rootWindow
        Tkinter.Frame.__init__( self, rootWindow )
        self.options = {}
        for option in Options : self.options[option] = Options[option]
        for option in options :
            if( option in [ 'xMin', 'xMax', 'yMin', 'yMax' ] ) :
                self.options[option] = options[option]
            elif( option in [ 'xLabel', 'yLabel', 'title' ] ) :
                if( type( options[option] ) != type( '' ) ) : raise Exception( 'option %s value must be a string' )
                self.options[option] = options[option]
            elif( option in [ 'xlog', 'ylog' ] ) :
                self.options[option] = 0
                if( options[option] ) : self.options[option] = 1
            else :
                raise Exception( 'Unsupported option = %s' % options )

        menubar = Tkinter.Menu( rootWindow )
        self.frame = Tkinter.LabelFrame( self, text = "Plot Options" )
        self.frame.grid( row = 0, sticky = Tkinter.W + Tkinter.E )
        self.frame.grid( )
        plotItemsFrame = Tkinter.Frame( self )
        plotItemsFrame.grid( row = 1, sticky = Tkinter.W + Tkinter.E )
        self.grid( )
        self.savedTitle = rootWindow.title( )

        self.plotItemEditting = fudgeNdMultiPlotMisc.fudgeNdMultiPlotItemsDialog( plotItemsFrame, self.redraw, activeMenubar = menubar )
        self.g = Gnuplot.Gnuplot( )
        self.g( 'set ticslevel 0.' )

        optionsMenu = Tkinter.Menu( menubar, tearoff = 0 )
        self.legendOnVar = Tkinter.IntVar( )
        optionsMenu.add_checkbutton( label = "Legend", variable = self.legendOnVar, command = self.legendCallback )
        menubar.insert_cascade( 0, label = "Options", menu = optionsMenu )

        dataMenu = Tkinter.Menu( menubar, tearoff = 0 )
        dataMenu.add_command( label = "Read file", command = self.readFile )
        menubar.insert_cascade( 0, label = "Data", menu = dataMenu )

        fileMenu = Tkinter.Menu( menubar, tearoff = 0 )

        savePlotAs = Tkinter.Menu( fileMenu, tearoff = 0 )
        savePlotAs.add_command( label = "pdf", command = self.plotSaveAsPDF )
        savePlotAs.add_command( label = "png", command = self.plotSaveAsPNG )
        savePlotAs.add_command( label = "svg", command = self.plotSaveAsSVG )
        savePlotAs.add_command( label = "eps", command = self.plotSaveAsEps )
        fileMenu.add_cascade( label = "Save plot as", menu = savePlotAs )


        fileMenu.add_command( label = "Print", command = self.plotPrint )
        fileMenu.add_separator( )
        fileMenu.add_command( label = "Exit", command = rootWindow.quit )
        menubar.insert_cascade( 0, label = "File", menu = fileMenu )
        rootWindow.config( menu = menubar )

        def entry_Callback( self, e ) :
            """Called whenever the <Return> key is pressed in certain entry window."""
            self.replot( )

        Row = 0
        if( not gnuplot_old ) :
            fontWidthLabel = Tkinter.Label( self.frame, text = "font size = " )
            fontWidthLabel.grid( row = Row, column = 0, sticky = Tkinter.E )
            self.fontSizeVar = Tkinter.IntVar( )
            self.fontSizeVar.set( 12 )
            fontSizeEntry = Tkinter.Entry( self.frame, textvariable = self.fontSizeVar )
            fontSizeEntry.grid( row = Row, column = 1, sticky = Tkinter.W + Tkinter.E )
            fontSizeEntry.bind( "<Return>", lambda e, self = self  : entry_Callback( self, e ) )
            self.legendBoxVar = Tkinter.IntVar( )
            self.legendBoxVar.set( 0 )
            b = logButton( self, "Legend box", Row, 2, self.legendBoxVar )
            Row += 1
        ticsSize = Tkinter.Label( self.frame, text = "tics size = " )
        ticsSize.grid( row = Row, column = 0, sticky = Tkinter.E )
        self.ticsVar = Tkinter.DoubleVar( )
        self.ticsVar.set( 1.5 )
        ticsSizeEntry = Tkinter.Entry( self.frame, textvariable = self.ticsVar )
        ticsSizeEntry.grid( row = Row, column = 1, sticky = Tkinter.W + Tkinter.E )
        ticsSizeEntry.bind( "<Return>", lambda e, self = self  : entry_Callback( self, e ) )

        Row += 1
        l = Tkinter.Label( self.frame, text = "legend x,y =" )
        l.grid( row = Row, column = 0, sticky = Tkinter.E )
        self.legendX = Tkinter.DoubleVar( )
        self.legendX.set( .98 )
        e = Tkinter.Entry( self.frame, textvariable = self.legendX )
        e.grid( row = Row, column = 1, sticky = Tkinter.W + Tkinter.E )
        e.bind( "<Return>", lambda e, self = self  : entry_Callback( self, e ) )
        self.legendY = Tkinter.DoubleVar( )
        self.legendY.set( .95 )
        e = Tkinter.Entry( self.frame, textvariable = self.legendY )
        e.grid( row = Row, column = 2, sticky = Tkinter.W + Tkinter.E )
        e.bind( "<Return>", lambda e, self = self  : entry_Callback( self, e ) )

        Row += 1
        xAxisLabel = Tkinter.Label( self.frame, text = "x label = " )
        xAxisLabel.grid( row = Row, column = 0, sticky = Tkinter.E )
        self.xAxisVar = Tkinter.StringVar( )
        self.xAxisVar.set( self.options["xLabel"] )
        xAxisEntry = Tkinter.Entry( self.frame, textvariable = self.xAxisVar )
        xAxisEntry.grid( row = Row, column = 1, columnspan = 3, sticky = Tkinter.W + Tkinter.E )
        xAxisEntry.bind( "<Return>", lambda e, self = self  : entry_Callback( self, e ) )

        Row += 1
        yAxisLabel = Tkinter.Label( self.frame, text = "y label = " )
        yAxisLabel.grid( row = Row, column = 0, sticky = Tkinter.E )
        self.yAxisVar = Tkinter.StringVar( )
        self.yAxisVar.set( self.options["yLabel"] )
        yAxisEntry = Tkinter.Entry( self.frame, textvariable = self.yAxisVar )
        yAxisEntry.grid( row = Row, column = 1, columnspan = 3, sticky = Tkinter.W + Tkinter.E )
        yAxisEntry.bind( "<Return>", lambda e, self = self  : entry_Callback( self, e ) )

        Row += 1
        TitleLabel = Tkinter.Label( self.frame, text = "title = " )
        TitleLabel.grid( row = Row, column = 0, sticky = Tkinter.E )
        self.TitleVar = Tkinter.StringVar( )
        self.TitleVar.set( self.options["title"] )
        TitleEntry = Tkinter.Entry( self.frame, textvariable = self.TitleVar )
        TitleEntry.grid( row = Row, column = 1, columnspan = 3, sticky = Tkinter.W + Tkinter.E )
        TitleEntry.bind( "<Return>", lambda e, self = self  : entry_Callback( self, e ) )

        Row += 1
        xlabel = Tkinter.Label( self.frame, text = "xrange" )
        xlabel.grid( row = Row, column = 0, columnspan = 1, sticky = Tkinter.E, pady = 4 )

        def xrange_Callback( self, e ) :
            """Called when the <Return> key is pressed in the xrange entry window."""
            s = "set xrange [ %s ] noreverse" % self.xrange.get( )
            self.replot( s )

        self.xrange = Tkinter.Entry( self.frame )
        RangeStr = "%s : %s" % ( self.options["xMin"], self.options["xMax"] )
        self.g( "set xrange [ %s ] noreverse" % RangeStr )
        self.xrange.insert( 0, RangeStr )
        self.xrange.grid( row = Row, column = 1, columnspan = 2, sticky = Tkinter.W + Tkinter.E )
        self.xrange.bind( "<Return>", lambda e, self = self : xrange_Callback( self, e ) )
        self.xlogVar = Tkinter.IntVar( )
        self.xlogVar.set( self.options['xlog'] )
        xLogButton = logButton( self, "xlog", Row, 3, self.xlogVar )

        Row += 1
        ylabel = Tkinter.Label( self.frame, text = "yrange" )
        ylabel.grid( row = Row, column = 0, columnspan = 1, sticky = Tkinter.E, pady = 4 )

        def yrange_Callback( self, e ) :
            """Called when the <Return> key is pressed in the yrange entry window."""
            s = "set yrange [ %s ] noreverse" % self.yrange.get( )
            self.replot( s )

        self.yrange = Tkinter.Entry( self.frame )
        RangeStr = "%s : %s" % ( self.options["yMin"], self.options["yMax"] )
        self.g( "set yrange [ %s ] noreverse" % RangeStr )
        self.yrange.insert( 0, RangeStr )
        self.yrange.grid( row = Row, column = 1, columnspan = 2, sticky = Tkinter.W + Tkinter.E )
        self.yrange.bind( "<Return>", lambda e, self = self : yrange_Callback( self, e ) )
        self.ylogVar = Tkinter.IntVar( )
        self.ylogVar.set( self.options['ylog'] )
        yLogButton = logButton( self, "ylog", Row, 3, self.ylogVar )

        try:
            def gnuplotCommand_Callback2( self, e ) :
                """Called when the <Return> key is pressed in the gnuplotCommand entry window."""
                self.replot( self.gnuplotCommandComboBoxVar.get( ) )

            Row += 1
            self.gnuplotCommandComboBox = Tix.ComboBox( self.frame )
            self.gnuplotCommandComboBox.grid( row = Row, column = 0, columnspan = 4, sticky = Tkinter.W + Tkinter.E )
            self.gnuplotCommandComboBox['editable'] = True
            self.gnuplotCommandComboBox['prunehistory'] = False
            self.gnuplotCommandComboBox['label'] = 'Plot command'
            self.gnuplotCommandComboBox['selectmode'] = 'browse'
            self.gnuplotCommandComboBox['command'] = lambda e, self = self : gnuplotCommand_Callback2( self, e )
            self.gnuplotCommandComboBoxVar = Tkinter.StringVar( )
            self.gnuplotCommandComboBox['variable'] = self.gnuplotCommandComboBoxVar
            self.gnuplotCommandComboBox['history'] = True
            self.gnuplotCommandComboBox['state'] = 'normal'
            self.gnuplotCommandComboBox.entry['state'] = 'normal'
        except TclError: pass

    def redraw( self ) :

        self.replot( )

    def replot( self, extraGnuplotCommand = None ) :
        """This function is called whenever the plot needs to be redrawn (e.g., axis is changed)."""

        title = self.TitleVar.get( )
        self.g( 'set title "%s"' % title )
        if( title == '' ) : title = self.savedTitle
        dialogTitle = title + ' -- dialog'
        plotTitle = title + ' -- plot'
        if( self.lastTitle != plotTitle ) : self.g( 'set terminal x11 title "%s"' % plotTitle )
        self.lastTitle = plotTitle
        self.rootWindow.title( dialogTitle )

        legendBox = "no"
        if( not gnuplot_old ) :
           if( self.legendBoxVar.get( ) ) : legendBox = ""
           self.g( "set tics scale %e" % self.ticsVar.get( ) )
           fs_ = self.fontSizeVar.get( )
           fs = max( 8, min( fs_, 72 ) )
           self.g( 'set terminal X11 font "times,%d"' % fs )

        self.g( "set nokey" )
        if( self.legendOnVar.get( ) ) : self.g( "set key %sbox right top at graph %e, %e" % ( legendBox, self.legendX.get( ), self.legendY.get( ) ) )

        self.g( 'set xlabel "%s"' % self.xAxisVar.get( ) )
        self.g( 'set ylabel "%s"' % self.yAxisVar.get( ) )
        self.g( "set nologscale xy" )
        if( self.xlogVar.get( ) ) : self.g( "set logscale x" )
        if( self.ylogVar.get( ) ) : self.g( "set logscale y" )

        self.g( 'set format x "%.10g"' )
        self.g( 'set format y "%.10g"' )

        activeCurves = self.plotItemEditting.getActiveCurves( )
        if( extraGnuplotCommand is not None ) : self.g( extraGnuplotCommand )
        if( len( activeCurves ) ) : apply( self.g.plot, activeCurves )

    def add2dData( self, data, title, active = True, replot = True ) :

        fMPF = fudgeNdMultiPlotMisc.fudgeNdMultiPlotData( data, title, dimension = 2, active = active, replot = self.replot )
        self.plotItemEditting.addFudgeNdMultiPlotFile( fMPF )
        if( replot ) : self.replot( )

    def add2dFile( self, fileName, title, xColumn = 1, yColumn = 2, active = True, replot = True ) :

        fMPF = fudgeNdMultiPlotMisc.fudgeNdMultiPlotFile( fileName, 2, title = title, active = active, xColumn = xColumn, yColumn = yColumn )
        self.plotItemEditting.addFudgeNdMultiPlotFile( fMPF )
        if( replot ) : self.replot( )
        
    def readFile( self ) :

        fileName, title, columns = fudgeNdMultiPlotMisc.fudgeNdMultiPlotFileRead( self.rootWindow, 2 )
        if( columns is not None ) : self.add2dFile( fileName, title, xColumn = columns[0], yColumn = columns[1] )

    def plotSaveAs_terminal( self, terminal, suffix = None, color = True ) :
        """Called when "File -> SaveAs eps" menu is selected. Uses a FileDialog to get the output 
        file name and produces an eps file of the current plot."""

        global saveAsCounter

        if( suffix is None ) : suffix = terminal
        fn = "Fudge%d.%s" % ( saveAsCounter, suffix )
        FileName = tkFileDialog.asksaveasfilename( initialfile = fn )
        if( ( type( FileName ) == type( "" ) ) and ( FileName != "" ) ) :
            if( color ) :
                self.g.hardcopy( filename = FileName, color = 1, terminal = terminal )
            else :
                self.g.hardcopy( filename = FileName, terminal = terminal )
            saveAsCounter += 1

    def plotSaveAsEps( self ) :
        """Called when "File -> SaveAs eps" menu is selected. Uses a FileDialog to get the output 
        file name and produces an eps file of the current plot."""

        self.plotSaveAs_terminal( 'postscript', suffix = 'ps' )

    def plotSaveAsPDF( self ) :
        """Called when "File -> SaveAs pdf" menu is selected. Uses a FileDialog to get the output 
        file name and produces an eps file of the current plot."""

        self.plotSaveAs_terminal( 'pdf' )

    def plotSaveAsPNG( self ) :
        """Called when "File -> SaveAs png" menu is selected. Uses a FileDialog to get the output 
        file name and produces an eps file of the current plot."""

        self.plotSaveAs_terminal( 'png', color = False )

    def plotSaveAsSVG( self ) :
        """Called when "File -> SaveAs svg" menu is selected. Uses a FileDialog to get the output 
        file name and produces an eps file of the current plot."""

        self.plotSaveAs_terminal( 'svg', color = False )

    def plotPrint( self ) :
        """Called when "File -> Print" menu is selected. Prints the current plot."""

        self.g.hardcopy( )

    def legendCallback( self ) :
        """Called whenever the "Options -> Legend" menu state is changed."""

        self.replot( )
