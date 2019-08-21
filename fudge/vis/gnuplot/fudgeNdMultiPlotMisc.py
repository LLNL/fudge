# <<BEGIN-copyright>>
# <<END-copyright>>

import Gnuplot
import Tkinter
import tkSimpleDialog
import tkFileDialog
import tkMessageBox

__metaclass__ = type

# Things to add line type, line width, symbol type, symbol width, symbol fill color

lineTypes = { 'none' : { 'value' : -2 }, 'system' : { 'value' : -1 }, 'solid' : { 'value' : 1 } }
defaultLineType = 'none'

colors = { 'system' : { 'value' : -1 }, 'red' : { 'value' : 1 }, 'green' : { 'value' : 2 }, 'blue' : { 'value' : 3 } }
defaultColor = 'system'

symbolTypes = { 'none' : { 'value' : -2 }, 'system' : { 'value' : -1 }, 'cross' : { 'value' : 1 }, 'x' : { 'value' : 2 }, 'square' : { 'value' : 4},
    'filled square' : { 'value' : 5 }, 'circle' : { 'value' : 6 }, 'filled circle' : { 'value' : 7 } }
defaultSymbolType = 'none'

class fudgeNdMultiPlotFileReadGetInfo( tkSimpleDialog.Dialog ):

#
# This class currently only works for dimension = 2.
#
    def __init__( self, parent, dimension, title, fileName, lines ) :

        self.values = None
        self.lineCount = len( lines )
        self.lines = ''.join( lines )
        self.fileName = fileName
        self.dimension = dimension
        tkSimpleDialog.Dialog.__init__( self, parent, title = title )

    def body( self, frame ) :

        Row = 0
        l = Tkinter.Label( frame, text = 'Please enter x, y column numbers' )
        l.grid( row = Row, column = 0, sticky = Tkinter.W )

        Row += 1
        cFrame = Tkinter.Frame( frame )
        cFrame.grid( row = Row, column = 0, sticky = Tkinter.W )

        cRow = 0
        l = Tkinter.Label( cFrame, text = 'xColumn = ' )
        l.grid( row = cRow, column = 0, sticky = Tkinter.W )
        self.xColumn = Tkinter.IntVar( )
        self.xColumn.set( 1 )
        ex = Tkinter.Entry( cFrame, textvariable = self.xColumn )
        ex.grid( row = cRow, column = 1, sticky = Tkinter.W )

        cRow += 1
        l = Tkinter.Label( cFrame, text = 'yColumn = ' )
        l.grid( row = cRow, column = 0, sticky = Tkinter.W )
        self.yColumn = Tkinter.IntVar( )
        self.yColumn.set( 2 )
        ey = Tkinter.Entry( cFrame, textvariable = self.yColumn )
        ey.grid( row = cRow, column = 1, sticky = Tkinter.W + Tkinter.E )

        cRow += 1
        l = Tkinter.Label( cFrame, text = 'label = ' )
        l.grid( row = cRow, column = 0, sticky = Tkinter.W )
        self.title = Tkinter.StringVar( )
        self.title.set( self.fileName )
        el = Tkinter.Entry( cFrame, textvariable = self.title )
        el.grid( row = cRow, column = 1, sticky = Tkinter.W + Tkinter.E )

        Row += 1
        l = Tkinter.Label( frame, text = 'First %d lines from file %s' % ( self.lineCount, self.fileName ), justify = Tkinter.LEFT, anchor = Tkinter.W )
        l.grid( row = Row, column = 0, sticky = Tkinter.W )

        Row += 1
        f = Tkinter.Frame( frame )
        f.grid( row = Row, column = 0 )
        t = Tkinter.Text( f )
        t.insert( Tkinter.END, self.lines )
        t.config( state = Tkinter.DISABLED )
        t.grid( row = Row, column = 0, sticky = Tkinter.W )

        frame.grid( )

    def validate( self ) :

        s = None
        try :
            x = self.xColumn.get( )
            try :
                y = self.yColumn.get( )
            except :
                s = 'y'
        except :
            s = 'x'
        if( s == None ) :
            if( x != y ) : 
                self.values = ( x, y, None, None )
                return( True )
            tkMessageBox.showerror( "Read file", "column numbers must differ" )
        else :
            tkMessageBox.showerror( "Read file", "%s column number not a valid integer" % s )
        return( False )

def fudgeNdMultiPlotFileRead( parent, dimension ) :

    fileName = tkFileDialog.askopenfilename( )
    if( type( fileName ) != type( "" ) ) : return( None, None, None )
    try :
        f = open( fileName ) 
    except :
        return( None, None, None )
    ls = []
    for i in xrange( 20 ) :
        l = f.readline( )
        if( l == '' ) : break
        ls.append( l )
    f.close( )

    fileInfo = fudgeNdMultiPlotFileReadGetInfo( parent, 2, 'input for 2d file', fileName, ls )
    return( fileName, fileInfo.title.get( ), fileInfo.values )

class fudgeNdMultiPlotItem :

    def __init__( self, title, dimension, active = True, replot = None, xColumn = 1, yColumn = 2, zColumn = 3, tColumn = None ) :

        self.parent = None
        self.fileName = None

        self.title = title
        self.dimension = dimension
        self.active = Tkinter.IntVar( )
        self.active.set( active )
        self.replot = replot
        self.xColumn = xColumn
        self.yColumn = yColumn
        self.zColumn = zColumn
        self.tColumn = tColumn

        self.gnuPlotItem = None
        self.color = 'system'
        self.lineType = 'system'
        self.lineWidth = 'system'
        self.symbolType = 'system' 
        self.symbolSize = 'system'

    def getXColumn( self ) :

        return( self.xColumn )

    def setXColumn( self, xColumn ) :

        self.xColumn = xColumn

    def getYColumn( self ) :

        return( self.yColumn )

    def setYColumn( self, yColumn ) :

        self.yColumn = zColumn

    def getZColumn( self ) :

        return( self.zColumn )

    def setZColumn( self, zColumn ) :

        self.zColumn = zColumn

    def getTitle( self ) :

        return( self.title.get( ) )

    def setTitle( self, title ) :

        self.title.set( title )

    def setActive( self, active ) :

        self.active.set( active )

    def isActive( self ) :

        return( self.active.get( ) )

    def setColor( self, color ) :

        self.color = color
        if( self.replot != None ) : self.replot( )

    def getColor( self ) :

        return( self.color )

    def setToUpdate( self ) :

        self.parent.setUpdatingPlotParameters( self )

class fudgeNdMultiPlotSystem( fudgeNdMultiPlotItem ) :

    def __init__( self ) :

        fudgeNdMultiPlotItem.__init__( self, 'system', 0, active = False )
        self.lineType = 'solid'
        self.lineWidth = 1
        self.symbolType = 'none'
        self.symbolSize = 2

class fudgeNdMultiPlotData( fudgeNdMultiPlotItem ) :

    def __init__( self, data, title, dimension = None, active = True, replot = None ) :

        if( len( data ) == 0 ) : raise Exception( 'data with title = %s is empty' % title )
        dimension_ = 1
        if( ( type( data[0] ) == type( [] ) ) or ( type( data[0] ) == type( () ) ) ) : dimension_ = len( data[0] )
        if( dimension != dimension_ ) : raise Exception( "dimension = %d does not agree with data's dimension of %d for title = %s" % \
            ( dimension, dimension_, title ) )
        fudgeNdMultiPlotItem.__init__( self, title, dimension, active = active )
        self.gnuPlotItem = Gnuplot.Data( data, title = title )

class fudgeNdMultiPlotFile( fudgeNdMultiPlotItem ) :

    def __init__( self, fileName, dimension, title = "", active = True, replot = None, xColumn = 1, yColumn = 2, zColumn = 3, tColumn = None ) :

        if( title == "" ) : title = fileName
        fudgeNdMultiPlotItem.__init__( self, title, dimension, active = active, xColumn = xColumn, yColumn = yColumn, zColumn = zColumn, tColumn = yColumn )
        self.gnuPlotItem = fileName
        self.gnuPlotItem = Gnuplot.File( fileName, title = title )

class fudgeNdMultiPlotItemsDialog :

    def __init__( self, frame, redraw, activeMenubar = None ) :

        self.redraw = redraw
        self.frame = frame
        self.systemPlotParameters = fudgeNdMultiPlotSystem( )
        self.nextPlotNumber = 0
        self.updatingPlot = None
        self.fudgeNdMultiPlotFile = []
        self.activeMenubar = activeMenubar

        if( self.activeMenubar != None ) :
            self.activeMenu = Tkinter.Menu( self.activeMenubar, tearoff = 1 )
            self.activeMenu.add_command( label = "all on", command = self.activeAllOnCallback )
            self.activeMenu.add_command( label = "all off", command = self.activeAllOffCallback )
            self.activeMenu.add_command( label = "step", command = self.activeStepCallback )
            self.activeMenubar.add_cascade( label = "Active", menu = self.activeMenu )

        # Define containment frame for all dataset options
        setOptionsFrame = Tkinter.LabelFrame( frame, text="Dataset Options" )
        setOptionsFrame.grid( row=0, column=0, sticky = Tkinter.S )

        # Define header frame of dataset options
        headerFrame = Tkinter.Frame( setOptionsFrame )
        headerFrame.grid( row = 0, column = 0, sticky = Tkinter.W + Tkinter.E )

        # Set selection box
        self.selectedMenuButton = Tkinter.Menubutton( headerFrame, text = 'Pick a set', relief = Tkinter.RAISED )
        self.selectedMenuButton.grid( row = 0, column = 0, sticky = Tkinter.W )
        self.selectedMenu = Tkinter.Menu( self.selectedMenuButton, tearoff = False )
        self.selectedMenuButton['menu'] = self.selectedMenu
        self.titleLabel = self.addLabel( headerFrame, text = '   Set title: ', row = 0, column = 1 )
        self.titleEntry, self.titleVar = self.addStringEntry( headerFrame, '', row = 0, column = 2, sticky = Tkinter.W + Tkinter.E )
        self.titleEntry['width'] = 40

        # Define body frame for dataset options
        brow = 0
        bcolumn = 0
        bodyFrame = Tkinter.Frame( setOptionsFrame )
        bodyFrame.grid( row = 1, column = 0, sticky = Tkinter.W + Tkinter.E )

        # Remove set button        
        self.removeButton = Tkinter.Button( bodyFrame,  text = 'Remove\nset?', command = self.removeCallback )
        self.removeButton.grid( row = brow, column = bcolumn )

        # Set color box
        bcolumn += 1
        colors['system']['command'] = self.setColorSystem
        colors['red']['command'] = self.setColorRed
        colors['green']['command'] = self.setColorGreen
        colors['blue']['command'] = self.setColorBlue
        colorFrame = Tkinter.LabelFrame( bodyFrame, text="Color" )
        colorFrame.grid( row = brow, column = bcolumn, rowspan = 3 )
        #self.lineColor = Tkinter.StringVar()
        #self.lineColor.set( 'red' )
        #self.colorMenu = Tkinter.OptionMenu( colorFrame, self.lineColor, *colors.keys() )
        #self.colorMenu.grid( row = 0, column = 0 )
        self.lineColor, self.colorMenu = self.addButtonWithMenu( colorFrame, text = 'Color', row = 0, column = 0, options = colors, default = 'red' )
        self.colorLabel = self.addLabel( colorFrame, text = 'system', row = 1, column = 0 )
        self.colorLabel['width'] = len( 'system' )

        # Set linestyle box
        bcolumn += 1
        lineFrame = Tkinter.LabelFrame( bodyFrame, text = "Linestyle" )
        lineFrame['borderwidth'] = 2
        lineFrame['relief'] = Tkinter.GROOVE
        lineFrame.grid( row = brow, column = bcolumn, rowspan = 3 )
        lineTypes['none']['command'] = self.setLineTypeNone
        lineTypes['system']['command'] = self.setLineTypeSystem
        lineTypes['solid']['command'] = self.setLineTypeSolid
        self.lineType, self.lineTypeMenu =\
            self.addButtonWithMenu( lineFrame, text = 'Type', row = 0, column = 0, options = lineTypes, default = 'system' )
        self.lineTypeLabel = self.addLabel( lineFrame, text = 'system', row = 1, column = 0 )
        self.lineTypeLabel['width'] = len( 'system' )
        self.lineWidth = self.addLabel( lineFrame, text = 'Width', row = 0, column = 1 )
        self.lineWidthVar = Tkinter.IntVar( )
        self.lineWidthSpinBox = Tkinter.Spinbox( lineFrame, from_ = -1, to = 50, increment = 1, width = 5, textvariable = self.lineWidthVar )
        self.lineWidthSpinBox.grid( row = 1, column = 1 )

        # Set symbol box
        bcolumn += 1
        symbolFrame = Tkinter.LabelFrame( bodyFrame, text="Symbol" )
        symbolFrame['borderwidth'] = 2
        symbolFrame['relief'] = Tkinter.GROOVE
        symbolFrame.grid( row = brow, column = bcolumn, rowspan = 3 )
        symbolTypes['none']['command'] = self.setSymbolTypeNone
        symbolTypes['system']['command'] = self.setSymbolTypeSystem
        symbolTypes['cross']['command'] = self.setSymbolTypeCross
        symbolTypes['x']['command'] = self.setSymbolX
        symbolTypes['square']['command'] = self.setSymbolSquare
        symbolTypes['filled square']['command'] = self.setSymbolFilledSquare
        symbolTypes['circle']['command'] = self.setSymbolCircle
        symbolTypes['filled circle']['command'] = self.setSymbolFilledCircle
        self.symbolType, self.symbolTypeMenu =\
            self.addButtonWithMenu( symbolFrame, text = 'Type', row = 0, column = 0, options = symbolTypes, default = 'none' )
        self.symbolTypeLabel = self.addLabel( symbolFrame, text = 'none', row = 1, column = 0 )
        self.symbolTypeLabel['width'] = len( 'filled circle' )
        self.symbolSize = self.addLabel( symbolFrame, text = 'Size', row = 0, column = 1 )
        self.symbolSizeVar = Tkinter.IntVar( )
        self.symbolSizeSpinBox = Tkinter.Spinbox( symbolFrame, from_ = -1, to = 50, increment = 1, width = 5, textvariable = self.symbolSizeVar )
        self.symbolSizeSpinBox.grid( row = 1, column = 1 )

        # Update and redraw buttons, these live outside the dataset options frame
        updateFrame = Tkinter.Frame( frame )
        updateFrame['borderwidth'] = 2
        updateFrame['relief'] = Tkinter.GROOVE
        updateFrame.grid( row = 1, column = 0, sticky = Tkinter.E + Tkinter.W )
        self.UpdateAndRedraw = Tkinter.Button( updateFrame, text = 'Update and redraw', command = self.updateAndRedrawCallback  )
        self.UpdateAndRedraw.grid( row = 0, column = 1, sticky = Tkinter.E )
        self.UpdateOnly = Tkinter.Button( updateFrame, text = 'Update only', command = self.updateOnlyCallback  )
        self.UpdateOnly.grid( row = 0, column = 0, sticky = Tkinter.W )
        
        # Final configuration
        self.frame.columnconfigure( bcolumn, weight = 10 )
        self.frame.pack( )
        self.addFudgeNdMultiPlotFile( self.systemPlotParameters )
        self.setUpdatingPlotParameters( self.systemPlotParameters )

    def addLabel( self, frame, text, row, column, columnspan = 1 ) :

        l = Tkinter.Label( frame, text = text )
        l.grid( row = row, column = column, columnspan = columnspan )
        return( l )

    def addButtonWithMenu( self, frame, text, row, column, options, default ) :

        mb = Tkinter.Menubutton( frame, text = text, relief = Tkinter.RAISED )
        mb.grid( row = row, column = column )
        m = self.addMenu( mb, options, default )
        mb['menu'] = m
        return( mb, m )

    def addMenu( self, menuButton, list, default ) :

        m = Tkinter.Menu( menuButton, tearoff = False )
        newlist = []
        for key in list : newlist.append( [ list[key]['value'] , key ] )
        newlist.sort( )
        for value, key in newlist : m.add_command( label = key, command = list[key]['command'] )
        return( m )

    def addIntEntry( self, frame, value, row, column, width = 5 ) :

        v = Tkinter.IntVar( )
        v.set( value )
        e = Tkinter.Entry( frame, textvariable = v, width = width )
        e.grid( row = row, column = column )
        return( e, v )

    def addStringEntry( self, frame, string, row, column, columnspan = 1, sticky = None ) :

        v = Tkinter.StringVar( )
        v.set( string )
        e = Tkinter.Entry( frame, textvariable = v )
        e.grid( row = row, column = column, columnspan = columnspan, sticky = sticky )
        return( e, v )

    def addFudgeNdMultiPlotFile( self, gnuPlotItem ) :

        gnuPlotItem.plotNumber = self.nextPlotNumber
        self.nextPlotNumber += 1
        self.fudgeNdMultiPlotFile.append( gnuPlotItem )
        self.selectedMenu.add_command( label = gnuPlotItem.title, command = gnuPlotItem.setToUpdate )
        gnuPlotItem.parent = self
        if( ( self.activeMenubar != None ) and ( gnuPlotItem.gnuPlotItem != None ) ) :
            self.activeMenu.add_checkbutton( label = gnuPlotItem.title, variable = gnuPlotItem.active, command = self.activeChanged )
            self.updateColumnBreaks( )

    def updateColumnBreaks( self ) :

        n = self.activeMenu.index( Tkinter.END )
        for i in xrange( 3, n ) :
            c = 0
            if( ( ( i - 1 ) % 30 ) == 0 ) : c = 1
            self.activeMenu.entryconfigure( i, columnbreak = c )
        n = self.selectedMenu.index( Tkinter.END )
        for i in xrange( 3, n ) :
            c = 0
            if( ( ( i - 1 ) % 30 ) == 0 ) : c = 1
            self.selectedMenu.entryconfigure( i, columnbreak = c )
            
#
# Color
#
    def setColorSystem( self ) :

        self.colorLabel['text'] = 'system'

    def setColorRed( self ) :

        self.colorLabel['text'] = 'red'

    def setColorGreen( self ) :

        self.colorLabel['text'] = 'green'

    def setColorBlue( self ) :

        self.colorLabel['text'] = 'blue'
#
# line type
#
    def setLineType( self, type ) :

        self.lineTypeLabel['text'] = type

    def setLineTypeNone( self ) :

        self.setLineType( 'none' )

    def setLineTypeSystem( self ) :

        self.setLineType( 'system' )

    def setLineTypeSolid( self ) :

        self.setLineType( 'solid' )
#
# symbol type
#
    def setSymbolType( self, type ) :

        self.symbolTypeLabel['text'] = type

    def setSymbolTypeNone( self ) :

        self.setSymbolType( 'none' )

    def setSymbolTypeSystem( self ) :

        self.setSymbolType( 'system' )

    def setSymbolTypeCross( self ) :

        self.setSymbolType( 'cross' )

    def setSymbolX( self ) : 

        self.setSymbolType( 'x' )

    def setSymbolSquare( self ) : 

        self.setSymbolType( 'square' )

    def setSymbolFilledSquare( self ) : 

        self.setSymbolType( 'filled square' )

    def setSymbolCircle( self ) : 

        self.setSymbolType( 'circle' )

    def setSymbolFilledCircle( self ) : 

        self.setSymbolType( 'filled circle' )
#
#
#
    def setUpdatingPlotParameters( self, gnuPlotItem ) :

        def setMenuItemState( menu, label, state ) :

            i = 0
            while( True ) :
                index = menu.index( i )
                if( index == None ) : break
                if( menu.entrycget( index, 'label' ) == label ) :
                    menu.entryconfigure( index, state = state )
                    break
                i += 1
                if( i > 20 ) : break                            # Just in case.

        removeState = Tkinter.NORMAL
        self.gnuPlotItem = gnuPlotItem
        if( gnuPlotItem.gnuPlotItem == None ) : removeState = Tkinter.DISABLED
        self.removeButton['state'] = removeState
        setMenuItemState( self.lineTypeMenu, 'system', removeState )
        self.titleEntry['state'] = removeState
        self.titleVar.set( gnuPlotItem.title )
        self.colorLabel['text'] = gnuPlotItem.color
        self.lineTypeLabel['text'] = gnuPlotItem.lineType
        lineWidth = gnuPlotItem.lineWidth
        if( lineWidth == 'system' ) : lineWidth = -1
        self.lineWidthVar.set( lineWidth )
        self.symbolTypeLabel['text'] = gnuPlotItem.symbolType
        symbolSize = gnuPlotItem.symbolSize
        if( symbolSize == 'system' ) : symbolSize = -1
        self.symbolSizeVar.set( symbolSize )
        self.updatingPlot = gnuPlotItem

    def updateAndRedrawCallback( self ) :

        self.updateOnlyCallback( )        
        self.redraw( )

    def updateOnlyCallback( self ) :

        title = self.titleVar.get( )
        if( title != self.gnuPlotItem.title ) :
            index, listIndex, activeIndex = self.getIndicesForLabel( self.gnuPlotItem.title )
            self.gnuPlotItem.title = title
            self.selectedMenu.entryconfigure( listIndex, label = title )
            self.activeMenu.entryconfigure( activeIndex, label = title )
        self.gnuPlotItem.color = self.colorLabel['text']
        self.gnuPlotItem.lineType = self.lineTypeLabel['text']
        self.gnuPlotItem.lineWidth = self.lineWidthVar.get( )
        self.gnuPlotItem.symbolType = self.symbolTypeLabel['text']
        self.gnuPlotItem.symbolSize = self.symbolSizeVar.get( )

    def removeCallback( self ) :

        index, listIndex, activeIndex = self.getIndicesForLabel( self.gnuPlotItem.title )
        if( index != None ) :
            self.selectedMenu.delete( listIndex )
            if( activeIndex != None ) : self.activeMenu.delete( activeIndex )
            self.setUpdatingPlotParameters( self.fudgeNdMultiPlotFile[index-1] )
            del self.fudgeNdMultiPlotFile[index]
            self.updateColumnBreaks( )
            self.redraw( )

    def getIndicesForLabel( self, label ) :

        n = len( self.fudgeNdMultiPlotFile )
        index = None
        listIndex = None
        activeIndex = None
        for i in xrange( 1, n ) :
            j = self.selectedMenu.index( i )
            if( self.selectedMenu.entrycget( j, 'label' ) == label ) :
                index = i
                listIndex = j
                break
        if( self.activeMenubar != None ) :
            for i in xrange( 3, n + 3 ) :
                j = self.activeMenu.index( i )
                if( self.activeMenu.entrycget( j, 'label' ) == label ) :
                    activeIndex = j
                    break
        return( index, listIndex, activeIndex )

    def redraw( self ) :

        self.redraw( )

    def activeChanged( self ) :

        self.redraw( )

    def getActiveCurves( self ) :

        activeCurves = []
        n = len( self.fudgeNdMultiPlotFile )
        systemGnuPlotItem = self.fudgeNdMultiPlotFile[0]
        for gnuPlotItem in self.fudgeNdMultiPlotFile :
            if( gnuPlotItem.gnuPlotItem == None ) :
                continue
            elif( gnuPlotItem.isActive( ) ) :
                gnuPlotItem.gnuPlotItem.set_option( title = gnuPlotItem.title )
                gnuPlotItem.gnuPlotItem.set_option( using = "%d:%d" % ( gnuPlotItem.xColumn, gnuPlotItem.yColumn ) )
                symbolType = gnuPlotItem.symbolType
                symbolSize = gnuPlotItem.symbolSize
                if( symbolType == 'system' ) :
                    symbolType = systemGnuPlotItem.symbolType
                    symbolSize = systemGnuPlotItem.symbolSize
                if( symbolSize < 1 ) : symbolSize = systemGnuPlotItem.symbolSize
                symbolTypeValue = symbolTypes[symbolType]['value']
                lineType = gnuPlotItem.lineType
                lineWidth = gnuPlotItem.lineWidth
                if( lineType == 'system' ) :
                    lineType = systemGnuPlotItem.lineType
                    lineWidth = systemGnuPlotItem.lineWidth
                if( lineWidth < 1 ) : lineWidth = systemGnuPlotItem.lineWidth
                lineTypeValue = lineTypes[lineType]['value']
                color = gnuPlotItem.color
                colorSpecs = ''
                if( color == 'system' ) : color = systemGnuPlotItem.color
                if( color != 'system' ) : colorSpecs = ' lt %d' % colors[color]['value']
                linespoints = ''
                lineSpecs = ''
                if( lineTypeValue > 0 ) :
                    linespoints = 'lines'
                    lineSpecs = ' lw %d' % lineWidth
                pointsSpecs = ''
                if( symbolTypeValue > 0 ) :
                    linespoints += 'points'
                    pointsSpecs = ' pt %d ps %d' % ( symbolTypeValue, symbolSize )
                elif( symbolType == 'system' ) :
                    linespoints += 'points'
                    pointsSpecs = ' ps %d' % symbolSize
                if( linespoints != '' ) :
                    _with = linespoints + colorSpecs + lineSpecs + pointsSpecs
                    gnuPlotItem.gnuPlotItem._options['with'] = ( _with, 'with ' + _with )
#                    gnuPlotItem.gnuPlotItem.set_option( with = linespoints + colorSpecs + lineSpecs + pointsSpecs )
                    activeCurves.append( gnuPlotItem.gnuPlotItem )
        return( activeCurves )

    def activeAllOnCallback( self ) :
        """Called when the "Active -> all on" menu item is selected. Sets all file's active state to on."""

        self.activeCallback( True )

    def activeAllOffCallback( self ) :
        """Called when the "Active -> all off" menu item is selected. Sets all file's active state to off."""

        self.activeCallback( False )

    def activeCallback( self, state ) :
        """Called from activeAllOnCallback and activeAllOffCallback."""

        for gnuPlotItem in self.fudgeNdMultiPlotFile :
            if( gnuPlotItem.gnuPlotItem == None ) : continue
            gnuPlotItem.active.set( state )
        self.redraw( )

    def activeStepCallback( self ) :
        """Called when the "Active -> step" menu item is selected."""

        index = 1
        counter = 0
        i = 1
        for gnuPlotItem in self.fudgeNdMultiPlotFile :
            if( gnuPlotItem.gnuPlotItem == None ) : continue
            if( gnuPlotItem.isActive( ) ) :
                index = i
                counter += 1
            i += 1
        if( counter == 1 ) :
            self.fudgeNdMultiPlotFile[index].setActive( False )
            index += 1
            if( index >= len( self.fudgeNdMultiPlotFile ) ) : index = 1
        elif( counter > 1 ) :
            index = 1
            for gnuPlotItem in self.fudgeNdMultiPlotFile : gnuPlotItem.setActive( False )
        if( len( self.fudgeNdMultiPlotFile ) > 1 ) : self.fudgeNdMultiPlotFile[index].setActive( True )
        self.redraw( )
