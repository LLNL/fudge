# <<BEGIN-copyright>>
# <<END-copyright>>

from fudge.core.utilities.brb import uniquify
from pqu.physicalQuantityWithUncertainty import PhysicalQuantityWithUncertainty as PQU
from pqu.physicalQuantityWithUncertainty import valueOrPQ
from fudge.core.math.endl2dmathClasses import endl2dmath
from fudge.core.math.endl3dmathClasses import endl3dmath
from fudge.core.math.xData.XYs import XYs
from fudge.core.math.xData.W_XYs import W_XYs

__metaclass__ = type

__docformat__ = "restructuredtext en"

INCLASSAXISSETTING = False


#--------------------------------------------------------
#
# global data
#
#--------------------------------------------------------
defaultSymbols = [ 'o', 'v', '^', '<', '>', '1', '2', '3', '4', 's', 'p', '*', 'h', 
                   'H', '+', 'x', 'D', 'd', '|', '_' ]
nSymbols = len( defaultSymbols )
    
defaultColors = [ 'k','r','g','b','y','brown','gray','violet','cyan',
                  'magenta','orange','indigo','maroon','turquoise',
                  'darkgreen','0.20','0.40','0.60','0.80' ]
nColors = len( defaultColors )

#--------------------------------------------------------
#
# Classes for new, matplotlib based plotting
#
#--------------------------------------------------------
class TextSettings:
    '''
    Base class for text objects in plots
    
    Member data (set these with __init__ function):
    
    * **text** the text itself
    * **font** the font (e.g. Helvetica) (optional)
    * **fontSize** font size in points (optional, defaults to 14)
    * **fontWeight** font weight (e.g. bold) (optional)
    
    '''
    def __init__(self, text=None, font=None, fontSize=12, fontWeight='bold' ):
        self.text = text
        self.font = font
        self.fontSize = fontSize
        self.fontWeight = fontWeight


	
class AxisSettings( TextSettings ):
    '''
    Member data (set these with __init__ function) ::
    
        - `axisMin` :  min value of the axis 
            This is optional and defaults to `None` causing autoscaling.  
            If given as float, we assume in units of `unit`.  
            If given as a physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty instance, we'll 
            attempt to convert it to `unit`.
        - `axisMax` : max value of the axis.  See axisMin for usage
        - `isLog` : (optional, defaults to False) False == linear, True == semilog
        - `unit` : (optional, defaults to None), must be None or something that is usable in a 
             physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty instance
        - `gridOn` : (optional)
    
    '''
    def __init__( self, axisMin=None, axisMax=None, autoscale=True, isLog=False, unit=None,
            gridOn=False, label=None, font=None, fontSize=20, fontWeight='normal' ):
        TextSettings.__init__(self, text=label, font=font, fontSize=fontSize, fontWeight=fontWeight )
        if isinstance( axisMin, PQU ): 
            if unit!=None: self.axisMin = axisMin.getValueAs( unit )
            else: self.axisMin = axisMin.getValue(  )
        else: self.axisMin = axisMin
        if isinstance( axisMax, PQU ): 
            if unit!=None: self.axisMax = axisMax.getValueAs( unit )
            else: self.axisMax = axisMax.getValue(  )
        else: self.axisMax = axisMax
        self.isLog = isLog
        self.gridOn = gridOn
        self.autoscale = autoscale
        if axisMin is not None and axisMax is not None:
            self.autoscale = False
        self.unit = unit
        self.label = self.text
        if self.unit != None:
            if self.unit == '': self.label += ' (dimensionless)'
            else: self.label += ' (' + self.unit + ')'
            
    def setupAxis( self, theAxis ):
        theAxis.grid( self.gridOn, linewidth=1.5 )
        if INCLASSAXISSETTING: # this doesn't' work, should I remove it?
            if self.isLog: theAxis.set_scale( 'log', linewidth=2 )
            else:          theAxis.set_scale( 'linear', linewidth=2 )
            if not self.autoscale: 
                theAxis.set_data_interval( self.axisMin, self.axisMax )
                theAxis.set_view_interval( self.axisMin, self.axisMax )


class DataSetParameters:
    '''
    Base class for all DataSet classes (DataSet2d, DataSet3d, ...)
    
    Member data (set these with __init__ function):
    
    * legend (optional, defaults to '')
    * symbol
    * symbolColor
    * symbolSize
    * lineStyle
    * lineColor
    * lineWidth

    Commands which take color arguments can use several formats to specify the colors. 
    For the basic builtin colors, you can use a single letter
    
      =========   ====================
      character   description
      =========   ====================
              b   blue
              g   green
              r   red
              c   cyan
              m   magenta
              y   yellow
              k   black
              w   white
      =========   ====================
    
    Gray shades can be given as a string encoding a float in the 0-1 range, e.g.:
    
    >>> color = '0.75'
    
    For a greater range of colors, you have two options. You can specify the color 
    using an html hex string, as in:
    
    >>> color = '#eeefff'
    
    or you can pass an R , G , B tuple, where each of R , G , B are in the range [0,1].
    Finally, legal html names for colors, like "red", "burlywood" and "chartreuse" are 
    supported.

    The following format string characters are accepted to control the line style or 
    marker:
    
        =========   ====================
        character   description
        =========   ====================
        '-'         solid line style
        '--'        dashed line style
        '-.'        dash-dot line style
        ':'         dotted line style
        '.'         point marker
        ','         pixel marker
        'o'         circle marker
        'v'         triangle_down marker
        '^'         triangle_up marker
        '<'         triangle_left marker
        '>'         triangle_right marker
        '1'         tri_down marker
        '2'         tri_up marker
        '3'         tri_left marker
        '4'         tri_right marker
        's'         square marker
        'p'         pentagon marker
        '*'         star marker
        'h'         hexagon1 marker
        'H'         hexagon2 marker
        '+'         plus marker
        'x'         x marker
        'D'         diamond marker
        'd'         thin_diamond marker
        '|'         vline marker
        '_'         hline marker
        =========   ====================
    '''
    def __init__( self, legend='', symbol=None, symbolSize=5, color=None, lineStyle=None, 
            lineWidth=None, errorbarCapSize=3, errorbarLineWidth=1.5, errorbarColor=None, 
            formatString=None, dataType=None ):
        self.legend            = legend
        self.symbol            = symbol
        self.color             = color
        self.symbolSize        = symbolSize # not yet used
        self.lineStyle         = lineStyle
        self.lineWidth         = lineWidth
        self.errorbarCapSize   = errorbarCapSize
        self.errorbarLineWidth = errorbarLineWidth
        self.errorbarColor     = errorbarColor
        self.formatString      = formatString
        
    def getFormatString( self ): 
        '''Attempt to make simple format string for maptplotlib use.'''
        if self.formatString != None: return self.formatString
        r = ''
        if self.color != None and len( self.color ) == 1: r += self.color
        if self.symbol and len( self.symbol ) == 1 != None: r+= self.symbol
        if self.lineStyle and len( self.lineStyle ) == 1 != None: r+= self.lineStyle
        return r
        
    def getFormatMap( self, doErrors = False ):
        def updateMap( theMap, key, theValue, defaultValue ):
            if theValue != None: theMap[ key ] = theValue
            elif defaultValue != None: theMap[ key ] = defaultValue
        kw = {}
        updateMap( kw, 'color', self.color, None )
        updateMap( kw, 'marker', self.symbol, None )
        if self.symbol in [ '+', 'x', '1', '2', '3', '4', '|', '_' ]: kw ['markeredgewidth'] = 2 
        else: kw ['markeredgewidth'] = 1
        updateMap( kw, 'markersize', self.symbolSize, 5 )
        updateMap( kw, 'linestyle', self.lineStyle, '-' )
        updateMap( kw, 'linewidth', self.lineWidth, 1 )
        updateMap( kw, 'label', self.legend, None )
        if doErrors:
            updateMap( kw, 'elinewidth', self.errorbarLineWidth, 1.5 )
            updateMap( kw, 'capsize', self.errorbarCapSize, 3 )
            updateMap( kw, 'ecolor', self.errorbarColor, self.color )
        return kw

    def addToPlot( self, thePlot ): raise NotImplementedError( "Implement in derived classes" )


class DataSet2d( DataSetParameters ):
    '''
    Member data (set these with __init__ function):
    
    * data: the data as either a list of [ x, y ] pairs or as fudge XYs of legacy 
      endl2dmath objects
    * uncertainty: the uncertainty on the data as either a list of [ x, y ] pairs or as 
      fudge XYs of legacy endl2dmath objects
    * dataType ('scatter', 'line', or None): by default endl2dmath or fudge XYs objects 
      are rendered as 'line' and lists of points as 'scatter'.  
      Override the default behavior with this keyword.
    * xUnits: Use this option to set the x units because neither endl2dmath 
      objects nor lists of pairs have units associated with them (on the other hand, fudge 
      XYs do have x and y units associated with them).  The default is None.  
      When the units are given as None and the data/uncertainty are not fudge XYs, the 
      units are assumed to match those of the plot axis.
    * yUnits: Use this option to set the y units because neither endl2dmath 
      objects nor lists of pairs have units associated with them (on the other hand, fudge 
      XYs do have x and y units associated with them).  The default is None.  
      When the units are given as None and the data/uncertainty are not fudge XYs, the 
      units are assumed to match those of the plot axis.
    
    type help(DataSetParameters) to find out about how to set the line properties, etc.  
    All of this information is passed on though the __init__ function of DataSet2d to 
    DataSetParameters.
    '''
    
    def __init__( self, data, uncertainty = None, xUnit = None, yUnit = None, dataType = None, **kw ):
        
        self.xerr = None
        self.yerr = None
        self.xUnit = xUnit
        self.yUnit = yUnit

        if type(data) == list:
            self.x = [ p[0] for p in data ]
            self.y = [ p[1] for p in data ]
            self.dataType = 'scatter'
            if 'lineStyle' not in kw: kw[ 'lineStyle' ] = '-'
            if type(uncertainty) == list:
                try: self.xerr = [ p[0] for p in uncertainty ]
                except IndexError: pass
                try: self.yerr = [ p[1] for p in uncertainty ]
                except IndexError: pass
        elif isinstance( data, endl2dmath ):
            data = data.toInterpolation( 0, 0.01 )
            self.x = [ p[0] for p in data.data ]
            if True:
                self.x.sort()
                self.y = [ data.getValue( xx ) for xx in self.x ]
            else:
                self.y = [ p[1] for p in data.data ]
            self.dataType = 'line'
            if isinstance( uncertainty, endl2dmath ): self.yerr = [ uncertainty.getValue( xx ) for xx in self.x ]
        elif isinstance( data, XYs ):
            self.x, self.y = data.toPointwiseLinear( 0, 0.00001 ).copyDataToXsAndYs()
            self.xUnit = data.axes[0].getUnit()
            self.yUnit = data.axes[1].getUnit()
            self.dataType = 'line'
            if isinstance( uncertainty, XYs ): self.yerr = [ uncertainty.getValue( xx ) for xx in self.x ]
        else: raise TypeError( "Unsupported type: " + str( type( data ) ) )
        
        if dataType != None: self.dataType = dataType
        DataSetParameters.__init__( self, **kw )
    
    def convertUnits( self, xUnit = None, yUnit = None ):
        '''Convert self's units to xUnit and yUnit'''
        if xUnit != None and self.xUnit != None and self.xUnit != xUnit:
            try:  conversionFactor = valueOrPQ( 1.0, unitFrom=self.xUnit, unitTo=xUnit, asPQU=False )
            except TypeError, err: 
                raise TypeError( err.message + ': '+self.xUnit+' -> '+xUnit+' for set '+str( self.legend ) )
            except NameError, err:
                raise NameError( err.message + ': '+self.xUnit+' -> '+xUnit+' for set '+str( self.legend ) )
            if self.xerr != None: self.xerr = [ xx*conversionFactor for xx in self.xerr ]
            self.x = [ xx*conversionFactor for xx in self.x ]
            self.xUnit = xUnit
        if yUnit != None and self.yUnit != None and self.yUnit != yUnit:        
            try:  conversionFactor = valueOrPQ( 1.0, unitFrom=self.yUnit, unitTo=yUnit, asPQU=False )
            except TypeError, err: 
                raise TypeError( err.message + ': '+self.yUnit+' -> '+yUnit+' for set '+str( self.legend ) )
            except NameError, err:
                raise NameError( err.message + ': '+self.yUnit+' -> '+yUnit+' for set '+str( self.legend ) )
            if self.yerr != None: self.yerr = [ yy*conversionFactor for yy in self.yerr ]
            self.y = [ yy*conversionFactor for yy in self.y ]
            self.yUnit = yUnit
        
    def addToPlot( self, thePlot, logY = False, minY = 1e-12, verbose=False ):      
        if self.dataType == 'scatter':
            thePlot.errorbar( self.x, self.y, yerr=self.yerr, xerr=self.xerr, 
                **self.getFormatMap( doErrors = self.xerr != None or self.yerr != None ) )
        elif self.dataType == 'line':
            nRange = range( len( self.y ) )
            if not logY:  
                thePlot.plot( self.x, self.y, **self.getFormatMap(  ) )
                if self.yerr != None:
                    lowerbound = [ self.y[i] - abs( self.yerr[i] ) for i in nRange ]
                    upperbound = [ self.y[i] + abs( self.yerr[i] ) for i in nRange ]
                    thePlot.fill_between( self.x, lowerbound, upperbound, facecolor = self.color, alpha = 0.25 )    
            else:
                theYs = []
                if self.yerr != None: 
                    theYLowBounds = []
                    theYUpBounds = []
                for i in range( len( self.y ) ):
                    if self.y[i] <= 0.0: 
                        if verbose: print "WARNING: y value <= 0.0 (",self.y[i],") at index",i,"on a log plot"
                        theYs.append( minY )
                    else: theYs.append( self.y[i] )
                    if self.yerr != None:
                        lb = self.y[i] - abs( self.yerr[i] )
                        if lb <= 0.0: 
                            if verbose: print "WARNING: y-yerr value <= 0.0 (",lb,") at index",i,"on a log plot"
                            theYLowBounds.append( minY )
                        else: theYLowBounds.append( lb )
                        ub = self.y[i] + abs( self.yerr[i] )
                        if ub <= 0.0: 
                            if verbose: print "WARNING: y+yerr value <= 0.0 (",ub,") at index",i,"on a log plot"
                            theYUpBounds.append( minY )
                        else: theYUpBounds.append( ub )
                thePlot.plot( self.x, theYs, **self.getFormatMap(  ) )
                if self.yerr != None: thePlot.fill_between( self.x, theYLowBounds, theYUpBounds, facecolor = self.color, alpha = 0.25 )    
        else: raise TypeError( 'dataType must be "scatter" or "line", found ' + str( self.dataType ) )
        
        

class DataSet3d( DataSetParameters ):
    '''
    Member data (set these with __init__ function):
    
    * data
    * dataType (None)
    
    type help(DataSetParameters) to find out about how to set the line properties, etc.  
    All of this information is passed on though the __init__ function of DataSet2d to 
    DataSetParameters.
    '''
    
    def __init__( self, data, dataType = None, **kw ):
        if isinstance( data, endl3dmath ):
            self.x = [ p[0] for p in data.data ]
            tmp = []
            for p in data.data: tmp += [ pp[0] for pp in p[1] ]
            self.y = sorted( uniquify( tmp ) )
            self.z = [ ]
            for yy in self.y:
                self.z.append( [ data.getValue( xx, yy ) for xx in self.x ] )
            self.dataType = 'line'
        elif isinstance( data, W_XYs ):
            newData = data.toPointwiseLinear( 0, 0 )
            self.x, self.y, self.z = newData.copyDataToGridWsAndXsAndYs()
            self.xUnit = data.axes[0].getUnit()
            self.yUnit = data.axes[1].getUnit()
            self.zUnit = data.axes[2].getUnit()
            self.dataType = 'line'
        else: raise TypeError( "Unsupported type: " + str( type( data ) ) )
        if dataType != None: self.dataType = dataType
        DataSetParameters.__init__( self, **kw )
        
    def convertUnits( self, xUnit = None, yUnit = None, zUnit = None ):
        '''Convert self's units to xUnit and yUnit and zUnit'''
        if xUnit != None and self.xUnit != None and self.xUnit != xUnit:
            try:  conversionFactor = valueOrPQ( 1.0, unitFrom=self.xUnit, unitTo=xUnit, asPQU=False )
            except TypeError, err: 
                raise TypeError( err.message + ': '+self.xUnit+' -> '+xUnit+' for set '+str( self.legend ) )
            except NameError, err:
                raise NameError( err.message + ': '+self.xUnit+' -> '+xUnit+' for set '+str( self.legend ) )
            if self.xerr != None: self.xerr = [ xx*conversionFactor for xx in self.xerr ]
            self.x = [ xx*conversionFactor for xx in self.x ]
            self.xUnit = xUnit
        if yUnit != None and self.yUnit != None and self.yUnit != yUnit:        
            try:  conversionFactor = valueOrPQ( 1.0, unitFrom=self.yUnit, unitTo=yUnit, asPQU=False )
            except TypeError, err: 
                raise TypeError( err.message + ': '+self.yUnit+' -> '+yUnit+' for set '+str( self.legend ) )
            except NameError, err:
                raise NameError( err.message + ': '+self.yUnit+' -> '+yUnit+' for set '+str( self.legend ) )
            if self.yerr != None: self.yerr = [ yy*conversionFactor for yy in self.yerr ]
            self.y = [ yy*conversionFactor for yy in self.y ]
            self.yUnit = yUnit
        if zUnit != None and self.zUnit != None and self.zUnit != zUnit:        
            try:  conversionFactor = valueOrPQ( 1.0, unitFrom=self.zUnit, unitTo=zUnit, asPQU=False )
            except TypeError, err: 
                raise TypeError( err.message + ': '+self.zUnit+' -> '+zUnit+' for set '+str( self.legend ) )
            except NameError, err:
                raise NameError( err.message + ': '+self.zUnit+' -> '+zUnit+' for set '+str( self.legend ) )
            if self.zerr != None: self.zerr = [ zz*conversionFactor for zz in self.zerr ]
            self.z = [ zz*conversionFactor for zz in self.z ]
            self.zUnit = zUnit
        
    def addToPlot( self, thePlot, plotType = 'contour', numContours = 10, logX = False, logY = False ):
        '''Default plotType set to 'contour' so that can use the __makePlot2d function to drive both the regular 2d plotting and contour plotting'''
        if self.dataType == 'line':
            if plotType == 'contour': 
                plt.clabel( thePlot.contour( self.x, self.y, self.z, numContours, **self.getFormatMap(  ) ), inline=1, fontsize=10 )
            else: raise TypeError( 'plotType ' + plotType + ' not supported' )
        else: raise TypeError( 'dataType must be "scatter" or "line", found ' + str( self.dataType ) )
        
    def getEndl3dmath( self ): 
        table = []
        for ix, x in enumerate( self.x ):
            table.append( [ x, [] ] )
            for iy, y in enumerate( self.y ):
                table[-1][1].append( [ y, self.z[ iy ][ ix ] ] )
        return endl3dmath( data = table )



def __makePlot2d( theSets, xAxisSettings=None, yAxisSettings=None, theTitle=None, 
        legendOn=False, legendXY=( 0.05, 0.95 ), figsize=( 20, 10 ), outFile=None, 
        thePlot=None, **kw ):
    '''
    
    Main driver routine for all 2d plots (regular, contour, slices, ...)
    
    '''
    if outFile == None: defaultBackend = 'TkAgg'
    else:  defaultBackend = 'Agg'
    backendMap = { 'png':'AGG', 'ps':'PS', 'eps':'PS', 'pdf':'PDF', 'svg':'SVG' }
    import matplotlib
    if matplotlib.get_backend() != defaultBackend: matplotlib.use( defaultBackend )
    import matplotlib.pyplot as plt

    # Create a plotting instance in a figure if one is needed.  
    # If the user passed in a plot, we'll return that.
    # Otherwise, all control of the plot will be handled in here.
    if thePlot == None: 
        returnPlotInstance = False
        # Set up the plot canvas
        fig = plt.figure( figsize=figsize )
        # Add the plot title
        if theTitle != None: fig.suptitle( theTitle.text, fontsize=theTitle.fontSize, fontweight=theTitle.fontWeight )
        # Set up the plot region
        thePlot = fig.add_subplot( 111 )
        fig.subplots_adjust( top=0.85 )    
    else: returnPlotInstance = True

    # Add the data to the plot
    for set in theSets: 
        try:
            set.convertUnits( xUnit = xAxisSettings.unit, yUnit = yAxisSettings.unit )
            set.addToPlot( thePlot, logY = yAxisSettings.isLog, **kw )
        except Exception, err: 
            print "WARNING: could not add set "+str(set.legend)+", found error:"
            print err.__class__, err
    # Put on the plot legend
    if legendOn:
        thePlot.legend( bbox_to_anchor = legendXY, loc = 2, borderaxespad = 0., markerscale = 1.0, ncol = max( 1, len( theSets )/20 ) ) #, scatterpoints = 1, frameon = False ) # these last 3 options don't work?

    # Set up the axes
    xAxisSettings.setupAxis( thePlot.get_xaxis() )
    yAxisSettings.setupAxis( thePlot.get_yaxis() )

    if not INCLASSAXISSETTING:
        # Set the scales using the SubPlot calls (direct calls to XAxis and YAxis don't work because the scale information is 
        # contained in the Axes class instance that contains both the XAxis and YAxis instances.  The Axis.set_limit functions
        # actually call the stuff in the Axes instance and there is no equivalent set_scale functions).
        if xAxisSettings.isLog: thePlot.set_xscale( 'log' )
        else:                   thePlot.set_xscale( 'linear' )
        if yAxisSettings.isLog: thePlot.set_yscale( 'log' )
        else:                   thePlot.set_yscale( 'linear' )
        if xAxisSettings.autoscale: thePlot.set_autoscalex_on(True)
        else:                       
            thePlot.set_xlim( ( xAxisSettings.axisMin, xAxisSettings.axisMax ) )
            thePlot.set_xbound( ( xAxisSettings.axisMin, xAxisSettings.axisMax ) )
        if yAxisSettings.autoscale: thePlot.set_autoscaley_on(True)
        else:                       
            thePlot.set_ylim( ( yAxisSettings.axisMin, yAxisSettings.axisMax ) )
            thePlot.set_ybound( ( yAxisSettings.axisMin, yAxisSettings.axisMax ) )
        thePlot.set_xlabel( xAxisSettings.label, size = xAxisSettings.fontSize, weight = xAxisSettings.fontWeight, family = xAxisSettings.font )
        thePlot.set_ylabel( yAxisSettings.label, size = yAxisSettings.fontSize, weight = yAxisSettings.fontWeight, family = yAxisSettings.font )
        for label in thePlot.get_xticklabels():
            label.set_size( xAxisSettings.fontSize )
            label.set_weight( xAxisSettings.fontWeight )
            label.set_family( xAxisSettings.font )
        for label in thePlot.get_yticklabels():
            label.set_size( yAxisSettings.fontSize )
            label.set_weight( yAxisSettings.fontWeight )
            label.set_family( yAxisSettings.font )
        # the following worked, but when you zoom in on the plot, the axes don't go too:
        #thePlot.set_xticklabels( thePlot.get_xticks(), size = xAxisSettings.fontSize, weight = xAxisSettings.fontWeight, family = xAxisSettings.font )
        #thePlot.set_yticklabels( thePlot.get_yticks(), size = yAxisSettings.fontSize, weight = yAxisSettings.fontWeight, family = yAxisSettings.font )
        thePlot.minorticks_on(  )


    if returnPlotInstance: return thePlot
    else:
        # Generate the output
        if outFile != None: plt.savefig(outFile)
        else:               plt.show()
        


def makePlot2d( sets, xAxisSettings = None, yAxisSettings = None, title = '', legendOn = False, outFile = None, legendXY = ( 0.05, 0.95 ), figsize=( 20, 10 ) ):
    '''
    Plain, vanilla, 2d plots
    
    * sets (required), list of DataSet objects to plot
    * xAxisSettings (optional) 
    * yAxisSettings (optional)
    * title (optional, defaults to '' )
    * legendOn (optional, defaults to True)
    
    '''
    # Process the function arguments
    if isinstance( title, TextSettings ): theTitle = title
    else:                                 theTitle = TextSettings( title, fontSize = 24 )
    if not isinstance( xAxisSettings, AxisSettings ): xAxisSettings = AxisSettings( )
    if not isinstance( yAxisSettings, AxisSettings ): yAxisSettings = AxisSettings( )
    if not isinstance( sets, list ): sets = [ sets ]
    theSets = [ ]
    for set in sets:
        if isinstance( set, DataSet2d ):    theSets.append( set )
        elif isinstance( set, endl2dmath ): theSets.append( DataSet2d( set ) )
        elif isinstance( set, XYs ):        theSets.append( DataSet2d( set ) )
        else:                               theSets.append( DataSet2d( set ) )
    __makePlot2d( theSets, xAxisSettings, yAxisSettings, theTitle = theTitle, legendOn = legendOn, legendXY = legendXY, figsize = figsize, outFile = outFile )


def makePlot2dContour( sets, xAxisSettings = None, yAxisSettings = None, numContours = 10, title = '', figsize=( 20, 10 ), legendOn = False, outFile = None ):
    '''
    Contour plots:
    
    * sets (required): list of DataSet objects to plot, one works better than many for a contour plot
    * xAxisSettings (optional) 
    * yAxisSettings (optional)
    * numContours (optional, defaults to 10): number of contours in the plot assuming linear spacing between intervals from the min to max data values
    * title (optional, defaults to '' )
    * legendOn (optional, defaults to True)
    * outFile (optional)
    
    '''
    # Process the function arguments
    if isinstance( title, TextSettings ): theTitle = title
    else:                                 theTitle = TextSettings( title, fontSize = 24 )
    if not isinstance( xAxisSettings, AxisSettings ): xAxisSettings = AxisSettings( )
    if not isinstance( yAxisSettings, AxisSettings ): yAxisSettings = AxisSettings( )
    if not isinstance( sets, list ): sets = [ sets ]
    theSets = [ ]
    for set in sets:
        if   isinstance( set, DataSet3d ):  theSets.append( set )
        elif isinstance( set, endl3dmath ): theSets.append( DataSet3d( set ) )
        else:                               theSets.append( DataSet3d( set ) )
    if len( sets ) != 1: print 'Warning:  Contour plots with more than one set may not be very readable'
    __makePlot2d( theSets, xAxisSettings, yAxisSettings, theTitle = theTitle, legendOn = legendOn, outFile = outFile, figsize = figsize, numContours = numContours )



def makePlot2dSlice( sets, xyAxisSettings = None, zAxisSettings = None, sliceXCoords = None, sliceYCoords = None, sliceUnits = '', title = '', legendOn = True, legendXY = ( 0.05, 0.95 ), figsize=( 20, 10 ), outFile = None ):
    '''
        Make a slice of a set of 3d data objects along a fixed x line or a fixed y line.  
        
        Any 2d objects get plotted as if they are xy data in the correct plane for plotting.
        
        In other words, at fixed sliceXCoord = 2, all 3d functions F( x, y ) get plotted 
        as y vs. F( 2.0, y ) and all 2d functions G( x ) are plotted as x vs. G( x ) even 
        if you meant something different.
        
        * sets (required): list of DataSet objects to plot
        * xyAxisSettings (optional) 
        * zAxisSettings (optional)
        * sliceXCoords: a list or value
        * sliceYCoords: a list or value
        * title (optional, defaults to '' )
        * legendOn (optional, defaults to True)
        * outFile (optional)
        
    '''
    # Process the function arguments
    if isinstance( title, TextSettings ): theTitle = title
    else:                                 theTitle = TextSettings( title, fontSize = 24 )
    if not isinstance( xyAxisSettings, AxisSettings ): xyAxisSettings = AxisSettings( )
    if not isinstance( zAxisSettings, AxisSettings ): zAxisSettings = AxisSettings( )
    if sliceXCoords == None and sliceYCoords == None: 
        raise ValueError( "sliceXCoords or sliceYCoords must have a list or value, otherwise where do I slice?" )
    if sliceXCoords != None and sliceYCoords != None: 
        raise ValueError( "sliceXCoords or sliceYCoords both have a list or value, which one do I slice?" )
    if not isinstance( sets, list ): sets = [ sets ]
    theSets = [ ]
    for set in sets:
        if isinstance( set, DataSet2d ):    theSets.append( set )
        elif isinstance( set, endl2dmath ): theSets.append( DataSet2d( set ) )
        elif isinstance( set, DataSet3d ) or isinstance( set, endl3dmath ): 
            if isinstance( set, endl3dmath ):
                asDataSet3d = DataSet3d( set )
                asEndl3dmath = set
            else:
                asDataSet3d = set
                asEndl3dmath = set.getEndl3dmath()
            if sliceXCoords != None:
                for x in sliceXCoords:
                    theSets.append( DataSet2d( endl2dmath( data = [ [ y, asEndl3dmath.getValue( x, y ) ] for y in asDataSet3d.y ] ), legend = asDataSet3d.legend + ' @ ' + str( x ) + ' ' + sliceUnits ) )
            if sliceYCoords != None:
                for y in sliceYCoords:
                    theSets.append( DataSet2d( endl2dmath( data = [ [ x, asEndl3dmath.getValue( x, y ) ] for x in asDataSet3d.x ] ), legend = asDataSet3d.legend + ' @ ' + str( y ) + ' ' + sliceUnits ) )
        else: raise TypeError( 'Can only put 2d or 3d objects in slice plots' )
    __makePlot2d( theSets, xyAxisSettings, zAxisSettings, theTitle = theTitle, legendOn = legendOn, outFile = outFile, legendXY = legendXY, figsize = figsize )





#--------------------------------------------------------
#
#  testing
#
#--------------------------------------------------------
def plotTests( tests = 11*[ False ] ):
    testData = '''
        #  Authors:   T.S.Suzuki, Y.Nagai, T.Shima, T.Kikuchi, H.Sato, T.Kii, M.Igashira
        #  Title:     First Measurement Of A P(N,Gamma)D Reaction Cross Section Between 10 And 80 Kev
        #  Year:      1995
        #  Institute: Tokyo Inst.of Technology, Tokyo
        #  Reference: Astrophysical Journal 439, (L), 59 (1995)
        #  Subent:    22310002
        #  Reaction:  Cross section for 1H(n,gamma)2H 
        # Note: the d(Energy) errorbars are fake and are used for plot testing
        #   Energy        Data          d(Energy)     d(Data)       
        #   MeV           barns         MeV           barns         
            0.02          0.000353      0.001         2.6e-05       
            0.02          0.000329      0.001         2.6e-05       
            0.02          0.000287      0.0           2.2e-05       
            0.02          0.000304      0.0           1.8e-05       
            0.04          0.00023       0.001         1.5e-05       
            0.04          0.000198      0.001         1.2e-05      
            0.04          0.000177      0.0015        1e-05         
            0.04          0.000207      0.0           1.3e-05       
            0.064         0.000156      0.0           1.1e-05       
            0.064         0.00015       0.0           7e-06         
            0.064         0.000158      0.0           1.1e-05       
            0.064         0.00014       0.0           9e-06
    '''
    from fudge.legacy.endl.endlProject import endlProject
    from fudge import __path__
    e=endlProject( __path__[0] + "/legacy/endl/test/testdb" )
    za=e.readZA(1001)
    za.read( )
    xAxis = AxisSettings( isLog = True, label = '$E_n$ (MeV)', axisMin = 0.5e-2, axisMax = 1.5e-1, gridOn = True, autoscale = False )
    yAxis = AxisSettings( isLog = True, label = '$\sigma(E_n)$ (barns)', gridOn = True )
    d = []
    u = []
    for line in testData.split( '\n' ):
        if line.strip().startswith( '#' ): continue
        sline = map( float, line.split() )
        if len( sline ) != 0:
            d.append( sline[0:2] )
            u.append( sline[2:4] )

    # Simple test, here we make a plot of one set
    if tests[0]:
        makePlot2d( [ za.findData( I = 0, C = 46 ) ], xAxisSettings = xAxis, yAxisSettings = yAxis, title = '$^1$H$(n,\gamma)$ Cross Section', outFile = None )        

    # Plot all the cross section data in the fudge2 test library, but unthemed!
    if tests[1]:
        makePlot2d( za.findDatas( I = 0 ), outFile = None )
            
    # Plot all the cross section data in the fudge2 test library
    if tests[2]:
        makePlot2d( za.findDatas( I = 0 ), xAxisSettings = xAxis, yAxisSettings = yAxis, title = '$^1$H$(n,*)$ Cross Sections', outFile = None )   

    # Fancy test, here we make a plot of one set (the testData above)
    if tests[3]:
        endfData = za.findData( I=0, C=46 ) #.slicex(xMin=1e-2,xMax=1e-1)
        endfUnc = 0.1 * endfData # 10% error bars
        makePlot2d( [ \
            DataSet2d( data = endfData, uncertainty = endfUnc, legend = 'ENDF/B-VII.0', color='g', lineStyle = '-' ),
            DataSet2d( data = d, uncertainty = u, legend = 'T.S.Suzuki, et al. (1995) EXFOR entry # 22310002', color='g', symbol = 'o' ),
        ], xAxisSettings = xAxis, yAxisSettings = yAxis, title = '$^1$H$(n,\gamma)$ Cross Section', legendOn = True, outFile = None )

    # Contour test
    if tests[4]:
        # Contour plots are meaningful if there is only one set plotted, how do we enforce this?
        endfData = za.findData( yo = 1, I = 1, C = 10 )
        EAxis  = AxisSettings( isLog = False, label = '$E_n$ (MeV)', axisMin = 1.0, axisMax = 20.0, gridOn = True, autoscale = False )
        muAxis = AxisSettings( isLog = False, label = '$\mu$ = cos( $\\theta$ )', axisMin = -1.0, axisMax = 1.0, autoscale = False, gridOn = True )
        makePlot2dContour( DataSet3d( data = endfData, legend = 'ENDF/B-VII.0' ), xAxisSettings = EAxis, yAxisSettings = muAxis, title = '$^1$H$(n,el)$ Angular Distribution', outFile = None )

    # Slice tests
    if tests[5]:
        endfDataI0 = za.findData( yo = 0, I = 0, C = 10 ).slicex(1.0,20.0) # simplify the plot
        endfDataI1 = za.findData( yo = 1, I = 1, C = 10 )
        Es = [ p[0] for p in endfDataI0.data ]  # in [MeV]
        mus = [ ] # in [mu]
        for t in endfDataI1.data: mus += [ p[0] for p in t[1] ]
        mus = sorted( uniquify( mus ) )
        table = []
        for E in Es:
            muDist = []
            for mu in mus: muDist.append( [ mu, endfDataI0.getValue( E ) * endfDataI1.getValue( E, mu ) ] )
            table.append( [ E, muDist ] )
        endfDataCombined = endl3dmath( data = table )
        EAxis  = AxisSettings( isLog = True, label = '$E_n$ (MeV)', axisMin = 1.0, axisMax = 20.0, gridOn = True, autoscale = False )
        muAxis = AxisSettings( isLog = False, label = '$\mu = \cos{( \\theta )}$', axisMin = -1.0, axisMax = 1.0, autoscale = False, gridOn = True )
        zAxis = AxisSettings( isLog = True, label = '$d\sigma(E)/d\mu$ (barns)', axisMin = 0.0, axisMax = 5.0, autoscale = False, gridOn = True )
        
        # unthemed slice tests
        if tests[6]:makePlot2dSlice( endfDataCombined, xyAxisSettings = EAxis, zAxisSettings = zAxis, sliceXCoords = None, sliceYCoords = [ -1.0, -0.75, -0.5, -.25,  0.0, .25, .5, .75, 1.0 ], title = '', outFile = None )
        if tests[7]:makePlot2dSlice( endfDataCombined, xyAxisSettings = muAxis, zAxisSettings = zAxis, sliceXCoords = [ 1.0, 5.0, 10.0, 15.0, 20.0 ], sliceYCoords = None, title = '', outFile = None )

        # themed slice test and contour test
        if tests[8]:makePlot2dContour( DataSet3d( data = endfDataCombined, legend = 'ENDF/B-VII.0' ), xAxisSettings = EAxis, yAxisSettings = muAxis, numContours = 10, title = '$d\sigma(E)/d\mu$ for $^1$H$(n,el)$', outFile = None )
        if tests[9]: makePlot2dSlice( DataSet3d( data = endfDataCombined, legend = 'ENDF/B-VII.0' ), xyAxisSettings = EAxis, zAxisSettings = zAxis, sliceXCoords = None, sliceYCoords = [ -1.0, -0.75, -0.5, -.25,  0.0, .25, .5, .75, 1.0 ], title = '', outFile = None )

        # themed slice test, with test 2d experimental data
        if tests[10]: 
            makePlot2dSlice( [ \
                DataSet3d( data = endfDataCombined, legend = 'ENDF/B-VII.0' ), \
                DataSet2d( data = d, uncertainty = u, legend = 'T.S.Suzuki, et al. (1995) EXFOR entry # 22310002', color='g', symbol = 'o' ),\
            ], xyAxisSettings = muAxis, zAxisSettings = zAxis, sliceXCoords = [ 1.0, 5.0, 10.0, 15.0, 20.0 ], sliceYCoords = None, sliceUnits = 'MeV', title = 'Slice test', outFile = None )


#--------------------------------------------------------
#
# main!
#
#--------------------------------------------------------
if __name__ == "__main__":
    plotTests( tests = [ False, False, False, False, False, True, False, False, False, False, True ]  )
