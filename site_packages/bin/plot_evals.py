#! /usr/bin/env python

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
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# <<END-copyright>>

import collections, os
try:                import argparse
except ImportError: from fudge.core.utilities import argparse

from BNL import data_io as io
from BNL import utils, endf_utils, gnd_utils, plot_utils

import fudge
from fudge.vis.matplotlib import plot2d, DataSet2d, DataSet3d, defaultSymbols, nSymbols, defaultColors, nColors
from fudge.core.utilities.brb import uniquify

#---------------------------------------------------
# Set up the command line parser
#---------------------------------------------------
def parseArgs():
    parser = argparse.ArgumentParser(description='Plot nuclear data from an ENDF or GND file')
    
    # Required things so we know what to plot
    parser.add_argument('mt',       metavar='mt', type=int, help='MT of the cross section to plot.  If set to 0, will try to make all plots for all open channels (requires -o option too)' )
    parser.add_argument('endf',     metavar='endf', type=str, nargs='+', help='ENDF file(s) whose cross section you want to plot.  Use "None" for no input file.  AMPX BOF files and tables of X,Y(,dY) pairs may work too.' )
    
    # Plot output file controls
    parser.add_argument('-o',       dest='outFile', default=None, type=str, help='Output file for plot (disables interactive plotting)' )

    # Override the target/projectile/product of interest
    parser.add_argument('--target',     default=None, type=str, help="The target nucleus, given in GND notation, e.g. 'Pu239' (Default is None which means to take it from the first ENDF file)" )
    parser.add_argument('--projectile', default=None, type=str, help="The projectile particles, given in GND notation, e.g. 'n' (Default is None which means to take it from the first ENDF file)" )
    parser.add_argument('--product',    default=None, type=str, help="The product particle of interest, given in GND notation, e.g. 'n' (Default is None which means to take the first emitted particle for this observable)" )

    # Resolved resonance region reconstruction controls
    parser.add_argument('--enableRRAngDist', default=False, action='store_true', help='Reconstruct the angular distribution from the resolved resonance parameters (default: False)' )
    parser.add_argument('--noReconstruct',   dest='doResonanceReconstruction', default=True, action='store_false', help="Don't reconstruct resonances (default: True)'" )
    parser.add_argument('--showURRCloud',    default=False, action='store_true', help='Overlay a contour plot of the PDF for the cross section in the URR (default: False)' )

    # Uncertainty related options
    parser.add_argument('--uncRatio',       default=False, action='store_true', help='Plot a ratio of the uncertainty over the data (default: False)' )
    parser.add_argument('--noUnc',          default=False, action='store_true', help='Do not plot uncertainties' )

    # Experimental data sources & controls for cross comparison
    parser.add_argument('--c4File',         default=None,  type=str, help="Optional C4 file to pull data from instead of (or in addition to) EXFOR" )
    parser.add_argument('--noX4',           default=False, action='store_true', help='Do not plot EXFOR data' )
    parser.add_argument('--showX4Evals',    default=True, action='store_false', help='Plot evaluations found in EXFOR library (with "V" SUBENT).  The default behavior is not to plot them.' )

    # Overall plot style controls.  Overrides defaults in DEFAULT_STYLE_DICT
    parser.add_argument('--style',          default=None, type=str, help='JSON file with plot style overrides, see "plot_defaults.json" in source distribution for examples')

    # Plot legend controls
    parser.add_argument('--noX4Legend',     default=False, action='store_true', help='Do not put legend for EXFOR data on plot' )
    parser.add_argument('--showX4Legend',   default=True,  action='store_true', help='Put legend for EXFOR data on plot, even if it would be automatically supressed' )
    parser.add_argument('--simpleLegend',   default=False, action='store_true', help='Do not put entries for the components of requested MT in the plot legend' )

    # Detailed control of the observable to plot
    parser.add_argument('--observable',  choices=plot_utils.DEFAULT_STYLE_DICT["plotStyles"], type=str, default="crossSection", help='The observable to plot (default is "crossSection")' )
    parser.add_argument('--showParts',   default=False, action='store_true', help='Show all of the channels that comprise the one requested (if applicable)' )
    parser.add_argument('--L',           type=int, default=0, help='Legendre moment to show (if applicable)' )
    parser.add_argument('--referenceFrame', choices=['com','lab'], type=str, default=None, help='The frame of the observable to plot.' )
    
    return parser.parse_args()


#---------------------------------------------------
# ENDF-specific plotting helpers
#---------------------------------------------------  

# Simplified MT - EXFOR reaction mapping
def getEXFORRxn( MT, projectile='N' ):
    if MT == 2: return 'EL'
    if MT in range( 50, 92 ) or MT==4:    
        if projectile == 'N': return 'INEL'
        else: return 'N'
    if MT in range( 600, 650 ) or MT==103:  
        if projectile == 'P': return 'INEL'
        else: return 'P'
    if MT in range( 650, 700 ) or MT==104:  
        if projectile == 'D': return 'INEL'
        else: return 'D'
    if MT in range( 700, 750 ) or MT==105:  
        if projectile == 'T': return 'INEL'
        else: return 'T'
    if MT in range( 750, 800 ) or MT==106:  
        if projectile == 'HE3': return 'INEL'
        else: return 'HE3'
    if MT in range( 800, 850 ) or MT==107:  
        if projectile == 'A': return 'INEL'
        else: return 'A'
    if MT in range( 875, 891 ):  return '2N'
    if MT in [ 19, 20, 21, 38 ]: return 'F'
    return { \
        1:'TOT', 2:'EL', 3:'NON', 4:'N', 5:'X', 11:'2N+D', 16:'2N', 17:'3N', 18:'F', \
        22:'N+A', 23:'N+3A', 24:'2N+A', 25:'3N+A', 27:'ABS', 28:'N+P', 29:'N+2A', \
        30:'2N+2A', 32:'N+D', 33:'N+T', 34:'N+HE3', 35:'N+D+2A', 36:'N+T+2A', 37:'4N', \
        41:'2N+P', 42:'3N+P', 44:'N+2P', 45:'N+P+A', 102:'G', 103:'P', 104:'D', 105:'T', \
        106:'HE3', 107:'A', 108:'2A', 109:'3A', 111:'2P', 112:'P+A', 113:'T+2A', \
        114:'D+2A', 115:'P+D', 116:'P+T', 117:'D+A' }.get( MT, "MT="+str(MT))


def getSuggestTitle( target, projectile, reaction, mt ):
     if reaction is None: return target.capitalize() + '(' + projectile+',X)'
     if reaction.lower() == 'inel': 
         if   mt in range(51,91):    reactionString = "n["+str( mt%50 )+"]"
         elif mt == 91:              reactionString = "n[c]"
         elif mt in range(600,649):  reactionString = "p["+str( mt%50 )+"]"
         elif mt == 649:             reactionString = "p[c]"
         elif mt in range(650,691):  reactionString = "d["+str( mt%50 )+"]"
         elif mt == 691:             reactionString = "d[c]"
         elif mt in range(700,749):  reactionString = "t["+str( mt%50 )+"]"
         elif mt == 749:             reactionString = "t[c]"
         elif mt in range(750,791):  reactionString = "3He["+str( mt%50 )+"]"
         elif mt == 791:             reactionString = "3He[c]"
         elif mt in range(800,849):  reactionString = "a["+str( mt%50 )+"]"
         elif mt == 849:             reactionString = "a[c]"
         elif mt in range(850,891):  reactionString = "2n["+str( mt%50 )+"]"
         elif mt == 891:             reactionString = "2n[c]"
         else:                       reactionString = 'inel'
     else: reactionString = reaction.lower()
     return target.capitalize() + '(' + projectile+','+reactionString +')'

 
def generatePlot( observable, dataSets, xyData=[], xydyData=[], xdxydyData=[], plotStyle={}, suggestTitle='', suggestXLog=False, suggestYLog=False, suggestFrame=None ):

    # Set up the XY data
    xyDataSets=[]
    if len( xyData ) != 0:  xyDataSets=getXYDataSets( xyData, plotStyle )

    # Set up the XYdY data
    xydyDataSets=[]
    if len( xydyData ) != 0:  xydyDataSets=getXYdYDataSets( xydyData, plotStyle )

    # Set up the XdXYdY data
    xdxydyDataSets=[]
    if len( xdxydyData ) != 0:  xdxydyDataSets=getXdXYdYDataSets( xdxydyData, plotStyle )

    # Set up this plot styles
    thisPlotStyle = plot_utils.getThisPlotStyle( plotStyle, observable )
    
    # Whether to do lin-lin or log-log based on suggestLog flag
    if thisPlotStyle["xAxis"]["log"] == None: 
        if suggestXLog: thisPlotStyle["xAxis"]["log"] = True
        else: thisPlotStyle["xAxis"]["log"] = False
    if thisPlotStyle["yAxis"]["log"] == None:
        if suggestYLog: thisPlotStyle["yAxis"]["log"] = True
        else: thisPlotStyle["yAxis"]["log"] = False
        
    # Check the reference frame
    if thisPlotStyle['referenceFrame'] != None:
        if suggestFrame != thisPlotStyle['referenceFrame']: 
            raise ValueError( "Suggested frame from data of %s does not match required plot style %s" % (suggestFrame, thisPlotStyle['referenceFrame']) )

    # load the axis style information into the plotting widget class instances
    xAxisSettings = plot2d.AxisSettings( axisMin = thisPlotStyle["xAxis"]["min"], \
                                         axisMax = thisPlotStyle["xAxis"]["max"],\
                                         label   = thisPlotStyle["xAxis"]["label"], \
                                         isLog   = thisPlotStyle["xAxis"]["log"], \
                                         unit    = thisPlotStyle["xAxis"]["unit"] )
    yAxisSettings = plot2d.AxisSettings( axisMin = thisPlotStyle["yAxis"]["min"], \
                                         axisMax = thisPlotStyle["yAxis"]["max"],\
                                         label = thisPlotStyle["yAxis"]["label"], \
                                         isLog = thisPlotStyle["yAxis"]["log"], \
                                         unit  = thisPlotStyle["yAxis"]["unit"] )
    # Set the plot title, if needed
    if thisPlotStyle['title']==None: thisPlotStyle['title'] = suggestTitle

    # Finally make the plot
    plot2d.makePlot2d( \
                      dataSets, \
                      xAxisSettings = xAxisSettings, \
                      yAxisSettings = yAxisSettings,\
                      title = thisPlotStyle['title'],
                      legendOn = True, \
                      outFile = outFile,\
                      legendXY=(thisPlotStyle['legendX'],thisPlotStyle['legendY']) )


#---------------------------------------------------
# stuff for interacting with various data sources
#---------------------------------------------------  
    
def getEXFORSets( sym, A, reaction = None, quantity = "SIG", nox4evals = True, nox4legend = False, forceLegend = False, plotSyle={} ): 
    exforData = []
    try:
        from x4i import exfor_manager, exfor_entry
    except ImportError:
        print 'WARNING: x4i not successfully imported (check your PYTHONPATH?), so EXFOR data not plotted'
        return exforData
    def barnsConverter( x ):
        if x == 'barns': return 'b'
        else: return x
    db = exfor_manager.X4DBManagerPlainFS( )  
    fileList = []
    i = 0
    subents = db.retrieve( target = sym+'-'+A, reaction = reaction, quantity = quantity )
    print utils.bigBanner( "Preparing EXFOR data for " + sym+'-'+A +'(' +reaction.upper()+')'+', '+ quantity )
    print 'Retrieving entries:'+str( subents.keys() )
    suppressEXFORLegend = ( nox4legend or len( subents.keys() ) > 20 ) and not forceLegend
    print 'Search parameters:', sym+'-'+A, reaction, quantity
    print 'Found following (sub)entries:'
    for e in subents:
        if nox4evals and e.startswith( 'V' ): continue
        print '    Entry:',e
        try:        
            # New version of x4i returns an X4Entry directly
            if isinstance( subents[ e ], exfor_entry.X4Entry ): ds = subents[ e ].getSimplifiedDataSets( makeAllColumns = True )
            # Old versions of x4i return a string that we have to convert to an X4Entry
            else: ds = exfor_entry.X4Entry( subents[ e ] ).getSimplifiedDataSets( makeAllColumns = True )
        except KeyError: continue
        for d in ds:
            # TO DO: skip if data ratio data or SPA or MXW ave'd
            legend = ds[d].legend()
            print '       ',d, legend
            dat = []
            unc = []
            for line in ds[d].data:
                if len( line ) != 4: continue
                if line[0] == None or line[1] == None: continue
                dx = 0.0
                dy = 0.0
                if line[2] != None: dx = line[2]
                if line[3] != None: dy = line[3]
                dat.append( [ line[0], line[1] ] )
                unc.append( [ dx, dy ] )
                if None in unc[-1]: unc[-1][ unc[-1].index( None ) ] = 0.0
            if len( dat ) > 0:
                if e.startswith( 'V' ) or not suppressEXFORLegend: theLegend = legend+' ('+str(d[0])+'.'+str(d[1][-3:])+')'
                else: theLegend = '_noLegend_'
                exforData.append( DataSet2d( data = dat, uncertainty = unc, xUnit=ds[d].units[0], yUnit=barnsConverter( ds[d].units[1] ), legend = theLegend, lineStyle = ' ', symbol = plot_utils.getPlotSymbol( i ), color=plot_utils.getPlotColor(theLegend,False) ) )
                i += 1
    print
    return exforData
    

def getC4DataSets( c4File, mt, mf, target, projectile, plotSyle={} ):
    try:
        from empire import c4
    except ImportError:
        try: 
            import c4
        except ImportError:
            print 'WARNING: c4 not successfully imported (check your PYTHONPATH?), so X4TOC4 data not plotted'
            return []
    print utils.smallBanner( 'Retrieving C4 data from: '+str( c4File ) )
    flist = filter( lambda x: x.MT == mt and x.MF == mf and x.target == gnd_utils.getZAFromGNDName( target ) and x.projectile == gnd_utils.getZAFromGNDName( projectile ), c4.readC4File( open( c4File ).readlines(), asPointList=True ) )
    i = -1
    dat = []
    unc = []
    c4Data = []
    theLegend = None
    lastSet = ( None, None )
    thisSet = ( None, None )
    for p in flist:
        thisSet = ( p.reference, p.exforEntry+str(p.exforSubEntry).zfill(3) )
        if lastSet != thisSet: 
            if lastSet != ( None, None ): 
                print '    Set:',lastSet
                theLegend = lastSet[0]+' ('+lastSet[1]+')'
                if mf == 3:   c4Data.append( DataSet2d( data = dat, uncertainty = unc, xUnit='eV', yUnit='b',           legend = theLegend, lineStyle = ' ', symbol = plot_utils.getPlotSymbol( i ), color=plot_utils.getPlotColor(theLegend,False) ) )
                elif mf == 4: c4Data.append( DataSet3d( data = dat, uncertainty = unc, xUnit='eV', yUnit='', zUnit='b', legend = theLegend, lineStyle = ' ', symbol = plot_utils.getPlotSymbol( i ), color=plot_utils.getPlotColor(theLegend,False) ) )
                dat = []
                unc = []
            i+=1
            lastSet = thisSet
        dat.append( [ p.energy, p.data ] )
        if p.dEnergy!=None: unc.append( [ p.dEnergy ] )
        else:               unc.append( [ 0.0 ] )
        if p.dData!=None:   unc[-1].append( p.dData )
        else:               unc[-1].append( 0.0 )
        if mf == 4:
            dat[-1].append( p.cosMuOrLegendreOrder )
            if p.dCosMuOrLegendreOrder!=None:   unc[-1].append( p.dCosMuOrLegendreOrder )
            else:                               unc[-1].append( 0.0 )
            cmFlag = p.cmFlag
    lastSet = thisSet
    if lastSet == (None, None): return []
    print '    Set:',lastSet
    theLegend = lastSet[0]+' ('+lastSet[1]+')'
    if mf == 3:   c4Data.append( DataSet2d( data = dat, uncertainty = unc, xUnit='eV', yUnit='b',           legend = theLegend, lineStyle = ' ', symbol = plot_utils.getPlotSymbol( i ), color=plot_utils.getPlotColor(theLegend,False) ) )
    elif mf == 4: c4Data.append( DataSet3d( data = dat, uncertainty = unc, xUnit='eV', yUnit='', zUnit='b', legend = theLegend, lineStyle = ' ', symbol = plot_utils.getPlotSymbol( i ), color=plot_utils.getPlotColor(theLegend,False) ) )
    print
    return c4Data


def getXYDataSets( xyData, plotStyle ):
    xyDataSets=[]
    for k in xyData:
        thisSetStyle=plot_utils.getThisSetStyle(plotStyle,'evaluation',k)
        xyDataSets.append( DataSet2d( data = xyData[k],\
                                      legend    = thisSetStyle[ 'legend' ], \
                                      lineWidth = thisSetStyle[ 'lineWidth' ], \
                                      lineStyle = thisSetStyle[ 'lineStyle' ], \
                                      color     = thisSetStyle[ 'lineColor' ],\
                                      symbol    = thisSetStyle[ 'symbol' ],\
                                      dataType  = thisSetStyle[ 'dataType' ],\
                                      xUnit     = thisSetStyle[ 'xUnit' ],\
                                      yUnit     = thisSetStyle[ 'yUnit' ] ) )
    return xyDataSets


def getXYdYDataSets( xydyData, plotStyle ):
    # Set up theXYdY data
    xydyDataSets=[]
    for k in xydyData:
        thisSetStyle=plot_utils.getThisSetStyle(plotStyle,'evaluation',k)
        d = [ [x[0],x[1]] for x in xydyData[k] ]
        u = [ [0.0 ,x[2]] for x in xydyData[k] ]
        xydyDataSets.append( DataSet2d( data = d, \
                                        uncertainty = u,\
                                        legend        = thisSetStyle[ 'legend' ], \
                                        lineWidth     = thisSetStyle[ 'lineWidth' ], \
                                        lineStyle     = thisSetStyle[ 'lineStyle' ], \
                                        color         = thisSetStyle[ 'lineColor' ], \
                                        symbol        = thisSetStyle[ 'symbol' ],\
                                        dataType      = thisSetStyle[ 'dataType' ],\
                                        xUnit         = thisSetStyle[ 'xUnit' ],\
                                        yUnit         = thisSetStyle[ 'yUnit' ],\
                                        errorbarColor = thisSetStyle[ 'errorColor' ] ) )
    return xydyDataSets


def getXdXYdYDataSets( xdxydyData, plotStyle ):
    # Set up the XdXYdY data
    xdxydyDataSets=[]
    for k in xdxydyData:
        thisSetStyle=plot_utils.getThisSetStyle(plotStyle,'evaluation',k)
        d = [ [x[0],x[2]] for x in xdxydyData[k] ]
        u = [ [x[1],x[3]] for x in xdxydyData[k] ]
        xdxydyDataSets.append( DataSet2d( data = d, \
                                        uncertainty = u,\
                                        legend        = thisSetStyle[ 'legend' ], \
                                        lineWidth     = thisSetStyle[ 'lineWidth' ], \
                                        lineStyle     = thisSetStyle[ 'lineStyle' ], \
                                        color         = thisSetStyle[ 'lineColor' ], \
                                        symbol        = thisSetStyle[ 'symbol' ],\
                                        dataType      = thisSetStyle[ 'dataType' ],\
                                        xUnit         = thisSetStyle[ 'xUnit' ],\
                                        yUnit         = thisSetStyle[ 'yUnit' ],\
                                        errorbarColor = thisSetStyle[ 'errorColor' ] ) )
    return xdxydyDataSets


#---------------------------------------------------
# cross section plots
#---------------------------------------------------
def makeCrossSectionPlot( gndMap={}, xyData={}, xydyData={}, xdxydyData={}, mt=0, projectile='n', target='1H', 
                          nounc=False, nox4=False, showparts=False, nox4evals=True, nox4legend=False, doResonanceReconstruction=True, heatToTemp=None,  
                          outFile=None, plotStyle={} ):

    # Declare the main reaction type
    try: reaction = getEXFORRxn( mt )
    except KeyError: 
        if mt in range( 850, 871, 1 ): print( "Got lumped covariance for MT = "+ str( mt ) )
        else: raise KeyError( "Unknown MT: "+str( mt ) )
    
    # Get the ENDF data for plotting
    endfData=[]
    rawEndfData=[]
    rawEndfUnc=[]
    for endf in gndMap:
        print utils.smallBanner( "Preparing data for "+endf )
        reactionSuite, covarianceSuite = gndMap[endf]

        #try:
        if True:
        
            # Set up MT list.  
            # This is the list of all MT that we will search for data to sum into the plot.  
            # Usually, it is one element long, but there are some special cases: summed channels that were not specified in the ENDF file and lumped channels....
            mtList = [ mt ]
            if not endf_utils.hasMT( reactionSuite, mt ):
                if   mt == 3: pass # (n,non-el), use for catch-all "C=55" gammas
                elif mt == 5: pass  # (n,X), no sumrule specified.  Treat as plain old reaction. 
                elif mt == 1:          
                    # (n,tot)
                    # Need better solution here!!!!
                    print( "WARNING: fudge-2.0 considers (n,tot) to be redundant so will attempt to compute it from summing parts" )
                    if 1 in mtList: del( mtList[ mtList.index(1) ] )
                    mtList += endf_utils.getEvaluationMTs( reactionSuite, mtFilter = [ 2 ] + range( 5, 120, 1 ) )
                    mtList += endf_utils.getEvaluationMTs( reactionSuite, mtFilter = range( 51, 92, 1 ) )
                    mtList += endf_utils.getEvaluationMTs( reactionSuite, mtFilter = range( 600, 650, 1 ) )
                    mtList += endf_utils.getEvaluationMTs( reactionSuite, mtFilter = range( 650, 700, 1 ) )
                    mtList += endf_utils.getEvaluationMTs( reactionSuite, mtFilter = range( 700, 750, 1 ) )
                    mtList += endf_utils.getEvaluationMTs( reactionSuite, mtFilter = range( 750, 800, 1 ) )
                    mtList += endf_utils.getEvaluationMTs( reactionSuite, mtFilter = range( 800, 850, 1 ) )
                    mtList += endf_utils.getEvaluationMTs( reactionSuite, mtFilter = range( 875, 892, 1 ) )
                    if endf_utils.hasMT( reactionSuite, 18 ): mtList += [ 18 ]
                    else: mtList += endf_utils.getEvaluationMTs( reactionSuite, mtFilter = [ 19, 20, 21, 38  ] )
                elif mt == 4:          
                    # (n,n'), sum: 51-91
                    print "Summing (n,n') reactions MT=51-91"
                    if 4 in mtList: del( mtList[ mtList.index(4) ] )
                    mtList += endf_utils.getEvaluationMTs( reactionSuite, mtFilter = range( 51, 92, 1 ) )
                elif mt == 103:        
                    # (n,p), sum: 600-649
                    print "Summing (n,p) reactions MT=600-649"
                    if 103 in mtList: del( mtList[ mtList.index(103) ] )
                    mtList += endf_utils.getEvaluationMTs( reactionSuite, mtFilter = range( 600, 650, 1 ) )
                elif mt == 104:        
                    # (n,d), sum: 650-699
                    mtList += endf_utils.getEvaluationMTs( reactionSuite, mtFilter = range( 650, 700, 1 ) )
                elif mt == 105:        
                    # (n,t), sum: 700-749
                    mtList += endf_utils.getEvaluationMTs( reactionSuite, mtFilter = range( 700, 750, 1 ) )
                elif mt == 106:        
                    # (n,3He), sum: 750-799
                    mtList += endf_utils.getEvaluationMTs( reactionSuite, mtFilter = range( 750, 800, 1 ) )
                elif mt == 107:        
                    # (n,a), sum: 800-849
                    mtList += endf_utils.getEvaluationMTs( reactionSuite, mtFilter = range( 800, 850, 1 ) )
                elif mt == 16:         
                    # (n,2n), sum 875-891
                    mtList += endf_utils.getEvaluationMTs( reactionSuite, mtFilter = range( 875, 892, 1 ) )
                elif mt == 18:
                    # (n,f), sum: 19-21,38
                    mtList += endf_utils.getEvaluationMTs( reactionSuite, mtFilter = [ 19, 20, 21, 38 ] ) 
            elif mt in range( 851, 871, 1 ): 
                # lumped covariance MT's are 851-870
                # must extract the MT's from the covariance file
                for lump in covarianceSuite.lumpedChannels:
                    if lump.ENDF_MFMT == '33,'+str(mt): break
                mtList += endf_utils.getEvaluationMTs( reactionSuite, mtFilter = [ int( channel.attributes['ENDF_MFMT'].split(',')[1] ) for channel in lump.channels] )
                reaction = 'N,' + getEXFORRxn( mtList[-1] )
            else: pass # plain old reaction, do nothing special
            mtList = uniquify( mtList )

            # Do resonance reconstruction if needed
            if doResonanceReconstruction:
                for MT in mtList: 
                    rxnList = endf_utils.getReactions( reactionSuite, MT )
                    if rxnList != [] and rxnList[0].crossSection.nativeData == 'resonancesWithBackground':
                        print "Attempting to reconstruct the resonances for",MT
                        try: reactionSuite.reconstructResonances(styleName='reconstructed')
                        except Exception, err: print 'WARNING: when reconstructing resonances for '+endf+' got '+str(err)
                        break
                        
            # Heat if needed
            if heatToTemp!=None: raise NotImplementedError( "Don't know how to heat cross sections yet" )
            
            print endf, mtList
            
            # collect complete list of matching reactions
            reactionList = []
            for MT in mtList: 
                try:
                    reactionList += endf_utils.getReactions( reactionSuite, MT )
                except: pass
            print "Exit Channels:", ', '.join( map( str, reactionList ) )
            
            # get the reaction corresponding to the requested mainMT
            if len( reactionList ) == 1: mainReaction = 0
            else:                        mainReaction = 'sum'
            
            # now plot the main reaction, this one gets the covariance
            if mainReaction != None:
                if True: #try:
                    
                    if reactionList == []: continue
                    
                    # Extract the cross section data from the endf file
                    csData = endf_utils.getPointwiseCrossSection( reactionList[0] )
                    if ( not nounc ) and mt not in [ 1, 2, 18, 102 ] :
                        try:
                            csData = csData.thicken( sectionSubdivideMax = 20 ) # only thicken cross sections that don't have enough points in them to resolve uncertainty steps
                        except Exception,err:
                            print "WARNING: "+str(err)
                    if mt != 18: # element 0 already was sum
                        for endfRxn in reactionList[1:]:
                            csData += endf_utils.getPointwiseCrossSection( endfRxn ) # this will need try/except block sooner or later
                    rawEndfData.append( csData )
                    
                    # Extract the uncertainty on the cross section
                    if nounc or covarianceSuite == None: rawEndfUnc.append( None )
                    else:
                        cov = [ c for c in covarianceSuite.sections if hasattr(c,'rowData') and c.rowData.attributes['ENDF_MFMT'] == '33,%i' % mt ]
                        if not cov: 
                            print ( "MT = %i (%s) covariance not present in the file" % ( mt, reaction ) )
                            rawEndfUnc.append( None )
                        else:
                            rawEndfUnc.append( endf_utils.getUncertainty( cov[0], rawEndfData[-1] ) )
        
                    thisSetStyle=plot_utils.getThisSetStyle(plotStyle,'evaluation',endf)
                    endfData.append( \
                        DataSet2d( \
                            rawEndfData[-1], \
                            uncertainty = rawEndfUnc[-1],\
                            legend        = thisSetStyle['legend'], \
                            lineWidth     = thisSetStyle['lineWidth'], \
                            lineStyle     = thisSetStyle['lineStyle'], \
                            color         = thisSetStyle['lineColor'], \
                            errorbarColor = thisSetStyle['errorColor'] ) )   

                #except Exception, err: print "Adding data from "+endf+" failed with error "+str(err)     

            # and plot all the rest, they don't get covariances...
            if showparts :
                for i, endfRxn in enumerate( reactionList ):
                    if i == mainReaction: continue # already did it
                    csData = endf_utils.getPointwiseCrossSection( endfRxn )
                    if csData != None:    rawEndfData.append( csData )
                    else:                 continue
                    if args.simplelegend:       endfData.append( DataSet2d( rawEndfData[-1] ) )
                    elif len( gndMap ) == 1:    endfData.append( DataSet2d( rawEndfData[-1], legend = str( endfRxn ) ) )
                    else:                       endfData.append( DataSet2d( rawEndfData[-1], legend = endf+': '+ str( endfRxn ) ) )
            print

#        except Exception, err: print "WARNING: could not add endf data in "+endf+" because got error: "+str(err)
        
    # Get the C4 data for plotting
    c4Data = []
    if ( args.c4File != None ): 
        c4Data = getC4DataSets( args.c4File, mt, 3, target, projectile, plotSyle=plotStyle )

    # Get the EXFOR data for plotting
    exforData = []
    if ( not nox4 ) and ( not reaction is None ): 
        sym, A, m = gnd_utils.getSymAFromGNDName( target )
        exforData = getEXFORSets( sym, A, reaction = projectile+','+reaction, quantity = "SIG", nox4evals = nox4evals, nox4legend = nox4legend, plotSyle=plotStyle  )

    # Actually make the plot
    if endfData + exforData + c4Data != [] and xyData != [] and xydyData != [] and xdxydyData != []:
 
        generatePlot( observable = 'crossSection', \
                      dataSets = endfData + exforData + c4Data, \
                      xyData = xyData, \
                      xydyData = xydyData, \
                      xdxydyData = xdxydyData, \
                      plotStyle = plotStyle, \
                      suggestTitle = getSuggestTitle( target, projectile, reaction, mt ), \
                      suggestXLog = len( set(mtList).intersection( [ 1, 2, 18, 102 ] + range(501, 574) ) ) != 0,\
                      suggestYLog = len( set(mtList).intersection( [ 1, 2, 18, 102 ] + range(501, 574) ) ) != 0 )
    


def makeCrossSectionUncertaintyPlot( gndMap={}, xyData={}, xydyData={}, xdxydyData={}, mt=0, projectile='n', target='1H', outFile=None, plotStyle={} ):
    
    raise NotImplementedError("Must refactor")

    # Assemble the cross section data into something matplotlib can plot
    if not cov: raise TypeError( "File has no covariance" ) 
    ratioDomain = ( max(rawEndfUnc[-1].domainMin(), rawEndfData[-1].domainMin()),
            min(rawEndfUnc[-1].domainMax(), rawEndfData[-1].domainMax()) )
    ratio = rawEndfUnc[-1].domainSlice(*ratioDomain) / rawEndfData[-1].domainSlice(*ratioDomain)
    endfData.append( \
        DataSet2d( ratio, \
            uncertainty = None, \
            legend = getStyle( plotStyle, 'evaluations', endf, 'legend', default=endf + ': (' + projectile+','+reaction.lower() + ')'), \
            lineWidth=getStyle( plotStyle, 'evaluations', endf, 'lineWidth' ), \
            lineStyle=getStyle( plotStyle, 'evaluations', endf, 'lineStyle' ), \
            color=getStyle( plotStyle, 'evaluations', endf, 'lineColor' ), \
            errorbarColor=getStyle( plotStyle, 'evaluations', endf, 'errorColor' ) ) )       


def makeCrossSectionIntegralsPlot( gndMap={}, xyData={}, xydyData={}, xdxydyData={}, mt=0, projectile='n', target='1H', outFile=None, plotStyle={} ): 
    
    raise NotImplementedError()


#---------------------------------------------------
# special atomic data plots
#---------------------------------------------------
def makeEnergyTransferPlot( gndMap={}, xyData={}, xydyData={}, xdxydyData={}, mt=0, projectile='e-', target='1H', product='e-', outFile=None, plotStyle={} ): 
    
    raise NotImplementedError("Ask Bret how to get at the energy transfer for LAW=8 data in MF=26, MT=528 or 527")
    
    # Get the ENDF data for plotting
    endfData=[]
    for endf in gndMap:
        print utils.smallBanner( "Preparing data for "+endf )
        reactionSuite, covarianceSuite = gndMap[endf]

        # collect complete list of matching reactions
        reactionList = []
        for MT in mtList: 
            try: reactionList += endf_utils.getReactions( reactionSuite, mt )
            except: pass
        print "Exit Channels:", ', '.join( map( str, reactionList ) )

        # Extract the cross section data from the endf file
        dist = reactionList[0].outputChannel.getProductWithName(product).distributions.energyTransfer.getNativeData().toPointwise_withLinearXYs(1e-8, 0, 1e-8)
        thisSetStyle=plot_utils.getThisSetStyle(plotStyle,'evaluation',endf)
        endfData.append( \
                         DataSet2d( \
                                    dist, \
                                    legend        = thisSetStyle['legend'], \
                                    lineWidth     = thisSetStyle['lineWidth'], \
                                    lineStyle     = thisSetStyle['lineStyle'], \
                                    color         = thisSetStyle['lineColor'], \
                                    errorbarColor = thisSetStyle['errorColor'] ) )   

    # Actually make the plot
    if endfData != [] and xyData != [] and xydyData != [] and xdxydyData != []:
 
        generatePlot( observable = 'energyTransfer', \
                      dataSets = endfData, \
                      xyData = xyData, \
                      xydyData = xydyData, \
                      xdxydyData = xdxydyData, \
                      plotStyle = plotStyle, \
                      suggestTitle = projectile + '+' + target + '->' + str(reactionList[0]), \
                      suggestXLog = len( set(mtList).intersection( range(501, 574) ) ) != 0,\
                      suggestYLog = len( set(mtList).intersection( range(501, 574) ) ) != 0 )


def makeFormFactorPlot( gndMap={}, xyData={}, xydyData={}, xdxydyData={}, mt=0, projectile='gamma', target='1H', product='gamma', outFile=None, plotStyle={} ): 

    # collect complete list of matching reactions
    if mt not in [ 502, 504 ]: mtList = [ 502, 504 ]
    else: mtList = [mt]
    mtPrefix = {502:"coherent", 504:"incoherent"}
    mtLineStyle = {502:"-", 504:":"}

    # Get the ENDF data for plotting
    endfData=[]
    for endf in gndMap:
        print utils.smallBanner( "Preparing data for "+endf )
        reactionSuite, covarianceSuite = gndMap[endf]

        # Extract the cross section data from the endf file
        for MT in mtList:
            dist = endf_utils.getReactions( reactionSuite, MT )[0].outputChannel.getProductWithName(product).distributions.getNativeData()
            if hasattr(dist,'formFactor'): dist=dist.formFactor
            dist = dist.toPointwise_withLinearXYs(1e-8, 0, 1e-8)
            thisSetStyle=plot_utils.getThisSetStyle(plotStyle,'evaluation',endf)
            endfData.append( \
                         DataSet2d( \
                                    dist, \
                                    legend        = thisSetStyle['legend']+" ("+mtPrefix[MT]+")", \
                                    lineWidth     = thisSetStyle['lineWidth'], \
                                    lineStyle     = mtLineStyle[MT],\
                                    #thisSetStyle['lineStyle'], \
                                    color         = thisSetStyle['lineColor'], \
                                    errorbarColor = thisSetStyle['errorColor'] ) )   

    # Actually make the plot
    if endfData != [] and xyData != [] and xydyData != [] and xdxydyData != []:
 
        generatePlot( observable = 'formFactor', \
                      dataSets = endfData, \
                      xyData = xyData, \
                      xydyData = xydyData, \
                      xdxydyData = xdxydyData, \
                      plotStyle = plotStyle, \
                      suggestTitle = projectile + '+' + target + ' elastic scattering form factors', \
                      suggestXLog = len( set(mtList).intersection( range(501, 574) ) ) != 0,\
                      suggestYLog = len( set(mtList).intersection( range(501, 574) ) ) != 0 )


def makeAnomolousScatteringFactorPlot( gndMap={}, xyData={}, xydyData={}, xdxydyData={}, mt=502, projectile='gamma', target='1H', product='gamma', outFile=None, plotStyle={} ):
    
    # Get the ENDF data for plotting
    endfData=[]
    for endf in gndMap:
        print utils.smallBanner( "Preparing data for "+endf )
        reactionSuite, covarianceSuite = gndMap[endf]

        # Extract the cross section data from the endf file
        dist = endf_utils.getReactions( reactionSuite, mt )[0].outputChannel.getProductWithName(product).distributions.getNativeData()
        for x in [ 'anomalousScatteringFactor_imaginaryPart', 'anomalousScatteringFactor_realPart' ]:
            if not hasattr(dist,x): continue
            thisSetStyle=plot_utils.getThisSetStyle(plotStyle,'evaluation',endf)
            endfData.append( \
                             DataSet2d( \
                                        getattr(dist,x), \
                                        legend        = thisSetStyle['legend']+', '+x.split('_')[-1], \
                                        lineWidth     = thisSetStyle['lineWidth'], \
                                        lineStyle     = {'imaginaryPart':':', 'realPart':'-'}[x.split('_')[-1]], \
                                        color         = thisSetStyle['lineColor'], \
                                        errorbarColor = thisSetStyle['errorColor'] ) )   

    # Actually make the plot
    if endfData != [] and xyData != [] and xydyData != [] and xdxydyData != []:
 
        generatePlot( observable = 'anomolousScatteringFactor', \
                      dataSets = endfData, \
                      xyData = xyData, \
                      xydyData = xydyData, \
                      xdxydyData = xdxydyData, \
                      plotStyle = plotStyle, \
                      suggestTitle = projectile + '+' + target  + ' coherent elastic anomolous scattering factors', \
                      suggestXLog = len( set(mtList).intersection( range(501, 574) ) ) != 0,\
                      suggestYLog = len( set(mtList).intersection( range(501, 574) ) ) != 0 )



#---------------------------------------------------
# multiplicity plots
#---------------------------------------------------
def makeMultiplicityPlot( gndMap={}, xyData={}, xydyData={}, xdxydyData={}, mt=0, projectile='n', target='1H', product='n', outFile=None, plotStyle={} ): 

    raise NotImplementedError()


#---------------------------------------------------
# energy/momentum deposition plots
#---------------------------------------------------
def makeEnergyDepositPlot( gndMap={}, xyData={}, xydyData={}, xdxydyData={}, mt=0, projectile='n', target='1H', product='n', outFile=None, plotStyle={}, observable=None ): 
    from fudge.gnd.productData import energyDeposition
    import fudge.core.math.xData.XYs as XYs

    endfData = []
    
    # Get the ENDF data for plotting
    for endf in gndMap:
        print utils.smallBanner( "Preparing data for "+endf )
        reactionSuite, covarianceSuite = gndMap[endf]
        
        # Set up MT list.  
        # This is the list of all MT that we will search for data to sum into the plot.  
        # Usually, it is one element long, but there are some special cases: 
        #     summed channels that were not specified in the ENDF file and lumped channels....    
        if not endf_utils.hasMT( reactionSuite, mt ): continue

        # We may need resonance reconstruction
        if mt in [ 1, 2, 18, 102 ] and args.doResonanceReconstruction:  
            reactionSuite.reconstructResonances(styleName='reconstructed')
        
        # We need to compute energy depositions
        reactionSuite.calculateDepositionData({'verbosity':0})
        
        # Get the right reaction
        reaction = endf_utils.getReactions( reactionSuite, mt )[0]       
        print "Exit Channel:", str(reaction)
        
        # Compute the energy deposits and available energy
        Q = reaction.outputChannel.Q.toPointwise_withLinearXYs( 1e-8, 0 )
        energyDep = [\
                     [prod.getLabel(), prod.data[ energyDeposition.component.genre ].getNativeData()] \
                     for prod in reaction.outputChannel.particles if energyDeposition.component.genre in prod.data ]
        if energyDep:
            totalEDep = energyDep[0][1].copy()
            for idx in range(1,len(energyDep)):
                if totalEDep.domain() != energyDep[idx][1].domain():
                    totalEDep, energyDep[idx][1] = totalEDep.mutualify(1e-8,0,0, energyDep[idx][1], 1e-8,0,0)
                totalEDep += energyDep[idx][1]
        else: totalEDep = []
        expectedTotalEnergy = XYs.XYs( axes_=Q.axes, data=[ [ x[0],x[1]+x[0] ] for x in Q.copyDataToXYs() ], accuracy=1e-8 )

        # Add sets to plot list       
        thisSetStyle=plot_utils.getThisSetStyle(plotStyle,'evaluation',endf,verbose=False)
        if observable == "energyDeposit":

            for i,p in enumerate(energyDep):
                endfData.append( DataSet2d(\
                                        p[1], \
                                        legend = 'Energy for product %s (%s)' % (str(p[0]),thisSetStyle['legend']), \
                                        lineWidth     = max(thisSetStyle['lineWidth']-2,2), \
                                        lineStyle     = plot_utils.getLineStyle(i+1), \
                                        color         = thisSetStyle['lineColor'] ) )
            endfData.append( DataSet2d( \
                                        totalEDep,\
                                        legend = 'Actual total energy (%s)' % thisSetStyle['legend'], \
                                        lineWidth     = thisSetStyle['lineWidth']+2, \
                                        lineStyle     = thisSetStyle['lineStyle'], \
                                        color         = thisSetStyle['lineColor'] ) )         
            endfData.append( DataSet2d( \
                                        expectedTotalEnergy,\
                                        legend = 'Expected total energy (%s)' % thisSetStyle['legend'], \
                                        lineWidth     = thisSetStyle['lineWidth']+5, \
                                        lineStyle     = thisSetStyle['lineStyle'], \
                                        color         = '#cccccc' ) )
      
        else: raise NotImplementedError( observable )
    
    # Actually make the plot
    if endfData != [] and xyData != [] and xydyData != [] and xdxydyData != []:
 
        generatePlot( observable = observable, \
                      dataSets = endfData, \
                      xyData = xyData, \
                      xydyData = xydyData, \
                      xdxydyData = xdxydyData, \
                      plotStyle = plotStyle, \
                      suggestTitle = getSuggestTitle( target, projectile, str(reaction), mt ), \
                      suggestXLog = len( set(mtList).intersection( [ 1, 2, 4, 18, 102, 103, 105, 106, 107 ] + range(501, 574) ) ) != 0,\
                      suggestYLog = len( set(mtList).intersection( [ 1, 2, 4, 18, 102, 103, 105, 106, 107 ] + range(501, 574) ) ) != 0,\
                      suggestFrame = 'lab' )


def makeMomentumDepositPlot( gndMap={}, xyData={}, xydyData={}, xdxydyData={}, mt=0, projectile='n', target='1H', product='n', outFile=None, plotStyle={}, observable=None ): 
    from fudge.gnd.productData import momentumDeposition
    import fudge.core.math.xData.XYs as XYs
    import math
    
    projectileMass = None
    endfData = []
    
    # Get the ENDF data for plotting
    for endf in gndMap:
        print utils.smallBanner( "Preparing data for "+endf )
        reactionSuite, covarianceSuite = gndMap[endf]
        
        if projectileMass == None:projectileMass = gndMap[endf][0].projectile.getMass('eV/c/c')
        
        # Set up MT list.  
        # This is the list of all MT that we will search for data to sum into the plot.  
        # Usually, it is one element long, but there are some special cases: 
        #     summed channels that were not specified in the ENDF file and lumped channels....    
        if not endf_utils.hasMT( reactionSuite, mt ): continue

        # We may need resonance reconstruction
        if mt in [ 1, 2, 18, 102 ] and args.doResonanceReconstruction:  
            reactionSuite.reconstructResonances(styleName='reconstructed')
        
        # We need to compute energy depositions
        reactionSuite.calculateDepositionData({'verbosity':0})
        
        # Get the right reaction
        reaction = endf_utils.getReactions( reactionSuite, mt )[0]       
        print "Exit Channel:", str(reaction)
        
        # Compute the energy deposits and available energy
#        Q = reaction.outputChannel.Q.toPointwise_withLinearXYs( 1e-8, 0 )
        momDep = [\
                     [prod.getLabel(), prod.data[ momentumDeposition.component.genre ].getNativeData()] \
                     for prod in reaction.outputChannel.particles if momentumDeposition.component.genre in prod.data ]
        if momDep:
            totalMomDep = momDep[0][1].copy()
            for idx in range(1,len(momDep)):
                if totalMomDep.domain() != momDep[idx][1].domain():
                    totalMomDep, momDep[idx][1] = totalMomDep.mutualify(1e-8,0,0, momDep[idx][1], 1e-8,0,0)
                totalMomDep += momDep[idx][1]
        else: totalMomDep = []
        expectedTotalMom = XYs.XYs( axes_=totalMomDep.axes, data=[ [ x[0], math.sqrt(2.0*x[0]*projectileMass) ] for x in totalMomDep.copyDataToXYs() ], accuracy=1e-8 )

        # Add sets to plot list       
        thisSetStyle=plot_utils.getThisSetStyle(plotStyle,'evaluation',endf)
        if observable == "momentumDeposit":

            for i,p in enumerate(momDep):
                endfData.append( DataSet2d(\
                                        p[1], \
                                        legend = 'Forward momentum of product %s (%s)' % (str(p[0]),thisSetStyle['legend']), \
                                        lineWidth     = max(thisSetStyle['lineWidth']-2,2), \
                                        lineStyle     = plot_utils.getLineStyle(i+1), \
                                        color         = thisSetStyle['lineColor'] ) )
            endfData.append( DataSet2d( \
                                        totalMomDep,\
                                        legend = 'Actual total forward momentum (%s)' % thisSetStyle['legend'], \
                                        lineWidth     = thisSetStyle['lineWidth']+2, \
                                        lineStyle     = thisSetStyle['lineStyle'], \
                                        color         = thisSetStyle['lineColor'] ) )         
            endfData.append( DataSet2d( \
                                        expectedTotalMom,\
                                        legend = 'Expected total forward momentum (%s)' % thisSetStyle['legend'], \
                                        lineWidth     = thisSetStyle['lineWidth']+5, \
                                        lineStyle     = thisSetStyle['lineStyle'], \
                                        color         = '#cccccc' ) )
      
        else: raise NotImplementedError( observable )
    
    # Actually make the plot
    if endfData != [] and xyData != [] and xydyData != [] and xdxydyData != []:
 
        generatePlot( observable = observable, \
                      dataSets = endfData, \
                      xyData = xyData, \
                      xydyData = xydyData, \
                      xdxydyData = xdxydyData, \
                      plotStyle = plotStyle, \
                      suggestTitle = getSuggestTitle( target, projectile, str(reaction), mt ), \
                      suggestXLog = len( set(mtList).intersection( [ 1, 2, 4, 18, 102, 103, 105, 106, 107 ] + range(501, 574) ) ) != 0,\
                      suggestYLog = len( set(mtList).intersection( [ 1, 2, 4, 18, 102, 103, 105, 106, 107 ] + range(501, 574) ) ) != 0,\
                      suggestFrame = 'lab' )


def makeFissionEnergyReleasePlot( gndMap={}, xyData={}, xydyData={}, xdxydyData={}, mt=0, projectile='n', target='1H', product='n', outFile=None, plotStyle={} ): 

    raise NotImplementedError()

#---------------------------------------------------
# angular distribution plots
#---------------------------------------------------
def makeAngDistMubarPlot( gndMap={}, xyData={}, xydyData={}, xdxydyData={}, mt=0, projectile='n', target='1H', product='n', outFile=None, plotStyle={} ): 

    # Declare the main reaction type
    try: reaction = getEXFORRxn( mt )
    except KeyError: 
        if mt in range( 850, 871, 1 ): 
            print( "Skipping MT = "+ str( mt )+", it is lumped covariance" )
            return
        else: raise KeyError( "Unknown MT: "+str( mt ) )
    
   # Get the ENDF data for plotting
    endfData = []
    for endf in gndMap:
        print utils.smallBanner( "Preparing data for "+endf )
        reactionSuite, covarianceSuite = gndMap[endf]
        reactionList = endf_utils.getReactions( reactionSuite, mt )
        for r in reactionList:
            print "Retrieving data for",r
            csData = endf_utils.getPointwiseCrossSection( r )
            angDistData = endf_utils.getAngularDistribution( r, product )
            for a in angDistData:
                suggestFrame = a.productFrame
                legend = product+' given as '+a.name+" in "+a.productFrame+' frame ('+endf+')' 
                xUnit = a.axes[0].unit
                endfData.append( 
                                DataSet2d( 
                                          [ [ E, a.muAverageAtEnergy(E) ] for E in a.getEnergyArray() ], \
                                          xUnit=a.axes[0].unit, \
                                          yUnit='', \
                                          lineWidth=3, \
                                          legend=legend) )
        print
        
    # Actually make the plot
    if endfData != [] and xyData != [] and xydyData != [] and xdxydyData != []:
 
        generatePlot( observable = 'mubar', \
                      dataSets = endfData, \
                      xyData = xyData, \
                      xydyData = xydyData, \
                      xdxydyData = xdxydyData, \
                      plotStyle = plotStyle, \
                      suggestTitle = getSuggestTitle( target, projectile, reaction, mt ), \
                      suggestXLog = len( set(mtList).intersection( [ 1, 2, 18, 102 ] + range(501, 574) ) ) != 0,\
                      suggestYLog = len( set(mtList).intersection( [ 1, 2, 18, 102 ] + range(501, 574) ) ) != 0,\
                      suggestFrame=suggestFrame )
        
        
        
def makeAngDistLegendreMomentPlot( gndMap={}, xyData={}, xydyData={}, xdxydyData={}, mt=0, L=1, projectile='n', target='1H', product='n', outFile=None, plotStyle={} ): 
    
   # Get the ENDF data for plotting
    endfData = []
    for endf in gndMap:
        print utils.smallBanner( "Preparing data for "+endf )
        reactionSuite, covarianceSuite = gndMap[endf]
        reactionList = endf_utils.getReactions( reactionSuite, mt )
        for r in reactionList:
            print "Retrieving data for",r
            angDistData = endf_utils.getAngularDistribution( r, product )
            for a in angDistData:
                suggestFrame = a.productFrame
                legend = 'L=%i Legendre moment for %s given in %s frame (%s)' % (L, product, a.productFrame, endf) 
                xUnit = a.axes[0].unit
                endfData.append( 
                                DataSet2d( 
                                          a.getCoefficientSafely( L ), \
                                          xUnit=a.axes[0].unit, \
                                          yUnit='', \
                                          lineWidth=3, \
                                          legend=legend) )
        print
        
    # Actually make the plot
    if endfData != [] and xyData != [] and xydyData != [] and xdxydyData != []:
 
        generatePlot( observable = 'LegendreMoment', \
                      dataSets = endfData, \
                      xyData = xyData, \
                      xydyData = xydyData, \
                      xdxydyData = xdxydyData, \
                      plotStyle = plotStyle, \
                      suggestTitle = getSuggestTitle( target, projectile, reaction, mt ), \
                      suggestXLog = len( set(mtList).intersection( [ 1, 2, 18, 102 ] + range(501, 574) ) ) != 0,\
                      suggestYLog = len( set(mtList).intersection( [ 1, 2, 18, 102 ] + range(501, 574) ) ) != 0,\
                      suggestFrame=suggestFrame )



#---------------------------------------------------
# main routine!
#---------------------------------------------------
if __name__ == "__main__": 

    args = parseArgs()

    # Declare the various particles involved in the plot
    target = args.target
    projectile = args.projectile
    product = args.product

    # Read in the plot style information
    userDefs = plot_utils.readUserStyles( args.style )
    if False:  # for debugging
        import json
        print json.dumps(userDefs,indent=2)
        exit()
    
    # Read the ENDF (& other user provided) data
    xyData = {}
    xydyData = {}
    xdxydyData = {}
    gndMap = collections.OrderedDict()
    mtMap = {}
    print utils.bigBanner( "Reading evaluation files" )
    for endf in args.endf:
        if endf == "None": continue
        print utils.smallBanner( "Reading "+endf )

        # Try plain XY data
        if endf.endswith( '.xy.dat'):       xyData[endf]   = io.readXYData( endf )

        # Try plain XYdY data
        elif endf.endswith( '.xydy.dat'):   xydyData[endf] = io.readXYdYData( endf )

        # Try plain XdXYdY data
        elif endf.endswith( '.xdxydy.dat'):   xdxydyData[endf] = io.readXdXYdYData( endf )

        # Must be an evaluation
        else:                               gndMap[endf]   = io.readEvaluation( endf )
                                
        if not endf in gndMap: continue
                
        mtMap[endf] = endf_utils.getEvaluationMTs( gndMap[endf][0] )


    # MT is specified by the user.  Just plot that.
    if args.mt != 0: mtList = [ args.mt ]

    # MT not specified, so assemble a list of mt's.
    # For discrete level excitations, we need to add the summed up version.
    # Also, put them in order so the first is always the sum.
    else: 
        mtList = []
        for mt in mtMap[ args.endf[0] ]: # MT Map has all the MT's in every evaluation.  We just want highlights.
            if   mt >= 50  and mt <= 91 :
                if 4 not in mtList:   mtList.append( 4 )   # (*,n)
            elif mt >= 600 and mt <= 649:
                if 103 not in mtList: mtList.append( 103 ) # (*,p)
            elif mt >= 650 and mt <= 699:
                if 104 not in mtList: mtList.append( 104 ) # (*,d)
            elif mt >= 700 and mt <= 749:
                if 105 not in mtList: mtList.append( 105 ) # (*,t)
            elif mt >= 750 and mt <= 799:
                if 106 not in mtList: mtList.append( 106 ) # (*,3He)
            elif mt >= 800 and mt <= 849:
                if 107 not in mtList: mtList.append( 107 ) # (*,a)
            elif mt >= 875 and mt <= 891:
                if 16 not in mtList:  mtList.append( 16 )  # (*,2n)
            elif mt in [ 19, 20, 21, 38 ]:
                if 18 not in mtList:  mtList.append( 18 )  # (*,f)
            elif mt not in mtList: mtList.append( mt )
            try: reaction = getEXFORRxn( mt )
            except KeyError: 
                if mt in range( 850, 871, 1 ): print( "Got lumped covariance for MT = "+ str( mt ) )
                else: raise KeyError( "Unknown MT: "+str( mt ) )
        mtList.sort()
      
      
    # Determine the names of the particles involved
    print utils.bigBanner( "Plot details" )
    if target == None:     target = str(gndMap.items()[0][1][0].target.name)
    if projectile == None: projectile = str(gndMap.items()[0][1][0].projectile.name)
    if product == None and args.observable not in ['crossSection','crossSectionIntegrals',"energyDeposit",'energyBalance',"momentumDeposit",'momentumBalance','fissionEnergyRelease']: 
        if args.observable in ["formFactor",'anomolousScatteringFactor']: 
            print "Assuming product is a 'gamma' for observable %s"%args.observable
            product = 'gamma'
        elif args.observable in ["energyTransfer"]:
            print "Assuming product is an 'e-' for observable %s"%args.observable
            product = 'e-'
        elif args.observable in ["nubar"]:
            print "Assuming product is an 'n' for observable %s"%args.observable
            product = 'n'
        else:
            raise ValueError( "For observable='"+str(args.observable)+"', need to declare a product using --product command line argument" )
    print "Projectile is "+projectile
    print "Target is "+target
    print "List of MT's to plot:",mtList

    
    # Determine the observable type
    if product ==None: print "Observable is "+args.observable
    else: print "Observable is "+args.observable+' for '+str(product)
    print

    # Loop over reactions
    for mt in mtList:    
    
        # What is the proper name of the reaction?
        reaction = projectile.capitalize() +','+ getEXFORRxn( mt )
        
        # Compute the output filename
        if args.outFile != None and not args.outFile.endswith('.png'):
            outFile = args.outFile+'_mt'+str(mt)+'_mf'+str(args.mf)+'.png'
        else: outFile = None
        
        # Determine the kind of plot and set up the axes
        print utils.bigBanner( "Generating plots for "+target+"("+reaction.lower()+')' )
        if args.observable == "mubar":                      
            makeAngDistMubarPlot( \
                gndMap, xyData, xydyData, xdxydyData, \
                mt=mt, \
                projectile=projectile, target=target, product=product, \
                outFile=outFile, \
                plotStyle=userDefs )
        elif args.observable == "formFactor":               
            if mt not in [502,504]: 
                print "in ENDF, formFactor data only exists for (in)coherent photon-atom scattering (MT=502, 504), so only showing MT=502 and 504"
            makeFormFactorPlot(
                gndMap, xyData, xydyData, xdxydyData, \
                projectile=projectile, target=target, product=product, \
                outFile=outFile, 
                plotStyle=userDefs )
            break
        elif args.observable == 'anomolousScatteringFactor':
            if mt != 502: 
                print "in ENDF, anomolousScatteringFactor data only exists for coherent photon-atom scattering (MT=502), so only showing MT=502"
            makeAnomolousScatteringFactorPlot(
                gndMap, xyData, xydyData, xdxydyData, \
                projectile=projectile, target=target, product=product, \
                outFile=outFile, 
                plotStyle=userDefs )
            break
        elif args.observable == 'energyTransfer':
            if mt not in [527, 528]: 
                print "in ENDF, energyTransfer data only exists for bremstrahlung reactions (MT=527, 528)"
                continue
            makeEnergyTransferPlot(                
                gndMap, xyData, xydyData, xdxydyData, \
                mt=mt, \
                projectile=projectile, target=target, product=product, \
                outFile=outFile, 
                plotStyle=userDefs )
        elif args.observable in ["energyDeposit",'energyBalance'] :            
            makeEnergyDepositPlot( \
                gndMap, xyData, xydyData, xdxydyData, \
                mt=mt, \
                projectile=projectile, target=target, product=product, \
                outFile=outFile, \
                plotStyle=userDefs,\
                observable =args.observable )
        elif args.observable == "fissionEnergyRelease":                    
            makeFissionEnergyReleasePlot( \
                gndMap, xyData, xydyData, xdxydyData, \
                mt=mt, \
                projectile=projectile, target=target, product=product, \
                outFile=outFile, \
                plotStyle=userDefs )
        elif args.observable == "nubar":                    
            makeMultiplicityPlot( \
                gndMap, xyData, xydyData, xdxydyData, \
                mt=mt, \
                projectile=projectile, target=target, product=product, \
                outFile=outFile, \
                plotStyle=userDefs ) 
        elif args.observable == "crossSectionIntegrals":    
            makeCrossSectionIntegralsPlot( \
                gndMap, xyData, xydyData, xdxydyData, \
                mt=mt, \
                projectile=projectile, target=target, product=product, \
                outFile=outFile, \
                plotStyle=userDefs )
        elif args.observable == "crossSection":
            if  args.uncRatio:
                makeCrossSectionUncertaintyPlot( \
                    gndMap, xyData, xydyData, xdxydyData, \
                    mt=0, \
                    projectile='n', target='1H', \
                    outFile=None, \
                    plotStyle={} )
            else:
                makeCrossSectionPlot( \
                    gndMap, xyData, xydyData, xdxydyData, \
                    mt=mt, \
                    projectile=projectile, target=target, \
                    nounc=args.noUnc,\
                    nox4=args.noX4, \
                    showparts=args.showParts, \
                    nox4evals=not args.showX4Evals, \
                    nox4legend=args.noX4Legend, \
                    doResonanceReconstruction=args.doResonanceReconstruction, \
                    heatToTemp=None, \
                    outFile=outFile, \
                    plotStyle=userDefs )
        elif args.observable in ["momentumDeposit", 'momentumBalance']:          
            makeMomentumDepositPlot(  \
                gndMap, xyData, xydyData, xdxydyData, \
                mt=mt, \
                projectile=projectile, target=target, product=product, \
                outFile=outFile, \
                plotStyle=userDefs,\
                observable =args.observable )
        elif args.observable == "LegendreMoment":           
            makeAngDistLegendreMomentPlot( \
                gndMap, xyData, xydyData, xdxydyData, \
                mt=mt, \
                L=args.L, \
                projectile=projectile, target=target, product=product, \
                outFile=outFile, \
                plotStyle=userDefs )
        else: raise ValueError( "Unsupported plot observable: %s" % args.observable )
        
            



