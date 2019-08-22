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

import collections, os, fudge, json, fnmatch
from fudge.vis.matplotlib import plot2d, DataSet2d, DataSet3d, defaultSymbols, nSymbols, defaultColors, nColors
from fudge.legacy.converting import endfFileToGND
from fudge.core.utilities import argparse
from fudge.core.utilities.brb import uniquify
from fudge.gnd.covariances.base import covarianceMatrix
from fudge.gnd.covariances.mixed import mixedForm
from fudge.particles import nuclear

 
#---------------------------------------------------
# Set up the command line parser
#---------------------------------------------------
def parseArgs():
    parser = argparse.ArgumentParser(description='Plot nuclear data from an ENDF or GND file')
    parser.add_argument('mt',       metavar='mt', type=int, help='MT of the cross section to plot.  If set to 0, will try to make all plots for all open channels (requires -o option too)' )
    parser.add_argument('endf',     metavar='endf', type=str, nargs='+', help='ENDF file(s) whose cross section you want to plot.  Use "None" for no input file.' )
    parser.add_argument('-o',       dest='outFile', default=None, type=str, help='Output file for plot (disables interactive plotting)' )
    parser.add_argument('--mf',     type=int, default=3, help='MF of original data to plot [default is 3, cross section data, use 4 for angular, 5 for energy and 6 for energy-angle distributions and 12 for multiplicity]' )
    
    parser.add_argument('--c4File',         default=None,  type=str, help="Optional C4 file to pull data from instead of (or in addition to) EXFOR" )
    parser.add_argument('--nox4',           default=False, action='store_true', help='Do not plot EXFOR data' )
    parser.add_argument('--showx4evals',    dest='nox4evals', default=True, action='store_false', help='Plot evaluations found in EXFOR library (with "V" SUBENT)' )
    parser.add_argument('--nox4legend',     default=False, action='store_true', help='Do not put legend for EXFOR data on plot' )
    parser.add_argument('--showx4legend',   default=True,  action='store_true', help='Put legend for EXFOR data on plot, even if it would be automatically supressed' )
    
    parser.add_argument('--uncratio',       default=False, action='store_true', help='Plot a ratio of the uncertainty over the data' )
    parser.add_argument('--angMubar',       default=False, action='store_true', help='Make plot of mubar for an angular distribution' )
    parser.add_argument('--enableRRAngDist',default=False, action='store_true', help='Reconstruct the angular distribution from the resolved resonance parameters' )
    parser.add_argument('--angContour',     default=False, action='store_true', help='Make contour plot of an angular distribution' )
    parser.add_argument('--Pzbar',          default=False, action='store_true', help='Make average forward scattering momentum <pz> plot from angular distribution' )
    parser.add_argument('--enContour',      default=False, action='store_true', help='Make contour plot of an energy distribution' )
    parser.add_argument('--Ebar',           default=False, action='store_true', help="Make average emitted/deposited energy <E'> plot of an energy distribution" )
    parser.add_argument('--nounc',          default=False, action='store_true', help='Do not plot uncertainties' )
    parser.add_argument('--simplelegend',   default=False, action='store_true', help='Do not put entries for the components of requested MT in the plot legend' )
    parser.add_argument('--showparts',      default=False, action='store_true', help='Show all of the channels that comprise the one requested (if applicable)' )
    
    parser.add_argument('--logx', default=None, action='store_true', help='Make x scale logrithmic (Default is linear except for MT in [1,2,18,102])' )
    parser.add_argument('--logy', default=None, action='store_true', help='Make y scale logrithmic (Default is linear except for MT in [1,2,18,102])' )
    parser.add_argument('--logz', default=None, action='store_true', help='Make z scale logrithmic' )
    
    parser.add_argument('--nologx', dest='logx', default=None, action='store_false', help='Make x scale linear (Default is linear except for MT in [1,2,18,102])' )
    parser.add_argument('--nology', dest='logy', default=None, action='store_false', help='Make y scale linear (Default is linear except for MT in [1,2,18,102])' )
    parser.add_argument('--nologz', dest='logz', default=None, action='store_false', help='Make z scale linear' )
    
    parser.add_argument('--legendx', default=0.05, type=float, help='X coordinate of upper left corner of plot legend (Default 0.05)' ) 
    parser.add_argument('--legendy', default=0.95, type=float, help='Y coordinate of upper left corner of plot legend (Default 0.95)' ) 
    
    parser.add_argument('--xUnit', default=None, type=str, help='Use units of xUnit for the x axis (Defaults to eV for MF=3)' )
    #parser.add_argument('--xMin', default=None, type=float, help='Min. x in plots, in units of xUnit.  Setting this turns off autoscaling for xmin.' )
    #parser.add_argument('--xMax', default=None, type=float, help='Max. x in plots, in units of xUnit.  Setting this turns off autoscaling for xmax.' )
    
    parser.add_argument('--yUnit', default=None, type=str, help='Use units of yUnit for the y axis (Defaults to b for MF=3)' )
    #parser.add_argument('--yMin', default=None, type=float, help='Min. y in plots, in units of yUnit.  Setting this turns off autoscaling for ymin.' )
    #parser.add_argument('--yMax', default=None, type=float, help='Max. y in plots, in units of yUnit.  Setting this turns off autoscaling for ymax.' )
    
    parser.add_argument('--zUnit', default=None, type=str, help='Use units of zUnit for the z axis' )
    #parser.add_argument('--zMin', default=None, type=float, help='Min. z in plots, in units of zUnit.  Setting this turns off autoscaling for zmin.' )
    #parser.add_argument('--zMax', default=None, type=float, help='Min. z in plots, in units of zUnit.  Setting this turns off autoscaling for zmax.' )
    
    parser.add_argument('--noReconstruct', dest='doResonanceReconstruction', default=True, action='store_false', help="Don't reconstruct resonances" )
    parser.add_argument('--target',     default=None, type=str, help="The target nucleus, given in GND notation, e.g. Pu239 (Default is None which means to take it from the first ENDF file)" )
    parser.add_argument('--projectile', default=None, type=str, help="The projectile particles, given in GND notation, e.g. n (Default is None which means to take it from the first ENDF file)" )
    parser.add_argument('--product',    default=None, type=str, help="The product particle of interest, given in GND notation, e.g. n (Default is None which means to take the first emitted particle for this MT)" )
    parser.add_argument('--style',      default=None, type=str, help='JSON file with plot style overrides, see "plot_defaults.json" in source distribution for examples')
    return parser.parse_args()


#---------------------------------------------------
# utilities
#---------------------------------------------------
observableNameMap = { 3:"Cross Section ( Sigma(E) )", 23:"Cross Section ( Sigma(E) )", 4:"Angular Distribution ( dSigma(E)/d(cos Theta) )", 5:"Energy Distribution ( dSigma(E)/d(E') )", 6:"Energy-Angle Distribution ( dSigma(E)/d(E')d(cos Theta) )" }

ENDFStandardsReactions = ( ('H1',2), ('C0',2), ('AU197',102), ('U235',18), ('U238',18), ('U238',102), ('PU239',18) )

def isStandards( projectile, target, mt ):
    if projectile != 'n': return False
    if mt not in [ 2, 102, 18 ]: return False
    if target not in [ 'H1','C0','AU197','U235', 'U238','PU239' ]: return False
    for x in ENDFStandardsReactions:
        if ( target, mt ) == x: return True
    return False

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
    return { 1:'TOT', 2:'EL', 3:'NON', 4:'N', 5:'X', 11:'2N+D', 16:'2N', 17:'3N', 18:'F', 22:'N+A', 23:'N+3A', 24:'2N+A', 25:'3N+A', 27:'ABS', 28:'N+P', 29:'N+2A', 30:'2N+2A', 32:'N+D', 33:'N+T', 34:'N+HE3', 35:'N+D+2A', 36:'N+T+2A', 37:'4N', 41:'2N+P', 42:'3N+P', 44:'N+2P', 45:'N+P+A', 102:'G', 103:'P', 104:'D', 105:'T', 106:'HE3', 107:'A', 108:'2A', 109:'3A', 111:'2P', 112:'P+A', 113:'T+2A', 114:'D+2A', 115:'P+D', 116:'P+T', 117:'D+A' }.get( MT, "MT="+str(MT))


# Plot symbol & color
def getPlotSymbol( i ): return defaultSymbols[ i % nSymbols ]
def getPlotColor( i ): return defaultColors[ i % nColors ]


# Big banner, with flowerbox
def bigBanner( x ): 
    sx = x.split( '\n' )
    l = max( map( len, sx ) )+4
    return '\n'.join( [ '+'+l*'-'+'+' ] + [ '|'+y.center( l )+'|' for y in sx ] + [ '+'+l*'-'+'+' ] )


# Small banner, with wings
def smallBanner( x, wingsize=10 ): 
    return wingsize*'*'+' '+x.replace( '\n', '; ' )+' '+wingsize*'*'

# Crack the isotope name to get the ZA
def getZAFromGNDName( name ):
    sym, A, m = getSymAFromGNDName( name )   
    Z = nuclear.elementZFromSymbol( sym )
    return Z*1000+int(A)

# Crack the isotope name to get the A & symbol
def getSymAFromGNDName( name ):
    sym = ''
    A = ''
    if '_' in name: m = name.split('_')[1]
    else: m = None
    for c in name.split('_')[0]:
        if c.isalpha(): sym += c
        else: A += c
    if sym == 'n': return sym, 1, None
    if sym == 'g': return sym, 0, None
    if m =='natural': A = '0'
    return sym, A, m


#---------------------------------------------------
# Uncertainty management functions
#---------------------------------------------------

def getUncertainty( theCovariance, theData ):
    """ extract absolute uncertainty vector (ie, in units of barn) from covariance """
    theUncert = None
    if isinstance( theCovariance.forms['eval'], (covarianceMatrix, mixedForm ) ):
        theUncert = theCovariance.forms[ 'eval' ].getUncertaintyVector( theData, relative=False )
    if theUncert == None:
        print( 'Cannot plot uncertainty in any of these forms: ' +str( [ f for f in theCovariance.forms ] ) )
        return
    return theUncert 


#---------------------------------------------------
# stuff for interacting with various data sources
#---------------------------------------------------

def hasMT( reactionSuite, mt ): 
    return mt in getEvaluationMTs( reactionSuite )


def getPlottableReactions( reactionSuite ):
    plottables = list(reactionSuite.reactions) + list(reactionSuite.sums) + list(reactionSuite.fissionComponents) # + list(reactionSuite.productions)
    return [plt for plt in plottables if hasattr(plt, 'crossSection')]


def getEvaluationMTs( reactionSuite, mtFilter = None ):
    reactionList = getPlottableReactions( reactionSuite )
    if mtFilter == None: return [ r.ENDF_MT for r in reactionList ]
    else:                return [ r.ENDF_MT for r in reactionList if r.ENDF_MT in mtFilter ]


def getReactions( reactionSuite, MT ): 
    return [ r for r in getPlottableReactions( reactionSuite ) if r.ENDF_MT == MT ]


def getPointwiseCrossSection( reac ):
    """
    :param reac: reaction instance.  Look in reaction.crossSection to find
    :return:
    """
    if   'reconstructed'     in reac.crossSection.forms: return reac.crossSection['reconstructed']
    elif 'linearized'        in reac.crossSection.forms: return reac.crossSection['linearized']
    elif 'eval'              in reac.crossSection.forms:
        xsc = reac.crossSection.forms['eval']
        if xsc.moniker == 'resonancesWithBackground':   # in case --noReconstruct option is used:
            return xsc.tabulatedData.toPointwise_withLinearXYs(1e-8, 0)
        return xsc.toPointwise_withLinearXYs(1e-8, 0)
    return None


def getAngularDistribution( reac, product ):
    results = []
    for p in reac.outputChannel.particles:
        if p.label != product: print "    Skipping distributions for",p
        else:
            print "    Processing distributions for",p,"    ("+str(p.distributions)+")"
            nativeComponant = p.distributions.nativeData
            for component in p.distributions.components.values():
                compName = str( component ).split('/')[-1]
                if compName != "angular": 
                    print '        Skipping ' + compName + ' distribution component     ('+str(component)+')'
                    if compName == 'uncorrelated': print '*** Check uncorrelated data -- we might be able to plot it too'
                    continue
                else:
                    print '        Processing ' + compName + ' distribution component     ('+str(component)+')'
                    print '        Native form is',component.nativeData
                    results.append( component.forms[ component.nativeData ] )
                    break
    return results
    
    
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
    print bigBanner( "Preparing EXFOR data for " + sym+'-'+A +'(' +reaction.upper()+')'+', '+ quantity )
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
                if e.startswith( 'V' ) or not suppressEXFORLegend: theLegend = legend
                else: theLegend = '_noLegend_'
                exforData.append( DataSet2d( data = dat, uncertainty = unc, xUnit=ds[d].units[0], yUnit=barnsConverter( ds[d].units[1] ), legend = theLegend, lineStyle = ' ', symbol = getPlotSymbol( i ) ) )
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
    print smallBanner( 'Retrieving C4 data from: '+str( c4File ) )
    flist = filter( lambda x: x.MT == mt and x.MF == mf and x.target == getZAFromGNDName( target ) and x.projectile == getZAFromGNDName( projectile ), c4.readC4File( open( c4File ).readlines(), asPointList=True ) )
    i = -1
    dat = []
    unc = []
    c4Data = []
    theLegend = None
    lastSet = ( None, None )
    for p in flist:
        thisSet = ( p.reference, p.exforEntry+str(p.exforSubEntry).zfill(3) )
        if lastSet != thisSet: 
            if lastSet != ( None, None ): 
                print '    Set:',lastSet
                if mf == 3:   c4Data.append( DataSet2d( data = dat, uncertainty = unc, xUnit='eV', yUnit='b', legend = lastSet[0]+' ('+lastSet[1]+')', lineStyle = ' ', symbol = getPlotSymbol( i ) ) )
                elif mf == 4: c4Data.append( DataSet3d( data = dat, uncertainty = unc, xUnit='eV', yUnit='', zUnit='b', legend = lastSet[0]+' ('+lastSet[1]+')', lineStyle = ' ', symbol = getPlotSymbol( i ) ) )
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
    print '    Set:',lastSet
    if mf == 3:   c4Data.append( DataSet2d( data = dat, uncertainty = unc, xUnit='eV', yUnit='b', legend = lastSet[0]+' ('+lastSet[1]+')', lineStyle = ' ', symbol = getPlotSymbol( i ) ) )
    elif mf == 4: c4Data.append( DataSet3d( data = dat, uncertainty = unc, xUnit='eV', yUnit='', zUnit='b', legend = lastSet[0]+' ('+lastSet[1]+')', lineStyle = ' ', symbol = getPlotSymbol( i ) ) )
    print
    return c4Data

def readXYData(filename):
    if not os.path.exists( filename ): raise IOError( "XY data file %s not found"%filename)
    results = []
    for line in open(filename).readlines():
        try:
            results.append( map( float, line.split()[0:2] ) )
        except ValueError: pass
    return results 

def readXYdYData(filename):
    if not os.path.exists( filename ): raise IOError( "XYdY data file %s not found"%filename)
    results = []
    for line in open(filename).readlines():
        try:
            results.append( map( float, line.split()[0:3] ) )
        except ValueError: pass
    return results 

#---------------------------------------------------
# plot style resolution
#---------------------------------------------------
allowedDataKinds = [u'xydyCurves', u'xyCurves', u'C4Sets', u'evaluations', u'EXFORSets']
allowedStyleKeys = [ "filePattern","legend","symbol","lineStyle","lineColor","lineWidth","errorColor","xUnit","yUnit","zUnit","dataType"]
        
def getStyle( styleDict, dataKind, dataName, styleKey, default=None ):
    '''
    Use this to work out what style to use for a data set
    
    styleDict: should be main one constructed by __main__ routine
    dataKind: in allowedDataKinds list
    dataName: name to key off of when attempting to resolve the style information.  Usually the data file name, but could be the EXFOR/C4 entry.
    styleKey: the style thing itself, like 'symbol'.  used by matplotlib
    '''
    if dataKind not in styleDict: raise KeyError( "Cannot resolve data kind %s, so cannot look up style information"%dataKind )
    if styleKey not in allowedStyleKeys: raise ValueError( "%s is not a valid plot style element for data" )
    
    if default != None: return default
    
    result = None
 
    # Set the defaults first
    result = styleDict[dataKind]['default'].get(styleKey,result)
    
    # Now search through for wildcard overrides
    for key in styleDict[dataKind]:
        if key == 'default': continue # did it already
        if fnmatch.fnmatch(dataName,styleDict[dataKind][key]['filePattern']): 
            result = styleDict[dataKind][key].get(styleKey,result)
 
    # Final pass for specific overrides
    for key in styleDict[dataKind]: 
        if dataName == styleDict[dataKind][key]['filePattern']:
            result = styleDict[dataKind][dataName][key].get(styleKey,result)

    return result

#---------------------------------------------------
# cross section plots
#---------------------------------------------------
def makeCSPlot( gndMap, mt, projectile, target, xyData={}, xydyData={}, uncratio=False, nounc=False, nox4=False, showparts=False, nox4evals=True, nox4legend=False, xLog=None, xUnit='eV', yLog=None, yUnit='b', outFile=None, plotStyle={}, legendXY=(0.05, 0.95), doResonanceReconstruction=True, heatToTemp=None ):

    endfData=[]

    # Declare the main reaction type
    try: reaction = getEXFORRxn( mt )
    except KeyError: 
        if mt in range( 850, 871, 1 ): print( "Got lumped covariance for MT = "+ str( mt ) )
        else: raise KeyError( "Unknown MT: "+str( mt ) )
    
    # Set default units as needed
    if xUnit == None: xUnit = 'eV'
    if yUnit == None: 
        if uncratio: yUnit = ""
        else: yUnit = 'b'
        
    # Get the ENDF data for plotting
    for endf in gndMap:
        print smallBanner( "Preparing data for "+endf )
        reactionSuite, covarianceSuite = gndMap[endf]

        #try:
        if True:
        
            # Set up MT list.  
            # This is the list of all MT that we will search for data to sum into the plot.  
            # Usually, it is one element long, but there are some special cases: summed channels that were not specified in the ENDF file and lumped channels....
            mtList = [ mt ]
            if not hasMT( reactionSuite, mt ):
                if   mt == 3: pass # (n,non-el), use for catch-all "C=55" gammas
                elif mt == 5: pass  # (n,X), no sumrule specified.  Treat as plain old reaction. 
                elif mt == 1:          
                    # (n,tot)
                    # Need better solution here!!!!
                    print( "WARNING: fudge-2.0 considers (n,tot) to be redundant so will attempt to compute it from summing parts" )
                    if 1 in mtList: del( mtList[ mtList.index(1) ] )
                    mtList += getEvaluationMTs( reactionSuite, mtFilter = [ 2 ] + range( 5, 120, 1 ) )
                    mtList += getEvaluationMTs( reactionSuite, mtFilter = range( 51, 92, 1 ) )
                    mtList += getEvaluationMTs( reactionSuite, mtFilter = range( 600, 650, 1 ) )
                    mtList += getEvaluationMTs( reactionSuite, mtFilter = range( 650, 700, 1 ) )
                    mtList += getEvaluationMTs( reactionSuite, mtFilter = range( 700, 750, 1 ) )
                    mtList += getEvaluationMTs( reactionSuite, mtFilter = range( 750, 800, 1 ) )
                    mtList += getEvaluationMTs( reactionSuite, mtFilter = range( 800, 850, 1 ) )
                    mtList += getEvaluationMTs( reactionSuite, mtFilter = range( 875, 892, 1 ) )
                    if hasMT( reactionSuite, 18 ): mtList += [ 18 ]
                    else: mtList += getEvaluationMTs( reactionSuite, mtFilter = [ 19, 20, 21, 38  ] )
                elif mt == 4:          
                    # (n,n'), sum: 51-91
                    print "Summing (n,n') reactions MT=51-91"
                    if 4 in mtList: del( mtList[ mtList.index(4) ] )
                    mtList += getEvaluationMTs( reactionSuite, mtFilter = range( 51, 92, 1 ) )
                elif mt == 103:        
                    # (n,p), sum: 600-649
                    print "Summing (n,p) reactions MT=600-649"
                    if 103 in mtList: del( mtList[ mtList.index(103) ] )
                    mtList += getEvaluationMTs( reactionSuite, mtFilter = range( 600, 650, 1 ) )
                elif mt == 104:        
                    # (n,d), sum: 650-699
                    mtList += getEvaluationMTs( reactionSuite, mtFilter = range( 650, 700, 1 ) )
                elif mt == 105:        
                    # (n,t), sum: 700-749
                    mtList += getEvaluationMTs( reactionSuite, mtFilter = range( 700, 750, 1 ) )
                elif mt == 106:        
                    # (n,3He), sum: 750-799
                    mtList += getEvaluationMTs( reactionSuite, mtFilter = range( 750, 800, 1 ) )
                elif mt == 107:        
                    # (n,a), sum: 800-849
                    mtList += getEvaluationMTs( reactionSuite, mtFilter = range( 800, 850, 1 ) )
                elif mt == 16:         
                    # (n,2n), sum 875-891
                    mtList += getEvaluationMTs( reactionSuite, mtFilter = range( 875, 892, 1 ) )
                elif mt == 18:
                    # (n,f), sum: 19-21,38
                    mtList += getEvaluationMTs( reactionSuite, mtFilter = [ 19, 20, 21, 38 ] ) 
            elif mt in range( 851, 871, 1 ): 
                # lumped covariance MT's are 851-870
                # must extract the MT's from the covariance file
                for lump in covarianceSuite.lumpedChannels:
                    if lump.ENDF_MFMT == '33,'+str(mt): break
                mtList += getEvaluationMTs( reactionSuite, mtFilter = [ int( channel.attributes['ENDF_MFMT'].split(',')[1] ) for channel in lump.channels] )
                reaction = 'N,' + getEXFORRxn( mtList[-1] )
            else: pass # plain old reaction, do nothing special
            mtList = uniquify( mtList )

            # Do resonance reconstruction if needed
            if doResonanceReconstruction:
                for MT in mtList: 
                    rxnList = getReactions( reactionSuite, MT )
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
                    reactionList += getReactions( reactionSuite, MT )
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
                    csData = getPointwiseCrossSection( reactionList[0] )
                    if ( not nounc ) and mt not in [ 1, 2, 18, 102 ] :
                        try:
                            csData = csData.thicken( sectionSubdivideMax = 20 ) # only thicken cross sections that don't have enough points in them to resolve uncertainty steps
                        except Exception,err:
                            print "WARNING: "+str(err)
                    if mt != 18: # element 0 already was sum
                        for endfRxn in reactionList[1:]:
                            csData += getPointwiseCrossSection( endfRxn ) # this will need try/except block sooner or later
                    rawEndfData.append( csData )
                    
                    # Extract the uncertainty on the cross section
                    if nounc or covarianceSuite == None: rawEndfUnc.append( None )
                    else:
                        cov = [ c for c in covarianceSuite.sections if hasattr(c,'rowData') and c.rowData.attributes['ENDF_MFMT'] == '33,%i' % mt ]
                        if not cov: 
                            print ( "MT = %i (%s) covariance not present in the file" % ( mt, reaction ) )
                            rawEndfUnc.append( None )
                        else:
                            rawEndfUnc.append( getUncertainty( cov[0], rawEndfData[-1] ) )
        
                    # Assemble the cross section data into something matplotlib can plot
                    if uncratio: # do it as a ratio
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
                    else: # just a straight up cross section
                        endfData.append( \
                            DataSet2d( \
                                rawEndfData[-1], \
                                uncertainty = rawEndfUnc[-1],\
                                legend = getStyle( plotStyle, 'evaluations', endf, 'legend', default=endf + ': (' + projectile+','+reaction.lower() + ')' ), \
                                lineWidth=getStyle( plotStyle, 'evaluations', endf, 'lineWidth' ), \
                                lineStyle=getStyle( plotStyle, 'evaluations', endf, 'lineStyle' ), \
                                color=getStyle( plotStyle, 'evaluations', endf, 'lineColor' ), \
                                errorbarColor=getStyle( plotStyle, 'evaluations', endf, 'errorColor' ) ) )   

                #except Exception, err: print "Adding data from "+endf+" failed with error "+str(err)     

            # and plot all the rest, they don't get covariances...
            if showparts and not uncratio:
                for i, endfRxn in enumerate( reactionList ):
                    if i == mainReaction: continue # already did it
                    csData = getPointwiseCrossSection( endfRxn )
                    if csData != None:    rawEndfData.append( csData )
                    else:                 continue
                    if args.simplelegend:       endfData.append( DataSet2d( rawEndfData[-1] ) )
                    elif len( gndMap ) == 1:    endfData.append( DataSet2d( rawEndfData[-1], legend = str( endfRxn ) ) )
                    else:                       endfData.append( DataSet2d( rawEndfData[-1], legend = endf+': '+ str( endfRxn ) ) )
            print

#        except Exception, err: print "WARNING: could not add endf data in "+endf+" because got error: "+str(err)
        
    # Whether to do lin-lin or log-log based on the MT
    if xLog == None: 
        if set(mtList).intersection( [ 1, 2, 18, 102 ] ): xLog = True
        else: xLog = False
    if yLog == None:
        if set(mtList).intersection( [ 1, 2, 18, 102 ] ): yLog = True
        else: yLog = False

    # Set up the axis labels
    if uncratio: yAxisLabel = "$\Delta\sigma(E)/\sigma(E)$"
    else:        yAxisLabel = "Cross Section"
    xAxisLabel = '$E$'

    xAxisSettings = plot2d.AxisSettings( label=xAxisLabel, isLog=xLog, unit=xUnit )
    yAxisSettings = plot2d.AxisSettings( label=yAxisLabel, isLog=yLog, unit=yUnit )

    # Set plot limits
    try:
        lows, highs = zip( *[ d.domain() for d in rawEndfData ] )
        xAxisSettings.axisMin = min(lows)
        xAxisSettings.axisMax = max(highs)
    except ValueError:
        xAxisSettings.axisMin = None
        xAxisSettings.axisMax = None

    # Get the C4 data for plotting
    c4Data = []
    if ( args.c4File != None ): c4Data = getC4DataSets( args.c4File, mt, 3, target, projectile, plotSyle=plotStyle )

    # Get the EXFOR data for plotting
    exforData = []
    if ( not nox4 ) and ( not reaction is None ) and ( not uncratio ): 
        sym, A, m = getSymAFromGNDName( target )
        exforData = getEXFORSets( sym, A, reaction = projectile+','+reaction, quantity = "SIG", nox4evals = nox4evals, nox4legend = nox4legend, plotSyle=plotStyle  )

    # Set up the XY and XYdY data
    xyDataSets=[]
    xydyDataSets=[]
    for k in xyData:
        xyDataSets.append( DataSet2d( data = xyData[k],\
                                      legend = getStyle( plotStyle, 'xyCurves', k, 'legend' ), \
                                      lineWidth=getStyle( plotStyle, 'xyCurves', k, 'lineWidth' ), \
                                      lineStyle=getStyle( plotStyle, 'xyCurves', k, 'lineStyle' ), \
                                      color=getStyle( plotStyle, 'xyCurves', k, 'lineColor' ),\
                                      symbol=getStyle( plotStyle, 'xyCurves', k, 'symbol' ),\
                                      dataType=getStyle( plotStyle, 'xyCurves', k, 'dataType' ),\
                                      xUnit=getStyle( plotStyle, 'xyCurves', k, 'xUnit' ),\
                                      yUnit=getStyle( plotStyle, 'xyCurves', k, 'yUnit' ) ) )

    for k in xydyData:
        d = [ [x[0],x[1]] for x in xydyData[k] ]
        u = [ [x[0],x[2]] for x in xydyData[k] ]
        xydyDataSets.append( DataSet2d( data = d, \
                                        uncertainty = u,\
                                        legend = getStyle( plotStyle, 'xydyCurves', k, 'legend' ), \
                                        lineWidth=getStyle( plotStyle, 'xydyCurves', k, 'lineWidth' ), \
                                        lineStyle=getStyle( plotStyle, 'xydyCurves', k, 'lineStyle' ), \
                                        color=getStyle( plotStyle, 'xydyCurves', k, 'lineColor' ), \
                                        symbol=getStyle( plotStyle, 'xydyCurves', k, 'symbol' ),\
                                        dataType=getStyle( plotStyle, 'xydyCurves', k, 'dataType' ),\
                                        xUnit=getStyle( plotStyle, 'xydyCurves', k, 'xUnit' ),\
                                        yUnit=getStyle( plotStyle, 'xydyCurves', k, 'yUnit' ),\
                                        errorbarColor=getStyle( plotStyle, 'xydyCurves', k, 'errorColor' ) ) )

    # Actually make the plot
    if endfData + exforData + c4Data + xyDataSets + xydyDataSets != []:
        if False:
            print 'num. endfData:', len(endfData)
            print 'num. exforData:', len(exforData)
            print 'num. c4Data:', len(c4Data)
            print 'num. xyDataSets:', len(xyDataSets)
            print 'num. xydyDataSets:', len(xydyDataSets)
        if reaction.lower() == 'inel': 
            if mt in [ 91, 799, 749, 699, 649, 849, 891 ]:reactionString = product+"[c]"
            else: reactionString = str(product)+"["+str( mt%50 )+"]"
        else: reactionString = reaction.lower()
        plot2d.makePlot2d( exforData + c4Data + endfData + xyDataSets + xydyDataSets, xAxisSettings = xAxisSettings, 
            yAxisSettings = yAxisSettings, title = target.capitalize() + '(' + projectile+','+reactionString +')', legendOn = True, outFile = outFile, legendXY=legendXY )


#---------------------------------------------------
# multiplicity plots
#---------------------------------------------------
def makeMultiplicityPlot( gndMap, mt, projectile, target, product, xLog=None, xUnit='MeV', yLog=None, yUnit='', yFrame='centerOfMass', zLog=None, zUnit='b', outFile=None, plotStyle={}, legendXY=(0.05, 0.95) ): 

    raise NotImplementedError()


#---------------------------------------------------
# energy distribution plots
#---------------------------------------------------
def makeEnDistSlicesPlot( gndMap, mt, projectile, target, product, xLog=None, xUnit='MeV', yLog=None, yUnit='1/MeV', yFrame='centerOfMass', zLog=None, zUnit='b', outFile=None, plotStyle={}, legendXY=(0.05, 0.95) ): 

    raise NotImplementedError()


def makeEnDistContourPlot( gndMap, mt, projectile, target, product, xLog=None, xUnit='MeV', yLog=None, yUnit='1/MeV', yFrame='centerOfMass', zLog=None, zUnit='b', outFile=None, plotStyle={}, legendXY=(0.05, 0.95) ): 

    raise NotImplementedError()


#---------------------------------------------------
# energy/momentum deposition plots
#---------------------------------------------------
def makeEnDepPlot( gndMap, mt, projectile, target, product, xLog=None, xUnit='MeV', yLog=None, yUnit='MeV', yFrame='centerOfMass', zLog=None, zUnit='b', outFile=None, plotStyle={}, legendXY=(0.05, 0.95) ): 
    from fudge.gnd.productData import energyDeposition

    # Declare the main reaction type
    try: reaction = getEXFORRxn( mt )
    except KeyError: 
        if mt in range( 850, 871, 1 ): KeyError( "Cannot do lumped covariance in MT = "+ str( mt ) )
        else: raise KeyError( "Unknown MT: "+str( mt ) )
    
    # Set default units as needed
    if xUnit == None: xUnit = 'eV'
    if yUnit == None: yUnit = 'MeV'
        
    mtList = [ mt ]
    endfData = []
    
    # Get the ENDF data for plotting
    for endf in gndMap:
        print smallBanner( "Preparing data for "+endf )
        reactionSuite, covarianceSuite, plotSettings = gndMap[endf][:2]
        
        # Set up MT list.  
        # This is the list of all MT that we will search for data to sum into the plot.  
        # Usually, it is one element long, but there are some special cases: 
        #     summed channels that were not specified in the ENDF file and lumped channels....
        if not hasMT( reactionSuite, mt ):
            if   mt == 3:          
                # (n,non-el), use for catch-all "C=55" gammas
                print( "WARNING: fudge-2.0 considers (n,non-el) to be redundant and will not process it, so it is being removed from your plot list" )
                del( mtList[ mtList.index(3) ] )
            elif mt == 5: pass  # (n,X), no sumrule specified.  Treat as plain old reaction. 
            elif mt == 1: raise KeyError( "Cannot compute average energies for (n,tot)" )
            elif mt == 4:          
                # (n,n'), sum: 51-91
                mtList += getEvaluationMTs( reactionSuite, mtFilter = range( 51, 92, 1 ) )
            elif mt == 103:        
                # (n,p), sum: 600-649
                mtList += getEvaluationMTs( reactionSuite, mtFilter = range( 600, 650, 1 ) )
            elif mt == 104:        
                # (n,d), sum: 650-699
                mtList += getEvaluationMTs( reactionSuite, mtFilter = range( 650, 700, 1 ) )
            elif mt == 105:        
                # (n,t), sum: 700-749
                mtList += getEvaluationMTs( reactionSuite, mtFilter = range( 700, 750, 1 ) )
            elif mt == 106:        
                # (n,3He), sum: 750-799
                mtList += getEvaluationMTs( reactionSuite, mtFilter = range( 750, 800, 1 ) )
            elif mt == 107:        
                # (n,a), sum: 800-849
                mtList += getEvaluationMTs( reactionSuite, mtFilter = range( 800, 850, 1 ) )
            elif mt == 16:         
                # (n,2n), sum 875-891
                mtList += getEvaluationMTs( reactionSuite, mtFilter = range( 875, 892, 1 ) )
            elif mt == 18:
                # (n,f), sum: 19-21,38
                mtList += getEvaluationMTs( reactionSuite, mtFilter = [ 19, 20, 21, 38 ] ) 
        elif mt in range( 851, 871, 1 ): 
            raise KeyError( "Cannot compute energy balance for MT's in range 851-870, these are for lumped covariance" )
        else: pass # plain old reaction, do nothing special
        mtList = uniquify( mtList )

        # We may need resonance reconstruction
        if set(mtList).intersection( [ 1, 2, 18, 102 ] ) and args.doResonanceReconstruction:  
            reactionSuite.reconstructResonances(styleName='reconstructed')
        
        # We need to compute energy depositions
        reactionSuite.calculateDepositionData({'verbosity':0})
        
        # Collect complete list of matching reactions
        reactionList = []
        for MT in mtList: reactionList += getReactions( reactionSuite, MT )
        
        print "Exit Channels:", ', '.join( map( str, reactionList ) )
        
        # Get the reaction corresponding to the requested mainMT
        if len( reactionList ) == 1: mainReaction = 0
        else:                        mainReaction = 'sum'
        
        # Now plot the main reaction
        if mainReaction != None:
            Q = reactionList[ mainReaction ].outputChannel.Q
            cs = reactionList[ mainReaction ].crossSection
            csTable = cs[cs.nativeData]
            #availableEnergy = [ [ csTable[0][0], Q + csTable[0][0] ], [ csTable[-1][0], Q + csTable[-1][0] ] ]
            for particle in reactionList[ mainReaction ].outputChannel.particles:
                try:
                    dat = particle.data[ energyDeposition.component.genre ]
                    endfData.append( DataSet2d( dat.forms[ dat.nativeData ],\
                        legend = str(particle), lineWidth = 3 ) )
                except KeyError: pass
                
        # Whether to do lin-lin or log-log based on the MT
        if xLog == None: 
            if set(mtList).intersection( [ 1, 2, 18, 102, 3, 105 ] ): xLog = True
            else: xLog = False
        if yLog == None: yLog = False
    
        # Set up the axis labels
        yAxisLabel = "$<E'>$"
        xAxisLabel = '$E$'
    
        xAxisSettings = plot2d.AxisSettings( label=xAxisLabel, isLog=xLog, unit=xUnit )
        yAxisSettings = plot2d.AxisSettings( label=yAxisLabel, isLog=yLog, unit=yUnit )

        # Set plot limits
        try:
            lows, highs = zip( *[ d.domain() for d in rawEndfData ] )
            xAxisSettings.axisMin = min(lows)
            xAxisSettings.axisMax = max(highs)
        except ValueError:
            xAxisSettings.axisMin = None
            xAxisSettings.axisMax = None

        # Actually make the plot
        if endfData != []:
            if reaction.lower() == 'inel': 
                if product == None:reactionString = reaction.lower()
                else:
                    if mt in [ 91, 799, 749, 699, 649, 849, 891 ]: reactionString = product+"[c]"
                    elif mt in [ 4, 103, 104, 105, 106, 107 ]: reactionString = reaction.lower()
                    else: reactionString = product+"["+str( mt%50 )+"]"
            else: reactionString = reaction.lower()
            plot2d.makePlot2d( endfData, xAxisSettings = xAxisSettings, yAxisSettings = yAxisSettings, \
                title = target.capitalize() + '(' + projectile+','+reactionString +')', \
                legendOn = True, outFile = outFile )


def makeMomDepPlot( gndMap, mt, projectile, target, product, xLog=None, xUnit='MeV', yLog=None, yUnit='MeV/c', yFrame='centerOfMass', zLog=None, zUnit='b', outFile=None, plotStyle={}, legendXY=(0.05, 0.95) ): 

    raise NotImplementedError()


#---------------------------------------------------
# angular distribution plots
#---------------------------------------------------
def makeAngDistContourPlot( gndMap, mt, projectile, target, product, xLog=None, xUnit='eV', yLog=None, yUnit='', yFrame='centerOfMass', zLog=None, zUnit='b', outFile=None, plotStyle={}, legendXY=(0.05, 0.95) ): 

    # Declare the main reaction type
    try: reaction = getEXFORRxn( mt )
    except KeyError: 
        if mt in range( 850, 871, 1 ): 
            print( "Skipping MT = "+ str( mt )+", it is lumped covariance" )
            return
        else: raise KeyError( "Unknown MT: "+str( mt ) )
    
    # Set default units as needed
    if xUnit == None: xUnit = "eV"
    if yUnit == None: yUnit = ''
    if zUnit == None: zUnit = 'b'

    # Get the EXFOR data for plotting
    exforData = []
    if True:
        sym, A, m = getSymAFromGNDName( target )
        exforData = getEXFORSets( sym, A, reaction = projectile+','+reaction, quantity = "DA", nox4evals = True, nox4legend = False )

    # Get the C4 data for plotting
    c4Data = []
    if False: c4Data = getC4DataSets( args.c4File, mt, 4, target, projectile )

    # Get the ENDF data for plotting
    endfData = []
    for endf in gndMap:
        print smallBanner( "Preparing data for "+endf )
        reactionSuite, covarianceSuite, plotSettings = gndMap[endf]
        reactionList = getReactions( reactionSuite, mt )
        for r in reactionList:
            print "Retrieving data for",r
            csData = getPointwiseCrossSection( r )
            angDistData = getAngularDistribution( r, product )
            for a in angDistData:
                endfData.append( a.toPointwise_withLinearXYs( 1e-6 ) )
                
    # Set up the axes
    xAxisSettings = plot2d.AxisSettings( label="$E'$", isLog=True, unit=xUnit )
    yAxisSettings = plot2d.AxisSettings( label="$\\mu = \\cos( \\theta )$", isLog=yLog, unit=yUnit )
    zAxisSettings = plot2d.AxisSettings( label="$d\\sigma(E)/d\\mu $", isLog=yLog, unit=yUnit )

    # Actually make the plot
    if endfData != []:
        if reaction.lower() == 'inel': reactionString = product+"["+str( mt%50 )+"]"
        else: reactionString = reaction.lower()
        plot2d.makePlot2dContour( endfData, xAxisSettings = xAxisSettings, yAxisSettings = yAxisSettings, 
            title = "$d^2\\sigma(E)/dE'd{\\mu}$ for emitted "+product+"'s from "+target.capitalize() + '(' + projectile+', '+reactionString +')', legendOn = False, outFile = outFile, legendXY=legendXY )
            
            
def makeAngDistSlicesPlot( gndMap, mt, projectile, target, product, xLog=None, xUnit='MeV', yLog=None, yUnit='', yFrame='centerOfMass', zLog=None, zUnit='b', outFile=None, plotStyle={}, legendXY=(0.05, 0.95) ): 

    raise NotImplementedError()


def makeAngDistMubarPlot( gndMap, mt, projectile, target, product, xLog=None, xUnit='eV', yLog=None, yUnit='', yFrame='centerOfMass', outFile=None, plotStyle={}, legendXY=(0.05,0.95) ): 

    # Declare the main reaction type
    try: reaction = getEXFORRxn( mt )
    except KeyError: 
        if mt in range( 850, 871, 1 ): 
            print( "Skipping MT = "+ str( mt )+", it is lumped covariance" )
            return
        else: raise KeyError( "Unknown MT: "+str( mt ) )
    
    # Set default units as needed
    if xUnit == None: xUnit = 'eV'
    if yUnit == None: yUnit = ""

    # Get the ENDF data for plotting
    endfData = []
    for endf in gndMap:
        print smallBanner( "Preparing data for "+endf )
        reactionSuite, covarianceSuite, plotSettings = gndMap[endf]
        reactionList = getReactions( reactionSuite, mt )
        for r in reactionList:
            print "Retrieving data for",r
            csData = getPointwiseCrossSection( r )
            angDistData = getAngularDistribution( r, product )
            for a in angDistData:
                xUnit = a.axes[0].unit
                if yFrame != a.productFrame: raise ValueError( "Found frame ",a.productFrame, "you requested frame", yFrame )
                endfData.append( DataSet2d( [ [ E, a.muAverageAtEnergy(E) ] for E in a.getEnergyArray() ], xUnit='eV', yUnit='', lineWidth=3, legend=product+' given as '+a.name+' ('+endf+')' ) )
        print
        
    # Set up the axes
    if yFrame == 'centerOfMass': yFrameStr = 'cm'
    else: yFrameStr = 'lab'
    xAxisSettings = plot2d.AxisSettings( label="$E$", isLog=xLog, unit=xUnit )
    yAxisSettings = plot2d.AxisSettings( label="$<\mu_{"+yFrameStr+"}>$", isLog=yLog, unit=yUnit )

    # Actually make the plot
    if endfData != []:
        if reaction.lower() == 'inel': reactionString = product+"["+str( mt%50 )+"]"
        else: reactionString = reaction.lower()
        plot2d.makePlot2d( endfData, xAxisSettings = xAxisSettings, yAxisSettings = yAxisSettings, 
            title = "$<\mu_{"+yFrameStr+"}>$ for emitted "+product+"'s from "+target.capitalize() + '(' + projectile+', '+reactionString +')', legendOn = True, legendXY=legendXY, outFile = outFile )
        
            
#---------------------------------------------------
# main routine!
#---------------------------------------------------
if __name__ == "__main__": 

    args = parseArgs()

    # Declare the various particles involved in the plot
    target = args.target
    projectile = args.projectile
    product = args.product

    # Declare the lists of data from various sources
    rawEndfData = []
    rawEndfUnc = []
    endfData = []
    
    # Get the target & projectile from the command line.
    # If it's not set, we'll have to extract it from one of the evaluation files.
    targ = args.target 
    proj = args.projectile
    
    # Read in the plot style information
    plotStyle = json.loads(open(os.path.abspath(__file__).replace('plotEvals.py','plot_defaults.json')).read())
    if args.style != None:
        if not os.path.exists( args.style ): raise IOError("Cannot file plot style file %s" % args.style)
        overrides = json.loads( open( args.style ).read() )
        for k in overrides:
            plotStyle[k].update(overrides[k])
    if False:  # for debugging
        print json.dumps(plotStyle,indent=2)
        exit()
    
    # Read the ENDF (& other user provided) data
    xyData = {}
    xydyData = {}
    gndMap = collections.OrderedDict()
    mtMap = {}
    print bigBanner( "Reading evaluation files" )
    for endf in args.endf:
        if endf == "None": continue
        print smallBanner( "Reading "+endf )

        # Is the file a GND file?
        if open(endf).readline().startswith( "<?xml" ):
            RS = fudge.gnd.reactionSuite.readXML( endf )
            try: 
                CS = fudge.gnd.covariances.covarianceSuite.readXML( endf.replace( '.gnd.', '.gndCov.' ) )
            except: 
                CS = fudge.gnd.covariances.covarianceSuite()
            gndMap[endf] = [ RS, CS ]

        # Maybe its an ENDF file?
        elif open(endf).readline().endswith(' 0  0    0\n') or endf.endswith('.endf'):
            results = endfFileToGND.endfFileToGND( endf, toStdOut = False, skipBadData = True )
            if type( results ) == dict:
                gndMap[endf] = [ results['reactionSuite'], results['covarianceSuite'] ]
            elif type( results ) == tuple:
                gndMap[endf] = [ results[0], results[1] ]
            else:
                raise TypeError( "endfFileToGND.endfFileToGND() returned a "+str(type(results))+", I don't know what to do with it" )

        # Try plain XY data
        elif endf.endswith( '.xy.dat'):
            xyData[endf]=readXYData( endf )
        
        # Try plain XYdY data
        elif endf.endswith( '.xydy.dat'):
            xydyData[endf]=readXYdYData( endf )

        # OK, try AMPX
        elif fnmatch.fnmatch(endf,'*.ampx*'): 
            try:
                import ampx2fudge, ampx
                ampx_za, bounds = ampx.readEvaluation( endf, str(targ), str(proj) )
                gndMap[endf] = [ 
                    ampx2fudge.convertAmpxNuclideToFudgeReactionSuite(ampx_za, bounds), 
                    fudge.gnd.covariances.covarianceSuite() ]
            except ImportError:
                print "WARNING: Could not load AMPX module.  Is it in your path?"
        
        # Failed!
        else: print "WARNING: Unknown file type, not reading %s"% endf
                                
        if not endf in gndMap: continue
        
        if targ == None:    targ = gndMap[endf][0].target
        if proj == None:    proj = gndMap[endf][0].projectile
        
        mtMap[endf] = getEvaluationMTs( gndMap[endf][0] )


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
    print bigBanner( "Plot details" )
    if target == None:     target = gndMap.items()[0][1][0].target.name
    if projectile == None: projectile = gndMap.items()[0][1][0].projectile.name
    if product == None and args.mf not in [3,23]: 
        raise ValueError( "For MF="+str(args.mf)+", need to declare a product using --product command line argument" )
    print "Projectile is "+projectile
    print "Target is "+target
    print "List of MT's to plot:",mtList

    
    # Determine the observable type
    if args.Ebar: print "Observable is average energy deposition"
    elif args.Pzbar: print "Observable is average forward momentum"
    elif args.mf in [ 3,23 ]: print "Observable is "+observableNameMap[args.mf]
    else: print "Observable is "+observableNameMap[args.mf]+' for '+product
    print

    # Loop over reactions
    for mt in mtList:    
    
        # What is the proper name of the reaction?
        reaction = projectile.capitalize() +','+ getEXFORRxn( mt )
        
        # Compute the output filename
        if args.outFile != None:
            if not args.outFile.endswith('.png'):
                outFile = args.outFile+'_mt'+str(mt)+'_mf'+str(args.mf)+'.png'
            else: outFile = args.outFile
        else: outFile = None
        
        # Determine the kind of plot and set up the axes
        print bigBanner( "Generating plots for "+target+"("+reaction.lower()+')' )
        if args.Ebar:
            makeEnDepPlot( gndMap, mt=mt, projectile=projectile, target=target, xLog=args.logx, xUnit=args.xUnit, yLog=args.logy, yUnit=args.yUnit, outFile=outFile, plotStyle=plotStyle )
        elif args.Pzbar:pass
        else:
            if args.mf in [3,23]: 
                makeCSPlot( \
                    gndMap, mt=mt, xyData=xyData, xydyData=xydyData,\
                    projectile=projectile, target=target, \
                    uncratio=args.uncratio, nounc=args.nounc,\
                    nox4=args.nox4, showparts=args.showparts, nox4evals=args.nox4evals, \
                    nox4legend=args.nox4legend, \
                    xLog=args.logx, xUnit=args.xUnit, yLog=args.logy, yUnit=args.yUnit, \
                    outFile=outFile, plotStyle=plotStyle, \
                    legendXY=( args.legendx, args.legendy ),\
                    doResonanceReconstruction=args.doResonanceReconstruction, heatToTemp=None )
            elif args.mf == 4: 
                if args.angContour: 
                    makeAngDistContourPlot( \
                        gndMap, mt=mt, \
                        projectile=projectile, target=target, product=product, \
                        xLog=args.logx, xUnit=args.xUnit, yLog=args.logy, yUnit=args.yUnit, zLog=args.logz, zUnit=args.zUnit, \
                        outFile=outFile, plotStyle=plotStyle, \
                        legendXY=( args.legendx, args.legendy )  )
                elif args.angMubar: 
                    makeAngDistMubarPlot( \
                        gndMap, mt=mt, \
                        projectile=projectile, target=target, product=product, \
                        xLog=args.logx, xUnit=args.xUnit, yLog=args.logy, yUnit=args.yUnit, \
                        outFile=outFile, plotStyle=plotStyle, \
                        legendXY=( args.legendx, args.legendy ) )
                elif args.Pzbar:    
                    makeMomDepPlot( \
                        gndMap, mt=mt, \
                        projectile=projectile, target=target, product=product, \
                        xLog=args.logx, xUnit=args.xUnit, yLog=args.logy, yUnit=args.yUnit, \
                        outFile=outFile, plotStyle=plotStyle, \
                        legendXY=( args.legendx, args.legendy ) )
                else:               
                    makeAngDistSlicesPlot( \
                        gndMap, mt=mt, \
                        projectile=projectile, target=target, product=product, \
                        xLog=args.logx, xUnit=args.xUnit, yLog=args.logy, yUnit=args.yUnit, zLog=args.logz, zUnit=args.zUnit, \
                        outFile=args.outFile, plotStyle=plotStyle, \
                        legendXY=( args.legendx, args.legendy ) )
            elif args.mf == 5: 
                if args.enContour:  
                    makeEnDistContourPlot( \
                        gndMap, mt=mt, \
                        projectile=projectile, target=target, product=product, \
                        xLog=args.logx, xUnit=args.xUnit, yLog=args.logy, yUnit=args.yUnit, zLog=args.logz, zUnit=args.zUnit, \
                        outFile=outFile, plotStyle=plotStyle, \
                        legendXY=( args.legendx, args.legendy ) )
                elif args.Ebar:     
                    makeEnDepPlot( \
                        gndMap, mt=mt, \
                        projectile=projectile, target=target, product=product, \
                        xLog=args.logx, xUnit=args.xUnit, yLog=args.logy, yUnit=args.yUnit, \
                        outFile=outFile, plotStyle=plotStyle, \
                        legendXY=( args.legendx, args.legendy ) )
                else:               
                    makeEnDistSlicesPlot( \
                        gndMap, mt=mt, \
                        projectile=projectile, target=target, product=product, \
                        xLog=args.logx, xUnit=args.xUnit, yLog=args.logy, yUnit=args.yUnit, zLog=args.logz, zUnit=args.zUnit, \
                        outFile=outFile, plotStyle=plotStyle, \
                        legendXY=( args.legendx, args.legendy ) )
            else:                   
                raise NotImplementedError( "Only MF=3 (cross section data), MF=4 (angular distributions) & MF=5 (energy distributions)  currently implemented, you requested MF="+str(args.mf) )

    
