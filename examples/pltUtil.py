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
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# <<END-copyright>>

import os, sys, math
sys.path.append( os.sep.join( [ os.environ[ 'HOME' ], 'apps', 'fudge2' ] ) )
#sys.path.append( os.sep.join( [ os.environ[ 'HOME' ], 'Projects', 'Current', 'x4i-1.0' ] ) )
sys.path.append( os.sep.join( [ os.environ[ 'HOME' ], 'apps', 'x4i-1.0' ] ) )
from fudge.vis.matplotlib import plot2d, DataSet2d, DataSet3d, defaultSymbols, nSymbols, defaultColors, nColors
from fudge.legacy.converting import endfFileToGND
from fudge.core.math.xData import XYs, axes
from pqu.physicalQuantityWithUncertainty import PhysicalQuantityWithUncertainty as PQU
from fudge.core.utilities import argparse, fudgeZA
from fudge.core.utilities.brb import uniquify
from fudge.gnd.covariances import tokens as covTokens 



#---------------------------------------------------
# Set up the command line parser
#---------------------------------------------------
parser = argparse.ArgumentParser(description='Plot cross section data')
parser.add_argument('mt',       metavar='mt', type=int, nargs=1, help='MT of the cross section to plot' )
parser.add_argument('endf',     metavar='endf', type=str, nargs='+', help='ENDF file(s) whose cross section you want to plot.  Use "None" for no input file.' )
parser.add_argument('--isGND',  default=False, action='store_true', help='Input file is not an ENDF formatted file, it is a GND formatted file' )
parser.add_argument('-o',       dest='outFile', default=None, type=str, help='Output file for plot (disables interactive plotting)' )
parser.add_argument('--mf',     type=int, default=3, help='MF of original data to plot [default is 3, cross section data, use 4 for angular, 5 for energy and 6 for energy-angle distributions and 12 for multiplicity]' )

parser.add_argument('--c4File',         default=None,  type=str, help="Optional C4 file to pull data from instead of (or in addition to) EXFOR" )
parser.add_argument('--nox4',           default=False, action='store_true', help='Do not plot EXFOR data' )
parser.add_argument('--nox4evals',      default=False, action='store_true', help='Do not plot evaluations foudn in EXFOR library (with "V" SUBENT)' )
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

parser.add_argument('--logx', default=None, action='store_true', help='Make x scale logrithmic (Default is False except for MT in [1,2,18,102])' )
parser.add_argument('--logy', default=None, action='store_true', help='Make y scale logrithmic (Default is False except for MT in [1,2,18,102])' )
parser.add_argument('--logz', default=None, action='store_true', help='Make z scale logrithmic' )

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
args = parser.parse_args()

#---------------------------------------------------
# utilities
#---------------------------------------------------
observableNameMap = { 3:"Cross Section ( Sigma(E) )", 4:"Angular Distribution ( dSigma(E)/d(cos Theta) )", 5:"Energy Distribution ( dSigma(E)/d(E') )", 6:"Energy-Angle Distribution ( dSigma(E)/d(E')d(cos Theta) )" }


# Simplified MT - EXFOR reaction mapping
def getEXFORRxn( MT ):
    if MT in range( 50, 92 ): return 'INEL'
    if MT in range( 600, 650 ): return 'P'
    if MT in range( 650, 700 ): return 'D'
    if MT in range( 700, 750 ): return 'T'
    if MT in range( 750, 800 ): return 'HE3'
    if MT in range( 800, 850 ): return 'A'
    return { 1:'TOT', 2:'EL', 3:'NON', 4:'INEL', 5:'X', 11:'2N+D', 16:'2N', 17:'3N', 18:'F', 22:'N+A', 102:'G', 103:'P', 104:'D', 105:'T', 107:'A', 111:'2P', 28:'N+P', 32:'N+D', 33:'N+T', 37:'4N', 41:'2N+P', 42:'3N+P' }[ MT ]


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
    Z = fudgeZA.SymbolToZ( sym )
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
    for form in ( covTokens.covarianceFormToken, covTokens.mixedFormToken ):
        if form in theCovariance.forms:
            theUncert = theCovariance.forms[ form ].getUncertaintyVector( theData, relative=False )
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
    return reactionSuite.reactions + reactionSuite.summedReactions + reactionSuite.fissionComponents # + reactionSuite.productions


def getEvaluationMTs( reactionSuite, mtFilter = None ):
    reactionList = getPlottableReactions( reactionSuite )
    if mtFilter == None: return [ int( r.attributes['ENDF_MT' ] ) for r in reactionList ]
    else:                return [ int( r.attributes['ENDF_MT' ] ) for r in reactionList if int( r.attributes['ENDF_MT' ] ) in mtFilter ]


def getReactions( reactionSuite, MT ): 
    return [ r for r in getPlottableReactions( reactionSuite ) if r.attributes['ENDF_MT' ] == str(MT) ]


def getPointwiseCrossSection( reac ):
    if   'linear'            in reac.crossSection.forms: return reac.crossSection['linear']
    elif 'pointwise'         in reac.crossSection.forms: return reac.crossSection['pointwise']
    elif 'piecewise'         in reac.crossSection.forms: return reac.crossSection['piecewise'].toPointwise_withLinearXYs( 0, 1e-8 )
    elif 'weightedPointwise' in reac.crossSection.forms: return reac.crossSection['weightedPointwise'].toPointwise_withLinearXYs( 0, 1e-8 )
    elif 'resonancesWithBackground' in reac.crossSection.forms: # in case --noReconstruct option is used:
        return reac.crossSection['resonancesWithBackground'].tabulatedData.toPointwise_withLinearXYs( 0, 1e-8 ) 
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
    
    
def getEXFORSets( sym, A, reaction = None, quantity = "SIG", nox4evals = True, nox4legend = False, forceLegend = False ): 
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
    print sym+'-'+A, reaction, quantity
    if len( subents.keys() ) > 15: print bigBanner( 'Retrieving %i entries' % ( len( subents.keys() ) ) )
    else: print smallBanner( 'Retrieving entries:'+str( subents.keys() ) )
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
    

def getC4DataSets( c4File, mt, mf, target, projectile ):
    try:
        from empire import c4
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



#---------------------------------------------------
# cross section plots
#---------------------------------------------------
def makeCSPlot( gndMap, mt, projectile, target, uncratio=False, nounc=False, nox4=False, showparts=False, nox4evals=True, nox4legend=False, xLog=None, xUnit='eV', yLog=None, yUnit='b', outFile=args.outFile, legendXY=(0.05, 0.95)):

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
        
        # Set up MT list.  
        # This is the list of all MT that we will search for data to sum into the plot.  
        # Usually, it is one element long, but there are some special cases: summed channels that were not specified in the ENDF file and lumped channels....
        mtList = [ mt ]
        if not hasMT( reactionSuite, mt ):
            if   mt == 3:          
                # (n,non-el), use for catch-all "C=55" gammas
                print( "WARNING: fudge-2.0 considers (n,non-el) to be redundant and will not process it, so it is being removed from your plot list" )
                del( mtList[ mtList.index(3) ] )
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
            # lumped covariance MT's are 851-870
            # must extract the MT's from the covariance file
            for lump in covarianceSuite.lumpedChannels:
                if lump.ENDF_MFMT == '33,'+str(mt): break
            mtList += getEvaluationMTs( reactionSuite, mtFilter = [ int( channel.attributes['ENDF_MFMT'].split(',')[1] ) for channel in lump.channels] )
            reaction = 'N,' + getEXFORRxn( mtList[-1] )
        else: pass # plain old reaction, do nothing special
        mtList = uniquify( mtList )

        if set(mtList).intersection( [ 1, 2, 18, 102 ] ) and args.doResonanceReconstruction:  # need resonance reconstruction
            reactionSuite.reconstructResonances()
        
        # collect complete list of matching reactions
        reactionList = []
        for MT in mtList: reactionList += getReactions( reactionSuite, MT )
        
        print "Exit Channels:", ', '.join( map( str, reactionList ) )
        
        # get the reaction corresponding to the requested mainMT
        if len( reactionList ) == 1: mainReaction = 0
        else:                        mainReaction = 'sum'
        
        # now plot the main reaction, this one gets the covariance
        if mainReaction != None:
            csData = getPointwiseCrossSection( reactionList[0] )
            if ( not nounc ) and mt not in [ 1, 2, 18, 102 ] :
                csData = csData.thicken( sectionSubdivideMax = 20 ) # only thicken cross sections that don't have enough points in them to resolve uncertainty steps
            if mt != 18: # element 0 already was sum
                for endfRxn in reactionList[1:]:
                    csData += getPointwiseCrossSection( endfRxn ) # this will need try/except block sooner or later
            rawEndfData.append( csData )
            if args.nounc: rawEndfUnc.append( None )
            else:
                cov = [ c for c in covarianceSuite.sections if hasattr(c,'rowData') and c.rowData.attributes['ENDF_MFMT'] == '33,%i' % mt ]
                if not cov: 
                    print ( "MT = %i (%s) covariance not present in the file" % ( mt, reaction ) )
                    rawEndfUnc.append( None )
                else:
                    rawEndfUnc.append( getUncertainty( cov[0], rawEndfData[-1] ) )
            if uncratio:
                if not cov: raise TypeError( "File has no covariance" ) 
                ratioDomain = ( max(rawEndfUnc[-1].domainMin(), rawEndfData[-1].domainMin()),
                        min(rawEndfUnc[-1].domainMax(), rawEndfData[-1].domainMax()) )
                ratio = rawEndfUnc[-1].xSlice(*ratioDomain) / rawEndfData[-1].xSlice(*ratioDomain)
                endfData.append( DataSet2d( ratio, uncertainty = None, legend = endf + ': (' + projectile+','+reaction.lower() + ')', lineWidth = 3 ) )
            else:
                endfData.append( DataSet2d( rawEndfData[-1], uncertainty = rawEndfUnc[-1],
                    legend = endf + ': (' + projectile+','+reaction.lower() + ')', lineWidth = 3 ) )
        
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
        lows, highs = zip( *[ d.getDomain() for d in rawEndfData ] )
        xAxisSettings.axisMin = min(lows)
        xAxisSettings.axisMax = max(highs)
    except ValueError:
        xAxisSettings.axisMin = None
        xAxisSettings.axisMax = None

    # Get the C4 data for plotting
    c4Data = []
    if ( args.c4File != None ): c4Data = getC4DataSets( args.c4File, mt, 3, target, projectile )

    # Get the EXFOR data for plotting
    exforData = []
    if ( not nox4 ) and ( not reaction is None ) and ( not uncratio ): 
        sym, A, m = getSymAFromGNDName( target )
        exforData = getEXFORSets( sym, A, reaction = projectile+','+reaction, quantity = "SIG", nox4evals = nox4evals, nox4legend = nox4legend )

    # Actually make the plot
    if endfData + exforData +c4Data != []:
        if reaction.lower() == 'inel': 
            if mt in [ 91, 799, 749, 699, 649, 849, 891 ]:reactionString = product+"[c]"
            else: reactionString = product+"["+str( mt%50 )+"]"
        else: reactionString = reaction.lower()
        plot2d.makePlot2d( exforData + c4Data + endfData, xAxisSettings = xAxisSettings, 
            yAxisSettings = yAxisSettings, title = target.capitalize() + '(' + projectile+','+reactionString +')', legendOn = True, outFile = outFile, legendXY=legendXY )


#---------------------------------------------------
# multiplicity plots
#---------------------------------------------------
def makeMultiplicityPlot( gndMap, mt, projectile, target, product, xLog=None, xUnit=args.xUnit, yLog=None, yUnit=args.yUnit, yFrame='centerOfMass', zLog=None, zUnit='b', outFile=args.outFile, legendXY=(0.05, 0.95) ): 

    raise NotImplementedError()


#---------------------------------------------------
# energy distribution plots
#---------------------------------------------------
def makeEnDistSlicesPlot( gndMap, mt, projectile, target, product, xLog=None, xUnit=args.xUnit, yLog=None, yUnit=args.yUnit, yFrame='centerOfMass', zLog=None, zUnit='b', outFile=args.outFile, legendXY=(0.05, 0.95) ): 

    raise NotImplementedError()


def makeEnDistContourPlot( gndMap, mt, projectile, target, product, xLog=None, xUnit=args.xUnit, yLog=None, yUnit=args.yUnit, yFrame='centerOfMass', zLog=None, zUnit='b', outFile=args.outFile, legendXY=(0.05, 0.95) ): 

    raise NotImplementedError()


#---------------------------------------------------
# energy/momentum deposition plots
#---------------------------------------------------
def makeEnDepPlot( gndMap, mt, projectile, target, product, xLog=None, xUnit=args.xUnit, yLog=None, yUnit=args.yUnit, yFrame='centerOfMass', zLog=None, zUnit='b', outFile=args.outFile, legendXY=(0.05, 0.95) ): 

    raise NotImplementedError()


def makeMomDepPlot( gndMap, mt, projectile, target, product, xLog=None, xUnit=args.xUnit, yLog=None, yUnit=args.yUnit, yFrame='centerOfMass', zLog=None, zUnit='b', outFile=args.outFile, legendXY=(0.05, 0.95) ): 

    raise NotImplementedError()


#---------------------------------------------------
# angular distribution plots
#---------------------------------------------------
def makeAngDistContourPlot( gndMap, mt, projectile, target, product, xLog=None, xUnit=args.xUnit, yLog=None, yUnit=args.yUnit, yFrame='centerOfMass', zLog=None, zUnit='b', outFile=args.outFile, legendXY=(0.05, 0.95) ): 

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
        reactionSuite, covarianceSuite = gndMap[endf]
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
            
            
def makeAngDistSlicesPlot( gndMap, mt, projectile, target, product, xLog=None, xUnit=args.xUnit, yLog=None, yUnit=args.yUnit, yFrame='centerOfMass', zLog=None, zUnit='b', outFile=args.outFile, legendXY=(0.05, 0.95) ): 

    raise NotImplementedError()


def makeAngDistMubarPlot( gndMap, mt, projectile, target, product, xLog=None, xUnit=args.xUnit, yLog=None, yUnit=args.yUnit, yFrame='centerOfMass', outFile=args.outFile, legendXY=(0.50,0.2) ): 

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
        reactionSuite, covarianceSuite = gndMap[endf]
        reactionList = getReactions( reactionSuite, mt )
        for r in reactionList:
            print "Retrieving data for",r
            csData = getPointwiseCrossSection( r )
            angDistData = getAngularDistribution( r, product )
            for a in angDistData:
                xUnit = a.axes[0].unit
                if yFrame != a.productFrame: raise ValueError( "Found frame ",a.productFrame, "you requested frame", yFrame )
                endfData.append( DataSet2d( [ [ E, a.muAverageAtEnergy(E) ] for E in a.getEnergyArray() ], xUnit='eV', yUnit='', lineWidth=3, legend=product+' given as '+a.getName()+' ('+endf+')' ) )
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
    
    # Declare the various particles involved in the plot
    target = args.target
    projectile = args.projectile
    product = args.product
    
    # Declare the lists of data from various sources
    rawEndfData = []
    rawEndfUnc = []
    endfData = []

    # Check the reaction name
    if args.mt == []: exit() # Nothing to plot!
    try: reaction = getEXFORRxn( args.mt[0] )
    except KeyError: 
        if args.mt[0] in range( 850, 871, 1 ): print( "Got lumped covariance for MT = "+ str( args.mt[0] ) )
        else: raise KeyError( "Unknown MT: "+str( args.mt[0] ) )
    
    # Read the ENDF data
    gndMap = {}
    print bigBanner( "Reading ENDF files" )
    for endf in args.endf:
        if endf == "None": continue
        print smallBanner( "Reading "+endf )
        try: rce = endfFileToGND.endfFileToGND( endf, toStdOut = False, skipBadData = True )
        except:
            from fudge.gnd import reactionSuite, covariances
            rce = { 'reactionSuite': reactionSuite.readXML( endf), 'covarianceSuite':covariances.covarianceSuite() }
            endfCovars = endf.replace('.xml','-cov.xml')
            if os.path.exists(endfCovars) and endfCovars!=endf:
                rce['covarianceSuite'] = covariances.readXML( endfCovars, reactionSuite=rce['reactionSuite'] )
        gndMap[endf] = ( rce['reactionSuite'], rce['covarianceSuite'] )
        print
          
    # Determine the names of the particles involved
    print bigBanner( "Plot details" )
    if target == None:     target = gndMap[endf][0].target.name
    if projectile == None: projectile = gndMap[endf][0].projectile.name
    reaction = projectile.capitalize() +','+ reaction
    if product == None and args.mf != 3: raise ValueError( "For MF="+str(args.mf)+", need to declare a product using --product command line argument" )
    print "Projectile is "+projectile
    print "Target is "+target
    print "Reaction is "+target+"("+reaction.lower()+')'
    if args.mf in [ 3 ]:
        print "Observable is "+observableNameMap[args.mf]
    else:
        print "Observable is "+observableNameMap[args.mf]+' for '+product
    print
    
    # Determine the kind of plot and set up the axes
    print bigBanner( "Generating plot" )
    if args.mf == 3:        
        makeCSPlot( \
            gndMap, mt=args.mt[0], \
            projectile=projectile, target=target, \
            uncratio=args.uncratio, nounc=args.nounc, \
            nox4=args.nox4, showparts=args.showparts, nox4evals=args.nox4evals, \
            nox4legend=args.nox4legend, \
            xLog=args.logx, xUnit=args.xUnit, yLog=args.logy, yUnit=args.yUnit, \
            outFile=args.outFile, \
            legendXY=( args.legendx, args.legendy ) )
    elif args.mf == 4: 
        if args.angContour: 
            makeAngDistContourPlot( \
                gndMap, mt=args.mt[0], \
                projectile=projectile, target=target, product=product, \
                xLog=args.logx, xUnit=args.xUnit, yLog=args.logy, yUnit=args.yUnit, zLog=args.logz, zUnit=args.zUnit, \
                outFile=args.outFile, \
                legendXY=( args.legendx, args.legendy )  )
        elif args.angMubar: 
            makeAngDistMubarPlot( \
                gndMap, mt=args.mt[0], \
                projectile=projectile, target=target, product=product, \
                xLog=args.logx, xUnit=args.xUnit, yLog=args.logy, yUnit=args.yUnit, \
                outFile=args.outFile, \
                legendXY=( args.legendx, args.legendy ) )
        elif args.Pzbar:    
            makeMomDepPlot( \
                gndMap, mt=args.mt[0], \
                projectile=projectile, target=target, product=product, \
                xLog=args.logx, xUnit=args.xUnit, yLog=args.logy, yUnit=args.yUnit, \
                outFile=args.outFile, \
                legendXY=( args.legendx, args.legendy ) )
        else:               
            makeAngDistSlicesPlot( \
                gndMap, mt=args.mt[0], \
                projectile=projectile, target=target, product=product, \
                xLog=args.logx, xUnit=args.xUnit, yLog=args.logy, yUnit=args.yUnit, zLog=args.logz, zUnit=args.zUnit, \
                outFile=args.outFile, \
                legendXY=( args.legendx, args.legendy ) )
    elif args.mf == 5: 
        if args.enContour:  
            makeEnDistContourPlot( \
                gndMap, mt=args.mt[0], \
                projectile=projectile, target=target, product=product, \
                xLog=args.logx, xUnit=args.xUnit, yLog=args.logy, yUnit=args.yUnit, zLog=args.logz, zUnit=args.zUnit, \
                outFile=args.outFile, \
                legendXY=( args.legendx, args.legendy ) )
        elif args.Ebar:     
            makeEnDepPlot( \
                gndMap, mt=args.mt[0], \
                projectile=projectile, target=target, product=product, \
                xLog=args.logx, xUnit=args.xUnit, yLog=args.logy, yUnit=args.yUnit, \
                outFile=args.outFile, \
                legendXY=( args.legendx, args.legendy ) )
        else:               
            makeEnDistSlicesPlot( \
                gndMap, mt=args.mt[0], \
                projectile=projectile, target=target, product=product, \
                xLog=args.logx, xUnit=args.xUnit, yLog=args.logy, yUnit=args.yUnit, zLog=args.logz, zUnit=args.zUnit, \
                outFile=args.outFile, \
                legendXY=( args.legendx, args.legendy ) )
    else:                   
        raise NotImplementedError( "Only MF=3 (cross section data), MF=4 (angular distributions) & MF=5 (energy distributions)  currently implemented, you requested MF="+str(args.mf) )

    
