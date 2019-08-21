#! /usr/bin/env python

# <<BEGIN-copyright>>
# <<END-copyright>>

#
# This routine calls resolutionTestFile.py.
#
# USAGE:
#
#   resolutionTestPlotDataset.py ENDL_file temperature dataset_number_within_ENDL_file
#
# OPTIONS:
#
#   ENDL_file                               ENDL cross section file to heat.
#   temperature                             Temperature to heat to in MeV.
#   dataset_number_within_ENDL_file         some ENDL files have more than 1 cross section (e.g., S = 1), this number is the index 
#                                               (0 based) of the cross section dataset to plot.
#
import sys, os, glob
from fudge.core.math.xData import axes, XYs
from fudge.vis.gnuplot import fudgeMultiPlots

if( len( sys.argv ) != 4 ) : raise Exception( 'need ENDL_file, temperature and dataset_number_within_ENDL_file' )

endlFile = sys.argv[1]
endlFileName = os.path.split( endlFile )[1]
status = os.system( 'resolutionTestFile.py %s %s' % ( endlFile, sys.argv[2] ) )
if( status != 0 ) : raise

nDataset = int( sys.argv[3] )
axes_ = axes.defaultAxes( labelsUnits = { 0 : ( 'energy_in', 'eV' ), 1 : ( 'crossSection' , 'b' ) } )

def getDataset( n, file ) :

    f = open( file )
    ls = f.readlines( )
    f.close( )
    for i in xrange( n + 1 ) :
        if( len( ls ) < 5 ) : raise Exception( 'file %s only has %d datasets' % ( file, i ) )
        ls = ls[2:]
        j = 0
        for l in ls :
            if( l == '                                                                       1\n' ) : break
            j += 1
        if( i == n ) :
            data = []
            for k in xrange( j ) : data.append( map( float, ls[k].split( ) ) )
            return( data )
        else :
            ls = ls[j+1:]

endlData = XYs.XYs( axes_, getDataset( nDataset, endlFile ), 1e-3 )
endlData.label = endlFileName
lowData = XYs.XYs( axes_, getDataset( nDataset, glob.glob( '%s_T*eV_low' % endlFileName )[0] ), 1e-3 )
lowData.label = 'low res.'
highData = XYs.XYs( axes_, getDataset( nDataset, glob.glob( '%s_T*eV_high' % endlFileName )[0] ), 1e-3 )
highData.label = 'high res.'

fudgeMultiPlots.multiPlot( [ endlData, lowData, highData ] )
