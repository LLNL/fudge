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

import argparse, fudge
import fudge.legacy.converting.endfFileToGND
import fudge.gnd.covariances.covarianceSuite
import fudge.processing.resonances.reconstructResonances


def parseArgs():
    '''Parse the command line arguments'''
    parser = argparse.ArgumentParser(description='Plot nuclear data from an ENDF or GND file')
    parser.add_argument('--mt',     metavar='mt', type=int, help='MT of the angular distribution to work on.' )
    parser.add_argument('--check',  default=False, action='store_true', help="Check the resulting evaluation")
    parser.add_argument('--format', choices=['gnd','endf','table'], default='gnd', help="Format to save the results in outFile, with results from --strategy switch")
    parser.add_argument('--strategy', choices=['merge','replace','dryrun'], default='dryrun', help="How to merge with existing fast angular distribution")
    parser.add_argument('endf',     metavar='endf', type=str, help='ENDF file' )
    parser.add_argument('-o',       dest='outFile', default=None, type=str, help='Output file' )
    return parser.parse_args()

def smallBanner( x, wingsize=10 ): 
    '''Small banner, with wings'''
    return wingsize*'*'+' '+x.replace( '\n', '; ' )+' '+wingsize*'*'


def getPlotTable( nativeData ):
    '''Make a table that we can feed to a spreadsheet from an instance of fudge.gnd.productData.distributions.angular.LegendrePointwise'''
    if not isinstance( nativeData, fudge.gnd.productData.distributions.angular.LegendrePointwise ): raise TypeError( "Must be of type fudge.gnd.productData.distributions.angular.LegendrePointwise")
    table = []
    for x in nativeData: 
        table.append([0.0]*(2+nativeData.maxLegendreOrder()))
        table[-1][0]=x.value
        for i,y in enumerate(x): table[-1][i+1]=y
    return [ '\t'.join( map(str,row) ) for row in table ]


if __name__ == "__main__":

    args = parseArgs()

    gndMap = {}

    print smallBanner( "Reading %s" % args.endf )

    # Is the file a GND file?
    if open(args.endf).readline().startswith( "<?xml" ):
        RS = fudge.gnd.reactionSuite.readXML( args.endf )
        try:    CS = fudge.gnd.covariances.covarianceSuite.readXML( args.endf.replace( '.gnd.', '.gndCov.' ) )
        except: CS = fudge.gnd.covariances.covarianceSuite()
        gndMap[args.endf] = [ RS, CS ]

    # Maybe its an ENDF file?
    elif open(args.endf).readline().endswith(' 0  0    0\n'): 
        results = fudge.legacy.converting.endfFileToGND.endfFileToGND( args.endf, toStdOut = False, skipBadData = True ) 
        if type( results ) == dict:    gndMap[args.endf] = [ results['reactionSuite'], results['covarianceSuite'] ]
        elif type( results ) == tuple: gndMap[args.endf] = [ results[0], results[1] ]
        else: 
            raise TypeError( "endfFileToGND.endfFileToGND() returned a "+str(type(results))+", I don't know what to do with it" )

    # Failed!
    else: print "WARNING: Unknown file type, not reading %s"% args.endf
                                
    # Reconstruct the angular distributions
    # The thing returned as the value in the angDists map should be a valid 
    # component that we can just attach to the product (possibly after matching onto an existing table)
    if args.strategy == 'dryrun': exit()
    print smallBanner( "Reconstructing resonance cross sections" )
    gndMap[args.endf][0].reconstructResonances()
    print smallBanner( "Reconstructing resonance angular distributions" )
    newAngDists = fudge.processing.resonances.reconstructResonances.reconstructAngularDistributions( gndMap[args.endf][0] )
    
    
    # If we reconstructed the angular distributions too, put them in the appropriate 
    # reaction/outputChannel/product/distributions/angular/LegendrePointwise
    # Right now only elastic is supported
    product = gndMap[args.endf][0].getReaction( 'elastic' ).outputChannel.getProductWithName( 'n' )
    if args.strategy == 'replace': product.addDistributionComponent( angDists['elastic'] )
    elif args.strategy =='merge':  
        print smallBanner( "Merging angular distributions" )
        originalNativeData = product.getDistributionComponentByToken( product.getDistributionNativeData(  ) ).getNativeData()
        newNativeData = newAngDists['elastic'].getNativeData()
        Emax = newNativeData.domainMax()
        if str( originalNativeData ).endswith( 'LegendrePointwise' ):
            for x in originalNativeData:
                if x.value > Emax: newNativeData.append(x)
#            print 'DEBUG: new',len(newNativeData)
            product.addDistributionComponent( newAngDists['elastic'] )
#            print 'DEBUG: old',len(originalNativeData)
#            print 'DEBUG: should  be new',len(product.getDistributionComponentByToken( product.getDistributionNativeData(  ) ).getNativeData())
        else: raise ValueError( "Don't know how to merge LegendrePointwise from reconstructed resonances with angular distribution of type %s"%str( nativeData ).split('/')[-1])
    else: pass    

    # Check evaluation
    if args.check:
        print smallBanner( "Checking final evaluation" )
        print gndMap[args.endf][0].check()

    # Save the results
    if args.outFile is not None:
        print smallBanner( "Saving results as a "+args.format+" in file %s"%args.outFile )
        if args.format == 'table': 
            open( args.outFile, mode='w' ).write( 
                '\n'.join( 
                    getPlotTable( product.getDistributionComponentByToken( product.getDistributionNativeData(  ) ).getNativeData() )))
        elif args.format == 'endf':
            open( args.outFile, mode='w' ).write( gndMap[args.endf][0].toENDF6(flags={'verbosity':0},covarianceSuite=gndMap[args.endf][1]) )
        else: 
            gndMap[args.endf][0].saveToFile( args.outFile )
    
    
    
    
