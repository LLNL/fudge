# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse, fudge
import fudge.covariances.covarianceSuite
import fudge.processing.resonances.reconstructResonances

import brownies.legacy.converting.endfFileToGNDS


def parseArgs():
    '''Parse the command line arguments'''
    parser = argparse.ArgumentParser(description='Plot nuclear data from an ENDF or GNDS file')
    parser.add_argument('--mt',     metavar='mt', type=int, help='MT of the angular distribution to work on.' )
    parser.add_argument('--check',  default=False, action='store_true', help="Check the resulting evaluation")
    parser.add_argument('--format', choices=['gnds','endf','table'], default='gnds', help="Format to save the results in outFile, with results from --strategy switch")
    parser.add_argument('--strategy', choices=['merge','replace','dryrun'], default='dryrun', help="How to merge with existing fast angular distribution")
    parser.add_argument('endf',     metavar='endf', type=str, help='ENDF file' )
    parser.add_argument('-o',       dest='outFile', default=None, type=str, help='Output file' )
    return parser.parse_args()

def smallBanner( x, wingsize=10 ): 
    '''Small banner, with wings'''
    return wingsize*'*'+' '+x.replace( '\n', '; ' )+' '+wingsize*'*'


def getPlotTable( nativeData ):
    '''Make a table that we can feed to a spreadsheet from an instance of fudge.productData.distributions.angular.LegendrePointwise'''
    if not isinstance( nativeData, fudge.productData.distributions.angular.LegendrePointwise): raise TypeError("Must be of type fudge.productData.distributions.angular.LegendrePointwise")
    table = []
    for x in nativeData: 
        table.append([0.0]*(2+nativeData.maxLegendreOrder()))
        table[-1][0]=x.value
        for i,y in enumerate(x): table[-1][i+1]=y
    return [ '\t'.join( map(str,row) ) for row in table ]


if __name__ == "__main__":

    args = parseArgs()

    gndsMap = {}

    print( smallBanner("Reading %s" % args.endf) )

    # Is the file a GNDS file?
    if open(args.endf).readline().startswith( "<?xml" ):
        RS = fudge.reactionSuite.readXML( args.endf)
        try:    CS = fudge.covariances.covarianceSuite.readXML( args.endf.replace('.gnds.', '.gndsCov.'))
        except: CS = fudge.covariances.covarianceSuite()
        gndsMap[args.endf] = [ RS, CS ]

    # Maybe its an ENDF file?
    elif open(args.endf).readline().endswith(' 0  0    0\n'): 
        results = brownies.legacy.converting.endfFileToGNDS.endfFileToGNDS(args.endf, toStdOut = False, skipBadData = True)
        if type( results ) == dict:    gndsMap[args.endf] = [ results['reactionSuite'], results['covarianceSuite'] ]
        elif type( results ) == tuple: gndsMap[args.endf] = [ results[0], results[1] ]
        else: 
            raise TypeError( "endfFileToGNDS.endfFileToGNDS() returned a "+str(type(results))+", I don't know what to do with it" )

    # Failed!
    else: print( "WARNING: Unknown file type, not reading %s"% args.endf )
                                
    # Reconstruct the angular distributions
    # The thing returned as the value in the angDists map should be a valid 
    # component that we can just attach to the product (possibly after matching onto an existing table)
    if args.strategy == 'dryrun': exit()
    print( smallBanner( "Reconstructing resonance cross sections" ) )
    gndsMap[args.endf][0].reconstructResonances()
    print( smallBanner( "Reconstructing resonance angular distributions" ) )
    newAngDists = fudge.processing.resonances.reconstructResonances.reconstructAngularDistributions( gndsMap[args.endf][0] )
    
    
    # If we reconstructed the angular distributions too, put them in the appropriate 
    # reaction/outputChannel/product/distributions/angular/LegendrePointwise
    # Right now only elastic is supported
    product = gndsMap[args.endf][0].getReaction( 'elastic' ).outputChannel.getProductWithName( 'n' )
    if args.strategy == 'replace': product.addDistributionComponent( angDists['elastic'] )
    elif args.strategy =='merge':  
        print( smallBanner( "Merging angular distributions" ) )
        originalNativeData = product.getDistributionComponentByToken( product.getDistributionNativeData(  ) ).getNativeData()
        newNativeData = newAngDists['elastic'].getNativeData()
        Emax = newNativeData.domainMax()
        if str( originalNativeData ).endswith( 'LegendrePointwise' ):
            for x in originalNativeData:
                if x.value > Emax: newNativeData.append(x)
#            print('DEBUG: new',len(newNativeData) )
            product.addDistributionComponent( newAngDists['elastic'] )
#            print('DEBUG: old',len(originalNativeData) )
#            print('DEBUG: should  be new',len(product.getDistributionComponentByToken( product.getDistributionNativeData(  ) ).getNativeData()) )
        else: raise ValueError( "Don't know how to merge LegendrePointwise from reconstructed resonances with angular distribution of type %s"%str( nativeData ).split('/')[-1])
    else: pass    

    # Check evaluation
    if args.check:
        print( smallBanner( "Checking final evaluation" ) )
        print( gndsMap[args.endf][0].check() )

    # Save the results
    if args.outFile is not None:
        print( smallBanner( "Saving results as a "+args.format+" in file %s"%args.outFile ) )
        if args.format == 'table': 
            open( args.outFile, mode='w' ).write( 
                '\n'.join( 
                    getPlotTable( product.getDistributionComponentByToken( product.getDistributionNativeData(  ) ).getNativeData() )))
        elif args.format == 'endf':
            open( args.outFile, mode='w' ).write( gndsMap[args.endf][0].toENDF6(flags={'verbosity':0},covarianceSuite=gndsMap[args.endf][1]) )
        else: 
            gndsMap[args.endf][0].saveToFile( args.outFile )
    
    
    
    
