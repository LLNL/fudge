#! /usr/bin/env python
import argparse
from . import reaclibreaction as reaclib
from fudge.core.utilities.brb import banner

testString='''4
         n    p    d                       ks02     -5.06470E+00          
 3.746377e+00 1.391141e-02-2.672281e+00 1.048889e+01                      
-7.198877e-01 4.782074e-02-3.543090e+00                                   
4
         p    p    d                       bet+w     1.44206e+00          
-3.478630e+01 0.000000e+00-3.511930e+00 3.100860e+00                      
-1.983140e-01 1.262510e-02-1.025170e+00                                   
4
         p    d  he3                       de04      5.49300e+00          
 8.935250e+00 0.000000e+00-3.720800e+00 1.986540e-01                      
 0.000000e+00 0.000000e+00 3.333330e-01                                   
4
         n    d    t                       ks02      6.25730E+00          
 3.301700e+00 7.998334e-05-1.021690e-01 3.713208e+00                      
 1.671925e-01-3.591191e-02-3.749098e-01                                   
4
         n fe56 fe57                       ks02      7.64630E+00          
 3.583399e+01-6.008732e-02 9.794135e+00-3.373665e+01                      
 2.971734e+00-2.364979e-01 1.187107e+01                                   
4
         nau197au198                       ks02      6.51230E+00          
 1.150735e+01 3.266589e-02-4.055318e+00 1.082196e+01                      
 1.578143e-01-2.636530e-01-4.502907e+00                                   
'''

def parse_args():
    parser = argparse.ArgumentParser( description="Monte Carlo Resonances (McRes)" )
    parser.add_argument('--endfFile', default=None, type=str, help='An ENDF file to extract reactions from')
    parser.add_argument('--verbose', dest='verbose', default=False, action='store_true', help='Run verbosely')
    parser.add_argument('--quiet', dest='verbose', default=False, action='store_false', help='Run quietly')
    parser.add_argument('--list', default=False, action="store_true", help="List generated resonances")
    parser.add_argument('--skynetlib', default=False, action='store_true', help="Read in SkyNet's reactlib data")
    parser.add_argument('--reaclibFile', default=None, type=str, help="Read in the specified reaclib file")
    return parser.parse_args()


if __name__=="__main__":
    args=parse_args()
    
    if args.skynetlib:
        print(banner("Reading SkyNetReaclib/reaclib"))
        skynetLib = reaclib.ReaclibLibrary(open('SkyNetReaclib/reaclib').read())
        
    if args.reaclibFile is not None:
        print(banner("Reading user specified reactlib file "+args.reaclibFile))
        userreacLib = reaclib.ReaclibLibrary(open(args.reaclibFile).read())
    else:
        userreacLib = None
  
    if args.endfFile is not None:
        print(banner("Reading user specified ENDF file "+args.endfFile))
        endfReacLib = reaclib.endfFile_to_ReaclibLibrary(args.endfFile)
    else:
        endfReacLib = None
        
    
    
    
    userreacLib = reaclib.ReaclibLibrary(testString)
        
    if True:
        # n+1H -> 2H+g
        print( userreacLib.chapters[4].chapter )
        endfReacLib= reaclib.endfFile_to_ReaclibLibrary('/Users/dbrown/Desktop/endfLibrary/neutrons/n-001_H_001.endf')
        print( userreacLib.chapters[4].reacts[0] )
        reaclib.plot_rates({'ks02':userreacLib.chapters[4].reacts[0], 'endf':endfReacLib[0].reacts[0]}, title="n+1H -> 2H+g")
    
    if False:
        # p+1H->nu+2H
        pass
        # p+1H->nu+2H really goes by p+1H->2He*, one of the p's does this: p->nu + e+ + n, left with 2H
        # we can't do it without including decay library.  
    
    if True:
        # p+2H -> 3He+g
        endfReacLib= reaclib.endfFile_to_ReaclibLibrary('/Users/dbrown/Desktop/endfLibrary/protons/p-001_H_002.endf')
        print( userreacLib.chapters[4].reacts[3] )
        reaclib.plot_rates({'de04':userreacLib.chapters[4].reacts[2], 'endf':endfReacLib[0].reacts[0]}, title="p+2H -> 3He+g")

    if True:
        # n+2H -> 3H+g
        endfReacLib= reaclib.endfFile_to_ReaclibLibrary('/Users/dbrown/Desktop/endfLibrary/neutrons/n-001_H_002.endf')
        print( userreacLib.chapters[4].reacts[3] )
        reaclib.plot_rates({'ks02':userreacLib.chapters[4].reacts[3], 'endf':endfReacLib[0].reacts[0]}, title="n+2H -> 3H+g")
        
    if True:
        # n+56Fe -> 57Fe+g
        endfReacLib= reaclib.endfFile_to_ReaclibLibrary('/Users/dbrown/Desktop/endfLibrary/neutrons/n-026_Fe_056.endf')
        print( userreacLib.chapters[4].reacts[4] )
        reaclib.plot_rates({'ks02':userreacLib.chapters[4].reacts[4], 'endf':endfReacLib[1].reacts[0]}, title="n+56Fe -> 57Fe+g")

    if True:
        # n+197Au -> 198Au+g
        endfReacLib= reaclib.endfFile_to_ReaclibLibrary('/Users/dbrown/Desktop/endfLibrary/neutrons/n-079_Au_197.endf')
        print( userreacLib.chapters[4].reacts[5] )
        reaclib.plot_rates({'ks02':userreacLib.chapters[4].reacts[5], 'endf':endfReacLib[0].reacts[0]}, title="n+197Au -> 198Au+g")
