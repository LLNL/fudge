#! /usr/bin/env python

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
import os
from fudge.legacy.endl import endlProject
import site_packages.legacy.toENDF6.toENDF6 # adds 'toENDF6' methods to GND classes

import argparse
parser = argparse.ArgumentParser( description = "Translate ENDL evaluations into GND, and optionally also ENDF-6." )
parser.add_argument( 'ZA', type = str, help = "za of target to translate" )
parser.add_argument( '-l', '--library', default = 'endl2011.0',
            help = "ENDL library (default=endl2011.0). May be full path or identifier" )
parser.add_argument( '-p', '--projectile', default = 'n',
            type = lambda val: int(val) if val.isdigit() else val,
            help = "Choose projectile (default='n'). May be string (one of n,p,d,t,He3,a,g) or yi number" )
parser.add_argument( '--includeAverageProductData', action='store_true',
            help = "Include I=10 and 13 (average outgoing energy and momentum) data" )
parser.add_argument( '-6', '--toENDF6', action = 'store_true',
            help = "After creating GND, also translate to ENDF-6" )
parser.add_argument( '-o', '--output', default='endl2gnd', help='prefix for resulting .endf and .xml files' )
parser.add_argument( '-v', '--version', default = '1.0.0',
            help = "Evaluation version number. Should be of form 'major.minor.patchlevel' (i.e. 2011.0.1)" )

if( __name__ == '__main__' ) :

    args = parser.parse_args( )

    e = endlProject( args.library, projectile = args.projectile, readOnly = True )

    if args.ZA[-1].isalpha():   # isomer, e.g. 95242m
        za = e.readZA( args.ZA[:-1], suffix = args.ZA[-1] )
    else:
        za = e.readZA( args.ZA )
    za.read ()
    if not args.includeAverageProductData:
        za.removeFile(I=10)
        za.removeFile(I=13)

    r = za.toGND( evaluationLibrary = os.path.basename( args.library ), evaluationVersion = args.version )
    r.saveToFile( args.output + ".xml" )

    if( args.toENDF6 ) :
        r.convertUnits( {'MeV':'eV'} )
        with open( args.output + ".endf", "w" ) as fout :
            fout.write( r.toENDF6( style = "eval", flags = { 'verbosity' : 32 } ) )
