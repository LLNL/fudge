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
import argparse

from fudge.lib import GNDSType as GNDSTypeModule
from PoPs import database as databaseModule

description1 = """Read one or more GNDS files into Fudge, then write back to GNDS/xml.  Intent is to test 
    for errors during reading or writing.  Sample use: python gnds2gnds.py n-001_H_001.xml n-001_H_002.xml ...
    If file n-001_H_001-cov.xml (or -covar.xml) exists, covariances will automatically be read and re-written.
    Resulting files will be in the local directory, and will have extension '.g2g' appended.
"""

__doc__ = description1

parser = argparse.ArgumentParser( description1 )
parser.add_argument( 'gnds', nargs = '+',                                       help = 'GNDS and/or PoPs file(s) to translate' )
parser.add_argument( '-o', '--outline', default = False, action = 'store_true', help = 'The outputted GNDS files are written in outline mode' )
parser.add_argument( '--energyUnit', type = str, default = None,                help = 'Convert all energies in the gnds file to this unit.' )

if __name__ == '__main__' :
    args = parser.parse_args()

    for fileName in args.gnds :

        covariance = None
        name, dummy = GNDSTypeModule.type( fileName )
        if( name == databaseModule.database.moniker ) :
            gnds = GNDSTypeModule.read( fileName )
        else :
            covarianceFileName = None
            for extension in ( '-cov.xml', '-covar.xml' ) :
                tmp = fileName.replace( '.xml', extension )
                if( os.path.exists( tmp ) ) : covarianceFileName = tmp

            gnds = GNDSTypeModule.read( fileName )
            if( covarianceFileName is not None ) :
                covariance = GNDSTypeModule.read( covarianceFileName, reactionSuite = gnds )

        if( args.energyUnit is not None ) :
            gnds.convertUnits( { 'MeV' : args.energyUnit, 'eV' : args.energyUnit } )
            if( covariance is not None ) : covariance.convertUnits( { 'MeV' : args.energyUnit, 'eV' : args.energyUnit } )
        gnds.saveToFile( os.path.basename( fileName ) + '.g2g', outline = args.outline, xs_pdf_cdf1d_singleLine = True )
        if( covariance is not None ) : covariance.saveToFile( os.path.basename( covarianceFileName ) + '.g2g' )
