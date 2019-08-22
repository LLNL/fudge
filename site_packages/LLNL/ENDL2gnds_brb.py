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
import glob
from argparse import ArgumentParser

from fudge.legacy.endl import endlZA as endlZAClass
from fudge import fudgeParameters as fudgeParametersModule
fudgeParametersModule.VerboseMode = 0

defaultEvaluationName = "Test"
defaultEvaluationVersion = "1.0"

usage = """
Converts the data in an endlZA directory (e.g., './yi01/za008016') to GNDS. The path to the endlZA directory
must end with a projectile-directory (e.g., 'yi01') followed by a ZA-directory (e.g., 'za008016').
"""

parser = ArgumentParser( description = usage )
parser.add_argument( 'ZAPath', type = str,                                          help = 'path to an endl ZA.' )
parser.add_argument( '-i', '--includeAverageProductData', action = 'store_true',    help = 'Include ENDL average product data (i.e., I 10, 11, and 13 data) in conversion.' )
parser.add_argument( '-e', '--evaluationName', action = 'store', default = defaultEvaluationName,
                                                                                    help = 'Evaluation name. Default is "%s".' % defaultEvaluationName )
parser.add_argument( '-o', '--output', action = 'store', default = None,            help = 'Output file full path. If not entered, ZA-directory with "xml" appended.' )
parser.add_argument( '-V', '--version', action = 'store', default = defaultEvaluationVersion,
                                                                                    help = 'Evaluation version. Default is "%s".' % defaultEvaluationVersion )

args = parser.parse_args( )

crossSectionFiles = glob.glob( args.ZAPath + os.sep + 'yo00c??i000s00?' )
if( len( crossSectionFiles ) == 0 ) : raise Exception( 'No ENDL cross section file in endl ZA directory given.' )

yiPath, ZA = os.path.split( args.ZAPath )
database, yi = os.path.split( yiPath )

endlZA = endlZAClass( ZA, yi, database, readOnly = True )
endlZA.read( )
GNDS = endlZA.toGNDS( args.evaluationName, args.version, excludeAverageProductData = not( args.includeAverageProductData ) )

output = args.output
if( output is None ) : output = ZA + '.xml'
GNDS.saveToFile( output )
