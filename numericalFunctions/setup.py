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

import os, glob, shutil
from distutils.core import setup, Extension

extra_compile_args = [ ]
# Option need for some MACs
# extra_compile_args = [ '-Wno-error=unused-command-line-argument-hard-error-in-future' ]

statusMessageReportingRoot = '../statusMessageReporting'

sep = os.path.sep

statusMessageReporting_c = glob.glob( statusMessageReportingRoot + sep + 'Src' + sep + '*.c' )
ptwC_c = glob.glob( 'ptwC' + sep + 'Src' + sep + '*.c' )
ptwX_c = glob.glob( 'ptwX' + sep + 'Src' + sep + '*.c' )
ptwX_Py_c = glob.glob( 'ptwX' + sep + 'Python' + sep + 'Src' + sep + '*.c' )
ptwXY_c = glob.glob( 'ptwXY' + sep + 'Src' + sep + '*.c' )
ptwXY_Py_c = glob.glob( 'ptwXY' + sep + 'Python' + sep + 'Src' + sep + '*.c' )
nf_Legendre_c = glob.glob( 'nf_Legendre' + sep + 'Src' + sep + '*.c' )
nf_Legendre_Py_c = glob.glob( 'nf_Legendre' + sep + 'Python' + sep + 'Src' + sep + '*.c' )
nf_specialFunctions_c = glob.glob( 'nf_specialFunctions' + sep + 'Src' + sep + '*.c' )
nf_specialFunctions_Py_c = glob.glob( 'nf_specialFunctions' + sep + 'Python' + sep + 'Src' + sep + '*.c' )
nf_integration_c = glob.glob( 'nf_integration' + sep + 'Src' + sep + '*.c' )
nf_integration_Py_c = glob.glob( 'nf_integration' + sep + 'Python' + sep + 'Src' + sep + '*.c' )

statusMessageReporting_hDir = statusMessageReportingRoot + sep + 'Src'
ptwC_hDir = 'ptwC' + sep + 'Src'
ptwX_hDir = 'ptwX' + sep + 'Src'
ptwXY_hDir = 'ptwXY' + sep + 'Src'
ptwXY_Py_hDir = 'ptwXY' + sep + 'Python' + sep + 'Src'
nf_Legendre_hDir = 'nf_Legendre' + sep + 'Src'
nf_specialFunctions_hDir = 'nf_specialFunctions' + sep + 'Src'
nf_specialFunctions_Py_hDir = 'nf_specialFunctions' + sep + 'Python' + sep + 'Src'
nf_integration_hDir = 'nf_integration' + sep + 'Src'

#
# Stuff to build listOfDoubles.so.
#
libs = glob.glob( 'build' + sep + 'lib*' + sep + 'listOfDoubles.so' )
for lib in libs : os.remove( lib )

listOfDoubles_C = Extension( 'listOfDoubles_C', \
    extra_compile_args = extra_compile_args, \
    sources = statusMessageReporting_c + ptwC_c + ptwX_c + ptwX_Py_c, \
    include_dirs = [ statusMessageReporting_hDir, ptwC_hDir, ptwX_hDir ] )

setup( name = 'listOfDoubles_C', 
    version = '1.0', \
    description = 'This module contains the listOfDoubles_C class and support routines.', \
    ext_modules = [ listOfDoubles_C ] )

libs = glob.glob( 'build' + sep + 'lib*' + sep + 'listOfDoubles_C.so' )
if( len( libs ) > 0 ) : shutil.copyfile( libs[0], os.path.join( 'lib', 'listOfDoubles_C.so' ) )

#
# Stuff to build pointwiseXY_C.so.
#
libs = glob.glob( 'build' + sep + 'lib*' + sep + 'pointwiseXY_C.so' )
for lib in libs : os.remove( lib )

pointwiseXY_C = Extension( 'pointwiseXY_C', \
    extra_compile_args = extra_compile_args, \
    sources = statusMessageReporting_c + ptwC_c + ptwX_c + nf_Legendre_c + nf_integration_c + ptwXY_c + ptwXY_Py_c, \
    include_dirs = [ statusMessageReporting_hDir, ptwC_hDir, ptwX_hDir, nf_Legendre_hDir, nf_integration_hDir, 
            ptwXY_hDir, ptwXY_Py_hDir ] )

setup( name = 'pointwiseXY_C', 
    version = '1.0', \
    description = 'This module contains the pointwiseXY_C class and support routines.', \
    ext_modules = [ pointwiseXY_C ] )

libs = glob.glob( 'build' + sep + 'lib*' + sep + 'pointwiseXY_C.so' )
if( len( libs ) > 0 ) : shutil.copyfile( libs[0], os.path.join( 'lib', 'pointwiseXY_C.so' ) )

#
# Stuff to build Legendre.so.
#
libs = glob.glob( 'build' + sep + 'lib*' + sep + 'Legendre.so' )
for lib in libs : os.remove( lib )

nf_Legendre_C = Extension( 'Legendre', \
    extra_compile_args = extra_compile_args, \
    sources = statusMessageReporting_c + ptwC_c + ptwX_c + ptwXY_c + nf_integration_c + nf_Legendre_c + nf_Legendre_Py_c, \
    include_dirs = [ statusMessageReporting_hDir, ptwC_hDir, ptwX_hDir, ptwXY_hDir, ptwXY_Py_hDir, 
            nf_integration_hDir, nf_Legendre_hDir ] )

setup( name = 'Legendre', 
    version = '1.0', \
    description = 'This module contains the Legendre class and support routines.', \
    ext_modules = [ nf_Legendre_C ] )

libs = glob.glob( 'build' + sep + 'lib*' + sep + 'Legendre.so' )
if( len( libs ) > 0 ) : shutil.copyfile( libs[0], os.path.join( 'lib', 'Legendre.so' ) )

#
# Stuff to build specialFunctions.so
#
libs = glob.glob( 'build' + sep + 'lib*' + sep + 'specialFunctions.so' )
for lib in libs : os.remove( lib )

specialFunctions = Extension( 'specialFunctions', \
    extra_compile_args = extra_compile_args, \
    sources = statusMessageReporting_c + ptwC_c + nf_specialFunctions_c + nf_specialFunctions_Py_c, \
    include_dirs = [ statusMessageReporting_hDir, ptwC_hDir, nf_specialFunctions_hDir, nf_specialFunctions_Py_hDir ] )

setup( name = 'specialFunctions', 
    version = '1.0', \
    description = 'This module contains some special math functions not in the python math module.', \
    ext_modules = [ specialFunctions ] )

libs = glob.glob( 'build' + sep + 'lib*' + sep + 'specialFunctions.so' )
if( len( libs ) > 0 ) : shutil.copyfile( libs[0], os.path.join( 'lib', 'specialFunctions.so' ) )

#
# Stuff to build integration.so
#
libs = glob.glob( 'build' + sep + 'lib*' + sep + 'integration.so' )
for lib in libs : os.remove( lib )

integration = Extension( 'integration', \
    extra_compile_args = extra_compile_args, \
    sources = statusMessageReporting_c + ptwC_c + ptwX_c + ptwXY_c + nf_Legendre_c + nf_integration_c + nf_integration_Py_c, \
    include_dirs = [ statusMessageReporting_hDir, ptwC_hDir, ptwX_hDir, ptwXY_hDir, nf_Legendre_hDir, nf_integration_hDir ] )

setup( name = 'integration', 
    version = '1.0', \
    description = 'This module contains functions for integrating a function representing an integrand.', \
    ext_modules = [ integration ] )

libs = glob.glob( 'build' + sep + 'lib*' + sep + 'integration.so' )
if( len( libs ) > 0 ) : shutil.copyfile( libs[0], os.path.join( 'lib', 'integration.so' ) )
