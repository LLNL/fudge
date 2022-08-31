#! /usr/bin/env python

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

#
# This routine loops over each file in the target given by the first argument, and
# calls resolutionTestFile.py on that file.
#
# USAGE:
#
#  resolutionTestTarget.py ENDL_target_path temperature
#
# OPTIONS:
#
#   ENDL_target_path    Path to an ENDL ZA target directory.
#   temperature         Temperature to heat to in MeV.
#
import sys, os, glob

if( len( sys.argv ) != 3 ) : raise Exception( 'need target and temperature' )

files = glob.glob( '%s/yo[0-9][0-9]c[0-9][0-9]i000s[0-9][0-9][0-9]' % sys.argv[1] )
if( len( files ) == 0 ) : raise Exception( 'no cross section files found in target = %s' % sys.argv[1] )

for file in files : 
    status = os.system( 'resolutionTestFile.py %s %s' % ( file, sys.argv[2] ) )
    if( status != 0 ) : raise
