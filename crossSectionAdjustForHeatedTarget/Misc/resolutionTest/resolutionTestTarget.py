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
