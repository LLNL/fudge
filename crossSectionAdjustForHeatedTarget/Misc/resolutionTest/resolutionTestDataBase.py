#! /usr/bin/env python

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

#
# This routine loops over each target in the directory given by the first argument, and
# calls resolutionTestTarget.py for target.
#
# USAGE:
#
#   resolutionTestDataBase.py ENDL_projectile_path temperature
#
# OPTIONS:
#
#   ENDL_projectile_path    Path to an ENDL yi projectile directory.
#   temperature             Temperature to heat to in MeV.
#

import sys, os, glob, time

if( len( sys.argv ) != 3 ) : raise Exception( 'need database and temperature' )

files = glob.glob( '%s/za[0-9][0-9][0-9][0-9][0-9][0-9]*' % sys.argv[1] )
files.sort( )
if( len( files ) == 0 ) : raise Exception( 'no targets found in database = %s' % sys.argv[1] )
t0 = time.time( )
for file in files : 
    status = os.system( 'resolutionTestTarget.py %s %s' % ( file, sys.argv[2] ) )
    if( status != 0 ) : raise SystemError()
    os.system( 'make -s clean' )
    t1 = time.time( )
    print(file, '%.2f' % (t1 - t0))
    t0 = t1
