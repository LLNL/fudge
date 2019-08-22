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

# This script requires that groupTestAll1 be run with -v option as it need its *.dat files.

from fudge import *

def printData( data ) :

    for d in data : print "%18.11e" % d

groupBoundaries = endlmisc.read1dDataFile( 'groupBoundaries.dat' )

flux = endl2dmath( endlmisc.read2dDataFile( 'flux.dat' ) )
print 'flux'
print flux
print
flux_grouped = fudge2dGrouping.groupOneFunction( groupBoundaries, flux )
print 'flux grouped'
printData( flux_grouped )

crossSection = endl2dmath( endlmisc.read2dDataFile( 'crossSection.dat' ) )
print 'crossSection'
# print crossSection
print
crossSection_grouped = fudge2dGrouping.groupTwoFunctions( groupBoundaries, crossSection, flux )
print 'cross section grouped'
printData( crossSection_grouped )

multiplicity = endl2dmath( endlmisc.read2dDataFile( 'multiplicity.dat' ) )
print 'multiplicity'
# print multiplicity
print
multiplicity_grouped = fudge2dGrouping.groupThreeFunctions( groupBoundaries, crossSection, multiplicity, flux )
print 'multiplicity grouped'
printData( multiplicity_grouped )

