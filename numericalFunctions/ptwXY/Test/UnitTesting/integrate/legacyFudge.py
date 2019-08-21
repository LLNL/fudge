# <<BEGIN-copyright>>
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

