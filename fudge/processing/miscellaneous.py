# <<BEGIN-copyright>>
# <<END-copyright>>

def _groupFunctionsAndFluxInit( projectile, processInfo, f1 ) :

    groupBoundaries = processInfo.getParticleGroups( projectile )
    flux = processInfo['flux']['data'][0]
    xMin, xMax = f1.getDomain( )
    return( groupBoundaries, flux.xSlice( xMin, xMax ) )

def _mutualifyGrouping3Data( f1, f2, f3 ) :

    xMin1, xMax1 = f1.getDomain( )
    xMin2, xMax2 = f2.getDomain( )
    xMin3, xMax3 = f2.getDomain( )
    if( ( xMin1 != xMin2 ) or ( xMin1 != xMin3 ) or ( xMax1 != xMax2 ) or ( xMax1 != xMax3 ) ) :
        xMin, xMax = max( xMin1, xMin2, xMin3 ), min( xMax1, xMax2, xMax3 )
        print "WARNING: making domains mutual for grouping,", xMin1, xMin2, xMin3, xMax1, xMax2, xMax3
        f1, f2, f3 = f1.xSlice( xMin = xMin, xMax = xMax ), f2.xSlice( xMin = xMin, xMax = xMax ), f3.xSlice( xMin = xMin, xMax = xMax )
    return( f1, f2, f3 )

def groupOneFunctionsAndFlux( projectile, processInfo, f1, norm = None ) :

    groupBoundaries, flux_ = _groupFunctionsAndFluxInit( projectile, processInfo, f1 )
    return( f1.groupTwoFunctions( groupBoundaries, flux_, norm = norm ) )

def groupTwoFunctionsAndFlux( projectile, processInfo, f1, f2, norm = None ) :

    groupBoundaries, flux_ = _groupFunctionsAndFluxInit( projectile, processInfo, f1 )
    f1, f2, flux_ = _mutualifyGrouping3Data( f1, f2, flux_ )
    return( f1.groupThreeFunctions( groupBoundaries, flux_, f2, norm = norm ) )
