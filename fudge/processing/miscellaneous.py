# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

accuracy = 1e-6

from xData import enums as xDataEnumsModule
from xData import link as linkModule
from xData import axes as axesModule
from xData import XYs1d as XYs1dModule
from xData import regions as regionsModule

from . import group as groupModule

def _toLinear( func ) :

    if( isinstance( func, XYs1dModule.XYs1d ) ) :
        if func.interpolation == xDataEnumsModule.Interpolation.linlin:
            return func
        func = func.copy( )
        return( func.toPointwise_withLinearXYs( accuracy = accuracy, upperEps = 1e-8 ) )
    elif( isinstance( func, regionsModule.Regions1d ) ) :
        return( func.toPointwise_withLinearXYs( accuracy = accuracy, upperEps = 1e-8 ) )

    print(type(func))
    raise Exception( 'FIX ME' )

def _groupFunctionsAndFluxInit( style, tempInfo, f1 ) :

    reactionSuite = tempInfo['reactionSuite']
    groupBoundaries = style.transportables[reactionSuite.projectile].group.boundaries
    flux = style.flux.getFluxAtLegendreOrder( 0 )
    if( f1 is not None ) :
        domainMin, domainMax = f1.domainMin, f1.domainMax
        flux = flux.domainSlice( domainMin, domainMax )
    return( groupBoundaries, flux )

def _mutualifyGrouping3Data( f1, f2, f3, printMutualDomainWarning = False ) :

    domainMin1, domainMax1 = f1.domainMin, f1.domainMax
    domainMin2, domainMax2 = f2.domainMin, f2.domainMax
    if( len( f3 ) == 0 ) :
        domainMin3, domainMax3 = domainMin2, domainMax2
    else :
        domainMin3, domainMax3 = f3.domainMin, f3.domainMax
    if( ( domainMin1 != domainMin2 ) or ( domainMin1 != domainMin3 ) or ( domainMax1 != domainMax2 ) or ( domainMax1 != domainMax3 ) ) :
        domainMin, domainMax = max( domainMin1, domainMin2, domainMin3 ), min( domainMax1, domainMax2, domainMax3 )
        if( printMutualDomainWarning ) :
            print("WARNING: making domains mutual for grouping,", domainMin1, domainMin2, domainMin3, domainMax1, domainMax2, domainMax3)
        f1 = f1.domainSlice( domainMin = domainMin, domainMax = domainMax )
        f2 = f2.domainSlice( domainMin = domainMin, domainMax = domainMax )
        if( len( f3 ) > 0 ) : f3 = f3.domainSlice( domainMin = domainMin, domainMax = domainMax )
    return( f1, f2, f3 )

def groupOneFunctionAndFlux( style, tempInfo, f1, styleFilter = None ) :

    if( isinstance( f1, linkModule.Link ) ) : f1 = style.findFormMatchingDerivedStyle( f1.link, styleFilter = styleFilter )
    f1 = _toLinear( f1 )
    groupBoundaries, flux = _groupFunctionsAndFluxInit( style, tempInfo, f1 )
    return( f1.groupTwoFunctions( groupBoundaries, flux, norm = tempInfo['groupedFlux'] ) )

def groupTwoFunctionsAndFlux( style, tempInfo, f1, f2, norm = None, printMutualDomainWarning = False ) :

    f1 = _toLinear( f1 )
    f2 = _toLinear( f2 )
    groupBoundaries, flux = _groupFunctionsAndFluxInit( style, tempInfo, f1 )
    f1, f2, flux = _mutualifyGrouping3Data( f1, f2, flux, printMutualDomainWarning = printMutualDomainWarning )
    return( f1.groupThreeFunctions( groupBoundaries, flux, f2, norm = norm ) )

def groupFunctionCrossSectionAndFlux( cls, style, tempInfo, f1, printMutualDomainWarning = False ) :

    crossSection = tempInfo['crossSection']
    norm = tempInfo['groupedFlux']
    grouped = groupTwoFunctionsAndFlux( style, tempInfo, f1, crossSection, norm = norm, printMutualDomainWarning = printMutualDomainWarning )

    unit = crossSection.axes[0].unit
    if( f1.axes[0].unit != '' ) : unit = '%s * %s' % ( f1.axes[0].unit, unit )
    axis = axesModule.Axis( label = '%s * %s' % ( f1.axes[0].label, crossSection.axes[0].label ), unit = unit, index = 0 )
    axes = f1.axes.copy( )
    axes[0] = axis
    return( groupModule.toMultiGroup1d( cls, style, tempInfo, axes, grouped, zeroPerTNSL = tempInfo['zeroPerTNSL'] ) )
