# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import math

from xData import axes as axesModule
from xData import XYs1d as XYs1dModule

from fudge.productData.distributions import angular as angularModule
from fudge.productData.distributions import uncorrelated as uncorrelatedModule

from . import base as baseModule

def calculatePofMu( energy, index, SofE, epsilon ) :

    def addDeltaFunction( PofMu, energy, energy_i, S_i, S_im1, epsilon ) :

        dS = S_i - S_im1
        if( dS == 0 ) : return
        mu_i = 1 - 2 * energy_i / energy
        delta = abs( epsilon * mu_i )
        if( delta == 0 ) : delta = epsilon * 1e-3
        mu_lower = mu_i - delta
        mu_upper = mu_i + delta
        if( mu_lower < -1 ) :
            mu_lower = -1
            mu_i = mu_lower + delta
            mu_upper = mu_i + delta
        elif( mu_upper > 1 ) :
            mu_upper = 1
            mu_i = mu_upper - delta
            mu_lower = mu_i - delta

        if( len( PofMu ) > 0 ) :
            if( PofMu[-1][0] <= mu_lower ) :
                raise ValueError( 'Delta function overlap at energy = %s and energy_i = %s: %s and %s' % 
                        ( energy, energy_i, PofMu[-1][0], mu_i ) )
        PofMu.insert( 0, [ mu_upper, 0 ] )
        PofMu.insert( 0, [ mu_i,     dS / delta ] )
        PofMu.insert( 0, [ mu_lower, 0 ] )

    if( index == 0 ) : index = 1

    PofMu = []
    S_im1 = 0
    for i1 in range( index ) :
        energy_i, S_i = SofE[i1]
        addDeltaFunction( PofMu, energy, energy_i, S_i, S_im1, epsilon )
        S_im1 = S_i

    if( len( PofMu ) > 0 ) :
        if( PofMu[0][0]  > -1 ) : PofMu.insert( 0, [ -1, 0 ] )
        if( PofMu[-1][0] <  1 ) : PofMu.append(    [  1, 0 ] )
    return( PofMu )

def bisectionTest( crossSectionFine, x1, y1, x2, y2, accuracy, level = 0 ) :

    y2p = y1 * x1 / x2
    xMid = 0.5 * ( x1 + x2 )
    yMid = y1 * x1 / xMid
    yEstimate = 0.5 * ( y1 + y2p )
    if( abs( yEstimate - yMid ) <= accuracy * yMid ) : return
    bisectionTest( crossSectionFine, x1,   y1,   xMid, yMid, accuracy, level + 1 )
    crossSectionFine.append( [ xMid, yMid ] )
    bisectionTest( crossSectionFine, xMid, yMid, x2,   y2,   accuracy, level + 1 )

def process( self, label, energyMin, energyMax, temperature, kwargs ) :

    accuracy = min( 0.1, max( 1e-5, kwargs['accuracy'] ) )
    epsilon = min( 1e-3, max( 1e-15, kwargs['epsilon'] ) )
    momentumConstant = math.sqrt( 2.0 * kwargs['neutronMass'] )

    atomsPerMolecule = baseModule.getNumberOfAtomsPerMolecule(self)

    S_table = self.S_table
    gridded2d = S_table.gridded2d

    array = gridded2d.array.constructArray( )

    energyGrid = gridded2d.axes[1].copy()

    temperatureGrid = gridded2d.axes[2].copy()
    temperatureGrid.convertToUnit( kwargs['temperatureUnit'] )
    for index2, temperature2 in enumerate( temperatureGrid.values ) :
        if( temperature2 >= temperature ) : break

    def expandGrid(energyGrid, array):
        """
        Add energyMax to the energy if not already present. energyMin dealt with later
        """
        Xs = list(energyGrid.values)
        Ys = list(array)
        if Xs[0] > energyMin:
            if Xs[0] <= energyMin * (1+epsilon):
                Xs[0] = energyMin
            else:
                Xs.insert(0, energyMin)
                Ys.insert(0, 0)
        if Xs[-1] < energyMax:
            if Xs[-1] >= energyMax * (1-epsilon):
                Xs[-1] = energyMax
            else:
                Xs.append(energyMax)
                Ys.append(Ys[-1])
        return XYs1dModule.XYs1d((Xs, Ys), dataForm="XsAndYs", interpolation=energyGrid.interpolation, axes=axesModule.Axes(2))

    if( temperature > temperature2 ) :
        SofE = expandGrid(energyGrid, array[-1])
    elif( ( index2 == 0 ) or ( temperature2 == temperature ) ) :
        SofE = expandGrid(energyGrid, array[index2])
    else :
        index1 = index2 - 1
        temperature1 = temperatureGrid.values[index1]
        fraction = ( temperature - temperature1 ) / ( temperature2 - temperature1 )

        SofE1 = expandGrid(energyGrid, array[index1])
        SofE2 = expandGrid(energyGrid, array[index2])
        SofE = ( 1 - fraction ) * SofE1 + fraction * SofE2

    _SofE = SofE
    SofE = []
    for energy, S2 in _SofE :
        if( energy > energyMax ) : break
        SofE.append( [ energy, S2 ] )
    if( ( SofE[-1][0] < energyMax ) and ( SofE[-1][0] < _SofE[-1][0] ) ) :
        eMax = _SofE[-1][0]
        if( energyMax < eMax ) : eMax = energyMax
        SofE.append( [ eMax, SofE[-1][1] ] )

    S1 = 0
    crossSection = []
    averageProductEnergy = []
    averageProductMomentum = []
    angularData = []
    n_minus1 = len( SofE ) - 1
    for index, energyS in enumerate( SofE ) :
        energy, S2 = energyS

        energy1 = ( 1 - epsilon ) * energy
        if( index == 0 ) :
            S1 = S2
            energy1 = energy

        energy2 = ( 1 + epsilon ) * energy
        if( index == n_minus1 ) : energy2 = energy

        if( energy >= energyMax ) :
            energy = energyMax
            energy2 = energyMax
        if( energy2 > energyMax ) :
            energy2 = energyMax

        crossSection.append( [ energy1, S1 / energy ] )
        if( energy != energy1 ) : crossSection.append( [ energy2, S2 / energy ] )

        averageProductEnergy.append( [ energy1, S1 ] )
        if( energy != energy1 ) : averageProductEnergy.append( [ energy2, S2 ] )

        momentumFactor = momentumConstant / math.sqrt( energy )
        averageProductMomentum.append( [ energy1, momentumFactor * S1 ] )                              # This is not correct. Missing mu_i.
        if( energy != energy1 ) : averageProductMomentum.append( [ energy2, momentumFactor * S2 ] )    # This is not correct. Missing mu_i.

        PofMu = calculatePofMu( energy1, index,     SofE, epsilon )
        if( len( PofMu ) > 0 ) : angularData.append( [ energy1, PofMu ] )
        if( energy != energy1 ) :
            PofMu = calculatePofMu( energy2, index + 1, SofE, epsilon )
            if( len( PofMu ) > 0 ) : angularData.append( [ energy2, PofMu ] )

        S1 = S2
        if( energy == energyMax ) : break

    crossSectionFine = [ crossSection.pop( 0 ) ]
    x1, y1 = crossSectionFine[0]
    for x2, y2 in crossSection :
        if( ( x2 - x1 ) > ( 2 * epsilon * x1 ) ) : bisectionTest( crossSectionFine, x1, y1, x2, y2, accuracy )
        crossSectionFine.append( [ x2, y2 ] )
        x1, y1 = x2, y2

    crossSection, averageProductEnergy, averageProductMomentum = baseModule.processedToForms( label, crossSectionFine, averageProductEnergy, averageProductMomentum, kwargs )
    crossSection /= atomsPerMolecule

    incidentEnergyUnit = kwargs['incidentEnergyUnit']
    axes = angularModule.defaultAxes( incidentEnergyUnit )
    angular2d = angularModule.XYs2d( axes = axes )
    for energy, PofMu in angularData :
        angular2d.append( angularModule.XYs1d( data = PofMu, outerDomainValue = energy, axes = axes ).normalize( ) )
    angular2d = uncorrelatedModule.AngularSubform( angular2d )

    energy2d = baseModule.energyDelta2d( crossSection[0][0], crossSection[-1][0], epsilon, incidentEnergyUnit )
    energy2d = uncorrelatedModule.EnergySubform( energy2d )

    distribution = uncorrelatedModule.Form( label, self.productFrame, angular2d, energy2d )

    return( crossSection, averageProductEnergy, averageProductMomentum, distribution )
