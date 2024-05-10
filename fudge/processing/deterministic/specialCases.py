# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import math

def twoBodyPhotoNuclearAtThreshold(masses, Q, modifiedProductGroupBoundaries, crossSection, angularData, legendreMax, groupFlux, firstCall = True):
    """
    This function calculates the part of the tranfer matrix at the threshold energy for the outgoing particle (product)
    groups specified by **modifiedProductGroupBoundaries**. This function exists because Merced has a
    hard time calculating this part of the tranfer matrix that this function is designed to handle. This function
    only treats the photon relativistically, using Newtonian kinematics for the other particles. Note, as of the
    writing of this function, Merced treats all particles relativistically when one of the particles has zero mass.

    This functions ignores any flux dependency on energy and only calculates the Legendre order 0 term, all other terms a returned a 0.0's.

    :param masses:                          A python dictionary containing the masses of the two incoming and two outgoing particles.
    :param Q:                               This is the Q-value for the two-body reaction.
    :param modifiedProductGroupBoundaries:  A smaller set of the actual multi-group boundaries.
    :param crossSection:                    This is the pointwise linear cross section for the reaction.
    :param angularData:                     This is the angular distribution data for the product.
    :param legendreMax:                     This is the maximum Legendre order to calculate the product matrices.
    :param groupFlux:                       This is the multi-grouped flux for the group.
    :param firstCall:                       A boolean that indicates if the function is calling itself.

    :returns:                               The number and energy conserving product matices.
    """

    def K_boost(photonEnergy, m2, m3):
        """
        Returns the kinetic energy the product with its speed equal to the boost speed.

        :param photonEnergy:    The energy of the incident photon.
        :param m2:              The mass of the target.
        :param m3:              The mass of the product.

        :returns:               A python float.
        """

        v_boost = photonEnergy / ( m2 + photonEnergy )
        return 0.5 * m3 * v_boost**2

    def K_com(photonEnergy, m2, m3, m4, Q):
        """
        Returns the kinetic energy of the product in the center-of-mass frame.

        :param photonEnergy:    The energy of the incident photon.
        :param m2:              The mass of the target.
        :param m3:              The mass of the product.
        :param m4:              The mass of the residual.
        :param Q:               The Q-value for the reaction.

        :returns:               A python float.
        """

        return ( photonEnergy * math.sqrt( m2 / ( m2 + 2 * photonEnergy ) ) + K_boost(photonEnergy, m2, m2) + Q ) * m4 / ( m3 + m4 )

    def findEgForK3(photonEnergy, m2, m3, m4, K3, sign):
        """
        Returns the incident photon energy where the minimum outgoing particle kinetic energy is K3.

        :param photonEnergy:    The energy of the incident photon.
        :param m2:              The mass of the target.
        :param m3:              The mass of the product.
        :param m4:              The mass of the residual.
        :param K3:              The kinetic energy of the product.
        :param sign:            -1 if lower photon energy is to be returned and 1 if upper photon energy is to be returned.

        :returns:               A python float.
        """

        Eg = photonEnergy
        for i in range(20):
            EgPrior = Eg
            Qp = -Q
            if( K3 != 0.0 ): Qp += sign * math.sqrt( K3 * ( 2 * ( K_com(Eg, m2, m3, m4, Q) + K_boost(Eg, m2, m3) ) - K3 ) ) * ( m3 + m4 ) / m4
            Eg = Qp / ( math.sqrt( m2 / (m2 + 2 * Eg) ) + 0.5 * Eg / ( m2 + Eg )**2 * ( ( m2 - m3 ) * m4 - m3 * m3 ) / m4 )
            if abs(Eg - EgPrior) < 4e-16 * Eg: break

        return Eg

    if firstCall:
        firstBoundaries = [ 0.0, modifiedProductGroupBoundaries[0] ]
        TM1_firstCall, TME_firstCall = twoBodyPhotoNuclearAtThreshold(masses, Q, firstBoundaries, crossSection, 
                angularData, legendreMax, groupFlux, firstCall = False)

    projectileMass = masses['Projectile']
    targetMass = masses['Target']
    productMass = masses['Product']
    residualMass = targetMass - productMass - Q

    TM1 = []
    TME = []

    EgMin = findEgForK3(-Q, targetMass, productMass, residualMass, 0, 1)
    crossSectionAtMin = crossSection.evaluate(EgMin)
    if crossSectionAtMin is None: crossSectionAtMin = 0.0
    if crossSectionAtMin == 0.0:
        for K_3_upper in modifiedProductGroupBoundaries[:-1]:
            LegendreCoefficients = ( legendreMax + 1 ) * [ 0.0 ]
            TM1.append(LegendreCoefficients)
            LegendreCoefficients = ( legendreMax + 1 ) * [ 0.0 ]
            TME.append(LegendreCoefficients)
        return TM1, TME
        
    PofMuAtMin = angularData.evaluate(EgMin)

    K_3_lower = None
    numberOfEgs = 10
    for K_3_upper in modifiedProductGroupBoundaries:
        if K_3_lower is not None:
            EgLowerLeft  = findEgForK3(EgMin, targetMass, productMass, residualMass, K_3_lower, -1)
            EgLowerRight = findEgForK3(EgMin, targetMass, productMass, residualMass, K_3_lower,  1)

            dmu = 0.0
            for index in range(numberOfEgs + 1):
                Eg = ( index * EgLowerLeft + ( numberOfEgs - index ) * EgLowerRight ) / numberOfEgs
                sqrtTerm = math.sqrt(K_boost(Eg, targetMass, productMass) * K_com(Eg, targetMass, productMass, residualMass, Q))
                PofMu = 0.5                                          # 0.5 is P(mu) which should be gotten from PofMuAtMin.
                dmu += 0.5 * ( K_3_upper - K_3_lower ) / sqrtTerm * PofMu

            EgUpperLeft  = findEgForK3(EgMin, targetMass, productMass, residualMass, K_3_upper, -1)
            EgUpperRight = findEgForK3(EgMin, targetMass, productMass, residualMass, K_3_upper,  1)
            dEg = 0.5 * ( ( EgUpperRight - EgUpperLeft ) + ( EgLowerRight - EgLowerLeft ) )
            area = dEg * dmu / ( numberOfEgs + 1 )

            LegendreCoefficients = ( legendreMax + 1 ) * [ 0.0 ]
            LegendreCoefficients[0] = area * crossSectionAtMin / groupFlux
            if firstCall:
                for index, value in enumerate(TM1_firstCall[0]): LegendreCoefficients[index] += value
            TM1.append(LegendreCoefficients)

            LegendreCoefficients = ( legendreMax + 1 ) * [ 0.0 ]
            LegendreCoefficients[0] = area * crossSectionAtMin * EgMin / groupFlux
            if firstCall:
                for index, value in enumerate(TME_firstCall[0]): LegendreCoefficients[index] += value
                firstCall = False
            TME.append(LegendreCoefficients)

        K_3_lower = K_3_upper

    return TM1, TME
