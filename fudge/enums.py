# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from LUPY import enums as enumsModule
from fudge import GNDS_formatVersion as GNDS_formatVersionModule

class Interaction(enumsModule.Enum):
    '''Class that defines the all interaction enums.'''

    nuclear = enumsModule.auto()
    atomic = enumsModule.auto()
    TNSL = 'thermalNeutronScatteringLaw'
    LLNL_TNSL = 'LLNL_TNSL'
    legacyTNSL = 'TNSL'

    @staticmethod
    def getTNSL_interaction(formatVersion):
        '''
        Returns the interation string based on the GNDS version format.

        :param formatVersion:   The GNDS version format.

        :return:                The interation string based on the GNDS version format *formatVersion*.
        '''

        majorVersion, minorVersion, betaVersion = GNDS_formatVersionModule.versionComponents(formatVersion)
        if majorVersion < 2:
            return Interaction.legacyTNSL

        return Interaction.TNSL

    @staticmethod
    def allowedStrings():

        return tuple(str(enum) for enum in Interaction)

class Genre(enumsModule.Enum):

    twoBody = enumsModule.auto()
    NBody = enumsModule.auto()
    production = enumsModule.auto()
    sumOfRemainingOutputChannels = enumsModule.auto()

class FissionGenre(enumsModule.Enum):
    '''Enum class representing the allowed fission genres.'''

    none = enumsModule.auto()
    total = enumsModule.auto()
    firstChance = enumsModule.auto()
    secondChance = enumsModule.auto()
    thirdChance = enumsModule.auto()
    fourthChance = enumsModule.auto()

class Conserve(enumsModule.Enum):

    number = enumsModule.auto()
    energyOut = enumsModule.auto()
