# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

'''
Functions to determine if the interaction attribute for the reactionSuite, and map Protare and Import classes 
is missing; and if so, return a good guess of what its value should be.
'''

from PoPs import IDs as PoPsIDsModule
from PoPs.chemicalElements import misc as PoPsGroupMiscModule
from fudge import enums as enumsModule

names = ( 'HinCH2',                 'HinH2O',               'HinC5O2H8',                'HinIceIh',     'HinYH2',
          'HinZrH',                 'ortho-H',              'para-H',                   'DinD2O',       'ortho-D',
          'para-D',                 'Be-metal',             'BeinBeO',                  'CinSiC',       'graphite',
          'crystalline-graphite',   'reactor-graphite-10P', 'reactor-graphite-30P',     'NinUN',        'OinBeO',
          'OinD2O',                 'OinIceIh',             'OinUO2',                   'tnsl-Al27',    'SiinSiC',
          'tnsl-Fe56',              'YinYH2',               'ZrinZrH',                  'UinUN',        'UinUO2',
          'SiO2',                   'SiO2-alpha',           'SiO2-beta',                's-CH4',        'l-CH4',
          'benzine',                'benzene' )

def isTNSL_name(name):
    '''Returns True if the specified name is a FUDGE TNSL name and False otherwise.'''

    return name in names

def guessInteraction(interaction, projectileID, targetID):
    '''If the interaction argument is None, the interaction is guessed by looking at the targetID. If interaction
    is not None, (False, interaction) are returned.  If interaction is None, (True, interaction) are returned
    where the following logic is used to guess (i.e., set) interaction:

        -) if passing targetID to the function *isTNSL_name* returns True, then the interaction is set to TNSL.
        -) if targetID is in PoPsGroupMiscModule.ZFromSymbol, then the interaction is set to atomic.
        -) otherwise, the interaction is set to nuclear.

    For ENDL 99120 and 99125 ZAs, the projectileID is also needed to determine if the interaction is nuclear or photo-atomic .
    '''

    guessedInteraction = False
    if interaction is None:
        guessedInteraction = True
        if 'ENDL9912' in targetID:
            interaction = enumsModule.Interaction.nuclear if projectileID != PoPsIDsModule.photon else enumsModule.Interaction.atomic
        elif isTNSL_name(targetID):
            interaction = enumsModule.Interaction.TNSL
        elif targetID in PoPsGroupMiscModule.ZFromSymbol:
            interaction = enumsModule.Interaction.atomic
        else:
            interaction = enumsModule.Interaction.nuclear
    else:
        interaction = enumsModule.Interaction.checkEnumOrString(interaction)

    return guessedInteraction, interaction
