# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module adds the method toACE to the reactionSuite class.
"""

from fudge import enums as enumsModule
from fudge import styles as stylesModule
from fudge import reactionSuite as reactionSuiteModule

from . import gndsToACE as gndsToACEModule
from . import GNDS_toACE_TNSL as GNDS_toACE_TNSL_Module

def toACE(self, args, styleName, fileName, evaluationId, addAnnotation=False, verbose=0, skipURR=False, skipILF_logic=False):
    """
    Produce an ACE file with data from self.

    :param styleName: The name of the griddedCrossSection style whose data are to be translated to ACE.
    :param str fileName: path to save resulting ACE file.
    :param int evaluationId:  evaluation identifier, 2 digits max.
        That is the 'nn' in the HZ string 'ZZZAAA.nnC'
    :type temperature: PQU or string
    :param bool addAnnotation: if True, adds comments to help navigate the resulting ACE file
    """

    if( verbose > 0 ) : print( self.inputParticlesToReactionString( ) )

    cdf_style = self.styles[styleName].findDerivedFromStyle( stylesModule.MonteCarlo_cdf )

    if self.interaction == enumsModule.Interaction.TNSL:
        GNDS_toACE_TNSL_Module.toACE(self, args, styleName, cdf_style, fileName, evaluationId, addAnnotation)
    else:

        ACE_data = []
        delayedNeutronRateAndDatas = []
        for reaction in self.reactions : reaction.toACE( styleName, cdf_style, ACE_data, delayedNeutronRateAndDatas, verbose )
        gndsToACEModule.toACE(self, styleName, cdf_style, fileName, evaluationId, ACE_data, delayedNeutronRateAndDatas, 
                addAnnotation=addAnnotation, skipURR=skipURR, skipILF_logic=skipILF_logic)

reactionSuiteModule.ReactionSuite.toACE = toACE
