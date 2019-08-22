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
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# <<END-copyright>>

import fudge.gnd.tokens as tokensModule
import fudge.gnd.product as productModule
from fudge.gnd.productData.distributions import unspecified as unspecifiedModule 

def toENDF6( self, MT, endfMFList, flags, targetInfo, verbosityIndent = '' ) :

    def getPromptOrTotalNubar( self ) :

        return( self.multiplicity[targetInfo['style']] )

    targetInfo['product'] = self
    targetInfo['delayedNubarWeight'] = None
    if( 'emissionMode' in self.attributes ) :
        if( self.getAttribute( 'emissionMode' ) == tokensModule.delayedToken ) :
            MT = 455
            if( MT not in endfMFList[5] ) : endfMFList[5][MT] = [ ]
            targetInfo['delayedNubarWeight'] = self.ENDF6_delayedNubarWeights
        elif( self.getAttribute( 'emissionMode' ) == tokensModule.promptToken ) :
            targetInfo['promptNubar'] = getPromptOrTotalNubar( self )
        elif( self.getAttribute( 'emissionMode' ) == 'total' ) :
                targetInfo['totalNubar'] = getPromptOrTotalNubar( self )

    if( flags['verbosity'] >= 10 ) : print '%s%s: label = %s: to ENDF6:' % ( verbosityIndent, self.name, self.label )
    priorMF6flag = targetInfo['doMF4AsMF6']
    if( self.attributes.get( 'ENDFconversionFlag' ) == 'MF6' ) :
        targetInfo['doMF4AsMF6'] = True # flag was set in reaction.py, but may need to be overwritten
    if( len( self.distribution ) ) :        # This should now always be true.
        distribution = self.distribution[targetInfo['style']]
        if( not( isinstance( distribution, unspecifiedModule.form ) ) ) :
            targetInfo['zapID'] = self.particle.name
            if( hasattr( self.particle,'groundState' ) ) :  # ENDF wants ground state mass:
                targetInfo['particleMass'] = self.particle.groundState.getMass( 'eV/c**2' )
            else:
                targetInfo['particleMass'] = self.getMass( 'eV/c**2' )
            targetInfo['multiplicity'] = self.multiplicity
            self.distribution.toENDF6( MT, endfMFList, flags, targetInfo )
    if( not( self.decayChannel is None ) ) :
        priorIndex, priorToken, priorLabel = targetInfo['productIndex'], targetInfo['productToken'], targetInfo['productLabel']
        for index, product in enumerate( self.decayChannel ) :
            if( product.name == 'gamma' ) :
                targetInfo['gammas'].append( product )
                continue
            targetInfo['productIndex'] = "%s.%s" % ( priorIndex, index )
            targetInfo['productToken'] = product.name
            targetInfo['productLabel'] = product.label
            product.toENDF6( MT, endfMFList, flags, targetInfo, verbosityIndent = verbosityIndent + '    ' )
        targetInfo['productIndex'], targetInfo['productToken'], targetInfo['productLabel'] = priorIndex, priorToken, priorLabel
    targetInfo['doMF4AsMF6'] = priorMF6flag

productModule.product.toENDF6 = toENDF6
