# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
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
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
# <<END-copyright>>

from PoPs.groups import misc as chemicalElementMiscPoPsModule
from PoPs import IDs as IDsPoPsModule

from fudge.gnds import tokens as tokensModule
from fudge.gnds import product as productModule
from fudge.gnds.productData.distributions import unspecified as unspecifiedModule

from xData import standards as standardsModule

from . import gndsToENDF6 as gndsToENDF6Module
from . import endfFormats as endfFormatsModule

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

    if( flags['verbosity'] >= 10 ) : print( '%s%s: label = %s: to ENDF6:' % ( verbosityIndent, self.id, self.label ) )
    priorMF6flag = targetInfo['doMF4AsMF6']
    conversionFlag = targetInfo['ENDFconversionFlags'].get(self,"")
    if( conversionFlag == 'MF6' ) :
        targetInfo['doMF4AsMF6'] = True # flag was set in reaction.py, but may need to be overwritten

    gammasPresent = False
    if( self.outputChannel is not None ) :
        for product in self.outputChannel : gammasPresent = gammasPresent or ( product.id == IDsPoPsModule.photon )
    if( len( self.distribution ) or gammasPresent ) :        # First part should now always be true.
        targetInfo['zapID'] = self.id
        particle = targetInfo['reactionSuite'].PoPs[self.id]
        ZA = chemicalElementMiscPoPsModule.ZA( particle )
        try :
            targetInfo['particleMass'] = targetInfo['massTracker'].getMassAWR( ZA, asTarget = False )
        except :
            pass
        targetInfo['multiplicity'] = self.multiplicity
        if conversionFlag != 'implicitProduct':
            self.distribution.toENDF6( MT, endfMFList, flags, targetInfo )
        if MT in (527, 528) and targetInfo['style'] in self.energyDeposition:
            energyLoss(self, MT, endfMFList, flags, targetInfo)
    if( self.outputChannel is not None ) :
        priorIndex, priorToken, priorLabel = targetInfo['productIndex'], targetInfo['productToken'], targetInfo['productLabel']
        for index, product in enumerate( self.outputChannel ) :
            if( product.id == IDsPoPsModule.photon ) :
                targetInfo['gammas'].append( product )
                continue
            targetInfo['productIndex'] = "%s.%s" % ( priorIndex, index )
            targetInfo['productToken'] = product.id
            targetInfo['productLabel'] = product.label
            product.toENDF6( MT, endfMFList, flags, targetInfo, verbosityIndent = verbosityIndent + '    ' )
        targetInfo['productIndex'], targetInfo['productToken'], targetInfo['productLabel'] = priorIndex, priorToken, priorLabel
    targetInfo['doMF4AsMF6'] = priorMF6flag

productModule.product.toENDF6 = toENDF6

def energyLoss( self, MT, endfMFList, flags, targetInfo ) :

    energyLoss = self.energyDeposition[ targetInfo['style'] ]
    data = [ [ energyLoss.domainMin, energyLoss.domainMin ], [ energyLoss.domainMax, energyLoss.domainMax ] ]
    energyLoss = energyLoss.__class__( data = data, axes = energyLoss.axes ) - energyLoss
    data = []
    for xy in energyLoss.copyDataToXYs( ) : data += xy
    NE = len( energyLoss )
    EInInterpolation = gndsToENDF6Module.gndsToENDFInterpolationFlag( energyLoss.interpolation )
    ENDFDataList = [ endfFormatsModule.endfContLine( 0, 0, 0, 0, 1, NE ) ] + \
            endfFormatsModule.endfInterpolationList( [ NE, EInInterpolation ] )
    ENDFDataList += endfFormatsModule.endfDataList( data )
    frame = standardsModule.frames.labToken
    gndsToENDF6Module.toENDF6_MF6( MT, endfMFList, flags, targetInfo, 8, frame, ENDFDataList )
