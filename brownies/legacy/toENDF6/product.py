# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from PoPs.groups import misc as chemicalElementMiscPoPsModule
from PoPs import IDs as IDsPoPsModule

from fudge import product as productModule
from fudge.productData import multiplicity as multiplicityModule
from fudge.productData.distributions import unspecified as unspecifiedModule, angular as angularModule

from xData import standards as standardsModule

from . import gndsToENDF6 as gndsToENDF6Module
from . import endfFormats as endfFormatsModule

def toENDF6( self, MT, endfMFList, flags, targetInfo, verbosityIndent = '' ) :

    def getPromptOrTotalNubar( self ) :

        return( gndsToENDF6Module.getForm( targetInfo['style'], self.multiplicity ) )

    targetInfo['product'] = self
    targetInfo['delayedNubarWeight'] = None
    if( MT == 455 ) :
        if( MT not in endfMFList[5] ) : endfMFList[5][MT] = []
        targetInfo['delayedNubarWeight'] = self.ENDF6_delayedNubarWeights
    elif( targetInfo['isFission'] ) :
        if( targetInfo['totalFission'] ) :
            totalNubar = getPromptOrTotalNubar( self )
            if( not( isinstance( totalNubar, ( multiplicityModule.unspecified, multiplicityModule.reference ) ) ) ) : targetInfo['totalNubar'] = getPromptOrTotalNubar( self )
        else :
            targetInfo['promptNubar'] = getPromptOrTotalNubar( self )

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
        else:
            form = gndsToENDF6Module.getForm(targetInfo['style'], self.distribution)
            if isinstance(form, unspecifiedModule.form):
                pass
            elif isinstance(form, angularModule.twoBodyForm) and isinstance(form.angularSubform, angularModule.recoil):
                pass
            else:
                print("WARNING: product %s in reaction %s is listed as implicit but has a real distribution" %
                      (self.pid, targetInfo['reaction']))
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

    energyLoss = gndsToENDF6Module.getForm( targetInfo['style'], self.energyDeposition )
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
