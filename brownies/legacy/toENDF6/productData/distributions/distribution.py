# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from PoPs.chemicalElements import misc as chemicalElementMiscPoPsModule

from fudge import product as productModule
from fudge.productData.distributions import distribution as distributionModule

from ... import endfFormats as endfFormatsModule
from ... import gndsToENDF6 as gndsToENDF6Module

def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

    if( len( self ) == 0 ) :
        raise Exception( 'This needs to be tested' )     # In particular, the method endfMultiplicityList has been removed.
        if( MT in [ 452, 455, 456 ] ) : return
        targetInfo['MF6LCTs'].append( None )
        particle = targetInfo['reactionSuite'].PoPs[targetInfo['zapID']]
        ZAP = chemicalElementMiscPoPsModule.ZA( particle )
        nPoints, multiplicityList = targetInfo['multiplicity'].endfMultiplicityList( targetInfo )
        AWP = targetInfo['massTracker'].getMassAWR(ZAP, asTarget = False)
        endfMFList[6][MT] += [ endfFormatsModule.endfContLine( ZAP, AWP, 0, 0, 1, nPoints ) ]
        endfMFList[6][MT] += multiplicityList
        return
    form = gndsToENDF6Module.getForm( targetInfo['style'], self )
    if( hasattr( form, 'toENDF6' ) ) :
        try:
            form.toENDF6( MT, endfMFList, flags, targetInfo )
        except Exception as ex:
            print("\nError encountered while writing distribution %s to ENDF-6:" % self.toXLink())
            raise
    else :
        product = self.findClassInAncestry(productModule.Product)
        if form is None:
            print( 'WARNING: found no distribution of style "%s" for product "%s"' % (targetInfo['style'], product.toXLink()))
        else:
            print( 'WARNING: unsupported distribution (type %s) for product "%s"' % (form.__class__, product.toXLink()) )

distributionModule.Component.toENDF6 = toENDF6
