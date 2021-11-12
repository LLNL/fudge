# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from fudge.productData.distributions import Legendre as LLNLLegendreModule

from ... import gndsToENDF6 as gndsToENDF6Module
from ... import endfFormats as endfFormatsModule

#
# form
#
def toENDF6( self, MT, endfMFList, flags, targetInfo ) :    # FIXME appears to be broken code

    LAW, LEP = 1, 2
    EInInterpolation = gndsToENDF6Module.axisToEndfInterpolationFlag( self.axes[1] )
    E_ins = [ [ EpCl.outerDomainValue, {} ] for EpCl in self[0].EpP ]
    for l_EEpCl in self :
        for indexE, EpCl in enumerate( l_EEpCl ) :
            if( EpCl.outerDomainValue != E_ins[indexE][0] ) :
                raise Exception( "E_in = %s not in list E_ins" % EpCl.outerDomainValue )
            for Ep, Cl in EpCl :
                if( Ep not in E_ins[indexE][1] ) :
                    E_ins[indexE][1][Ep] = []
                    for l in xrange( l_EEpCl.l ) : E_ins[indexE][1][Ep].append( 0. )
                E_ins[indexE][1][Ep].append( Cl )
    ENDFDataList = [ endfFormatsModule.endfContLine( 0, 0, 1, LEP, 1, len( E_ins ) ) ]
    ENDFDataList += endfFormatsModule.endfInterpolationList( [ len( E_ins ), EInInterpolation ] )
    for Es in E_ins :
        NA, data = 0, []
        for key in sorted( Es[1] ) :
            LegendreSeries = Es[1][key]
            NA = max( len( LegendreSeries ) , NA )
            data += [ key ] + LegendreSeries
        ENDFDataList.append( endfFormatsModule.endfContLine( 0, Es[0], 0, NA - 1, len( data ), len( data ) / ( NA + 1 ) ) )
        ENDFDataList += endfFormatsModule.endfDataList( data )
    gndsToENDF6Module.toENDF6_MF6( MT, endfMFList, flags, targetInfo, LAW, self.productFrame, ENDFDataList )

LLNLLegendreModule.form.toENDF6 = toENDF6
