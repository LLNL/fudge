# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module adds the URR probability tables if they exists.
"""

from xData import XYs as XYs1dModule
from fudge import styles as stylesModule

def toACE_tables( reactionCrossSection, GNDS_URR_data, cumulatives ) :

    ACE_URR_data = {}
    data = GNDS_URR_data.data
    for function in data.functionals :
        pdf = XYs1dModule.XYs1d( [ function.xs.values, function.pdf.values ], dataForm = 'xsandys' )
        cdf = XYs1dModule.XYs1d( [ function.xs.values, function.cdf.values ], dataForm = 'xsandys' )
        inv_cdf = cdf.inverse( )

        ACE_URR_crossSections = []
        xs = [ inv_cdf.evaluate( cumulative ) for cumulative in cumulatives ] + [ inv_cdf[-1][1] ]
        for i1, x2 in enumerate( xs )  :
            if( i1 > 0 ) :
                ACE_URR_crossSections.append( float( pdf.integrateWithWeight_x( domainMin = x1, domainMax = x2 ) / pdf.integrate( domainMin = x1, domainMax = x2 ) ) )
            x1 = x2
        ACE_URR_data[function.outerDomainValue] = [ reactionCrossSection.evaluate( function.outerDomainValue ), ACE_URR_crossSections ]

    return( ACE_URR_data )

def URR_probabilityTable( self, numberOfProbabilities = 20 ) :

    URRPT = []
    URR_styles = self.styles.getStylesOfClass( stylesModule.URR_probabilityTables )
    if( len( URR_styles ) > 0 ) :
        if( len( URR_styles ) > 1 ) :
            print( '    WARNING: more than 1 URR style found. This is currently not supported hence no URR probability tables is being skipped.' )
        else :
            URR_style = URR_styles[0]
            heatedStyle = URR_style
            while( True ) :
                heatedStyle = self.styles[heatedStyle.derivedFrom]
                if( heatedStyle is None ) : raise Exception( 'Could not find heated style for URR data.' )
                if( isinstance( heatedStyle, stylesModule.heated ) ) : break

            cumulatives = [ i1 / float( numberOfProbabilities ) for i1 in range( numberOfProbabilities ) ]
            URR_label = URR_style.label
            MTs_URR = {}
            for reaction in self.reactions :
                if( URR_label in reaction.crossSection ) :
                    reactionCrossSection = reaction.crossSection[heatedStyle.label]
                    ACE_I_index = { 2 : 3, 18 : 4, 102 : 5 }[reaction.ENDF_MT]
                    MTs_URR[ACE_I_index] = toACE_tables( reactionCrossSection, reaction.crossSection[URR_label], cumulatives )

            energies = sorted( MTs_URR[3].keys( ) )
            URRPT = [ len( energies ), numberOfProbabilities, 2, 0, 0, 1 ] + energies

            zeros = numberOfProbabilities * [ 0.0 ]
            cumulatives = cumulatives + [ 1.0 ]
            probabilities = cumulatives[1:]
            for energy in energies :
                elastic = MTs_URR[3][energy]
                fission = [ 0.0, zeros ]
                if( 4 in MTs_URR ) : fission = MTs_URR[4][energy]
                capture = MTs_URR[5][energy]
                total = []
                totalCrossSection = elastic[0] + fission[0] + capture[0]
                for i1 in range( len( elastic[1] ) ) :
                    total.append( ( elastic[0] * elastic[1][i1] + fission[0] * fission[1][i1] + capture[0] * capture[1][i1] ) / totalCrossSection )

                URRPT += probabilities
                URRPT += total
                URRPT += elastic[1]
                URRPT += fission[1]
                URRPT += capture[1]
                URRPT += zeros

    return( URRPT )
