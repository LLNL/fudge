# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains miscellaneous functions used by the modules that convert a GNDS reactionSuite into an ACE file.
"""


floatFormatOthers = '%20.12E'
floatFormatEnergyGrid = '%21.13E'
floatFormatEnergyGrid = floatFormatOthers

def updateXSSInfo(label, annotates, XSS, data):

    annotates.append((len(XSS), label))
    XSS += data

def intArrayToRecords(intArray):

    stringData = ['%9d' % d for d in intArray]
    lines = []
    while(len(stringData)) :
        lines.append(''.join(stringData[:8]))
        stringData = stringData[8:]
    return lines

def XSSToStrings(annotates, XSS, addAnnotation, numberOfEnergiesInGrid):

    strData, record, i3 = [], [], 0
    annotates.append((len(XSS), ''))
    i2, label = annotates[i3]
    for i1, datum in enumerate(XSS):
        if( type( datum ) == int ) :
            record.append( "%20d" % datum )
        else :
            try :
                floatFormat = floatFormatOthers
                if( i1 < numberOfEnergiesInGrid ) : floatFormat = floatFormatEnergyGrid
                record.append( floatFormat % datum )
            except :
                print( i1, len(XSS), datum, annotates[i3] )
                raise
        if( len( record ) == 4 ) :
            annotate = ''
            if( i2 <= i1 ) :
                annotate = ' !'
                sep = ' '
                while( i2 <= i1 ) :
                    annotate += sep + label + ' (%s)' % ( i2 - i1 + 3 )
                    sep = ', '
                    i3 += 1
                    i2, label = annotates[i3]
            if( not addAnnotation ) :   annotate = ''
            strData.append( ''.join( record ) + annotate )
            record = []

    if len(record) > 0: strData.append(''.join(record))

    return strData
