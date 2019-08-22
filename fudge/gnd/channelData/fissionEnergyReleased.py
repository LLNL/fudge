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
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# <<END-copyright>>

from fudge.gnd import baseClasses
import base
from fudge.legacy.converting import endfFormats
from fudge.gnd import tokens

__metaclass__ = type

fissionEnergyReleasedForms = [ tokens.polynomialFormToken ]

#
# fissionEnergyReleased genre and forms
#
class component( baseClasses.componentBase ) :

    genre = base.fissionEnergyReleasedToken

    def __init__( self ) :

        baseClasses.componentBase.__init__( self, fissionEnergyReleasedForms )

    def check( self, info ) :

        from fudge.gnd import warning
        warnings = []
        for form in self.forms:
            fwarnings = self.forms[form].check( info )
            if fwarnings:
                warnings.append( warning.context('%s:' % form, fwarnings) )
        return warnings

    @staticmethod
    def parseXMLNode( FERelement, xPath=[], linkData={} ):
        """ parse <fissionEnergyReleased> from xml """

        xPath.append( FERelement.tag )
        fer = component()
        fer.nativeData = FERelement.get('nativeData')
        for form in FERelement:
            if form.tag==tokens.polynomialFormToken:
                fer.addForm( polynomial.parseXMLNode( form, xPath, linkData ) )
        xPath.pop()
        return fer

    def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

        self.forms[self.nativeData].toENDF6( MT, endfMFList, flags, targetInfo )

class polynomial( baseClasses.formBase ) :

    genre = component.genre
    form = tokens.polynomialFormToken

    labels = [ 'promptProductKE', 'promptNeutronKE', 'delayedNeutronKE', 'promptGammaEnergy', 'delayedGammaEnergy', 'delayedBetaEnergy', 'neutrinoEnergy',
                'nonNeutrinoEnergy', 'totalEnergy' ]

    def __init__( self, order, data, energyUnit, hasUncertainties = False ) :

        baseClasses.formBase.__init__( self )
        self.order = order
        self.data = data
        self.energyUnit = energyUnit                        # Need to use this in toENDF6.
        self.hasUncertainties = hasUncertainties

    def check( self, info ):

        from fudge.gnd import warning
        from pqu.physicalQuantityWithUncertainty import PhysicalQuantityWithUncertainty as PQU
        warnings = []
        domain = [d.getValueAs( self.energyUnit ) for d in info['crossSectionDomain']]
        for key in self.data:   # 'totalEnergy', 'nonNeutrinoEnergy', etc
            coefs = [a[0] for a in self.data[key]]
            poly = lambda e: sum( [coef * e**i for i,coef in enumerate(coefs)] )
            for ein in domain:
                if not 0 < poly(ein) < 4e+8: # 0 MeV -> 400 MeV considered acceptable
                    warnings.append( warning.badFissionEnergyRelease( key, PQU(ein,self.energyUnit), 
                        PQU(poly(ein),self.energyUnit) ) )
        return warnings

    @staticmethod
    def parseXMLNode( element, xPath=[], linkData={} ):
        """ translate <polynomial> element from xml """

        xPath.append( element.tag )
        order = int( element.get("order") )
        hasUncertainties = False
        dataDict = {}
        for data in element:
            coefs = map(float, data.text.split())
            if len(coefs) == (order+1)*2:
                hasUncertainties=True
                coefs = zip( coefs[::2], coefs[1::2] )
            assert len(coefs) == order+1
            dataDict[data.tag] = coefs
        Poly = polynomial( order, dataDict, element.get('energyUnit'), hasUncertainties )
        xPath.pop()
        return Poly

    def toXMLList( self, indent = "" ) :

        hasUncertainties, indent2 = '', indent + '  '
        if( self.hasUncertainties ) : hasUncertainties = ' hasUncertainties="true"'
        xmlString = [ '%s<%s order="%s" energyUnit="%s"%s>' % ( indent, self.form, self.order, self.energyUnit, hasUncertainties ) ]
        for label in self.labels :
            data = ""
            for d in self.data[label] :
                for v in d : data += ' %s' % v
            xmlString.append( '%s<%s>%s</%s>' % ( indent2, label, data, label ) )
        xmlString[-1] += '</%s>' % self.form
        return( xmlString )

    def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

        n, data = self.order + 1, {}
        for i in xrange( n ) :
            data[i] = []
            for index, label in enumerate( self.labels ) :
                d = self.data[label][i][0]
                u = 0
                if( self.hasUncertainties ) : u = self.data[label][i][1]
                data[i].append( d )
                data[i].append( u )
        list = []
        for i in xrange( n ) : list += data[i]

        endfMFList[MT][458] = [ endfFormats.endfContLine( targetInfo['ZA'], targetInfo['mass'], 0, 0, 0, 0 ) ]
        endfMFList[MT][458].append( endfFormats.endfContLine( 0, 0, 0, self.order, 18 * n, 9 * n ) )
        endfMFList[MT][458] += endfFormats.endfDataList( list )
        endfMFList[MT][458].append( endfFormats.endfSENDLineNumber( ) )
