# <<BEGIN-copyright>>
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
    def parseXMLNode( FERelement, linkData={} ):
        """ parse <fissionEnergyReleased> from xml """
        fer = component()
        fer.nativeData = FERelement.get('nativeData')
        for form in FERelement:
            if form.tag==tokens.polynomialFormToken:
                fer.addForm( polynomial.parseXMLNode( form, linkData ) )
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
    def parseXMLNode( polynomialElement, linkData={} ):
        """ translate <polynomial> element from xml """
        order = int( polynomialElement.get("order") )
        hasUncertainties = False
        dataDict = {}
        for data in polynomialElement:
            coefs = map(float, data.text.split())
            if len(coefs) == (order+1)*2:
                hasUncertainties=True
                coefs = zip( coefs[::2], coefs[1::2] )
            assert len(coefs) == order+1
            dataDict[data.tag] = coefs
        return polynomial( order, dataDict, polynomialElement.get('energyUnit'), hasUncertainties )

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
