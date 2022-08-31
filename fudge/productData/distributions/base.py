# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""Base classes for distributions."""

from LUPY import ancestry as ancestryModule
from xData import enums as xDataEnumsModule

from PoPs.chemicalElements import misc as chemicalElementMiscPoPsModule

from fudge import abstractClasses as abstractClassesModule

class Form( abstractClassesModule.Form ) :

    keyName = 'label'

    def __init__( self, label, productFrame, subForms ) :

        abstractClassesModule.Form.__init__( self )

        self.label = label
#        self.productFrame = productFrame                # This does not seem to work, hence the next line.
        self.__productFrame = xDataEnumsModule.Frame.checkEnumOrString(productFrame)

        for i1, subform_ in enumerate( subForms ) :
            setattr( self, self.subformAttributes[i1], subform_ )
            if( subform_ is not None ) : subform_.setAncestor( self )
        self.subforms = [ subform_ for subform_ in subForms ]

    @property
    def label( self ) :

        return( self.__label )

    @label.setter
    def label( self, value ) :

        if( value is not None ) :
            if( not( isinstance( value, str ) ) ) : raise TypeError( 'value must be a string' )
        self.__label = value

    @property
    def product( self ) :

        from fudge import product as productModule

        return self.findClassInAncestry(productModule.Product)

    @property
    def productFrame( self ) :

        return( self.__productFrame )

    @productFrame.setter
    def productFrame(self, frame):
        """Set the distribution's product frame."""

        self.__productFrame = xDataEnumsModule.Frame.checkEnumOrString(frame)

    def convertUnits( self, unitMap ) :
        """See documentation for reactionSuite.convertUnits."""

        if( hasattr( self, 'subforms' ) ) :
            for subform_ in self.subforms :
                if( subform_ is None ) : continue                    # FIXME, needed for photo-atomic coherentScattering data with no anomalousScattering subforms.
                subform_.convertUnits( unitMap )

    def findEntity( self, entityName, attribute = None, value = None ):
        """
        Find specific subform within the form. Overrides ancestry.findEntity.
        """

        for subform_ in self.subforms :
            if( subform_.moniker == entityName ) : return( subform_ )
        return( ancestryModule.Ancestry.findEntity( self, entityName, attribute, value ) )

    def isTwoBody( self ) :

        return( False )

    def isSpecified( self ) :

        return( True )

    def residualMass( self, PoPs, residualZA, massUnit, compoundMass, product = None, verbose = 1 ) :

        residual = None
        residualSymbol = chemicalElementMiscPoPsModule.symbolFromZ[residualZA//1000]
        residualID = chemicalElementMiscPoPsModule.isotopeSymbolFromChemicalElementIDAndA( residualSymbol, residualZA % 1000 )

        try :
            residual = PoPs[residualID]
            return( residual.getMass( massUnit ) )
        except :
            pass

        if( product is None ) : product = self.product
        _residualMass = compoundMass - product.getMass( massUnit )
        if( verbose > 0 ) : print('    Could not find residual in particle database: id = "%s", ZA = %d, using mass %s %s' % ( residualID, residualZA, _residualMass, massUnit ))
        return( _residualMass )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        return( self.subforms[0].toPointwise_withLinearXYs( **kwargs ) )

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attributeStr = ''
        if self.label is not None:
            attributeStr += ' label="%s"' % self.label
        if self.productFrame != xDataEnumsModule.Frame.none:
            attributeStr += ' productFrame="%s"' % self.productFrame

        xmlString = [ '%s<%s%s>' % ( indent, self.moniker, attributeStr ) ]
        for subform_ in self.subforms : xmlString += subform_.toXML_strList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker

        return( xmlString )

class Subform( ancestryModule.AncestryIO ) :

    def __init__( self ) :

        ancestryModule.AncestryIO.__init__( self )
        self.label = None

    def XMLStartTagString( self, indent = '', extraAttributesAsStrings = '', emptyTag = False ) :

        emptyTagStr = ''
        if( emptyTag ) : emptyTagStr = '/'
        attributeStr = ''
        if( self.label is not None ) : attributeStr += ' label="%s"' % self.label
        return( '%s<%s%s%s%s>' % ( indent, self.moniker, attributeStr, extraAttributesAsStrings, emptyTagStr ) )
