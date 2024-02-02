# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains abstract base classes for distribution forms and sub-forms.

This module contains the following classes:
        
    +---------------+-----------------------------------------------------------------------+
    | Class         | Description                                                           |
    +===============+=======================================================================+ 
    | Form          | Abstract base class for all forms.                                    |
    +---------------+-----------------------------------------------------------------------+
    | Subform       | Abstract base class for all sub-forms.                                |
    +---------------+-----------------------------------------------------------------------+
"""

from LUPY import ancestry as ancestryModule
from xData import enums as xDataEnumsModule

from PoPs.chemicalElements import misc as chemicalElementMiscPoPsModule

from fudge import abstractClasses as abstractClassesModule

class Form( abstractClassesModule.Form ) :
    """
    Abstract base class for all forms.

    The following table list the primary members of this class:

    +---------------+---------------------------------------------------------------+
    | Member        | Description                                                   |
    +===============+===============================================================+
    | label         | Unique label for the form within the distribution instance.   |
    +---------------+---------------------------------------------------------------+
    | productFrame  | The frame the product data are in.                            |
    +---------------+---------------------------------------------------------------+
    | subForms      | The list of sub-forms for the form.                           |
    +---------------+---------------------------------------------------------------+

    :param label:               The label for this form.
    :param productFrame:        The frame the product data are specified in.
    :param subForms:            The list of sub-forms for the form. 
    """

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
        """
        Returns the label for *self*.
        """

        return( self.__label )

    @label.setter
    def label( self, value ) :
        """
        Sets the label of *self* to *value*.

        :param value:   The value of the new label for *self*.
        """

        if( value is not None ) :
            if( not( isinstance( value, str ) ) ) : raise TypeError( 'value must be a string' )
        self.__label = value

    @property
    def product( self ) :
        """
        Returns the :py:class:`productModule.Product` instance that *self* is in.

        :return:    A :py:class:`productModule.Product` instance.
        """

        from fudge import product as productModule

        return self.findClassInAncestry(productModule.Product)

    @property
    def productFrame( self ) :
        """Returns the frame the product data are specified in."""

        return( self.__productFrame )

    @productFrame.setter
    def productFrame(self, frame):
        """
        Sets the frame the product data are specified in.

        :param frame:       The frame for the product.
        """

        self.__productFrame = xDataEnumsModule.Frame.checkEnumOrString(frame)

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        if( hasattr( self, 'subforms' ) ) :
            for subform_ in self.subforms :
                if( subform_ is None ) : continue                    # FIXME, needed for photo-atomic coherentScattering data with no anomalousScattering subforms.
                subform_.convertUnits( unitMap )

    def findEntity( self, entityName, attribute = None, value = None ):
        """
        Overrides :py:meth:`ancestryModule.Ancestry.findEntity`.

        :param entityName:  Name is the desired entry in *self*.
        :param attribute:   Name of an attribute in *self*.
        :param value:       Value of the named attribute.

        :return:            The entry matching *entityName* or *attribute*/*value*.
        """

        for subform_ in self.subforms :
            if( subform_.moniker == entityName ) : return( subform_ )
        return( ancestryModule.Ancestry.findEntity( self, entityName, attribute, value ) )

    def isTwoBody( self ) :
        """
        Returns True if the form represents a two-body interaction and False otherwise.

        :return:    Boolean.
        """

        return( False )

    def isSpecified( self ) :
        """
        This method returns *True*. However, this method is overwritten in the :py:class:`Unspecified` class and that
        method returns False.  That is, for all distribution classes except the :py:class:`Unspecified`, the 
        **isSpecified** method return True.

        :return:                Boolean.
        """

        return( True )

    def residualMass( self, PoPs, residualZA, massUnit, compoundMass, product = None, verbose = 1 ) :
        """
        Calculates the mass of the residual from the specified data. If the residual is not in the
        is not in *PoPs*, the residual mass is estimated.

        :param PoPs:            A GNDS PoPs instance.
        :param residualZA:      The ZA of the residual.
        :param massUnit:        The unit of the mass value returned.
        :param compoundMass:    The mass of the compound.
        :param product:         The :py:class:`productModule.Product` other than *self* is desired.
        :param verbose:         If greater than 0, a warning message if residual is not in *PoPs*.
        """

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
        """
        Returns a pointwise represent of *self*.

        :param kwargs:              A dictionary that contains data to control the way this method acts.

        :return:                    ?
        """

        return( self.subforms[0].toPointwise_withLinearXYs( **kwargs ) )

    def toXML_strList( self, indent = '', **kwargs ) :
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The amount of indentation for each line. Child nodes and text may be indented more.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

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
    """
    Abstract base class for all distribution sub-forms.
    """

    def __init__( self ) :

        ancestryModule.AncestryIO.__init__( self )
        self.label = None

    def XMLStartTagString( self, indent = '', extraAttributesAsStrings = '', emptyTag = False ) :
        """
        Returns the XML start tag of *self*.

        :return:        Python str instance for the XML start tag of *self*.
        """

        emptyTagStr = ''
        if( emptyTag ) : emptyTagStr = '/'
        attributeStr = ''
        if( self.label is not None ) : attributeStr += ' label="%s"' % self.label
        return( '%s<%s%s%s%s>' % ( indent, self.moniker, attributeStr, extraAttributesAsStrings, emptyTagStr ) )
