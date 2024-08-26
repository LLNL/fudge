# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains classes that are base classes for some of the double differential cross section classes.

This module contains the following classes:
        
    +---------------------------------------+-----------------------------------------------------------------------------------+
    | Class                                 | Description                                                                       |
    +=======================================+===================================================================================+
    | Form                                  | This class is a base class for the Form classes in the sub-directories.           |
    +---------------------------------------+-----------------------------------------------------------------------------------+
    | XYs1d                                 | This class is a base class for several XYs1d classes.                             |
    +---------------------------------------+-----------------------------------------------------------------------------------+
    | Regions1d                             | This class is a base class for several Regions1d classes.                         |
    +---------------------------------------+-----------------------------------------------------------------------------------+
"""

from xData import enums as xDataEnumsModule
from xData import axes as axesModule
from xData import XYs1d as XYs1dModule
from xData import regions as regionsModule

from fudge import abstractClasses as abstractClassesModule

class Form( abstractClassesModule.Form ) :
    """
    This class is a base class for the Form classes in the various sub-directories.

    The following table list the primary members of this class:

    +-----------------------+-------------------------------------------------------------------+
    | Member                | Description                                                       |
    +=======================+===================================================================+
    | pid                   | The GNDS PoPs id of the product.                                  |
    +-----------------------+-------------------------------------------------------------------+
    | label                 | The label for this form.                                          |
    +-----------------------+-------------------------------------------------------------------+
    | productFrame          | The frame the product data are specified in.                      |
    +-----------------------+-------------------------------------------------------------------+
    | identicalParticles    | If True, both products in a two-body reaction are identical,      |
    |                       | otherwise they are not.                                           |
    +-----------------------+-------------------------------------------------------------------+
    """

    def __init__( self, pid, label, productFrame, subForms, identicalParticles = False ) :
        """
        :param pid:                     The GNDS PoPs id of the product.
        :param label:                   The label for this form.
        :param productFrame:            The frame the product data are specified in.
        :param subForms:                A python list of sub-forms to load into *self*.
        :param identicalParticles:      If True, both products in a two-body reaction are identical, otherwise they are not.
        """

        abstractClassesModule.Form.__init__( self )

        if( not( isinstance( pid, str ) ) ) : raise TypeError( 'pid must be a str instance' )
        self.__pid = pid

        self.label = label

        self.__productFrame = xDataEnumsModule.Frame.checkEnumOrString(productFrame)

        self.__identicalParticles = identicalParticles

        for i1, subform in enumerate( subForms ) :
            setattr( self, self.subformAttributes[i1], subform )
            if( subform is not None ) : subform.setAncestor( self )
        self.subforms = [ subform for subform in subForms ]

    @property
    def label( self ) :
        """This method returns the label of *self*."""

        return( self.__label )

    @label.setter
    def label( self, value ) :
        """
        This method sets the label of *self* to *value*.

        :param value:       The value of the new label for *self*.
        """

        if( value is not None ) :
            if( not( isinstance( value, str ) ) ) : raise TypeError( 'value must be a string' )
        self.__label = value

    @property
    def pid( self ) :
        """This method returns the pid of *self*."""

        return( self.__pid )

    @property
    def productFrame( self ) :
        """This method returns the productFrame of *self*."""

        return( self.__productFrame )

    @property
    def identicalParticles( self ) :
        """This method returns the identicalParticles of *self*."""

        return( self.__identicalParticles )

    def isThermalNeutronScatteringLaw( self ) :
        """This method always returned False but will be overwritten by some derived classes."""

        return( False )

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        if( hasattr( self, 'subforms' ) ) :
            for subform in self.subforms :
                if( subform is None ) : continue           # FIXME, needed for photo-atomic coherentScattering data with no anomalousScattering subforms.
                subform.convertUnits( unitMap )

    def toXML_strList( self, indent = '', **kwargs ) :
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attributeStr = ''
        if( self.label is not None ) : attributeStr = ' label="%s"' % self.label
        attributeStr += ' pid="%s" productFrame="%s"' % ( self.__pid, self.__productFrame )
        if( self.identicalParticles ) : attributeStr += ' identicalParticles="true"'
        xmlString = [ '%s<%s%s>' % ( indent, self.moniker, attributeStr ) ]
        for subform in self.subforms : 
            if( subform is not None ) : xmlString += subform.toXML_strList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

class XYs1d( XYs1dModule.XYs1d ) :
    """This class is a base class for several XYs1d classes in modules below this module's location."""

    pass

class Regions1d( regionsModule.Regions1d ) :
    """This class is a base class for several Regions1d classes in modules below this module's location."""

    pass

def defaultAxes( factorLabel, energyUnit ) :
    """This function does not seem to be used and should be deleted."""

    axes = axesModule.Axes(2)
    axes[0] = axesModule.Axis( factorLabel, 0, "" )
    axes[1] = axesModule.Axis( 'energy_in', 1, energyUnit )
    return( axes )

def check( self, info ) :
    """This function does not seem to be used and should be deleted."""

    return( [] )
