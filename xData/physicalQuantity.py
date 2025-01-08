# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the :py:class:`PhysicalQuantity` class which is used to represents a number with a unit (e.g., '1 cm', '5.4 kg').

This module contains the following classes:

    +-------------------+-------------------------------------------------------------------------------------------------------+
    | Class             | Description                                                                                           |
    +===================+=======================================================================================================+
    | PhysicalQuantity  | This class represents a physical quantity which is a number with a unit (e.g., '1 cm', '5.4 kg').     |
    +-------------------+-------------------------------------------------------------------------------------------------------+
"""

from LUPY import ancestry as ancestryModule

from pqu import PQU as PQUModule

from .uncertainty.physicalQuantity import uncertainty as uncertaintyModule

class PhysicalQuantity( ancestryModule.AncestryIO ) :
    """
    This class is use to store a physical quantity which is a number with a unit (e.g., '1 cm', '5.4 kg'), and
    to support various operations on the physical quantity.

    The following table list the primary members of this class:

    +---------------+-------------------------------------------------------------------------------+
    | Member        | Description                                                                   |
    +===============+===============================================================================+
    | value         | This is the value for the physical quantity.                                  |
    +---------------+-------------------------------------------------------------------------------+
    | unit          | This is the unit for the physical quantity.                                   |
    +---------------+-------------------------------------------------------------------------------+
    | label         | This is the label member.                                                     |
    +---------------+-------------------------------------------------------------------------------+
    | uncertainty   | This is the uncertainy of value.                                              |
    +---------------+-------------------------------------------------------------------------------+
    | PQ            | A :py:class:`PQUModule.PQU` called with *self*' value and unit.               |
    +---------------+-------------------------------------------------------------------------------+
    """

    moniker = 'physicalQuantity'

    def __init__( self, value, unit, label = None ) :
        """
        :param value:   This is the value of the physical quantity.
        :param unit:    This is the unit of the physical quantity.
        :param label:   This is the label of the physical quantity.
        """

        ancestryModule.AncestryIO.__init__( self )
        self.__PQ = PQUModule.PQU( value, unit )

        if( label is not None ) :
            if( not( isinstance( label, str ) ) ) : raise TypeError( 'label must be a str instance.' )
        self.__label = label

        self.__uncertainty = None

    def __str__( self ) :
        """
        This method returns the value and unit of *self* as a str.

        :returns:       A python str.
        """

        return( self.toString( ) )

    def __float__( self ) :
        """
        This method returns *self*'s value as a float.

        :returns:       A python float.
        """

        return( float( self.__PQ ) )

    def __eq__( self, other ) :
        """
        This method returns True if *self* is equal to *other* and False otherwise.
        If *other* is not an instance of :py:class:`PhysicalQuantity`, False is returned. If *self*'s unit is not compatible
        with *other*'s unit, False is returned. Otherwise, returns True if :py:func:`compare` with epsilonFactor = 5 returns
        0 and False otherwise.

        :param other:   A :py:class:`PhysicalQuantity` to compare to *self*.

        :returns:       A python boolean.
        """

        if( not isinstance( other, PhysicalQuantity ) ) : return( False )
        if( not( self.__PQ.isCompatible( other.__PQ.unit ) ) ) : return( False )
        return( self.compare( other, 5 ) == 0 )

    def __ne__( self, other ) :
        """
        This method returns True if *self* is not equal to *other* and False otherwise.
        If *other* is not an instance of :py:class:`PhysicalQuantity`, False is returned. If *self*'s unit is not compatible
        with *other*'s unit, False is returned. Otherwise, returns False if :py:func:`compare` with epsilonFactor = 5 returns
        0 and True otherwise.

        :param other:   A :py:class:`PhysicalQuantity` to compare to *self*.

        :returns:       A python boolean.
        """

        if( not isinstance( other, PhysicalQuantity ) ) : return( True )
        if( not( self.__PQ.isCompatible( other.__PQ.unit ) ) ) : return( True )
        return( self.compare( other, 5 ) != 0 )

    def __lt__( self, other ) :
        """
        This method returns True if *self* is less than *other* and False otherwise.
        If *other* is not an instance of :py:class:`PhysicalQuantity`, False is returned. If *self*'s unit is not compatible
        with *other*'s unit, False is returned. Otherwise, returns True if :py:func:`compare` with epsilonFactor = 5 returns
        a negative number and False otherwise.

        :param other:   A :py:class:`PhysicalQuantity` to compare to *self*.

        :returns:       A python boolean.
        """


        return( self.compare( other, 5 ) < 0 )

    def __le__( self, other ) :
        """
        This method returns True if *self* is less than or equal to *other*, and False otherwise.
        If *other* is not an instance of :py:class:`PhysicalQuantity`, False is returned. If *self*'s unit is not compatible
        with *other*'s unit, False is returned. Otherwise, returns True if :py:func:`compare` with epsilonFactor = 5 returns
        a negative number or 0, and False otherwise.

        :param other:   A :py:class:`PhysicalQuantity` to compare to *self*.

        :returns:       A python boolean.
        """

        return( self.compare( other, 5 ) <= 0 )

    def __gt__( self, other ) :
        """
        This method returns True if *self* is greater than *other* and False otherwise.
        If *other* is not an instance of :py:class:`PhysicalQuantity`, False is returned. If *self*'s unit is not compatible
        with *other*'s unit, False is returned. Otherwise, returns True if :py:func:`compare` with epsilonFactor = 5 returns
        a positive number, and False otherwise.

        :param other:   A :py:class:`PhysicalQuantity` to compare to *self*.

        :returns:       A python boolean.
        """

        return( self.compare( other, 5 ) > 0 )

    def __ge__( self, other ) :
        """
        This method returns True if *self* is greater than or equal to *other* and False otherwise.
        If *other* is not an instance of :py:class:`PhysicalQuantity`, False is returned. If *self*'s unit is not compatible
        with *other*'s unit, False is returned. Otherwise, returns True if :py:func:`compare` with epsilonFactor = 5 returns
        a positive number or 0, and False otherwise.

        :param other:   A :py:class:`PhysicalQuantity` to compare to *self*.

        :returns:       A python boolean.
        """

        return( self.compare( other, 5 ) >= 0 )

    def compare( self, other, epsilonFactor = 0 ) :
        """
        This method calls :py:func:`PQUModule.PQU.compare` with epsilonFactor = 5 and returns the results. 
        *Other* must be an instance of :py:func:`PhysicalQuantity`.
        If *self* is less than *other* -1 is returned. If *self* is equal to *other* 0 is returned. Otherwise, 1 is returned.

        :param other:   A :py:class:`PhysicalQuantity` to compare to *self*.

        :returns:       A python in that is one of -1, 0 or 1.
        """

        if( not isinstance( other, PhysicalQuantity ) ) : raise TypeError( 'other not instance of PhysicalQuantity' )
        return( self.__PQ.compare( other.__PQ, epsilonFactor ) )

    @property
    def key( self ) :
        """
        This method returns the keyName for *self*.

        :returns:       A python str.
        """

        return( self.__label )

    @property
    def label( self ) :
        """
        This method returns the label for *self*.

        :returns:       A python str.
        """

        return( self.__label )

    @property
    def value( self ) :
        """
        This method returns *self*'s value as a float.

        :returns:       A python float.
        """

        return( float( self ) )

    @property
    def unit( self ) :
        """
        This method returns *self*'s unit.

        :returns:       A python str.
        """

        return( self.__PQ.unit )

    @property
    def uncertainty( self ) :
        """
        This method returns a reference to *self*'s uncertainty.

        :returns:       An instance of :py:class:`uncertaintyModule.Uncertainty` or None.
        """

        return( self.__uncertainty )

    @uncertainty.setter
    def uncertainty( self, _uncertainty ) :
        """
        This method sets *self*'s uncertainty to *_uncertainty*.

        :param _uncertainty:        An instance of :py:class:`uncertaintyModule.Uncertainty` or None.
        """

        if( _uncertainty is not None ) :
            if( not( isinstance( _uncertainty, uncertaintyModule.Uncertainty ) ) ) : raise TypeError( 'Invalid uncertainty instance.' )

        self.__uncertainty = _uncertainty
        if( self.__uncertainty is not None ) : self.__uncertainty.setAncestor( self )

    def convertToUnit( self, unit ) :
        """
        This method changes *self*'s unit to *unit* and updates its value accordingly.

        :param unit:    The new unit for *self*'s value.
        """

        self.__PQ.convertToUnit( unit )

    def convertUnits( self, unitMap ) :
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        unit, factor = PQUModule.convertUnits( self.unit, unitMap )
        self.__PQ.convertToUnit( unit )
        if( self.__uncertainty is not None ) : self.__uncertainty.parentConvertingUnits( [ factor ] )

    def copy( self ) :
        """
        This method returns a copy of *self*.
        """

        cls = self.__class__( self.value, self.unit, self.label )
        if( self.__uncertainty is not None ) : cls.uncertainty = self.__uncertainty.copy( )
        return( cls )

    __copy__ = copy

    def copyToUnit( self, unit ) :
        """
        This method returns a copy of *self* with new unit *unit*.

        :param unit:    The unit of the copy of *self*.

        :returns:       An instance of :py:class:`PhysicalQuantity`.
        """

        _copy = self.copy( )
        _copy.convertToUnit( unit )
        return( _copy )

    def getValueAs( self, unit ) :
        """
        This method returns a float that is *self*'s value converted to *unit*.

        :param uhit:        The unit of the returned value.

        :returns:           A python float.
        """

        return( self.__PQ.getValueAs( unit ) )

    def toString( self, significantDigits = None, keepPeriod = True ) :
        """
        This method returns the value and unit of *self* as a str.

        :returns:       A python str.
        """

        return( self.__PQ.toString( significantDigits = significantDigits, keepPeriod = keepPeriod ) )

    def toXML_strList(self, indent='', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get('incrementalIndent', '  ')
        moniker = kwargs.get('moniker', self.moniker)

        keyNameValue = ''                                       # The next few lines are a kludge until keyName is properly supported.
        keyName = self.keyName
        if keyName is None:
            keyName = 'label'
            keyValue = self.label
        else:
            keyValue = self.keyValue
        if keyValue is not None:
            keyNameValue = ' %s="%s"' % (keyName, keyValue)

        unit = ''
        if not self.unit.isDimensionless(): unit = ' unit="%s"' % self.unit

        ending = '>'
        if self.__uncertainty is None: ending = '/>'
        XML_strList = [ '%s<%s%s value="%s"%s%s' % ( indent, moniker, keyNameValue, PQUModule.floatToShortestString(self.value, 12), unit, ending ) ]

        if ending == '>':
            if self.__uncertainty is not None: XML_strList += self.__uncertainty.toXML_strList(indent = indent2, **kwargs)
            XML_strList[-1] += '</%s>' % moniker

        return XML_strList

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance of *cls*.

        :param cls:         Form class to return.
        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :returns:           An instance of *cls* representing *node*.
        """

        xPath.append(node.tag)

        value = node.get('value')
        unit = node.get('unit')
        if cls.keyName is None:
            keyValue = node.get('label', None)
        else:
            keyValue = node.get(cls.keyName, None)
        instance = cls(value, unit, keyValue)

        for child in node:
            if child.tag == uncertaintyModule.Uncertainty.moniker:
                instance.uncertainty = uncertaintyModule.Uncertainty.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            else:
                raise ValueError('Invalid child node with tag "%s".' % child.tag)

        xPath.pop() 

        return instance
