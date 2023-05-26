# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

'''
This module contains the classes needed to represent the nodes of an isotopic abundance library.
An isotopic abundance library is a list of chemical elements with each chemical element having 
a unique list of natural occurring isotopes for that chemical element.  Each isotope in the library
has an atom fraction and uncertainly. If the uncertainly for an isotope is not known, its 
value must be -1.0.
'''

from pqu import PQU as PQUModule

from LUPY import enums as enumsModule
from LUPY import ancestry as ancestryModule
from xData.Documentation import documentation as documentationModule
from fudge import suites as suitesModule

class Format(enumsModule.Enum):
    '''
    This enum class represents the allowed formats for the IsotopicAbundancesByChemicalElement class.
    '''

    version_2_0 = '2.0'

    @staticmethod
    def default():
        '''Returns the current default format.'''

        return Format.version_2_0

class Isotope(ancestryModule.AncestryIO):
    '''
    This class represents an iam *isotope* node and stores the atom fraction and its uncertainty for an isotope.

    The following table list the primary members of this class:

        +---------------+-----------------------------------------------------------+
        | Member        | Description                                               |
        +===============+===========================================================+
        | id            | The Pops id for the isotope (e.g., "Fe56").               |
        +---------------+-----------------------------------------------------------+
        | atomFraction  | The atom fraction of *self*.                              |
        +---------------+-----------------------------------------------------------+
        | uncertainty   | The uncertainty of the atomFraction.                      |
        +---------------+-----------------------------------------------------------+
    '''

    moniker = 'isotope'
    keyName = 'id'

    def __init__(self, id, atomFraction, uncertainty):
        '''
        Constructor.

        :param id:              The PoPs id for the isotope.
        :param atomFraction:    The atom fraction of the isotope.
        :param uncertainty:     The uncertainly of the *atomFraction*.
        '''

        self.__id = id
        self.__atomFraction = float(atomFraction)
        self.__uncertainty = float(uncertainty)

    def __str__(self):

        return '    %-6s %s %s' % (self.__id, self.__atomFraction, self.__uncertainty)

    @property
    def label(self):
        '''Returns the label (i.e., id) for *self*.'''

        return self.id

    @property
    def id(self):
        '''Returns the id for *self*.'''

        return self.__id

    @property
    def atomFraction(self):
        '''Returns the atom fraction for *self*.'''

        return self.__atomFraction

    @property
    def uncertainty(self):
        '''Returns the uncertainty for *self*.'''

        return self.__uncertainty

    def toXML_strList(self, indent='', **kwargs):
        '''
        Converts *this* to a list of XML strings.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that control how *self* is converted to a list of XML strings.
        '''

        incrementalIndent = kwargs.get('incrementalIndent', '  ')
        indent2 = indent + incrementalIndent

        atomFraction = PQUModule.floatToShortestString(self.__atomFraction)
        uncertainty = PQUModule.floatToShortestString(self.__uncertainty, significantDigits=4)
        XML_string = ['%s<%s id="%s" atomFraction="%s" uncertainty="%s"/>' % (indent, self.moniker, self.__id, atomFraction, uncertainty)]

        return XML_string

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        '''
        Class method that parses a XML representation of an isotope node.

        :param cls:             The class used to create the isotope instance.
        :param node:            The XML node to parse.
        :param xPath:           The xPath to the creation of this instance.
        :param linkData:        Internal dictionary of extra arguments that aids in parsing.
        :param kwargs:          External dictionary of extra arguments that aids in parsing.
        '''

        isotope = cls(node.get('id'), node.get('atomFraction'), node.get('uncertainty'))

        return isotope

class Isotopes(suitesModule.ExclusiveSuite):
    '''
    This class represents an isotopes instanse.
    '''

    moniker = 'isotopes'
    keyName = 'id'

    def __init__(self):
        '''Constructor.'''

        suitesModule.ExclusiveSuite.__init__(self, (Isotope,))

class ChemicalElement(ancestryModule.AncestryIO):
    '''
    This class represents a chemical element which has a list of isotopes with "natural" abundance data
    for this chemical element and a documentation member.

    The following table list the primary members of this class:

        +-------------------+-----------------------------------------------------------+
        | Member            | Description                                               |
        +===================+===========================================================+
        | symbol            | The PoPs symbol for the chemical element (e.g., "Fe").    |
        +-------------------+-----------------------------------------------------------+
        | documentation     | A documentation node for *self*.                          |
        +-------------------+-----------------------------------------------------------+
        | __isotopes        | The list of isotopes for *self*                           |
        +-------------------+-----------------------------------------------------------+
    '''

    # Missing member uncertainty.

    moniker = 'chemicalElement'
    keyName = 'symbol'

    def __init__(self, symbol):
        '''
        Constructor.

        :param symbol:  The PoPs symbol of the chemical element.
        '''

        self.__symbol = symbol

        self.__documentation = documentationModule.Documentation()
        self.__documentation.setAncestor(self)

        self.__isotopes = Isotopes()
        self.__isotopes.setAncestor(self)

    def __getitem__(self, id):
        '''Returns the isotope with whose id is *id*.'''

        return self.__isotopes[id]

    def __iter__(self):
        '''Iterates over the entries in *self.__isotopes*.'''

        for isotope in self.__isotopes:
            yield isotope

    @property
    def label(self):
        '''Returns the label (i.e., symbol) for *self*.'''

        return self.symbol

    @property
    def symbol(self):
        '''Returns the symbol for *self*.'''

        return self.__symbol

    @property
    def documentation(self):
        '''Returns a reference to *self*'s documentation.'''

        return self.__documentation

    @property
    def isotopes(self):
        '''Returns a reference to *self*'s list of isotopes.'''

        return self.__isotopes

    def atomFractionSum(self):
        '''
        Calculations the sum of the atom factions of the isotopes in *self*.

        :return:                The sum of the atom fractions.
        '''

        atomFractionSum = 0
        for isotope in self.__isotopes:
            atomFractionSum += isotope.atomFraction

        return atomFractionSum

    def mass(self, pops, unit='amu'):
        '''
        Calculates the element mass for *self* from the isotope atomFractions data in *self*
        and the isotope masses given in *pops*.

        :param pops:        A FUDGE PoPs instance from which masses are looked up.
        :param unit:        The requested unit of the returned mass.

        :return:            The calculated element mass as a float.
        '''

        mass = 0.0
        for isotope in self:
            mass += isotope.atomFraction * pops[isotope.id].getMass(unit)

        return mass

    def toXML_strList(self, indent='', **kwargs):
        '''
        Converts *this* to a list of XML strings.

        :param indent:  The minimum amount of indentation.
        :param kwargs:  A dictionary of extra arguments that control how *self* is converted to a list of XML strings.
        '''

        incrementalIndent = kwargs.get('incrementalIndent', '  ')
        indent2 = indent + incrementalIndent

        XML_string = ['%s<%s symbol="%s">' % (indent, self.moniker, self.__symbol)]
        XML_string += self.__documentation.toXML_strList(indent=indent2, **kwargs)
        XML_string += self.__isotopes.toXML_strList(indent=indent2, **kwargs)
        XML_string[-1] += '</%s>' % self.moniker

        return XML_string

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        '''
        Class method that parses a XML representation of a chemicalElement node.

        :param cls:             The class used to create the isotope instance.
        :param node:            The XML node to parse.
        :param xPath:           The xPath to the creation of this instance.
        :param linkData:        Internal dictionary of extra arguments that aids in parsing.
        :param kwargs:          External dictionary of extra arguments that aids in parsing.
        '''

        chemicalElement = ChemicalElement(node.get('symbol'))
        for child in node:
            if child.tag == chemicalElement.__documentation.moniker:
                chemicalElement.__documentation.parseNode(child, xPath, linkData, **kwargs)
            elif child.tag == chemicalElement.__isotopes.moniker:
                chemicalElement.__isotopes.parseNode(child, xPath, linkData, **kwargs)
            else:
                raise Exception('Invalid child with tag name "%s".' % child.tag)

        return chemicalElement

class ChemicalElements(suitesModule.ExclusiveSuite):
    '''
    This class represents an chemicalElements instanse which contains a list of **ChemicalElement** instances.
    '''

    moniker = 'chemicalElements'
    keyName = 'symbol'

    def __init__(self):
        '''Constructor.'''

        suitesModule.ExclusiveSuite.__init__(self, (ChemicalElement,))

class IsotopicAbundancesByChemicalElement(ancestryModule.AncestryIO):
    '''
    The classes represents an evaluation with a list of chemical elements evaluated.

    The following table list the primary members of this class:

        +-------------------+-----------------------------------------------------------+
        | Member            | Description                                               |
        +===================+===========================================================+
        | evaluation        | The evaluation string for *self*.                         |
        +-------------------+-----------------------------------------------------------+
        | format            | The format of the file the data were stored in.           |
        +-------------------+-----------------------------------------------------------+
        | documentation     | The documentation for *self*.                             |
        +-------------------+-----------------------------------------------------------+
        | chemicalElements  | The list of unique chemical elements in *self*.           |
        +-------------------+-----------------------------------------------------------+
    '''

    moniker = 'isotopicAbundancesByChemicalElement'
    keyName = 'evaluation'

    def __init__(self, evaluation, format=Format.default()):
        '''
        Constructor.

        :param evaluation:      The evaluation for the **isotopicAbundancesByChemicalElement** instance.
        :param format:          The format the output will be written to.
        '''

        ancestryModule.AncestryIO.__init__(self)

        self.__evaluation = evaluation

        self.__format = Format.checkEnumOrString(format)

        self.__documentation = documentationModule.Documentation()
        self.__documentation.setAncestor(self)

        self.__chemicalElements = ChemicalElements()
        self.__chemicalElements.setAncestor(self)

    def __getitem__(self, symbol):
        '''Returns the chemical element with symbol *symbol*.'''

        return self.__chemicalElements[symbol]

    def __iter__(self):
        '''Iterates over the entries in *self.__chemicalElements*.'''

        for chemicalElement in self.__chemicalElements:
            yield chemicalElement

    @property
    def evaluation(self):
        '''Returns a reference to *self*'s evaluation.'''

        return self.__evaluation

    @property
    def format(self):
        '''Returns a reference to *self*'s format.'''

        return(self.__format)

    @property
    def documentation(self):
        '''Returns a reference to *self*'s documentation.'''

        return(self.__documentation)

    @property
    def chemicalElements(self):
        '''Returns a reference to *self*'s chemicalElements.'''

        return(self.__chemicalElements)

    def getIsotope(self, id):

        for chemicalElement in self.__chemicalElements:
            for isotope in chemicalElement.isotopes:
                if isotope.id == id:
                    return(isotope)

        return None

    def symbols(self):
        '''Returns a list of available chemical elements in *self*.'''

        return [chemicalElement.symbol for chemicalElement in self.__chemicalElements]

    def toXML_strList(self, indent = '', **kwargs):
        '''
        Converts *this* to a list of XML strings.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that control how *self* is converted to a list of XML strings.
        '''

        incrementalIndent = kwargs.get('incrementalIndent', '  ')
        indent2 = indent + incrementalIndent

        XML_string = ['%s<%s evaluation="%s" format="%s">' % (indent, self.moniker, self.__evaluation, self.__format)]
        XML_string += self.__documentation.toXML_strList(indent=indent2, **kwargs)
        XML_string += self.__chemicalElements.toXML_strList(indent=indent2, **kwargs)
        XML_string[-1] += '</%s>' % self.moniker

        return XML_string

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        '''
        Class method that parses a XML representation of an isotopicAbundancesByChemicalElement node.

        :param cls:             The class used to create the isotope instance.
        :param node:            The XML node to parse.
        :param xPath:           The xPath to the creation of this instance.
        :param linkData:        Internal dictionary of extra arguments that aids in parsing.
        :param kwargs:          External dictionary of extra arguments that aids in parsing.
        '''

        isotopicAbundancesByChemicalElement = IsotopicAbundancesByChemicalElement(node.get('evaluation'), node.get('format'))
        for child in node:
            if child.tag == isotopicAbundancesByChemicalElement.__documentation.moniker:
                isotopicAbundancesByChemicalElement.__documentation.parseNode(child, xPath, linkData, **kwargs)
            elif child.tag == isotopicAbundancesByChemicalElement.__chemicalElements.moniker:
                isotopicAbundancesByChemicalElement.__chemicalElements.parseNode(child, xPath, linkData, **kwargs)
            else:
                raise Exception('Invalid child with tag name "%s".' % child.tag)

        return isotopicAbundancesByChemicalElement

    @staticmethod
    def read(fileName, **kwargs):
        '''
        Reads in the file named *fileName* and returns a **IsotopicAbundancesByChemicalElement** instance of it.

        :param fileName:        The name of the file to read.
        :param kwargs:          External dictionary of extra arguments that aids in parsing.
        '''

        return IsotopicAbundancesByChemicalElement.readXML_file(fileName, **kwargs)

def read(fileName, **kwargs):
    '''
    Reads in the file named *fileName* and returns a **IsotopicAbundancesByChemicalElement** instance of it.

    :param fileName:            The name of the file to read.
    :param kwargs:              External dictionary of extra arguments that aids in parsing.
    '''

    return IsotopicAbundancesByChemicalElement.read(fileName, **kwargs)
