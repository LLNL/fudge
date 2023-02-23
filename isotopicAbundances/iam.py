# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

'''
This module contains the classes needed to represent the nodes of an isotopic abundance map (iam) instance.
'''

import pathlib 

from LUPY import misc as LUPY_miscModule
from LUPY import ancestry as ancestryModule
from LUPY import checksums as checksumsModule

from . import isotopicAbundances as isotopicAbundancesModule

class Base(ancestryModule.AncestryIO):
    '''Base class used by all other iam related classes. Inherits class ancestryModule.AncestryIO, and adds checksum and algorithm properties.'''

    def __init__(self, checksum=None, algorithm=None):
        '''
        :param checksum:        The check sum of the file referrence by the path member.
        :param algorithm:       The algorithm  used to calculate *checksum*.
        '''

        ancestryModule.AncestryIO.__init__(self)

        self.checksum = checksum
        self.algorithm = algorithm

    @property
    def checksum(self):
        '''Returns self's checksum.'''

        return self.__checksum

    @checksum.setter
    def checksum(self, checksum):
        '''Sets self's checksum to *checksum*.'''

        if checksum is not None and not isinstance(checksum, str):
            raise TypeError('Checksum must be a string.')
        self.__checksum = checksum

    @property
    def algorithm(self):
        '''Returns the algorithm used to compute the checksum for this entry.'''

        if self.__algorithm:
            return self.__algorithm
        if isinstance(self, Iam):
            return None

        return self.ancestor.algorithm

    @algorithm.setter
    def algorithm(self, algorithm):
        '''Set self's algorithm to *algorithm*.'''

        if algorithm not in checksumsModule.supportedAlgorithms:
            raise ValueError('Unsupported checksum algorithm "%s".' % algorithm)
        self.__algorithm = algorithm

    @property
    def fileName(self):
        '''Returns the file name for the location of self. Unlike path, this will always be the absolute path.'''

        filename = self.path
        if self.ancestor is not None:
            if not pathlib.Path(filename).is_absolute():
                filename = str(pathlib.Path(self.ancestor.path).resolve().parent / filename)

        return filename

    def standardXML_attributes(self):
        '''
        Returns the XML attribute string for *checksum* and *algorithm*. For internal use.
        '''

        attributes = ''
        if self.checksum:
            attributes += ' checksum="%s"' % self.checksum

        algorithm = self.algorithm
        if self.ancestor is not None:
            if not isinstance(self, Iam) and algorithm == self.ancestor.algorithm:
                algorithm = None
        if algorithm is not None:
            attributes += ' algorithm="%s"' % algorithm

        return attributes

class Iam(Base):
    '''This class represents an Isotopic Abundance Map (IAM) instance.'''

    moniker = 'iam'

    def __init__(self, library, path, format=isotopicAbundancesModule.Format.default(), checksum=None, algorithm=None):
        ''' 
        Constructor for Iam class.

        :param library:         The name of the library.
        :param path:            The path to the iam file read in to construct *self*. If created anew, this can be the string "./".
        :param format:          The format version for the iam.
        :param checksum:        Hash for the iam (computed from the concatenation of all checksums inside the iam)
        :param algorithm:       Algorithm used to compute checksums.
        '''

        Base.__init__(self, checksum=checksum, algorithm=algorithm)

        self.__library = LUPY_miscModule.isString(library, 'Library must be a string.')
        self.__path = LUPY_miscModule.isString(path, 'Path must be a string.')
        self.__format = isotopicAbundancesModule.Format.checkEnumOrString(format)
        self.__entries = []

    def __getitem__(self, index):
        '''
        Returns the (index-1)^th item of self.

        :param index:           The index of the iter to return.
        '''

        return self.__entries[index]

    def __iter__(self):
        '''Iterators over each entry in self. Does not dive into import entries.'''

        for entry in self.__entries:
            yield entry

    def __len__(self):
        '''Returns the number of entries in self. Does not dive into import entries.'''

        return len(self.__entries)

    @property
    def path(self):
        '''Returns the path this iam was read from or may be written to. It may be absolote or relative, depending on how it was initialize.'''

        return self.__path

    @property
    def format(self):
        '''Returns to format for self.'''

        return self.__format

    @property
    def library(self):
        '''Returns the library name for self.'''

        return self.__library

    def updateAllChecksums(self, algorithm=checksumsModule.Sha1sum.algorithm, iamDirectory=None):
        '''
        Calls *updateChecksum* on all entries and then updates self's *checksum*.

        :param algorithm:       Algorithm that will be used to compute checksums.
        :param iamDirectory:    See documentation for the method **buildFileName** in class **EntryBase**.
        '''

        for entry in self:
            entry.updateChecksum(algorithm, iamDirectory=iamDirectory)
        self.updateChecksum(algorithm)

    def updateChecksum(self, algorithm=checksumsModule.Sha1sum.algorithm):
        '''
        Computes iam file's checksum: concatenate checksums for all entries into a string and compute the checksum
        of that string.

        :param algorithm:       Algorithm that will be used to compute checksums.
        '''

        s1 = ''.join([entry.checksum for entry in self if entry.checksum is not None])
        self.checksum = checksumsModule.checkers[algorithm].from_string(s1)
        self.algorithm = algorithm

    def append(self, entry):
        '''Appends *entry* to self.'''

        if not isinstance(entry, EntryBase):
            TypeError('Invalid entry.')

        self.__entries.append(entry)
        entry.setAncestor(self)

    def insert(self, index, entry):
        '''Inserts the entry into self at index.'''

        if not isinstance(entry, EntryBase):
            TypeError('Invalid entry of type "%s".' % type(entry))

        self.__entries.insert(index, entry)
        entry.setAncestor(self)

    def iterate(self):
        '''Iterates over all IsotopicAbundancesByChemicalElement in *self*. Dives into import entries.'''

        for entry in self.__entries:
            if isinstance(entry, Import):
                yield from entry.iterate()
            else:
                yield entry

    def evaluations(self):
        '''Returns a list of available evaluations.'''

        evaluations = []
        return [isotopicAbundancesByChemicalElement.evaluation for isotopicAbundancesByChemicalElement in self.iterate()]

    def find(self, symbol, evaluation=None):
        '''
        Returns the first entry matching a chemical elememnt's *symbol*. Searches each entry in the order they
        were appended and dives into imported iams.

        :param symbol:          The requested chemical element's PoPs symbol to match.
        :param evaluation:      The name of the evaluation to match. Can be None to match the first *symbol* found.

        :return:                Returns the found entry or None is no match was found.
        '''

        for entry in self.__entries:
            foundChemicalElement =  entry.find(symbol, evaluation)
            if foundChemicalElement is not None:
                return foundChemicalElement

        return None

    def findEvaluation(self, evaluation):
        '''Returns an IsotopicAbundancesByChemicalElement from from the module isotopicAbundances.py matching *evaluation* or None if no match is found.'''

        for isotopicAbundancesByChemicalElement in self.iterate():
            if isotopicAbundancesByChemicalElement.evaluation == evaluation:
                return isotopicAbundancesByChemicalElement.isotopicAbundancesByChemicalElement

        return None

    def toXML_strList(self, indent='', **kwargs):
        '''
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The amount of indentation for each line. Child nodes and text may be indented more.
        :param kwargs:          A keyword list.

        :return:                List of str instances representing the XML lines of self.
        '''

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        format = isotopicAbundancesModule.Format.checkEnumOrString(kwargs.get('format', kwargs.get('formatVersion', isotopicAbundancesModule.Format.default())))

        attrs = self.standardXML_attributes()

        XML_list = ['%s<%s library="%s" format="%s"%s>' % (indent, self.moniker, self.library, format, attrs)]
        for entry in self.__entries:
            XML_list += entry.toXML_strList(indent2, **kwargs)
        XML_list[-1] +=  '</%s>' % self.moniker

        return XML_list

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        '''
        Creates a Iam instance from an iam node.

        :param node:            node to parse.
        :param xPath:           Currently not used.
        :param linkData:        Currently not used.
        :param kwargs:          A keyword list.

        :return:                Iam instance.
        '''

        if node.tag != Iam.moniker:
            raise TypeError('Invalid node name "%s" for a iam instance.' % node.tag)

        sourcePath = kwargs['sourcePath']

        library = node.get('library', None)
        format = isotopicAbundancesModule.Format.checkEnumOrString(node.get('format'))

        iam = Iam(library, sourcePath, format=format, checksum=node.get('checksum'), algorithm=node.get('algorithm'))

        kwargs['format'] = format
        for child in node:
            if child.tag == Import.moniker:
                iam.append(Import.parseNodeUsingClass(child, xPath, linkData, **kwargs))
            elif child.tag == IsotopicAbundancesByChemicalElement.moniker:
                iam.append(IsotopicAbundancesByChemicalElement.parseNodeUsingClass(child, xPath, linkData, **kwargs))
            else:
                raise ValueError('Invalid child tag "%s" for iam file.' % child.tag)

        return iam

    @staticmethod
    def read(fileName, **kwargs):
        '''
        Reads in the file name *fileName* and returns a **Iam** instance.

        :param kwargs:          A keyword list.

        :return:                Iam instance.
        '''

        return Iam.readXML_file(fileName, **kwargs)

class EntryBase(Base):
    '''Base class for all ian entry classes.'''

    def __init__(self, path, checksum=None, algorithm=None):
        '''
        Base constructor for map entries.

        :param path:            Path attribute for the entry.
        :param checksum:        Checksum for this entry, computed using the indicated algorithm.
        :param algorithm:       Algorithm for computing the checksum. Only required if different from parent.
        '''

        Base.__init__(self, checksum=checksum, algorithm=algorithm)
        self.__path = LUPY_miscModule.isString(path, 'Path must be a string.')

    @property
    def path(self):
        '''Returns the path of the file. This may be absolute or relative.'''

        return self.__path

    def buildFileName(self, iamDirectory=None):
        '''
        This method is designed to aid in building an iam file when the iam file does not reside in its final resting place.
        The value of *iamDirectory* should be the final resting place (i.e., directory) of the iam file.
        Otherwise, this method returns the file name as the method *fileName*.

        :param iamDirectory:    Alternative directory where the iam file will be placed.
        '''

        if iamDirectory is None:
            return self.fileName

        return str(pathlib.Path(iamDirectory).revolve() / self.path)

    def updateChecksum(self, algorithm=checksumsModule.Sha1sum.algorithm, iamDirectory=None):
        '''
        Compute the checksum of the file specified by *path* member and store it into the *checksum* member.

        :param algorithm:       The algorithm to use for the checksum.
        :param iamDirectory:    See documentation for the method **buildFileName**.
        '''

        self.algorithm = algorithm
        self.checksum = checksumsModule.checkers[algorithm].from_file(self.buildFileName(iamDirectory))

    def standardXML_attributes(self):
        '''Returns the XML attribute string for *this*.'''

        return ' path="%s"%s' % (self.path, Base.standardXML_attributes(self))

class Import(EntryBase):
    '''This class represents an import instance within an Isotopic Abundance Map (IAM) instance.'''

    moniker = 'import'

    def __init__(self, path, checksum=None, algorithm=None):
        '''
        Constructor for the import entry.

        :param path:            Path to the iam file to import.
        :param checksum:        The check sum of the file referrence by the path member.
        :param algorithm:       The algorithm  used to calculate *checksum*.
        '''

        EntryBase.__init__(self, path, checksum=checksum, algorithm=algorithm)
        self.__iam = None

    def __str__(self):
        '''Returns a simple string representation of *self*.'''

        return '%s with path "%s".' % (self.moniker, self.path)

    @property
    def derivedPath(self):
        '''Returns the parent directory of the Iam instance containing *self*.'''

        ancestor = self.ancestor
        if isinstance(ancestor, Iam):
            return str(pathlib.Path(ancestor.path).parent)

        raise TypeError('Import not an entry of a Iam instance.')

    @property
    def iam(self):
        '''Returns a reference to self.__iam.'''

        self.readIam()

        return self.__iam

    def find(self, symbol, evaluation=None):
        '''
        Calls find on self.__iam and returns its results.

        :param symbol:          The requested chemical element's PoPs symbol to match.
        :param evaluation:      The name of the evaluation to match.

        :return:                Returns the found entry or None is no match was found.
        '''

        self.readIam()

        return self.__iam.find(symbol, evaluation)

    def iterate(self):
        '''Iterates over all entries of self's iam.'''

        self.readIam()
        yield from self.iam.iterate()

    def readIam(self):
        '''Reads in the iam file pointed to by self if not already read in. An import only reads its iam file when needed.'''

        if self.__iam is None:
            iamFilePath = self.path
            if not pathlib.Path(iamFilePath).is_absolute():
                iamFilePath = str(pathlib.Path(self.derivedPath) / iamFilePath)
            self.__iam = Iam.read(iamFilePath)
            self.__iam.setAncestor(self)

        return self.__iam

    def toXML_strList(self, indent = '', **kwargs):
        '''
        Returns a list of str instances representing the XML lines of self.

        :param indent:          The amount of indentation for each line. Child nodes and text may be indented more.
        :param kwargs:          A keyword list.

        :return:                List of str instances representing the XML lines of self.
        '''

        return ['%s<%s%s/>' % (indent, self.moniker, self.standardXML_attributes())]

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        '''
        Creates an Import instance from an import node.

        :param node:            The XML node to parse.
        :param xPath:           Currently not used.
        :param linkData:        Currently not used.
        :param kwargs:          A keyword list.

        :return:                Import instance.
        '''

        if node.tag != Import.moniker:
            raise TypeError('Invalid node name.')

        _import = Import(node.get('path'), checksum=node.get('checksum'), algorithm=node.get('algorithm'))

        return _import

class IsotopicAbundancesByChemicalElement(EntryBase):
    '''This class represents an isotopicAbundancesByChemicalElement instance within an Isotopic Abundance Map (IAM) instance.'''

    moniker = 'isotopicAbundancesByChemicalElement'

    def __init__(self, evaluation, path, checksum=None, algorithm=None):
        '''
        Construtor for IsotopicAbundancesByChemicalElement instance.

        :param evaluation:      Name for the evaluation.
        :param path:            Path to the file.
        :param checksum:        The checksum for the file referenced by *path*.
        :param algorithm:       The algorithm used to calculate the checksum.
        '''

        EntryBase.__init__(self, path, checksum=checksum, algorithm=algorithm)

        if not isinstance(evaluation, str):
            raise TypeError('Evaluation must be a string.')
        self.__evaluation = evaluation

        self.__isotopicAbundancesByChemicalElement = None

    def __str__(self):
        '''Returns a simple string representation of self.'''

        return '%s: evaluation "%s" and path "%s".' % (self.moniker, self.evaluation, self.path)

    @property
    def evaluation(self) :
        '''Returns the self's evaluation string.'''

        return self.__evaluation

    @property
    def derivedPath(self):
        '''Returns the parent directory of the instance containing *self*.'''

        ancestor = self.ancestor
        if pathlib.Path(self.path).is_absolute() or not isinstance(ancestor, Base):
            return self.path

        return str(pathlib.Path(ancestor.path).parent / self.path)

    @property
    def isotopicAbundancesByChemicalElement(self):
        '''Returns self.__isotopicAbundancesByChemicalElement.'''

        self.readIsotopicAbundancesByChemicalElement()
        return self.__isotopicAbundancesByChemicalElement

    def find(self, symbol, evaluation=None):
        '''
        Returns the entry matching a chemical elememnt's *symbol*. If non is found, None is returned.

        :param symbol:          The requested chemical element's PoPs symbol to match.
        :param evaluation:      The requested evaluation to match.

        :return:                Returns found chemical element or None.
        '''

        self.readIsotopicAbundancesByChemicalElement()

        for chemicalElement in self.__isotopicAbundancesByChemicalElement.chemicalElements:
            if self.adaptable(chemicalElement.symbol, symbol):
                 if self.adaptable(self.__evaluation, evaluation):
                    return chemicalElement

        return None

    def readIsotopicAbundancesByChemicalElement(self):
        '''Reads in the isotopicAbundancesByChemicalElement file pointed to by self if not already read in.'''

        if self.__isotopicAbundancesByChemicalElement is None:
            filePath = self.path
            if not pathlib.Path(filePath).is_absolute():
                filePath = pathlib.Path(self.derivedPath)
            self.__isotopicAbundancesByChemicalElement = isotopicAbundancesModule.read(filePath)
#            self.__isotopicAbundancesByChemicalElement.setAncestor(self)           # FIXME, is this needed?

        return self.__isotopicAbundancesByChemicalElement

    def toXML_strList(self, indent='', **kwargs):
        '''
        Returns a list of str instances representing the XML lines of self.

        :param indent:          The amount of indentation for each line. Child nodes and text may be indented more.
        :param kwargs:          A keyword list.

        :return:                List of str instances representing the XML lines of self.
        '''

        return ['%s<%s evaluation="%s"%s/>' % (indent, self.moniker, self.evaluation, self.standardXML_attributes())]

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        '''
        Creates a IsotopicAbundancesByChemicalElement instance from a isotopicAbundancesByChemicalElement node.

        :param node:            The node to parse.
        :param xPath:           Currently not used.
        :param linkData:        Currently not used.
        :param kwargs:          A keyword list.

        :return:                A iam's IsotopicAbundancesByChemicalElement instance.
        '''

        if node.tag != IsotopicAbundancesByChemicalElement.moniker:
            raise TypeError('Invalid node name.')

        kwargs = {attr: node.get(attr) for attr in ('evaluation', 'path', 'checksum', 'algorithm')}

        isotopicAbundancesByChemicalElement = IsotopicAbundancesByChemicalElement(**kwargs)

        return isotopicAbundancesByChemicalElement

    @staticmethod
    def adaptable(item1, item2):
        '''
        Returns True whenever item2 is None, or when item1 and item2 are the same. For internal use.

        :param item1:           Item to compare to item2.
        :param item2:           Item to compare to item1.
        '''

        if item2 is None:
            return True

        return item1 == item2

def read(fileName, **kwargs):
    '''
    Reads in the file name *fileName* and returns an **Ian** instance.

    :param fileName:            The name of the file to read.
    :param kwargs:              External dictionary of extra arguments that aids in parsing.
    '''

    return Iam.read(fileName, **kwargs)
