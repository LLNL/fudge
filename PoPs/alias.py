# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the alias classes.
"""

import abc

from fudge import GNDS_formatVersion as GNDS_formatVersionModule
from . import misc as miscModule
from . import suite as suiteModule

from .chemicalElements import misc as chemicalElementMiscModule

class BaseAlias(miscModule.ClassWithIDKey, abc.ABC):

    def __init__(self, id, pid):

        miscModule.ClassWithIDKey.__init__(self, id)

        if not isinstance(pid, str): raise TypeError('pid not str')
        self.__pid = pid

    @property
    def pid(self):

        return self.__pid

    def check(self, info):
        from . import warning as warningModule

        warnings = []
        if self.pid not in info['PoPs']:
            warnings.append(warningModule.AliasToNonExistentParticle(self.id, self.pid, self))
        return warnings

    def copy(self):

        return self.__class__(self.id, self.pid)

    def isAlias(self):

        return True

    def isMetaStable(self):

        return isinstance(self, MetaStable)

    def toXML_strList(self, indent='', **kwargs):

        formatVersion = kwargs.get('formatVersion', GNDS_formatVersionModule.default)
        name = self.moniker
        if formatVersion == GNDS_formatVersionModule.version_1_10: name = 'particle'

        return ['%s<%s id="%s" pid="%s"/>' % (indent, name, self.id, self.pid)]

class Alias(BaseAlias):
    moniker = 'alias'

    def intid(self, intidDB={}):
        '''
        Converts the particle id into a unique integer dubbed an INTeger ID (INTID).
        '''

        from . import database as databaseModule

        if self.id in ['d', 't', 'h', 'a']:
            return self.findClassInAncestry(databaseModule.Database)[self.pid].intid()

        raise ValueError('Alias "%s" for "%s" does not have a defined intid.' % (self.id, self.pid))

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        xPath.append(node.tag)

        self = cls(node.get('id'), node.get('pid'))

        xPath.pop()

        return self

class MetaStable(BaseAlias):
    moniker = 'metaStable'

    def __init__(self, id, pid, metaStableIndex):

        BaseAlias.__init__(self, id, pid)

        if not isinstance(metaStableIndex, int):
            raise TypeError('metaStableIndex must be an int: %s' % miscModule.toLimitedString(metaStableIndex))
        self.__metaStableIndex = metaStableIndex

    @property
    def metaStableIndex(self):

        return self.__metaStableIndex

    def toXML_strList(self, indent='', **kwargs):

        return (['%s<%s id="%s" pid="%s" metaStableIndex="%s"/>' % (
            indent, self.moniker, self.id, self.pid, self.metaStableIndex)])

    def copy(self):

        return self.__class__(self.id, self.pid, self.metaStableIndex)

    def intid(self, intidDB={}):
        '''
        Converts the particle id into a unique integer dubbed an INTeger ID (INTID).
        '''

        from . import database as databaseModule

        nuclide = self.findClassInAncestry(databaseModule.Database)[self.pid]
        intid1 = nuclide.intid()
        sign = -1 if intid1 < 0 else 1

        return sign * (1000000 * (self.metaStableIndex + 480) + abs(intid1) % 1000000)

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)

        instance = cls(element.get('id'), element.get('pid'), int(element.get('metaStableIndex')))

        xPath.pop()
        return instance

    @staticmethod
    def metaStableNameFromNuclearLevelNameAndMetaStableIndex(nuclideName, metaStableIndex):

        if not isinstance(metaStableIndex, int):
            raise TypeError('metaStableIndex must be an int: %s' % miscModule.toLimitedString(metaStableIndex))
        if metaStableIndex < 1: raise ValueError('metaStableIndex must be greater than 0 got "%s".' % metaStableIndex)

        baseID, chemicalElementID, A, levelID, isNucleus, anti, qualifier = chemicalElementMiscModule.chemicalElementALevelIDsAndAnti(
            nuclideName)

        if isNucleus: chemicalElementID = chemicalElementID[0].lower() + chemicalElementID[1:]
        return "%s_m%d" % (
            chemicalElementMiscModule.isotopeSymbolFromChemicalElementIDAndA(chemicalElementID, A), metaStableIndex)

    @staticmethod
    def nuclideNameAndMetaStableIndexFromName(name):
        """
        This function returns the nuclide name and meta-stable index from its argument name. If name does not appear to be a 
        meta-stable, ( name, 0 ) are returned. This function splits on the string '_m'. If the number of sub-strings returned 
        is not 2, the name is considered not to be a meta-stable and ( name, 0 ) are returned. For example, name = 'O16' will 
        return ( 'O16', 0 ), name = 'Am242_m1' will return ( 'Am242_m1', 1 ) and name = 'Am242_m1_m2' will return ( 'Am242_m1_m2', 0 ).
        """

        if name.count('_m') != 1: return name, 0
        nuclideName, metaStableIndex = name.split('_m')
        try:
            metaStableIndex = int(metaStableIndex)
        except:
            metaStableIndex = 0

        return nuclideName, metaStableIndex

class Suite(suiteModule.Suite):
    moniker = 'aliases'

    def __init__(self):

        suiteModule.Suite.__init__(self, (BaseAlias,))

    def has_pid(self, ParticleID):
        """Returns True if one of the aliases has pid equal to ParticleID and False otherwise."""

        for alias in self:
            if alias.pid == ParticleID: return True

        return False

    def parseNode(self, element, xPath, linkData, **kwargs):

        for child in element:
            if child.tag == MetaStable.moniker:
                self.add(MetaStable.parseNodeUsingClass(child, xPath, linkData, **kwargs))
            else:
                self.add(Alias.parseNodeUsingClass(child, xPath, linkData, **kwargs))

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)

        instance = cls()
        instance.parseNode(element, xPath, linkData, **kwargs)

        xPath.pop()
        return instance
