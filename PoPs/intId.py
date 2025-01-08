# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

'''
Module for all particle intid (INTeger ID) classes and functions.

The intid is a signed integer that uniquely maps to a GNDS PoPs (i.e., particle) id. 
An intid of -1 implies an unknown particle. For example, currently, all TNSL ids have an intid of -1.
The sign of the intid indicates whether the particle is a particle or anti-particle.
The particle intid is then defined by its base-10 digits for a signed 32-bit integer (s32).
The following tables definite the meaning of an intid base-10 digits. For a s32 the 9th base-10 digit can
only be 0, 1, 2. Firstly, the digits of the intid can be written as::

         base-10 digit     9876543210
         definition       PN?????????

where the P is ' ' or '+' for a particle and '-' for an anti-particle, and N indicates if the particle is 
a nuclear type particle or not.
If the value is 0, then the particle is a nuclear type particle (i.e., a nuclide or nucleus), and
if the value is 1, the particle is not a nuclear type particle; ergo, anything else. Currently N = 2 is not defined.

For a nuclear type particle, the base-10 digits are defined as::

         base-10 digit     9876543210
         definition       P0IIIZZZAAA

That is, the intid can be written as::

     intid = P ( 1000000 * III + 1000 * ZZZ + AAA )

where ZZZ is the number of protons of the particle (also called an atomic number), AAA is the number of protons 
plus neutrons (also called the atomic mass number) and III is the level index with the following meaning:

    *) For 0 <= III < 500 the particle is a nuclide and III is the nuclear level of its nucleus.
        For example, III = 0 is a nuclide with its nucleus in its nuclear ground state (e.g., "O16")
        and III = 1 is a nuclide with its nucleus in its first nuclear excited state (e.g., "O16_e1").

    *) For 500 <= III < 1000 the particle is a nucleus in nuclear excited state (III - 500).
        For example, III = 500 is a nucleus in its ground state (e.g., "o16") and
        III = 501 is a nucleus in its first excited nuclear state (e.g., "o16_e1").

For non-nuclear paricles, the base-10 digits are defined as::

         base-10 digit     9876543210
         definition       P1FFSSSSSSS

where P is as before, FF is the particle's family and SSSSSSS is the particle's indentifier within 
its family. The following families are currently defined (note a space is equivalent to the digit
0, for example 1 and 01 are equivalent):

    |-------------------+-----------------------|
    | FF value          | Particle's family     |
    |===================+=======================|
    |  0                | Gauge boson           |
    |-------------------+-----------------------|
    |  1                | Lepton                |
    |-------------------+-----------------------|
    |  2                | Baryon with 3 quarks  |
    |-------------------+-----------------------|
    | 50 - 59           | Nuclide meta-stables. |
    |-------------------+-----------------------|
    | 60 - 69           | Nuclear meta-stables. |
    |-------------------+-----------------------|
    | 98                | TNSL: this is         |
    |                   | currently not neeeded |
    |                   | hence not used.       |
    |-------------------+-----------------------|
    | 99                | ENDL fission product  |
    |-------------------+-----------------------|

For gauge bosons the SSSSSSS values are defined as:

    |-------------------+-----------------------|
    | SSSSSSS           | Particle              |
    |===================+=======================|
    |       0           | photon                |
    |-------------------+-----------------------|
    |       1           | gluon                 |
    |-------------------+-----------------------|
    |       2           | W boson               |
    |-------------------+-----------------------|
    |       3           | Z boson               |
    |-------------------+-----------------------|

For leptons the generation is designated by the right most S value as (the table only
defines the values under the capital S columns, not the ones under the small s column:

    |-------------------+-----------------------|
    | ssssssS           | Generation            |
    |===================+=======================|
    |       0           | Electron              |
    |-------------------+-----------------------|
    |       1           | Muon                  |
    |-------------------+-----------------------|
    |       2           | Tau                   |
    |-------------------+-----------------------|

The second S value designates if the lepton is a neutrino as:

    |-------------------+-----------------------|
    | sssssSs           | Is a neutrino?        |
    |===================+=======================|
    |      0            | False                 |
    |-------------------+-----------------------|
    |      1            | True                  |
    |-------------------+-----------------------|

The baryons are grouped with the group designated by the upper S value as:

    |-------------------+-----------------------|
    | Sssssss           | Group                 |
    |===================+=======================|
    | 0                 | Nucleon               |
    |-------------------+-----------------------|
    | 1                 | Delta                 |
    |-------------------+-----------------------|
    | 2                 | Lambda                |
    |-------------------+-----------------------|
    | 3                 | Sigma                 |
    |-------------------+-----------------------|
    | 4                 | Xi                    |
    |-------------------+-----------------------|
    | 5                 | Omega                 |
    |-------------------+-----------------------|

The particles in the nucleon group of the baryon family are designated by the lower 6 S values as:

    |-------------------+-----------------------|
    | sSSSSSS           | Particle              |
    |===================+=======================|
    |       0           | Neutron (n)           |
    |-------------------+-----------------------|
    |       1           | Proton (p)            |
    |-------------------+-----------------------|

Note, in addition to a neutron (i.e., id = 'n') having the intid of 1020000000 (and 'n_anti' of -1020000000), the 
parseIntid and idFromIntid functions also returns 'n' ('n_anti') when an intid of 1 (-1) is entered due to the way 
PoPs in FUDGE is implemented. However, this should not be considered as an alias intid but an idiosyncrasy of FUDGE.

For nuclide meta-stables the coding is::

         base-10 digit     9876543210
         definition       P15MMZZZAAA

where ZZZ and AAA are as for nuclear type particles, and MM is the meta-stable index. For example, the first and second
meta-stables for Zn61 are::

    |----------+----------------|
    | GNDS id  |  intid         |
    |==========+================|
    | Zn61_m1  |  1501030061    |
    | Zn61_m2  |  1502030061    |
    |----------+----------------|

Note, MM = 00 is not a valid meta-stable index.

For nucleus meta-stables the coding is the same as for nuclide meta-stables except the "5" in "P15" is replaced with a "6" as "P16".

The SSSSSSS value for TNSL are currently not needed or used, and are thus not defined.
The following aliases are mapped to other ids:

    |-------------------+-----------------------|
    | Alias id          | Id used for intid     |
    |===================+=======================|
    | d                 | h2                    |
    |-------------------+-----------------------|
    | t                 | h3                    |
    |-------------------+-----------------------|
    | h                 | he3                   |
    |-------------------+-----------------------|
    | a                 | he4                   |
    |-------------------+-----------------------|

Bit SSSSSSS values for ENDL fission products 99120 and 99125 are designated as:

    |-------------------+---------------------------|
    | SSSSSSS           | Fission product           |
    |-------------------+---------------------------|
    |   99120           | FissionProductENDL99120   |
    |-------------------+---------------------------|
    |   99125           | FissionProductENDL99125   |
    |-------------------+---------------------------|

Questions/comments:

    -) Is 2095242 (Am242_e2) this the same as 481095242 (Am242_m1)?
'''

import pathlib

from LUPY import enums as enumsModule

from . import IDs as IDsModule
from . import misc as miscModule
from . import specialNuclearParticleID as specialNuclearParticleIDModule
from .chemicalElements import misc as chemicalElementsMiscModule

class Family(enumsModule.Enum):
    '''
    Enum for the allowed intid families.
    '''

    nuclide = enumsModule.auto()
    nucleus = enumsModule.auto() 
    gaugeBoson = enumsModule.auto()
    lepton = enumsModule.auto()
    baryon = enumsModule.auto()
    nuclideMetaStable = enumsModule.auto()
    nucleusMetaStable = enumsModule.auto()
    TNSL = enumsModule.auto()
    ENDL_fissionProduct = enumsModule.auto()

def family2Integer(family):
    '''
    Returns an integer representing the *family* value. Family.nucleus and Family.nuclide do not have a
    family integer value and hence return -1.

    :param family:      Family enum value.

    :return:            Integer.
    '''

    return {Family.nucleus: -1, Family.nuclide: -2, Family.gaugeBoson: 0, Family.lepton: 1, Family.baryon: 2, 
            Family.nuclideMetaStable: 50, Family.nucleusMetaStable: 60, Family.TNSL: 98, Family.ENDL_fissionProduct: 99}[family]

def parseIntid(intid):
    '''
    Returns a tuple of the parts of intid. The number items returned depends on the family. Each return value will have
    the family enum, a boolean specifying if the particle is an anti particle. The number of values returned are family
    dependent. The following tables list the values returned by family:

    |-----------------------+-----------------------------------------------|
    | family                | values returned                               |
    |=======================+===============================================|
    | nuclide               | family enum, isAnti, III, ZZZ, AAA            |
    |-----------------------+-----------------------------------------------|
    | nucleus               | family enum, isAnti, III, ZZZ, AAA            |
    |-----------------------+-----------------------------------------------|
    | gauge boson           | family enum, isAnti, particle id              |
    |-----------------------------------------------------------------------|
    | lepton                | family enum, isAnti, generation, isNeutrino   |
    |-----------------------+-----------------------------------------------|
    | baryon                | family enum, isAnti, group, particle id       |
    |-----------------------+-----------------------------------------------|
    | ENDL fission product  | family enum, isAnti, particle id              |
    |-----------------------+-----------------------------------------------|
    | TNSL                  | family enum, isAnti, particle id              |
    |-----------------------+-----------------------------------------------|

    If intid is not a valid particle, then **None** is returned.

    :param intid:   Particle's intid.

    :return:        Python tuple whose length depends on the family.
    '''

    isAnti = intid < 0

    abs_intid = abs(intid)

    nuclearLike = abs_intid // 10**9                    # Digit can to 0 (nuclear type) or 1 (other).
    family = (abs_intid // 10**7) % 100                 # Only used for non-nuclear type.
    SSSSSSS = abs_intid % 10**7                         # Only used for non-nuclear type.

    if nuclearLike == 0:
        AAA = abs_intid % 1000
        ZZZ = abs_intid % 1000000 // 1000
        III = ( abs_intid % 1000000000 // 1000000 )
        family = Family.nuclide
        if III >= 500:
            family = Family.nucleus
        return family, isAnti, III, ZZZ, AAA

    if nuclearLike != 1:
        raise ValueError('Invalid intid = "%s".' % intid)

    if nuclearLike == 1:
        if family == 0:                                 # Gauge boson.
            return Family.gaugeBoson, isAnti, SSSSSSS
        elif family == 1:                               # Lepton.
            generation = SSSSSSS % 10
            isNeutrino = (SSSSSSS % 100) // 10 != 0
            return Family.lepton, isAnti, generation, isNeutrino
        elif family == 2:                               # Baryon with 3 quarks.
            group = SSSSSSS // 10**6
            baryon = SSSSSSS % 10**6
            return Family.baryon, isAnti, group, baryon
        elif 50 <= family < 70:                         # Nuclide or nucleus meta-stable.
            AAA = abs_intid % 1000
            ZZZ = abs_intid % 1000000 // 1000
            III = abs_intid % 1000000000 // 1000000
            MM = III % 100
            family = Family.nuclideMetaStable
            if III // 100 == 6:
                family = Family.nucleusMetaStable
            return family, isAnti, MM, ZZZ, AAA
        elif family == 98:                               # TNSL target.
            return Family.TNSL, isAnti, SSSSSSS
        elif family == 99:                               # ENDL fission product.
            return Family.ENDL_fissionProduct, isAnti, SSSSSSS

    return None

def idFromIntid(intid, intidDB=None, mode=specialNuclearParticleIDModule.Mode.familiar):
    '''
    Returns the GNDS particle id from the particle's intid (INTeger ID).

    :param intid:       The requested particle's intid.
    :param intidDB:     A intid database for particles that are not defined intid particles.
    :param mode:        Several particles have more than one possible id. This specified which of the possible ids to return.

    :return:            The particle's id.
    '''

    if intid == -1:
        return ''

    if intidDB is not None:
        if intid in intidDB.intids:
            return intidDB.intids[intid]

    family, isAnti, *others = parseIntid(intid)
    anti = miscModule.antiSuffix if isAnti else ''

    if family in (Family.nuclide, Family.nucleus, Family.nuclideMetaStable, Family.nucleusMetaStable):   # Nuclide/nuclear-like particle or thier meta-stable alias.
        III, ZZZ, AAA = others
        pid = chemicalElementsMiscModule.idFromZAndA(ZZZ, AAA, False)
        if family == Family.nucleus:
            III -= 500
        if family in (Family.nucleus, Family.nucleusMetaStable):
            pid = pid.lower()
        if III > 0:
            if family in (Family.nuclide, Family.nucleus):
                pid = '%s_e%d' % (pid, III)
            else:
                pid = '%s_m%d' % (pid, III)
    elif family == Family.gaugeBoson:                                   # Gauge boson.
        particleID = others[0]
        if particleID != 0:
            raise Exception('Oops, unsupported gauge boson: %s' % particleID)
        pid = IDsModule.photon
    elif family == Family.lepton:                                       # Lepton.
        generation, isNeutrino = others
        if generation != 0:
            raise Exception('Oops, unsupported lepton generation: %s' % generation)
        if isNeutrino:
            raise Exception('Oops, unsupported lepton: neutrino')
        pid = IDsModule.electron
    elif family == Family.baryon:                                       # Baryon with 3 quarks
        baryonGroup, baryon = others
        if baryonGroup == 0:
            if baryon == 0:
                pid = IDsModule.neutron
            elif baryon == 1:
                pid = IDsModule.proton
            else:
                raise Exception('Oops, unsupported baryon: %s' % baryon)
        else:
            raise Exception('Oops, unsupported baryon group: %s' % baryonGroup)
    elif family == Family.ENDL_fissionProduct:                          # ENDL fission product.
        particleID = others[0]
        if particleID == 99120:
            pid = IDsModule.FissionProductENDL99120
        elif particleID == 99125:
            pid = IDsModule.FissionProductENDL99125
        else:
            raise Exception('Oops, unsupported ENDL fission product: %s' % particleID)
    else:
        raise Exception('Oops, unsupported family: %s' % family)

    pid = specialNuclearParticleIDModule.specialNuclearParticleID(pid, mode)

    return pid + anti

class IntidDB:
    '''
    This class is currently not needed.
    Stores ids and there associated intid values for particles that do not have a standard intid definitions.
    This class acts like a python **dict** object except that it has an inverse dictionary that maps intids
    to there associated ids. Also, this class has a read method for reading data from a file. Each line of the
    read file can::

        -) the first line must contain a "format" key followed by a version number of the form "I.I" 
            (without the double quotes) where "I" is an integer (e.g., "1.0" without the double quotes).
        -) be blank (after its comment is removed),
        -) contain an id and its associated intid (e.g., "HinH2O 2751463432" without the double quotes),
        -) contain an import key followed by a file name to import (e.g., "import TNSL.intid").

    For each line, the comment character "#" and all characters after it are ignored.
    '''

    def __init__(self, path=None):
        '''
        :param path:            If not **None**, the path to the file to read intid data from.
        '''

        self.__ids = {}
        self.__intids = {}

        if path is not None:
            self.read(path)

    def __len__(self):
        '''
        Returns the number of entries (i.e., ids) in *self*.
        '''

        return len(self.__ids)

    def __getitem__(self, id):
        '''
        Returns the intid for key *id*.

        :param id:      The name of the id whose intid is returned.
        '''

        return self.__ids[id]

    def __iter__(self):
        '''
        Iterates over the entries (i.e., ids) in *self*.
        '''

        for id in self.__ids:
            yield id

    def __contains__(self, id):
        '''
        Returns **True** if *id* in *self* and **False** otherwise.

        :param id:      The name of the id whose intid is returned.
        '''

        return id in self.__ids

    def add(self, id, intid):
        '''
        Adds key *id* to dictionary with value *intid*. If *id* or *intid* already exist, a raise is executed.

        :param id:      The name of the key to add.
        :param value:   The value for key *id*.
        '''

        if id in self.__ids:
            raise KeyError('Key "%s" already exist.' % id)
        if intid in self.__intids:
            raise KeyError('intid "%s" already exist.' % intid)

        if not isinstance(id, str):
            raise TypeError('id must be a str instance.')
        if not isinstance(intid, int):
            raise TypeError('intid must be an int instance.')

        self.__ids[id] = intid
        self.__intids[intid] = id

    def get(self, id, default=None):
        '''
        Return the intid for *id* if *id* is in the dictionary, else default.

        :param id:          The *id* of the particle to lookup.
        :default:           Value to return if *id* is not in *self*.
        '''

        return self.__ids.get(id, default)

    def pop(self, id):
        '''
        Removes *id* from *self* and returns its value (i.e., intid).

        :param id:      id to remove from self.

        :return:        The intid value of the removed *id* removed.
        '''

        if id not in self:
            raise KeyError('id "%s" not in *self*.' % id)

        intid = self.__ids[id]
        del self.__ids[id]
        del self.__intids[intid]

        return intid

    def read(self, path):
        '''
        Reads in id/intid data from file *path*.

        :param path:    The path to the file to read in.
        '''

        path = pathlib.Path(path)
        with open(path) as fIn:
            lines = fIn.readlines()

        key, value = lines.pop(0).split('#')[0].strip().split()
        if key != 'format':
            raise KeyError('First line must specify format version.')
        if value != '1.0':
            raise ValueError('Unsupport intid format version "%s".' % value )

        for line in lines:
            line = line.split('#')[0].strip()
            if len(line) == 0:
                continue
            key, value =  line.split()
            if key == 'import':
                importPath = pathlib.Path(value)
                if importPath.is_absolute():
                    self.read(importPath)
                else:
                    self.read(pathlib.Path(path.parent) / importPath)
            else:
                self.add(key, int(value))

    def saveToFile(self, path):
        '''
        Writes the contents of *self* to file *path*.

        :param path:    Name of file to write data to.
        '''

        fmt = '%%-%ds' % max(map(len, self.__ids))
        with open(path, 'w') as fOut:
            print('format 1.0', file=fOut)
            for intid in sorted(self.__intids):
                print(fmt % self.__intids[intid], intid, file=fOut)

    @property
    def intids(self):
        '''Returns a reference to the __intids member of *self*.'''

        return self.__intids

def intidHelper(isAnti, family, SSSSSSS):
    '''
    Returns intid for requeted *family* and its *SSSSSSS* value.

    :param isAnti:          **True** is particle is an anti-particle and **False** otherwise.
    :param family:          FUDGE PoPs family or negative integer for special cases (e.g., 'FissionProductENDL9912[05]').
    :param SSSSSSS:         The particle's indentifier within its family.

    :return:                Integer with bits 24-31 of a particle's intid set.
    '''

    if isinstance(isAnti, str):
        isAnti = {'': False, miscModule.antiSuffix: True}[isAnti]
    sign = -1 if isAnti else 1
    

    intid = family2Integer(family)
    if intid < 0:
        raise ValueType('Family "%s" not supported.')
    intid *= 10**7
    intid += 10**9

    return sign * (intid + SSSSSSS)
