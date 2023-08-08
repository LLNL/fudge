# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

'''
Defines the function idFromIntid and the class IntidDB.
'''

import pathlib

from . import IDs as IDsModule
from . import misc as miscModule
from . import specialNuclearParticleID as specialNuclearParticleIDModule
from PoPs.chemicalElements import misc as chemicalElementsMiscModule

def idFromIntid(intid, intidDB, mode=specialNuclearParticleIDModule.Mode.familiar):
    '''
    Returns the GNDS particle id from the particles intid (INTeger ID). The following list the
    intid values given a family/particle.

        |---------------+-------------------------------------------------------------------|
        | Bits          | Description                                                       |
        |===============+===================================================================|
        | 31            | Particle/anti-particle bit. 0 is particle, 1 is anti-particle.    |
        |               | This is the sign bit so a negative intid is an anti-particle.     |
        |---------------+-------------------------------------------------------------------|
        | 30            | Nuclide/nuclear-like bit. 0 if Nuclide/nuclear-like and           |
        |               | 1 if not.                                                         |
        |---------------+-------------------------------------------------------------------|

    For nuclide/nuclear-like particle, the other bits are:

         IZA = 1000000 * III + 1000 * ZZZ + AAA

    where III is the level index with the following meaning:

        *) For 0 <= III < 481 the particle is a nuclide and III is the nuclear level of its nucleus.
            For example, III = 0 is a nuclide with its nucleus in its nuclear ground state (e.g., "O16")
            and III = 1 is a nuclide with its nucleus in its first nuclear excited state (e.g., "O16_e1").

        *) For 481 <= III < 500 the particle is a meta-stable for a nuclide with meta-stable index
            (III - 480). For example, the IZA for "Am242_m1" is 481095242.

        *) For 500 <= III < 980 the particle is a nucleus in nuclear excited state (III - 500).
            For example, III = 500 is a nucleus in its ground state (e.g., "o16") and
            III = 501 is a nucleus in its first excited nuclear state (e.g., "o16_e1").

        *) For 981 <= III < 1000 the particle is a meta-stable for a nucleus with meta-stable index
            (III - 980). For example, the IZA for "am242_m1" is 981095242.

        |---------------+---------------------------------------------------------------|
        | Bits          | Description                                                   |
        |===============+===============================================================|
        | 26-29         | The particle's family when bit 31 is 1.                       |
        |---------------+---------------------------------------------------------------|

        |-------------------+-----------------------|
        | Bit 26-29 value   | Particle's family     |
        | when bit 30 is 1  |                       |
        |===================+=======================|
        | 0                 | Gauge boson           |
        |-------------------+-----------------------|
        | 1                 | Lepton                |
        |-------------------+-----------------------|
        | 2                 | Baryon with 3 quarks  |
        |-------------------+-----------------------|
        | 8                 | 99120 and 99125       |
        |-------------------+-----------------------|
        | 9                 | TNSL                  |
        |-------------------+-----------------------|

    Gauge boson bits 0 to 25.

        |-------------------+-----------------------|
        | Bits 0-25 value   | Particle              |
        |===================+=======================|
        | 0                 | photon                |
        |-------------------+-----------------------|
        | 1                 | gluon                 |
        |-------------------+-----------------------|
        | 2                 | W boson               |
        |-------------------+-----------------------|
        | 3                 | Z boson               |
        |-------------------+-----------------------|

    Lepton generation is designated by bits:

        |-------------------+-----------------------|
        | Bits 0-1 value    | Generation            |
        |===================+=======================|
        | 0                 | Electron              |
        |-------------------+-----------------------|
        | 1                 | Muon                  |
        |-------------------+-----------------------|
        | 2                 | Tau                   |
        |-------------------+-----------------------|

    Bit 2 designates if the lepton is a neutrino as:

        |-------------------+-----------------------|
        | Bit 2 value       | Is a neutrino?        |
        |===================+=======================|
        | 0                 | False                 |
        |-------------------+-----------------------|
        | 1                 | True                  |
        |-------------------+-----------------------|

    Baryon group is designated by bits 20 to 25.

        |-------------------+-----------------------|
        | Bits 20-25 value  | Group                 |
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

    Baryon nucleon group particle designated by bits 0 to 22.

        |-------------------+-----------------------|
        | Bits 0-19 value   | Particle              |
        |===================+=======================|
        | 0                 | Neutron (n)           |
        |-------------------+-----------------------|
        | 1                 | Proton (p)            |
        |-------------------+-----------------------|


    That is, the neutron has bits 0 to 25 being "0" or integer 0, and the proton is "1" or integer 1.

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


    :param intid:       The requested particle's intid.
    :param intidDB:     A intid database for particles that are not defined intid particles.
    :param mode:        Several particles have more than one possible id. This specified which of the possible ids to return.

    :return:            The particle's id.
    '''

    def getBits(intid, begin, end):
        '''
        Returns the value of *intid* between bits begin to end.
        '''

        return (intid >> begin) % 2**(end-begin)

    if intid in intidDB.intids:
        return intidDB.intids[intid]

    anti = ''
    if intid < 1:
        anti = miscModule.antiSuffix

    abs_intid = abs(intid)
    if getBits(abs_intid, 30, 31) == 0:
        AAA = abs_intid % 1000
        ZZZ = abs_intid % 1000000 // 1000
        III = abs_intid % 1000000000 // 1000000
        try:
            pid = chemicalElementsMiscModule.idFromZAndA(ZZZ, AAA)
        except:
            print(intid, III, ZZZ, AAA)
            raise
        if III >= 500:
            pid = pid.lower()
            III -= 500
        if III > 0:
            if III <= 480:
                pid = '%s_e%d' % (pid, III)
            else:
                III -= 480
                pid = '%s_m%d' % (pid, III)
    else:
        family = getBits(abs_intid, 26, 30)
        if family == 0:
            pid = IDsModule.photon
        elif family == 1:
            pid = IDsModule.electron
        elif family == 2:
            baryonGroup = getBits(abs_intid, 20, 26)
            if baryonGroup == 0:
                baryon = getBits(abs_intid, 0, 20)
                if baryon == 0:
                    pid = IDsModule.neutron
                elif baryon == 1:
                    pid = IDsModule.proton
                else:
                    raise Exception('Oops, unsupported baryon : %s' % baryon)
            else:
                raise Exception('Oops, unsupported baryon group: %s' % baryonGroup)
        else:
            raise Exception('Oops, unsupported family: %s' % family)

    pid = specialNuclearParticleIDModule.specialNuclearParticleID(pid, mode)

    return pid + anti

class IntidDB(dict):
    '''
    Stores ids and there associated intid values for particles that do not have a standard intid definition.
    This classes acts like a python **dict** object except that it has an inverse dictionary that maps intids
    to there associated ids. Also, this class has a read method for reading data from a file. Each line of the
    read file can::

        -) be blank (after its comment is removed),
        -) contain an id and its associated intid (e.g., "HinH2O 2751463432"),
        -) contain an import key followed by a file name to import (e.g., "import TNSL.intid".

    For each line, all characters at and after a "#" character are ignored.

    :param fileName:        If not **None**, the path to the file to read intid data from.
    '''

    def __init__(self, fileName=None):

        self.__intids = {}

        if fileName is not None:
            self.read(fileName)

    def read(self, fileName):
        '''
        Reads in the data in file *fileName*.

        :param fileName:    The path to the file to read in.
        '''

        file = pathlib.Path(fileName)
        with open(fileName) as fIn:
            lines = fIn.readlines()

        for line in lines:
            line = line.split('#')[0].strip()
            if len(line) == 0:
                continue
            item1, item2 = line.split()
            if item1 == 'import':
                self.read(pathlib.Path(file.parent) / item2)
            else:
                intid = int(item2)
                self[item1] = intid
                self.__intids[intid] = item1

    @property
    def intids(self):
        '''Returns a reference to the __intids member.'''

        return self.__intids
