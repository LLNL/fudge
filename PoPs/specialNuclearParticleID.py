# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
The modules help convert a special nuclear particle id (i.e., name) between its familiar id, its nucleus id and its 
nuclide id.  The familiar ids are "p", "d", "t", "h" and "a". For example, the familiar id for the nucleus of a "He4" 
atom is "a" (for alpha particle). The nucleus id for "a" is "he4" and the nuclide id is "He4".
"""

from LUPY import enums as enumsModule

from . import IDs as IDsModule
from .families import nuclide as nuclideModule
from .families import nucleus as nucleusModule

class Mode(enumsModule.Enum):

    familiar = enumsModule.auto()
    nucleus = enumsModule.auto()
    nuclide = enumsModule.auto()

lightParticleIDs = [{IDsModule.familiarProton:   Mode.familiar, 'h1':  Mode.nucleus, 'H1':  Mode.nuclide},
                    {IDsModule.familiarDeuteron: Mode.familiar, 'h2':  Mode.nucleus, 'H2':  Mode.nuclide},
                    {IDsModule.familiarTriton:   Mode.familiar, 'h3':  Mode.nucleus, 'H3':  Mode.nuclide},
                    {IDsModule.familiarHelion:   Mode.familiar, 'he3': Mode.nucleus, 'He3': Mode.nuclide},
                    {IDsModule.familiarAlpha:    Mode.familiar, 'he4': Mode.nucleus, 'He4': Mode.nuclide}]

orderedLightParticles = (IDsModule.photon, IDsModule.neutron, IDsModule.proton, IDsModule.familiarDeuteron, 
                         IDsModule.familiarTriton, IDsModule.familiarHelion, IDsModule.familiarAlpha)

def sortLightParticle(ids):
    """
    From a list of particle *ids*, this function returns the ids sorted. The sorting order is 'n', 'p', 'd', 't', 'h', 'a', 'g', 
    (i.e., 'photon'). Other ids in the list are return after the ones mentioned and in the order they appear in the list.
    Note, an id will only appear once in the list even it appears more than once in *ids*.

    :param ids:     List of particle ids to sort.

    :return:        The sorted list of ids.
    """

    mapping = {}
    for id1 in set(ids): mapping[specialNuclearParticleID(id1, mode=Mode.familiar)] = id1
    photon = []
    if IDsModule.familiarPhoton in mapping: photon.append(mapping.pop(IDsModule.familiarPhoton))
    ids = list(mapping.keys())
    orderedIDs = []
    for orderedLightParticle in orderedLightParticles:
        if orderedLightParticle in ids:
            orderedIDs.append(orderedLightParticle)
            ids.pop(ids.index(orderedLightParticle))
    if len(photon) > 0: mapping[photon[0]] = photon[0]

    return [mapping[id1] for id1 in orderedIDs + photon + ids]

def findID(a_id):
    """
    For internal use. Returns the dictionary from *lightParticleIDs* that contain *a_id* as a key. 
    If no match is found **None** is returned.

    :param a_id:    The GNDS PoPs id of the particle to look up in the *lightParticleIDs* data set.

    :return:        The dictionary from *lightParticleIDs* that contain *a_id* as a key or **None**.
    """

    for lightPartcleID in lightParticleIDs:
        if a_id in lightPartcleID: return lightPartcleID

    return None
    
def familiarID(a_id):
    """
    Returns the "familiar" id for *a_id*. If no match is found, **None** is returned. Mainly for internal use.

    :param a_id:    The GNDS PoPs id of the particle to look up in the *lightParticleIDs* data set.

    :return:        Returns the familiar id for *a_id* or None if *a_id* does not have a familiar id.
    """

    lightPartcleID = findID(a_id)
    if lightPartcleID is None: return None
    for pid in lightPartcleID:
        if lightPartcleID[pid] == Mode.familiar: return pid

def nucleusID(a_id):
    """Returns the "nucleus" id for *a_id*. If no match is found, **None** is returned. Mainly for internal use.

    :param a_id:    The GNDS PoPs id of the particle to look up in the *lightParticleIDs* data set.

    :return:        Returns the nucleus id for *a_id* or None if *a_id* is not a familiar particle.
    """

    lightPartcleID = findID(a_id)
    if lightPartcleID is None: return None
    for pid in lightPartcleID:
        if lightPartcleID[pid] == Mode.nucleus: return pid

def nuclideID(a_id):
    """Returns the "nuclide" id for *a_id*. If no match is found, **None** is returned. Mainly for internal use.

    :param a_id:    The GNDS PoPs id of the particle to look up in the *lightParticleIDs* data set.

    :return:        Returns the nuclide id for *a_id* or None if *a_id* is not a familiar particle.
    """

    lightPartcleID = findID(a_id)
    if lightPartcleID is None: return None
    for pid in lightPartcleID:
        if lightPartcleID[pid] == Mode.nuclide: return pid

def specialNuclearParticleID(a_id, mode=Mode.familiar, pops=None):
    """
    Returns the special nuclear particle id per the specified *a_id* and *mode*. Note, if *pops* is not **None** 
    and *a_id* is the id for a nuclide or nucleus, then the id specified by *mode* is returned. For example, if
    *a_id* is "O16" and *mode* is **Mode.nucleus**, the "o16" is returned.

    :param a_id:    The GNDS PoPs id of the particle to look up in the *lightParticleIDs* data set.
    :param mode:    The **Mode** enum which specifying the special nuclear particle id to return for *a_id*.
    :param pops:    If not None, and *mode* is not **Mode.familiar**, then the *a_id* nuclide or nucleus id is return.

    :return:        Returns the id for *a_id* that match *mode*.
    """

    if type(mode) is str: mode = Mode(mode)

    if a_id == IDsModule.photon or a_id == IDsModule.familiarPhoton:
        if mode == Mode.familiar: return IDsModule.familiarPhoton
        return IDsModule.photon

    if   mode == Mode.familiar:
        pid = familiarID(a_id)
    elif mode == Mode.nucleus:
        pid = nucleusID(a_id)
    elif mode == Mode.nuclide:
        pid = nuclideID(a_id)

    if pid is None and mode != Mode.familiar:
        if pops is not None:
            try:
                particle = pops[a_id]
                if isinstance(particle, nuclideModule.Particle) and mode == Mode.nucleus:
                     pid = particle.nucleus.id
                if isinstance(particle, nucleusModule.Particle) and mode == Mode.nuclide:
                     pid = particle.nuclide.id
            except:
                pass
    if pid is None: pid = a_id

    return pid

def sameSpecialNuclearParticle(id1, id2):
    """
    Returns **True** is *id1* specified the same particle as *id2*. The function **specialNuclearParticleID** is called with **Mode.familiar**
    for path partcile.
    """

    return specialNuclearParticleID(id1, Mode.familiar) == specialNuclearParticleID(id2, Mode.familiar)
