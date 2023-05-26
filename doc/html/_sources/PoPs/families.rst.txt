.. _families:

Particle families
=================

Within a PoPs database, particles are organized into families.
Several particle families should be familiar to nuclear and particle physicists:

* gaugeBoson
    contains the photon, gluon, W and Z
* lepton
    contains electron, neutrino, etc.
* baryon
    contains the neutron and proton
* nucleus
    for nuclei in either the ground state or an excited state. A nucleus has a charge equal to Z (the proton number).
* nuclide
    particles made up of a nucleus plus Z electrons. A nuclide has charge = 0 since the number of protons and
    electrons balance out.

An additional particle family is defined for less-familiar types of particle:

* unorthodox
    contains particle-like objects. Examples include representative fission fragments (i.e. average particles
    meant to represent a range of isotopes produced by fission) or 'thermal scattering particles'
    representing the collective behavior of atoms inside a molecule.
    Unorthodox particles are sometimes used to support modeling nuclear data, and
    like normal particles they may be assigned basic particle properties like a mass, charge, etc.

Each particle family inherits from an abstract base particle class, defined in PoPs.families.particle.
A family may define additional particle properties that its members can contain.


Base particle
-------------

.. automodule:: PoPs.families.particle
    :members:
    :special-members:

gaugeBoson
----------

.. automodule:: PoPs.families.gaugeBoson
    :members:
    :special-members:
    :exclude-members: alias, suite

lepton
------

.. automodule:: PoPs.families.lepton
    :members:
    :special-members:
    :exclude-members: alias, suite

baryon
------

.. automodule:: PoPs.families.baryon
    :members:
    :special-members:
    :exclude-members: alias, suite

nuclide
-------

.. automodule:: PoPs.families.nuclide
    :members:
    :special-members:
    :exclude-members: alias, suite

nucleus
-------

.. automodule:: PoPs.families.nucleus
    :members:
    :special-members:
    :exclude-members: alias, suite

unorthodox
----------

.. automodule:: PoPs.families.unorthodox
    :members:
    :special-members:
    :exclude-members: alias, suite
