.. _groups:

Particle groups
===============

In addition to families, PoPs also organizes some particles into 'groups'.
Particle groups are collections of particles that share some common properties.
For example, all isotopes of a chemical element share the same proton number Z,
and are organized together inside PoPs.  Two types of particle groups are defined:

* chemicalElement
    Defines the proton number Z, along with the chemical name and symbol.  A chemicalElement also contains a list
    of isotopes

* isotope
    Inherits the proton number Z from its parent chemicalElement, and defines the total nucleon number A.
    The isotope contains a list of nuclides (corresponding to the ground state and excited states).

Particle groups are not particles, and do not generally define particle properties like mass, spin or charge.
They are not assigned ids in PoPs, but they do have unique symbols that can be used to look them up from a PoPs
database. For example, the unique symbol for the chemicalElement Iron is 'Fe'.

chemicalElements
----------------

.. automodule:: PoPs.chemicalElements.chemicalElements
    :members:
    :special-members:

chemicalElement
---------------

.. automodule:: PoPs.chemicalElements.chemicalElement
    :members:
    :special-members:

isotope
-------

.. automodule:: PoPs.chemicalElements.isotope
    :members:
    :special-members:

misc
----

.. automodule:: PoPs.chemicalElements.misc
    :members:
    :special-members:
