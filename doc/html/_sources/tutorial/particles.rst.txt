Particles
=========

Nuclei and their properties
---------------------------

GNDS evaluations store a list of all particles involved in each reaction, including the projectile and target,
outgoing products and residuals.  Some properties of reaction products (like the multiplicity and distribution data)
depend on the reaction that produced them.  Other properties (like the mass, spin and parity, charge, etc.) are
reaction-independent.  These properties are collected in the Properties of Particles (``PoPs``) section of the
``reactionSuite``.

TBD: expand description of PoPs (it is still changing rapidly to accommodate the ENDF decay sub-library).

Transportable particles
-----------------------

When GNDS files are processed, the ``styles`` section contains a list of the particles that can be transported
with this data. Current evaluations typically only support transporting neutrons, photons and light charged particles
(isotopes of Hydrogen and Helium), but GNDS places no limitation on what particles may be transported.
However, evaluators must ensure that outgoing product multiplicities and distributions are listed for all
transportable particles so that they are properly tracked in transport applications.
