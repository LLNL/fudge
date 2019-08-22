Fission
=======

While fission is treated as a single reaction in both ENDF and GND, in reality it is composed of many different reactions each of which produces a different set of outgoing particles. Fission poses a special challenge for GND since, unlike most other reactions, we cannot explicitly list all the outgoing products. Instead we list average values for outgoing neutrons and gammas, and ignore the fission fragments (except to give the average energy deposited in them).

A neutron-induced fission channel generally contains several outgoing products::

- A prompt neutron
- One or more delayed neutrons, and
- Prompt gammas

In each case the multiplicity is generally a function of incident energy.

First, second, third chance, ...
--------------------------------

For most evaluations in major libraries, total fission (MT=18) contains a full description of the cross section,
product multiplicities and distributions. Additional fission cross sections like first-chance, 2nd-chance, etc.
are often also listed, but do not contain any description of outgoing products.
GND (and Fudge) handle these with the ``fissionComponent`` container, which contains a cross section
and Q-value but not necessarily a list of products.

The cross section
-----------------

Fission cross sections behave like all other reaction cross sections. For non-threshold fissioners,
part of the cross section is typically described by resonance parameters, so the cross section requires
reconstruction.

nubar
-----

The average neutron multiplicity 'nubar' is stored as the multiplicity of outgoing fission neutrons.
Nubar is typically separated into two components: 'prompt', for neutrons produced during the fission process,
and 'delayed' for beta-delayed neutrons emitted by the fission fragments.
Evaluations typically contain multiple delayed fission neutrons, grouped by the average beta-decay lifetime of
the parent fission products.  The ENDF-6 format only supports storing a single delayed nubar, so in GND
we typically see multiple delayed neutrons with multiplicities that link back to the first delayed neutron using
the ``reference`` multiplicity class.

PFNS
----

The prompt fission neutron spectrum (PFNS) is stored within the ``distribution`` information for the prompt neutron
product. All existing evaluations make the assumption that the outgoing angles and energies are uncoupled (or only
weakly coupled), and can be represented using ``uncorrelated`` energy and angle distributions.
The angular distribution is typically ``isotropic`` in the laboratory frame.
The outgoing energy distribution can be stored in several different ways:

- XYs2d  (i.e. a point-wise P(E' | E) distribution)
- Watt spectrum
- Simple Maxwellian fission spectrum
- Madland-Nix model

Some evaluations (notably Th232 in ENDF-VII.1 and ENDF-VIII-beta4) use a double-differential energy-angle distribution
for the PFNS, with a Legendre expansion given for each combination of incident and outgoing energy.

Gammas
------

TBD: expand documentation

Energy release
--------------

Delayed neutrons
----------------

Fission product yields
----------------------
