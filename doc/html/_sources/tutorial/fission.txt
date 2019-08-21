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

The cross section
-----------------

PFNS
----

nubar
-----

Gammas
------

Energy release
--------------

Delayed neutrons
----------------

Fission product yields
----------------------
