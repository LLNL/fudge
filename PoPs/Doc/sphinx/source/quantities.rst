.. _quantities:

Particle properties
===================

PoPs includes several classes for storing different types of physical quantities.
These are all organized inside the PoPs.quantities module.

Each particle property container derives from the base class called 'quantity'.
Quantities may be stored as ints, doubles, fractions and strings.

PoPs supports storing more than one value for each quantity, for example when an
evaluator cannot make a firm spin/parity assignment to an excited state but can
make a list of possible assignments. Thus each quantity (like 'mass', 'spin', etc.) is
actually a suite of possible assignments. The first item in the suite should be the
recommended value, as it will be used as the default unless another assignment is explicitly chosen.

quantity
--------

.. automodule:: PoPs.quantities.quantity
    :members:
    :special-members:

mass
----

.. automodule:: PoPs.quantities.mass
    :members:
    :special-members:

charge
------

.. automodule:: PoPs.quantities.charge
    :members:
    :special-members:

spin
----

.. automodule:: PoPs.quantities.spin
    :members:
    :special-members:

parity
------

.. automodule:: PoPs.quantities.parity
    :members:
    :special-members:

halflife
--------

.. automodule:: PoPs.quantities.halflife
    :members:
    :special-members:

nuclearEnergyLevel
------------------

.. automodule:: PoPs.quantities.nuclearEnergyLevel
    :members:
    :special-members:
