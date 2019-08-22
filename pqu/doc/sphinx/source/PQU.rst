.. automodule:: pqu.PQU

Variables
---------

This sections list some of the "static" variables in the PQU module.

maxSignificantDigits
++++++++++++++++++++
.. py:data:: pqu.PQU.maxSignificantDigits

This is an interger value representing the maximum number of significant digits that the :py:class:`PQU` class will 
treat as having a significant meaning. Instances of :py:class:`PQU` with significantDigits greater than
:py:data:`maxSignificantDigits` are treated as a normal python float.

Classes
-------

This sections list the classes in the PQU module.

PQU
+++

.. autoclass:: pqu.PQU.PQU
    :members:
    :undoc-members:
    :show-inheritance:

pqu_float
+++++++++

.. autoclass:: pqu.PQU.pqu_float
    :members:
    :undoc-members:
    :show-inheritance:

pqu_uncertainty
+++++++++++++++

.. autoclass:: pqu.PQU.pqu_uncertainty
    :members:
    :undoc-members:
    :show-inheritance:

parsers
+++++++

.. autoclass:: pqu.PQU.parsers
    :members:
    :undoc-members:
    :show-inheritance:

PhysicalUnit
++++++++++++

.. autoclass:: pqu.PQU.PhysicalUnit
    :members:
    :undoc-members:
    :show-inheritance:

Functions
---------

This sections list some of the functions in the PQU module.

compare
+++++++

.. autofunction:: pqu.PQU.compare

floatToShortestString
+++++++++++++++++++++

.. autofunction:: pqu.PQU.floatToShortestString

valueOrPQ
+++++++++

.. autofunction:: pqu.PQU.valueOrPQ
