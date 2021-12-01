The Styles Module
=================

.. todo::
   - The python implementation for fudge.styles.multiGroup.convertUnits() is incomplete and consequently
     it's associated docstring.
   - The python implementation for fudge.styles.heatedMultiGroup.convertUnits() is incomplete and consequently
     it's associated docstring.

Styles Module: Ancestry diagrams
--------------------------------

The inheretence diagram for the classes in the ``styles`` module are shown below.

.. inheritance-diagram:: fudge.styles

Supporting Background Theory
----------------------------

The Singularity in the Coulomb Scattering Cross Section
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The differential cross section for Rutherford scattering into a solid angle 
:math:`d\Omega = 2\pi sin\theta d\theta` may be expressed as (Podgoršak, 2006; Eisberg and Resnick, 1985)

.. math::
   \frac{d\sigma}{d\Omega} = \left(\frac{D}{4(1 - cos\theta)}\right)^2,

where :math:`D` is the effective characteristic distance between the charge particle nucleus scattering, and
:math:`\theta` is the center of mass scattering angle.

Noting that :math:`\mu = cos\theta`, this may be rewritten as,

.. math::
   \frac{d\sigma}{d\mu} = -2\pi \left(\frac{D}{4(1 - \mu)}\right)^2,

which only has finite values for :math:`\mu \neq 1`.   
   


Styles Module: Python Source Documentation 
------------------------------------------

.. automodule:: fudge.styles
   :members:
   :undoc-members:
   :show-inheritance:


