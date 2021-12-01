The Neutron Transport Equation
==============================
In steady-state applications a discrete version of the linearized Boltzmann equation is solved 
(Cullen, 2010; Shultis & Faw, 2000) and the transfer matrix is used to approximate the kernel of the integral 
operator in this equation. This transport equation, in the absence of an external source, may be 
written as

.. math::
    \bm{\Omega} \cdot \nabla \psi \left( \bm{r}, E, \bm{\Omega} \right) + N \left( \bm{r} \right)
    \sigma_t \left( \bm{r}, E \right) \psi \left( \bm{r}, E, \bm{\Omega} \right) =
    \frac{N \left( \bm{r} \right)}{4\pi} \int_{E ^{\prime}} dE ^ \prime \int_{\bm{\Omega}^{\prime}} d\Omega^{\prime}
    \sigma \left(\bm{r}, E ^{\prime} \rightarrow E, 
    \bm{\Omega}^{\prime} \rightarrow \bm{\Omega} \right)
    \psi \left( \bm{r}, E^{\prime}, \bm{\Omega}^{\prime} \right)

where 
    - :math:`\psi \left( \bm{r}, E, \bm{\Omega} \right)` is the angular flux of neutrons with energy
      :math:`E`, location :math:`\bm{r}` and moving in direction :math:`\bm{\Omega}` ;
    - :math:`\sigma_t \left( \bm{r}, E \right)` is the total microscopic cross section at position
      :math:`\bm{r}` for a neutron with energy :math:`E` ;
    - :math:`\sigma \left(\bm{r}, E ^{\prime} \rightarrow E, \bm{\Omega}^{\prime} \rightarrow \bm{\Omega} \right)`
      is the microscopic differential cross section for a neutron with energy :math:`E ^{\prime}` and direction
      :math:`\bm{\Omega}^{\prime}` which results in secondary neutrons with energy :math:`E` and 
      direction :math:`\bm{\Omega}` at location :math:`\bm{r}`. This includes all reactions that produce secondary 
      neutrons (e.g. fission, :math:`(n,2n)`, etc.); and
    - :math:`N \left( \bm{r} \right)` is the material number density at position :math:`\bm{r}`.

The microscopic differential cross section or kernel is expressed in various forms in nuclear data processing and
these different forms are investigated in the next section. These various forms may be expressed in either the 
laboratory or the center of mass coordinate systems. In the above transport equation the energies 
:math:`E` and :math:`E ^{\prime}` and the directions :math:`\bm{\Omega}` and :math:`\bm{\Omega} ^{\prime}`
are in the laboratory coordinate system.


The differential cross section
------------------------------

Assuming an infinite,  homogeneous and isotropic medium, the angular dependance is reduced to the angle between the
initial and final neutron directions (i.e. :math:`\bm{\Omega}^{\prime}` and :math:`\bm{\Omega}`). 
This allow for the expression of the differential cross section in terms of the cosine of the scattering angle 
(:math:`\mu = \bm{\Omega}^{\prime} \cdot \bm{\Omega}`), i.e. 

.. math::
    \sigma \left(E ^{\prime} \rightarrow E, \bm{\Omega}^{\prime} \rightarrow \bm{\Omega} \right) =
    \sigma \left(E ^{\prime} \rightarrow E, \mu\right)

The differential cross section may be further expressed in terms of the contributions from the individual
reactions (Cullen, 2010). 

.. math::
    \sigma \left(E ^{\prime} \rightarrow E, \mu\right) = \sum_k \sigma_k \left(E ^{\prime} \rightarrow E, \mu\right)


These reaction kernels may be further reduced to their component factors, i.e.

.. math::
    \sigma_k \left(E ^{\prime} \rightarrow E, \mu\right) = M_k\left( E ^{\prime} \right)
    \sigma_k \left( E ^{\prime} \right) P_k \left(E, \mu \mid E ^{\prime} \right)

where
    - :math:`M_k\left( E ^{\prime} \right)` is the multiplicity;
    - :math:`\sigma_k \left( E ^{\prime} \right)` is the reaction cross section for process :math:`k`; and
    - :math:`P_k \left(E, \mu \mid E ^{\prime} \right)` is the joint probability distribution for neutrons with energy
      :math:`E` and angle :math:`\mu` given an incident energy of :math:`E ^{\prime}` for process
      :math:`k`. Note (Trkov et al., 2018) that 
      :math:`\int_{E, \mu} P_k \left(E, \mu \mid E ^{\prime} \right) = 1`.

The joint probability distribution :math:`P_k \left(E, \mu \mid E ^{\prime} \right)` is also referred to as the 
double differential probability density (Hedstrom et al., 2016)


Uncorrelated energy-angle probability densities
-----------------------------------------------
The simpliest form of the joint probability distribution in the previous section is for a uncorrelated dependence on
the outgoing energy :math:`E` and direction cosine :math:`\mu` (Hedstrom et al., 2016), i.e.

.. math::
   P_k \left(E, \mu \mid E ^{\prime} \right) = P_k \left(E\mid E ^{\prime} \right) \times 
   P_k \left(\mu \mid E ^{\prime} \right)

The energies and direction cosines may be either in the laboratory or center-of-mass coordinate systems.
   
Correlated energy-angle probability densities
-----------------------------------------------

.. todo::
   Expand description of correlated energy-angle probability densities and include Legendre expansions, 
   Kalbach-Mann model, representation for use in Monte Carlo codes, etc.

.. The latter probability distribution may be further reduced to
.. 
..     P_k \left(E, \mu \mid E ^{\prime} \right) = P_k \left(E \mid E ^{\prime} \right) \times
..     P_k \left(\mu \mid E, E ^{\prime}\right) = P_k \left(\mu \mid E ^{\prime} \right) \times
..     P_k \left(E \mid \mu, E ^{\prime}\right)
.. 
.. where
..     - :math:`P_k \left(E, E ^{\prime}\right)` is the probability for exit energy :math:`E` given
..       incident energy :math:`E ^{\prime}`;
..     - :math:`P_k \left(\mu, E ^{\prime}\right)` is the probability for angle :math:`\mu` given
..       incident energy :math:`E ^{\prime}`;
..     - :math:`P_k \left(\mu \mid E, E ^{\prime}\right)` is the probability for angle :math:`\mu` given
..       exit energy :math:`E` and incident energy :math:`E ^{\prime}`; and
..     - :math:`P_k \left(E \mid \mu, E ^{\prime}\right)` is the probability for exit angle :math:`E` given
..       angle :math:`\mu` and incident energy :math:`E ^{\prime}`.
..     - :math:`P_k \left(\mu \mid E ^{\prime} \right)` is the angular distribution as a function of the 
..       incident energy. This is usually given in the center-of-mass system for correlated distributions
..       and the laboratory system for uncorrelated distributions.
..    - :math:`P_k \left(E, \mu \mid E ^{\prime}\right)` is the energy distribution which correlates
..      :math:`\mu` to :math:`tem for uncorrelated distributions.
..    - :math:`P_k \left(\mu, E ^{\prime} \rightarrow E \right)` for correlated distributions as oppose to the 
..      actual energy distribution (independent of :math:`mu`) in the case of uncorrelated distributions.
 