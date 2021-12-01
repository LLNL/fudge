/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2015-05-20 -0800 (Wed, 20 May 2015) $
 * $Author: hedstrom $
 * $Id: relativistic.hpp 1 2015-05-20 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/

// header for the classes used in relativistic mechanics

#ifndef DEF_RELATIVISTIC
#define DEF_RELATIVISTIC

#include "mappings.hpp"  // for the Maps::particleInfo class

namespace Rel
{
// ----------- class Ep_vector -----------------
//! Class for an energy-momentum vector
class Ep_vector
{
private:

public:
  double Ep[ 4 ];

  //! constructor
  inline Ep_vector( ) {}

  //! destructor
  inline ~Ep_vector( ) {}

  //! gets the length of the momentum
  double get_p( );
};

// ----------- class relativistic_masses -----------------
//! Class for mass data used in relativistic boosts
class relativistic_masses
{
private:

public:
  Maps::particleInfo *rest_masses;
  double M_total;  // sum of all particle masses

  //! the threshold for the reaction
  double threshold;

  double Q_value;

  //! Default constructor
  inline relativistic_masses( ): Q_value(0.0) {}

  //! Default destructor
  inline ~relativistic_masses() {}

  //! Saves the rest masses; modified for consistency
  //! \param to_save the rest masses of the particles
  //! \param file_Q the Q value given in the input file
  void setup_masses( Maps::particleInfo *to_save, double *file_Q );

  //! Calculates the threshold energy
  void get_threshold( );

};

// ----------- class relativistic_map -----------------
//! Class for mapping used in relativistic 2-body reactions
class relativistic_map
{
private:
  //! the boost
  double cosh_chi;
  double sinh_chi;

public:
  double Tin_lab;  // lab kinetic energy of the incident particle
  double pin_lab;  // lab momentum of the incident particle
  double pin_cm;   // center-of-mass momentum of the incident particle
  double Minkowski;  // the Minkowski length, sqrt( E*E - p*p*c*c )
  double T_cm_out;  // cm kinetic energy of the outgoing particle
  double p_cm_out;  // cm momentum of the outgoing particle

  //! data for the reaction
  Rel::relativistic_masses *relMasses;

  //! Default constructor
  inline relativistic_map( ) {}
  
  //! Default destructor
  inline ~relativistic_map( ) {}

  //! Sets up the boost to the lab frame
  //! \param T_in_lab lab frame kinetic energy of the incident particle
  void set_boost( double T_in_lab );

  //! Boosts the ejected particle from the center-of-mass to the lab frame
  //! \param mu_cm direction cosine in the center-of-mass frame
  //! \param T_lab computed lab frame kinetic energy
  //! \param p_lab_parallel computed lab frame parallel component of the momentum
  void boost( double mu_cm, double *T_lab, double *p_lab_parallel );

  //! Calculates the center-of-mass energy and momentum for discrete 2-body reactions
  void get_p_cm_out( );

  //! Calculates the lab-frame kinetic energy for discrete 2-body reactions
  //! \param mu_cm center-of-mass frame direction cosine
  double get_T_lab_out( double mu_cm );

  //! Gets the value of mu_cm given the laboratory kinetic energy of outgoing particle.
  //! \param T_lab the desired lab frame kinetic energy
  double get_mu_cm( double T_lab );

  //! Gets the laboratory kinetic energy and cosine of outgoing particle
  //! \param mu_cm direction cosine in the center-of-mass frame
  //! \param Tout_lab computed lab frame kinetic energy
  //! \param mu_lab computed lab frame direction cosine
  void get_E_mu_lab( double mu_cm, double *Tout_lab, double *mu_lab );

  //! For exothermic reactions, returns the Newton incident energy for minimal outgoing lab energy
  double Newton_min_Eout( );

  //! Returns the incident energy for zero outgoing lab-frame energy, Newtonian
  double Newton_zero_Eout( );
};

// ----------- class relativistic_param -----------------
//! Class for the parameters in relativistic_F functions
class relativistic_param
{
private:
public:
  double mu_cm;   // for use in relativistic_F::T_out_lab

  Rel::relativistic_map relMap;

  //! Default constructor
  inline relativistic_param( ) {}

  //! Default destructor
  inline ~relativistic_param() {}

  //! Returns the incident kinetic energy for zero outgoing lab-frame kinetic energy, relativistic
  double zero_Eout( );

  //! For negative mu, gets the minimal lab frame outgoing kinetic energy and the
  //! corresponding incident kinetic energy
  //! \param mu the direction cosine in the center-of-mass frame
  //! \param guess_Ein an initial guess of the required incident kinetic energy
  //! \param Eout the computed minimal outgoing kinetic energy in the lab frame
  //! \param Ein the computed incident kinetic energy for this Eout
  void get_min_Eout( double mu, double guess_Ein, double *Eout, double *Ein );

  //! Finds the incident kinetic energies for given T_lab and mu_cm.
  //! Returns the incident kinetic energy for which T_lab is the outgoing kinetic energy
  //! \param T_lab lab-frame kinetic energy of ejected particle
  //! \param mu_cm center-of-mass direction cosine of ejected particle
  //! \param pair_0 ( Ein, Eout ) at lower incident kinetic energy
  //! \param pair_1 ( Ein, Eout ) at higher incident kinetic energy
  double find_hit( double T_lab, double mu_cm, const Ddvec::dd_entry &pair_0,
		   const Ddvec::dd_entry &pair_1 );

  //! For mu < 0 find the incident kinetic energy which minimizes Eout.
  //! Returns the location of the minimum and its outgoing kinetic energy.
  //! \param mu: the direction cosine in the center-of-mass frame
  //! \param pair_in_0: the lower bound and its outgoing kinetic energy
  //! \param pair_in_1: the upper bound and its outgoing kinetic energy
  //! \param tol: an error tolerance
  Ddvec::dd_entry find_bottom( double mu, const Ddvec::dd_entry &pair_in_0, const Ddvec::dd_entry &pair_in_1,
     double tol );

  //! Find the incident kinetic energy which minimizes Eout for mu = -1
  //  Ddvec::dd_entry find_lowest_bottom(  );
};

// ----------- class two_step_relativistic_masses -----------------
//! Class for mass data used in relativistic 2-step reactions
class two_step_relativistic_masses
{
private:

public:
  //! data for the reaction
  Rel::relativistic_masses step1Masses;
  Rel::relativistic_masses step2Masses;

  //! Default constructor
  inline two_step_relativistic_masses( ) {}
  
  //! Default destructor
  inline ~two_step_relativistic_masses( ) {}

  //! Sets up the map from center-of-mass to laboratory coordinates
  //! \param step1_particles, the particle masses for step 1
  //! \param file_Q1, the Q value for step 1 given in the file
  //! \param step2_particles, the particle masses for step 2
  //! \param file_Q2, the Q value for step 2 given in the file
  void set_masses( Maps::particleInfo *step1_particles, double file_Q1,
		   Maps::particleInfo *step2_particles, double file_Q2 );
  
  //! Returns the outgoing center-of-mass momentum for step 2
  double get_p_cm2( );
};

// ----------- class two_step_relativistic_map -----------------
//! Class for mapping used in relativistic 2-step reactions
class two_step_relativistic_map
{
private:
  //! the boost for step 1
  double cosh_chi1;
  double sinh_chi1;

  //! rotation for step 1
  double mu_lab1;  

  //! Initial steps for getting the laboratory kinetic energy of the outgoing particle for 2-step reactions
  //! \param T_in, kinetic energy of incident particle
  //! \param mucm_1, center-of-mass direction cosine of ejected particle
  //! \param mucm_2, the direction cosine for the second step
  //! \param w, the longitude for the second step, $-\pi/2 \le w \le \pi/2$.
  //! \param Ep_lab, computed lab-frame kinetic energy and momentun of ejected particle
  void two_step_start_E_lab( double T_in, double mucm_1,
			     double mucm_2, double w,
			     Rel::Ep_vector *Ep_lab );

public:
  //! data for the reaction
  Rel::two_step_relativistic_masses *twoStepMasses;

  Rel::relativistic_map step1Map;

  //! Default constructor
  inline two_step_relativistic_map( ) {}
  
  //! Default destructor
  inline ~two_step_relativistic_map( ) {}

  //! Sets up the boost and rotation for step 1
  //! \param T_in, kinetic energy of incident particle
  //! \param mucm1, center-of-mass direction cosine for first step
  void set_boost_Theta( double T_in, double mucm_1 );

  //! Does the boost to the lab frame, without the rotation
  //! \param Ep_cm, the energy-momentum vector in the center-of-mass frame
  //! \param Ep_lab, the energy-momentum vector in the lab frame
  void do_boost_Theta( const Rel::Ep_vector &Ep_cm, Rel::Ep_vector *Ep_lab );

  //! Does the rotation to coordinants parallel to the incident particle
  //! \param Ep_lab1, the energy-momentum vector in frame of step 1 outgoing particle
  //! \param Ep_lab, the energy-momentum vector in frame of incident particle
  void rotate_Theta( const Rel::Ep_vector &Ep_lab1, Rel::Ep_vector *Ep_lab );

  //! Gets the laboratory kinetic energy of the outgoing particle for 2-step reactions
  //! \param T_in, kinetic energy of incident particle
  //! \param mucm1, center-of-mass direction cosine of ejected particle
  //! \param mucm2, the direction cosine for the second step
  double two_step_get_E_lab( double T_in, double mucm1, double mucm2 );

  //! Gets the laboratory cosine of outgoing particle for 2-step reactions
  //! \param T_in, kinetic energy of incident particle
  //! \param mucm_1, center-of-mass direction cosine of ejected particle
  //! \param mucm_2, the direction cosine for the second step
  //! \param w, the longitude for the second step, $-\pi/2 \le w \le \pi/2$.
  double two_step_get_mu_lab( double T_in, double mucm_1, double mucm_2, double w );

  //! Gets the laboratory kinetic energy and cosine of outgoing particle for 2-step reactions
  //! \param T_in, kinetic energy of incident particle
  //! \param mucm_1, center-of-mass direction cosine of ejected particle
  //! \param mucm_2, the direction cosine for the second step
  //! \param w, the longitude for the second step, $-\pi/2 \le w \le \pi/2$.
  //! \param Tout_lab, computed lab-frame kinetic energy of ejected particle
  //! \param mu_lab, computed lab-frame direction cosine of ejected particle
  void two_step_get_E_mu_lab( double T_in, double mucm_1, double mucm_2, double w,
			      double *Tout_lab, double *mu_lab );
};

// ----------- class relativistic_2_step_param -----------------
//! Class for the parameters in relativistic_F functions
class relativistic_2_step_param
{
private:

public:
  // the particle masses and Qs
  Rel::two_step_relativistic_map relTwoStepMap;
  
  // for the function parameters
  double Tin_lab;
  double mu_cm_1;
  double mu_cm_2;

  //! Default constructor
  inline relativistic_2_step_param( ) {}

  //! Default destructor
  inline ~relativistic_2_step_param() {}

  //! Returns the direction cosine for the first step for given incident kinetic energy and desired outgoing kinetic energy
  //! \param T_in, the incident kinetic energy
  //! \param mucm2, direction cosine for step 2 ($\pm 1$)
  //! \param Tout_lab, the desired outgoing kinetic energy, lab frame
  double get_mucm1( double T_in, double mucm2, double Tout_lab );

  //! Returns the direction cosine for the second step for given translational energy and desired outgoing kinetic energy
  //! \param Tout_lab, the desired outgoing kinetic energy, lab frame
  double get_mucm2( double Tout_lab );

  //! Finds the incident kinetic energies for given T_lab and mu_cm.
  //! Returns the incident kinetic energy for which T_lab is the outgoing kinetic energy
  //! \param T_lab lab-frame kinetic energy of ejected particle
  //! \param pair_0 ( Ein, Tout ) at lower incident kinetic energy
  //! \param pair_1 ( Ein, Tout ) at higher incident kinetic energy
  double find_hit( double T_lab, const Ddvec::dd_entry &pair_0,
		   const Ddvec::dd_entry &pair_1 );

  //! Find the incident kinetic energy which minimizes Tout
  //! On exit pair_0 and pair_1 are tighter bounds.
  //! \param pair_0, (T_in, T_out) below the bottom
  //! \param pair_0, (T_in, T_out) above the bottom
  //! \param tol, when we are close enough
  Ddvec::dd_entry find_bottom( const Ddvec::dd_entry &pair_in_0,
			       const Ddvec::dd_entry &pair_in_1,
		 double tol );

  //! Find the incident kinetic energy which minimizes Tout for mucm_1 = mucm_2 = -1
  Ddvec::dd_entry find_lowest_bottom(  );
};
} // end of namespace Rel

// ************* functions *************
namespace relativistic_F
{
  // ------------------ p_from_T -----------------------
  //! Gets the momentum from the kinetic energy
  //! \param T particle kinetic energy
  //! \param E0 particle rest mass
  double p_from_T( double T, double E0 );

  // ------------------ T_from_p -----------------------
  //! Gets the kinetic energy from the momentum
  //! \param p particle momentum
  //! \param E0 particle rest mass
  double T_from_p( double p, double E0 );

  // ------------------ T_out_lab -----------------------
  //! Returns the kinetic energy of the emitted particle in the lab frame
  //! \param T_in_lab kinetic energy of the incident particle in the lab frame
  //! \param params the data for the relativistic boost
  double T_out_lab( double T_in_lab, void *params );

  // ------------------ p_out_lab -----------------------
  //! Returns the parallel momentum of the emitted particle in the lab frame.
  //! Routine is used to find when back scattering gives zero outgoing kinetic energy.
  //! \param T_in_lab kinetic energy of the incident particle in the lab frame
  //! \param params the data for the relativistic boost
  double p_out_lab( double T_in_lab, void *params );

  // ------------- two_step_T_out_lab --------------------
  //! For a 2-step reaction, returns the kinetic energy of the emitted particle in the lab frame as a function of T_in
  //! \param T_in kinetic energy of the incident particle in the lab frame
  //! \param params the data for the relativistic boost
  double two_step_T_out_lab( double T_in, void *params );
  
  // ------------------ mucm1_T_out_lab -----------------------
  //! For a 2-step reaction, returns the kinetic energy of the emitted particle in the lab frame as a function of mucm_1
  //! \param mucm_1, the direction cosine for step 1
  //! \param params the data for the relativistic boost
  double mucm1_T_out_lab( double mucm_1, void *params );

  // ------------------ mucm2_T_out_lab -----------------------
  //! For a 2-step reaction, returns the kinetic energy of the emitted particle in the lab frame as a function of mucm_2
  //! \param mucm_2, the direction cosine for step 2
  //! \param params the data for the relativistic boost
  double mucm2_T_out_lab( double mucm_2, void *params );
}

#endif
