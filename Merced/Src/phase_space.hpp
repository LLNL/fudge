/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2010-10-16 (Tue, Oct 16, 2010) $
 * $Author: hedstrom $
 * $Id: phase_space.cpp 1 2010-10-16 hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// declaration of the classes used to handle phase-space energy model

#ifndef PHASE_SPACE_CLASS
#define PHASE_SPACE_CLASS

#include "Ecm_Elab_geom.hpp"
#include "transfer.hpp"

class phase_space_Ein_param;  // forward reference

//! Class for parameters for the 2-d quadrature over cm cosine and Eout_cm
// ---------------- class phase_space_Ecm_param ------------------
class phase_space_Ecm_param : public Ecm_Elab_Ecm_param
{
public:
  // the mapping and phase-spcace parameters
  phase_space_map *PSmap;
  Quadrature_Method  mu_quad_method;  // quadrature method for outgoing cosine
  long int mu_F_count;  // number of calls to phase_space_F::mu_F

  inline phase_space_Ecm_param( ): mu_F_count( 0 ) {}
  inline ~phase_space_Ecm_param( ) {}

  //! Sets up the data for this incident energy
  //! \param E_in energy of the incident particle
  //! \param Eoutmin minimum lab-frame energy for this outgoing energy bin
  //! \param Eoutmax maximum lab-frame energy for this outgoing energy bin
  void setup( double E_in, double Eoutmin, double Eoutmax );

  //! Gets the center-of-mass energy probability density
  //! \param E_cm outgoing center-of-mass energy for this particle
  //! \param E_cm maximum outgoing center-of-mass energy for this data
  inline double get_Ecm_prob( double E_cm, double Ecm_max )
  { return PSmap->get_prob( E_cm, Ecm_max ); }
};

//! class for Phase_Space data
// ---------------- class Phase_Space ------------------
class phase_space
{
private:
  phase_space_map PSmap;

  int first_Ein;  // index of the left-hand end of the first significant energy bin
  int last_Ein;  // index of the right-hand end of the last significant energy bin

  //! Ensures that the input data is complete
  void check_input( );

  //! Sets up the map from center-of-mass to laboratory coordinates
  void setup_map( );

  //!  Gets the range of nontrivial incident energy bins; computes first_Ein and last_Ein
  //! returns true if the threshold is too high for the energy bins
  //! \param sigma the cross section data
  //! \param multiple the multiplicity of the outgoing particle
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param e_flux approximate flux used to weight the transfer matrix
  //! \param Ein_groups the boundaries of the incident energy groups
  bool get_Ein_range( const dd_vector& sigma, const dd_vector& multiple,
    const dd_vector& weight,
    const Flux_List& e_flux, const Energy_groups& Ein_groups );

  //! Loops over the E_out bins for a given cm_Eout bin
  //! \param transfer the computed transfer matrix
  //! \param Ein_param the quadrature parameters for intration over incident energy
  void Eout_ladder( T_matrix& transfer, phase_space_Ein_param *Ein_param );

  //! Does the integration for one E-E' box between 2 eta = const hyperbolas
  //! \param transfer the computed transfer matrix
  //! \param Eout_count identifies the matrix row to update
  //! \param Ein_param the quadrature parameters for intration over incident energy
  void one_Ebox( T_matrix& transfer, int Eout_count, phase_space_Ein_param *Ein_param );

  //! Adds to an element of transfer the integral between the intersections of 2 Eout_cm = const arcs with the Eout_lab box
  //! \param transfer the computed transfer matrix
  //! \param Eout_count identifies the matrix row to update
  //! \param Ein_param the quadrature parameters for intration over incident energy
  void update_T( T_matrix &transfer, int Eout_count, phase_space_Ein_param *Ein_param );

public:
  double mProj;
  double mTarg;
  double mEject;
  double totalMass;
  double Q_value;
  int numParticles;

  phase_space( ): mProj( -1.0 ), mTarg( -1.0 ), mEject( -1.0 ), 
     totalMass( -1.0 ), Q_value( -2.0e10 ), numParticles( 0 ) {}
  ~phase_space( ) {}

  //! Calculates the transfer matrix for this particle
  //! \param sigma the cross section data
  //! \param multiple the outgoing particle multiplicity data
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param transfer the transfer matrix
  void get_T( const dd_vector& sigma, const dd_vector& multiple,
	      const dd_vector& weight, T_matrix& transfer );

  //! Stores the masses
  //! \param particle_info the masses of the particles involved in the reaction
  void copy_masses( const particleInfo &particle_info );

};

//! Class for parameters for the 3-d quadrature over Ein, cm cosine, and Eout_cm
// ---------------- class phase_space_Ein_param ------------------
class phase_space_Ein_param : public Ecm_Elab_Ein_param
{
public:
  Vcm_Vlab_hit_list upper_hits;
  Vcm_Vlab_hit_list lower_hits;

  // parameters for integration over center-of-mass energy
  phase_space_Ecm_param Ecm_params;

  // the mapping and phase-spcace parameters
  phase_space_map *map;
  Quadrature_Method Eout_quad_method;  // quadrature method for outgoing energy
  Quadrature_Method mu_quad_method;  // quadrature method for outgoing cosine

  long int quad_count;  // number of 3-d quadratures
  long int Ein_F_count;  // number of calls to phase_space_F::Ein_F
  long int Ecm_F_count;  // number of calls to phase_space_F::Ecm_F
  long int mu_F_count;  // number of calls to phase_space_F::mu_F
  int Vcm_hit_count;     // number of calls to quad_F::integrate in phase_space_F::Ein_F

  phase_space_Ein_param( ): quad_count( 0 ), Ein_F_count( 0 ), Ecm_F_count( 0 ),
     mu_F_count( 0 ), Vcm_hit_count( 0 ) {}
  ~phase_space_Ein_param( ) {}

  //! Copies the data for mapping to the lab frame
  //! \param PSmap data for mapping to the lab frame
  void setup_map( phase_space_map *PSmap );

  // Sets up the parameters for integration over center-of-mass outgoing energy
  //! \param E_in energy of the incident particle
  void  set_Ecm_param( double E_in );

};

// ************* functions to integrate ******************
namespace phase_space_F
{
  // ---------------- phase_space_F::mu_F ------------------
  //! Function for the 1-d quadrature over cm cosine
  //! \param mu the center-of-mass direction cosine of the outgoing particle
  //! \param mu_quad_param the function parameters
  //! \param value the value of the integrand, a set of Legendre coefficients
  void mu_F( double mu, QuadParamBase *mu_quad_param, coef_vector *value );

  // ---------------- phase_space_F::Ecm_F ------------------
  //! Function for the 2-d quadrature over cm cosine and outgoing energy
  //! \param Eout_cm the center-of-mass energy of the outgoing particle
  //! \param Ecm_quad_param the function parameters
  //! \param value the value of the integrand, a set of Legendre coefficients
  void Ecm_F( double Eout_cm, QuadParamBase *Ecm_quad_param, coef_vector *value );

  // ---------------- phase_space_F::Ein_F ------------------
  //! Function for the 3-d quadrature over incident energy, cm cosine, and outgoing energy
  //! \param E_in the energy of the incident particle
  //! \param quad_param the function parameters
  //! \param value the value of the integrand, a set of Legendre coefficients
  void Ein_F( double E_in, QuadParamBase *quad_param, coef_vector *value );
}

#endif
