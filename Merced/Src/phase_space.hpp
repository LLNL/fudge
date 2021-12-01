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

namespace Phase
{
class phase_space_Ein_param;  // forward reference

//! Class for parameters for the 2-d quadrature over cm cosine and Eout_cm
// ---------------- class phase_space_Ecm_param ------------------
class phase_space_Ecm_param : public Egeom::Ecm_Elab_Ecm_param
{
public:
  // the mapping and phase-spcace parameters
  Maps::phase_space_map *PSmap;
  Qmeth::Quadrature_Rule  mu_quad_rule;  // quadrature rule for outgoing cosine
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
// ---------------- class phase_space ------------------
class phase_space
{
private:
  Maps::phase_space_map PSmap;

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
  bool get_Ein_range( const Ddvec::dd_vector& sigma, const Ddvec::dd_vector& multiple,
    const Ddvec::dd_vector& weight,
    const Lgdata::Flux_List& e_flux, const Egp::Energy_groups& Ein_groups );

  //! Loops over the E_out bins for a given cm_Eout bin
  //! \param transfer the computed transfer matrix
  //! \param Ein_param the quadrature parameters for intration over incident energy
  void Eout_ladder( Trf::T_matrix& transfer,
		    Phase::phase_space_Ein_param *Ein_param );

  //! Does the integration for one E-E' box between 2 eta = const hyperbolas
  //! \param transfer the computed transfer matrix
  //! \param Eout_count identifies the matrix row to update
  //! \param Ein_param the quadrature parameters for intration over incident energy
  void one_Ebox( Trf::T_matrix& transfer, int Eout_count,
		 Phase::phase_space_Ein_param *Ein_param );

  //! Adds to an element of transfer the integral between the intersections of 2 Eout_cm = const arcs with the Eout_lab box
  //! \param transfer the computed transfer matrix
  //! \param Eout_count identifies the matrix row to update
  //! \param Ein_param the quadrature parameters for intration over incident energy
  void update_T( Trf::T_matrix &transfer, int Eout_count,
		 Phase::phase_space_Ein_param *Ein_param );

public:
  double mProj;
  double mTarg;
  double mEject;
  double totalMass;
  double Q_value;
  double threshold;
  int numParticles;

  phase_space( ): mProj( -1.0 ), mTarg( -1.0 ), mEject( -1.0 ), 
     totalMass( -1.0 ), Q_value( -2.0e10 ), numParticles( 0 ) {}
  ~phase_space( ) {}

  //! Calculates the transfer matrix for this particle
  //! \param sigma the cross section data
  //! \param multiple the outgoing particle multiplicity data
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param transfer the transfer matrix
  void get_T( const Ddvec::dd_vector& sigma, const Ddvec::dd_vector& multiple,
	      const Ddvec::dd_vector& weight, Trf::T_matrix& transfer );

  //! Stores the masses
  //! \param particle_info the masses of the particles involved in the reaction
  void copy_masses( const Maps::particleInfo &particle_info );

};

//! Class for parameters for the 3-d quadrature over Ein, cm cosine, and Eout_cm
// ---------------- class phase_space_Ein_param ------------------
class phase_space_Ein_param : public Egeom::Ecm_Elab_Ein_param
{
public:
  Vhit::Vcm_Vlab_hit_list upper_hits;
  Vhit::Vcm_Vlab_hit_list lower_hits;

  // parameters for integration over center-of-mass energy
  Phase::phase_space_Ecm_param Ecm_params;

  // the mapping and phase-spcace parameters
  Maps::phase_space_map *map;
  Qmeth::Quadrature_Rule Eout_quad_rule;  // quadrature rule for outgoing energy
  Qmeth::Quadrature_Rule mu_quad_rule;  // quadrature rule for outgoing cosine

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
  void setup_map( Maps::phase_space_map *PSmap );

  // Sets up the parameters for integration over center-of-mass outgoing energy
  //! \param E_in energy of the incident particle
  void  set_Ecm_param( double E_in );

};
} // end of namespace Phase

// ************* functions to integrate ******************
namespace phase_space_F
{
  // ---------------- phase_space_F::mu_F ------------------
  //! Function for the 1-d quadrature over cm cosine
  //! Returns true if the interpolation is OK.
  //! \param mu the center-of-mass direction cosine of the outgoing particle
  //! \param mu_quad_param the function parameters
  //! \param value the value of the integrand, a set of Legendre coefficients
  bool mu_F( double mu, Qparam::QuadParamBase *mu_quad_param,
	     Coef::coef_vector *value );

  // ---------------- phase_space_F::Ecm_F ------------------
  //! Function for the 2-d quadrature over cm cosine and outgoing energy
  //! Returns true if the interpolation is OK.
  //! \param Eout_cm the center-of-mass energy of the outgoing particle
  //! \param Ecm_quad_param the function parameters
  //! \param value the value of the integrand, a set of Legendre coefficients
  bool Ecm_F( double Eout_cm, Qparam::QuadParamBase *Ecm_quad_param,
	      Coef::coef_vector *value );

  // ---------------- phase_space_F::Ecm_F_flip ------------------
  //! Function for the 2-d quadrature over cm cosine and outgoing energy with
  //! singularity sqrt( flip_Eout_cm )
  //! Returns true if the interpolation is OK.
  //! \param flip_Eout_cm the center-of-mass energy of the outgoing particle
  //! \param Ecm_quad_param the function parameters
  //! \param value the value of the integrand, a set of Legendre coefficients
  bool Ecm_F_flip( double flip_Eout_cm, Qparam::QuadParamBase *Ecm_quad_param,
		   Coef::coef_vector *value );

  // ---------------- phase_space_F::Ein_F ------------------
  //! Function for the 3-d quadrature over incident energy, cm cosine, and outgoing energy
  //! Returns true if the interpolation is OK.
  //! \param E_in the energy of the incident particle
  //! \param quad_param the function parameters
  //! \param value the value of the integrand, a set of Legendre coefficients
  bool Ein_F( double E_in, Qparam::QuadParamBase *quad_param,
	      Coef::coef_vector *value );
}

#endif
