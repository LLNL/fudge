/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2011-01-18 11:06:56 -0800 (Tue, 18 Jan 2011) $
 * $Author: hedstrom $
 * $Id: Legendre2Body.hpp 1 2011-01-18 11:06:56 -0800 hedstrom $
 *
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */
// define the classes used for ENDF Legendre expansions for discrete 2-body reactions

#ifndef LEGENDRE_2BODY_DEF
#define LEGENDRE_2BODY_DEF

// we use the Adist::angle_hit_list class
#include "angle_dist.hpp"

namespace Lg2b
{

class Legendre2d_param;  // forward declaration

//! Class for angular distributions as Legendre expansions
// ---------------- class Legendre_angle_dist ------------------
class Legendre_angle_dist : public std::list< Lgdata::Legendre_coefs >
{
private:
  double threshold;
  double threshold_out;  // outgoing energy at the threshold
  Ddvec::dd_entry flip;  // Ein for minimal Eout

  Maps::two_body_map map;
  Rel::relativistic_map relMap;
  
  //! The smallest incident energy for cross section, multiplicity,
  //! model weight, flux weight, and energy groups
  double E_first;

  int first_Ein;  // index of the left-hand end of the first significant energy bin
  int last_Ein;  // index of the right-hand end of the last significant energy bin

  //! Uses the mass difference to set the threshold
  void set_threshold( );

  //! Initializes the quadrature parameters
  //! \param mu_quad_rule rule for integration over direction cosine
  //! \param Ein_param parameters for quadrature over incident energy
  void initialize_param( Qmeth::Quadrature_Rule mu_quad_rule, 
     Legendre2d_param *Ein_param );

  //!  Gets the range of nontrivial incident energy bins; computes first_Ein and last_Ein
  //! returns true if the threshold is too high for the energy bins
  //! \param sigma the cross section data
  //! \param mult the multiplicity of the outgoing particle
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param e_flux approximate flux used to weight the transfer matrix
  //! \param Ein_groups the boundaries of the incident energy groups
  bool get_Ein_range( const Ddvec::dd_vector& sigma, const Ddvec::dd_vector& mult,
    const Ddvec::dd_vector& weight,
    const Lgdata::Flux_List& e_flux, const Egp::Energy_groups& Ein_groups );

  //! Sets up the map from center-of-mass to laboratory coordinates
  //! \param double top_E_in, the highest incident energy
  void setup_map( double top_E_in );

  //! Sets up the Ein_param->map from center-of-mass to laboratory coordinates
  //! \param Ein_param parameters for quadrature over incident energy
  void setup_param_map( Legendre2d_param *Ein_param );

  //! Sets the range of incident energies for this intergration
  //! \param Ein_param parameters for quadrature over incident energy
  void set_Ein_range( Legendre2d_param *Ein_param );

  //! Initializes the quadrature parameters
  //! \param Ein_param parameters for quadrature over incident energy
  void setup_data( Legendre2d_param *Ein_param );

  //! Adds to the transfer matrix for all E_out bins for a pair of incident energies.
  //! \param transfer the transfer matrix
  //! \param Ein_param parameters for quadrature over incident energy
  void Eout_ladder( Trf::T_matrix& transfer, Legendre2d_param *Ein_param );

  //! Does the integration for one E-E' box between 2 eta = const hyperbolas
  //! \param transfer the transfer matrix
  //! \param Eout_count count of the current outgoing energy bin
  //! \param Ein_param parameters for quadrature over incident energy
  void one_Ebox( Trf::T_matrix& transfer, int Eout_count, Legendre2d_param *Ein_param );

  //! Go to the next pair of incident energies.  Returns "true" when finished.
  //! \param E_in the next incident energy
  //! \param Ein_param parameters for quadrature over incident energy
  bool next_ladder( double E_in, Legendre2d_param *Ein_param );

  //! Adds to an element of transfer the integral between the intersections of 2 eta = const hyperbolas with the E-E' box
  //! \param transfer the transfer matrix
  //! \param Eout_count identify the row of the transfer matrix
  //! \param Ein_param parameters for quadrature over incident energy
  void update_T( Trf::T_matrix &transfer, int Eout_count, Legendre2d_param *Ein_param );

public:
  bool use_relativistic;
  Rel::relativistic_masses relMass;

  Terp::Interp_Type Ein_interp;  // interpolation between incident energies

  Legendre_angle_dist( ): use_relativistic( false ), Ein_interp( Terp::NOTSET ) {}
  ~Legendre_angle_dist( ) {}

  Maps::particleInfo particle_info;

  double Q;  // the reaction Q value

  //! Reads the data from Python
  //! \param input_file input file
  //! \param num_Ein number of incident energies for this reaction
  void read_data( Dpar::data_parser &input_file, int num_Ein );

  //! Calculates the transfer matrix for this particle
  //! \param sigma the cross section data
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param transfer the transfer matrix
  void get_T( const Ddvec::dd_vector& sigma, const Ddvec::dd_vector& weight,
	      Trf::T_matrix& transfer );

  // for debugging
  void print( );

};

//! Class for parameters for the 2-d quadrature
// ---------------- class Legendre2d_param ------------------
class Legendre2d_param : public Pbase::param_base
{
private:
public:
  double threshold;
  double threshold_out;  // outgoing energy at the threshold

  Maps::Newton_map_param Newton_map;
  Rel::relativistic_map relMap;

  Legendre_angle_dist::const_iterator left_data;        // 4 pointers to (eta, probability) data
  Legendre_angle_dist::const_iterator right_data;

  int num_negative;
  bool use_relativistic;

  // where the lower and upper eta values hit the quadrature box
  Adist::angle_hit_list lower_hits;
  Adist::angle_hit_list upper_hits;
 
  long int quad_count;  // number of 2-d quadratures
  long int Ein_F_count;  // number of calls to Legendre2Body_F::E_quad_F
  long int mu_F_count;  // number of calls to Legendre2Body_F::mu_cm_quad_F

  // quadrature rule for center-of-mass cosine
  Qmeth::Quadrature_Rule mu_quad_rule;
  Terp::Interp_Type Ein_interp;

  inline Legendre2d_param(): num_negative( 0 ), quad_count( 0 ),
			     Ein_F_count( 0 ), mu_F_count( 0 ) {}
  inline ~Legendre2d_param() {}

};

//! Class for parameters for the 1-d quadrature based on a set of Legendre coefficients
// ---------------- class Legendre_param ------------------
class Legendre_param : public Qparam::QuadParamBase
{
public:
  Lgdata::Legendre_coefs coefs;
  Maps::Newton_map_param *Newton_map;
  Rel::relativistic_map relMap;
  
  int num_negative;
  bool flag_set;
  bool use_relativistic;

  Terp::Interp_Type Ein_interp;
  inline Legendre_param( ): num_negative( 0 ), flag_set( false )  {}
  inline ~Legendre_param( ) {}

  //! Interpolates between two incident energies
  //! Returns false if there are problems
  //! \param Ein, the intermediate incident energy
  //! \param prev_coef, Legendre coefficients at the lower incident energy
  //! \param next_coef, Legendre coefficients at the higher incident energy
  bool interpolate( double Ein, Legendre_angle_dist::const_iterator prev_coef,
		    Legendre_angle_dist::const_iterator next_coef );
};

} // end of namespace Lg2b

// **************** functions to integrate **********
namespace Legendre2Body_F
{
  // ---------------- mu_cm_quad_F ------------------
  //! Function for the 1-d quadrature over mu_cm
  //! \param mu_cm center-of-mass direction cosine for outgoing particle
  //! \param mu_cm_quad_param parameters for this function
  //! \param value computed contribution to the transfer matrix
  bool mu_cm_quad_F( double mu_cm, Qparam::QuadParamBase *mu_cm_quad_param,
		     Coef::coef_vector *value );

  // ---------------- E_quad_F ------------------
  //! Function for the 2-d quadrature over (E_in, mu_cm )
  //! \param E_in inergy of the incident particle
  //! \param e_quad_param parameters for this function
  //! \param value computed contribution to the transfer matrix
  bool E_quad_F( double E_in, Qparam::QuadParamBase *e_quad_param,
		 Coef::coef_vector *value );
}

#endif
