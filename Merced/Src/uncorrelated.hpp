/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: uncorrelated.hpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// define the classes used for uncorrelated energy-angle distributions

#ifndef UNCORRELATED_CLASS
#define UNCORRELATED_CLASS

#include <iostream>
#include "angle_dist.hpp"
#include "Legendre2Body.hpp"
#include "energy_dist.hpp"

namespace Uncor
{
class uncorrelated_param;  // forward declaration

//! Class for uncorrelated energy-angle distributions
//--------------- class uncorrelated ----------------
class uncorrelated : public Ebase::energy_dist_base
{
private:
  // cumulative probabilities for cumulative points interpolation
  Cum::cumulative_prob_list cum_prob;

  //! The smallest incident energy for cross section, multiplicity,
  //! model weight, flux weight, and energy groups
  double E_first;

  int first_Ein;  // index of the left-hand end of the first significant energy bin
  int last_Ein;  // index of the right-hand end of the last significant energy bin

  Lg2b::Legendre_angle_dist Legendre_coef_data;   // mu data as Legendre coefficients

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

  //! Initializes the quadrature parameters
  //! \param Ein_param parameters for the 2d quadrature over energy
  void setup_param( Uncor::uncorrelated_param *Ein_param );

  //! Sets the range of incident energies for this intergration
  //! \param Ein_param parameters for the 2d quadrature over energy
  void set_Ein_range( Uncor::uncorrelated_param *Ein_param );

  //! Get the next incident energy range for tabular angular data
  //! \param  E_in the next incident energy
  //! \param Ein_param parameters for the 2d quadrature over energy
  bool next_tabular_ladder( double E_in, Uncor::uncorrelated_param *Ein_param );

  //! Get the next incident energy range for Legendre angular data
  //! \param  E_in the next incident energy
  //! \param Ein_param parameters for the 2d quadrature over energy
  bool next_Legendre_ladder( double E_in, Uncor::uncorrelated_param *Ein_param );

  //! Adds to the transfer matrix for all E_out bins for a pair of incident energies.
  //! \param transfer the computed transfer matrix
  //! \param Ein_param parameters for the 2d quadrature over energy
  void Eout_ladder( Trf::T_matrix& transfer,
		    Uncor::uncorrelated_param *Ein_param );

  //! Starts one staircase of the Eout data
  //! \param Ein_param the quadrature parameters
  void start_Eout( Uncor::uncorrelated_param *Ein_param );

  //! Go to the next set of (E_out, probability) pairs for unitbase
  //! interpolation in incident energy.
  //! Returns true if we are finished with the data.
  //! \param Ein_param parameters for the 2d quadrature over energy
  bool next_Eout_ubase( Uncor::uncorrelated_param *Ein_param );

  //! Go to the next set of (E_out, probability) pairs for histogram
  //! interpolation in incident energy.
  //! Returns true if we are finished with the data.
  //! \param Ein_param parameters for the 2d quadrature over energy
  bool next_Eout_flat( Uncor::uncorrelated_param *Ein_param );

  //! Go to the next pair of incident energies.  Returns "true" when finished.
  //! \param E_in the next incident energy
  bool next_ladder( double E_in, Uncor::uncorrelated_param *Ein_param );

  //! Integrates over one E-E' box
  //! \param transfer the computed transfer matrix
  //! \param Eout_count the current row of the transfer matrix
  //! \param Ein_param parameters for the 2d quadrature over energy
  void one_Ebox( Trf::T_matrix& transfer, int Eout_count,
		 Uncor::uncorrelated_param *Ein_param );

  //! Increments the transfer matrix
  //! \param transfer the transfer matrix to compute
  //! \param Eout_count identifies the matrix entry to update
  //! \param Ein_param parameters for the 2d quadrature over energy
  void update_T( Trf::T_matrix &transfer, int Eout_count,
		 Uncor::uncorrelated_param *Ein_param );

public:
  bool mu_table;    // angular data as a table if true; otherwise Legendre
  //  Terp::two_d_interp Ein_interp;  // for interpolation of incident energy
  //  Terp::Interp_Type Eout_interp;  // for interpolation of E_out probability
  Terp::Interp_Type mu_interp;
  Adist::angle_dist *angles;

  uncorrelated( );
  ~uncorrelated( );

  //! Reads the Legendre coefficients
  //! \param infile input file
  //! \param num_Ein number of incident energies for this reaction
  void read_Legendre( Dpar::data_parser& infile, int num_Ein );

  //! Reads the ENDL data
  //! \param infile input file
  //! \param num_I4 number of incident energies with pairs ( E_out, probability )
  //! \param ang the angular probability densities (already read)
  void read_data( Dpar::data_parser& infile, int num_I4, Adist::angle_dist *ang );

  // Calculates the transfer matrix for this particle.
  //! \param sigma the cross section data
  //! \param mult the outgoing particle multiplicity data
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param transfer the transfer matrix
  void get_T( const Ddvec::dd_vector& sigma, const Ddvec::dd_vector& mult,
	      const Ddvec::dd_vector& weight,
    Trf::T_matrix& transfer );
};

//! Class for parameters for the 2-d quadrature
// ---------------- class uncorrelated_param ------------------
class uncorrelated_param : public Pbase::param_base
{
private:

  //! Interpolates (Eout, probability) data to the higher common cumulative probability
  void setup_high_A( double higher_A );

  //! Interpolates the angular data
  //! Returns true if the interpolation is OK
  //! \param E_in, the incident energy to interpolate to
  bool interpolate_mu( double E_in );

public:
  //! the interpolated unit-base mapping parameters
  Ddvec::unit_base_map mid_ubase_map;

  Ddvec::dd_entry Eout_0_range; // values of Eout for left_data and next_left_data
  Ddvec::dd_entry Eout_1_range; // values of Eout for right_data and next_right_data

  // pointers to the current pair of energy distirbutions
  Uncor::uncorrelated::const_iterator this_E_dist;
  Uncor::uncorrelated::const_iterator next_E_dist;

  // copies of energy distributions for direct interpolation
  Ddvec::dd_vector this_Ein_direct;
  Ddvec::dd_vector next_Ein_direct;

  // pointers to (E', probability) data, not interpolated in Eout
  Ddvec::dd_vector::const_iterator left_data;      
  Ddvec::dd_vector::const_iterator next_left_data;
  Ddvec::dd_vector::const_iterator last_left_data;
  Ddvec::dd_vector::const_iterator right_data;
  Ddvec::dd_vector::const_iterator next_right_data;
  Ddvec::dd_vector::const_iterator last_right_data;

  Cum::cumulative_prob_list::const_iterator left_cum_prob;  // cumulative probability at lower Ein and lower Eout
  Cum::cumulative_prob_list::const_iterator next_left_cum_prob;  // cumulative probability at lower Ein and higher Eout
  Cum::cumulative_prob_list::const_iterator right_cum_prob;  // cumulative probability at higher Ein and lower Eout
  Cum::cumulative_prob_list::const_iterator next_right_cum_prob;  // cumulative probability at higher Ein and higher Eout

  // current unit-base outgoing energy range
  double current_prev_Eout;
  double current_next_Eout;

  // (E', probability) data interpolated in Eout to current_prev_Eout and current_next_Eout
  Ddvec::cum_points_pair Ein0_data;  // lower incident energy
  Ddvec::cum_points_pair Ein1_data;  // higher incident energy
  Ddvec::cum_points_pair mid_Eout_data;  // intermediate incident energy

  double left_Ein;      // the incident energies for the probability data
  double right_Ein;
  //  int order;             // The Legendre order of the problem
  int L_order;           // The Legendre order of this term

  Lgdata::Legendre_coefs mu_integral[ 2 ];   // integrals over mu
  Lgdata::Legendre_coefs mid_mu_integral;    // interpolated mu_intergral

  Adist::angle_dist::iterator this_mu_dist;  // for tabular angular data
  Adist::angle_dist::iterator next_mu_dist;
  int prev_mu;                         // 0 or 1, index in mu_intergral
  bool mu_table;    // angular data as a table if true; otherwise Legendre
  Lg2b::Legendre_angle_dist::const_iterator prev_L_coefs;  // for Legendre coefficients
  Lg2b::Legendre_angle_dist::const_iterator next_L_coefs;

  // where the lower and upper eta values hit the quadrature box
  Box::energy_hit_list lower_hits;
  Box::energy_hit_list upper_hits;
 
  long int quad_count;  // number of 2d quadratures
  long int Ein_F_count;  // number of calls to  uncorrelated_F::Ein_F
  long int mu_quad_count;  // number of integrations over mu
  long int mu_F_count;  // number of calls to uncorrelated_F::mu_F

  Terp::two_d_interp Ein_interp; // interpolation rule for incident energy
  Terp::Interp_Type Eout_interp; // interpolation rule for outgoing energy
  Terp::Interp_Type mu_interp; // interpolation rule for direction cosine

  inline uncorrelated_param( ): quad_count( 0 ), Ein_F_count( 0 ),
				mu_quad_count( 0 ), mu_F_count( 0 ) {}
  inline ~uncorrelated_param( ) {}

  //! Does truncation or extrapolation of data for direct linlin interpolation
  void setup_direct( );
  
  //! Sets up the data for interpolation with respect to incident energy
  void set_data( );

  //! Sets up the data for cumulative-points interpolation in incident energy
  void setup_Ein_cum_prob( );

  //! Interpolates to set up for the integration under unitbase interplation
  //! Returns true if the interpolation is OK
  //! \param E_in energy of incident particle
  bool interp_data_ubase( double E_in );

  //! Interpolates to set up for the integration under histogram interplation
  //! Returns true if the interpolation is OK
  //! \param E_in energy of incident particle
  bool interp_data_flat( double E_in );

  //! Interpolates to set up for the integration linlin direct interplation
  //! Returns true if the interpolation is OK
  //! \param E_in energy of incident particle
  bool interp_linlin_direct( double E_in );

  //! Go to the next set of (E_out, probability) pairs for cumulative
  //! points interpolation in incident energy.
  //! Returns true if we are finished with the data.
  bool next_Eout_cum_prob( );

};

} // end of namespace Uncor

namespace uncorrelated_F
{
  // -------------- functions to integrate -------------------

  // ---------------- uncorrelated_F::Ein_F ------------------
  //! Function for the quadrature over incident energy
  //! \param E_in the energy of the incident particle
  //! \param param the function parameters
  //! \param value the value of the integrand, a set of Legendre coefficients
  bool Ein_F( double E_in, Qparam::QuadParamBase *param, Coef::coef_vector *value );
  
}  // end of namespace uncorrelated_F

#endif
