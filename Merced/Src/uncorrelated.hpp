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

class uncorrelated_param;  // forward declaration

//! Class for parameters for the quadrature over mu
// ---------------- class mu_param ------------------
class mu_param: public QuadParamBase
{
public:
  dd_vector::const_iterator left_data; 
  dd_vector::const_iterator right_data;

  inline mu_param() {}
  inline ~mu_param() {}

  //! Evaluate by linear-linear interpolation
  //! \param mu an intermediate direction cosine
  double value( double mu );
};

// ---------------- class angular_moments ------------------
class angular_moments : public Legendre_coefs
{
public:
  inline angular_moments( ) {}
  inline ~angular_moments( ) {}

  //! Initializes the elements to zero
  void set_zero( );

  //! Addition operator
  //! \param to_add the Legendre moments to add
  angular_moments& operator+=( const coef_vector &to_add );

  //! Evaluates the Legendre moments
  //! Returns the number of function evaluations used in the quadrature
  //! \param this_mu_dist pairs ( direction cosine, probability density )
  int get_moment( const dd_vector& this_mu_dist );
};

//! Class for uncorrelated energy-angle distributions
//--------------- class uncorrelated ----------------
class uncorrelated : public energy_dist_base
{
private:
  // cumulative probabilities for cumulative points interpolation
  cumulative_prob_list cum_prob;

  //! The smallest incident energy for cross section, multiplicity,
  //! model weight, flux weight, and energy groups
  double E_first;

  int first_Ein;  // index of the left-hand end of the first significant energy bin
  int last_Ein;  // index of the right-hand end of the last significant energy bin

  Legendre_angle_dist Legendre_coef_data;   // mu data as Legendre coefficients

  //!  Gets the range of nontrivial incident energy bins; computes first_Ein and last_Ein
  //! returns true if the threshold is too high for the energy bins
  //! \param sigma the cross section data
  //! \param mult the multiplicity of the outgoing particle
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param e_flux approximate flux used to weight the transfer matrix
  //! \param Ein_groups the boundaries of the incident energy groups
  bool get_Ein_range( const dd_vector& sigma, const dd_vector& mult,
    const dd_vector& weight,
    const Flux_List& e_flux, const Energy_groups& Ein_groups );

  //! Initializes the quadrature parameters
  //! \param Ein_param parameters for the 2d quadrature over energy
  void setup_param( uncorrelated_param *Ein_param );

  //! Sets the range of incident energies for this intergration
  //! \param Ein_param parameters for the 2d quadrature over energy
  void set_Ein_range( uncorrelated_param *Ein_param );

  //! Get the next incident energy range for tabular angular data
  //! \param  E_in the next incident energy
  //! \param Ein_param parameters for the 2d quadrature over energy
  bool next_tabular_ladder( double E_in, uncorrelated_param *Ein_param );

  //! Get the next incident energy range for Legendre angular data
  //! \param  E_in the next incident energy
  //! \param Ein_param parameters for the 2d quadrature over energy
  bool next_Legendre_ladder( double E_in, uncorrelated_param *Ein_param );

  //! Adds to the transfer matrix for all E_out bins for a pair of incident energies.
  //! \param transfer the computed transfer matrix
  //! \param Ein_param parameters for the 2d quadrature over energy
  void Eout_ladder( T_matrix& transfer, uncorrelated_param *Ein_param );

  //! Starts one staircase of the Eout data
  //! \param Ein_param the quadrature parameters
  void start_Eout( uncorrelated_param *Ein_param );

  //! Go to the next set of (E_out, probability) pairs for unitbase
  //! interpolation in incident energy.
  //! Returns true if we are finished with the data.
  //! \param Ein_param parameters for the 2d quadrature over energy
  bool next_Eout_ubase( uncorrelated_param *Ein_param );

  //! Go to the next set of (E_out, probability) pairs for histogram
  //! interpolation in incident energy.
  //! Returns true if we are finished with the data.
  //! \param Ein_param parameters for the 2d quadrature over energy
  bool next_Eout_flat( uncorrelated_param *Ein_param );

  //! Go to the next pair of incident energies.  Returns "true" when finished.
  //! \param E_in the next incident energy
  bool next_ladder( double E_in, uncorrelated_param *Ein_param );

  //! Integrates over one E-E' box
  //! \param transfer the computed transfer matrix
  //! \param Eout_count the current row of the transfer matrix
  //! \param Ein_param parameters for the 2d quadrature over energy
  void one_Ebox( T_matrix& transfer, int Eout_count, uncorrelated_param *Ein_param );

  //! Increments the transfer matrix
  //! \param transfer the transfer matrix to compute
  //! \param Eout_count identifies the matrix entry to update
  //! \param Ein_param parameters for the 2d quadrature over energy
  void update_T( T_matrix &transfer, int Eout_count, uncorrelated_param *Ein_param );

public:
  bool mu_table;    // angular data as a table if true; otherwise Legendre
  //  two_d_interp Ein_interp;  // for interpolation of incident energy
  //  Interp_Type Eout_interp;  // for interpolation of E_out probability
  Interp_Type mu_interp;
  angle_dist *angles;

  uncorrelated( );
  ~uncorrelated( );

  //! Reads the Legendre coefficients
  //! \param infile input file
  //! \param num_Ein number of incident energies for this reaction
  void read_Legendre( data_parser& infile, int num_Ein );

  //! Reads the ENDL data
  //! \param infile input file
  //! \param num_I4 number of incident energies with pairs ( E_out, probability )
  //! \param ang the angular probability densities (already read)
  void read_data( data_parser& infile, int num_I4, angle_dist *ang );

  // Calculates the transfer matrix for this particle.
  //! \param sigma the cross section data
  //! \param mult the outgoing particle multiplicity data
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param transfer the transfer matrix
  void get_T( const dd_vector& sigma, const dd_vector& mult, const dd_vector& weight,
    T_matrix& transfer );
};

//! Class for parameters for the 2-d quadrature
// ---------------- class uncorrelated_param ------------------
class uncorrelated_param : public param_base
{
private:

  //! Interpolates (Eout, probability) data to the higher common cumulative probability
  void setup_high_A( double higher_A );

public:
  //! the interpolated unit-base mapping parameters
  unit_base_map mid_ubase_map;

  dd_entry Eout_0_range; // values of Eout for left_data and next_left_data
  dd_entry Eout_1_range; // values of Eout for right_data and next_right_data

  // pointers to the current pair of energy distirbutions
  uncorrelated::const_iterator this_E_dist;
  uncorrelated::const_iterator next_E_dist;

  // pointers to (E', probability) data, not interpolated in Eout
  dd_vector::const_iterator left_data;      
  dd_vector::const_iterator next_left_data;
  dd_vector::const_iterator right_data;
  dd_vector::const_iterator next_right_data;

  cumulative_prob_list::const_iterator left_cum_prob;  // cumulative probability at lower Ein and lower Eout
  cumulative_prob_list::const_iterator next_left_cum_prob;  // cumulative probability at lower Ein and higher Eout
  cumulative_prob_list::const_iterator right_cum_prob;  // cumulative probability at higher Ein and lower Eout
  cumulative_prob_list::const_iterator next_right_cum_prob;  // cumulative probability at higher Ein and higher Eout

  // current unit-base outgoing energy range
  double current_prev_Eout;
  double current_next_Eout;

  // (E', probability) data interpolated in Eout to current_prev_Eout and current_next_Eout
  cum_points_pair Ein0_data;  // lower incident energy
  cum_points_pair Ein1_data;  // higher incident energy
  cum_points_pair mid_Eout_data;  // intermediate incident energy

  double left_Ein;      // the incident energies for the probability data
  double right_Ein;
  //  int order;             // The Legendre order of the problem
  int L_order;           // The Legendre order of this term

  angular_moments mu_integral[ 2 ];   // integrals over mu
  Legendre_coefs mid_mu_integral;    // interpolated mu_intergral

  angle_dist::iterator this_mu_dist;  // for tabular angular data
  angle_dist::iterator next_mu_dist;
  int prev_mu;                         // 0 or 1, index in mu_intergral
  bool mu_table;    // angular data as a table if true; otherwise Legendre
  Legendre_angle_dist::const_iterator prev_L_coefs;  // for Legendre coefficients
  Legendre_angle_dist::const_iterator next_L_coefs;

  // where the lower and upper eta values hit the quadrature box
  energy_hit_list lower_hits;
  energy_hit_list upper_hits;
 
  long int quad_count;  // number of 2d quadratures
  long int Ein_F_count;  // number of calls to  uncorrelated_F::Ein_F
  long int mu_quad_count;  // number of integrations over mu
  long int mu_F_count;  // number of calls to uncorrelated_F::mu_F

  two_d_interp Ein_interp; // interpolation rule for incident energy
  Interp_Type Eout_interp; // interpolation rule for outgoing energy
  Interp_Type mu_interp; // interpolation rule for direction cosine

  inline uncorrelated_param( ): quad_count( 0 ), Ein_F_count( 0 ),
				mu_quad_count( 0 ), mu_F_count( 0 ) {}
  inline ~uncorrelated_param( ) {}

  //! Sets up the data for interpolation with respect to incident energy
  void set_data( );

  //! Sets up the data for cumulative-points interpolation in incident energy
  void setup_Ein_cum_prob( );

  //! Interpolates to set up for the integration under unitbase interplation
  //! \param E_in energy of incident particle
  void interp_data_ubase( double E_in );

  //! Interpolates to set up for the integration under histogram interplation
  //! \param E_in energy of incident particle
  void interp_data_flat( double E_in );

  //! Go to the next set of (E_out, probability) pairs for cumulative
  //! points interpolation in incident energy.
  //! Returns true if we are finished with the data.
  bool next_Eout_cum_prob( );

};

// -------------- functions to integrate -------------------
namespace uncorrelated_F
{
  // ---------------- mu_F ------------------
  //! Function for the quadrature over mu
  //! \param mu the direction cosine
  //! \param param the function parameters
  //! \param value the value of the integrand, a set of Legendre coefficients
  void mu_F( double mu_in, QuadParamBase *param, coef_vector *value );

  // ---------------- Ein_F ------------------
  //! Function for the quadrature over incident energy
  //! \param E_in the energy of the incident particle
  //! \param param the function parameters
  //! \param value the value of the integrand, a set of Legendre coefficients
  void Ein_F( double E_in, QuadParamBase *param, coef_vector *value );
}

#endif
