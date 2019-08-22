/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 381 $
 * $Date: 2013-06-14 (Fri, Jun 14, 2013) $
 * $Author: hedstrom $
 * $Id: cm_Legendre.hpp 1 2013-06-14 hedstrom $
 *
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// declaration of the classes used to handle Legendre expansions of energy 
// probability density given in the center-of-mass frame

#ifndef CM_LEGENDRE_ENERGY_DIST
#define CM_LEGENDRE_ENERGY_DIST

#include "Legendre_data.hpp"
#include "Ecm_Elab_geom.hpp"
#include "transfer.hpp"
#include "cumulative_points.hpp"

class cm_Legendre_Ein_param;  // forward declaration

//! Class for parameters for the 2-d quadrature over cm cosine and Eout_cm
// ---------------- class cm_Legendre_Ecm_param ------------------
class cm_Legendre_Ecm_param : public Ecm_Elab_Ecm_param
{
public:
  Quadrature_Method mu_quad_method;
  Interp_Type Eout_interp;
  long int mu_F_count;  // number of calls to cm_Legendre_F::mu_F

  // the data entries for this incident energy
  Legendre_data_range *this_data;

  inline cm_Legendre_Ecm_param( ): mu_F_count ( 0 ) {}
  inline ~cm_Legendre_Ecm_param( ) {}
};

//! Class for parameters for the 1-d quadrature over cm cosine
// ---------------- class cm_Legendre_mu_param ------------------
class cm_Legendre_mu_param : public Ecm_Elab_mu_param
{
public:
  Legendre_coefs this_data;  // Legendre data for this outgoing energy

  inline cm_Legendre_mu_param( ) {}
  inline ~cm_Legendre_mu_param( ) {}

  //! Sets up the data for this incident energy and outgoing energy
  //! \param Ein energy of incident particle
  //! \param Eoutcm center-of-mass energy of outgoing particle
  //! \param Ecm_param parameters for integration over outgoing energy
  void setup( double Ein, double Eoutcm, const cm_Legendre_Ecm_param& Ecm_param );
};

//! Class for one energy distribution
//--------------- class cm_Legendre_vector ----------------
class cm_Legendre_vector : public Legendre_list_base
{
private:

public:
  // cumulative probabilities for cumulative points interpolation
  cumulative_prob_list cum_prob;

  inline cm_Legendre_vector( ) {}

  inline ~cm_Legendre_vector( ) {}

  //! Forms the list of cumulative probabilities
  void form_cum_prob( );

  //! Reads the Legendre coefficients of the energy probability density
  //! \param infile input file
  //! \param num_Eout number of outgoing energies for this incident energy
  void read_coef( data_parser& infile, int num_Eout );

  //! Reads the isotropic energy probability density
  //! \param infile input file
  //! \param num_Eout number of outgoing energies for this incident energy
  void read_order_zero( data_parser& infile, int num_Eout );

};

//! Class for energy distributions
//--------------- class cm_Legendre ----------------
class cm_Legendre : public list< cm_Legendre_vector >
{
private:
  //! The smallest incident energy for cross section, multiplicity,
  //! model weight, flux weight, and energy groups
  double E_first;

  int first_Ein;  // index of the left-hand end of the first significant energy bin
  int last_Ein;  // index of the right-hand end of the last significant energy bin

  //! Sets up the map from center-of-mass to laboratory coordinates
  void setup_map( );

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

  //! Loops through the center-of-mass data for a pair of incident energies
  //! \param transfer: the entries in the transfer matrix get updated
  //! \param Ein_param the quadrature parameters
  void Eout_data_ladder( T_matrix& transfer, cm_Legendre_Ein_param *Ein_param );

  //! Loops through the outgoing energy bins for given cm_Eout data
  //! \param transfer: the entries in the transfer matrix get updated
  //! \param Ein_param the quadrature parameters
  void lab_Eout_ladder( T_matrix& transfer, cm_Legendre_Ein_param *Ein_param );

  //! Integrates over one E-E' box
  //! \param transfer: the entries in the transfer matrix get updated
  //! \param Eout_count identifies the matrix entry to update
  //! \param Ein_param the quadrature parameters
  void one_Ebox( T_matrix& transfer, int Eout_count, cm_Legendre_Ein_param *Ein_param );

  //! Go to the next pair of incident energies.  Returns "true" when finished.
  //! \param  E_in the next incident energy
  //! \param Ein_param the quadrature parameters
  bool next_ladder( double E_in, cm_Legendre_Ein_param *Ein_param );

  //! Increments the transfer matrix
  //! \param transfer: the entries in the transfer matrix get updated
  //! \param Eout_count identifies the matrix entry to update
  //! \param Ein_param the quadrature parameters
  void update_T( T_matrix &transfer, int Eout_count, cm_Legendre_Ein_param *Ein_param );

  //! Initializes the quadrature parameters
  //! \param Ein_param the quadrature parameters
  void setup_param( cm_Legendre_Ein_param *Ein_param );

  //! Sets the range of incident energies for this intergration
  //! \param Ein_param the quadrature parameters
  void set_Ein_range( cm_Legendre_Ein_param *Ein_param );

public:
  particleInfo particles;
  map_cm_lab map;

  two_d_interp Ein_interp;  // interpolation rule for incident energy
  Interp_Type Eout_interp; // interpolation rule for outgoing energy
  int order;  // the Legendre order of the output

  cm_Legendre( ): Eout_interp( NOTSET ) {}

  ~cm_Legendre( ) {}

  //! Reads the Legendre data
  //! \param infile input file
  //! \param num_Ein number of incident energies for this reaction
  void read_data( data_parser& infile, int num_Ein );

  //! Reads isotropic data
  //! \param infile input file
  //! \param num_Ein number of incident energies for this reaction
  void read_isotropic( data_parser& infile, int num_Ein );

  //! Computes the transfer matrix
  //! \param sigma the cross section data
  //! \param multiple the outgoing particle multiplicity data
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param transfer the transfer matrix
  void get_T( const dd_vector& sigma, const dd_vector& multiple, 
    const dd_vector& weight, T_matrix& transfer );

  // Prints the lists for debugging
  void print( );

};

// ---------------- class cm_Legendre_Ein_param ------------------
class cm_Legendre_Ein_param : public Ecm_Elab_Ein_param
{
private:
  //! Sets up the data for unit-base interpolation in incident energy
  void setup_Ein_ubase( );

  //! Sets up the data for cumulative points interpolation in incident energy
  void setup_Ein_cum_prob( );

  //! Sets (Eout, probability) data to the lower zero cumulative probability
  void setup_low_A( );

  //! Interpolates (Eout, probability) data to the higher common cumulative probability
  //! \param higher_A the higher cumulative probability value
  void setup_high_A( double higher_A );

public:
  // where the lower and upper eta values hit the quadrature box
  Vcm_Vlab_hit_list lower_hits;
  Vcm_Vlab_hit_list upper_hits;
 
  // Where we are in the incident energy data
  cm_Legendre::const_iterator left_Ein_data;
  cm_Legendre::const_iterator right_Ein_data;

  // Where we are in the outgoing energy data (unit base)
  cm_Legendre_vector::const_iterator left_Eout_data;
  cm_Legendre_vector::const_iterator right_Eout_data;
  cm_Legendre_vector::const_iterator next_left_Eout_data;
  cm_Legendre_vector::const_iterator next_right_Eout_data;

  // current unit-base outgoing energy range
  double prev_data_Eout;
  double next_data_Eout;

  cumulative_prob_list::const_iterator left_cum_prob;  // cumulative probability at lower Ein and lower Eout
  cumulative_prob_list::const_iterator next_left_cum_prob;  // cumulative probability at lower Ein and higher Eout
  cumulative_prob_list::const_iterator right_cum_prob;  // cumulative probability at higher Ein and lower Eout
  cumulative_prob_list::const_iterator next_right_cum_prob;  // cumulative probability at higher Ein and higher Eout

  Legendre_data_range Ein0_data;  // Legendre data for incident energy Ein0
  Legendre_data_range Ein1_data;  // Legendre data for incident energy Ein1
  Legendre_data_range mid_data;         // Legendre data for interpolated incident energy

  cm_Legendre_Ecm_param Ecm_params;   // parameters for integration over center-of-mass energy

  long int quad_count;  // number of 3d quadratures
  long int Ein_F_count;  // number of calls to cm_Legendre_F::Ein_F
  long int Ecm_F_count;  // number of calls to cm_Legendre_F::Ecm_F
  int Vcm_hit_count;     // number of calls to quad_F::integrate in cm_Legendre_F::Ein_F

  Quadrature_Method Eout_quad_method;  // for integration over outgoing energy
  Quadrature_Method mu_quad_method;   // for integration over cosine
  two_d_interp Ein_interp;
  Interp_Type Eout_interp;

  cm_Legendre_Ein_param( ): quad_count( 0 ), Ein_F_count ( 0 ), Ecm_F_count ( 0 ),
			    Vcm_hit_count( 0 ) {}
  ~cm_Legendre_Ein_param( ) {}

  //! Sets up the unit-base Legendre data for interpolation with respect to incident energy
  void set_data( );

//! Class for parameters for the 3-d quadrature over Ein, cm cosine, and Eout_cm
  //! Sets up the parameters for integration over center-of-mass outgoing energy
  void set_Ecm_param( );

  //! Do unit-base interpolation and return the physical Legendre coefficients
  //! \param E_in energy of incident particle
  void ubase_interpolate( double E_in );

  //! Starts one staircase of the Eout_cm data
  //! \param Ein_param the quadrature parameters
  void start_Eout_cm( );

  //! Gets the next interval of unit-base outgoing energy
  //! \param Ein_param the quadrature parameters
  bool next_Ecm( );

  //! Gets the next interval of outgoing energy for cumulative points interpolation
  bool CP_next_Ecm( );
};

// **************** Functions to integrate *********************
namespace cm_Legendre_F
{
  //! Function for the 1-d quadrature over cm cosine
  //! \param mu the center-of-mass direction cosine of the outgoing particle
  //! \param mu_quad_param the function parameters
  //! \param value the value of the integrand, a set of Legendre coefficients
  void mu_F( double mu, QuadParamBase *mu_quad_param, coef_vector *value );

  //! Function for the 2-d quadrature over cm cosine and Eout_cm
  //! \param Eout_cm the center-of-mass energy of the outgoing particle
  //! \param Ecm_quad_param the function parameters
  //! \param value the value of the integrand, a set of Legendre coefficients
  void Ecm_F( double Eout_cm, QuadParamBase *Ecm_quad_param, coef_vector *value );

  //! Function for the 3-d quadrature over E_in, and Eout_cm and cm cosine
  //! The value of Ein_F is itself an integral over Eout_cm and cm cosine.
  //! \param E_in the energy of the incident particle
  //! \param Ein_quad_param the function parameters
  //! \param value the value of the integrand, a set of Legendre coefficients
  void Ein_F( double E_in, QuadParamBase *Ein_quad_param, coef_vector *value );
}
#endif
