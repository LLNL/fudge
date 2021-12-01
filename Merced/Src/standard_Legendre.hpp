/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2011-03-07 (Mon, Mar 7, 2011) $
 * $Author: hedstrom $
 * $Id: standard_Legendre.hpp 1 2011-03-07 hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// declaration of the classes used to handle Legendre expansions of energy probability density



#ifndef standard_LEGENDRE_ENERGY_DIST
#define standard_LEGENDRE_ENERGY_DIST

#include "param_base.hpp"
#include "math_util.hpp"
#include "transfer.hpp"
#include "box_geom.hpp"  // for Box::energy_hit_list
#include "cumulative_points.hpp"

namespace StdLg
{
class standard_Legendre_param;

//! Class for one energy distribution
//--------------- class standard_Legendre_vector ----------------
class standard_Legendre_vector : public Lgdata::Legendre_list_base
{
private:
  //! Appends a copy of the data to the list
  void append_data( double E_out, const Lgdata::Legendre_coefs &to_copy );

public:
  // cumulative probabilities for cumulative points interpolation
  Cum::cumulative_prob_list cum_prob;

  inline standard_Legendre_vector( ) {}

  inline ~standard_Legendre_vector( ) {}

  //! Reads the Legendre coefficients of the energy probability density
  //! \param infile input file
  //! \param num_Eout number of outgoing energies for this data
  //! \param max_order themaximum Legendre order in the data
  void read_coef( Dpar::data_parser& infile, int num_Eout, int max_order );

  //! Forms the list of cumulative probabilities
  void form_cum_prob( );

  //! Copies a vector
  //! \param vector_from the vector to copy
  void copy( const StdLg::standard_Legendre_vector& vector_from );

  //! Copies a vector with truncation
  //! \param vector_from the vector to copy
  //! \param min_E if necessary, truncate to this lowest energy
  //! \param max_E if necessary, truncate to this highest energy
  bool truncate_copy( const StdLg::standard_Legendre_vector& vector_from,
     double min_E, double max_E );

  //! Copies a vector with extrapolation
  //! \param vector_from the vector to copy
  //! \param min_E if necessary, extrapolate as zero to this lowest energy
  //! \param max_E if necessary, extrapolate as zero to this highest energy
  void extrapolate_copy( const StdLg::standard_Legendre_vector& vector_from,
     double min_E, double max_E );
};

//! Class for energy distributions
//--------------- class standard_Legendre ----------------
class standard_Legendre : public std::list< StdLg::standard_Legendre_vector >
{
private:
  //! The smallest incident energy for cross section, multiplicity,
  //! model weight, flux weight, and energy groups
  double E_first;

  //! the range of nonzero incident energy bins
  int first_Ein;
  int last_Ein;

  //! Gets the range of nontrivial incident energy bins; computes E_first, first_Ein and last_Ein
  //! returns true if the threshold is too high for the energy bins
  //! \param sigma_ the cross section data
  //! \param mult_ the outgoing particle multiplicity data
  //! \param weight_ the weighting to apply to the transfer matrix entries
  //! \param e_flux_ the initial approximation to apply to the particle flux
  //! \param Ein_groups the boundaries of the incident energy groups
  bool get_Ein_range( const Ddvec::dd_vector& sigma_, const Ddvec::dd_vector& mult_,
                      const Ddvec::dd_vector& weight_,
		      const Lgdata::Flux_List& e_flux_, const Egp::Energy_groups& Ein_groups );

  //! Initializes the quadrature parameters
  void setup_data( StdLg::standard_Legendre_param *Ein_param );

  //! Loops through the energy data
  //! \param transfer the transfer matrix to compute
  //! \param Ein_param the quadrature parameters
  void E_data_ladder( Trf::T_matrix& transfer,
		      StdLg::standard_Legendre_param *Ein_param );

  //! Loops through the outgoing energy bins
  //! \param transfer the transfer matrix to compute
  void Eout_ladder( Trf::T_matrix& transfer,
		    StdLg::standard_Legendre_param *Ein_param );

  //! Integrates over one E-E' box
  //! \param transfer the transfer matrix to compute
  //! \param Eout_count identifies the matrix entry to update
  //! \param Ein_param the quadrature parameters
  void one_Ebox( Trf::T_matrix& transfer, int Eout_count,
		 StdLg::standard_Legendre_param *Ein_param );

  //! Go to the next pair of incident energies.  Returns "true" when finished.
  //! \param  E_in the next incident energy
  //! \param Ein_param the quadrature parameters
  bool next_ladder( double E_in, StdLg::standard_Legendre_param *Ein_param );

  //! Increments the transfer matrix
  //! \param transfer the transfer matrix to compute
  //! \param Eout_count identifies the matrix entry to update
  //! \param Ein_param the quadrature parameters
  void update_T( Trf::T_matrix &transfer, int Eout_count,
		 StdLg::standard_Legendre_param *Ein_param );

  //! Sets the range of incident energies for this intergration
  //! \param Ein_bin the incident energy bin
  //! \param Ein_param the quadrature parameters
  void set_Ein_range( int Ein_bin, StdLg::standard_Legendre_param &Ein_param );

public:
  Terp::two_d_interp Ein_interp;  // interpolation rule for incident energy
  Terp::Interp_Type Eout_interp; // interpolation rule for outgoing energy
  int order;  // the Legendre order of the output

  standard_Legendre( );

  ~standard_Legendre( );

  //! Reads the Legendre data
  //! \param infile input file
  //! \param num_Ein number of incident energies
  void read_data( Dpar::data_parser& infile, int num_Ein );

  //! Computes the transfer matrix
  //! \param sigma the cross section data
  //! \param multiople the outgoing particle multiplicity data
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param transfer the transfer matrix
  void get_T( const Ddvec::dd_vector& sigma, const Ddvec::dd_vector& multiple, 
    const Ddvec::dd_vector& weight, Trf::T_matrix& transfer );

  // Prints the lists for debugging
  void print( );

};

//! Class for parameters for 2d integration
//--------------- class standard_Legendre_param ----------------
class standard_Legendre_param : public Pbase::param_base
{
private:
  //! Extrapolated data used with linlin interpolation of incident energy
  StdLg::standard_Legendre_vector low_linlin;  // (Eout, probability) at lower incident energy
  StdLg::standard_Legendre_vector high_linlin;  // (Eout, probability) at higher incident energy

  //! The current data
  double left_Ein;  // lower incident energy
  double right_Ein;  // higher incident energy

  //! Interpolates (Eout, probability) data to the lower common Eout value
  //! \param lower_Eout the lower end of the Eout range
  void common_low_Eout( double lower_Eout );

  //! Interpolates (Eout, probability) data to the higher common Eout value
  //! \param higher_Eout the upper end of the Eout range
  void common_high_Eout( double higher_Eout );

  //! Sets (Eout, probability) data to the lower zero cumulative probability
  void setup_low_A( );

  //! Interpolates (Eout, probability) data to the higher common cumulative probability
  void setup_high_A( double higher_A );

  //! Sets up the data for unit-base interpolation in incident energy
  void setup_Ein_ubase( );

  //! Sets up the data for cumulative-points interpolation in incident energy
  void setup_Ein_cum_prob( );

  //! Sets up the data for linlin interpolation in incident energy
  void setup_Ein_linlin( );

  //! Sets up the data for histogram interpolation in incident energy
  void setup_Ein_flat( );

public:
  // (Eout_cm, probability) data interpolated to common Eout_cm values
  Lgdata::Legendre_data_range Ein0_data;  // Legendre data for lower incident energy
  Lgdata::Legendre_data_range Ein1_data;  // Legendre data for higher incident energy

  Ddvec::unit_base_map mid_ubase_map;  // for interpolated data

  // where the lower and upper eta values hit the quadrature box
  Box::energy_hit_list lower_hits;
  Box::energy_hit_list upper_hits;
 
  // where we are in the data
  StdLg::standard_Legendre::iterator this_Ein;
  StdLg::standard_Legendre::iterator next_Ein;

  // pointers to the data entries
  StdLg::standard_Legendre_vector::const_iterator left_ptr;  //(lower Eout, probability), low Ein
  StdLg::standard_Legendre_vector::const_iterator next_left_ptr; //(higher Eout, probability), low Ein
  StdLg::standard_Legendre_vector::const_iterator last_left_ptr; // end of (Eout, probability),low Ein

  StdLg::standard_Legendre_vector::const_iterator right_ptr;  //(lower Eout, probability), high Ein
  StdLg::standard_Legendre_vector::const_iterator next_right_ptr; //(higher Eout, probability), high Ein
  StdLg::standard_Legendre_vector::const_iterator last_right_ptr; //end of (Eout, probability), high Ein

  Cum::cumulative_prob_list::const_iterator left_cum_prob;  // cumulative probability at lower Ein and lower Eout
  Cum::cumulative_prob_list::const_iterator next_left_cum_prob;  // cumulative probability at lower Ein and higher Eout
  Cum::cumulative_prob_list::const_iterator right_cum_prob;  // cumulative probability at higher Ein and lower Eout
  Cum::cumulative_prob_list::const_iterator next_right_cum_prob;  // cumulative probability at higher Ein and higher Eout

  //! Interpolated data at intermediate incident energy
  Lgdata::Legendre_coefs mid_lower_Eout;  // at lower unit-base outgoing energy for the data
  Lgdata::Legendre_coefs mid_upper_Eout;  // at higher unit-base outgoing energy for the data
  Lgdata::Legendre_coefs use_prev_Eout;  // at lower unit-base outgoing energy for integration range
  Lgdata::Legendre_coefs use_next_Eout;  // at higher unit-base outgoing energy for integration range

  Ddvec::dd_entry Eout_0_range; // values of physical Eout for left_data and next_left_data
  Ddvec::dd_entry Eout_1_range; // values of physical Eout for right_data and next_right_data
  Ddvec::dd_entry Eout_range; // values of physical Eout for interpolated data
  Terp::two_d_interp Ein_interp; // interpolation rule for incident energy
  Terp::Interp_Type Eout_interp; // interpolation rule for outgoing energy

  long int quad_count;  // number of 2-d quadratures
  long int Ein_F_count;  // number of calls to standard_Legendre_F::Ein_F

  //! Constructor
  inline standard_Legendre_param( ): quad_count(0), Ein_F_count(0) {}

  //! Destructor
  inline ~standard_Legendre_param( ) {}

  //! Allocates space
  //! \param Order the Legendre order of the output transfer matrix
  void initialize( int Order );

  //! Initializes the data pointers for one incident energy range
  void reset_start( );

  //! Increments the data pointers for one incident energy range
  bool get_next_Eout( );

  //! Interpolates unitbase between two incident energies
  //! Returns true if the interpolation is OK
  //! \param E_in the energy of the incident particle
  bool unitbase_interpolate( double Ein );

  //! Interpolates direct linlin between two incident energies
  //! \param E_in the energy of the incident particle
  void direct_linlin_interpolate( double Ein );

  //! Interpolates histogram between two incident energies
  void flat_interpolate( );
};

} // end of namespace StdLg

namespace standard_Legendre_F
{
  // ---------------- standard_Legendre_F::Ein_F ------------------
  //! Function for the quadrature over incident energy
  //! Returns true if the interpolation is OK
  //! \param E_in the energy of the incident particle
  //! \param e_quad_param the function parameters
  //! \param value the value of the integrand, a set of Legendre coefficients
  bool Ein_F( double E_in, Qparam::QuadParamBase *e_quad_param,
    Coef::coef_vector *value );
}

#endif
