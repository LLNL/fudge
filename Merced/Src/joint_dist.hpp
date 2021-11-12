/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2014-04-22 (Tue, Apr 22, 2014) $
 * $Author: hedstrom $
 * $Id: joint_dist.hpp 1 2014-04-22Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// header for the classes used on joint energy-angle distributions

#ifndef JOINT_DIST_CLASS
#define JOINT_DIST_CLASS

#include "joint_dist_data.hpp"
#include "param_base.hpp"
#include "transfer.hpp"
#include "angle_dist.hpp"
#include "box_geom.hpp"
#include "cumulative_points.hpp"

namespace Jdist
{
class joint_dist_param;  // forward declaration

// ----------- class joint_dist_hits -----------------
//! Class for the list of intersections of a unit-base curve with an integration box
class joint_dist_hits : public Box::energy_hit_list
{
private:

public:
  // ! where we are in the data
  double Ein_0;  // lower incident energy
  double Ein_1;  // higher incident energy
  Jdata::E_mu_P_data *Ein0_data;  // lower Ein
  Jdata::E_mu_P_data *Ein1_data;  // higher Ein

  //! Constructor
  inline joint_dist_hits( )
  {}

  //! Default destructor
  inline ~joint_dist_hits( )
  {}

  //! Sets the range of integration over incident energy
  //! \param A the lower incident energy for this integral
  //! \param B the higher incident energy for this integral
  void set_incident_range( double A, double B );

  //! Gets the physical outgoing energies at the ends of the incident energy bin
  void get_phys_Eout( );
};

// ----------- class joint_Eout_param -----------------
//! parameters for 1-d integration in outgoing energy
class joint_Eout_param: public Pbase::param_base
{
public:
  //! Rule for interpolation in outgoing energy
  Terp::Interp_Type Eout_interp;

  //! (Eout, probability) data at the lower outgoing energy
  Jdata::E_mu_P_data Eout0_data;

  //! (Eout, probability) data at the higher outgoing energy
  Jdata::E_mu_P_data Eout1_data;

  //! Default constructor
  inline joint_Eout_param( )
  {}

  //! Default destructor
  inline ~joint_Eout_param( )
  {}

};

// ----------- class joint_mu_param -----------------
//! parameters for 2-d integration over mu and E_out
class joint_mu_param: public Pbase::param_base
{
private:
 
public:
  //! (Eout, probability) data at this incident energy
  Jdata::current_data this_data;

  //! parameters for integration over outgoing energy
  Jdist::joint_Eout_param Eout_params;

  //! unit-base maps of the outgoing energy
  Ddvec::unit_base_map mu0_ubase_map;  // outgoing energy range at lower_mu
  Ddvec::unit_base_map mu1_ubase_map;  // outgoing energy range at upper_mu
  Ddvec::unit_base_map mid_ubase_map;  // outgoing energy range at the current mu

  Qmeth::Quadrature_Rule Eout_quad_rule;  // quadrature rule for outgoing energy
  long int Eout_F_count;  // number of calls to joint_dist_F::Eout_F

  //! the outgoing energy bin boundaries
  std::vector< double >::const_iterator Eout_bottom;
  std::vector< double >::const_iterator Eout_top;

  //! Default constructor
  inline joint_mu_param( ): Eout_F_count ( 0 ) 
  {}

  //! Default destructor
  inline ~joint_mu_param( )
  {}

  //! determines the geometry for 2d integration over outgoing cosine and energy
  //! returns true if the geometry makes sense
  //! \param this_data (Eout, probability) data at the current incident energy
  bool geometry( const Jdata::current_data &this_data );
};

// ----------- class one_mu -----------------
//! Class for the energy distribution at given incident energy and angle
class one_mu : public Ddvec::dd_vector
{
private:

public:
  // for cumulative-points interpolation
  Cum::cumulative_prob_list cum_prob;

  // for unit-base interpolation
  Ddvec::unit_base_map ubase_map;

  //! probability of the cosine for ENDL data
  double mu_Prob;

  //! Default constructor
  inline one_mu( )
  {}

  //! Default destructor
  inline ~one_mu( )
  {}

  //! Use the direction cosine, mu, as the tag
  inline double get_mu( ) const
  {
    return get_tag( );
  }

  //! Sets the direction cosine tag
  //! \param mu the direction cosine of the outgoing particle
  inline void set_mu( double mu )
  {
    set_tag( mu );
  }

  //! We need to make a copy for the one_joint_dist_quad::copy routine
  //! \param to_copy the vector of pairs ( outgoing energy, probability ) to copy
  void copy( const Jdist::one_mu &to_copy );

  //! Sets up cumulative probabilities for outgoing energy
  void form_cum_prob( );

  //! Sets up Jdata::E_mu_P_data at unit-base outgoing energy UB_Eout
  //! \param UB_Eout an intermediate unit-base outgoing energy
  //! \param phys_mu the physical direction cosine
  //! \param prev_data ( unit-base outgoing energy, probability ) at lower E_out
  //! \param next_data ( unit-base outgoing energy, probability ) at higher E_out
  //! \param mid_data Jdata::E_mu_P_data at unit-base outgoing energy UB_Eout
  void set_E_mu_P_data( double UB_Eout, double phys_mu,
    Jdist::one_mu::const_iterator prev_data,
    Jdist::one_mu::const_iterator next_data, Jdata::E_mu_P_data *mid_data );
};

// ----------- class one_joint_dist -----------------
//! Class for the joint energy-angle distribution at given incident energy
class one_joint_dist : public std::list< Jdist::one_mu >
{
private:
  double tag_;

public:
  //! for mapping direction cosines to 0 <= mu <= 1
  Ddvec::unit_base_map mu_ubase_map;

  Terp::two_d_interp Ein_interp;
  Terp::Interp_Type Eout_interp;
  Terp::two_d_interp mu_interp;

  //! is this ENDL data?
  bool ENDL_data;

  //! Default constructor
  inline one_joint_dist( )
  {}

  //! Default destructor
  inline ~one_joint_dist( )
  {}

  //! Use E_in as the tag
  inline double get_E_in( ) const
  {
    return tag_;
  }

  //! Sets the tag
  //! \param E_in the energy of the incident particle
  inline void set_E_in( double E_in )
  {
    tag_ = E_in;
  }

  //! Makes a copy at the threshold, used by the ENDL_kludge routine
  //! \param to_copy the vectors of pairs ( outgoing energy, probability ) to copy
  void copy( const Jdist::one_joint_dist &to_copy );

  //! Sets mu_Prob for ENDL double-differential data
  //! \param angles vector of pairs ( outgoing direction cosine, probability )
  void set_mu_Prob( Adist::angle_dist::iterator &angles );

  //! Converts ENDF double-differential data to our format
  void convert_ENDF( );

  //! Converts ENDL double-differential data to our format
  //! \param angles vector of pairs ( outgoing direction cosine, probability )
  void convert_ENDL( Adist::angle_dist::iterator &angles );

  //! Maps the direction cosines to 0 <= mu <= 1
  //! \param angles vector of pairs ( outgoing direction cosine, probability )
  void mu_to_unit_base( Adist::angle_dist::iterator &angles );

  //! Maps the outgoing energies to 0 <= Eout <= 1
  void Eout_to_unit_base( );

  //! Sets up cumulative probabilities for outgoing energy
  void form_Eout_cum_prob( );

  //! Insert energy probability densities for an intermediate mu
  //! This routine is used when ENDL data has (mu, probability) data
  //! at a mu value which is missing from the (mu, Eout, probability) table
  //! This version is for unit-base separate interpolation of
  //! $P( E' \mid E, \mu )$ and $P( \mu \mid E )$.
  //! \param angle_ptr, the intermediate direction cosine and its probability density
  //! \param prev_joint_mu (Eout, probability) data at a lower mu value
  //! \param joint_mu (Eout, probability) data at a higher mu value
  void interpolate_mu_UB( Ddvec::dd_vector::const_iterator &angle_ptr,
    Jdist::one_joint_dist::iterator &prev_joint_mu,
    Jdist::one_joint_dist::iterator &joint_mu );

  //! Insert energy probability densities for an intermediate mu
  //! This routine is used when ENDL data has (mu, probability) data
  //! at a mu value which is missing from the (mu, Eout, probability) table
  //! This version is for cumulative-points separate interpolation of
  //! $P( E' \mid E, \mu )$ and $P( \mu \mid E )$.
  //! \param angle_ptr, the intermediate direction cosine and its probability density
  //! \param prev_joint_mu (Eout, probability) data at a lower mu value
  //! \param next_joint_mu (Eout, probability) data at a higher mu value
  void interpolate_mu_CP( Ddvec::dd_vector::const_iterator &angle_ptr,
    Jdist::one_joint_dist::iterator &prev_joint_mu,
    Jdist::one_joint_dist::iterator &next_joint_mu );

  //! Sets the mu data for cumulative-points interpolation
  //! \param this_mu, interpolate data to this mu
  //! \param low_mu, the outgoing energy data at the lower mu value
  //! \param high_mu, the outgoing energy data at the higher mu value
  //! \param Ein_data, set the mu values in this current data
  //! \param level, 0 for the lower mu value, 1 for the higher
  void set_mu_data( double this_mu,
		 Jdist::one_joint_dist::const_iterator low_mu,
		  Jdist::one_joint_dist::const_iterator high_mu,
		      Jdata::current_data *Ein_data, int mu_level );
};

// ----------- class joint_dist -----------------
//! Class for joint energy-angle distributions
class joint_dist : public std::list< Jdist::one_joint_dist >
{
private:

  int first_Ein;  // index of the left-hand end of the first significant energy bin
  int last_Ein;  // index of the right-hand end of the last significant energy bin

  //! ENDL has the convention that the distribution is assumed independent of energy at low incident energies.
  void ENDL_kludge( );

  //! Converts ENDL double-differential data to our format
  void convert_ENDL( );

  //! Converts ENDF double-differential data to our format
  void convert_ENDF( );

  //! Maps the direction cosines to 0 <= mu <= 1
  void mu_to_unit_base( );

  //! Maps the outgoing energies to 0 <= Eout <= 1
  void Eout_to_unit_base( );

  //! Sets up cumulative probabilities for outgoing energy
  void form_Eout_cum_prob( );

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
  //! \param Ein_param the qudrature parameters
  void setup_param( Jdist::joint_dist_param *Ein_param );

  //! Sets the range of incident energies for this intergration
  //! \param Ein_param the qudrature parameters
  void set_Ein_range( Jdist::joint_dist_param *Ein_param );

  //! Handles the ( cosine, Eout, probability ) data for one pair of incident energies
  //! \param transfer the transfer matrix
  //! \param Ein_param the qudrature parameters
  void mu_data_ladder( Trf::T_matrix& transfer,
		       Jdist::joint_dist_param *Ein_param );

  //! Go to the next pair of incident energies.  Returns "true" when finished.
  //! \param E_in energy of the incident particle
  //! \param Ein_param the qudrature parameters
  bool next_Ein_pair( double E_in, Jdist::joint_dist_param *Ein_param );

  //! Initializes the pointers to the ( Eout, probability ) data for current cosines
  //! \param Ein_param the qudrature parameters
  void start_mu_data( Jdist::joint_dist_param *Ein_param );

  //! Go to the next pairs of direction cosines.  Returns "true" when finished.
  //! \param Ein_param the qudrature parameters
  bool next_mu_pairs( Jdist::joint_dist_param *Ein_param );

  //! Adds to the transfer matrix for all E_out bins for given Ein and mu ranges
  //! \param transfer the transfer matrix
  //! \param Ein_param the qudrature parameters
  void Eout_data_ladder( Trf::T_matrix& transfer,
			 Jdist::joint_dist_param *Ein_param );

  //! Starts one staircase of the Eout data for unit-base interpolation
  //! \param Ein_param the qudrature parameters
  void start_Eout_data_UB( Jdist::joint_dist_param *Ein_param );

  //! Starts one staircase of the Eout data for cumulative-points interpolation
  //! \param Ein_param the qudrature parameters
  void start_Eout_data_CP( Jdist::joint_dist_param *Ein_param );

  //! Adds to the transfer matrix for the current data
  //! \param transfer the transfer matrix
  //! \param Eout_count the current row of the transfer matrix
  //! \param Ein_param the qudrature parameters
  void one_box( Trf::T_matrix& transfer, int Eout_count,
		Jdist::joint_dist_param *Ein_param );

  //! Increments the transfer matrix
  //! \param transfer the transfer matrix
  //! \param Eout_count the current row of the transfer matrix
  //! \param Ein0_orig the lower incident energy determined by the probability data
  //! \param Ein1_orig the higher incident energy determined by the probability data
  //! \param Ein_param the qudrature parameters
  void update_T( Trf::T_matrix& transfer, int Eout_count,
                 double Ein0_orig, double Ein1_orig,
		 Jdist::joint_dist_param *Ein_param );


public:
  //! The I=1 angular probability densities
  Adist::angle_dist angle_data;

  Terp::two_d_interp Ein_interp; // interpolation rule for incident energy
  Terp::Interp_Type Eout_interp; // interpolation rule for outgoing energy
  Terp::two_d_interp mu_interp; // interpolation rule for direction cosine

  bool ENDL_data;

  //! Default constructor
  joint_dist( ): Eout_interp( Terp::LINLIN ), ENDL_data( false )
  {}

  //! Default destructor
  inline ~joint_dist( )
  {}

  //! Reads the python data
  //! \param infile input file
  //! \param num_Ein number of incident energies for this reaction
  void read_data( Dpar::data_parser &inFile, int num_Ein );

  //! Calculates the transfer matrix for this particle.
  //! \param sigma the cross section data
  //! \param mult the outgoing particle multiplicity data
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param transfer the transfer matrix
  void get_T( const Ddvec::dd_vector& sigma, const Ddvec::dd_vector& mult, 
	      const Ddvec::dd_vector& weight, Trf::T_matrix& transfer );
};

// ----------- class joint_dist_param -----------------
//! parameters for 3-d integration
class joint_dist_param: public Pbase::param_base
{
private:
  //! Increments cumulative probability data for one incident energy and one mu
  //! Returns true of we run out of data
  //! \param low_cum, the current lower cumulative probability
  //! \param prev_CP, current CP data at lower outgoing energy
  //! \param next_CP, current CP data at higher outgoing energy
  //! \param end_CP, end of CP data
  bool next_cum_data( double low_cum,
		      Cum::cumulative_prob_list::iterator &prev_CP,
		      Cum::cumulative_prob_list::iterator &next_CP,
		      Cum::cumulative_prob_list::iterator end_CP );
  
  //! Sets the value of upper_E_cum
  void set_upper_cum( );
  
public:
  // where the lower and upper physical Eout values hit the quadrature box
  Jdist::joint_dist_hits lower_mu0_hits;  // lower mu, lower Eout
  Jdist::joint_dist_hits lower_mu1_hits;  // higher mu, lower Eout
  Jdist::joint_dist_hits upper_mu0_hits;  // lower mu, higher Eout
  Jdist::joint_dist_hits upper_mu1_hits;  // higher mu, higher Eout

  //! where we are in the list, which incident energies
  Jdist::joint_dist::iterator this_Ein_ptr;
  Jdist::joint_dist::iterator next_Ein_ptr;

  //! where we are in the data, which mu values
  Jdist::one_joint_dist::iterator this_Ein_this_mu;  // lower Ein, lower mu
  Jdist::one_joint_dist::iterator this_Ein_next_mu;  // lower Ein, higher mu
  Jdist::one_joint_dist::iterator next_Ein_this_mu;  // higher Ein, lower mu
  Jdist::one_joint_dist::iterator next_Ein_next_mu;  // higher Ein, higher mu

  //! where we are in the (Eout, probability) data
  Jdist::one_mu::const_iterator Ein0_mu0_Eout0;  // lower Ein, lower mu, lower Eout
  Jdist::one_mu::const_iterator Ein0_mu0_Eout1;  // lower Ein, lower mu, higher Eout
  Jdist::one_mu::const_iterator Ein0_mu1_Eout0;  // lower Ein, higher mu, lower Eout
  Jdist::one_mu::const_iterator Ein0_mu1_Eout1;  // lower Ein, higher mu, higher Eout
  Jdist::one_mu::const_iterator Ein1_mu0_Eout0;  // higher Ein, lower mu, lower Eout
  Jdist::one_mu::const_iterator Ein1_mu0_Eout1;  // higher Ein, lower mu, higher Eout
  Jdist::one_mu::const_iterator Ein1_mu1_Eout0;  // higher Ein, higher mu, lower Eout
  Jdist::one_mu::const_iterator Ein1_mu1_Eout1;  // higher Ein, higher mu, higher Eout

  //! The corners of the parallelpiped for 3-d interpolation of data
  //! The range of incident energies is as in the data
  double lower_mu;  // the lower direction cosine
  double upper_mu;  // the higher direction cosine
  double lower_Eout;  // the lower outgoing energy
  double upper_Eout;  // the higher outgoing energy
  double lower_E_cum;  // lower cumulative probability for energy
  double upper_E_cum;  // upper cumulative probability for energy

  Cum::cumulative_prob_list::iterator Ein0_mu0_Eout0_CP;  // cumulative probability at lower Ein, lower mu, and lower Eout
  Cum::cumulative_prob_list::iterator Ein0_mu0_Eout1_CP;  // cumulative probability at lower Ein, lower mu, and higher Eout
  Cum::cumulative_prob_list::iterator Ein0_mu1_Eout0_CP;  // cumulative probability at lower Ein, higher mu, and lower Eout
  Cum::cumulative_prob_list::iterator Ein0_mu1_Eout1_CP;  // cumulative probability at lower Ein, higher mu, and higher Eout
  Cum::cumulative_prob_list::iterator Ein1_mu0_Eout0_CP;  // cumulative probability at higher Ein, lower mu, and lower Eout
  Cum::cumulative_prob_list::iterator Ein1_mu0_Eout1_CP;  // cumulative probability at higher Ein, lower mu, and higher Eout
  Cum::cumulative_prob_list::iterator Ein1_mu1_Eout0_CP;  // cumulative probability at higher Ein, higher mu, and lower Eout
  Cum::cumulative_prob_list::iterator Ein1_mu1_Eout1_CP;  // cumulative probability at higher Ein, higher mu, and higher Eout

  // parameters for integration over outgoing cosine and energy
  Jdist::joint_mu_param mu_params;

  //! (Eout, probability) data at the lower incident energy
  Jdata::current_data Ein0_data;

  //! (Eout, probability) data at the higher incident energy
  Jdata::current_data Ein1_data;

  //! How the physical lower outgoing energy hits the quadrature box
  Box::energy_hit_list lower_2d_hits;

  //! How the physical higher outgoing energy hits the quadrature box
  Box::energy_hit_list upper_2d_hits;
 
  Qmeth::Quadrature_Rule Eout_quad_rule;  // quadrature rule for outgoing energy
  Qmeth::Quadrature_Rule mu_quad_rule;  // quadrature rule for outgoing cosine

  long int quad_count;  // number of 3-d quadratures
  long int Ein_F_count;  // number of calls to joint_dist_F::Ein_F
  long int Eout_F_count;  // number of calls to joint_dist_F::Eout_F
  long int mu_F_count;  // number of calls to joint_dist_F::mu_F

  //! constructor
  joint_dist_param( );

  //! Default destructor
  inline ~joint_dist_param( )
  {}

  //! determines the geometry for 2d integration over outgoing cosine and energy
  //! returns true if the geometry makes sense
  //! \param this_data (Eout, probability) data at the current incident energy
  bool geometry( const Jdata::current_data &this_data );

  //! Interpolates data to common values of unit-base mu and outgoing energy
  //! \param UB_Eout the desired value of unit-base outgoing energy
  //! \param Ein0_mu0_data computed  Jdata::E_mu_P_data at lower incident energy, lower_mu
  //! \param Ein0_mu1_data computed  Jdata::E_mu_P_data at lower incident energy, upper_mu
  //! \param Ein1_mu0_data computed  Jdata::E_mu_P_data at higher incident energy, lower_mu
  //! \param Ein1_mu1_data computed  Jdata::E_mu_P_data at higher incident energy, upper_mu
  void common_mu_Eout( double UB_Eout, Jdata::E_mu_P_data *Ein0_mu0_data,
    Jdata::E_mu_P_data *Ein0_mu1_data, Jdata::E_mu_P_data *Ein1_mu0_data,
    Jdata::E_mu_P_data *Ein1_mu1_data );

  //! go to the next sets of (E_out, probability) pairs, unit base
  bool next_Eout_UB( );

  //! go to the next sets of (E_out, probability) pairs, cumulative points
  bool next_Eout_CP( );

  //! Sets the data for cumulative-points interpolation for one incident energy
  //! Returns 0 if unit-base maps are OK,
  //!         1 if the unit-base map for mu0 is bad,
  //!         2 if the unit-base map for mu1 is bad,
  //!         3 if the unit-base map for mu0 and mu1 are bad
  //! \param mu0_Eout0_CP, cumulative-points data at lower_mu
  //! \param mu1_Eout0_CP, cumulative-points data at higher_mu
  //! \param Ein_data, cumulative-points data interpolated to common mu and Eout
  int set_cum_data( Cum::cumulative_prob_list::const_iterator mu0_Eout0_CP,
     Cum::cumulative_prob_list::const_iterator mu1_Eout0_CP,
		     Jdata::current_data *Ein_data );
};
} // end of namespace Jdist

// ************* functions to integrate ******************
namespace joint_dist_F
{
  // ---------------- joint_dist_F::Eout_F ------------------
  //! Function for the 1-d quadrature over outgoing energy
  //! Returns true if the interpolation is OK
  //! \param Eout the lab-frame energy of the outgoing particle
  //! \param Eout_quad_param the function parameters
  //! \param value the value of the integrand, a set of Legendre coefficients
  bool Eout_F( double Eout, Qparam::QuadParamBase *Eout_quad_param,
	       Coef::coef_vector *value );

  // ---------------- joint_dist_F::mu_F ------------------
  //! Function for the 2-d quadrature over outgoing energy and cosine
  //! Returns true if the interpolation is OK
  //! \param mu the lab-frame direction cosine of the outgoing particle
  //! \param mu_quad_param the function parameters
  //! \param value the value of the integrand, a set of Legendre coefficients
  bool mu_F( double mu, Qparam::QuadParamBase *mu_quad_param,
	     Coef::coef_vector *value );

  // ---------------- joint_dist_F::Ein_F ------------------
  //! Function for the 3-d quadrature over incident energy, cosine, and outgoing energy
  //! Returns true if the interpolation is OK
  //! \param E_in the energy of the incident particle
  //! \param Ein_quad_param the function parameters
  //! \param value the value of the integrand, a set of Legendre coefficients
  bool Ein_F( double E_in, Qparam::QuadParamBase *Ein_quad_param,
	      Coef::coef_vector *value );
}

#endif
