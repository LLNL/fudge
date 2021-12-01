/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: angle_dist.hpp 1 2006-01 19:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */
// define the classes used for ENDL i=1 angular distributions

#ifndef ANGLE_DIST_DEF
#define ANGLE_DIST_DEF

#include "mappings.hpp"
#include "relativistic.hpp"
#include "transfer.hpp"
#include "math_util.hpp"
#include "box_geom.hpp"
#include "param_base.hpp"

namespace Adist
{
class two_body_Ein_param;  // forward declaration

// ---------------- class min_Eout_info ------------------
//! This class holds
//! the negative direction cosine, the minimal outgoing energy in the lab
//! frame, and the corresponding incident energy.
class min_Eout_info
{
public:
  double mu;  // direction cosine in center-of-mass frame (negative)
  double Ein;  // incident energy at the minimum
  double Eout; // minimum outgoing energy in lab frame

  min_Eout_info( ) {}
  ~min_Eout_info( ) {}
};

// ---------------- class min_Eout_info_list ------------------
//! This structure holds a list of
//! the negative direction cosine, the minimal outgoing energy in the lab
//! frame, and the corresponding incident energy.
class min_Eout_info_list : public std::list< Adist::min_Eout_info >
{
public:
  min_Eout_info_list( ) {}
  ~min_Eout_info_list( ) {}
};

//! Class for parameters for the 1-d quadrature based on a pair of (mucm, probability) values.
// ---------------- class two_body_mucm_param ------------------
class two_body_mucm_param : public Ddvec::dd_pair, public Qparam::QuadParamBase
{
public:
  bool use_relativistic;

  Maps::Newton_map_param *Newton_map;
  Rel::relativistic_map relMap;

  inline two_body_mucm_param( ): use_relativistic( false ) {}
  inline ~two_body_mucm_param( ) {}
};

//! Class for the list of intersections of a curve with an integration box
//! In this class the hit_list eta parameter is the center-of-mass direction cosine.
// ---------------- class angle_hit_list ------------------
class angle_hit_list : public Box::hit_list
{
private:
  // Implement the virtual functions
  //! Calculates E_out from E_in, for testing the side of a box
  //! \param E_in energy of the incident particle
  double get_Eout( double E_in );

  //! Finds an intersection with the bottom or top of a box
  //! Returns incident energy corresponding to E_out
  //! \param E_out the level to hit
  //! \param pair_0 ( Ein, Eout ) at lower incident energy
  //! \param pair_1 ( Ein, Eout ) at higher incident energy
  double find_hit( double E_out, const Ddvec::dd_entry &pair_0,
		   const Ddvec::dd_entry &pair_1 );

  //! Finds the intersections with the bottom of a box
  //! \param E_out bottom energy of the outgoing energy bin
  //! \param Ein_hits incident energies which give outgoing energy E_out
  void find_bottom_hits( double E_out, std::vector< Box::Ein_Eta_Hit > *Ein_hits );

  //! Finds the intersections with the top of a box
  //! \param E_out top energy of the outgoing energy bin
  //! \param Ein_hits incident energies which give outgoing energy E_out
  void find_top_hits( double E_out, std::vector< Box::Ein_Eta_Hit > *Ein_hits );

public:
  //! data for determination of how this outgoing energy curve hits the quadrature box
  Ddvec::dd_entry flip;  // Ein for minimal Eout, minimal Eout for this mucm 
  Ddvec::dd_entry left_Ein_Eout;  // ( Ein, Eout ) for this mucm and for lower Ein
  Ddvec::dd_entry right_Ein_Eout;  // ( Ein, Eout ) for this mucm and for higher Ein
  bool use_relativistic;

  Maps::Newton_map_param Newton_map;
  Rel::relativistic_param relParam;

  angle_hit_list( ) : use_relativistic( false ) {}
  ~angle_hit_list( ) {}
};

//! Class for angular distributions
// ---------------- class angle_dist ------------------
class angle_dist : public std::list< Ddvec::dd_vector >
{
private:
  Maps::two_body_map map;
  Rel::relativistic_map relMap;

  Adist::min_Eout_info_list min_Eout_list;  // minimal outgoing energies

  //! The smallest incident energy for cross section, multiplicity,
  //! model weight, flux weight, and energy groups
  double E_first;

  int first_Ein;  // index of the left-hand end of the first significant energy bin
  int last_Ein;  // index of the right-hand end of the last significant energy bin

  //! Accounts for the ENDL convention of isotropic emission at low energies
  void ENDL_kludge( );

  //! Returns the common negative cosines for the relativistic treatment of endothermic reactions
  std::list< double > negative_mu_list( );

  //! Initializes min_Eout_list, used in the treatment of endothermic reactions
  void init_min_Eout_list( );

  //! Computes min_Eout_list
  //! \param top_E_in: the highest incident energy
  void make_min_Eout_list( double top_E_in );

  //! Uses the mass difference to set the threshold
  void set_threshold( );

  //! Initializes the quadrature parameters
  //! \param mu_quad_rule rule of integration over cosine
  //! \param Ein_param the quadrature parameters for integration over incident energy
  void initialize_param( Qmeth::Quadrature_Rule mu_quad_rule, Adist::two_body_Ein_param
    *Ein_param );

  //! Adds to the transfer matrix for all E_out bins for a pair of incident energies.
  //! \param transfer the computed transfer matrix
  //! \param Ein_param the quadrature parameters for integration over incident energy
  void mucm_ladder( Trf::T_matrix& transfer, Adist::two_body_Ein_param *Ein_param );

  //! Sets up the map from center-of-mass to laboratory coordinates
  void setup_map( );

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
  //! \param Ein_param the quadrature parameters for integration over incident energy
  void setup_data( Adist::two_body_Ein_param *Ein_param );

  //! Sets the range of incident energies for this intergration
  //! \param Ein_bin the number of this incident energy bin
  //! \param Ein_param the quadrature parameters for integration over incident energy
  void set_Ein_range( int Ein_bin, Adist::two_body_Ein_param &Ein_param );

  //! Does the integration for one E-E' box between 2 mucm = const hyperbolas
  //! \param transfer the transfer matrix
  //! \param Eout_count count of the current outgoing energy bin
  //! \param Ein_param the quadrature parameters for integration over incident energy
  void one_Ebox( Trf::T_matrix& transfer, int Eout_count, Adist::two_body_Ein_param *Ein_param );

  //! Go to the next pair of incident energies.  Returns "true" when finished.
  //! \param E_in the next incident energy
  //! \param Ein_param the quadrature parameters for integration over incident energy
  bool next_ladder( double E_in, Adist::two_body_Ein_param *Ein_param );

  //! Adds to an element of transfer the integral between the intersections of 2 mucm = const hyperbolas with the E-E' box
  //! \param transfer the transfer matrix
  //! \param Eout_count count of the current outgoing energy bin
  //! \param Ein_param the quadrature parameters for integration over incident energy
  void update_T( Trf::T_matrix &transfer, int Eout_count, Adist::two_body_Ein_param *Ein_param );

public:

  double threshold;  // the computed threshold
  double threshold_out;  // outgoing energy at the computed threshold
  Terp::two_d_interp Ein_interp;  // interpolation between incident energies
  Terp::Interp_qualifier Ein_qualifier;  // for interpolation between incident energies
  Terp::Interp_Type mu_interp;

  bool use_relativistic;
  Rel::relativistic_masses relMass;

  //! Constructor
  inline angle_dist( ) : threshold( -1.0 ),
    	 mu_interp( Terp::LINLIN ), use_relativistic( false ) {}

  //! Destructor
  inline ~angle_dist( ) {}

  Maps::particleInfo particle_info;

  double Q;  // the reaction Q value

  //! Reads the data from Python
  //! \param input_file input file
  //! \param num_Ein number of incident energies for this reaction
  void read_data( Dpar::data_parser &input_file, int num_Ein );

  //! Calculates the transfer matrix for this particle
  //! \param sigma the cross section data
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param transfer the computed transfer matrix
  void get_T( const Ddvec::dd_vector& sigma, Ddvec::dd_vector& weight,
	      Trf::T_matrix& transfer );

  //! Checks whether an angular probability density is isotropic
  bool isotropic( );

  // for debugging
  void print( );

};

//! Class for parameters for the 2-d quadrature
// ---------------- class two_body_Ein_param ------------------
class two_body_Ein_param : public Pbase::param_base
{
private:
  //! Interpolates (mu_cm, probability) data to the lower common mu_cm value
  //! Returns true if the interpolation is OK
  //! \param lower_mu the lower end of the mu_cm range
  bool common_low_mucm( double lower_mu );

  //! Interpolates (mu_cm, probability) data to the higher common mu_cm value
  //! Returns true if the interpolation is OK
  //! \param higher_mu the upper end of the mu_cm range
  bool common_high_mucm( double higher_mu );

  //! Returns the Adist::min_Eout_info for given center-of-mass outgoing cosine
  //! \param mu_cm the center-of-mass outgoing cosine
  Adist::min_Eout_info_list::const_iterator get_min_Eout_info( double mu_cm ) const;

  //! Initializes upper_hits for the value of mucm
  //! \param mucm the center-of-mass direction cosine
  void set_upper_hits( double mucm );

public:
  bool use_relativistic;

  Maps::Newton_map_param Newton_map;
  Rel::relativistic_map relMap;
  
  double threshold;  // the computed threshold
  double threshold_out;  // outgoing energy at the computed threshold

  double left_data_Ein;                       // incident energy for left data
  double right_data_Ein;                      // incident energy for right data
  // pointers to the current pair of angular distirbutions
  Adist::angle_dist::const_iterator this_mucm_dist;
  Adist::angle_dist::const_iterator next_mucm_dist;
 
  // (mu_cm, probability) data interpolated to common mu_cm values
  Ddvec::dd_entry left_data;        // (lower mucm, probability) at low Ein
  Ddvec::dd_entry next_left_data;   // (higher mucm, probability) at low Ein
  Ddvec::dd_entry right_data;       // (lower mucm, probability) at high Ein
  Ddvec::dd_entry next_right_data;  // (higher mucm, probability) at high Ein
  long int quad_count;  // number of 2-d quadratures
  long int Ein_F_count;  // number of calls to angle_dist_F::E_quad_F
  long int mu_F_count;  // number of calls to angle_dist_F::mu_cm_quad_F
  Qmeth::Quadrature_Rule mu_quad_rule;  // rule of integration over cosine

  // where the lower and upper mucm values hit the quadrature box
  Adist::angle_hit_list lower_hits;
  Adist::angle_hit_list upper_hits;

  // pointers to the (mu_cm, probability) data
  Ddvec::dd_vector::const_iterator left_ptr;
  Ddvec::dd_vector::const_iterator next_left_ptr;
  Ddvec::dd_vector::const_iterator right_ptr;
  Ddvec::dd_vector::const_iterator next_right_ptr;

  // the range of outgoing energies for this data
  double upper_data_max_Eout;
  double lower_data_min_Eout;

  Adist::min_Eout_info_list *min_Eout_list;  // minimal outgoing energies

  inline two_body_Ein_param(): use_relativistic( false ), quad_count( 0 ),
      Ein_F_count( 0 ), mu_F_count( 0 ) {}
  inline ~two_body_Ein_param() {}

  //! Copies the relativistic mapping from center-of-mass to lab frame
  //! \param map_ the mapping data to copy
  void set_rel_map( Rel::relativistic_masses *relMass );

  //! Copies the Newtonian mapping from center-of-mass to lab frame
  //! \param map_ the mapping data to copy
  void set_Newton_map( Maps::two_body_map *map_ );

  //! Initializes the pointers to the angular probabilities for this E_in
  void reset_start( );

  //! Sets up the next interval of mu_cm values
  //! returns true if there is no more data
  bool next_mucm( );

  //! Initializes lower_hits and upper_hits for the values of mucm
  void reset_hits( );
};

} // end of namespace Adist

// **************** functions to integrate **********
namespace angle_dist_F
{
  // ---------------- mu_cm_quad_F ------------------
  //! Function for the 1-d quadrature
  //! Returns true if the interpolation is OK
  //! \param mucm center-of-mass direction cosine for outgoing particle
  //! \param mucm_quad_param parameters for this function
  //! \param value computed contribution to the transfer matrix
  bool mu_cm_quad_F( double mucm, Qparam::QuadParamBase *mucm_quad_param,
		     Coef::coef_vector *value );

  // ---------------- E_quad_F ------------------
  //! Function for the 2-d quadrature
  //! Returns true if the interpolation is OK
  //! \param E_in inergy of the incident particle
  //! \param mucm_quad_param parameters for this function
  //! \param value computed contribution to the transfer matrix
  bool E_quad_F( double E_in, Qparam::QuadParamBase *e_quad_param,
		 Coef::coef_vector *value );
}

#endif
