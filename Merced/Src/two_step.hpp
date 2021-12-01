/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2018-10-23 $
 * $Author: hedstrom $
 * $Id: two_step.hpp 1 2018-10-23Z hedstrom $
 *
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */
// define the classes used for 2-step 2-body reactions

#ifndef TWO_STEP_DEF
#define TWO_STEP_DEF

#include "coef_vector.hpp"
#include "param_base.hpp"
#include "transfer.hpp"
#include "quad_methods.hpp"
#include "two_step_hit.hpp"
#include "relativistic.hpp"

namespace Tstep
{

class two_step_param;  // forward declaration

//! Class for angular distributions as Legendre expansions
// ---------------- class two_step ------------------
class two_step : public std::list< Lgdata::Legendre_coefs >
{
private:
  double threshold;
  double threshold_out;  // outgoing energy at the threshold
  Ddvec::dd_entry flip;  // Ein for minimal Eout

  Maps::two_step_map Newton_map;
  Rel::two_step_relativistic_masses twoStepMasses;
  
  //! The smallest incident energy for cross section, multiplicity,
  //! model weight, flux weight, and energy groups
  double E_first;

  int first_Ein;  // index of the left-hand end of the first significant energy bin
  int last_Ein;  // index of the right-hand end of the last significant energy bin

  //! Uses the mass difference to set the threshold
  void set_threshold( );

  //! Initializes the quadrature parameters
  //! \param mucm1_quad_rule rule for integration over step 1 direction cosine
  //! \param mucm2_quad_rule rule for integration over step 2 direction cosine
  //! \param w_quad_rule rule for integration over step 2 w parameter
  //! \param Ein_param parameters for quadrature over incident energy
  void initialize_param( Qmeth::Quadrature_Rule mucm1_quad_rule,
     Qmeth::Quadrature_Rule mucm2_quad_rule, Qmeth::Quadrature_Rule w_quad_rule, 
     Tstep::two_step_param *Ein_param );

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
  void setup_map( );

  //! Sets up the Ein_param->map from center-of-mass to laboratory coordinates
  //! \param Ein_param parameters for quadrature over incident energy
  void setup_param_map( Tstep::two_step_param *Ein_param );

  //! Sets the range of incident energies for this intergration
  //! \param Ein_param parameters for quadrature over incident energy
  void set_Ein_range( Tstep::two_step_param *Ein_param );

  //! Initializes the quadrature parameters
  //! \param Ein_param parameters for quadrature over incident energy
  void setup_data( Tstep::two_step_param *Ein_param );

  //! Adds to the transfer matrix for all E_out bins for a pair of incident energies.
  //! \param transfer the transfer matrix
  //! \param Ein_param parameters for quadrature over incident energy
  void Eout_ladder( Trf::T_matrix& transfer, Tstep::two_step_param *Ein_param );

  //! Does the integration for one E-E' box between 2 eta = const hyperbolas
  //! \param transfer the transfer matrix
  //! \param Eout_count count of the current outgoing energy bin
  //! \param Ein_param parameters for quadrature over incident energy
  void one_Ebox( Trf::T_matrix& transfer, int Eout_count,
		 Tstep::two_step_param *Ein_param );

  //! Go to the next pair of incident energies.  Returns "true" when finished.
  //! \param E_in the next incident energy
  //! \param Ein_param parameters for quadrature over incident energy
  bool next_ladder( double E_in, Tstep::two_step_param *Ein_param );

  //! Adds to an element of transfer the integral between the intersections of 2 eta = const hyperbolas with the E-E' box
  //! \param transfer the transfer matrix
  //! \param Eout_count identify the row of the transfer matrix
  //! \param Ein_param parameters for quadrature over incident energy
  void update_T( Trf::T_matrix &transfer, int Eout_count,
		 Tstep::two_step_param *Ein_param );

public:
  // the particle masses and Qs
  Maps::particleInfo step1_particles;  // particles for the first step
  Maps::particleInfo step2_particles;  // particles for the second step

  double first_Q;  // the first reaction Q value
  double second_Q;  // the second reaction Q value

  bool use_relativistic;
  
  Terp::Interp_Type Ein_interp;  // interpolation between incident energies

  two_step( ): first_Q( -1.0 ), second_Q( -1.0 ), use_relativistic( false ),
	       Ein_interp( Terp::NOTSET ) {}
  ~two_step( ) {}

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

//! Class for parameters for the 4-d quadrature based on a set of Legendre coefficients
// ---------------- class two_step_param ------------------
class two_step_param : public Pbase::param_base
{
public:
  Lgdata::Legendre_coefs coefs;
  Maps::two_step_map_param twoStepMap;  // Newtonian 
  Rel::two_step_relativistic_map relTwoStepMap;

  int num_negative;
  
  bool use_relativistic;

  Terp::Interp_Type Ein_interp;
  
  int order;   // the Legendre order of the transfer matrix
  Coef::Conserve conserve;  // what we calculate in the transfer matrix

  double threshold;
  double threshold_out;  // outgoing energy at the threshold

  Tstep::two_step::const_iterator left_data;   // pointers to (eta, probability) data
  Tstep::two_step::const_iterator right_data;

  // where the fore and aft emissions hit the quadrature box
  std::vector< double >::const_iterator Eout_ptr;  // bottom of the energy bin
  Twohit::two_step_hit_list_array fore_aft_hits;
 
  long int quad_count;  // number of 4-d quadratures
  long int Ein_F_count;  // number of calls to twoStep_F::Ein_F
  long int mucm1_F_count;  // number of calls to twoStep_F::mucm_1_F
  long int mucm2_F_count;  // number of calls to twoStep_F::mucm_2_F
  long int w_F_count;  // number of calls to twoStep_F::step2_w_F

  // quadrature rule for center-of-mass cosine
  Qmeth::Quadrature_Rule mucm1_quad_rule;
  Qmeth::Quadrature_Rule mucm2_quad_rule;
  Qmeth::Quadrature_Rule w_quad_rule;

  // for the range of integration over mucm_1
  bool use_mucm1_max;  // for mucm1_max = 1
  bool use_mucm1_min;  // for mucm1_min = -1
  
  two_step_param( );
  inline ~two_step_param( ) {}

  //! Interpolates between two incident energies
  //! Returns false if there are problems
  //! \param Ein, the intermediate incident energy
  //! \param prev_coef, Legendre coefficients at the lower incident energy
  //! \param next_coef, Legendre coefficients at the higher incident energy
  bool interpolate( double Ein, Tstep::two_step::const_iterator prev_coef,
		    Tstep::two_step::const_iterator next_coef );
};

//! Class for 3-d quadrature over mucm_1, mucm_2, and w
// ------------------- class mucm_1_param ---------------
class mucm_1_param :  public Qparam::QuadParamBase
{
private:

  //! Gets the range of mucm_2 cosines for this sector
  //! \param mucm1, the direction cosine for the first reaction
  //! \param mucm2_aft, the minimal mucm2 (computed)
  //! \param mucm2_fore, the maximal mucm2 (computed)
  void get_mucm2_range( double mucm1, 
			double *mucm2_aft, double *mucm2_fore );
 
public:
  double E_in;  // incident energy

  std::vector< double >::const_iterator Eout_ptr;  // bottom of the energy bin
  double Eout_min;  // bottom of the outgoing energy bin
  double Eout_max;  // top of the outgoing energy bin

  double Etrans2;   // energy component due to motion of final center of mass

  long int mucm2_F_count;  // number of calls to twoStep_F::mucm_2_F
  long int w_F_count;  // number of calls to twoStep_F::step2_w_F

  Qmeth::Quadrature_Rule mucm2_quad_rule;
  Qmeth::Quadrature_Rule w_quad_rule;

  int order;   // the Legendre order of the transfer matrix
  Coef::Conserve conserve;  // what we calculate in the transfer matrix

  Lgdata::Legendre_coefs coefs;  // interpolated Legendre coefficients
  Terp::Interp_Type Ein_interp;  // interpolation between incident energies

  int num_negative;
  
  Maps::two_step_map *map;
  Rel::two_step_relativistic_map relTwoStepMap;

  bool use_relativistic;

  inline mucm_1_param( ): mucm2_F_count( 0 ), w_F_count( 0 ),
			  num_negative( 0 ), use_relativistic( false ) {}
  inline ~mucm_1_param( ) {}

  //! Sets up the data for step 2 cosine mucm2
  //! \param Ein, incident energy
  //! \param Ein_param parameters for the 4-d quadrature
  void setup( double Ein, Tstep::two_step_param *Ein_param );

  //! Does the integration over mucm_2 and w
  //! Returns true if there are no problems
  //! \param mucm1, direction cosine for the first reaction
  //! \param value, the computed integral over mucm_2 and w
  bool integrate( double mucm1, Coef::coef_vector *value );
  
  //! Finds how the geometry depends on mucm_1
  //! \param mucm1, the direction cosine for step 1
  void get_geometry( double mucm1, std::vector< double >::const_iterator Eout_ptr );
 
  //! Interpolates between two incident energies
  //! Returns false if there are problems
  //! \param Ein, the intermediate incident energy
  //! \param prev_coef, Legendre coefficients at the lower incident energy
  //! \param next_coef, Legendre coefficients at the higher incident energy
  bool interpolate( double Ein, Tstep::two_step::const_iterator prev_coef,
                    Tstep::two_step::const_iterator next_coef );
};

//! Class for 2-d quadrature over mucm_2 and w for step 2 of the reaction
// ------------------- class mucm_2_param ---------------
class mucm_2_param :  public Qparam::QuadParamBase
{
public:
  long int w_F_count;  // number of calls to twoStep_F::step2_w_F

  Qmeth::Quadrature_Rule w_quad_rule;

  double E_in;  // incident energy
  double mucm_1; // direction cosine for step 1
  double mucm_2; // direction cosine for step 2
  double Etrans2;   // energy component due to motion of final center of mass
  double Eout_2; // outgoing energy for particle of step 2

  Maps::two_step_map *map;
  Rel::two_step_relativistic_map relTwoStepMap;

  bool use_relativistic;

  inline mucm_2_param( ): w_F_count( 0 ) {}
  inline ~mucm_2_param( ) {}

  //! Sets up the data for step 2 cosine mucm2
  //! \param mucm1, direction cosine for step 1
  //! \param mucm1_param parameters for quadrature over cosine for step 1
  void setup( double mucm1, const Tstep::mucm_1_param *mucm1_param );

};

//! Class for 1-d quadrature over w for step 2 of the reaction
// ------------------- class step2_w_param ---------------
class step2_w_param :  public Qparam::QuadParamBase
{
public:
  double E_in;  // incident energy
  double mucm_1; // direction cosine for step 1
  double mucm_2; // direction cosine for step 2
  double Etrans2;   // energy component due to motion of final center of mass
  double Eout_2; // outgoing energy for particle of step 2

  Maps::two_step_map *map;
  Rel::two_step_relativistic_map relTwoStepMap;

  bool use_relativistic;

  inline step2_w_param( ) {}
  inline ~step2_w_param( ) {}

  //! Sets up the data for step 2 cosine mucm2
  //! \param mucm2, direction cosine for step 2
  //! \param Eout2, outgoing energy for step 2
  //! \param mucm2_param parameters for quadrature over cosine for step 2
  void setup( double mucm2, double Eout2,
	      const Tstep::mucm_2_param *mucm2_param );

};

} // end of namespace Tstep

// **************** functions to integrate **********
namespace twoStep_F
{
  // ---------------- step2_w_F ------------------
  //! Function for the 1-d quadrature over w
  //! \param w, parameter for direction of outgoing particle in step 2
  //! \param step2w_param parameters for this function
  //! \param value computed contribution to the transfer matrix
  bool step2_w_F( double w, Qparam::QuadParamBase *step2w_param,
		  Coef::coef_vector *value );

  // ---------------- mucm_2_F ------------------
  //! Function for the 2-d quadrature over mucm_2, w
  //! \param mucm_2 direction cosine for outgoing particle in step 2
  //! \param mucm2_param parameters for this function
  //! \param value computed contribution to the transfer matrix
  bool mucm_2_F( double mucm_2, Qparam::QuadParamBase *mucm2_param,
		 Coef::coef_vector *value );

  // ---------------- mucm_1_F ------------------
  //! Function for the 3-d quadrature over mucm_1, mucm_2, w
  //! \param mucm_1 center-of-mass direction cosine for outgoing particle in step 1
  //! \param mucm1_param parameters for this function
  //! \param value computed contribution to the transfer matrix
  bool mucm_1_F( double mucm_1, Qparam::QuadParamBase *mucm1_param,
		 Coef::coef_vector *value );

  // ---------------- Ein_F ------------------
  //! Function for the 4-d quadrature over (E_in, mucm_1, mucm_2, w )
  //! \param E_in inergy of the incident particle
  //! \param e_quad_param parameters for this function
  //! \param value computed contribution to the transfer matrix
  bool Ein_F( double E_in, Qparam::QuadParamBase *e_quad_param,
	      Coef::coef_vector *value );
}

#endif
