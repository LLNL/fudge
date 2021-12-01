/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2018-10-23 $
 * $Author: hedstrom $
 * $Id: two_step.cpp 1 2018-10-23Z hedstrom $
 *
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */
// Implement the classes used for 2-step 2-body reactions

#ifdef _OPENMP
 #include <omp.h>
#endif

#include <iostream>

#include "two_step.hpp"
#include "adapt_quad.hpp"
#include "global_params.hpp"
#include "math_util.hpp"
#include "messaging.hpp"

// **************** Tstep::step2_w_param ****************
// -------------- Tstep::step2_w_param::setup ----------------------
// Sets up the data for step 2 cosine mucm2
void Tstep::step2_w_param::setup( double mucm2, double Eout2,
			   const Tstep::mucm_2_param *mucm2_param )
{
  // copy the data
  map = mucm2_param->map;
  E_in = mucm2_param->E_in;
  mucm_1 = mucm2_param->mucm_1;
  mucm_2 = mucm2;
  Etrans2 = mucm2_param->Etrans2;
  Eout_2 = Eout2;

  use_relativistic = mucm2_param->use_relativistic;
  relTwoStepMap.twoStepMasses = mucm2_param->relTwoStepMap.twoStepMasses;
}

// ************* Tstep::mucm_2_param *************
// ---------------- Tstep::mucm_2_param::setup ------------------
void Tstep::mucm_2_param::setup( double mucm1, const Tstep::mucm_1_param *mucm1_param )
{
  w_quad_rule = mucm1_param->w_quad_rule;

  mucm_1 = mucm1;
  E_in = mucm1_param->E_in;
  Etrans2 = mucm1_param->Etrans2;
  
  map = mucm1_param->map;
  use_relativistic = mucm1_param->use_relativistic;
  relTwoStepMap.twoStepMasses = mucm1_param->relTwoStepMap.twoStepMasses;
}

// ************* Tstep::mucm_1_param *************
// ---------------- Tstep::mucm_1_param::setup ------------------
void Tstep::mucm_1_param::setup( double Ein, Tstep::two_step_param *Ein_param )
{

  Ein_interp = Ein_param->Ein_interp;
  
  mucm2_quad_rule = Ein_param->mucm2_quad_rule;
  w_quad_rule = Ein_param->w_quad_rule;

  E_in = Ein;
  Eout_ptr = Ein_param->Eout_ptr;
  Eout_min = Ein_param->Eout_min;
  Eout_max = Ein_param->Eout_max;
  map = Ein_param->twoStepMap.map;

  use_relativistic = Ein_param->use_relativistic;
  relTwoStepMap.twoStepMasses = Ein_param->relTwoStepMap.twoStepMasses;
}

// ---------------- Tstep::mucm_1_param::integrate ------------------
// Does the integration over mucm_2 and w
bool Tstep::mucm_1_param::integrate( double mucm_1, Coef::coef_vector *value )
{
  static double tol = Global.Value( "quad_tol" );

  // parameters for integration over mucm_2 and w
  Tstep::mucm_2_param mucm2_param;
  mucm2_param.setup( mucm_1, this );
  Qparam::QuadParamBase *params = static_cast< Qparam::QuadParamBase* >( &mucm2_param );
  
  // the range of integration
  double mucm2_aft;
  double mucm2_fore;

  get_mucm2_range( mucm_1, &mucm2_aft, &mucm2_fore );

  bool is_OK = quad_F::integrate( twoStep_F::mucm_2_F, mucm2_quad_rule,
	      mucm2_aft, mucm2_fore, params, tol, value );
    
  // increment the function counts
  mucm2_F_count += mucm2_param.func_count;
  w_F_count += mucm2_param.w_F_count;

  return is_OK;
}
// ---------------- Tstep::mucm_1_param::get_mucm2_range ------------------
// Gets the range of mucm_2 cosines for this sector
void Tstep::mucm_1_param::get_mucm2_range( double mucm1, double *mucm2_aft,
				    double *mucm2_fore )
{
  // parameters for get_mucm2
  Rel::relativistic_2_step_param step_2_param;
  step_2_param.relTwoStepMap.twoStepMasses = relTwoStepMap.twoStepMasses;

  step_2_param.Tin_lab = E_in;
  step_2_param.mu_cm_1 = mucm1;

  double aft_Eout;
  double fore_Eout;

  if( use_relativistic )
  {
    aft_Eout = step_2_param.relTwoStepMap.two_step_get_E_lab( E_in, mucm1, -1.0 );
    fore_Eout = step_2_param.relTwoStepMap.two_step_get_E_lab( E_in, mucm1, 1.0 );
  }
  else
  {
    aft_Eout = map->two_step_get_E_lab( E_in, mucm1, -1.0 );
    fore_Eout = map->two_step_get_E_lab( E_in, mucm1, 1.0 );
  }
  
  if( fore_Eout <= Eout_min )
  {
    Msg::Warning( "Tstep::mucm_1_param::get_mucm2_range", "data below bin" );
    *mucm2_aft = -1.0;
    *mucm2_fore = -1.0;
  }
  else if( fore_Eout <= Eout_max )
  {
    *mucm2_fore = 1.0;
    if( aft_Eout < Eout_min )
    {
      if( use_relativistic )
      {
        *mucm2_aft = step_2_param.get_mucm2( Eout_min );
      }
      else
      {
        *mucm2_aft = map->get_mucm2( Etrans2, Eout_min );
      }
    }
    else
    {
      *mucm2_aft = -1.0;
    }
  }
  // from here on, fore__Eout > Eout_max
  else if( aft_Eout < Eout_min )
  {
    if( use_relativistic )
    {
      *mucm2_fore = step_2_param.get_mucm2( Eout_max );
      *mucm2_aft = step_2_param.get_mucm2( Eout_min );
    }
    else
    {
      *mucm2_fore = map->get_mucm2( Etrans2, Eout_max );
      *mucm2_aft = map->get_mucm2( Etrans2, Eout_min );
    }
  }
  else if( aft_Eout < Eout_max )
  {
    if( use_relativistic )
    {
      *mucm2_fore = step_2_param.get_mucm2( Eout_max );
    }
    else
    {
      *mucm2_fore = map->get_mucm2( Etrans2, Eout_max );
    }
    *mucm2_aft = -1.0;
  }
  else
  {
    Msg::Warning( "Tstep::mucm_1_param::get_mucm2_range", "data above bin" );
    *mucm2_fore = -1.0;
    *mucm2_aft = -1.0;
  }
}
// ---------------- Tstep::mucm_1_param::interpolate ------------------
// Interpolates between two incident energies
bool Tstep::mucm_1_param::interpolate( double Ein,
   Tstep::two_step::const_iterator prev_coefs,
   Tstep::two_step::const_iterator next_coefs )
{
  bool interp_OK;
  if( Ein_interp == Terp::LINLIN )
  {
    interp_OK = coefs.linlin_interp( Ein, *prev_coefs, *next_coefs );
  }
  else
  {
    interp_OK = coefs.linlog_interp( Ein, *prev_coefs, *next_coefs );
  }
  return interp_OK;
}

// ************* Tstep::two_step_param *************
// ---------------- Tstep::two_step_param::two_step_param ------------------
// constructor
Tstep::two_step_param::two_step_param( )
{
  num_negative = 0;
  quad_count = 0;
  Ein_F_count = 0;
  mucm1_F_count = 0;
  mucm2_F_count = 0;
  w_F_count = 0;
}
// ---------------- Tstep::two_step_param::interpolate ------------------
// Interpolates between two incident energies
bool Tstep::two_step_param::interpolate( double Ein,
   Tstep::two_step::const_iterator prev_coefs,
   Tstep::two_step::const_iterator next_coefs )
{
  bool interp_OK;
  if( Ein_interp == Terp::LINLIN )
  {
    interp_OK = coefs.linlin_interp( Ein, *prev_coefs, *next_coefs );
  }
  else
  {
    interp_OK = coefs.linlog_interp( Ein, *prev_coefs, *next_coefs );
  }
  return interp_OK;
}

// ************* Tstep::two_step *************
// ---------------- Tstep::two_step::setup_map ------------------
void Tstep::two_step::setup_map( )
// set up the map from center-of-mass to laboratory coordinates
{
  // function parameter
  void *params;

  if( use_relativistic )
  {
    twoStepMasses.set_masses( &step1_particles, first_Q, &step2_particles,
				  second_Q );

    twoStepMasses.step1Masses.get_threshold( );
    threshold = twoStepMasses.step1Masses.threshold;

    Rel::relativistic_2_step_param relTwoStepParam;
    relTwoStepParam.relTwoStepMap.twoStepMasses = &twoStepMasses;
    
    // we need to set the direction cosines
    relTwoStepParam.mu_cm_1 = -1.0;
    relTwoStepParam.mu_cm_2 = -1.0;
    params = static_cast< void * >( &relTwoStepParam );
    threshold_out = relativistic_F::two_step_T_out_lab( threshold, params );
    flip = relTwoStepParam.find_lowest_bottom( );
  }
  else
  {
    Newton_map.set_map( step1_particles, first_Q, step2_particles );
    Newton_map.get_threshold( );
    threshold = Newton_map.threshold;

    // parameters for the Newtonian_F functions
    Maps::two_step_map_param Newton_param;
    Newton_param.mucm_1 = -1.0;
    Newton_param.mucm_2 = -1.0;
    Newton_param.map = &Newton_map;
    params = static_cast< void * >( &Newton_param );
    threshold_out = Newtonian_F::two_step_T_out_lab( threshold, params );
    flip = Newton_param.find_lowest_bottom( );
  }
}
// ---------------- Tstep::two_step::setup_param_map ------------------
void Tstep::two_step::setup_param_map( Tstep::two_step_param *Ein_param )
// set up the Ein_param->map from center-of-mass to laboratory coordinates
{
  if( use_relativistic )
  {
    Ein_param->relTwoStepMap.twoStepMasses = &twoStepMasses;
    
    for( int list_count = 0; list_count < 4; ++list_count )
    {
      Ein_param->fore_aft_hits.hit_lists[ list_count ].relTwoStepParam.relTwoStepMap.twoStepMasses =
	&twoStepMasses;
    }
  }
  else
  {
    Ein_param->twoStepMap.map = &Newton_map;
    for( int list_count = 0; list_count < 4; ++list_count )
    {
      Ein_param->fore_aft_hits.hit_lists[ list_count ].twoStepMap.map = &Newton_map;
    }
  }
  
  // set the direction cosines
  Ein_param->fore_aft_hits.hit_lists[ 0 ].set_mucm_1( -1.0 );  // aft, aft
  Ein_param->fore_aft_hits.hit_lists[ 0 ].set_mucm_2( -1.0 );  // aft, aft

  Ein_param->fore_aft_hits.hit_lists[ 1 ].set_mucm_1( -1.0 );  // aft, fore
  Ein_param->fore_aft_hits.hit_lists[ 1 ].set_mucm_2( 1.0 );  // aft, fore

  Ein_param->fore_aft_hits.hit_lists[ 2 ].set_mucm_1( 1.0 );  // fore, aft
  Ein_param->fore_aft_hits.hit_lists[ 2 ].set_mucm_2( -1.0 );  // fore, aft

  Ein_param->fore_aft_hits.hit_lists[ 3 ].set_mucm_1( 1.0 );  // fore, fore
  Ein_param->fore_aft_hits.hit_lists[ 3 ].set_mucm_2( 1.0 );  // fore, fore  
}
// ----------- Tstep::two_step::set_threshold --------------
// Uses the mass difference to set the threshold
void Tstep::two_step::set_threshold( )
{
  if( use_relativistic )
  {
    threshold = twoStepMasses.step1Masses.threshold;
  }
  else
  {
    threshold = Newton_map.threshold;
  }

  // adjust the data if necessary
  Tstep::two_step::iterator data_ptr = begin( );
  double first_Ein = data_ptr->get_E_in( );
  if( first_Ein < threshold )
  {
    data_ptr->set_E_in( threshold );
  }
}
// ----------- Tstep::two_step::initialize_param --------------
// Initializes the quadrature parameters
void Tstep::two_step::initialize_param( Qmeth::Quadrature_Rule mucm1_quad_rule, 
     Qmeth::Quadrature_Rule mucm2_quad_rule, Qmeth::Quadrature_Rule w_quad_rule, 
     Tstep::two_step_param *Ein_param )
{
  Ein_param->use_relativistic = use_relativistic;
  
  Ein_param->fore_aft_hits.get_array( 4 );  // allocate space
  for( int list_count = 0; list_count < 4; ++list_count )
  {
    Ein_param->fore_aft_hits.hit_lists[ list_count ].use_relativistic = use_relativistic;
  }

  Ein_param->threshold = threshold;
  Ein_param->threshold_out = threshold_out;
  for( int list_count = 0; list_count < 2; ++list_count )
  {
    Ein_param->fore_aft_hits.hit_lists[ list_count ].flip.x = flip.x;  // for mucm_1 = -1
    Ein_param->fore_aft_hits.hit_lists[ list_count ].flip.y = flip.y;  // for mucm_1 = -1
  }
  for( int list_count = 2; list_count < 4; ++list_count )
  {
    Ein_param->fore_aft_hits.hit_lists[ list_count ].flip.x = threshold;  // for mucm_1 = 1
    Ein_param->fore_aft_hits.hit_lists[ list_count ].flip.y = threshold_out;  // for mucm_1 = 1
  }

  Ein_param->Ein_interp = Ein_interp;
  Ein_param->mucm1_quad_rule = mucm1_quad_rule;
  Ein_param->mucm2_quad_rule = mucm2_quad_rule;
  Ein_param->w_quad_rule = w_quad_rule;
}
// ----------- Tstep::two_step::read_data --------------
void Tstep::two_step::read_data( Dpar::data_parser &input_file, int num_Ein )
{
  Ein_interp = interp_flag_F::read_1d_interpolation( input_file );

  Lgdata::Legendre_coefs new_angle_dist;  // angular distribution for one E_in
  Tstep::two_step::iterator new_angle_ptr;

  // read the data
  for( int Ein_count = 0; Ein_count < num_Ein; ++Ein_count )
  {
    // insert a new angular distribution
    new_angle_ptr = insert( end( ), Lgdata::Legendre_coefs( ) );
    // get the incident energy and the data pairs
    new_angle_ptr->set_E_in( input_file.get_next_double( ) );
    int Order = input_file.get_next_int( ) - 1;
    new_angle_ptr->initialize( Order );
    for( int coef_count = 0; coef_count <= Order; ++coef_count )
    {
      (*new_angle_ptr)[ coef_count ] = input_file.get_next_double( );
    }
  }
  new_angle_ptr->truncate_zeros( );
}
// ----------- Tstep::two_step::get_T --------------
// Calculates the transfer matrix for this particle.
// sigma is the cross section.
void Tstep::two_step::get_T( const Ddvec::dd_vector& sigma,
		      const Ddvec::dd_vector& weight,
  Trf::T_matrix& transfer )
{
  // set the flag for relativistic mechanics
  if( Global.Flag( "kinetics" ) == "newtonian" )
  {
    use_relativistic = false;
  }
  else
  {
    use_relativistic = true;
  }

  if( ( Ein_interp != Terp::LINLIN ) && ( Ein_interp != Terp::LINLOG ) )
  {
    Msg::FatalError( "Tstep::two_step::get_T",
      "Incident energy interpolation not implemented" );
  }

  // the multiplicity is one
  Ddvec::dd_vector multiple;
  multiple.make_flat( sigma, 1.0 );
  bool done = get_Ein_range( sigma, multiple, weight, transfer.e_flux,
    transfer.in_groups );
  if( done )
  {
    transfer.zero_transfer( );
  }

  // for center-of-mass data
  setup_map( );
  // use the threshold computed from the kinetics
  set_threshold( ); 
  transfer.threshold = threshold;
  
  // the computed threshold may be above 20 MeV
  Ddvec::dd_vector::const_iterator last_sigma = sigma.end( );
  --last_sigma;
  if( last_sigma->x < threshold )
  {
    Msg::Info( "Tstep::two_step::get_T",
	       "computed threshold outside the data range" );
    transfer.zero_transfer( );
  }
  
  int num_negative = 0; // number of negative sums of Legendre series
  long int quad_count = 0;  // number of 4-d quadratures
  long int Ein_F_count= 0;  // number of calls to twoStep_F::Ein_F
  long int mucm1_F_count = 0;  // number of calls to twoStep_F::mucm_1_F
  long int mucm2_F_count = 0;  // number of calls to twoStep_F::mucm_2_F
  long int w_F_count= 0;  // number of calls to twoStep_F::w_F

  // now do the integrals incident bin by incident bin
#pragma omp parallel for schedule( dynamic, 1 ) default( none )            \
  shared( sigma, multiple, weight, transfer )				\
  reduction( +: quad_count ) reduction( +: Ein_F_count ) \
  reduction( +: mucm1_F_count ) reduction( +: num_negative ) \
  reduction( +: mucm2_F_count ) reduction( +: w_F_count )
  for( int Ein_bin = first_Ein; Ein_bin < last_Ein; ++Ein_bin )
  {
    Tstep::two_step_param Ein_param;
    initialize_param( transfer.mu_quad_rule, transfer.mucm2_quad_rule,
		      transfer.w_quad_rule, &Ein_param );
    setup_param_map( &Ein_param );

    // set up the data range for this bin
    Ein_param.setup_bin( Ein_bin, sigma, multiple, weight, transfer.e_flux,
                         transfer.in_groups );
    setup_data( &Ein_param );
    for( ; ; )
    {
      set_Ein_range( &Ein_param );   // get the incident energy interval
      Eout_ladder( transfer, &Ein_param );
      bool done = next_ladder( Ein_param.data_E_1, &Ein_param );   // go to the next interval
      if( done )
      {
        break;
      
      }
    }
    num_negative += Ein_param.num_negative;
    quad_count += Ein_param.quad_count;
    Ein_F_count += Ein_param.Ein_F_count;
    mucm1_F_count += Ein_param.mucm1_F_count;
    mucm2_F_count += Ein_param.mucm2_F_count;
    w_F_count += Ein_param.w_F_count;

  }
  if( num_negative > 0 )
  {
    Msg::Info( "Tstep::two_step::get_T", Msg::pastenum( "got ", num_negative ) +
	  " negative Legendre sums." );
  }  // end of parallel loop

  // print the counts of function evaluations
  std::cout << "4d quadratures: " << quad_count << std::endl;
  std::cout << "twoStep_F::Ein_F calls: " << Ein_F_count << std::endl;
  std::cout << "twoStep_F::mucm_1_F calls: " << mucm1_F_count << std::endl;
  std::cout << "twoStep_F::mucm_2_F calls: " << mucm2_F_count << std::endl;
  std::cout << "twoStep_F:w_F calls: " << w_F_count << std::endl;

  std::cout << "average twoStep_F::Ein_F calls: " << 1.0*Ein_F_count/quad_count << std::endl;
  std::cout << "average twoStep_F::mucm_1_F calls: " << 1.0*mucm1_F_count/Ein_F_count << std::endl;
  std::cout << "average twoStep_F::mucm_2_F calls: " << 1.0*mucm2_F_count/mucm1_F_count << std::endl;
  std::cout << "average twoStep_F::w_1_F calls: " << 1.0*w_F_count/mucm2_F_count << std::endl;
}
// ----------- Tstep::two_step::get_Ein_range --------------
//  Gets the range of nontrivial incident energy bins; computes first_Ein and last_Ein
// returns true if the threshold is too high for the energy bins
bool Tstep::two_step::get_Ein_range( const Ddvec::dd_vector& sigma, const Ddvec::dd_vector& mult,
    const Ddvec::dd_vector& weight,
    const Lgdata::Flux_List& e_flux, const Egp::Energy_groups& Ein_groups )
{
  double E_last;

  Tstep::two_step_param initial_param;
  bool done = initial_param.get_Ein_range( sigma, mult, weight, e_flux,
                                         Ein_groups, &E_first, &E_last );
  if( done ) return true;

  // check the range of incident energies for the probability data
  Tstep::two_step::const_iterator Ein_data_ptr = begin( );
  double E_data = Ein_data_ptr->get_E_in( );
  if( E_data > E_first )
  {
    E_first = E_data;
  }
  first_Ein = Ein_groups.first_bin_ID( E_first );

  Ein_data_ptr = end( );
  --Ein_data_ptr;
  E_data = Ein_data_ptr->get_E_in( );
  if( E_data < E_last )
  {
    E_last = E_data;
  }
  last_Ein = Ein_groups.last_bin_ID( E_last );

  return false;
}
// ----------- Tstep::two_step::setup_data --------------
// Initializes the quadrature parameters
void Tstep::two_step::setup_data( Tstep::two_step_param *Ein_param )
{
  static double skip_tol = Global.Value( "tight_tol" );
  
  Ein_param->Ein_interp = Ein_interp;
    
  Ein_param->left_data = begin( );
  Ein_param->right_data = Ein_param->left_data;
  ++Ein_param->right_data;

  while( ( Ein_param->right_data->get_E_in( ) < E_first * ( 1.0 + skip_tol ) ) ||
	 ( Ein_param->right_data->get_E_in( ) < (*Ein_param->Ein_ptr) *
           ( 1.0 + skip_tol ) ) )
  {
    Ein_param->left_data = Ein_param->right_data;
    ++Ein_param->right_data;
  }

  double first_E = Ein_param->left_data->get_E_in( );
  if( first_E > Ein_param->data_E_0 )
  {
    Ein_param->data_E_0 = first_E;
    bool data_bad = Ein_param->update_pointers( first_E );
    if( data_bad )
    {
      Msg::FatalError( "Tstep::two_step::setup_data", "energies inconsistent" );
    }
  }
}
// ----------- Tstep::two_step::next_ladder --------------
bool Tstep::two_step::next_ladder( double E_in, Tstep::two_step_param *Ein_param )
{
  bool done = Ein_param->update_bin_pointers( E_in );
  static double etol = Global.Value( "tight_tol" );
  if( !done )
  {
    double E_tol = E_in * etol;
    if( E_in + E_tol >= Ein_param->right_data->get_E_in( ) )
    {
      while( E_in + E_tol >= Ein_param->right_data->get_E_in( ) )
      {
        // get the next angular data
        Ein_param->left_data = Ein_param->right_data;
        ++Ein_param->right_data;
        if( Ein_param->right_data == end ( ) )
        {
          return true;
        }
      }
    }
  }
  return done;
}

// ----------- Tstep::two_step::set_Ein_range --------------
// Sets the range of incident energies for this intergration
void Tstep::two_step::set_Ein_range( Tstep::two_step_param *Ein_param )
{
  Ein_param->set_Ein_range( );
  double this_E = Ein_param->left_data->get_E_in( );
  if( this_E > Ein_param->data_E_0 ) Ein_param->data_E_0 = this_E;
  this_E = Ein_param->right_data->get_E_in( );
  if( this_E < Ein_param->data_E_1 ) Ein_param->data_E_1 = this_E;

  // the data may be below the threshold
  if( Ein_param->data_E_1 < threshold )
  {
    Ein_param->data_E_0 = threshold;
    Ein_param->data_E_1 = threshold;
    // Msg::Warning( "Tstep::two_step::set_Ein_range",
    //    "ignoring 2 data below threshold" );
  }
  else if( Ein_param->data_E_0 < threshold )
  {
    Ein_param->data_E_0 = threshold;
    // Msg::Warning( "Tstep::two_step::set_Ein_range",
    //    "ignoring data below threshold" );
  }

  if( Ein_param->data_E_1 < Ein_param->data_E_0 )
  {
    Msg::FatalError( "Tstep::two_step::set_Ein_range",
		     "check the incident energies" );
  }
  Ein_param->set_sigma_range( );

  // set the energies for Ein_param->fore_aft_hits
  double left_data_Ein = Ein_param->left_data->get_E_in( );
  for( int list_count = 0; list_count < 4; ++list_count )
  {
    Ein_param->fore_aft_hits.hit_lists[ list_count ].left_Ein_Eout.x = left_data_Ein;
  }

  double right_data_Ein = Ein_param->right_data->get_E_in( );
  for( int list_count = 0; list_count < 4; ++list_count )
  {
    Ein_param->fore_aft_hits.hit_lists[ list_count ].right_Ein_Eout.x = right_data_Ein;
  }

  if( use_relativistic )
  {
    Rel::two_step_relativistic_map relTwoStepMap;
    relTwoStepMap.twoStepMasses = &twoStepMasses;

    for( int list_count = 0; list_count < 4; ++list_count )
    {
      Ein_param->fore_aft_hits.hit_lists[ list_count ].left_Ein_Eout.y =
	relTwoStepMap.two_step_get_E_lab( left_data_Ein,
		    Ein_param->fore_aft_hits.hit_lists[ list_count ].relTwoStepParam.mu_cm_1,
		    Ein_param->fore_aft_hits.hit_lists[ list_count ].relTwoStepParam.mu_cm_2 );
      
      Ein_param->fore_aft_hits.hit_lists[ list_count ].right_Ein_Eout.y =
	relTwoStepMap.two_step_get_E_lab( right_data_Ein,
		    Ein_param->fore_aft_hits.hit_lists[ list_count ].relTwoStepParam.mu_cm_1,
		    Ein_param->fore_aft_hits.hit_lists[ list_count ].relTwoStepParam.mu_cm_2 );
    }
  }
  else
  {
    for( int list_count = 0; list_count < 4; ++list_count )
    {
      Ein_param->fore_aft_hits.hit_lists[ list_count ].left_Ein_Eout.y =
	Newton_map.two_step_get_E_lab( left_data_Ein,
		    Ein_param->fore_aft_hits.hit_lists[ list_count ].twoStepMap.mucm_1,
		    Ein_param->fore_aft_hits.hit_lists[ list_count ].twoStepMap.mucm_2 );
      
      Ein_param->fore_aft_hits.hit_lists[ list_count ].right_Ein_Eout.y =
	Newton_map.two_step_get_E_lab( right_data_Ein,
		    Ein_param->fore_aft_hits.hit_lists[ list_count ].twoStepMap.mucm_1,
		    Ein_param->fore_aft_hits.hit_lists[ list_count ].twoStepMap.mucm_2 );		  
     }
  }
}
// ----------- Tstep::two_step::Eout_ladder --------------
// This routine uses the this angular Legendre expansion and the
// next to calculate the contribution to the E_out boxes of the
// transfer matrix between incident energies Ein_param->data_E_0 and
// Ein_param->data_E_1.
void Tstep::two_step::Eout_ladder( Trf::T_matrix& transfer,
   Tstep::two_step_param *Ein_param )
{
  bool geom_OK;  // for checking the consistency of the geometry

  // loop through the outgoing energies (column of transfer)
  for( int Eout_count = 0; Eout_count < transfer.num_Eout_bins;
    ++Eout_count )
  {
    std::vector< double >::const_iterator Eout_ptr = transfer.out_groups.begin( )
      + Eout_count;
    Ein_param->Eout_ptr = Eout_ptr;
    
    // how does the mucm_1 = mucm_2 = -1 curve meet this E-E' box?
    geom_OK = Ein_param->fore_aft_hits.hit_lists[ 0 ].hit_box( -1.0, Eout_ptr,
      Ein_param->data_E_0, Ein_param->data_E_1 );
    if( !geom_OK )
    {
      Msg::Warning( "Tstep::two_step::Eout_ladder",
		    "Check the coding, mucm_1 = mucm_2 = -1" );
      continue;
    }
    if( ( Eout_count < transfer.num_Eout_bins - 1 ) &&
        ( Ein_param->fore_aft_hits.hit_lists[ 0 ].is_above( ) ) )
    {
      // go on to the next E-E' box
      continue;
    }
    // how does the mucm_1 = mucm_2 = 1 curve meet this E-E' box?
    geom_OK = Ein_param->fore_aft_hits.hit_lists[ 3 ].hit_box( 1.0, Eout_ptr,
      Ein_param->data_E_0, Ein_param->data_E_1 );
    if( !geom_OK )
    {
      Msg::Warning( "Tstep::two_step::Eout_ladder",
		    "Check the coding, mucm_1 = mucm_2 = 1" );
      continue;
    }
    if( ( Eout_count > 0 ) && ( Ein_param->fore_aft_hits.hit_lists[ 3 ].is_below( ) ) )
    {
      // we are done with this pair of eta values
      break;
    }
    
    // how does the mucm_1 = -1, mucm_2 = 1 curve meet this E-E' box?
    geom_OK = Ein_param->fore_aft_hits.hit_lists[ 1 ].hit_box( -1.0, Eout_ptr,
      Ein_param->data_E_0, Ein_param->data_E_1 );
    if( !geom_OK )
    {
      Msg::Warning( "Tstep::two_step::Eout_ladder",
		    "Check the coding, mucm_1 = -1, mucm_2 = 1" );
      continue;
    }
    // how does the mucm_1 = 1, mucm_2 = -1 curve meet this E-E' box?
    geom_OK = Ein_param->fore_aft_hits.hit_lists[ 2 ].hit_box( -1.0, Eout_ptr,
      Ein_param->data_E_0, Ein_param->data_E_1 );
    if( !geom_OK )
    {
      Msg::Warning( "Tstep::two_step::Eout_ladder",
		    "Check the coding, mucm_1 = 1, mucm_2 = -1" );
      continue;
    }
    
    // integrate over this E-E' box
    one_Ebox( transfer, Eout_count, Ein_param );
  }
}
// ----------- Tstep::two_step::one_Ebox --------------
// Integrate over one E-E' box
void Tstep::two_step::one_Ebox( Trf::T_matrix& transfer, int Eout_count,
   Tstep::two_step_param *Ein_param )
{
  // the E' energy range
  Ein_param->Eout_ptr = transfer.out_groups.begin( ) + Eout_count;
  Ein_param->Eout_min = transfer.out_groups[ Eout_count ];
  Ein_param->Eout_max = transfer.out_groups[ Eout_count + 1 ];

  // set up common incident energies
  Ein_param->fore_aft_hits.common_hits( );

  bool debug = false;
  if( debug )
  {
    std::cout << "mucm_1 = mucm_2 = -1" << std::endl;
    Ein_param->fore_aft_hits.hit_lists[0].print( );
    std::cout << "mucm_1 =  -1, mucm_2 = 1" << std::endl;
    Ein_param->fore_aft_hits.hit_lists[1].print( );
    std::cout << "mucm_1 = 1, mucm_2 = -1" << std::endl;
    Ein_param->fore_aft_hits.hit_lists[2].print( );
    std::cout << "mucm_1 = mucm_2 = 1" << std::endl;
    Ein_param->fore_aft_hits.hit_lists[3].print( );
  }
  
  // integrate depending on how the curves mucm = const meet the box
  Twohit::two_step_hit_list::iterator hit_ptr = Ein_param->fore_aft_hits.hit_lists[ 0 ].begin( );
  Twohit::two_step_hit_list::iterator next_ptr = hit_ptr;
  ++next_ptr;
  for( ; next_ptr != Ein_param->fore_aft_hits.hit_lists[ 0 ].end( );
       hit_ptr = next_ptr, ++next_ptr )
  {
    double E_mid = 0.5 * ( hit_ptr->E_in + next_ptr->E_in );
    
    double aft_aft_E_out;
    double aft_fore_E_out;
    double fore_aft_E_out;
    double fore_fore_E_out;

    if( use_relativistic )
    {
      Rel::two_step_relativistic_map relTwoStepMap;
      relTwoStepMap.twoStepMasses = &twoStepMasses;

      aft_aft_E_out = relTwoStepMap.two_step_get_E_lab( E_mid, -1.0, -1.0 );
      aft_fore_E_out = relTwoStepMap.two_step_get_E_lab( E_mid, -1.0, 1.0 );
      fore_aft_E_out = relTwoStepMap.two_step_get_E_lab( E_mid, 1.0, -1.0 );
      fore_fore_E_out = relTwoStepMap.two_step_get_E_lab( E_mid, 1.0, 1.0 );
    }
    else
    {
      aft_aft_E_out = Newton_map.two_step_get_E_lab( E_mid, -1.0, -1.0 );
      aft_fore_E_out = Newton_map.two_step_get_E_lab( E_mid, -1.0, 1.0 );
      fore_aft_E_out = Newton_map.two_step_get_E_lab( E_mid, 1.0, -1.0 );
      fore_fore_E_out = Newton_map.two_step_get_E_lab( E_mid, 1.0, 1.0 );
    }

    if( ( aft_aft_E_out >= Ein_param->Eout_max ) ||
	( fore_fore_E_out <= Ein_param->Eout_min ) )
    {
      // do nothing---we are above or below the E-E' box
      continue;
    }

    // What is the range of mucm_1?
    if( fore_aft_E_out <= Ein_param->Eout_max )
    {
      // the lower mucm_1 = 1 curve is inside the E-E' box
      Ein_param->use_mucm1_max = true;
    }
    else
    {
      Ein_param->use_mucm1_max = false;
    }
    
    if( aft_fore_E_out >= Ein_param->Eout_min )
    {
      // the lower mucm_1 = -1 curve is inside the E-E' box
      Ein_param->use_mucm1_min = true;
    }
    else
    {
      Ein_param->use_mucm1_min = false;
    }
       
    // the range of integration in incident energy
    Ein_param->Ein_0 = hit_ptr->E_in;
    Ein_param->Ein_1 = next_ptr->E_in;
    update_T( transfer, Eout_count, Ein_param );
  }
}
// ----------- Tstep::two_step::update_T --------------
void Tstep::two_step::update_T( Trf::T_matrix &transfer, int Eout_count,
   Tstep::two_step_param *Ein_param )
{
  static double tol = Global.Value( "quad_tol" );
  // a vector to store the integrals
  Coef::coef_vector value( transfer.order, transfer.conserve );
  value.set_zero( );

  // parameters for the integration
  Qparam::QuadParamBase *params = static_cast< Qparam::QuadParamBase* >( Ein_param );

  // loop over the cross section data
  Ein_param->this_sigma = Ein_param->first_ladder_sigma;
  Ein_param->next_sigma = Ein_param->this_sigma;
  ++Ein_param->next_sigma;
  // Ein_param->Ein_0 may be past Ein_param->next_sigma
  while( ( Ein_param->this_sigma != Ein_param->last_ladder_sigma ) &&
         ( Ein_param->next_sigma->x < Ein_param->Ein_0 ) )
  {
    Ein_param->this_sigma = Ein_param->next_sigma;
    ++Ein_param->next_sigma;
  }
  for( ; ( Ein_param->this_sigma != Ein_param->last_ladder_sigma ) &&
         ( Ein_param->this_sigma->x <  Ein_param->Ein_1 );
       Ein_param->this_sigma = Ein_param->next_sigma, ++Ein_param->next_sigma )
  {
    double left_E = ( Ein_param->this_sigma->x < Ein_param->Ein_0 ) ?
      Ein_param->Ein_0 : Ein_param->this_sigma->x;
    double right_E = ( Ein_param->next_sigma->x > Ein_param->Ein_1 ) ?
      Ein_param->Ein_1 : Ein_param->next_sigma->x;
    // evaluate the integral
    quad_F::integrate( twoStep_F::Ein_F, transfer.Ein_quad_rule,
       left_E, right_E, params, tol, &value );

    // Add this integral
    transfer( Ein_param->Ein_count, Eout_count ) += value;
    // increment the function counts
    Ein_param->Ein_F_count += Ein_param->func_count;
    ++Ein_param->quad_count;
  }
}

// **************** functions to integrate **********
// Function for the 1-d quadrature
// ---------------- twoStep_F::step2_w_F ------------------
bool twoStep_F::step2_w_F( double w,
   Qparam::QuadParamBase *step2w_param, Coef::coef_vector *value )
{
  // the parameters are really Tstep::step2_w_param
  Tstep::step2_w_param *params =
    static_cast<Tstep::step2_w_param*>( step2w_param );
  params->func_count += 1;

  double mu_lab;
  if( params->use_relativistic )
  {
    mu_lab = params->relTwoStepMap.two_step_get_mu_lab( params->E_in,
	     params->mucm_1, params->mucm_2, w );
  }
  else
  {
    mu_lab = params->map->two_step_get_mu_lab( params->E_in,
      params->mucm_1, params->mucm_2, params->Etrans2, params->Eout_2, w );
  }
  
  // the Legendre polynomials
  math_F::Legendre( mu_lab, value );

  // do the energy weighting if necessary
  if( ( value->conserve == Coef::ENERGY ) || ( value->conserve == Coef::BOTH ) )
  {
    value->scale_E( params->Eout_2 );
  }

  return true;
}
// Function for the 2-d quadrature over mucm2 and w
// ---------------- twoStep_F::mucm_2_F ------------------
bool twoStep_F::mucm_2_F( double mucm2,
   Qparam::QuadParamBase *mucm2_param, Coef::coef_vector *value )
{
  // the parameters are really Tstep::mucm_2_param
  Tstep::mucm_2_param *params =
    static_cast<Tstep::mucm_2_param*>( mucm2_param );
  params->func_count += 1;

  // set up the parameters for the integration over w
  Tstep::step2_w_param w_param;

  double Eout_2;
  if( params->use_relativistic )
  {
    Eout_2 = params->relTwoStepMap.two_step_get_E_lab( params->E_in,
		  params->mucm_1, mucm2 );
  }
  else
  {
    Eout_2 = params->map->two_step_get_E_lab( params->E_in,
		  params->mucm_1, mucm2 );
  }
  
  w_param.setup( mucm2, Eout_2, params );

  // evaluate the integral over w
  Qparam::QuadParamBase *WParam = static_cast< Qparam::QuadParamBase* >( &w_param );
  static double tol = Global.Value( "quad_tol" );

  bool is_OK = quad_F::integrate( twoStep_F::step2_w_F, params->w_quad_rule,
	         -0.5*M_PI, 0.5*M_PI, WParam, tol, value );

  // scale by the area of the hemisphere
  *value *= 1.0 / ( 2 * M_PI );
  params->w_F_count += w_param.func_count;

  return is_OK;
}
// Function for the 3-d quadrature over mucm1, mucm2, and w
// ---------------- twoStep_F::mucm_1_F ------------------
bool twoStep_F::mucm_1_F( double mucm1,
    Qparam::QuadParamBase *mucm1Param, Coef::coef_vector *value )
{
  // the parameters are really Tstep::mucm_1_param
  Tstep::mucm_1_param *mucm1_param = static_cast<Tstep::mucm_1_param*>( mucm1Param );
  mucm1_param->func_count += 1;
  //  if( mucm1_param->func_count % 100 == 0 )
  //  {
  //    Msg::Info( "twoStep_F::mucm_1_F",
  //          Msg::pastenum( "got ", mucm1_param->func_count ) + " evaluations");
  //  }
  // get Eout_lab and mu_lab
  double Eout_lab;
  double mu_lab;

  if( mucm1_param->use_relativistic )
  {
    /*
    We could get cosh(chi_2) and sinh(chi_2) here.
    */
  }
  else
  {
    double E_in = mucm1_param->E_in;
    mucm1_param->map->two_body_get_E_mu_lab( E_in, mucm1,
			       &Eout_lab, &mu_lab );

    mucm1_param->Etrans2 = mucm1_param->map->two_body_get_Etrans2(
				     Eout_lab );
  }


  // evaluate the integral over mucm2 and w
  bool is_OK = mucm1_param->integrate( mucm1, value );

  static bool negative_flag = false;

  // the energy-angle probability density
  double Prob = mucm1_param->coefs.sum_Legendre( mucm1 );
  if(( Prob < 0.0 ) && !negative_flag )
  {
    Msg::Warning( "twoStep_F::mucm_1_F::mu_F",
		  Msg::pastenum( "probability negative for mu: ", mu_lab ) +
        Msg::pastenum( "  mucm1: ", mucm1 ) +
		  Msg::pastenum( " Ein: ", mucm1_param->E_in) );
    negative_flag = true;
  }
  *value *= Prob;

  return is_OK;
}
// ---------------- twoStep_F::Ein_F ------------------
bool twoStep_F::Ein_F( double E_in, Qparam::QuadParamBase *e_quad_param,
   Coef::coef_vector *value )
// Function for the 4-d quadrature
{
  // the parameters are really Tstep::two_step_param *
  Tstep::two_step_param *e_params =
    static_cast<Tstep::two_step_param *>( e_quad_param );
  e_params->func_count += 1;
  //  if( e_params->func_count % 100 == 0 )
  //  {
  //    Msg::Info( "twoStep_F::Ein_F", Msg::pastenum( "got ",
  //      e_params->func_count ) + " evaluations");
  //  }

  // The value of Ein_F is itself an integral over mucm_1, mucm_2, and w.
  // *value comes in as 0.  

  // parameters for the integration over mucm_1, mucm_2, and w.
  Tstep::mucm_1_param mucm1_params;
  mucm1_params.setup( E_in, e_params );
  mucm1_params.coefs.set_E_in( E_in );

  // interpolate the (mucm, probability) with respect to incident energy
  int max_order = ( e_params->left_data->order > e_params->right_data->order ) ?
    e_params->left_data->order : e_params->right_data->order;
  mucm1_params.coefs.initialize( max_order );
  bool is_OK = mucm1_params.interpolate( E_in, e_params->left_data, e_params->right_data );
  if( !is_OK )
  {
    return false;
  }
  
  // the range of integration
  double mu_cm_0;
  double mu_cm_1;
  if( e_params->use_relativistic )
  {
    // parameters for get_mucm1
    Rel::relativistic_2_step_param rel_2_step_param;
    rel_2_step_param.relTwoStepMap.twoStepMasses =
      e_params->relTwoStepMap.twoStepMasses;
    
    mu_cm_0 = ( e_params->use_mucm1_min ) ? -1.0 :
      rel_2_step_param.get_mucm1( E_in, 1, e_params->Eout_min );
    mu_cm_1 = ( e_params->use_mucm1_max ) ? 1.0 :
      rel_2_step_param.get_mucm1( E_in, -1, e_params->Eout_max );
  }
  else
  {
    mu_cm_0 = ( e_params->use_mucm1_min ) ? -1.0 :
      mucm1_params.map->get_mucm1( E_in, 1, e_params->Eout_min );
    mu_cm_1 = ( e_params->use_mucm1_max ) ? 1.0 :
      mucm1_params.map->get_mucm1( E_in, -1, e_params->Eout_max );
  }

  // evaluate the integral over mucm
  Qparam::QuadParamBase *params = static_cast< Qparam::QuadParamBase* >( &mucm1_params );
  static double tol = Global.Value( "quad_tol" );

  is_OK = quad_F::integrate( twoStep_F::mucm_1_F, e_params->mucm1_quad_rule,
         mu_cm_0, mu_cm_1, params, tol, value );

  e_params->num_negative += mucm1_params.num_negative;
  e_params->mucm1_F_count += mucm1_params.func_count;
  e_params->mucm2_F_count += mucm1_params.mucm2_F_count;
  e_params->w_F_count += mucm1_params.w_F_count;
  
  // weight it by flux * cross section
  e_params->set_weight( E_in );
  *value *= e_params->current_weight;
  //  std::cout << "E_in: " << E_in << " mu_cm_0: " << mu_cm_0 << " mu_cm_1: " <<
  //    mu_cm_1 << std::endl;
  //  value->print( );

  return is_OK;
}
