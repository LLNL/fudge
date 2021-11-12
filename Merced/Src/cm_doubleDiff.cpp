/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 601 $
 * $Date: 2017-12-15 $
 * $Author: hedstrom $
 * $Id: cm_doubleDiff.cpp 601 2017-12-15Z hedstrom $
 *1
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// implementation of the classes used to handle pointwise energy-angle
// probability density given in the center-of-mass frame

#ifdef _OPENMP
 #include <omp.h>
#endif

#include "cm_doubleDiff.hpp"
#include "adapt_quad.hpp"
#include "messaging.hpp"
#include "global_params.hpp"

// ************* class cmDD::joint_mu_param **********************
// ----------- cmDD::joint_mu_param::setup --------------
void cmDD::joint_mu_param::setup( double Ein, double Eoutcm,
  cmDD::joint_Eout_param& Ecm_param )
{
  Egeom::Ecm_Elab_mu_param::setup( Ein, Eoutcm, Ecm_param );
  mu_quad_rule = Ecm_param.mu_quad_rule;

  // pointers to the data
  Ein0_Eout0 = Ecm_param.Ein0_Eout0;
  Ein0_Eout1 = Ecm_param.Ein0_Eout1;
  Ein1_Eout0 = Ecm_param.Ein1_Eout0;
  Ein1_Eout1 = Ecm_param.Ein1_Eout1;

  Ein0_Eout0_mu0 = Ein0_Eout0->begin( );
  Ein0_Eout0_mu1 = Ein0_Eout0_mu0;
  ++Ein0_Eout0_mu1;

  Ein0_Eout1_mu0 = Ein0_Eout1->begin( );
  Ein0_Eout1_mu1 = Ein0_Eout1_mu0;
  ++Ein0_Eout1_mu1;

  Ein1_Eout0_mu0 = Ein1_Eout0->begin( );
  Ein1_Eout0_mu1 = Ein1_Eout0_mu0;
  ++Ein1_Eout0_mu1;

  Ein1_Eout1_mu0 = Ein1_Eout1->begin( );
  Ein1_Eout1_mu1 = Ein1_Eout1_mu0;
  ++Ein1_Eout1_mu1;

  E_in0 = Ecm_param.E_in0;
  E_in1 = Ecm_param.E_in1;
  mid_ubase_map = &Ecm_param.mid_ubase_map;
}
// ----------- cmDD::joint_mu_param::scan_data --------------
// Performs the integral over cm mu
bool cmDD::joint_mu_param::scan_data( Coef::coef_vector *value )
{
  bool is_OK = true;
  
  int evaluations = 0;
  start_mu_data( );
  bool done = false;
  value->set_zero( );
  Coef::coef_vector one_arc( value->order, value->conserve );
  
  // loop over the data
  Qparam::QuadParamBase *params = static_cast< Qparam::QuadParamBase* >( this );
  static double tol = Global.Value( "quad_tol" );
  while( !done )
  {
    one_arc.set_zero( );
    bool one_OK = quad_F::integrate( cmDD_F::mu_F, mu_quad_rule, current_data.first.x,
                     current_data.second.x, params, tol, &one_arc );
    if( !one_OK )
    {
      is_OK = false;
    }
    *value += one_arc;
    evaluations += func_count;
    
    done = next_mu_data( );
  }
  func_count = evaluations;
  return is_OK;
}
// ----------- cmDD::joint_mu_param::start_mu_data --------------
// Interpolates to the initial mu data
void cmDD::joint_mu_param::start_mu_data( )
{
  // locate the data for the lower mu value
  find_mu( mu_cm_min );

  // get the probability density at the lower mu value
  double lower_P = get_prob( mu_cm_min );
  current_data.first.x = mu_cm_min;
  current_data.first.y = lower_P;

  // get the probability density at the higher mu value
  set_high_mu( );
}
// ----------- cmDD::joint_mu_param::next_mu_data --------------
// Interpolates to the next set of mu data
// Returns true when the are no more data
bool cmDD::joint_mu_param::next_mu_data( )
{
  bool done;
  current_data.first = current_data.second;

  static double mu_tol = Global.Value( "tight_tol" );
  double prev_mu = current_data.first.x;
  if( prev_mu > mu_cm_max - mu_tol )
  {
    done = true;
  }
  else
  {
    // increment data if necessary
    if( Ein0_Eout0_mu1->x < prev_mu + mu_tol )
    {
      Ein0_Eout0_mu0 = Ein0_Eout0_mu1;
      ++Ein0_Eout0_mu1;
    }
    if( Ein0_Eout1_mu1->x < prev_mu + mu_tol )
    {
      Ein0_Eout1_mu0 = Ein0_Eout1_mu1;
      ++Ein0_Eout1_mu1;
    }
    if( Ein1_Eout0_mu1->x < prev_mu + mu_tol )
    {
      Ein1_Eout0_mu0 = Ein1_Eout0_mu1;
      ++Ein1_Eout0_mu1;
    }
    if( Ein1_Eout1_mu1->x < prev_mu + mu_tol )
    {
      Ein1_Eout1_mu0 = Ein1_Eout1_mu1;
      ++Ein1_Eout1_mu1;
    }
    set_high_mu( );
    done = false;
  }
  return done;
}
// ----------- cmDD::joint_mu_param::find_mu --------------
// Finds the pointers to the data for the current mu value
void cmDD::joint_mu_param::find_mu( double Mu )
{
  // go through the data
  Ein0_Eout0->locate_x( Mu, &Ein0_Eout0_mu0, &Ein0_Eout0_mu1 );
  Ein0_Eout1->locate_x( Mu, &Ein0_Eout1_mu0, &Ein0_Eout1_mu1 );
  Ein1_Eout0->locate_x( Mu, &Ein1_Eout0_mu0, &Ein1_Eout0_mu1 );
  Ein1_Eout1->locate_x( Mu, &Ein1_Eout1_mu0, &Ein1_Eout1_mu1 );
}
// ----------- cmDD::joint_mu_param::get_prob --------------
// Finds the interpolated probability at this mu
double cmDD::joint_mu_param::get_prob( double Mu )
{
  bool interp_OK = true;
  bool one_OK;
  
  // interpolate in mu at the lower incident energy
  Ddvec::dd_pair Ein0_Eout_pair;
  Ein0_Eout_pair.first.x = Ein0_Eout0->get_tag( );
  Ein0_Eout_pair.first.y = Ein0_Eout0_mu0->linlin_interp( Mu,
			  *Ein0_Eout0_mu1, &one_OK );
  if( !one_OK ) interp_OK = false;
  
  Ein0_Eout_pair.second.x = Ein0_Eout1->get_tag( );
  Ein0_Eout_pair.second.y = Ein0_Eout1_mu0->linlin_interp( Mu,
			  *Ein0_Eout1_mu1, &one_OK );
  if( !one_OK ) interp_OK = false;

  // interpolate in mu at the higher incident energy
  Ddvec::dd_pair Ein1_Eout_pair;
  Ein1_Eout_pair.first.x = Ein1_Eout0->get_tag( );
  Ein1_Eout_pair.first.y = Ein1_Eout0_mu0->linlin_interp( Mu,
			  *Ein1_Eout0_mu1, &one_OK );
  if( !one_OK ) interp_OK = false;

  Ein1_Eout_pair.second.x = Ein1_Eout1->get_tag( );
  Ein1_Eout_pair.second.y = Ein1_Eout1_mu0->linlin_interp( Mu,
			  *Ein1_Eout1_mu1, &one_OK );
  if( !one_OK ) interp_OK = false;

  // interpolate in Eout at the lower incident energy
  double Eout_cm_UB = mid_ubase_map->to_unit_base( Eout_cm, &one_OK );
  if( !one_OK ) interp_OK = false;
  Ddvec::dd_pair Ein_pair;
  Ein_pair.first.x = E_in0;
  Ein_pair.first.y = Ein0_Eout_pair.value( Eout_cm_UB, &one_OK );
  if( !one_OK ) interp_OK = false;

  // interpolate in Eout at the higher incident energy
  Ein_pair.second.x = E_in1;
  Ein_pair.second.y = Ein1_Eout_pair.value( Eout_cm_UB, &one_OK );
  if( !one_OK ) interp_OK = false;

  // interpolate in incident energy
  double probability = Ein_pair.value( E_in, &one_OK );
  if( !one_OK ) interp_OK = false;

  if( !interp_OK )
  {
    Msg::FatalError( "cmDD::joint_mu_param::get_prob",
		     "bad interpolation" );
  }
  return probability;
}
// ----------- cmDD::joint_mu_param::set_high_mu --------------
//! Sets the interpolated probability density at the higher mu
void cmDD::joint_mu_param::set_high_mu( )
{
  // determine the next mu
  double next_mu = mu_cm_max;
  if( Ein0_Eout0_mu1->x < next_mu ) next_mu = Ein0_Eout0_mu1->x;
  if( Ein0_Eout1_mu1->x < next_mu ) next_mu = Ein0_Eout1_mu1->x;
  if( Ein1_Eout0_mu1->x < next_mu ) next_mu = Ein1_Eout0_mu1->x;
  if( Ein1_Eout1_mu1->x < next_mu ) next_mu = Ein1_Eout1_mu1->x;

  // interpolate data
  current_data.second.x = next_mu;
  current_data.second.y = get_prob( next_mu );
}

// ************* class cmDD::joint_dist_param ************************
// ----------- cmDD::joint_dist_param::set_Ecm_param --------------
// Sets up the parameters for integration over center-of-mass outgoing energy
void cmDD::joint_dist_param::set_Ecm_param( double E_in )
{
  Ecm_params.map = map;
  Ecm_params.mu_quad_rule = mu_quad_rule;  // rule for integration over cosine
  Ecm_params.Eout_interp = Eout_interp;
  Ecm_params.current_E_in = E_in;
  Ecm_params.E_in0 = Ein0_data->get_E_in( );
  Ecm_params.E_in1 = Ein1_data->get_E_in( );

  double alpha = ( E_in - Ecm_params.E_in0 ) /
    ( Ecm_params.E_in1 - Ecm_params.E_in0 );
  Ecm_params.mid_ubase_map.interpolate( alpha, Ein0_data->Eout_ubase_map,
					Ein1_data->Eout_ubase_map );

  // the energy range of the current data
  double phys_Eout_min = Ecm_params.mid_ubase_map.un_unit_base( prev_data_Eout );
  double phys_Eout_max = Ecm_params.mid_ubase_map.un_unit_base( next_data_Eout );
  Ecm_params.setup( E_in, Eout_min, Eout_max, phys_Eout_min, phys_Eout_max );
  Ecm_params.V_lab_sectors( );
}
// ----------- cmDD::joint_dist_param::start_Eout_cm --------------
// Starts one staircase of the Eout_cm data
void cmDD::joint_dist_param::start_Eout_cm( )
{
  //  Ein0_data.set_E_in( Ein0_data->get_E_in( ) );
  //  Ein1_data.set_E_in( Ein1_data->get_E_in( ) );

  Ecm_params.Ein0_Eout0 = Ein0_data->begin( );
  Ecm_params.Ein0_Eout1 = Ecm_params.Ein0_Eout0;
  ++Ecm_params.Ein0_Eout1;

  Ecm_params.Ein1_Eout0 = Ein1_data->begin( );
  Ecm_params.Ein1_Eout1 = Ecm_params.Ein1_Eout0;
  ++Ecm_params.Ein1_Eout1;

  setup_Ein_ubase( );

}
// ----------- cmDD::joint_dist_param::setup_Ein_ubase --------------
// Sets up the data for unit-base interpolation in incident energy
void cmDD::joint_dist_param::setup_Ein_ubase( )
{
  // The following coding is safe, because the unit-base outgoing energies are 0 <= E <= 1.
  prev_data_Eout = 0.0;
  double left_Eout = Ecm_params.Ein0_Eout1->get_Eout( );
  double right_Eout = Ecm_params.Ein1_Eout1->get_Eout( );
  static double etol = Global.Value( "tight_tol" );
  if( left_Eout < right_Eout*(1 + etol ) )
  {
    next_data_Eout = left_Eout;
  }
  else
  {
    next_data_Eout = right_Eout;
  }
}
// ----------- cmDD::joint_dist_param::next_Ecm --------------
// Gets the next interval of unit-base outgoing energy
bool cmDD::joint_dist_param::next_Ecm( )
{
  // which outgoing intervals do we increment?
  prev_data_Eout = next_data_Eout;
  double left_Eout = Ecm_params.Ein0_Eout1->get_Eout( );
  static double etol = Global.Value( "tight_tol" );
  if( left_Eout < prev_data_Eout*(1 + etol ) )
  {
    Ecm_params.Ein0_Eout0 = Ecm_params.Ein0_Eout1;
    ++Ecm_params.Ein0_Eout1;
    if( Ecm_params.Ein0_Eout1 == Ein0_data->end( ) )
    {
      return( true );
    }
  }
  double right_Eout = Ecm_params.Ein1_Eout1->get_Eout( );
  if( right_Eout < prev_data_Eout*(1 + etol ) )
  {
    Ecm_params.Ein1_Eout0 = Ecm_params.Ein1_Eout1;
    ++Ecm_params.Ein1_Eout1;
    if( Ecm_params.Ein1_Eout1 == Ein1_data->end( ) )
    {
      return( true );
    }
  }

  // find the upper common outgoing energy
  left_Eout = Ecm_params.Ein0_Eout1->get_Eout( );
  right_Eout = Ecm_params.Ein1_Eout1->get_Eout( );
  if( left_Eout < right_Eout*(1 + etol ) )
  {
    next_data_Eout = left_Eout;
  }
  else
  {
    next_data_Eout = right_Eout;
  }

  return( false );
}

// ************* class cmDD::joint_dist *****************************
// ----------- cmDD::joint_dist::setup_map ------------------
// Sets up the map from center-of-mass to laboratory coordinates
void cmDD::joint_dist::setup_map( )
{
  particles.check_data( );
  map.setup_ratios( particles.mProj, particles.mTarg, particles.mProd );
}
// ----------- cmDD::joint_dist::get_Ein_range --------------
//  Gets the range of nontrivial incident energy bins; computes E_first, first_Ein and last_Ein
// returns true if the threshold is too high for the energy bins
bool cmDD::joint_dist::get_Ein_range( const Ddvec::dd_vector& sigma, const Ddvec::dd_vector& mult,
    const Ddvec::dd_vector& weight,
    const Lgdata::Flux_List& e_flux, const Egp::Energy_groups& Ein_groups )
{
  double E_last;

  cmDD::joint_dist_param initial_param;
  bool done = initial_param.get_Ein_range( sigma, mult, weight, e_flux,
                                         Ein_groups, &E_first, &E_last );
  if( done ) return true;

  // check the range of incident energies for the probability data
  cmDD::joint_dist::const_iterator Ein_data_ptr = begin( );
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
// ----------- cmDD::joint_dist::get_T --------------
// Calculates the transfer matrix for this particle.
void cmDD::joint_dist::get_T( const Ddvec::dd_vector& sigma,
			      const Ddvec::dd_vector& multiple,
  const Ddvec::dd_vector& weight, Trf::T_matrix& transfer )
{
  if( ( Ein_interp.qualifier != Terp::UNITBASE ) || ( Ein_interp.flag != Terp::LINLIN ) )
  {
    Msg::FatalError( "cmDD::joint_dist::get_T",
		     "Ein_interp not implemented" );
  }
  if( mu_interp != Terp::LINLIN )
  {
    Msg::FatalError( "cmDD::joint_dist::get_T",
		     "mu_interp not implemented" );
  }
  if( ( Eout_interp.flag != Terp::LINLIN ) && ( Eout_interp.flag != Terp::HISTOGRAM ) )
  {
    Msg::FatalError( "cmDD::joint_dist::get_T",
		     "Eout_interp not implemented" );
  }

  setup_map( );
  bool done = get_Ein_range( sigma, multiple, weight, transfer.e_flux,
    transfer.in_groups );
  if( done )
  {
    transfer.zero_transfer( );
  }
  // This model is used for many different reactions
  transfer.threshold = sigma.begin( )->x;

  long int quad_count = 0;  // number of 2-d quadratures
  long int Ein_F_count= 0;  // number of calls to joint_dist_F::Ein_F
  long int Eout_F_count= 0;  // number of calls to joint_dist_F::Eout_F
  long int mu_F_count = 0;  // number of calls to joint_dist_F::mu_F

  // now do the integrals bin by bin
#pragma omp parallel for schedule( dynamic, 1 ) default( none ) \
  shared( sigma, multiple, weight, transfer ) \
  reduction( +: quad_count ) reduction( +: Ein_F_count ) \
  reduction( +: Eout_F_count ) reduction( +: mu_F_count )
  for( int Ein_bin = first_Ein; Ein_bin < last_Ein; ++Ein_bin )
  {
    cmDD::joint_dist_param Ein_param;
    // set up the data range for this bin
    Ein_param.setup_bin( Ein_bin, sigma, multiple, weight, transfer.e_flux,
                         transfer.in_groups );
    Ein_param.Eout_quad_rule = transfer.Eout_quad_rule;
    Ein_param.mu_quad_rule = transfer.mu_quad_rule;
    setup_param( &Ein_param );
    for( ; ; )
    {
      set_Ein_range( &Ein_param );   // get the incident energy interval
      Eout_data_ladder( transfer, &Ein_param );  // loop over the direction cosines
      bool Done = next_ladder( Ein_param.data_E_1, &Ein_param );   // go to the next interval
      if( Done )
      {
        break;
      }
    }
    quad_count += Ein_param.quad_count;
    Ein_F_count += Ein_param.Ein_F_count;
    Eout_F_count += Ein_param.Eout_F_count;
    mu_F_count += Ein_param.Ecm_params.mu_F_count;
  } // end of parallel loop

  // print the counts of function evaluations
  std::cout << "3d quadratures: " << quad_count << std::endl;
  std::cout << "joint_dist_F::Ein_F calls: " << Ein_F_count << std::endl;
  std::cout << "joint_dist_F::Eout_F calls: " << Eout_F_count << std::endl;
  std::cout << "joint_dist_F::mu_F calls: " << mu_F_count << std::endl;
  std::cout << "average joint_dist_F::Ein_F_count: " << 1.0*Ein_F_count/quad_count << std::endl;
  std::cout << "average joint_dist_F::Eout_F_count: " << 1.0*Eout_F_count/Ein_F_count << std::endl;
  std::cout << "average joint_dist_F::mu_F_count: " << 1.0*mu_F_count/Eout_F_count << std::endl;
}
// ----------- cmDD::joint_dist::setup_param ------------------
// Initializes the quadrature parameters
void cmDD::joint_dist::setup_param( cmDD::joint_dist_param *Ein_param )
{
  static double skip_tol = Global.Value( "tight_tol" );

  Ein_param->Ein_interp = Ein_interp;
  Ein_param->Eout_interp = Eout_interp;
  Ein_param->Ecm_params.Eout_interp = Eout_interp;

  Ein_param->map = &map;

  Ein_param->Ein0_data = begin( );
  Ein_param->Ein1_data = Ein_param->Ein0_data;
  ++Ein_param->Ein1_data;

  while( Ein_param->Ein1_data->get_E_in( ) < Ein_param->data_E_0 *
         ( 1.0 + skip_tol ) )
  {
    Ein_param->Ein0_data = Ein_param->Ein1_data;
    ++Ein_param->Ein1_data;
  }
  double first_Ein = Ein_param->Ein0_data->get_E_in( );
  if( first_Ein > Ein_param->data_E_0 )
  {
    Ein_param->data_E_0 = first_Ein;
    bool data_bad = Ein_param->update_pointers( first_Ein );
    if( data_bad )
    {
      Msg::FatalError( "cmDD::joint_dist::setup_param",
		       "energies inconsistent" );
    }
  }
  // the Vhit::Vcm_Vlab_hit_list objects need the gamma for the energy of translation of the center of mass
  Ein_param->lower_hits.G0_data.gamma = map.gamma;
  Ein_param->lower_hits.G1_data.gamma = map.gamma;
  Ein_param->upper_hits.G0_data.gamma = map.gamma;
  Ein_param->upper_hits.G1_data.gamma = map.gamma;
}
// ----------- cmDD::joint_dist::set_Ein_range ------------------
// Sets the range of incident energies for this intergration
void cmDD::joint_dist::set_Ein_range( cmDD::joint_dist_param *Ein_param )
{
  Ein_param->set_Ein_range( );
  double this_E = Ein_param->Ein0_data->get_E_in( );
  if( this_E > Ein_param->data_E_0 ) Ein_param->data_E_0 = this_E;
  this_E = Ein_param->Ein1_data->get_E_in( );
  if( this_E < Ein_param->data_E_1 ) Ein_param->data_E_1 = this_E;

  if( Ein_param->data_E_1 < Ein_param->data_E_0 )
  {
    Msg::FatalError( "cmDD::joint_dist::set_Ein_range",
		     "check the incident energies" );
  }
}
// ----------- cmDD::joint_dist::next_ladder ------------------
// Go to the next pair of incident energies.  Returns "true" when finished.
bool cmDD::joint_dist::next_ladder( double E_in, cmDD::joint_dist_param *Ein_param )
{
  bool done = Ein_param->update_bin_pointers( E_in );
  static double etol = Global.Value( "tight_tol" );
  double E_tol = E_in * etol;
  if( !done )
  {
    if( E_in + E_tol >= Ein_param->Ein0_data->get_E_in( ) )
    {
      while( E_in + E_tol >= Ein_param->Ein1_data->get_E_in( ) )
      {
        // get the next E_in cmDD::joint_dist data
        Ein_param->Ein0_data = Ein_param->Ein1_data;
        ++Ein_param->Ein1_data;
        if( Ein_param->Ein1_data == end( ) )
        {
          return true;
        }
      }
    }
  }
  return done;
}
// ----------- cmDD::joint_dist::Eout_data_ladder --------------
// Loops through the center-of-mass data for a pair of incident energies
void cmDD::joint_dist::Eout_data_ladder( Trf::T_matrix& transfer,
   cmDD::joint_dist_param *Ein_param )
{
  Ein_param->start_Eout_cm( );
  bool done = false;
  double Ein0 = Ein_param->Ein0_data->get_E_in( );
  double Ein1 = Ein_param->Ein1_data->get_E_in( );
  // loop through the center-of-mass outgoing data
  while( !done )
  {
    double phys_Ecm = Ein_param->Ein0_data->Eout_ubase_map.un_unit_base(
       Ein_param->prev_data_Eout );
    Ein_param->lower_hits.G0_data.set_energies( Ein0, phys_Ecm );

    phys_Ecm = Ein_param->Ein1_data->Eout_ubase_map.un_unit_base(
       Ein_param->prev_data_Eout );
    Ein_param->lower_hits.G1_data.set_energies( Ein1, phys_Ecm );


    phys_Ecm = Ein_param->Ein0_data->Eout_ubase_map.un_unit_base(
        Ein_param->next_data_Eout );
    Ein_param->upper_hits.G0_data.set_energies( Ein0, phys_Ecm );

    phys_Ecm = Ein_param->Ein1_data->Eout_ubase_map.un_unit_base(
       Ein_param->next_data_Eout );
    Ein_param->upper_hits.G1_data.set_energies( Ein1, phys_Ecm );

    lab_Eout_ladder( transfer, Ein_param );

    done = Ein_param->next_Ecm( );
  }
}
// ----------- cmDD::joint_dist::lab_Eout_ladder --------------
// Loops through the outgoing energy bins given cm_Eout data
void cmDD::joint_dist::lab_Eout_ladder( Trf::T_matrix& transfer,
					cmDD::joint_dist_param *Ein_param)
{
  //  bool check_geometry = true;
  bool check_geometry = false;
  bool geom_OK;  // for checking the consistency of the geometry
  bool upper_hits_set = false;
  bool lower_hits_set = false;
  Vhit::Vcm_Vlab_hit_list test_hits;
  test_hits.G0_data.gamma = map.gamma;
  test_hits.G1_data.gamma = map.gamma;
  double left_Ein = Ein_param->Ein0_data->get_E_in( );
  double right_Ein = Ein_param->Ein1_data->get_E_in( );
  double dummy = 0.0;
  double phys_Ecm;
  int Eout_count = 0;
  std::vector< double >::const_iterator Eout_ptr = transfer.out_groups.begin( );
  std::vector< double >::const_iterator next_Eout = Eout_ptr;
  ++next_Eout;
  // Check for only forward emission
  if( Ein_param->upper_hits.G1_data.E_cm <
       Ein_param->upper_hits.G1_data.get_Etrans( ) )
  {
    geom_OK = Ein_param->upper_hits.hit_box( dummy, Eout_ptr,
         Ein_param->data_E_0, Ein_param->data_E_1 );
    upper_hits_set = true;
    if( check_geometry )
    {
      std::cout << "Forward with next_Eout: " << *next_Eout << std::endl;
      Ein_param->upper_hits.print( );
    }
    if( !geom_OK )
    {
      phys_Ecm = Ein_param->Ein0_data->Eout_ubase_map.un_unit_base(
          Ein_param->next_data_Eout );
      test_hits.G0_data.set_energies( left_Ein, phys_Ecm );
      phys_Ecm = Ein_param->Ein1_data->Eout_ubase_map.un_unit_base(
          Ein_param->next_data_Eout );
      test_hits.G1_data.set_energies( right_Ein, phys_Ecm );
      geom_OK = test_hits.hit_box( dummy, Eout_ptr, Ein_param->data_E_0,
          Ein_param->data_E_1 );
      test_hits.print( );
      Msg::FatalError( "cmDD::joint_dist::lab_Eout_ladder",
		       "Check the coding, 1" );
    }
    while( Ein_param->upper_hits.is_above( ) )
    {
      Eout_ptr = next_Eout;
      ++next_Eout;
      if( next_Eout == transfer.out_groups.end( ) )
      {
        return;
      }
      ++Eout_count;
      geom_OK = Ein_param->upper_hits.hit_box( dummy, Eout_ptr,
        Ein_param->data_E_0, Ein_param->data_E_1 );
      if( check_geometry )
      {
        std::cout << "next Forward with next_Eout: " << *next_Eout << std::endl;
        Ein_param->upper_hits.print( );
      }
      if( !geom_OK )
      {
        phys_Ecm = Ein_param->Ein0_data->Eout_ubase_map.un_unit_base(
            Ein_param->next_data_Eout );
        test_hits.G0_data.set_energies( left_Ein, phys_Ecm );
        phys_Ecm = Ein_param->Ein1_data->Eout_ubase_map.un_unit_base(
            Ein_param->next_data_Eout );
        test_hits.G1_data.set_energies( right_Ein, phys_Ecm );
        geom_OK = test_hits.hit_box( dummy, Eout_ptr, Ein_param->data_E_0,
             Ein_param->data_E_1 );
        test_hits.print( );
        Msg::FatalError( "cmDD::joint_dist::lab_Eout_ladder",
			 "Check the coding, 2" );
      }
    }
  }
  else if( Ein_param->lower_hits.G1_data.E_cm >
       Ein_param->lower_hits.G1_data.get_Etrans( ) )
  {
    // Check whether all emission is above the lab energy bin
    geom_OK = Ein_param->lower_hits.hit_box( dummy, Eout_ptr,
        Ein_param->data_E_0, Ein_param->data_E_1 );
    lower_hits_set = true;
    if( check_geometry )
    {
      std::cout << "Backward with E_out: " << *Eout_ptr << std::endl;
      Ein_param->lower_hits.print( );
    }
    if( !geom_OK )
    {
      phys_Ecm = Ein_param->Ein0_data->Eout_ubase_map.un_unit_base(
         Ein_param->prev_data_Eout );
      test_hits.G0_data.set_energies( left_Ein, phys_Ecm );
      phys_Ecm = Ein_param->Ein1_data->Eout_ubase_map.un_unit_base(
           Ein_param->prev_data_Eout );
      test_hits.G1_data.set_energies( right_Ein, phys_Ecm );
      geom_OK = test_hits.hit_box( dummy, Eout_ptr, Ein_param->data_E_0,
           Ein_param->data_E_1 );
      test_hits.print( );
      Msg::FatalError( "cmDD::joint_dist::lab_Eout_ladder",
		       "Check the coding, 3" );
    }
    while( Ein_param->lower_hits.is_below( ) )
    {
      Eout_ptr = next_Eout;
      ++next_Eout;
      if( next_Eout == transfer.out_groups.end( ) )
      {
        return;
      }
      ++Eout_count;
      geom_OK = Ein_param->lower_hits.hit_box( dummy, Eout_ptr,
          Ein_param->data_E_0, Ein_param->data_E_1 );
      if( check_geometry )
      {
        std::cout << "backward with E_out: " << *Eout_ptr << std::endl;
        Ein_param->lower_hits.print( );
      }
      if( !geom_OK )
      {
        phys_Ecm = Ein_param->Ein0_data->Eout_ubase_map.un_unit_base(
           Ein_param->prev_data_Eout );
        test_hits.G0_data.set_energies( left_Ein, phys_Ecm );
        phys_Ecm = Ein_param->Ein1_data->Eout_ubase_map.un_unit_base(
             Ein_param->prev_data_Eout );
        test_hits.G1_data.set_energies( right_Ein, phys_Ecm );
        geom_OK = test_hits.hit_box( dummy, Eout_ptr, Ein_param->data_E_0,
            Ein_param->data_E_1 );
        test_hits.print( );
        Msg::FatalError( "cmDD::joint_dist::lab_Eout_ladder",
			 "Check the coding, 4" );
      }
    }
  }
  // Now, compute integrals until the lab energy bin is above the E_cm data
  for( ; Eout_count < transfer.num_Eout_bins;
       ++Eout_count, Eout_ptr = next_Eout, ++next_Eout )
  {
    if( !upper_hits_set )
    {
      geom_OK = Ein_param->upper_hits.hit_box( dummy, Eout_ptr,
           Ein_param->data_E_0, Ein_param->data_E_1 );
      if( check_geometry )
      {
        std::cout << "upper_hits for Eout: " << *Eout_ptr << std::endl;
        Ein_param->upper_hits.print( );
      }
      if( !geom_OK )
      {
        phys_Ecm = Ein_param->Ein0_data->Eout_ubase_map.un_unit_base(
            Ein_param->next_data_Eout );
        test_hits.G0_data.set_energies( left_Ein, phys_Ecm );
        phys_Ecm = Ein_param->Ein1_data->Eout_ubase_map.un_unit_base(
           Ein_param->next_data_Eout );
        test_hits.G1_data.set_energies( right_Ein, phys_Ecm );
        geom_OK = test_hits.hit_box( dummy, Eout_ptr, Ein_param->data_E_0,
           Ein_param->data_E_1 );
        test_hits.print( );
        Msg::FatalError( "cmDD::joint_dist::lab_Eout_ladder",
			 "Check the coding, 5" );
      }
    }
    if( Ein_param->upper_hits.is_below( ) )
    {
      break;  // we are done
    }
    if( !lower_hits_set )
    {
      geom_OK = Ein_param->lower_hits.hit_box( dummy, Eout_ptr,
          Ein_param->data_E_0, Ein_param->data_E_1 );
      if( check_geometry )
      {
        std::cout << "lower_hits for Eout: " << *Eout_ptr << std::endl;
        Ein_param->lower_hits.print( );
      }
      if( !geom_OK )
      {
        phys_Ecm = Ein_param->Ein0_data->Eout_ubase_map.un_unit_base(
            Ein_param->prev_data_Eout );
        test_hits.G0_data.set_energies( left_Ein, phys_Ecm );
        phys_Ecm = Ein_param->Ein1_data->Eout_ubase_map.un_unit_base(
           Ein_param->prev_data_Eout );
        test_hits.G1_data.set_energies( right_Ein, phys_Ecm );
        geom_OK = test_hits.hit_box( dummy, Eout_ptr, Ein_param->data_E_0,
           Ein_param->data_E_1 );
        test_hits.print( );
        Msg::FatalError( "cmDD::joint_dist::lab_Eout_ladder",
			 "Check the coding, 6" );
      }
    }
    // handle the overlap between this cm outgoing energy range and this 
    // lab energy bin
    one_Ebox( transfer, Eout_count, Ein_param );
    upper_hits_set = false;
    lower_hits_set = false;
  }
}
// ----------- cmDD::joint_dist::one_Ebox ------------------
// Does the integration for one Eout_lab annulus between a pair of incident energies
void cmDD::joint_dist::one_Ebox( Trf::T_matrix& transfer, int Eout_count,
  cmDD::joint_dist_param *Ein_param )
{ 
  // the E' energy range
  Ein_param->Eout_min = transfer.out_groups[ Eout_count ];
  Ein_param->Eout_max = transfer.out_groups[ Eout_count + 1 ];
  //  std::cout << Ein_param->Eout_min << " < E_out < " << Ein_param->Eout_max << std::endl;
  // set up common incident energies 
  Ein_param->lower_hits.common_hits( Ein_param->upper_hits );

  // integrate depending on how the arcs E_cm = const meet the box
  Vhit::Vcm_Vlab_hit_list::iterator low_hit_ptr = Ein_param->lower_hits.begin( );
  Vhit::Vcm_Vlab_hit_list::iterator next_low_ptr = low_hit_ptr;
  ++next_low_ptr;
  Vhit::Vcm_Vlab_hit_list::iterator high_hit_ptr = Ein_param->upper_hits.begin( );
  Vhit::Vcm_Vlab_hit_list::iterator next_high_ptr = high_hit_ptr;
  ++next_high_ptr;
  for( ; ( next_low_ptr != Ein_param->lower_hits.end( ) ) &&
         ( next_high_ptr != Ein_param->upper_hits.end( ) );
       low_hit_ptr = next_low_ptr, ++next_low_ptr,
         high_hit_ptr = next_high_ptr, ++next_high_ptr )
  { 
    if( ( ( low_hit_ptr->hit_edge == Box::ABOVE ) && 
          ( high_hit_ptr->hit_edge == Box::ABOVE ) ) ||
        ( ( next_low_ptr->hit_edge == Box::ABOVE ) && 
          ( next_high_ptr->hit_edge == Box::ABOVE ) ) ||
        ( ( low_hit_ptr->hit_edge == Box::ABOVE_FORWARD ) && 
          ( high_hit_ptr->hit_edge == Box::ABOVE_FORWARD ) ) ||
        ( ( next_low_ptr->hit_edge == Box::ABOVE_FORWARD ) && 
          ( next_high_ptr->hit_edge == Box::ABOVE_FORWARD ) ) ||
        ( ( low_hit_ptr->hit_edge == Box::BELOW ) && 
          ( high_hit_ptr->hit_edge == Box::BELOW ) ) ||
        ( ( next_low_ptr->hit_edge == Box::BELOW ) && 
          ( next_high_ptr->hit_edge == Box::BELOW ) ) )
    { 
      continue;
    }
    // the range of integration in incident energy
    Ein_param->Ein_0 = low_hit_ptr->E_in;
    Ein_param->Ein_1 = next_low_ptr->E_in;
    update_T( transfer, Eout_count, Ein_param );
  }
}
// ----------- cmDD::joint_dist::update_T ------------------
// Adds to an element of transfer the integral between the intersections of 2 Eout_cm = const arcs with the Eout_lab box
void cmDD::joint_dist::update_T( Trf::T_matrix &transfer, int Eout_count,
   cmDD::joint_dist_param *Ein_param )
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
    double left_E = ( Ein_param->this_sigma->x < Ein_param->Ein_0 ) ? Ein_param->Ein_0 :
      Ein_param->this_sigma->x;
    double right_E = ( Ein_param->next_sigma->x > Ein_param->Ein_1 ) ? Ein_param->Ein_1 :
      Ein_param->next_sigma->x;
    // evaluate the integral
    quad_F::integrate( cmDD_F::Ein_F, transfer.Ein_quad_rule,
                       left_E, right_E, params, tol, &value );

    // add this integral
    transfer( Ein_param->Ein_count, Eout_count ) += value;
    // increment the function counts
    Ein_param->Ein_F_count += Ein_param->func_count;
    Ein_param->quad_count += Ein_param->Vcm_hit_count;
  }
}

// **************** Functions to integrate *********************
// ------------------- cmDD_F::mu_F ------------------
// Function for the 1-d quadrature over cm cosine
bool cmDD_F::mu_F( double mu, Qparam::QuadParamBase *mu_quad_param,
		   Coef::coef_vector *value )
{
  // the parameters are really cmDD::joint_mu_param
  cmDD::joint_mu_param *mu_params =
    static_cast< cmDD::joint_mu_param* >( mu_quad_param );
  mu_params->func_count += 1;
  //   if( mu_params->func_count % 100 == 0 )
  //   {
  //     Msg::Info( "cmDD_F::mu_F",
  //          Msg::pastenum( "got ", mu_params->func_count ) + " evaluations");
  //   }

  double Eout_lab;
  double mu_lab;
  mu_params->map->get_E_mu_lab( mu_params->E_in, mu_params->Eout_cm, mu,
    &Eout_lab, &mu_lab );

  // the Legendre polynomials
  math_F::Legendre( mu_lab, value );

  // the unit-base energy-angle probability density
  bool is_OK;
  double Prob = mu_params->current_data.value( mu, &is_OK );
  if( !is_OK )
  {
    return false;
  }
  // physical energy-angle probability density
  Prob /= ( mu_params->mid_ubase_map->Eout_max -
	    mu_params->mid_ubase_map->Eout_min );

  *value *= Prob;

  // do the energy weighting if necessary
  if( ( value->conserve == Coef::ENERGY ) || ( value->conserve == Coef::BOTH ) )
  {
    value->scale_E( Eout_lab );
  }

  return true;
}
// ------------------- cmDD_F::Ecm_F ------------------
// Function for the 2-d quadrature over cm cosine and Eout_cm
bool cmDD_F::Ecm_F( double Eout_cm, Qparam::QuadParamBase *Ecm_quad_param,
		    Coef::coef_vector *value )
{
  // the parameters are really cmDD::joint_Eout_param *
  cmDD::joint_Eout_param *Ecm_param =
    static_cast<cmDD::joint_Eout_param *>( Ecm_quad_param );
  Ecm_param->func_count += 1;

  //   if( Ecm_param->func_count % 100 == 0 )
  //   {
  //     Msg::Info( "cmDD_F::Ecm_F",
  //        Msg::pastenum( "got ", Ecm_param->func_count ) + " evaluations");
  //   }

  // The value of cmDD_F::Ecm_F is itself an integral over mu
  // *value comes in as 0.  

  // set up the parameters for the integration over cm cosine
  cmDD::joint_mu_param mu_param;
  mu_param.setup( Ecm_param->current_E_in, Eout_cm, *Ecm_param );

  // evaluate the integral over mu
  bool is_OK = mu_param.scan_data( value );
  Ecm_param->mu_F_count += mu_param.func_count;

  return is_OK;
}
// ------------------- cmDD_F::Ein_F ------------------
// Function for the 3-d quadrature over E_in, and Eout_cm and cm cosine
// The value of cmDD_F::Ein_F is itself an integral over Eout_cm and cm cosine.
bool cmDD_F::Ein_F( double E_in, Qparam::QuadParamBase *Ein_quad_param,
		    Coef::coef_vector *value )
{
  bool is_OK = true;
  
  value->set_zero( );  // initialize to zero
  // the parameters are really cmDD::joint_dist_param *
  cmDD::joint_dist_param *Ein_param =
    static_cast<cmDD::joint_dist_param *>( Ein_quad_param );

  //   if( Ein_param->func_count % 100 == 0 )
  //   {
  //     Msg::Info( "cmDD_F::Ein_F",
  //        Msg::pastenum( "got ", Ein_param->func_count ) + " evaluations");
  //   }
  Ein_param->Vcm_hit_count = 0;   // number of local calls to quad_F::integrate

  // set up the parameters for the integration over Eout_cm and cm cosine
  Ein_param->set_Ecm_param( E_in );
  Qparam::QuadParamBase *params =
    static_cast< Qparam::QuadParamBase* >( &Ein_param->Ecm_params );

  // Integrate over sectors of ( Eout_lab, mu_cm ) space
  Coef::coef_vector one_value( value->order, value->conserve );

  std::list< Vhit::Vcm_quadBox_Hit >::const_iterator this_V_hit =
    Ein_param->Ecm_params.V_cm_limits.begin( );
  std::list< Vhit::Vcm_quadBox_Hit >::const_iterator next_V_hit = this_V_hit;
  ++next_V_hit;
  for( ; next_V_hit != Ein_param->Ecm_params.V_cm_limits.end( );
         this_V_hit = next_V_hit, ++next_V_hit )
  {
    if( ( next_V_hit->V_cm <= Ein_param->Ecm_params.min_V_cm ) ||
        ( this_V_hit->hit_corner == Vhit::V_BELOW ) )
    {
      continue;  // current V_cm values are below
    }
    else if( ( this_V_hit->V_cm >= Ein_param->Ecm_params.max_V_cm ) ||
             ( next_V_hit->hit_corner == Vhit::V_ABOVE ) )
    {
      break;  // all remaining V_cm values are above
    }
    else
    {
      Ein_param->Vcm_hit_min = *this_V_hit;
      Ein_param->Vcm_hit_max = *next_V_hit;
    }

    Ein_param->Ecm_params.min_hit_corner = Ein_param->Vcm_hit_min.hit_corner;
    Ein_param->Ecm_params.max_hit_corner = Ein_param->Vcm_hit_max.hit_corner;
    double tol = Ein_param->Ecm_params.Ecm_range( );

    bool one_OK = quad_F::integrate( cmDD_F::Ecm_F, Ein_param->Eout_quad_rule,
                       Ein_param->Ecm_params.Ecm_min,
                       Ein_param->Ecm_params.Ecm_max, params, tol, &one_value );
    if( !one_OK )
    {
      is_OK = false;
    }
    *value += one_value;
    // we actually want to count the number of 3d integrals
    Ein_param->Vcm_hit_count += 1;
    Ein_param->func_count += 1;
    Ein_param->Eout_F_count += Ein_param->Ecm_params.func_count;
  }
  // weight it by flux * cross section
  Ein_param->set_weight( E_in );
  *value *= Ein_param->current_weight;
  //  std::cout << "E_in: " << E_in << " eta_0: " << eta_0 << " eta_1: " <<
  //    eta_1 << std::endl;
  //  value->print( );

  return is_OK;
}
