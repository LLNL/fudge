/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: uncorrelated.cpp 1 2006-02-02 03:06:56Z hedstrom $
 *
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// implement the classes used for uncorrelated energy-angle distributions

#ifdef _OPENMP
 #include <omp.h>
#endif

#include "uncorrelated.hpp"
#include "adapt_quad.hpp"
#include "messaging.hpp"
#include "global_params.hpp"

// **************** class Uncor::uncorrelated_param *****************
// ------------- Uncor::uncorrelated_param::setup_direct ------------
// Does truncation or extrapolation of data for direct linlin interpolation
void Uncor::uncorrelated_param::setup_direct( )
{
  // remove previous data
  if( !this_Ein_direct.empty( ) )
  {
    this_Ein_direct.erase( this_Ein_direct.begin( ), this_Ein_direct.end( ) );
    next_Ein_direct.erase( next_Ein_direct.begin( ), next_Ein_direct.end( ) );
  }

  // trancate or extrapolate data?
  static int truncate = Global.Value( "truncate_direct" );
  bool use_truncate = ( truncate > 0 );

  // get the outgoing energy range
  left_data = this_E_dist->begin( );
  right_data = next_E_dist->begin( );

  double Eout_min;
  double left_Eout = left_data->x;
  double right_Eout = right_data->x;
  if( use_truncate )
  {
    Eout_min = ( left_Eout > right_Eout ) ? left_Eout : right_Eout;
  }
  else
  {
    Eout_min = ( left_Eout < right_Eout ) ? left_Eout : right_Eout;
  }

  left_data = this_E_dist->end( );
  --left_data;
  right_data = next_E_dist->end( );
  --right_data;
  double Eout_max;
  left_Eout = left_data->x;
  right_Eout = right_data->x;
  if( use_truncate )
  {
    Eout_max = ( left_Eout < right_Eout ) ? left_Eout : right_Eout;
  }
  else
  {
    Eout_max = ( left_Eout > right_Eout ) ? left_Eout : right_Eout;
  }

  // make copies, truncated or extrapolated
  if( use_truncate )
  {
    // renormalize the copy
    bool low_OK = this_Ein_direct.truncate_copy( *this_E_dist, Eout_min, Eout_max, true );
    bool high_OK = next_Ein_direct.truncate_copy( *next_E_dist, Eout_min, Eout_max, true );
    if( !low_OK || !high_OK )
    {
      Msg::Warning( "standard_Legendre_param::setup_Ein_linlin",
               "truncation gave norm 0, using histogram" );
      // we got norm zero; use histogram in incident energy
      if( !this_Ein_direct.empty( ) )
      {
        this_Ein_direct.erase( this_Ein_direct.begin( ), this_Ein_direct.end( ) );
        next_Ein_direct.erase( next_Ein_direct.begin( ), next_Ein_direct.end( ) );
      }

      this_Ein_direct.copy( *this_E_dist );
      next_Ein_direct.copy( *this_E_dist );
      next_Ein_direct.set_E_in( next_E_dist->get_E_in( ) );
    }
  }
  else
  {
    this_Ein_direct.extrapolate_copy( *this_E_dist, Eout_min, Eout_max );
    next_Ein_direct.extrapolate_copy( *next_E_dist, Eout_min, Eout_max );
  }

  // set pointers at lower incident energy
  left_data = this_Ein_direct.begin( );
  next_left_data = left_data;
  ++next_left_data;
  last_left_data = this_Ein_direct.end( );

  // higher incident energy
  right_data = next_Ein_direct.begin( );
  next_right_data = right_data;
  ++next_right_data;
  last_right_data = next_Ein_direct.end( );
}
// ---------------- Uncor::uncorrelated_param::set_data ------------------
// Sets up the data for interpolation with respect to incident energy
void Uncor::uncorrelated_param:: set_data( )
{
  Ein0_data.set_data( *left_data, *next_left_data, current_prev_Eout,
		      current_next_Eout );
  if( ( Ein_interp.qualifier == Terp::UNITBASE ) ||
      ( ( Ein_interp.qualifier == Terp::DIRECT ) && ( Ein_interp.flag == Terp::LINLIN ) ) )
  {
    Ein1_data.set_data( *right_data, *next_right_data, current_prev_Eout,
		      current_next_Eout );
  }
  else  // DIRECT HISTOGRAM
  {
    Ein1_data.set_pair( Ein0_data.first, Ein0_data.second );
  }
}
// ---------------- Uncor::uncorrelated_param::setup_Ein_cum_prob ------------------
// Sets up the data for cumulative points interpolation in incident energy
void Uncor::uncorrelated_param::setup_Ein_cum_prob( )
{
  left_cum_prob = this_E_dist->cum_prob.begin( );
  next_left_cum_prob = left_cum_prob;
  ++next_left_cum_prob;

  // skip intervals with zero probability
  while( ( left_cum_prob->Prob == 0.0 ) &&
	 ( left_cum_prob->slope == 0.0 ) )
  {
    left_cum_prob = next_left_cum_prob;
    ++next_left_cum_prob;
    left_data = next_left_data;
    ++next_left_data;
  }

  right_cum_prob = next_E_dist->cum_prob.begin( );
  next_right_cum_prob = right_cum_prob;
  ++next_right_cum_prob;

  // skip intervals with zero probability
  while( ( right_cum_prob->Prob == 0.0 ) &&
	 ( right_cum_prob->slope == 0.0 ) )
  {
    right_cum_prob = next_right_cum_prob;
    ++next_right_cum_prob;
    right_data = next_right_data;
    ++next_right_data;
  }

  Ein0_data.first = *left_data;
  Ein1_data.first = *right_data;

  double higher_A = ( next_left_cum_prob->cum_prob < next_right_cum_prob->cum_prob )?
      next_left_cum_prob->cum_prob : next_right_cum_prob->cum_prob;

  setup_high_A( higher_A );

  if( Ein0_data.too_short( ) )
  {
    Ein0_data.short_to_unit_base( higher_A );
  }
  else
  {
    Ein0_data.to_unit_base( );
  }
  
  if( Ein1_data.too_short( ) )
  {
    Ein1_data.short_to_unit_base( higher_A );
  }
  else
  {
    Ein1_data.to_unit_base( );
  }
}
// ---------------- Uncor::uncorrelated_param::setup_high_A ------------------
// Interpolates (Eout, probability) data to the higher common cumulative probability
void Uncor::uncorrelated_param::setup_high_A( double higher_A )
{
  // Do not check the interpolations here
  bool is_OK;
  
  double higher_Eout;

  if( next_left_cum_prob->cum_prob == higher_A )
  {
    Ein0_data.second = *next_left_data;
  }
  else
  {
    higher_Eout = left_cum_prob->get_cum_inv( higher_A );
    Ein0_data.second.x = higher_Eout;
    if( Eout_interp == Terp::HISTOGRAM )
    {
      Ein0_data.second.y = left_data->y;
    }
    else
    {
      Ein0_data.second.y =left_data->linlin_interp( higher_Eout, *next_left_data,
						    &is_OK );
    }
  }

  if( next_right_cum_prob->cum_prob == higher_A )
  {
    Ein1_data.second = *next_right_data;
  }
  else
  {
    higher_Eout = right_cum_prob->get_cum_inv( higher_A );
    Ein1_data.second.x = higher_Eout;
    if( Eout_interp == Terp::HISTOGRAM )
    {
      Ein1_data.second.y = right_data->y;
    }
    else
    {
      Ein1_data.second.y =right_data->linlin_interp( higher_Eout, *next_right_data,
						     &is_OK );
    }
  }
}
// ---------------- Uncor::uncorrelated_param::interp_data_ubase ------------------
// Interpolates to set up for the integration under unitbase interplation over E_in
bool Uncor::uncorrelated_param::interp_data_ubase( double E_in )
{
  bool is_OK = true;  // interpolation is OK
  
  // interpolate the unit-base map and the data to E_in
  double alpha;
  Ddvec::dd_pair mid_Ein_data;
  if( Ein_interp.flag == Terp::LINLIN )
  {
    is_OK = mid_Ein_data.linlin_interp( E_in, Ein0_data, Ein1_data );
    if( !is_OK ) return false;
  }
  else  // Terp::LINLOG
  {
    is_OK = mid_Ein_data.linlog_interp( E_in, Ein0_data, Ein1_data );
    if( !is_OK ) return false;
  }
  // The energy range interpolates lin-lin
  alpha = ( E_in - Ein0_data.get_E_in( ) )/
      ( Ein1_data.get_E_in( ) - Ein0_data.get_E_in( ) );
  if( Ein_interp.qualifier == Terp::UNITBASE )
  {
    mid_ubase_map.interpolate( alpha, this_E_dist->ubase_map,
       next_E_dist->ubase_map );
  }
  else // cumulative points
  {
    mid_ubase_map.interpolate( alpha, Ein0_data.ubase_map,
       Ein1_data.ubase_map );
  }

  // the range of integration over E_out
  double Eout_min_ubase;
  double Eout_max_ubase;
  if( Eout_interp == Terp::HISTOGRAM )
  {
    if( use_Eout_min )
    {
      Eout_min_ubase = mid_ubase_map.to_unit_base( Eout_min, &is_OK );
      if( !is_OK ) return false;
      mid_Eout_data.first.x = Eout_min_ubase;
      mid_Eout_data.first.y = mid_Ein_data.first.y;
    }
    else
    {
      mid_Eout_data.first.x = mid_Ein_data.first.x;
      mid_Eout_data.first.y = mid_Ein_data.first.y;
    }
    if( use_Eout_max )
    {
      Eout_max_ubase = mid_ubase_map.to_unit_base( Eout_max, &is_OK );
      if( !is_OK ) return false;
      mid_Eout_data.second.x = Eout_max_ubase;
      mid_Eout_data.second.y = mid_Ein_data.first.y;
    }
    else
    {
      mid_Eout_data.second.x = mid_Ein_data.second.x;
      mid_Eout_data.second.y = mid_Ein_data.first.y;
    }
  }
  else  // Terp::LINLIN
  {
    if( use_Eout_min )
    {
      Eout_min_ubase = mid_ubase_map.to_unit_base( Eout_min, &is_OK );
      if( !is_OK ) return false;
      mid_Eout_data.first.x = Eout_min_ubase;
      mid_Eout_data.first.y = mid_Ein_data.value( Eout_min_ubase, &is_OK );
      if( !is_OK ) return false;
    }
    else
    {
      mid_Eout_data.first.x = mid_Ein_data.first.x;
      mid_Eout_data.first.y = mid_Ein_data.first.y;
    }
    if( use_Eout_max )
    {
      Eout_max_ubase = mid_ubase_map.to_unit_base( Eout_max, &is_OK );
      if( !is_OK ) return false;
      mid_Eout_data.second.x = Eout_max_ubase;
      mid_Eout_data.second.y = mid_Ein_data.value( Eout_max_ubase, &is_OK );
      if( !is_OK ) return false;
    }
    else
    {
      mid_Eout_data.second.x = mid_Ein_data.second.x;
      mid_Eout_data.second.y = mid_Ein_data.second.y;
    }
  }
  // interpolate mu
  is_OK = interpolate_mu( E_in );
  return is_OK;
}
// ---------------- Uncor::uncorrelated_param::interp_data_flat ------------------
// Interpolates to set up for the integration under histogram interplation over E_in
bool Uncor::uncorrelated_param::interp_data_flat( double E_in )
{
  bool is_OK = true;  // interpolation is OK

  // interpolate the data in E_in
  
  // the range of integration over E_out
  if( Eout_interp == Terp::HISTOGRAM )
  {
    mid_Eout_data.first.y = Ein0_data.first.y;
    mid_Eout_data.second.y = Ein0_data.second.y; // not used
    
    if( use_Eout_min )
    {
      mid_Eout_data.first.x = Eout_min;
    }
    else
    {
      mid_Eout_data.first.x = Ein0_data.first.x;
    }

    if( use_Eout_max )
    {
      mid_Eout_data.second.x = Eout_max;
    }
    else
    {
      mid_Eout_data.second.x = Ein0_data.second.x;
    }
  }
  else  // Terp::LINLIN
  {
    if( use_Eout_min )
    {
      mid_Eout_data.first.x = Eout_min;
      mid_Eout_data.first.y = Ein0_data.value( Eout_min, &is_OK );
      if( !is_OK ) return false;
    }
    else
    {
      mid_Eout_data.first.x = Ein0_data.first.x;
      mid_Eout_data.first.y = Ein0_data.first.y;
    }
    if( use_Eout_max )
    {
      mid_Eout_data.second.x = Eout_max;
      mid_Eout_data.second.y = Ein0_data.value( Eout_max, &is_OK );
      if( !is_OK ) return false;
    }
    else
    {
      mid_Eout_data.second.x = Ein0_data.second.x;
      mid_Eout_data.second.y = Ein0_data.second.y;
    }
  }

  mid_mu_integral.set_E_in( E_in );

  // interpolate mu
  is_OK = interpolate_mu( E_in );
  return is_OK;
}
// ---------------- Uncor::uncorrelated_param::interp_linlin_direct -----------------
// Interpolates to set up for the integration under linlin direct interplation
bool Uncor::uncorrelated_param::interp_linlin_direct( double E_in )
{
  bool is_OK = true;
  
  // interpolate the data in E_in
  Ddvec::dd_pair mid_Ein_data;
  if( Ein_interp.flag == Terp::LINLIN )
  {
    is_OK = mid_Ein_data.linlin_interp( E_in, Ein0_data, Ein1_data );
    if( !is_OK ) return false;
  }
  else  // Terp::LINLOG
  {
    is_OK = mid_Ein_data.linlog_interp( E_in, Ein0_data, Ein1_data );
    if( !is_OK ) return false;
  }

  // the range of integration in E_out
  if( Eout_interp == Terp::HISTOGRAM )
  {
    if( use_Eout_min )
    {
      mid_Eout_data.first.x = Eout_min;
      mid_Eout_data.first.y = mid_Ein_data.first.y;
    }
    else
    {
      mid_Eout_data.first.x = mid_Ein_data.first.x;
      mid_Eout_data.first.y = mid_Ein_data.first.y;
    }
    if( use_Eout_max )
    {
      mid_Eout_data.second.x = Eout_max;
      mid_Eout_data.second.y = mid_Ein_data.first.y;
    }
    else
    {
      mid_Eout_data.second.x = mid_Ein_data.second.x;
      mid_Eout_data.second.y = mid_Ein_data.first.y;
    }
  }
  else
  {
    if( use_Eout_min )
    {
      mid_Eout_data.first.x = Eout_min;
      mid_Eout_data.first.y = mid_Ein_data.value( Eout_min, &is_OK );
      if( !is_OK ) return false;
    }
    else
    {
      mid_Eout_data.first.x = mid_Ein_data.first.x;
      mid_Eout_data.first.y = mid_Ein_data.first.y;
    }
    if( use_Eout_max )
    {
      mid_Eout_data.second.x = Eout_max;
      mid_Eout_data.second.y = mid_Ein_data.value( Eout_max, &is_OK );
      if( !is_OK ) return false;
    }
    else
    {
      mid_Eout_data.second.x = mid_Ein_data.second.x;
      mid_Eout_data.second.y = mid_Ein_data.second.y;
    }
  }
  // interpolate mu
  is_OK = interpolate_mu( E_in );
  return is_OK;
}
// ----------- Uncor::uncorrelated_param::interpolate_mu ------------
// Interpolates the angular data
bool Uncor::uncorrelated_param::interpolate_mu( double E_in )
{
  bool is_OK = true;
  
  if( mu_table )
  {
    int next_mu = 1 - prev_mu;
    if( mu_interp == Terp::LINLIN )
    {
      is_OK = mid_mu_integral.linlin_interp( E_in, mu_integral[ prev_mu ],
        mu_integral[ next_mu ] );
      if( !is_OK ) return false;
    }
    else  // Terp::HISTOGRAM
    {
      mid_mu_integral.only_copy_coef( mu_integral[ prev_mu ] );
    }
  }
  else
  {
    if( mu_interp == Terp::LINLIN )
    {
      is_OK = mid_mu_integral.linlin_interp( E_in, *prev_L_coefs, *next_L_coefs );
      if( !is_OK ) return false;
    }
    else  //  Terp::HISTOGRAM
    {
      mid_mu_integral.only_copy_coef( *prev_L_coefs );
    }
  }

  return true;
}
// ------------- Uncor::uncorrelated_param::next_Eout_cum_prob ----
// Go to the next set of (E_out, probability) pairs for cumulative
// points interpolation in incident energy.
bool Uncor::uncorrelated_param::next_Eout_cum_prob( )
{
  // undo the unit-base map before testing data
  Ein0_data.un_unit_base( );
  Ein1_data.un_unit_base( );

  Ein0_data.first = Ein0_data.second;
  Ein1_data.first = Ein1_data.second;

  // ignore intervals with probability less than skip_tol
  static double skip_tol = Global.Value( "tight_tol" );

  // update the pointers
  if( next_left_cum_prob->E_out <= Ein0_data.first.x )
  {
    left_cum_prob = next_left_cum_prob;
    ++next_left_cum_prob;
    left_data = next_left_data;
    ++next_left_data;
    if( next_left_data == this_E_dist->end( ) )
    {
      return true;
    }
    // skip intervals with essentially zero probability
    while( next_left_cum_prob->cum_prob - left_cum_prob->cum_prob <= skip_tol )
    {
      left_cum_prob = next_left_cum_prob;
      ++next_left_cum_prob;
      left_data = next_left_data;
      ++next_left_data;
      if( next_left_data == this_E_dist->end( ) )
      {
        return true;
      }
      Ein0_data.first = *left_data;
    }
  }

  if( next_right_cum_prob->E_out <= Ein1_data.first.x )
  {
    right_cum_prob = next_right_cum_prob;
    ++next_right_cum_prob;
    right_data = next_right_data;
    ++next_right_data;
    if( next_right_data == next_E_dist->end( ) )
    {
      return true;
    }
    // skip intervals with essentially zero probability
    while( next_right_cum_prob->cum_prob - right_cum_prob->cum_prob <= skip_tol )
    {
      right_cum_prob = next_right_cum_prob;
      ++next_right_cum_prob;
      right_data = next_right_data;
      ++next_right_data;
      if( next_right_data == next_E_dist->end( ) )
      {
        return true;
      }
      Ein1_data.first = *right_data;
    }
  }

  double lower_A = ( left_cum_prob->cum_prob > right_cum_prob->cum_prob )?
      left_cum_prob->cum_prob : right_cum_prob->cum_prob;

  double higher_A = ( next_left_cum_prob->cum_prob < next_right_cum_prob->cum_prob )?
      next_left_cum_prob->cum_prob : next_right_cum_prob->cum_prob;

  setup_high_A( higher_A );
  double dA = higher_A - lower_A;

  if( Ein0_data.too_short( ) )
  {
    Ein0_data.short_to_unit_base( dA );
  }
  else
  {
    Ein0_data.to_unit_base( );
  }
  
  if( Ein1_data.too_short( ) )
  {
    Ein1_data.short_to_unit_base( dA );
  }
  else
  {
    Ein1_data.to_unit_base( );
  }

  // Reset the physical E_out ranges
  Eout_0_range.x = Ein0_data.ubase_map.un_unit_base( Ein0_data.first.x );
  Eout_1_range.x = Ein1_data.ubase_map.un_unit_base( Ein1_data.first.x );
  Eout_0_range.y = Ein0_data.ubase_map.un_unit_base( Ein0_data.second.x );
  Eout_1_range.y = Ein1_data.ubase_map.un_unit_base( Ein1_data.second.x );
  return false;
}

//*************** class Uncor::uncorrelated ****************
// ----------- Uncor::uncorrelated::constructor --------------
Uncor::uncorrelated::uncorrelated( )
{
  Ein_interp.qualifier = Terp::UNITBASE;
  Ein_interp.flag = Terp::LINLIN;
  Eout_interp = Terp::LINLIN;
  mu_interp = Terp::LINLIN;
}
// ----------- Uncor::uncorrelated::destructor --------------
Uncor::uncorrelated::~uncorrelated( )
{
}
// ------------------ Uncor::uncorrelated::read_data --------------
// Reads the ENDL data
void Uncor::uncorrelated::read_data( Dpar::data_parser& infile, int num_I4,
			      Adist::angle_dist *ang )
{
  // the angular data
  angles = ang;
  // do some checking
  if( mu_table && ( angles->threshold < 0.0 ) )
  {
    Msg::FatalError( "Uncor::uncorrelated::read_data",
		     "There is no angular data." );
  }

  Uncor::uncorrelated::iterator new_energy_ptr;
  for( int Ein_count = 0; Ein_count < num_I4; ++Ein_count )
  {
    // make a new energy distribution
    new_energy_ptr = insert( end( ), Ebase::Eprob_vector( ) );
    new_energy_ptr->set_E_in( infile.get_next_double( ) );
    new_energy_ptr->interp_type = Eout_interp;
    // read the (energy, probability density) pairs
    int num_Eout = infile.get_next_int( );
    for( int Eout_count = 0; Eout_count < num_Eout; ++Eout_count )
    {
      double E_out = infile.get_next_double( );
      double Prob = infile.get_next_double( );
      new_energy_ptr->add_entry( E_out, Prob );
    }
    // ensure that the norm is 1
    new_energy_ptr->renorm( false );

    // the map to unit base is done in get_T
    if( Ein_interp.qualifier == Terp::CUMULATIVE_POINTS )
    {
      new_energy_ptr->form_cum_prob( );
    }
  }
}
// ------------------ Uncor::uncorrelated::read_Legendre --------------
// Reads the Legendre coefficients
void Uncor::uncorrelated::read_Legendre( Dpar::data_parser& infile, int num_Ein )
{
  Legendre_coef_data.read_data( infile, num_Ein );
}
// ----------- Uncor::uncorrelated::get_Ein_range --------------
//  Gets the range of nontrivial incident energy bins; computes E_first, first_Ein and last_Ein
// returns true if the threshold is too high for the energy bins
bool Uncor::uncorrelated::get_Ein_range( const Ddvec::dd_vector& sigma, const Ddvec::dd_vector& mult,
    const Ddvec::dd_vector& weight,
    const Lgdata::Flux_List& e_flux, const Egp::Energy_groups& Ein_groups )
{
  double E_last;

  Uncor::uncorrelated_param initial_param;
  bool done = initial_param.get_Ein_range( sigma, mult, weight, e_flux,
                                         Ein_groups, &E_first, &E_last );
  if( done ) return true;

  // check the range of incident energies for the probability data
  Uncor::uncorrelated::const_iterator Ein_data_ptr = begin( );
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
// ----------- Uncor::uncorrelated::setup_param --------------
void Uncor::uncorrelated::setup_param( Uncor::uncorrelated_param *Ein_param )
{
  Ein_param->L_order = 0;  // Legendre order of energy data

  for( int j = 0; j < 2; ++j )
  {
    Ein_param->mu_integral[ j ].initialize( Ein_param->order );
  }
  Ein_param->mid_mu_integral.initialize( Ein_param->order );

  Ein_param->Ein_interp = Ein_interp;
  Ein_param->Eout_interp = Eout_interp;
  Ein_param->mu_interp = mu_interp;
  Ein_param->mid_Eout_data.Eout_interp = Eout_interp;
  Ein_param->mu_table = mu_table;  // is the angular data a table or Legendre?

  Ein_param->this_E_dist = begin( );
  Ein_param->next_E_dist = Ein_param->this_E_dist;
  ++Ein_param->next_E_dist;
  while( Ein_param->next_E_dist->get_E_in( ) <= Ein_param->data_E_0 )
  {
    Ein_param->this_E_dist = Ein_param->next_E_dist;
    ++Ein_param->next_E_dist;
    if( Ein_param->next_E_dist == end( ) )
    {
      Msg::FatalError( "Uncor::uncorrelated::setup_param",
             "incident energy ranges inconsistent" );
    }
  }
  double first_Ein = Ein_param->this_E_dist->get_E_in( );
  if( first_Ein > Ein_param->data_E_0 )
  {
    Ein_param->data_E_0 = first_Ein;
    bool data_bad = Ein_param->update_pointers( first_Ein );
    if( data_bad )
    {
      Msg::FatalError( "Uncor::uncorrelated::setup_param",
		       "energies inconsistent" );
    }
  }

  if( mu_table )  // tabular angular data
  {
    Ein_param->this_mu_dist = angles->begin( );
    Ein_param->next_mu_dist = Ein_param->this_mu_dist;
    ++Ein_param->next_mu_dist;
    while( Ein_param->next_mu_dist->get_E_in( ) <= Ein_param->data_E_0 )
    {
      Ein_param->this_mu_dist = Ein_param->next_mu_dist;
      ++Ein_param->next_mu_dist;
      if( Ein_param->next_mu_dist == angles->end( ) )
      {
        Msg::FatalError( "Uncor::uncorrelated::setup_param",
             "angular data range inconsistent" );
      }
    }
    first_Ein = Ein_param->this_mu_dist->get_E_in( );
  }
  else  // Legendre coefficients
  {
    Ein_param->prev_L_coefs = Legendre_coef_data.begin( );
    Ein_param->next_L_coefs = Ein_param->prev_L_coefs;
    ++Ein_param->next_L_coefs;
    while( Ein_param->next_L_coefs->get_E_in( ) <= Ein_param->data_E_0 )
    {
      Ein_param->prev_L_coefs = Ein_param->next_L_coefs;
      ++Ein_param->next_L_coefs;
      if( Ein_param->next_L_coefs == Legendre_coef_data.end( ) )
      {
        Msg::FatalError( "Uncor::uncorrelated::setup_param",
             "Legendre data range inconsistent" );
      }
    }
    first_Ein = Ein_param->prev_L_coefs->get_E_in( );
  }
  if( first_Ein > Ein_param->data_E_0 )
  {
    Ein_param->data_E_0 = first_Ein;
  }

  if( mu_table )
  {
    Ein_param->prev_mu = 0;
    // do the integrals over mu for the current data
    Ein_param->mu_F_count +=
      to_Legendre_F::from_table( *Ein_param->this_mu_dist,
				 &Ein_param->mu_integral[ 0 ] );
    Ein_param->mu_F_count +=
      to_Legendre_F::from_table( *Ein_param->next_mu_dist,
				 &Ein_param->mu_integral[ 1 ] );

    Ein_param->prev_mu = 0;
    Ein_param->mu_quad_count += Ein_param->this_mu_dist->size( ) +
      Ein_param->next_mu_dist->size( ) - 2;
  }
}
// ----------- Uncor::uncorrelated::set_Ein_range --------------
// Sets the range of incident energies for this integration
void Uncor::uncorrelated::set_Ein_range( Uncor::uncorrelated_param *Ein_param )
{
  // Don't check the interpolation in this routine
  
  Ein_param->set_Ein_range( );
  static double E_tol = Global.Value( "tight_tol" );
  double this_E = Ein_param->this_E_dist->get_E_in( );
  if( this_E > Ein_param->data_E_0 * ( 1.0 + E_tol ) )
  {
    Ein_param->data_E_0 = this_E;
  }
  this_E = Ein_param->next_E_dist->get_E_in( );
  if( this_E < Ein_param->data_E_1 * ( 1.0 - E_tol ) ) Ein_param->data_E_1 = this_E;
  if( Ein_param->data_E_1 < Ein_param->data_E_0 )
  {
    Msg::FatalError( "Uncor::uncorrelated::set_Ein_range",
		     "check the I=4 incident energies" );
  }

  if( mu_table )
  {
    this_E = Ein_param->this_mu_dist->get_E_in( );
  }
  else
  {
    this_E = Ein_param->prev_L_coefs->get_E_in( );
  }
  if( this_E > Ein_param->data_E_0 * ( 1.0 + E_tol ) )
  {
    Ein_param->data_E_0 = this_E;
  }
  if( mu_table )
  {
    this_E = Ein_param->next_mu_dist->get_E_in( );
  }
  else
  {
    this_E = Ein_param->next_L_coefs->get_E_in( );
  }
  if( this_E < Ein_param->data_E_1 * ( 1.0 - E_tol ) )
  {
     Ein_param->data_E_1 = this_E;
  }
  if( Ein_param->data_E_1 < Ein_param->data_E_0 )
  {
    Msg::FatalError( "Uncor::uncorrelated::set_Ein_range",
		     "check the angular incident energies" );
  }
  //  Ein_param->set_sigma_range( );
}
// ----------- Uncor::uncorrelated::Eout_ladder --------------
// This routine uses the energy distributions Ein_param->this_E_dist and the
// next to calculate the contribution to the E_out boxes of the
// transfer matrix between incident energies E_0 and E_1
void Uncor::uncorrelated::Eout_ladder( Trf::T_matrix& transfer,
				       Uncor::uncorrelated_param *Ein_param )
{
  // do truncation or extrapolation of data for direct linlin interpolation
  if( ( Ein_interp.qualifier == Terp::DIRECT ) && ( Ein_interp.flag ==Terp:: LINLIN ) )
  {
    Ein_param->setup_direct( );
  }
  start_Eout( Ein_param );
  
  // Save the physical E_out ranges
  if( Ein_interp.qualifier == Terp::UNITBASE )
  {
    Ein_param->Eout_0_range.x = Ein_param->this_E_dist->ubase_map.Eout_min;
    Ein_param->Eout_1_range.x = Ein_param->next_E_dist->ubase_map.Eout_min;
    Ein_param->Eout_0_range.y =
      Ein_param->this_E_dist->ubase_map.un_unit_base( Ein_param->current_next_Eout );
    Ein_param->Eout_1_range.y =
      Ein_param->next_E_dist->ubase_map.un_unit_base( Ein_param->current_next_Eout );
  }
  else if( Ein_interp.qualifier == Terp::CUMULATIVE_POINTS )
  {
    Ein_param->Eout_0_range.x = Ein_param->Ein0_data.ubase_map.Eout_min;
    Ein_param->Eout_0_range.y = Ein_param->Ein0_data.ubase_map.Eout_max;
    Ein_param->Eout_1_range.x = Ein_param->Ein1_data.ubase_map.Eout_min;
    Ein_param->Eout_1_range.y = Ein_param->Ein1_data.ubase_map.Eout_max;
  }
  else // Terp::DIRECT
  {
    Ein_param->Eout_0_range.x = Ein_param->Ein0_data.first.x;
    Ein_param->Eout_0_range.y = Ein_param->Ein0_data.second.x;
    Ein_param->Eout_1_range.x = Ein_param->Ein1_data.first.x;
    Ein_param->Eout_1_range.y = Ein_param->Ein1_data.second.x;
  }

  // loop through the energy data
  for( ; ; )
  {
    // the physical (E, E') values for the data
    Ein_param->lower_hits.E_Eout.first.x = Ein_param->this_E_dist->get_E_in( );     // E
    Ein_param->lower_hits.E_Eout.first.y = Ein_param->Eout_0_range.x;   // E'
    Ein_param->lower_hits.E_Eout.second.x = Ein_param->next_E_dist->get_E_in( );    // E
    Ein_param->lower_hits.E_Eout.second.y = Ein_param->Eout_1_range.x;  // E'
    Ein_param->upper_hits.E_Eout.first.x = Ein_param->this_E_dist->get_E_in( );     // E
    Ein_param->upper_hits.E_Eout.first.y = Ein_param->Eout_0_range.y;   // E'
    Ein_param->upper_hits.E_Eout.second.x = Ein_param->next_E_dist->get_E_in( );    // E
    Ein_param->upper_hits.E_Eout.second.y = Ein_param->Eout_1_range.y;  // E'

    // loop through the outgoing energies (column of transfer)
    for( int Eout_count = 0; Eout_count < transfer.num_Eout_bins;
      ++Eout_count )
    {
      std::vector< double >::const_iterator Eout_ptr = transfer.out_groups.begin( )
        + Eout_count;
      // how does the lowest unit-base interpolation line meet this E-E' box?
      Ein_param->lower_hits.hit_box( Ein_param->Ein0_data.first.x, Eout_ptr,
        Ein_param->data_E_0, Ein_param->data_E_1 );
      if( ( Eout_count < transfer.num_Eout_bins - 1 ) &&
	  ( Ein_param->lower_hits.is_above( ) ) )
      {
        // go on to the next E-E' box
        continue;
      }
      // how does the next unit-base interpolation line meet this E-E' box?
      Ein_param->upper_hits.hit_box( Ein_param->Ein0_data.second.x, Eout_ptr,
        Ein_param->data_E_0, Ein_param->data_E_1 );
      if( ( Eout_count > 0 ) && ( Ein_param->upper_hits.is_below( ) ) )
      {
        // we are done with this pair of E values
        break;
      }
      // integrate over this E-E' box
      one_Ebox( transfer, Eout_count, Ein_param );
    }

    // go to the next pairs of (E_out, probability)
    bool done;
    if( ( Ein_interp.qualifier == Terp::UNITBASE ) ||
	( ( Ein_interp.qualifier == Terp::DIRECT ) && ( Ein_interp.flag == Terp::LINLIN ) ) )
    {
      done = next_Eout_ubase( Ein_param );
    }
    else if( Ein_interp.qualifier == Terp::CUMULATIVE_POINTS )
    {
      done = Ein_param->next_Eout_cum_prob( );
    }
     else  // DIRECT HISTOGRAM
    {
      done = next_Eout_flat( Ein_param );
    }
    if( done )
    {
      break;
    }
  }
}
// ------------------ Uncor::uncorrelated::start_Eout --------------
// Initializes the pointers to the energy probabilities for this E_in range
void Uncor::uncorrelated::start_Eout( Uncor::uncorrelated_param *Ein_param )
{
  if( ( Ein_interp.qualifier == Terp::DIRECT ) && ( Ein_interp.flag == Terp::LINLIN ) )
  {
    Ein_param->left_data = Ein_param->this_Ein_direct.begin( );
    Ein_param->next_left_data = Ein_param->left_data;
    ++Ein_param->next_left_data;
    Ein_param->last_left_data = Ein_param->this_Ein_direct.end( );

    Ein_param->right_data = Ein_param->next_Ein_direct.begin( );
    Ein_param->next_right_data = Ein_param->right_data;
    ++Ein_param->next_right_data;
    Ein_param->last_right_data = Ein_param->next_Ein_direct.end( );
  }
  else
  {
    Ein_param->left_data = Ein_param->this_E_dist->begin( );
    Ein_param->next_left_data = Ein_param->left_data;
    ++Ein_param->next_left_data;
    Ein_param->last_left_data = Ein_param->this_E_dist->end( );

    Ein_param->right_data = Ein_param->next_E_dist->begin( );
    Ein_param->next_right_data = Ein_param->right_data;
    ++Ein_param->next_right_data;
    Ein_param->last_right_data = Ein_param->next_E_dist->end( );
  }

  if( Ein_interp.qualifier == Terp::UNITBASE )
  {
    // The following coding is safe, because the unit-base outgoing energies are 0 <= E <= 1.
    Ein_param->current_prev_Eout = 0.0;
    double left_Eout = Ein_param->next_left_data->x;
    double right_Eout = Ein_param->next_right_data->x;
    static double etol = Global.Value( "tight_tol" );
    if( left_Eout < right_Eout*(1 + etol ) )
    {
      Ein_param->current_next_Eout = left_Eout;
    }
    else
    {
      Ein_param->current_next_Eout = right_Eout;
    }
  }
  else if( Ein_interp.qualifier == Terp::CUMULATIVE_POINTS )
  {
    Ein_param->setup_Ein_cum_prob( );
  }
  else // DIRECT
  {
    Ein_param->current_prev_Eout = Ein_param->left_data->x;
    Ein_param->current_next_Eout = Ein_param->next_left_data->x;
  }

  // set up the parameters
  Ein_param->Ein0_data.set_E_in( Ein_param->this_E_dist->get_E_in( ) );
  Ein_param->Ein1_data.set_E_in( Ein_param->next_E_dist->get_E_in( ) );
  Ein_param->Ein0_data.Eout_interp = Eout_interp;
  Ein_param->Ein1_data.Eout_interp = Eout_interp;

  if ( Ein_interp.qualifier != Terp::CUMULATIVE_POINTS )
  {
    Ein_param->set_data( );
  }
}
// ----------- Uncor::uncorrelated::next_ladder --------------
bool Uncor::uncorrelated::next_ladder( double E_in, Uncor::uncorrelated_param *Ein_param )
{
  bool done = Ein_param->update_bin_pointers( E_in );
  static double etol = Global.Value( "tight_tol" );
  if( !done )
  {
    double E_tol = E_in * etol;
    if( E_in + E_tol >= Ein_param->next_E_dist->get_E_in( ) )
    {
      while( E_in + E_tol >= Ein_param->next_E_dist->get_E_in( ) )
      {
        // get the next energy data
        Ein_param->this_E_dist = Ein_param->next_E_dist;
        ++Ein_param->next_E_dist;
        if( Ein_param->next_E_dist == end ( ) )
        {
          return true;
        }
      }
    }
  }
  return done;
}
// ----------- Uncor::uncorrelated::next_tabular_ladder --------------
// Get the next incident energy range for tabular angular data
bool Uncor::uncorrelated::next_tabular_ladder( double E_in, Uncor::uncorrelated_param *Ein_param )
{
  bool done = next_ladder( E_in, Ein_param );
  if( !done )
  {
    static double etol = Global.Value( "tight_tol" );
    double E_tol = E_in * etol;
    //    double E_tol = 0.0;
    // go to the next angular distribution?
    if( E_in + E_tol >= Ein_param->next_mu_dist->get_E_in( ) )
    {
      int mu_count = 0;
      while( E_in + E_tol >= Ein_param->next_mu_dist->get_E_in( ) )
      {
        // get the next angular data
        Ein_param->this_mu_dist = Ein_param->next_mu_dist;
        ++Ein_param->next_mu_dist;
        ++mu_count;
        if( Ein_param->next_mu_dist == angles->end ( ) )
        {
          return true;
        }
      }
      if( mu_count == 1 )
      {
        // get the new integral over mu
	Ein_param->mu_F_count +=
          to_Legendre_F::from_table( *Ein_param->next_mu_dist,
				 &Ein_param->mu_integral[ Ein_param->prev_mu ] );

	Ein_param->prev_mu = 1 - Ein_param->prev_mu;
        Ein_param->mu_quad_count += Ein_param->next_mu_dist->size( ) - 1;
      }
      else
      {
        // get both new integrals over mu
	Ein_param->mu_F_count +=
          to_Legendre_F::from_table( *Ein_param->this_mu_dist,
				 &Ein_param->mu_integral[ 0 ] );
        Ein_param->mu_F_count +=
          to_Legendre_F::from_table( *Ein_param->next_mu_dist,
				 &Ein_param->mu_integral[ 1 ] );

        Ein_param->prev_mu = 0;
        Ein_param->mu_quad_count += Ein_param->this_mu_dist->size( ) +
	  Ein_param->next_mu_dist->size( ) - 2;
      }
    }
  }
  return done;
}
// ----------- Uncor::uncorrelated::next_Legendre_ladder --------------
// Get the next incident energy range for Legendre coefficient angular data
bool Uncor::uncorrelated::next_Legendre_ladder( double E_in, Uncor::uncorrelated_param *Ein_param )
{
  bool done = next_ladder( E_in, Ein_param );
  if( !done )
  {
    static double etol = Global.Value( "tight_tol" );
    double E_tol = E_in * etol;
    // go to the next angular distribution?
    if( E_in + E_tol >= Ein_param->next_L_coefs->get_E_in( ) )
    {
      while( E_in + E_tol >= Ein_param->next_L_coefs->get_E_in( ) )
      {
        // get the Legendre coefficients for the next incident energy
        Ein_param->prev_L_coefs = Ein_param->next_L_coefs;
        ++Ein_param->next_L_coefs;
        if( Ein_param->next_L_coefs == Legendre_coef_data.end ( ) )
        {
          return true;
        }
      }
    }
  }
  return done;
}
// ----------- Uncor::uncorrelated::next_Eout_ubase --------------
// go to the next set of (E_out, probability) pairs for unit-base interpolation
bool Uncor::uncorrelated::next_Eout_ubase( Uncor::uncorrelated_param *Ein_param )
{
  Ein_param->current_prev_Eout = Ein_param->current_next_Eout;

  // which outgoing intervals do we increment?
  double left_Eout = Ein_param->next_left_data->x;
  static double etol = Global.Value( "tight_tol" );
  if( left_Eout < Ein_param->current_prev_Eout*(1 + etol ) )
  {
    Ein_param->left_data = Ein_param->next_left_data;
    ++Ein_param->next_left_data;
    if( Ein_param->next_left_data == Ein_param->last_left_data )
    {
      return( true );
    }
  }
  double right_Eout = Ein_param->next_right_data->x;
  if( right_Eout < Ein_param->current_prev_Eout*(1 + etol ) )
  {
    Ein_param->right_data = Ein_param->next_right_data;
    ++Ein_param->next_right_data;
    if( Ein_param->next_right_data == Ein_param->last_right_data )
    {
      return( true );
    }
  }

  // find the upper common outgoing energy
  left_Eout = Ein_param->next_left_data->x;
  right_Eout = Ein_param->next_right_data->x;
  if( left_Eout < right_Eout*(1 + etol ) )
  {
    Ein_param->current_next_Eout = left_Eout;
  }
  else
  {
    Ein_param->current_next_Eout = right_Eout;
  }

  Ein_param->set_data( );

  // Reset the physical E_out ranges
  Ein_param->Eout_0_range.x = Ein_param->Eout_0_range.y;
  Ein_param->Eout_1_range.x = Ein_param->Eout_1_range.y;
  if( Ein_interp.qualifier == Terp::UNITBASE )
  {
    Ein_param->Eout_0_range.y =
      Ein_param->this_E_dist->ubase_map.un_unit_base( Ein_param->current_next_Eout );
    Ein_param->Eout_1_range.y =
      Ein_param->next_E_dist->ubase_map.un_unit_base( Ein_param->current_next_Eout );
  }
  else
  {
    Ein_param->Eout_0_range.y = Ein_param->current_next_Eout;
    Ein_param->Eout_1_range.y = Ein_param->current_next_Eout;
  }
  return false;
}
// ----------- Uncor::uncorrelated::next_Eout_flat --------------
// go to the next set of (E_out, probability) pairs for histograms in Ein
bool Uncor::uncorrelated::next_Eout_flat( Uncor::uncorrelated_param *Ein_param )
{
  Ein_param->current_prev_Eout = Ein_param->current_next_Eout;

  // increment the left data
  Ein_param->left_data = Ein_param->next_left_data;
  ++Ein_param->next_left_data;
  if( Ein_param->next_left_data == Ein_param->this_E_dist->end( ) )
  {
    return( true );
  }

  // find the upper outgoing energy
  double left_Eout = Ein_param->next_left_data->x;
  Ein_param->current_next_Eout = left_Eout;

  Ein_param->set_data( );

  // Reset the physical E_out ranges
  Ein_param->Eout_0_range.x = Ein_param->Eout_0_range.y;
  Ein_param->Eout_1_range.x = Ein_param->Eout_1_range.y;
  Ein_param->Eout_0_range.y = Ein_param->current_next_Eout;
  Ein_param->Eout_1_range.y = Ein_param->current_next_Eout;

  return false;
}
// ------------------ Uncor::uncorrelated::get_T --------------
// Calculates the transfer matrix for this particle.
void Uncor::uncorrelated::get_T( const Ddvec::dd_vector& sigma,
				 const Ddvec::dd_vector& mult,
			       const Ddvec::dd_vector& weight,
				 Trf::T_matrix& transfer )
{
  bool interp_OK = ( ( Ein_interp.qualifier == Terp::UNITBASE ) &&
		     ( ( Ein_interp.flag == Terp::LINLIN ) ||
		       ( Ein_interp.flag == Terp::LINLOG ) ) ) ||
    ( ( Ein_interp.qualifier == Terp::CUMULATIVE_POINTS ) &&
      ( ( Ein_interp.flag == Terp::LINLIN ) ||
        ( Ein_interp.flag == Terp::HISTOGRAM ) ) ) ||
    ( ( Ein_interp.qualifier == Terp::DIRECT ) &&
      ( ( Ein_interp.flag == Terp::LINLIN ) ||
        ( Ein_interp.flag == Terp::HISTOGRAM ) ) );

  if( !interp_OK )
  {
    Msg::FatalError( "Uncor::uncorrelated::get_T",
      "Incident energy interpolation_type not implemented" );
  }
  interp_OK = ( Eout_interp == Terp::LINLIN ) || ( Eout_interp == Terp::HISTOGRAM );
  if( !interp_OK )
  {
    Msg::FatalError( "Uncor::uncorrelated::get_T",
      "Outgoing energy interpolation not implemented not implemented" );
  }
  interp_OK = ( mu_interp == Terp::LINLIN ) || ( mu_interp == Terp::HISTOGRAM );
  if( !interp_OK )
  {
    Msg::FatalError( "Uncor::uncorrelated::get_T",
		     "cosine interpolation not implemented" );
  }

  // first, do a unit-base transformation
  if( Ein_interp.qualifier == Terp::UNITBASE )
  {
    unit_base( 0 );
  }
  bool done = get_Ein_range( sigma, mult, weight, transfer.e_flux, transfer.in_groups );
  if( done )
  {
    transfer.zero_transfer( );
  }
  transfer.threshold = sigma.begin( )->x;

  long int quad_count = 0;  // number of 2-d quadratures
  long int Ein_F_count= 0;  // number of calls to uncorrelated_F::Ein_F
  long int mu_quad_count = 0;  // number of integrations over mu
  long int mu_F_count = 0;  // number of calls to uncorrelated_F::mu_F


  // now do the integrals bin by bin
#pragma omp parallel for schedule( dynamic, 1 ) default( none )	\
  shared( sigma, mult, weight, transfer ) \
  reduction( +: quad_count ) reduction( +: Ein_F_count ) \
  reduction( +: mu_quad_count ) reduction( +: mu_F_count )
  for( int Ein_bin = first_Ein; Ein_bin < last_Ein; ++Ein_bin )
  {
    Uncor::uncorrelated_param Ein_param;
    Ein_param.order = transfer.order;
    // set up the data range for this bin
    Ein_param.setup_bin( Ein_bin, sigma, mult, weight, transfer.e_flux,
                         transfer.in_groups );
    setup_param( &Ein_param );
    for( ; ; )
    {
      set_Ein_range( &Ein_param );   // get the incident energy interval
      Eout_ladder( transfer, &Ein_param );
      bool Done;
      if( mu_table )
      {
        Done = next_tabular_ladder( Ein_param.data_E_1, &Ein_param );   // go to the next interval
      }
      else
      {
        Done = next_Legendre_ladder( Ein_param.data_E_1, &Ein_param );   // go to the next interval
      }
      if( Done )
      {
        break;
      }
    }
    quad_count += Ein_param.quad_count;
    Ein_F_count += Ein_param.Ein_F_count;
    mu_quad_count += Ein_param.mu_quad_count;
    mu_F_count += Ein_param.mu_F_count;
  } // end of parallel loop

  // print the counts of function evaluations
  std::cout << "2d quadratures: " << quad_count << std::endl;
  std::cout << "uncorrelated_F::Ein_F calls: " << Ein_F_count << std::endl;
  std::cout << "average uncorrelated_F::Ein_F calls: " << 1.0*Ein_F_count/quad_count << std::endl;
  std::cout << "quadratures over cosine: " << mu_quad_count << std::endl;
  if( mu_quad_count > 0 )
  {
    std::cout << "uncorrelated_F::mu_F calls: " << mu_F_count << std::endl;
    std::cout << "average uncorrelated_F::mu_F calls: " << 1.0*mu_F_count/mu_quad_count << std::endl;
  }
}
// ----------- Uncor::uncorrelated::one_Ebox --------------
// Integrate over one E-E' box
void Uncor::uncorrelated::one_Ebox( Trf::T_matrix& transfer, int Eout_count,
    Uncor::uncorrelated_param *Ein_param )
{
  // the E' energy range
  Ein_param->Eout_min = transfer.out_groups[ Eout_count ];
  Ein_param->Eout_max = transfer.out_groups[ Eout_count + 1 ];

  // set up common incident energies
  Ein_param->lower_hits.common_hits( Ein_param->upper_hits );

  // integrate depending on how the hyperbolas eta = const meet the box
  Box::energy_hit_list::iterator low_hit_ptr = Ein_param->lower_hits.begin( );
  Box::energy_hit_list::iterator next_low_ptr = low_hit_ptr;
  ++next_low_ptr;
  Box::energy_hit_list::iterator high_hit_ptr = Ein_param->upper_hits.begin( );
  Box::energy_hit_list::iterator next_high_ptr = high_hit_ptr;
  ++next_high_ptr;
  for( ; ( next_low_ptr != Ein_param->lower_hits.end( ) ) &&
         ( next_high_ptr != Ein_param->upper_hits.end( ) );
       low_hit_ptr = next_low_ptr, ++next_low_ptr,
         high_hit_ptr = next_high_ptr, ++next_high_ptr )
  {
    if( ( low_hit_ptr->hit_edge == Box::ABOVE ) ||
        ( low_hit_ptr->hit_edge == Box::TOP_OUT ) )
    {
      // do nothing---we are above the E-E' box
      continue;
    }
    else if( ( low_hit_ptr->hit_edge == Box::TOP_IN ) ||
             ( low_hit_ptr->hit_edge == Box::BOTTOM_IN ) ||
             ( low_hit_ptr->hit_edge == Box::INSIDE ) )
    {
      // the lower eta = const hyperbola is inside the E-E' box
      Ein_param->use_Eout_min = false;
      // where is the upper hyperbola?
      if( ( high_hit_ptr->hit_edge == Box::ABOVE ) ||
          ( high_hit_ptr->hit_edge == Box::TOP_OUT ) )
      {
        // integrate up to the top of the E-E' bin
        Ein_param->use_Eout_max = true;
      }
      else
      {
        // integrate up to the next eta = const hyperbola
        Ein_param->use_Eout_max = false;
      }
    }
    else
    {
      // the lower eta = const hyperbola is below the E-E' box;
      // integrate from Eout_min
      Ein_param->use_Eout_min = true;
      // where is the upper eta = const hyperbola?
      if( ( high_hit_ptr->hit_edge == Box::BOTTOM_OUT ) ||
          ( high_hit_ptr->hit_edge == Box::BELOW ) )
      {
        // do nothing---we are below the E-E' box
        continue;
      }
      else if( ( high_hit_ptr->hit_edge == Box::TOP_IN ) ||
               ( high_hit_ptr->hit_edge == Box::BOTTOM_IN ) ||
               ( high_hit_ptr->hit_edge == Box::INSIDE ) )
      {
        // the upper eta = const hyperbola is inside the E-E' box
        Ein_param->use_Eout_max = false;
      }
      else
      {
        // the upper eta = const hyperbola is above the E-E' box
        Ein_param->use_Eout_max = true;
      }
    }
    // the range of integration in incident energy
    Ein_param->Ein_0 = low_hit_ptr->E_in;
    Ein_param->Ein_1 = next_low_ptr->E_in;
    update_T( transfer, Eout_count, Ein_param );
  }
}
// ----------- Uncor::uncorrelated::update_T --------------
void Uncor::uncorrelated::update_T( Trf::T_matrix &transfer, int Eout_count,
   Uncor::uncorrelated_param *Ein_param )
{
  // a vector to store the integrals, one Legendre order
  Coef::coef_vector value( transfer.order, transfer.conserve );
  value.set_zero( );
  // parameters for the integration
  Qparam::QuadParamBase *params = static_cast< Qparam::QuadParamBase* >( Ein_param );

  double Ein_0 = Ein_param->Ein_0;
  double Ein_1 = Ein_param->Ein_1;
  // loop over the cross section data
  Ein_param->this_sigma = Ein_param->first_ladder_sigma;
  Ein_param->next_sigma = Ein_param->this_sigma;
  ++Ein_param->next_sigma;
  // Ein_param->Ein_0 may be past Ein_param->next_sigma
  while( ( Ein_param->this_sigma != Ein_param->last_ladder_sigma ) &&
         ( Ein_param->next_sigma->x < Ein_0 ) )
  {
    Ein_param->this_sigma = Ein_param->next_sigma;
    ++Ein_param->next_sigma;
  }
  for( ; ( Ein_param->this_sigma != Ein_param->last_ladder_sigma ) &&
         ( Ein_param->this_sigma->x < Ein_1 );
       Ein_param->this_sigma = Ein_param->next_sigma, ++Ein_param->next_sigma )
  {
    Ein_param->Ein_0 = ( Ein_param->this_sigma->x < Ein_0 ) ? Ein_0 :
      Ein_param->this_sigma->x;
    Ein_param->Ein_1 = ( Ein_param->next_sigma->x > Ein_1 ) ? Ein_1 :
      Ein_param->next_sigma->x;

    // evaluate the integral
    static double quad_tol = Global.Value( "quad_tol" );
    quad_F::integrate( uncorrelated_F::Ein_F, transfer.Ein_quad_rule, Ein_param->Ein_0,
		       Ein_param->Ein_1, params, quad_tol, &value );

    if( !std::isfinite( value.weight_1[ 0 ] ) )
    {
      std::cout << "energy_dist::update_T: bad Legendre" << std::endl;
    }

    // add this integral
    transfer( Ein_param->Ein_count, Eout_count ) += value;
    //  std::cout << "E_0: " << E_0 << " E_1: " << E_1 << std::endl;
    //  value.print( );
    // increment the function counts
    Ein_param->Ein_F_count += Ein_param->func_count;
    ++Ein_param->quad_count;
  }
}

// **************** uncorrelated_F::Ein_F ******************
// Function for the quadrature over incident energy
bool uncorrelated_F::Ein_F( double E_in, Qparam::QuadParamBase *void_param,
   Coef::coef_vector *value )
{
  bool is_OK = true;
  
  // the parameters are really Uncor::uncorrelated_param
  Uncor::uncorrelated_param *e_params =
    static_cast<Uncor::uncorrelated_param*>( void_param );
  
  e_params->func_count += 1;
  //  if( e_params->func_count % 500 == 0 )
  //  {
  //    Info( "uncorrelated_F::Ein_F", pastenum( "got ", e_params->func_count ) + " evaluations");
  //  }

  // do the data interpolation
  if( ( e_params->Ein_interp.qualifier == Terp::UNITBASE ) ||
      ( e_params->Ein_interp.qualifier == Terp::CUMULATIVE_POINTS ) )
  {
    is_OK = e_params->interp_data_ubase( E_in );
    if( !is_OK ) return false;
  }
  else  // DIRECT
  {
    if( e_params->Ein_interp.flag == Terp::HISTOGRAM )
    {
      e_params->interp_data_flat( E_in );
    }
    else
    {
      e_params->interp_linlin_direct( E_in );
    }
  }

  double dE_out = e_params->mid_Eout_data.second.x - e_params->mid_Eout_data.first.x;
  if( dE_out <= 0.0 )
  {
    value->set_zero( );
    return false;
  }

  double av_number;
  if( e_params->Eout_interp == Terp::LINLIN )
  {
    av_number = 0.5*dE_out*
     ( e_params->mid_Eout_data.second.y + e_params->mid_Eout_data.first.y );
  }
  else
  {
    av_number = dE_out*e_params->mid_Eout_data.first.y;  // histogram
  }
  if( ( value->conserve == Coef::NUMBER ) || ( value->conserve == Coef::BOTH ) )
  {
    for( int L_count = 0; L_count <= value->order; ++L_count )
    {
      value->weight_1[ L_count ] = av_number*e_params->mid_mu_integral.value( L_count );
    }
  }
  if( ( value->conserve == Coef::ENERGY ) || ( value->conserve == Coef::BOTH ) )
  {
    double phys_Eout_0;
    double phys_Eout_1;
    if( ( e_params->Ein_interp.qualifier == Terp::UNITBASE ) ||
        ( e_params->Ein_interp.qualifier == Terp::CUMULATIVE_POINTS ) )
    {
      phys_Eout_0 = e_params->mid_ubase_map.un_unit_base( e_params->mid_Eout_data.first.x, &is_OK );
      if( !is_OK ) return false;
      phys_Eout_1 = e_params->mid_ubase_map.un_unit_base( e_params->mid_Eout_data.second.x, &is_OK  );
      if( !is_OK ) return false;
    }
    else // Terp::DIRECT
    {
      phys_Eout_0 = e_params->mid_Eout_data.first.x;
      phys_Eout_1 = e_params->mid_Eout_data.second.x;
    }

    double av_E;
    if( e_params->Eout_interp == Terp::LINLIN )
    {
      double d_phys_Eout = ( phys_Eout_1 - phys_Eout_0 )/dE_out;
      double d_prob = ( e_params->mid_Eout_data.second.y - e_params->mid_Eout_data.first.y )/dE_out;
      av_E = 0.5*av_number*( phys_Eout_0 + phys_Eout_1 ) +
        dE_out*dE_out*dE_out*d_phys_Eout*d_prob/12.0;
    }
    else
    {
      av_E = 0.5*av_number*( phys_Eout_0 + phys_Eout_1 );
    }
    for( int L_count = 0; L_count <= value->order; ++L_count )
    {
      value->weight_E[ L_count ] = av_E*e_params->mid_mu_integral.value( L_count );
    }
  }
  e_params->set_weight( E_in );
  *value *= e_params->current_weight;

  return is_OK;
}
