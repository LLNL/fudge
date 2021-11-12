/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2011-01-28 (Mon, Oct 31, 2011) $
 * $Author: hedstrom $
 * $Id: joint_dist.cpp 1 2011-10-30 hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// implementation of the classes used to handle ENDL energy-angle probability density

#include <cmath>
#include <cstdlib>
#include <cfloat>
#ifdef _OPENMP
 #include <omp.h>
#endif

#include "joint_dist.hpp"
#include "adapt_quad.hpp"
#include "messaging.hpp"
#include "global_params.hpp"

// ********* class Jdist::joint_dist_hits *********
// ----------- Jdist::joint_dist_hits::set_incident_range --------------
// Sets the range of integration over incident energy
void Jdist::joint_dist_hits::set_incident_range( double E0, double E1 )
{
  E_Eout.first.x = E0;
  E_Eout.second.x = E1;
}
// ----------- Jdist::joint_dist_hits::get_phys_Eout --------------
// Gets the physical outgoing energies at the ends of the working energy interval
void Jdist::joint_dist_hits::get_phys_Eout( )
{
  // check the incident energies
  double Ein_diff = Ein_1 - Ein_0;
  if( Ein_diff <= 0.0 )
  {
    Msg::FatalError( "Jdist::joint_dist_hits::get_phys_Eout",
		"incident energies out of order" );
  }
  // interpolate in incident energy
  double alpha = ( E_Eout.first.x - Ein_0 )/Ein_diff;
  E_Eout.first.y = ( 1.0 - alpha )*Ein0_data->phys_Eout +
    alpha*Ein1_data->phys_Eout;

  alpha = ( E_Eout.second.x - Ein_0 )/Ein_diff;
  E_Eout.second.y = ( 1.0 - alpha )*Ein0_data->phys_Eout +
    alpha*Ein1_data->phys_Eout;
}

// ********* class Jdist::joint_dist_param *********
// ----------- Jdist::joint_dist_param::joint_dist_param --------------
// Constructor
Jdist::joint_dist_param::joint_dist_param( )
{
  quad_count = 0;
  Ein_F_count = 0;
  Eout_F_count = 0;
  mu_F_count = 0;

  // pointers to the data for computing intersections with the quadrature box
  lower_mu0_hits.Ein0_data = &Ein0_data.mu0_Eout0;  // lower Ein, lower mu
  lower_mu1_hits.Ein0_data = &Ein0_data.mu1_Eout0;  // lower Ein, higher mu
  lower_mu0_hits.Ein1_data = &Ein1_data.mu0_Eout0;  // higher Ein, lower mu
  lower_mu1_hits.Ein1_data = &Ein1_data.mu1_Eout0;  // higher Ein, higher mu

  upper_mu0_hits.Ein0_data = &Ein0_data.mu0_Eout1;  // lower Ein, lower mu
  upper_mu1_hits.Ein0_data = &Ein0_data.mu1_Eout1;  // lower Ein, higher mu
  upper_mu0_hits.Ein1_data = &Ein1_data.mu0_Eout1;  // higher Ein, lower mu
  upper_mu1_hits.Ein1_data = &Ein1_data.mu1_Eout1;  // higher Ein, higher mu
}
// ----------- Jdist::joint_dist_param::geometry --------------
// determines the geometry for 2d integration over outgoing cosine and energy
// returns true if the geometry makes sense
bool Jdist::joint_dist_param::geometry( const Jdata::current_data &this_data )
{
  bool geom_OK = true;

  // the physical (mu, E_out) values for the data
  lower_2d_hits.E_Eout.first.x = this_data.mu0_Eout0.UB_mu;     // unit-base mu
  lower_2d_hits.E_Eout.first.y = this_data.mu0_Eout0.phys_Eout;   // physical E_out
  lower_2d_hits.E_Eout.second.x = this_data.mu1_Eout0.UB_mu;    // unit-base mu
  lower_2d_hits.E_Eout.second.y = this_data.mu1_Eout0.phys_Eout;  // physical E_out
  upper_2d_hits.E_Eout.first.x = this_data.mu0_Eout1.UB_mu;     // unit-base mu
  upper_2d_hits.E_Eout.first.y = this_data.mu0_Eout1.phys_Eout;   // physical E_out
  upper_2d_hits.E_Eout.second.x = this_data.mu1_Eout1.UB_mu;    // unit-base mu
  upper_2d_hits.E_Eout.second.y = this_data.mu1_Eout1.phys_Eout;  // physical E_out

  // parameter not used
  double dummy = 0.0;
  // how do the lower physical outgoing energies meet this mu-E_out box?
  lower_2d_hits.hit_box( dummy, mu_params.Eout_bottom, this_data.mu0_Eout0.UB_mu,
		      this_data.mu1_Eout0.UB_mu ); 
  static bool no_lower_warning = true;
  if( lower_2d_hits.is_above( ) && no_lower_warning )
  {
    // something fishy
    Msg::Warning( "Jdist::joint_dist_param::geometry",
		  "all outgoing energies too high" );
    geom_OK = false;
    no_lower_warning = false;
  }

  // how do the higher physical outgoing energies meet this mu-E_out box?
  upper_2d_hits.hit_box( dummy, mu_params.Eout_bottom, this_data.mu0_Eout0.UB_mu,
		      this_data.mu1_Eout0.UB_mu ); 
  static bool no_upper_warning = true;
  if( upper_2d_hits.is_below( ) && no_upper_warning ) 
  {
    // something fishy
    Msg::Warning( "Jdist::joint_dist_param::geometry",
		  "all outgoing energies too low" );
    geom_OK = false;
    no_upper_warning = false;
  }

  // set up common mu values
  lower_2d_hits.common_hits( upper_2d_hits );

  return geom_OK;
}
// ----------- Jdist::joint_dist_param::common_mu_Eout --------------
// Interpolates data to common values of unit-base mu and outgoing energy
void Jdist::joint_dist_param::common_mu_Eout( double UB_Eout,
				       Jdata::E_mu_P_data *Ein0_mu0_data,
  Jdata::E_mu_P_data *Ein0_mu1_data, Jdata::E_mu_P_data *Ein1_mu0_data,
				       Jdata::E_mu_P_data *Ein1_mu1_data )
{
  // set up the Jdata::E_mu_P_data at given incident energy, mu and outgoing energy
  // temporary storage for interpolation in unit-base mu
  Jdata::E_mu_P_data mu0_ENDF;  // for UB_Eout data at lower ENDF mu
  Jdata::E_mu_P_data mu1_ENDF;  // for UB_Eout data at higher ENDF mu

  // handle the first incident energy
  // interpolate in unit-base outgoing energy first
  double phys_mu = this_Ein_ptr->mu_ubase_map.un_unit_base( this_Ein_this_mu->get_mu( ) );
  this_Ein_this_mu->set_E_mu_P_data( UB_Eout, phys_mu, Ein0_mu0_Eout0,
				     Ein0_mu0_Eout1, &mu0_ENDF );
  phys_mu = this_Ein_ptr->mu_ubase_map.un_unit_base( this_Ein_next_mu->get_mu( ) );
  this_Ein_next_mu->set_E_mu_P_data( UB_Eout, phys_mu, Ein0_mu1_Eout0,
				     Ein0_mu1_Eout1, &mu1_ENDF );
 
  // interpolate in mu
  if( lower_mu == mu0_ENDF.UB_mu )
  {
    Ein0_mu0_data->copy( mu0_ENDF );
  }
  else
  {
    // we should not need the returned bool
    mu0_ENDF.mu_interp( lower_mu, mu1_ENDF, Ein0_mu0_data );
  }
  
  if( upper_mu == mu1_ENDF.UB_mu )
  {
    Ein0_mu1_data->copy( mu1_ENDF );
  }
  else
  {
    // we should not need the returned bool
    mu0_ENDF.mu_interp( upper_mu, mu1_ENDF, Ein0_mu1_data );
  }

  // repeat at the higher incident energy
  // interpolate in unit-base outgoing energy first
  phys_mu = next_Ein_ptr->mu_ubase_map.un_unit_base( next_Ein_this_mu->get_mu( ) );
  next_Ein_this_mu->set_E_mu_P_data( UB_Eout, phys_mu, Ein1_mu0_Eout0,
				     Ein1_mu0_Eout1, &mu0_ENDF );
  phys_mu = next_Ein_ptr->mu_ubase_map.un_unit_base( next_Ein_next_mu->get_mu( ) );
  next_Ein_next_mu->set_E_mu_P_data( UB_Eout, phys_mu, Ein1_mu1_Eout0,
				     Ein1_mu1_Eout1, &mu1_ENDF );
 
  // interpolate in mu
  if( lower_mu == mu0_ENDF.UB_mu )
  {
    Ein1_mu0_data->copy( mu0_ENDF );
  }
  else
  {
    // we should not need the returned bool
    mu0_ENDF.mu_interp( lower_mu, mu1_ENDF, Ein1_mu0_data );
  }
  
  if( upper_mu == mu1_ENDF.UB_mu )
  {
    Ein1_mu1_data->copy( mu1_ENDF );
  }
  else
  {
    // we should not need the returned bool
    mu0_ENDF.mu_interp( upper_mu, mu1_ENDF, Ein1_mu1_data );
  }
}
// ----------- Jdist::joint_dist_param::next_Eout_UB --------------
// go to the next sets of (E_out, probability) pairs, unit base
bool Jdist::joint_dist_param::next_Eout_UB( )
{
  bool done = false;
  double this_Eout;

  // the new lower limit for the integral over Eout is the former upper limit
  lower_Eout = upper_Eout;
  Ein0_data.mu0_Eout0.copy( Ein0_data.mu0_Eout1 );
  Ein0_data.mu1_Eout0.copy( Ein0_data.mu1_Eout1 );
  Ein1_data.mu0_Eout0.copy( Ein1_data.mu0_Eout1 );
  Ein1_data.mu1_Eout0.copy( Ein1_data.mu1_Eout1 );

  static double e_tol = Global.Value( "tight_tol" );
  if( Ein0_mu0_Eout1->x <= lower_Eout*(1 + e_tol ) )
  {
    Ein0_mu0_Eout0 = Ein0_mu0_Eout1;
    ++Ein0_mu0_Eout1;
    if( Ein0_mu0_Eout1 == this_Ein_this_mu->end( ) )
    {
      return true;
    }
  }
  upper_Eout = Ein0_mu0_Eout1->x;

  if( Ein0_mu1_Eout1->x <= lower_Eout*(1 + e_tol ) )
  {
    Ein0_mu1_Eout0 = Ein0_mu1_Eout1;
    ++Ein0_mu1_Eout1;
    if( Ein0_mu1_Eout1 == this_Ein_next_mu->end( ) )
    {
      return true;
    }
  }
  this_Eout = Ein0_mu1_Eout1->x;
  if( upper_Eout > this_Eout )
  {
    upper_Eout = this_Eout;
  }

  if( Ein1_mu0_Eout1->x <= lower_Eout*(1 + e_tol ) )
  {
    Ein1_mu0_Eout0 = Ein1_mu0_Eout1;
    ++Ein1_mu0_Eout1;
    if( Ein1_mu0_Eout1 == next_Ein_this_mu->end( ) )
    {
      return true;
    }
  }
  this_Eout = Ein1_mu0_Eout1->x;
  if( upper_Eout > this_Eout )
  {
    upper_Eout = this_Eout;
  }


  if( Ein1_mu1_Eout1->x <= lower_Eout*(1 + e_tol ) )
  {
    Ein1_mu1_Eout0 = Ein1_mu1_Eout1;
    ++Ein1_mu1_Eout1;
    if( Ein1_mu1_Eout1 == next_Ein_next_mu->end( ) )
    {
      return true;
    }
  }
  this_Eout = Ein1_mu1_Eout1->x;
  if( upper_Eout > this_Eout )
  {
    upper_Eout = this_Eout;
  }

  // reset the Jdata::E_mu_P_data for the higher Eout
  common_mu_Eout( upper_Eout,
                  &Ein0_data.mu0_Eout1,
		  &Ein0_data.mu1_Eout1, &Ein1_data.mu0_Eout1,
		  &Ein1_data.mu1_Eout1 );

  // Reset the parameters for the intersection of the physical enegies with the quadrature box
  lower_mu0_hits.get_phys_Eout( );
  lower_mu1_hits.get_phys_Eout( );
  upper_mu0_hits.get_phys_Eout( );
  upper_mu1_hits.get_phys_Eout( );

  return done;
}
// ----------- Jdist::joint_dist_param::next_Eout_CP --------------
// go to the next sets of (E_out, probability) pairs, cumulative points
bool Jdist::joint_dist_param::next_Eout_CP( )
{
  bool done = false;

  // the new lower limit for the integral over Eout is the former upper limit
  lower_E_cum = upper_E_cum;

  if( next_cum_data( lower_E_cum, Ein0_mu0_Eout0_CP, Ein0_mu0_Eout1_CP,
		     this_Ein_this_mu->cum_prob.end( ) ) )
  {
    return true;
  } 

  if( next_cum_data( lower_E_cum, Ein0_mu1_Eout0_CP, Ein0_mu1_Eout1_CP,
		     this_Ein_next_mu->cum_prob.end( ) ) )
  {
    return true;
  } 

  if( next_cum_data( lower_E_cum, Ein1_mu0_Eout0_CP, Ein1_mu0_Eout1_CP,
		     next_Ein_this_mu->cum_prob.end( ) ) )
  {
    return true;
  } 

  if( next_cum_data( lower_E_cum, Ein1_mu1_Eout0_CP, Ein1_mu1_Eout1_CP,
		     next_Ein_next_mu->cum_prob.end( ) ) )
  {
    return true;
  } 

  set_upper_cum( );
      
  // initialize Ein0_data
  int svar = set_cum_data( Ein0_mu0_Eout0_CP,
			   Ein0_mu1_Eout0_CP,
			   &Ein0_data );

  // initialize Ein1_data
  svar += set_cum_data( Ein1_mu0_Eout0_CP,
			   Ein1_mu1_Eout0_CP,
			   &Ein1_data );

  if( svar >= 6 )
  {
    done = next_Eout_CP( );
  }
  
  // Reset the parameters for the intersection of the physical enegies with the quadrature box
  lower_mu0_hits.get_phys_Eout( );
  lower_mu1_hits.get_phys_Eout( );
  upper_mu0_hits.get_phys_Eout( );
  upper_mu1_hits.get_phys_Eout( );

  return done;
}
// ----------- Jdist::joint_dist_param::next_cum_data --------------
// Increments cumulative probability data for one incident energy and one mu
// Returns true of we run out of data
bool Jdist::joint_dist_param::next_cum_data( double low_cum,
				      Cum::cumulative_prob_list::iterator &prev_CP,
				      Cum::cumulative_prob_list::iterator &next_CP,
				      Cum::cumulative_prob_list::iterator end_CP )
{
  // skip intervals of small probability
  static double CP_skip = Global.Value( "cum_prob_skip" );

  // this test forces a skip of intervals with small probability
  while( next_CP->cum_prob <= low_cum + CP_skip )
  {
    prev_CP = next_CP;
    ++next_CP;
    if( next_CP == end_CP )
    {
      return true;
    }
  }

  return false;
}
// ----------- Jdist::joint_dist_param::set_upper_cum --------------
// Sets the value of upper_E_cum
void Jdist::joint_dist_param::set_upper_cum( )
{
  double Ein0_cum = ( Ein0_mu0_Eout1_CP->cum_prob <
		      Ein0_mu1_Eout1_CP->cum_prob ) ?
	Ein0_mu0_Eout1_CP->cum_prob : Ein0_mu1_Eout1_CP->cum_prob;
      
  double Ein1_cum = ( Ein1_mu0_Eout1_CP->cum_prob <
		      Ein1_mu1_Eout1_CP->cum_prob ) ?
	Ein1_mu0_Eout1_CP->cum_prob : Ein1_mu1_Eout1_CP->cum_prob;
      
  upper_E_cum = ( Ein0_cum < Ein1_cum ) ?
	Ein0_cum : Ein1_cum;
}
// ----------- Jdist::joint_dist_param::set_cum_data --------------
// Sets the data for cumulative-points interpolation for one incident energy.
// returns 3 if both intervals have length zero
int Jdist::joint_dist_param::set_cum_data( Cum::cumulative_prob_list::const_iterator mu0_Eout0_CP,
     Cum::cumulative_prob_list::const_iterator mu1_Eout0_CP,
				     Jdata::current_data *Ein_data )
{
  int svar = 0;
  
  // for interpolation in outgoing energy
  Ddvec::cum_points_pair low_cum_prob;  // entry at lower mu
  Ddvec::cum_points_pair high_cum_prob;  // entry at higher mu

  double low_data_mu = this_Ein_this_mu->get_mu( );
  double high_data_mu = this_Ein_next_mu->get_mu( );

  // set up low_data_mu
  low_cum_prob.set_mu( low_data_mu );
  low_cum_prob.first.x = mu0_Eout0_CP->get_cum_inv( lower_E_cum );
  low_cum_prob.first.y = mu0_Eout0_CP->get_prob( low_cum_prob.first.x );
  low_cum_prob.second.x = mu0_Eout0_CP->get_cum_inv( upper_E_cum );
  low_cum_prob.second.y = mu0_Eout0_CP->get_prob( low_cum_prob.second.x );
  low_cum_prob.ubase_map.Eout_min = low_cum_prob.first.x;
  low_cum_prob.ubase_map.Eout_max = low_cum_prob.second.x;
  if( low_cum_prob.ubase_map.too_short( ) )
  {
    svar += 1;
    low_cum_prob.short_to_unit_base( upper_E_cum - lower_E_cum );
  }
  else
  {
    low_cum_prob.to_unit_base( );
  }
  
  // set up high_data_mu
  high_cum_prob.set_mu( high_data_mu );
  high_cum_prob.first.x = mu1_Eout0_CP->get_cum_inv( lower_E_cum );
  high_cum_prob.first.y = mu1_Eout0_CP->get_prob( high_cum_prob.first.x );
  high_cum_prob.second.x = mu1_Eout0_CP->get_cum_inv( upper_E_cum );
  high_cum_prob.second.y = mu1_Eout0_CP->get_prob( high_cum_prob.second.x );
  high_cum_prob.ubase_map.Eout_min = high_cum_prob.first.x;
  high_cum_prob.ubase_map.Eout_max = high_cum_prob.second.x;
  if( high_cum_prob.ubase_map.too_short( ) )
  {
    svar += 2;
    high_cum_prob.short_to_unit_base( upper_E_cum - lower_E_cum );
  }
  else
  {
    high_cum_prob.to_unit_base( );
  }
  
  // for interpolation in mu
  double denom = high_data_mu - low_data_mu;
  if( denom <= 0.0 )
  {
    Msg::FatalError( "Jdist::joint_dist_param::set_cum_data",
		  "mu values out of order 0" );
  }

  // for lower common mu
  double alpha = ( lower_mu - low_data_mu ) / denom;
  Ein_data->mu0_Eout0.UB_Eout = 0.0;
  Ein_data->mu0_Eout0.Eout_Prob = ( 1.0 - alpha ) * low_cum_prob.first.y +
    alpha * high_cum_prob.first.y;

  Ein_data->mu0_Eout1.UB_Eout = 1.0;
  Ein_data->mu0_Eout1.Eout_Prob = ( 1.0 - alpha ) * low_cum_prob.second.y +
    alpha * high_cum_prob.second.y;

  Ein_data->mu0_ubase_map.interpolate( alpha,
    low_cum_prob.ubase_map, high_cum_prob.ubase_map );

  Ein_data->mu0_Eout0.phys_Eout = Ein_data->mu0_ubase_map.un_unit_base( 0.0 );
  Ein_data->mu0_Eout1.phys_Eout = Ein_data->mu0_ubase_map.un_unit_base( 1.0 );
  
  // for higher common mu
  alpha = ( upper_mu - low_data_mu ) / denom;
  Ein_data->mu1_Eout0.UB_Eout = 0.0;
  Ein_data->mu1_Eout0.Eout_Prob = ( 1.0 - alpha ) * low_cum_prob.first.y +
    alpha * high_cum_prob.first.y;

  Ein_data->mu1_Eout1.UB_Eout = 1.0;
  Ein_data->mu1_Eout1.Eout_Prob = ( 1.0 - alpha ) * low_cum_prob.second.y +
    alpha * high_cum_prob.second.y;

  Ein_data->mu1_ubase_map.interpolate( alpha,
    low_cum_prob.ubase_map, high_cum_prob.ubase_map );

  Ein_data->mu1_Eout0.phys_Eout = Ein_data->mu1_ubase_map.un_unit_base( 0.0 );
  Ein_data->mu1_Eout1.phys_Eout = Ein_data->mu1_ubase_map.un_unit_base( 1.0 );

  return svar;
}

// ********* class Jdist::one_mu *********
// ----------- Jdist::one_mu::copy --------------
// Make a copy for the ENDL_kludge routine
void Jdist::one_mu::copy( const Jdist::one_mu &to_copy )
{
  // copy the basic information
  set_mu( to_copy.get_mu( ) );
  interp_type = to_copy.interp_type;
  ubase_map.copy( to_copy.ubase_map );

  // copy the entries
  for( Jdist::one_mu::const_iterator this_entry = to_copy.begin( );
         this_entry != to_copy.end( ); ++this_entry )
  {
    add_entry( this_entry->x, this_entry->y );
  }
}
// ----------- Jdist::one_mu::form_cum_prob --------------
// Forms the list of cumulative probabilities
void Jdist::one_mu::form_cum_prob( )
{
  // copy the data
  cum_prob.Eout_interp = interp_type;
  for( Jdist::one_mu::const_iterator Eout_ptr = begin( );
       Eout_ptr != end( ); ++Eout_ptr )
  {
    Cum::cumulative_prob_list::iterator cum_prob_ptr = cum_prob.insert(
      cum_prob.end( ), Cum::cumulative_prob_entry( ) );
    cum_prob_ptr->E_out = Eout_ptr->x;
    cum_prob_ptr->Prob = Eout_ptr->y;
  }
  // now form the slopes and cumulative probabilities
  if( interp_type == Terp::HISTOGRAM )
  {
    cum_prob.get_cum_prob_flat( );
  }
  else // lin-lin
  {
    cum_prob.get_cum_prob_linlin( );
  }
}
// ----------- Jdist::one_mu::set_E_mu_P_data --------------
// Sets up Jdata::E_mu_P_data at unit-base outgoing energy UB_Eout
void Jdist::one_mu::set_E_mu_P_data( double UB_Eout, double phys_mu,
    Jdist::one_mu::const_iterator prev_data,
    Jdist::one_mu::const_iterator next_data, Jdata::E_mu_P_data *mid_data )
{
  mid_data->UB_mu = get_mu( );
  mid_data->phys_mu = phys_mu;
  mid_data->mu_Prob = mu_Prob;
  mid_data->UB_Eout = UB_Eout;
  mid_data->phys_Eout = ubase_map.un_unit_base( UB_Eout );
  static double etol = Global.Value( "tight_tol" );
  bool is_OK;
  
  if( interp_type == Terp::HISTOGRAM )
  {
    if( UB_Eout > ( 1 - etol )*next_data->x )
    {
      mid_data->Eout_Prob = next_data->y;
    }
    else
    {
      mid_data->Eout_Prob = prev_data->y;
    }
  }
  else
  {
    mid_data->Eout_Prob = prev_data->linlin_interp( UB_Eout, *next_data, &is_OK );
  }

}

// ********* class Jdist::one_joint_dist *********
// ----------- Jdist::one_joint_dist::copy --------------
void Jdist::one_joint_dist::copy( const Jdist::one_joint_dist &to_copy )
// make a copy at the threshold, used by the ENDL_kludge routine
{
  for( Jdist::one_joint_dist::const_iterator copy_ptr = to_copy.begin( );
       copy_ptr != to_copy.end( ); ++copy_ptr )
  {
    // make a new energy distribution
    Jdist::one_joint_dist::iterator new_e_dist_ptr = insert( end( ), Jdist::one_mu( ) );
    // copy into it
    new_e_dist_ptr->copy( *copy_ptr );
  }
}
// ----------- Jdist::one_joint_dist::set_mu_Prob --------------
// Sets mu_Prob for ENDL double-differential data
void Jdist::one_joint_dist::set_mu_Prob( Adist::angle_dist::iterator &angles )
{
  static double abs_tol = Global.Value( "tight_tol" );
  
  Ddvec::dd_vector::const_iterator prev_mu = angles->begin( );
  Ddvec::dd_vector::const_iterator next_mu = prev_mu;
  ++next_mu;

  Jdist::one_joint_dist::iterator joint_ptr = begin( );
  if( std::abs( prev_mu->x - joint_ptr->get_mu( ) ) >
      abs_tol * std::abs( prev_mu->x + joint_ptr->get_mu( ) ) )
  {
    Msg::FatalError( "Jdist::one_joint_dist::set_mu_Prob",
		"starting mu values different" );
  }
  
  for( ; joint_ptr != end( );
       ++joint_ptr )
  {
    // increment the angular data
    while( ( next_mu != angles->end( ) ) &&
           ( next_mu->x <= joint_ptr->get_mu( ) +
	       abs_tol * std::abs( joint_ptr->get_mu( ) ) ) )
    {
      prev_mu = next_mu;
      ++next_mu;
    }

    if( joint_ptr->get_mu( ) <= prev_mu->x + abs_tol *
	std::abs( prev_mu->x ) )
    {
      joint_ptr->mu_Prob = prev_mu->y;
    }
    else
    {
      double alpha = ( joint_ptr->get_mu( ) - prev_mu->x )/
	( next_mu->x - prev_mu->x );
      joint_ptr->mu_Prob = ( 1.0 - alpha ) * prev_mu->y +
	alpha * next_mu->y;
    }
  }
}
// ----------- Jdist::one_joint_dist::convert_ENDF --------------
// Converts ENDF double-differential data to our format
void Jdist::one_joint_dist::convert_ENDF( )
{
  Adist::angle_dist angleData;
  Adist::angle_dist::iterator angle_ptr = angleData.insert( angleData.begin( ),
						     Ddvec::dd_vector( ) );
  angle_ptr->interp_type = mu_interp.flag;
  angle_ptr->set_E_in( get_E_in( ) );
  
  // get the angular probabilities
  for( Jdist::one_joint_dist::iterator joint_ptr = begin( );
       joint_ptr != end( ); ++joint_ptr )
  {
    double Prob = joint_ptr->get_norm( );
    angle_ptr->add_entry( joint_ptr->get_mu( ), Prob );

    // normalize the energy distribution
    if( Prob == 0.0 )
    {
      Prob = 1.0;
    }
    *joint_ptr *= 1.0 / Prob;
  }

  angle_ptr->renorm( false );

  if( Ein_interp.qualifier == Terp::UNITBASE )
  {
    mu_to_unit_base( angle_ptr );
  }
  
  Adist::angle_dist::iterator these_angles = angleData.begin( );
  convert_ENDL( these_angles );

  // We no longer need angleData
  angleData.erase( angleData.begin( ), angleData.end( ) );
}
// ----------- Jdist::one_joint_dist::convert_ENDL --------------
// Converts ENDL double-differential data to our format
void Jdist::one_joint_dist::convert_ENDL( Adist::angle_dist::iterator &angles )
{
  // first check the incident energy
  if( get_E_in( ) != angles->get_E_in( ) )
  {
    Msg::FatalError("Jdist::one_joint_dist::convert_ENDL",
      Msg::pastenum("E_in for angles: ", angles->get_E_in( ) ) +
      Msg::pastenum(" different from E_in for joint_dist: ", get_E_in( ) ) );
  }

  // Set the mu_Prob values
  set_mu_Prob( angles );
  
  Jdist::one_joint_dist::iterator joint_mu = begin( );
  Jdist::one_joint_dist::iterator prev_joint_mu = begin( );
  Ddvec::dd_vector::const_iterator angle_mu = angles->begin( );

  // First, check to see whether we need to interpolate the table to intermediate mu
  for( ; ( joint_mu != end( ) ) && ( angle_mu != angles->end( ) ); ++angle_mu )
  {
    if( ( joint_mu->get_mu( ) < angle_mu->x ) ||
        ( ( joint_mu->get_mu( ) > angle_mu->x ) && ( joint_mu == prev_joint_mu ) ) )
    {
      Msg::FatalError("Jdist::one_joint_dist::convert_ENDL ",
        Msg::pastenum(" mu for angles: ", angle_mu->x ) +
        Msg::pastenum(" different from for joint_mu: ", joint_mu->get_mu( ) ) );
    }
    else if( joint_mu->get_mu( ) > angle_mu->x )
    {
      // insert energy probability densities for an intermediate mu
      // ignore the returned bool
      if( mu_interp.qualifier == Terp::UNITBASE )
      {
        interpolate_mu_UB( angle_mu, prev_joint_mu, joint_mu );
        ++prev_joint_mu;  // points to the inserted intermediate list
      }
      else if( mu_interp.qualifier == Terp::CUMULATIVE_POINTS )
      {
        interpolate_mu_CP( angle_mu, prev_joint_mu, joint_mu );
        ++prev_joint_mu;  // points to the inserted intermediate list
      }
    }
    else
    {
      // they match
      prev_joint_mu = joint_mu;  // go to the next mu
      ++joint_mu;
    }
  }
}
// ----------- Jdist::one_joint_dist::interpolate_mu_UB --------------
// insert energy probability densities for an intermediate mu
// In this version the angular distribution and energy distribution
// are interpolated separately and then multiplied.
// Use unit-base interpolation.
void Jdist::one_joint_dist::interpolate_mu_UB( Ddvec::dd_vector::const_iterator &angle_ptr,
  Jdist::one_joint_dist::iterator &prev_joint_mu,
  Jdist::one_joint_dist::iterator &next_joint_mu )
{
  Jdist::one_joint_dist::iterator new_mu_ptr = insert( next_joint_mu, Jdist::one_mu( ) );
  new_mu_ptr->mu_Prob = angle_ptr->y;
  new_mu_ptr->interpolate( angle_ptr->x, *prev_joint_mu, *next_joint_mu );

  // Set up the unit-base map.
  double mu_diff = next_joint_mu->get_mu( ) - prev_joint_mu->get_mu( );
  if( mu_diff <= 0.0 )
  {
    Msg::FatalError( "Jdist::one_joint_dist::interpolate_mu_UB",
		"mu values out of order" );
  }
  double alpha = ( angle_ptr->x - prev_joint_mu->get_mu( ) )/mu_diff;
  new_mu_ptr->ubase_map.interpolate( alpha, prev_joint_mu->ubase_map,
				       next_joint_mu->ubase_map );
}
// ----------- Jdist::one_joint_dist::interpolate_mu_CP --------------
// insert energy probability densities for an intermediate mu
// In this version the angular distribution and energy distribution
// are interpolated separately and then multiplied.
// Use cumulative-points interpolation.
void Jdist::one_joint_dist::interpolate_mu_CP( Ddvec::dd_vector::const_iterator &angle_ptr,
  Jdist::one_joint_dist::iterator &prev_joint_mu,
  Jdist::one_joint_dist::iterator &next_joint_mu )
{
  Jdist::one_joint_dist::iterator new_mu_ptr = insert( next_joint_mu, Jdist::one_mu( ) );
  new_mu_ptr->set_mu( angle_ptr->x );
  new_mu_ptr->mu_Prob = angle_ptr->y;

  // Set up the interpolation weight
  double mu_diff = next_joint_mu->get_mu( ) - prev_joint_mu->get_mu( );
  if( mu_diff <= 0.0 )
  {
    Msg::FatalError( "Jdist::one_joint_dist::interpolate_mu_CP",
		"mu values out of order" );
  }
  double alpha = ( angle_ptr->x - prev_joint_mu->get_mu( ) )/mu_diff;

  Cum::cumulative_prob_list::const_iterator prev_CP_low = prev_joint_mu->cum_prob.begin( );
  Cum::cumulative_prob_list::const_iterator next_CP_low = next_joint_mu->cum_prob.begin( );

  Cum::cumulative_prob_list::const_iterator prev_CP_high = prev_CP_low;
  ++prev_CP_high;
  Cum::cumulative_prob_list::const_iterator next_CP_high = next_CP_low;
  ++next_CP_high;

  static double abs_tol = Global.Value( "tight_tol" );

  double low_A = 0.0;
  
  while( ( prev_CP_high != prev_joint_mu->cum_prob.end( ) ) &&
	 ( next_CP_high != next_joint_mu->cum_prob.end( ) ) )
  {
    // append a new cumulative_prob_entry
    Cum::cumulative_prob_list::iterator new_CP_ptr =
      new_mu_ptr->cum_prob.insert( new_mu_ptr->cum_prob.end( ),
				   Cum::cumulative_prob_entry( ) );

    // interpolate
    new_CP_ptr->CP_interpolate( alpha, low_A, *prev_CP_low, *next_CP_low );

    // the next cumulative probability
    double high_A = ( prev_CP_high->cum_prob <= next_CP_high->cum_prob ) ?
      prev_CP_high->cum_prob : next_CP_high->cum_prob;

    // special for the last link
    if( high_A > 1.0 - abs_tol )
    {
      new_CP_ptr =
        new_mu_ptr->cum_prob.insert( new_mu_ptr->cum_prob.end( ),
                                   Cum::cumulative_prob_entry( ) );
      new_CP_ptr->CP_interpolate( alpha, 1.0, *prev_CP_high, *next_CP_high );
    }

    if( prev_CP_high->cum_prob < ( 1.0 + abs_tol ) * high_A )
    {
      prev_CP_low = prev_CP_high;
      ++prev_CP_high;
    }
    if( next_CP_high->cum_prob < ( 1.0 + abs_tol ) * high_A )
    {
      next_CP_low = next_CP_high;
      ++next_CP_high;
    }

    low_A = high_A;
  }

  // set the slopes
  new_mu_ptr->cum_prob.set_slopes( );
  
  // Set up the unit-base map.
  new_mu_ptr->ubase_map.interpolate( alpha, prev_joint_mu->ubase_map,
				     next_joint_mu->ubase_map );
}
// ----------- Jdist::one_joint_dist::mu_to_unit_base --------------
// Maps the direction cosines to 0 <= mu <= 1
void Jdist::one_joint_dist::mu_to_unit_base( Adist::angle_dist::iterator &angles )
{
  // map the angular data to unit base; there is already a consistency check.
  Ddvec::unit_base_map ubase_Map;
  angles->unit_base( true, &ubase_Map );
  
  Jdist::one_joint_dist::iterator this_mu = begin( );
  Jdist::one_joint_dist::iterator last_mu = end( );
  --last_mu;
  mu_ubase_map.Eout_min = this_mu->get_mu( );
  mu_ubase_map.Eout_max = last_mu->get_mu( );
  double scale_mu = mu_ubase_map.Eout_max - mu_ubase_map.Eout_min;

  for( ; this_mu != end( ); ++this_mu )
  {
    bool is_OK;
    double current_mu = mu_ubase_map.to_unit_base( this_mu->get_mu( ),
						   &is_OK );
    if( !is_OK )
    {
      Msg::FatalError( "Jdist::one_joint_dist::mu_to_unit_base",
		  "bad interpolation" );
    }
    this_mu->set_mu( current_mu );
    this_mu->mu_Prob *= scale_mu;
  }
}
// ----------- Jdist::one_joint_dist::Eout_to_unit_base --------------
// Maps the outgoing energies to 0 <= Eout <= 1
void Jdist::one_joint_dist::Eout_to_unit_base( )
{
  for( Jdist::one_joint_dist::iterator this_mu = begin( );
       this_mu !=end( ); ++this_mu )
  {
    this_mu->unit_base( false, &this_mu->ubase_map );
  }
}
// ----------- Jdist::one_joint_dist::form_Eout_cum_prob --------------
// Sets up cumulative probabilities for outgoing energy
void Jdist::one_joint_dist::form_Eout_cum_prob( )
{
  Jdist::one_joint_dist::iterator this_mu = begin( );
  for( ; this_mu != end( ); ++this_mu )
  {
    this_mu->form_cum_prob( );
  }
}
// ----------- Jdist::one_joint_dist::set_mu_data --------------
// Sets the mu data for cumulative-points interpolation 
void Jdist::one_joint_dist::set_mu_data( double this_mu,
				  Jdist::one_joint_dist::const_iterator low_mu,
				  Jdist::one_joint_dist::const_iterator high_mu,
		      Jdata::current_data *Ein_data, int mu_level )
{
  // for interpolation in mu
  double denom = high_mu->get_mu( ) - low_mu->get_mu( );
  if( denom <= 0.0 )
  {
    Msg::FatalError( "Jdist::one_joint_dist::set_mu_data",
		"mu values out of order" );
  }
  
  double alpha = ( this_mu - low_mu->get_mu( ) ) / denom;
  double mu_Prob = ( 1.0 - alpha ) * low_mu->mu_Prob +
    alpha * high_mu->mu_Prob;

  double phys_mu = mu_ubase_map.un_unit_base( this_mu );
  
  if( mu_level == 0 )  // mu_0
  {
    Ein_data->mu0_Eout0.UB_mu = this_mu;
    Ein_data->mu0_Eout0.mu_Prob = mu_Prob;
    Ein_data->mu0_Eout0.phys_mu = phys_mu;
    Ein_data->mu0_Eout1.UB_mu = this_mu;
    Ein_data->mu0_Eout1.mu_Prob = mu_Prob;
    Ein_data->mu0_Eout1.phys_mu = phys_mu;
  }
  else
  {
    Ein_data->mu1_Eout0.UB_mu = this_mu;
    Ein_data->mu1_Eout0.mu_Prob = mu_Prob;
    Ein_data->mu1_Eout0.phys_mu = phys_mu;
    Ein_data->mu1_Eout1.UB_mu = this_mu;
    Ein_data->mu1_Eout1.mu_Prob = mu_Prob;
    Ein_data->mu1_Eout1.phys_mu = phys_mu;
  }
}

// ********* class Jdist::joint_dist *********
// ----------- Jdist::joint_dist::read_data --------------
// Reads double-differential energy-angle data
void Jdist::joint_dist::read_data( Dpar::data_parser &inFile, int num_Ein )
{
  interp_flag_F::read_3d_interpolation_ENDL( inFile, &Ein_interp, &mu_interp,
                       &Eout_interp );

  Jdist::joint_dist::iterator new_Ein_ptr;
  Jdist::one_joint_dist::iterator new_mu_ptr;

  double E_out;
  double Prob;

  // read the data
  for( int Ein_count = 0; Ein_count < num_Ein; ++Ein_count )
  {
    // make a new link for this incident energy
    new_Ein_ptr = insert( end( ), Jdist::one_joint_dist( ) );
    new_Ein_ptr->set_E_in( inFile.get_next_double( ) );
    new_Ein_ptr->mu_interp = mu_interp;
    new_Ein_ptr->Ein_interp = Ein_interp;
    new_Ein_ptr->Eout_interp = Eout_interp;
    new_Ein_ptr->ENDL_data = ENDL_data;

    // loop over cosine
    int num_mu = inFile.get_next_int( );
    for( int mu_count = 0; mu_count < num_mu; ++mu_count )
    {
      // make a new energy distribution
      new_mu_ptr = new_Ein_ptr->insert( new_Ein_ptr->end( ), Jdist::one_mu( ) );
      new_mu_ptr->set_mu( inFile.get_next_double( ) );
      new_mu_ptr->interp_type = Eout_interp;
      // read the (energy, probability density) pairs
      int num_Eout = inFile.get_next_int( );
      for( int Eout_count = 0; Eout_count < num_Eout; ++Eout_count )
      {
        E_out = inFile.get_next_double( );
        Prob = inFile.get_next_double( );
        new_mu_ptr->add_entry( E_out, Prob );
      }
    }
  }
  
  if( ENDL_data )
  {
    // we may need to copy the first distribution at the threshold
    ENDL_kludge( );
  }
  
  if( mu_interp.qualifier == Terp::UNITBASE )
  {
    Eout_to_unit_base( );
  }
  else if( mu_interp.qualifier == Terp::CUMULATIVE_POINTS )
  {
    form_Eout_cum_prob( );
  }

  if( ENDL_data && ( Ein_interp.qualifier == Terp::UNITBASE ) )
  {
    mu_to_unit_base( );
  }

  if( ENDL_data )
  {
    // convert to our format
    convert_ENDL( );
  }
  else  // ENDF data format
  {
    // convert to our format
    convert_ENDF( );
  }

}
// ----------- Jdist::joint_dist::ENDL_kludge --------------
void Jdist::joint_dist::ENDL_kludge( )
// ENDL has the convention that the distribution is assumed
// independent of energy at low incident energies.
{
  Jdist::joint_dist::iterator first_dist = begin( );
  if( first_dist->get_E_in( ) > angle_data.threshold )
  {
    // prepend a copy at the threshold energy
    Jdist::joint_dist::iterator second_dist = begin( );
    first_dist = insert( begin( ), Jdist::one_joint_dist( ) );
    first_dist->copy( *second_dist );
    first_dist->set_E_in( angle_data.threshold );
  }
}
// ----------- Jdist::joint_dist::convert_ENDL --------------
// Converts ENDL double-differential data to our format
void Jdist::joint_dist::convert_ENDL( )
{
  // normalize the angular data
  Jdist::joint_dist::iterator joint_mu = begin( );
  Adist::angle_dist::iterator angle_mu = angle_data.begin( );

  for( ; ( joint_mu != end( ) ) && ( angle_mu != angle_data.end( ) );
       ++joint_mu, ++angle_mu )
  {
    angle_mu->renorm( false );
  
    joint_mu->convert_ENDL( angle_mu );
  }
}
// ----------- Jdist::joint_dist::convert_ENDF --------------
// Converts ENDF double-differential data to our format
void Jdist::joint_dist::convert_ENDF( )
{
  // There should be no angular data.
  if( !angle_data.empty( ) )
  {
    angle_data.erase( angle_data.begin( ), angle_data.end( ) );
  }

  Jdist::joint_dist::iterator joint_mu = begin( );

  for( ; joint_mu != end( ); ++joint_mu )
  {
    joint_mu->convert_ENDF( );
  }
}
// ----------- Jdist::joint_dist::mu_to_unit_base --------------
// Maps the direction cosines to 0 <= mu <= 1
void Jdist::joint_dist::mu_to_unit_base( )
{
  Adist::angle_dist::iterator angle_mu = angle_data.begin( );
  Jdist::joint_dist::iterator joint_mu = begin( );

  for( ; ( joint_mu != end( ) ) && ( angle_mu != angle_data.end( ) );
       ++joint_mu, ++angle_mu )
  {
    joint_mu->mu_to_unit_base( angle_mu );
  }
}
// ----------- Jdist::joint_dist::Eout_to_unit_base --------------
// Maps the outgoing energies to 0 <= Eout <= 1
void Jdist::joint_dist::Eout_to_unit_base( )
{
  Jdist::joint_dist::iterator joint_mu = begin( );
  for( ; joint_mu != end( ); ++joint_mu )
  {
    joint_mu->Eout_to_unit_base( );
  }
}
// ----------- Jdist::joint_dist::form_Eout_cum_prob --------------
// Sets up cumulative probabilities for outgoing energy
void Jdist::joint_dist::form_Eout_cum_prob( )
{
  Jdist::joint_dist::iterator joint_mu = begin( );
  for( ; joint_mu != end( ); ++joint_mu )
  {
    joint_mu->form_Eout_cum_prob( );
  }
}
// ----------- Jdist::joint_dist::get_T --------------
// Calculates the transfer matrix for this particle.
void Jdist::joint_dist::get_T( const Ddvec::dd_vector& sigma, const Ddvec::dd_vector& multiple,
  const Ddvec::dd_vector& weight, Trf::T_matrix& transfer )
{
  if( ( Ein_interp.qualifier != Terp::UNITBASE ) || ( Ein_interp.flag != Terp::LINLIN ) )
  {
    Msg::FatalError( "Jdist::joint_dist::get_T",
		     "Ein_interp not implemented" );
  }

  bool interp_OK = ( ( mu_interp.qualifier == Terp::UNITBASE ) ||
		     ( mu_interp.qualifier == Terp::CUMULATIVE_POINTS ) ) &&
    ( mu_interp.flag == Terp::LINLIN );
  if( !interp_OK )
  {
    Msg::FatalError( "Jdist::joint_dist::get_T",
		     "mu_interp not implemented" );
  }
  
  if( ( Eout_interp != Terp::LINLIN ) && ( Eout_interp != Terp::HISTOGRAM ) )
  {
    Msg::FatalError( "Jdist::joint_dist::get_T",
		     "Eout_interp not implemented" );
  }
  if( ENDL_data && ( angle_data.size( ) == 0 ) )
  {
    Msg::FatalError( "Jdist::joint_dist::get_T",
      "ENDL table missing angular data" );
  }
  if( ( angle_data.size( ) > 0 ) &&
      ( Ein_interp.qualifier != angle_data.Ein_interp.qualifier ) )
  {
    Msg::FatalError( "Jdist::joint_dist::get_T",
      "Incident energy interpolation not consistent with the angular data" );
  }

  bool done = get_Ein_range( sigma, multiple, weight, transfer.e_flux,
    transfer.in_groups );
  if( done )
  {
    transfer.zero_transfer( );
  }
  transfer.threshold = sigma.begin( )->x;

  long int quad_count = 0;  // number of 2-d quadratures
  long int Ein_F_count= 0;  // number of calls to joint_dist_F::Ein_F
  long int Eout_F_count= 0;  // number of calls to joint_dist_F::Eout_F
  long int mu_F_count = 0;  // number of calls to joint_dist_F::mu_F

  // now do the integrals bin by bin
#pragma omp parallel for schedule( dynamic, 1 ) default( none )	\
  shared( sigma, multiple, weight, transfer ) \
  reduction( +: quad_count ) reduction( +: Ein_F_count ) \
  reduction( +: Eout_F_count ) reduction( +: mu_F_count )
  for( int Ein_bin = first_Ein; Ein_bin < last_Ein; ++Ein_bin )
  {
    Jdist::joint_dist_param Ein_param;
    // set up the data range for this bin
    Ein_param.setup_bin( Ein_bin, sigma, multiple, weight, transfer.e_flux,
                         transfer.in_groups );
    Ein_param.Eout_quad_rule = transfer.Eout_quad_rule;
    Ein_param.mu_quad_rule = transfer.mu_quad_rule;
    setup_param( &Ein_param );
    for( ; ; )
    {
      set_Ein_range( &Ein_param );   // get the incident energy interval
      mu_data_ladder( transfer, &Ein_param );  // loop over the direction cosines
      bool Done = next_Ein_pair( Ein_param.data_E_1, &Ein_param );   // go to the next interval
      if( Done )
      {
        break;
      }
    }
    quad_count += Ein_param.quad_count;
    Ein_F_count += Ein_param.Ein_F_count;
    Eout_F_count += Ein_param.Eout_F_count;
    mu_F_count += Ein_param.mu_F_count;
  } // end of parallel loop

  // print the counts of function evaluations
  std::cout << "3d quadratures: " << quad_count << std::endl;
  std::cout << "joint_dist_F::Ein_F calls: " << Ein_F_count << std::endl;
  std::cout << "joint_dist_F::mu_F calls: " << mu_F_count << std::endl;
  std::cout << "joint_dist_F::Eout_F calls: " << Eout_F_count << std::endl;
  std::cout << "average joint_dist_F::Ein_F_count: " << 1.0*Ein_F_count/quad_count << std::endl;
  std::cout << "average joint_dist_F::mu_F_count: " << 1.0*mu_F_count/Ein_F_count << std::endl;
  std::cout << "average joint_dist_F::Eout_F_count: " << 1.0*Eout_F_count/mu_F_count << std::endl;
}
// ----------- Jdist::joint_dist::setup_param ------------------
// Initializes the quadrature parameters
void Jdist::joint_dist::setup_param( Jdist::joint_dist_param *Ein_param )
{
  static double skip_tol = Global.Value( "tight_tol" );

  Ein_param->mu_params.Eout_params.Eout_interp = Eout_interp;
  
  Ein_param->this_Ein_ptr = begin( );
  Ein_param->next_Ein_ptr = Ein_param->this_Ein_ptr;
  ++Ein_param->next_Ein_ptr;

  while( Ein_param->next_Ein_ptr->get_E_in( ) < Ein_param->data_E_0 *
	 ( 1.0 + skip_tol ) )
  {
    Ein_param->this_Ein_ptr = Ein_param->next_Ein_ptr;
    ++Ein_param->next_Ein_ptr;
  }
  double first_Ein = Ein_param->this_Ein_ptr->get_E_in( );
  if( first_Ein > Ein_param->data_E_0 )
  {
    Ein_param->data_E_0 = first_Ein;
    bool data_bad = Ein_param->update_pointers( first_Ein );
    if( data_bad )
    {
      Msg::FatalError( "Jdist::joint_dist::setup_param",
		       "energies inconsistent" );
    }
  }
}
// ----------- Jdist::joint_dist::get_Ein_range --------------
//  Gets the range of nontrivial incident energy bins; computes first_Ein and last_Ein
// returns true if the threshold is too high for the energy bins
bool Jdist::joint_dist::get_Ein_range( const Ddvec::dd_vector& sigma, const Ddvec::dd_vector& mult,
    const Ddvec::dd_vector& weight,
    const Lgdata::Flux_List& e_flux, const Egp::Energy_groups& Ein_groups )
{
  double E_first;
  double E_last;
  Jdist::joint_dist_param initial_param;
  bool done = initial_param.get_Ein_range( sigma, mult, weight, e_flux,
                                         Ein_groups, &E_first, &E_last );
  if( done ) return true;

  // check the range of incident energies for the probability data
  Jdist::joint_dist::const_iterator Ein_data_ptr = begin( );
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
// -----------  Jdist::joint_dist::set_Ein_range ------------------
// Sets the range of incident energies for this intergration
void Jdist::joint_dist::set_Ein_range( Jdist::joint_dist_param *Ein_param )
{
  Ein_param->set_Ein_range( );
  double this_E = Ein_param->this_Ein_ptr->get_E_in( );
  if( this_E > Ein_param->data_E_0 ) Ein_param->data_E_0 = this_E;
  this_E = Ein_param->next_Ein_ptr->get_E_in( );
  if( this_E < Ein_param->data_E_1 ) Ein_param->data_E_1 = this_E;

  if( Ein_param->data_E_1 < Ein_param->data_E_0 )
  {
    Msg::FatalError( "Jdist::joint_dist::set_Ein_range",
		     "check the incident energies" );
  }
}
// ----------- Jdist::joint_dist::one_box ------------------
// Integrate over one E-E' box
void Jdist::joint_dist::one_box( Trf::T_matrix& transfer, int Eout_count,
   Jdist::joint_dist_param *Ein_param )
{
  // pointers to the hit_lists
  Jdist::joint_dist_hits::iterator mu0_low_hit_ptr = Ein_param->lower_mu0_hits.begin( );
  Jdist::joint_dist_hits::iterator next_mu0_low_ptr = mu0_low_hit_ptr;
  ++next_mu0_low_ptr;
  Jdist::joint_dist_hits::iterator mu1_low_hit_ptr = Ein_param->lower_mu1_hits.begin( );
  Jdist::joint_dist_hits::iterator next_mu1_low_ptr = mu1_low_hit_ptr;
  ++next_mu1_low_ptr;
  Jdist::joint_dist_hits::iterator mu0_high_hit_ptr = Ein_param->upper_mu0_hits.begin( );
  Jdist::joint_dist_hits::iterator next_mu0_high_ptr = mu0_high_hit_ptr;
  ++next_mu0_high_ptr;
  Jdist::joint_dist_hits::iterator mu1_high_hit_ptr = Ein_param->upper_mu1_hits.begin( );
  Jdist::joint_dist_hits::iterator next_mu1_high_ptr = mu1_high_hit_ptr;
  ++next_mu1_high_ptr;

  // the range of incident energies
  double lower_Ein = mu0_low_hit_ptr->E_in;
  double upper_Ein = next_mu0_low_ptr->E_in;
  double this_Ein = next_mu1_low_ptr->E_in;
  if( this_Ein < upper_Ein )
  {
    upper_Ein = this_Ein;
  }
  this_Ein = next_mu0_high_ptr->E_in;
  if( this_Ein < upper_Ein )
  {
    upper_Ein = this_Ein;
  }
  this_Ein = next_mu1_high_ptr->E_in;
  if( this_Ein < upper_Ein )
  {
    upper_Ein = this_Ein;
  }

  for( ; ; )
  {
    if( ( ( mu0_low_hit_ptr->hit_edge == Box::ABOVE ) ||
          ( mu0_low_hit_ptr->hit_edge == Box::TOP_OUT ) ) &&
	( ( mu1_low_hit_ptr->hit_edge == Box::ABOVE ) ||
          ( mu1_low_hit_ptr->hit_edge == Box::TOP_OUT ) ) )
    {
      // do nothing---we are above the E-E' box
    }
    else if( ( ( mu0_high_hit_ptr->hit_edge == Box::BOTTOM_OUT ) ||
               (  mu0_high_hit_ptr->hit_edge == Box::BELOW ) ) &&
	     ( ( mu1_high_hit_ptr->hit_edge == Box::BOTTOM_OUT ) ||
               (  mu1_high_hit_ptr->hit_edge == Box::BELOW ) ) )
    {
      // do nothing---we are below the E-E' box
    }
    else
    {
      update_T( transfer, Eout_count, lower_Ein, upper_Ein, Ein_param );
    }
    // update the energy range
    lower_Ein = upper_Ein;
    static double e_tol = Global.Value( "looser_tol" );
    if( next_mu0_low_ptr->E_in <= lower_Ein*(1 + e_tol ) )
    {
      mu0_low_hit_ptr = next_mu0_low_ptr;
      ++next_mu0_low_ptr;
      if( next_mu0_low_ptr == Ein_param->lower_mu0_hits.end( ) )
      {
        break;
      }
    }
    upper_Ein =  next_mu0_low_ptr->E_in;

    if( next_mu1_low_ptr->E_in <= lower_Ein*(1 + e_tol ) )
    {
      mu1_low_hit_ptr = next_mu1_low_ptr;
      ++next_mu1_low_ptr;
      if( next_mu1_low_ptr == Ein_param->lower_mu1_hits.end( ) )
      {
        break;
      }
    }
    this_Ein = next_mu1_low_ptr->E_in;
    if( upper_Ein > this_Ein )
    {
      upper_Ein = this_Ein;
    }

    if( next_mu0_high_ptr->E_in <= lower_Ein*(1 + e_tol ) )
    {
      mu0_high_hit_ptr = next_mu0_high_ptr;
      ++next_mu0_high_ptr;
      if( next_mu0_high_ptr == Ein_param->upper_mu0_hits.end( ) )
      {
        break;
      }
    }
    this_Ein = next_mu0_high_ptr->E_in;
    if( upper_Ein > this_Ein )
    {
      upper_Ein = this_Ein;
    }

    if( next_mu1_high_ptr->E_in <= lower_Ein*(1 + e_tol ) )
    {
      mu1_high_hit_ptr = next_mu1_high_ptr;
      ++next_mu1_high_ptr;
      if( next_mu1_high_ptr == Ein_param->upper_mu1_hits.end( ) )
      {
        break;
      }
    }
    this_Ein = next_mu1_high_ptr->E_in;
    if( upper_Ein > this_Ein )
    {
      upper_Ein = this_Ein;
    }
  }
}
// ----------- Jdist::joint_dist::update_T --------------
// Increments the transfer matrix
void Jdist::joint_dist::update_T( Trf::T_matrix &transfer, int Eout_count,
       double Ein0_orig, double Ein1_orig, Jdist::joint_dist_param *Ein_param )
{
  static double tol = Global.Value( "quad_tol" );
  // differences of nearly-equal numbers can cause problems; when to skip an interval
  static double skip_tol = Global.Value( "tight_tol" );
    
  // a vector to store the integrals, one Legendre order
  Coef::coef_vector value( transfer.order, transfer.conserve );
  value.set_zero( );
  // parameters for the integration
  Qparam::QuadParamBase *params = static_cast< Qparam::QuadParamBase* >( Ein_param );

  // initial values for the lower and upper limits of integration
  double use_Ein0 = Ein0_orig;
  double use_Ein1 = Ein1_orig;

  // loop over the cross section data
  Ein_param->this_sigma = Ein_param->first_ladder_sigma;
  Ein_param->next_sigma = Ein_param->this_sigma;
  ++Ein_param->next_sigma;
  // use_Ein0 may be past Ein_param->next_sigma
  while( ( Ein_param->this_sigma != Ein_param->last_ladder_sigma ) &&
         ( Ein_param->next_sigma->x < use_Ein0 ) )
  {
    Ein_param->this_sigma = Ein_param->next_sigma;
    ++Ein_param->next_sigma;
  }
  for( ; ( Ein_param->this_sigma != Ein_param->last_ladder_sigma ) &&
         ( Ein_param->this_sigma->x < Ein1_orig );
       Ein_param->this_sigma = Ein_param->next_sigma, ++Ein_param->next_sigma )
  {
    use_Ein0 = ( Ein_param->this_sigma->x < Ein0_orig ) ? Ein0_orig :
      Ein_param->this_sigma->x;
    use_Ein1 = ( Ein_param->next_sigma->x > Ein1_orig ) ? Ein1_orig :
      Ein_param->next_sigma->x;
    // differences of nearly-equal numbers can cause problems
    if( use_Ein1 - use_Ein0 <= use_Ein1 * skip_tol )
    {
      continue;  // skip this interval
    }
    else
    {
      quad_F::integrate( joint_dist_F::Ein_F, transfer.Ein_quad_rule,
                         use_Ein0, use_Ein1, params, tol, &value );
    }
    // add this integral
    transfer( Ein_param->Ein_count, Eout_count ) += value;
    // increment the function counts
    Ein_param->Ein_F_count += Ein_param->func_count;
    ++Ein_param->quad_count;
  }
}
// ----------- Jdist::joint_dist::next_Ein_pair --------------
// Go to the next pair of incident energies.  Returns "true" when finished.
bool Jdist::joint_dist::next_Ein_pair( double E_in, Jdist::joint_dist_param *Ein_param )
{
  bool done = Ein_param->update_bin_pointers( E_in );
  if( !done )
  {
    static double etol = Global.Value( "tight_tol" );
    double E_tol = E_in * etol;
    if( E_in + E_tol >= Ein_param->next_Ein_ptr->get_E_in( ) )
    {
      while( E_in + E_tol >= Ein_param->next_Ein_ptr->get_E_in( ) )
      {
        // get the next energy data
        Ein_param->this_Ein_ptr = Ein_param->next_Ein_ptr;
        ++Ein_param->next_Ein_ptr;
        if( Ein_param->next_Ein_ptr == end ( ) )
        {
          return true;
        }
      }
    }
  }
  return done;
}
// ----------- Jdist::joint_dist::mu_data_ladder --------------
// Handles the ( cosine, Eout, probability ) data for one pair of incident energies
void Jdist::joint_dist::mu_data_ladder( Trf::T_matrix& transfer,
					Jdist::joint_dist_param *Ein_param )
{
  start_mu_data( Ein_param );
  for( ; ; )
  {
    Eout_data_ladder( transfer, Ein_param );
    bool done = next_mu_pairs( Ein_param );   // go to the next mu interval
    if( done )
    {
      break;
    }
  }
}
// ----------- Jdist::joint_dist::start_mu_data------------
// Initializes the pointers to the ( Eout, probability ) data for current cosines
void Jdist::joint_dist::start_mu_data( Jdist::joint_dist_param *Ein_param )
{
  // This code is OK for unit-base data
  Ein_param->this_Ein_this_mu = Ein_param->this_Ein_ptr->begin( );
  Ein_param->this_Ein_next_mu = Ein_param->this_Ein_this_mu;
  ++Ein_param->this_Ein_next_mu;
  Ein_param->next_Ein_this_mu = Ein_param->next_Ein_ptr->begin( );
  Ein_param->next_Ein_next_mu = Ein_param->next_Ein_this_mu;
  ++Ein_param->next_Ein_next_mu;

  // Ein_param->lower_mu and Ein_param->upper_mu define the range of interpolation in unit-base mu
  Ein_param->lower_mu = 0.0;

  Ein_param->upper_mu = ( Ein_param->this_Ein_next_mu->get_mu( ) <  Ein_param->next_Ein_next_mu->get_mu( ) ) ?
    Ein_param->this_Ein_next_mu->get_mu( ) : Ein_param->next_Ein_next_mu->get_mu( );

  // save the incident energies for checking the geometry
  Ein_param->lower_mu0_hits.Ein_0 = Ein_param->this_Ein_ptr->get_E_in( );
  Ein_param->lower_mu0_hits.Ein_1 = Ein_param->next_Ein_ptr->get_E_in( );
  Ein_param->upper_mu0_hits.Ein_0 = Ein_param->this_Ein_ptr->get_E_in( );
  Ein_param->upper_mu0_hits.Ein_1 = Ein_param->next_Ein_ptr->get_E_in( );
  Ein_param->lower_mu0_hits.set_incident_range( Ein_param->data_E_0, Ein_param->data_E_1 );
  Ein_param->upper_mu0_hits.set_incident_range( Ein_param->data_E_0, Ein_param->data_E_1 );
  Ein_param->lower_mu1_hits.Ein_0 = Ein_param->this_Ein_ptr->get_E_in( );
  Ein_param->lower_mu1_hits.Ein_1 = Ein_param->next_Ein_ptr->get_E_in( );
  Ein_param->upper_mu1_hits.Ein_0 = Ein_param->this_Ein_ptr->get_E_in( );
  Ein_param->upper_mu1_hits.Ein_1 = Ein_param->next_Ein_ptr->get_E_in( );
  Ein_param->lower_mu1_hits.set_incident_range( Ein_param->data_E_0, Ein_param->data_E_1 );
  Ein_param->upper_mu1_hits.set_incident_range( Ein_param->data_E_0, Ein_param->data_E_1 );
}
// ----------- Jdist::joint_dist::next_mu_pairs --------------
// Go to the next pairs of direction cosines.  Returns "true" when finished.
bool Jdist::joint_dist::next_mu_pairs( Jdist::joint_dist_param *Ein_param )
{
  bool done = false;
  // Ein_param->lower_mu and Ein_param->upper_mu define the range of interpolation in unit-base mu
  Ein_param->lower_mu = Ein_param->upper_mu;  // the new lower limit for the integral over mu
  static double mu_tol = Global.Value( "tight_tol" );
  if( Ein_param->this_Ein_next_mu->get_mu( ) <= Ein_param->lower_mu*(1 + mu_tol ) )
  {
    Ein_param->this_Ein_this_mu = Ein_param->this_Ein_next_mu;
    ++Ein_param->this_Ein_next_mu;
    if( Ein_param->this_Ein_next_mu == Ein_param->this_Ein_ptr->end( ) )
    {
      return true;
    }
  }
  if( Ein_param->next_Ein_next_mu->get_mu( ) <= Ein_param->lower_mu*(1 + mu_tol ) )
  {
    Ein_param->next_Ein_this_mu = Ein_param->next_Ein_next_mu;
    ++Ein_param->next_Ein_next_mu;
    if( Ein_param->this_Ein_next_mu == Ein_param->this_Ein_ptr->end( ) )
    {
      return true;
    }
  }
  Ein_param->upper_mu = ( Ein_param->this_Ein_next_mu->get_mu( ) <
	       Ein_param->next_Ein_next_mu->get_mu( ) ) ?
    Ein_param->this_Ein_next_mu->get_mu( ) :
    Ein_param->next_Ein_next_mu->get_mu( );

  return done;
}
// ----------- Jdist::joint_dist::Eout_data_ladder --------------
// Adds to the transfer matrix for all E_out bins for given Ein and mu ranges
void Jdist::joint_dist::Eout_data_ladder( Trf::T_matrix& transfer,
					  Jdist::joint_dist_param *Ein_param )
{
  double dummy = 0.0;  // needed for Ein_param->lower_mu0_hits.hit_box
  if( mu_interp.qualifier == Terp::UNITBASE )
  {
    start_Eout_data_UB( Ein_param );
  }
  else // mu_interp.qualifier == CUMULATIVE_POINTS
  {
    start_Eout_data_CP( Ein_param );
  }
  
  // loop through the energy data
  for( ; ; )
  {
    // loop through the outgoing energy bins
    for( int Eout_count = 0; Eout_count < transfer.num_Eout_bins;
      ++Eout_count )
    {
      std::vector< double >::const_iterator Eout_ptr = transfer.out_groups.begin( )
        + Eout_count;
      // how do the lowest unit-base interpolation lines meet this E-E' box?
      Ein_param->lower_mu0_hits.hit_box( dummy, Eout_ptr, Ein_param->data_E_0, Ein_param->data_E_1 );
      Ein_param->lower_mu1_hits.hit_box( dummy, Eout_ptr, Ein_param->data_E_0, Ein_param->data_E_1 );
      if( ( Eout_count < transfer.num_Eout_bins - 1 ) &&
          ( Ein_param->lower_mu0_hits.is_above( ) ) &&
          ( Ein_param->lower_mu1_hits.is_above( ) ) )
      {
        // go on to the next E-E' box
        continue;
      }
      // how do the next unit-base interpolation lines meet this E-E' box?
      Ein_param->upper_mu0_hits.hit_box( dummy, Eout_ptr, Ein_param->data_E_0, Ein_param->data_E_1 );
      Ein_param->upper_mu1_hits.hit_box( dummy, Eout_ptr, Ein_param->data_E_0, Ein_param->data_E_1 );
      if( ( Eout_count > 0 ) && ( Ein_param->upper_mu0_hits.is_below( ) ) &&
	  ( Ein_param->upper_mu1_hits.is_below( ) ) )
      {
        // we are done with this pair of E values
        break;
      }
      // integrate over this E-E' box
      Ein_param->mu_params.Eout_bottom = Eout_ptr;
      Ein_param->mu_params.Eout_top = Eout_ptr + 1;
      one_box( transfer, Eout_count, Ein_param );
    }
    // go to the next pairs of (E_out, probability)
    bool done = false;
    if( mu_interp.qualifier == Terp::UNITBASE )
    {
      done =  Ein_param->next_Eout_UB( );
    }
    else if( mu_interp.qualifier == Terp::CUMULATIVE_POINTS )
    {
      done =  Ein_param->next_Eout_CP( );
    }
    if( done ) break;
  }
}
// ----------- Jdist::joint_dist::start_Eout_data_UB --------------
// Starts one staircase of the Eout data
// Sets up Ein_param->Ein0_data and Ein_param->Ein1_data
void Jdist::joint_dist::start_Eout_data_UB( Jdist::joint_dist_param *Ein_param )
{
  // This code is OK for unit-base data
  Ein_param->lower_Eout = 0.0;

  // where we are in the data
  Ein_param->Ein0_mu0_Eout0 = Ein_param->this_Ein_this_mu->begin( );
  Ein_param->Ein0_mu0_Eout1 = Ein_param->Ein0_mu0_Eout0;
  ++Ein_param->Ein0_mu0_Eout1;
  Ein_param->upper_Eout = Ein_param->Ein0_mu0_Eout1->x;

  Ein_param->Ein0_mu1_Eout0 = Ein_param->this_Ein_next_mu->begin( );
  Ein_param->Ein0_mu1_Eout1 = Ein_param->Ein0_mu1_Eout0;
  ++Ein_param->Ein0_mu1_Eout1;
  double this_Eout = Ein_param->Ein0_mu1_Eout1->x;
  if( Ein_param->upper_Eout > this_Eout )
  {
    Ein_param->upper_Eout = this_Eout;
  }

  Ein_param->Ein1_mu0_Eout0 = Ein_param->next_Ein_this_mu->begin( );
  Ein_param->Ein1_mu0_Eout1 = Ein_param->Ein1_mu0_Eout0;
  ++Ein_param->Ein1_mu0_Eout1;
  this_Eout = Ein_param->Ein1_mu0_Eout1->x;
  if( Ein_param->upper_Eout > this_Eout )
  {
    Ein_param->upper_Eout = this_Eout;
  }

  Ein_param->Ein1_mu1_Eout0 = Ein_param->next_Ein_next_mu->begin( );
  Ein_param->Ein1_mu1_Eout1 = Ein_param->Ein1_mu1_Eout0;
  ++Ein_param->Ein1_mu1_Eout1;
  this_Eout = Ein_param->Ein1_mu1_Eout1->x;
  if( Ein_param->upper_Eout > this_Eout )
  {
    Ein_param->upper_Eout = this_Eout;
  }

  // Set the incident energies
  double E_in = Ein_param->this_Ein_ptr->get_E_in( );
  Ein_param->Ein0_data.set_E_in( E_in );
  E_in = Ein_param->next_Ein_ptr->get_E_in( );
  Ein_param->Ein1_data.set_E_in( E_in );

  // For 3-d interpolation set up the Jdata::E_mu_P_data at given incident energies for
  // UB_Eout = Ein_param->lower_Eout, Ein_param->upper_Eout and UB_mu = Ein_param->lower_mu, Ein_param->upper_mu
  Ein_param->common_mu_Eout( Ein_param->lower_Eout, &Ein_param->Ein0_data.mu0_Eout0,
		  &Ein_param->Ein0_data.mu1_Eout0, &Ein_param->Ein1_data.mu0_Eout0,
		  &Ein_param->Ein1_data.mu1_Eout0 );
  // For UB_Eout = Ein_param->upper_Eout
  Ein_param->common_mu_Eout( Ein_param->upper_Eout, &Ein_param->Ein0_data.mu0_Eout1,
		  &Ein_param->Ein0_data.mu1_Eout1, &Ein_param->Ein1_data.mu0_Eout1,
		  &Ein_param->Ein1_data.mu1_Eout1 );

  // Save the outgoing energy ranges for the unit-base maps
  double alpha;
  if( Ein_param->lower_mu == Ein_param->this_Ein_this_mu->get_mu( ) )
  {
    Ein_param->Ein0_data.mu0_ubase_map.copy( Ein_param->this_Ein_this_mu->ubase_map );
  }
  else
  {
    alpha = ( Ein_param->lower_mu - Ein_param->this_Ein_this_mu->get_mu( ) )/
      ( Ein_param->this_Ein_next_mu->get_mu( ) - Ein_param->this_Ein_this_mu->get_mu( ) );
    Ein_param->Ein0_data.mu0_ubase_map.interpolate( alpha,
      Ein_param->this_Ein_this_mu->ubase_map, Ein_param->this_Ein_next_mu->ubase_map );
  }
  
  if( Ein_param->upper_mu == Ein_param->this_Ein_next_mu->get_mu( ) )
  {
    Ein_param->Ein0_data.mu1_ubase_map.copy( Ein_param->this_Ein_next_mu->ubase_map );
  }
  else
  {
    alpha = ( Ein_param->upper_mu - Ein_param->this_Ein_this_mu->get_mu( ) )/
      ( Ein_param->this_Ein_next_mu->get_mu( ) - Ein_param->this_Ein_this_mu->get_mu( ) );
    Ein_param->Ein0_data.mu1_ubase_map.interpolate( alpha,
      Ein_param->this_Ein_this_mu->ubase_map, Ein_param->this_Ein_next_mu->ubase_map );
  }

  // repeat at the higher incident energy
  if( Ein_param->lower_mu == Ein_param->next_Ein_this_mu->get_mu( ) )
  {
    Ein_param->Ein1_data.mu0_ubase_map.copy( Ein_param->next_Ein_this_mu->ubase_map );
  }
  else
  {
    alpha = ( Ein_param->lower_mu - Ein_param->next_Ein_this_mu->get_mu( ) )/
      ( Ein_param->next_Ein_next_mu->get_mu( ) - Ein_param->next_Ein_this_mu->get_mu( ) );
    Ein_param->Ein1_data.mu0_ubase_map.interpolate( alpha,
      Ein_param->next_Ein_this_mu->ubase_map, Ein_param->next_Ein_next_mu->ubase_map );
  }
  
  if( Ein_param->upper_mu == Ein_param->next_Ein_next_mu->get_mu( ) )
  {
    Ein_param->Ein1_data.mu1_ubase_map.copy( Ein_param->next_Ein_next_mu->ubase_map );
  }
  else
  {
    alpha = ( Ein_param->upper_mu - Ein_param->next_Ein_this_mu->get_mu( ) )/
      ( Ein_param->next_Ein_next_mu->get_mu( ) - Ein_param->next_Ein_this_mu->get_mu( ) );
    Ein_param->Ein1_data.mu1_ubase_map.interpolate( alpha, 
      Ein_param->next_Ein_this_mu->ubase_map, Ein_param->next_Ein_next_mu->ubase_map );
  }

  // Initialize the parameters for the intersection of the physical enegies with the quadrature box
  Ein_param->lower_mu0_hits.get_phys_Eout( );
  Ein_param->lower_mu1_hits.get_phys_Eout( );
  Ein_param->upper_mu0_hits.get_phys_Eout( );
  Ein_param->upper_mu1_hits.get_phys_Eout( );
}
// ----------- Jdist::joint_dist::start_Eout_data_CP --------------
// Starts one staircase of the Eout data
// Sets up Ein_param->Ein0_data and Ein_param->Ein1_data
void Jdist::joint_dist::start_Eout_data_CP( Jdist::joint_dist_param *Ein_param )
{
  // where we are in the data
  // for lower incident energy and lower mu
  Ein_param->Ein0_mu0_Eout0_CP = Ein_param->this_Ein_this_mu->cum_prob.begin( );
  Ein_param->Ein0_mu0_Eout1_CP = Ein_param->Ein0_mu0_Eout0_CP;
  ++Ein_param->Ein0_mu0_Eout1_CP;

  // skip intervals with zero probability
  while( ( Ein_param->Ein0_mu0_Eout0_CP->Prob == 0.0 ) &&
         ( Ein_param->Ein0_mu0_Eout0_CP->slope == 0.0 ) )
  {
    Ein_param->Ein0_mu0_Eout0_CP = Ein_param->Ein0_mu0_Eout1_CP;
    ++Ein_param->Ein0_mu0_Eout1_CP;
  }

  // for lower incident energy and higher mu
  Ein_param->Ein0_mu1_Eout0_CP = Ein_param->this_Ein_next_mu->cum_prob.begin( );
  Ein_param->Ein0_mu1_Eout1_CP = Ein_param->Ein0_mu1_Eout0_CP;
  ++Ein_param->Ein0_mu1_Eout1_CP;

  // skip intervals with zero probability
  while( ( Ein_param->Ein0_mu1_Eout0_CP->Prob == 0.0 ) &&
         ( Ein_param->Ein0_mu1_Eout0_CP->slope == 0.0 ) )
  {
    Ein_param->Ein0_mu1_Eout0_CP = Ein_param->Ein0_mu1_Eout1_CP;
    ++Ein_param->Ein0_mu1_Eout1_CP;
  }

  // for higher incident energy and lower mu
  Ein_param->Ein1_mu0_Eout0_CP = Ein_param->next_Ein_this_mu->cum_prob.begin( );
  Ein_param->Ein1_mu0_Eout1_CP = Ein_param->Ein1_mu0_Eout0_CP;
  ++Ein_param->Ein1_mu0_Eout1_CP;

  // skip intervals with zero probability
  while( ( Ein_param->Ein1_mu0_Eout0_CP->Prob == 0.0 ) &&
         ( Ein_param->Ein1_mu0_Eout0_CP->slope == 0.0 ) )
  {
    Ein_param->Ein1_mu0_Eout0_CP = Ein_param->Ein1_mu0_Eout1_CP;
    ++Ein_param->Ein1_mu0_Eout1_CP;
  }

  // for higher incident energy and higher mu
  Ein_param->Ein1_mu1_Eout0_CP = Ein_param->next_Ein_next_mu->cum_prob.begin( );
  Ein_param->Ein1_mu1_Eout1_CP = Ein_param->Ein1_mu1_Eout0_CP;
  ++Ein_param->Ein1_mu1_Eout1_CP;

  // skip intervals with zero probability
  while( ( Ein_param->Ein1_mu1_Eout0_CP->Prob == 0.0 ) &&
         ( Ein_param->Ein1_mu1_Eout0_CP->slope == 0.0 ) )
  {
    Ein_param->Ein1_mu1_Eout0_CP = Ein_param->Ein1_mu1_Eout1_CP;
    ++Ein_param->Ein1_mu1_Eout1_CP;
  }

  // Set the incident energies
  double E_in = Ein_param->this_Ein_ptr->get_E_in( );
  Ein_param->Ein0_data.set_E_in( E_in );
  E_in = Ein_param->next_Ein_ptr->get_E_in( );
  Ein_param->Ein1_data.set_E_in( E_in );

  // the range of mu values
  Ein_param->lower_mu = ( Ein_param->this_Ein_this_mu->get_mu( ) >
			  Ein_param->next_Ein_this_mu->get_mu( ) ) ?
    Ein_param->this_Ein_this_mu->get_mu( ) :
    Ein_param->next_Ein_this_mu->get_mu( );
  
  Ein_param->upper_mu = ( Ein_param->this_Ein_next_mu->get_mu( ) <
			  Ein_param->next_Ein_next_mu->get_mu( ) ) ?
    Ein_param->this_Ein_next_mu->get_mu( ) :
    Ein_param->next_Ein_next_mu->get_mu( );

  // Set the mu data
  Ein_param->this_Ein_ptr->set_mu_data( Ein_param->lower_mu, Ein_param->this_Ein_this_mu,
      Ein_param->this_Ein_next_mu, &Ein_param->Ein0_data, 0 );
  
  Ein_param->this_Ein_ptr->set_mu_data( Ein_param->upper_mu, Ein_param->this_Ein_this_mu,
      Ein_param->this_Ein_next_mu, &Ein_param->Ein0_data, 1 );

  Ein_param->next_Ein_ptr->set_mu_data( Ein_param->lower_mu, Ein_param->next_Ein_this_mu,
      Ein_param->next_Ein_next_mu, &Ein_param->Ein1_data, 0 );
  
  Ein_param->next_Ein_ptr->set_mu_data( Ein_param->upper_mu, Ein_param->next_Ein_this_mu,
      Ein_param->next_Ein_next_mu, &Ein_param->Ein1_data, 1 );

  // the range of cumulative probabilities
  Ein_param->lower_E_cum = 0.0;
  
  double upper_E_cum0 = ( Ein_param->Ein0_mu0_Eout1_CP->cum_prob <
			  Ein_param->Ein0_mu1_Eout1_CP->cum_prob ) ?
    Ein_param->Ein0_mu0_Eout1_CP->cum_prob :
    Ein_param->Ein0_mu1_Eout1_CP->cum_prob;

  double upper_E_cum1 = ( Ein_param->Ein1_mu0_Eout1_CP->cum_prob <
			  Ein_param->Ein1_mu1_Eout1_CP->cum_prob ) ?
    Ein_param->Ein1_mu0_Eout1_CP->cum_prob :
    Ein_param->Ein1_mu1_Eout1_CP->cum_prob;

  Ein_param->upper_E_cum = ( upper_E_cum0 < upper_E_cum1 ) ?
    upper_E_cum0 : upper_E_cum1;

  // initialize Ein_param->Ein0_data
  Ein_param->set_cum_data( Ein_param->Ein0_mu0_Eout0_CP,
			   Ein_param->Ein0_mu1_Eout0_CP,
			   &Ein_param->Ein0_data );

  // initialize Ein_param->Ein1_data
  Ein_param->set_cum_data( Ein_param->Ein1_mu0_Eout0_CP,
			   Ein_param->Ein1_mu1_Eout0_CP,
			   &Ein_param->Ein1_data );
  
  // Initialize the parameters for the intersection of the physical enegies with the quadrature box
  Ein_param->lower_mu0_hits.get_phys_Eout( );
  Ein_param->lower_mu1_hits.get_phys_Eout( );
  Ein_param->upper_mu0_hits.get_phys_Eout( );
  Ein_param->upper_mu1_hits.get_phys_Eout( );
}

// **************** Functions to integrate *********************
// ---------------- joint_dist_F::Eout_F ------------------
// Function for the 1-d quadrature over outgoing energy
bool joint_dist_F::Eout_F( double Eout, Qparam::QuadParamBase *Eout_quad_param,
			   Coef::coef_vector *value )
{
  // the parameters are really Jdist::joint_Eout_param
  Jdist::joint_Eout_param *Eout_params =
    static_cast< Jdist::joint_Eout_param* >( Eout_quad_param );
  Eout_params->func_count += 1;

  // for interpolation
  Jdata::E_mu_P_data mid_data;
  if( Eout_params->Eout_interp ==Terp:: LINLIN )
  {
    bool is_OK = Eout_params->Eout0_data.UB_Eout_interp( Eout, Eout_params->Eout1_data,
					     &mid_data);
    if( !is_OK ) return false;
  }
  else
  {
    Eout_params->Eout0_data.UB_Eout_histogram( Eout, Eout_params->Eout1_data,
					     &mid_data);
  }

  value->weight_1[ 0 ] = mid_data.Eout_Prob * mid_data.mu_Prob;
  value->weight_E[ 0 ] = value->weight_1[ 0 ] * mid_data.phys_Eout;

  return true;
}
// ---------------- joint_dist_F::mu_F ------------------
// Function for the 2-d quadrature over cosine and outgoing energy
bool joint_dist_F::mu_F( double mu, Qparam::QuadParamBase *mu_quad_param,
			   Coef::coef_vector *value )
{
  // the parameters are really Jdist::joint_mu_param
  Jdist::joint_mu_param *mu_params =
    static_cast< Jdist::joint_mu_param* >( mu_quad_param );
  mu_params->func_count += 1;

  // interpolate data
  Ddvec::unit_base_map mid_ubase_map;
  bool is_OK = mu_params->this_data.mu_interpolate( mu, &mu_params->Eout_params.Eout0_data,
					&mu_params->Eout_params.Eout1_data,
					&mid_ubase_map );
  if( !is_OK ) return false;

  // the range of integration
  double Eout0;
  if( mu_params->Eout_params.Eout0_data.phys_Eout < *mu_params->Eout_bottom )
  {
    Eout0 = mid_ubase_map.to_unit_base( *mu_params->Eout_bottom, &is_OK );
    if( !is_OK ) return false;
  }
  else
  {
    Eout0 = mu_params->Eout_params.Eout0_data.UB_Eout;
  }

  double Eout1;
  if( mu_params->Eout_params.Eout1_data.phys_Eout > *mu_params->Eout_top )
  {
    Eout1 = mid_ubase_map.to_unit_base( *mu_params->Eout_top, &is_OK );
    if( !is_OK ) return false;
  }
  else
  {
    Eout1 = mu_params->Eout_params.Eout1_data.UB_Eout;
  }

  // for integration over outgoing energy
  Qparam::QuadParamBase *params = static_cast< Qparam::QuadParamBase* >( &mu_params->Eout_params );
  // the mu dependence is scalar, but we might need the weight by ougtoing energy
  Coef::coef_vector Eout_integral( 0, Coef::BOTH );
  static double tol = Global.Value( "quad_tol" );
  // 2nd-order Gaussian quadrature is exact
  Qmeth::Quadrature_Rule quad_rule;
  quad_rule.adaptive = false;
  quad_rule.quad_method = Qmeth::GAUSS2;
  is_OK = quad_F::integrate( joint_dist_F::Eout_F, quad_rule, Eout0, Eout1,
		     params, tol, &Eout_integral );

  // the Legendre polynomials
  math_F::Legendre( mu_params->Eout_params.Eout0_data.phys_mu, value );
  // scale by the integrals over outgoing energy
  double Prob = Eout_integral.weight_1[ 0 ];
  *value *= Prob;
  // do the energy weighting if necessary
  if( ( Prob > 0.0 ) &&
      ( ( value->conserve == Coef::ENERGY ) || ( value->conserve == Coef::BOTH ) ) )
  {
    value->scale_E( Eout_integral.weight_E[ 0 ]/Prob );
  }
  mu_params->Eout_F_count += mu_params->Eout_params.func_count;

  return is_OK;
}
// ----------- joint_dist_F::Ein_F ------------------
// Function for the 3-d quadrature over incident energy, cosine, and outgoing energy
bool joint_dist_F::Ein_F( double E_in, Qparam::QuadParamBase *Ein_quad_param,
  Coef::coef_vector *value )
{
  // the parameters are really Jdist::joint_dist_param
  Jdist::joint_dist_param *Ein_params =
    static_cast< Jdist::joint_dist_param* >( Ein_quad_param );
  Ein_params->func_count += 1;

  // Interpolate the data in incident energy
  bool is_OK = Ein_params->Ein0_data.Ein_interpolate( E_in, Ein_params->Ein1_data,
					 &Ein_params->mu_params.this_data );
  if( !is_OK ) return false;

  // parameters for the integration over mu
  Qparam::QuadParamBase *params = static_cast< Qparam::QuadParamBase* >( &Ein_params->mu_params );
  static double tol = Global.Value( "quad_tol" );

  // We are integrating a quartic in mu, so 4-th order Gauss is exact
  Qmeth::Quadrature_Rule mu_quad_rule;
  mu_quad_rule.adaptive = false;
  mu_quad_rule.quad_method = Qmeth::GAUSS4;

  // to hold the value of the integral over a mu subinterval
  Coef::coef_vector temp_value( value->order, value->conserve );

  // Set the value to zero in case the geometry is bad
  value->set_zero( );

  // Determine how the data meets the outgoing physical energy range
  if( Ein_params->geometry( Ein_params->mu_params.this_data ) )
  {
    Box::energy_hit_list::iterator low_hit_ptr = Ein_params->lower_2d_hits.begin( );
    Box::energy_hit_list::iterator next_low_ptr = low_hit_ptr;
    ++next_low_ptr;
    Box::energy_hit_list::iterator high_hit_ptr = Ein_params->upper_2d_hits.begin( );
    Box::energy_hit_list::iterator next_high_ptr = high_hit_ptr;
    ++next_high_ptr;
    for( ; ( next_low_ptr != Ein_params->lower_2d_hits.end( ) ) &&
           ( next_high_ptr != Ein_params->upper_2d_hits.end( ) );
         low_hit_ptr = next_low_ptr, ++next_low_ptr,
           high_hit_ptr = next_high_ptr, ++next_high_ptr )
    {
      if( ( low_hit_ptr->hit_edge == Box::ABOVE ) ||
          ( low_hit_ptr->hit_edge == Box::TOP_OUT ) )
      {
        // do nothing---we are above the E-E' box
        continue;
      }
      else if( ( high_hit_ptr->hit_edge == Box::BOTTOM_OUT ) ||
               ( high_hit_ptr->hit_edge == Box::BELOW ) )
      {
        // do nothing---we are below the E-E' box
        continue;
      }
      else
      {
        bool one_OK = quad_F::integrate( joint_dist_F::mu_F, mu_quad_rule,
                         low_hit_ptr->E_in, next_low_ptr->E_in, 
                         params, tol, &temp_value );
	if( !one_OK )
	{
	  is_OK = false;
	}
	*value += temp_value;
        Ein_params->mu_F_count += Ein_params->mu_params.func_count;
        Ein_params->Eout_F_count += Ein_params->mu_params.Eout_F_count;
	Ein_params->mu_params.Eout_F_count = 0;
      }
    }
  }
  // weight it by flux * cross section * multiplicity * model weight
  Ein_params->set_weight( E_in );
  *value *= Ein_params->current_weight;

  return is_OK;
}
