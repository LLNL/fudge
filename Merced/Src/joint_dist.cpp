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
#include "messaging.hpp"
#include "global_params.hpp"

// ********* class E_mu_P_data *********
// ----------- E_mu_P_data::UB_Eout_interp --------------
// Interpolates with respect to unit-base outgoing energy
void E_mu_P_data::UB_Eout_interp( double mid_UB_Eout,
   const E_mu_P_data &next_data,
		       E_mu_P_data *mid_data ) const
{
  // check the inputs
  double E_diff = next_data.UB_Eout - UB_Eout;
  if( E_diff <= 0.0 )
  {
    FatalError( "E_mu_P_data::UB_Eout_interp", "data out of order" );
  }
  double alpha = ( mid_UB_Eout - UB_Eout )/E_diff;
  static double etol = Global.Value( "E_tol" );
  static bool no_warning = true;
  if( ( ( alpha < 0.0 ) || ( alpha > 1.0 ) ) &&
      ( E_diff > etol ) && no_warning )
  {
    Warning( "E_mu_P_data::UB_Eout_interp", "extrapolation" );
    no_warning = false;
  }
  mid_data->UB_mu = UB_mu;
  mid_data->phys_mu = phys_mu;
  mid_data->UB_Eout = mid_UB_Eout;
  mid_data->phys_Eout = ( 1.0 - alpha )*phys_Eout +
    alpha*next_data.phys_Eout;
  mid_data->Prob = ( 1.0 - alpha )*Prob + alpha*next_data.Prob;
}
// ----------- E_mu_P_data::UB_Eout_histogram --------------
// Interpolates histogram data with respect to unit-base outgoing energy
void E_mu_P_data::UB_Eout_histogram( double mid_UB_Eout,
    const E_mu_P_data &next_data,
		       E_mu_P_data *mid_data ) const
{
  // check the inputs
  double E_diff = next_data.UB_Eout - UB_Eout;
  if( E_diff <= 0.0 )
  {
    FatalError( "E_mu_P_data::UB_Eout_histogram", "data out of order" );
  }
  double alpha = ( mid_UB_Eout - UB_Eout )/E_diff;
  static double etol = Global.Value( "E_tol" );
  static bool no_warning = true;
  if( ( ( alpha < 0.0 ) || ( alpha > 1.0 ) ) &&
      ( E_diff > etol ) && no_warning )
  {
    Warning( "E_mu_P_data::UB_Eout_histogram", "extrapolation" );
    no_warning = false;
  }
  mid_data->UB_mu = UB_mu;
  mid_data->phys_mu = phys_mu;
  mid_data->UB_Eout = mid_UB_Eout;
  mid_data->phys_Eout = ( 1.0 - alpha )*phys_Eout +
    alpha*next_data.phys_Eout;
  mid_data->Prob = Prob;
}
// ----------- E_mu_P_data::mu_interp --------------
// Interpolates with respect to direction cosine
void E_mu_P_data::mu_interp( double mid_mu, const E_mu_P_data &next_data,
		       E_mu_P_data *mid_data ) const
{
  // check the inputs
  double mu_diff = next_data.UB_mu - UB_mu;
  if( mu_diff <= 0.0 )
  {
    FatalError( "E_mu_P_data::mu_interp", "data out of order" );
  }
  double alpha = ( mid_mu - UB_mu )/mu_diff;
  static double mu_tol = Global.Value( "E_tol" );
  static bool no_warning = true;
  if( ( ( alpha <= -mu_tol ) || ( alpha >= 1.0 + mu_tol ) ) && no_warning )
  {
    Warning( "E_mu_P_data::mu_interp", "extrapolation" );
    no_warning = false;
  }
  mid_data->UB_mu = mid_mu;
  mid_data->phys_mu = ( 1.0 - alpha )*phys_mu + alpha*next_data.phys_mu;
  mid_data->UB_Eout = UB_Eout;
  mid_data->phys_Eout = ( 1.0 - alpha )*phys_Eout +
    alpha*next_data.phys_Eout;
  mid_data->Prob = ( 1.0 - alpha )*Prob + alpha*next_data.Prob;
}
// ----------- E_mu_P_data::Ein_interp --------------
// Interpolates with respect to unit-base outgoing energy
void E_mu_P_data::Ein_interp( double alpha, const E_mu_P_data &next_data,
		       E_mu_P_data *mid_data ) const
{
  // check the inputs
  if( next_data.UB_mu != UB_mu )
  {
    FatalError( "E_mu_P_data::Ein_interp",
		"interpolation valid only for idential unit-base mu values" );
  }

  static bool no_warning = true;
  if( ( ( alpha < 0.0 ) || ( alpha > 1.0 ) ) && no_warning )
  {
    Warning( "E_mu_P_data::Ein_interp", "extrapolation" );
    no_warning = false;
  }
  //  mid_data->UB_mu = ( 1.0 - alpha )*UB_mu + alpha*next_data.UB_mu;
  mid_data->UB_mu = UB_mu;
  mid_data->phys_mu = ( 1.0 - alpha )*phys_mu + alpha*next_data.phys_mu;
  mid_data->UB_Eout = UB_Eout;
  mid_data->phys_Eout = ( 1.0 - alpha )*phys_Eout +
    alpha*next_data.phys_Eout;
  mid_data->Prob = ( 1.0 - alpha )*Prob + alpha*next_data.Prob;
}
// ----------- E_mu_P_data::copy --------------
// Copies the data
void E_mu_P_data::copy( const E_mu_P_data &to_copy )
{
  UB_mu = to_copy.UB_mu;
  phys_mu = to_copy.phys_mu;
  UB_Eout = to_copy.UB_Eout;
  phys_Eout = to_copy.phys_Eout;
  Prob = to_copy.Prob;
}

// ********* class current_data *********
// ----------- current_data::Ein_interpolate --------------
// Interpolate in incident energy, returns mid_data
void current_data::Ein_interpolate( double mid_Ein,
                    const current_data &next_data,
                    current_data *mid_data ) const
{
  // check the inputs
  double E_diff = next_data.get_E_in( ) - get_E_in( );
  if( E_diff <= 0.0 )
  {
    FatalError( "current_data::Ein_interpolate", "data out of order" );
  }
  double alpha = ( mid_Ein - get_E_in( ) )/E_diff;
  static bool no_warning = true;
  if( ( ( alpha < 0.0 ) || ( alpha > 1.0 ) ) && no_warning )
  {
    Warning( "current_data::Ein_interpolate", "extrapolation" );
    no_warning = false;
  }
  mid_data->set_E_in( mid_Ein );
  mu0_Eout0.Ein_interp( alpha, next_data.mu0_Eout0, &( mid_data->mu0_Eout0 ) );
  mu0_Eout1.Ein_interp( alpha, next_data.mu0_Eout1, &( mid_data->mu0_Eout1 ) );
  mu1_Eout0.Ein_interp( alpha, next_data.mu1_Eout0, &( mid_data->mu1_Eout0 ) );
  mu1_Eout1.Ein_interp( alpha, next_data.mu1_Eout1, &( mid_data->mu1_Eout1 ) );

  // interpolate the unit-base map ranges of outgoing energy
  mid_data->mu0_ubase_map.interpolate( alpha, mu0_ubase_map,
				       next_data.mu0_ubase_map );
  mid_data->mu1_ubase_map.interpolate( alpha, mu1_ubase_map,
				       next_data.mu1_ubase_map );
}
// ----------- current_data::mu_interpolate --------------
// Interpolate in mu, returns Eout0_data, Eout1_data, and mid_ubase_map
void current_data::mu_interpolate( double mid_mu, E_mu_P_data *Eout0_data,
                    E_mu_P_data *Eout1_data, unit_base_map *mid_ubase_map ) const
{
  mu0_Eout0.mu_interp( mid_mu, mu1_Eout0, Eout0_data );
  mu0_Eout1.mu_interp( mid_mu, mu1_Eout1, Eout1_data );

  // interpolate the unit-base map
  double alpha = ( mid_mu - mu0_Eout0.UB_mu )/( mu1_Eout0.UB_mu - mu0_Eout0.UB_mu );
   mid_ubase_map->interpolate( alpha, mu0_ubase_map, mu1_ubase_map );
}

// ********* class joint_dist_hits *********
// ----------- joint_dist_hits::set_incident_range --------------
// Sets the range of integration over incident energy
void joint_dist_hits::set_incident_range( double E0, double E1 )
{
  E_Eout.first.x = E0;
  E_Eout.second.x = E1;
}
// ----------- joint_dist_hits::get_phys_Eout --------------
// Gets the physical outgoing energies at the ends of the working energy interval
void joint_dist_hits::get_phys_Eout( )
{
  // check the incident energies
  double Ein_diff = Ein_1 - Ein_0;
  if( Ein_diff <= 0.0 )
  {
    FatalError( "joint_dist_hits::get_phys_Eout",
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

// ********* class joint_dist_param *********
// ----------- joint_dist_param::joint_dist_param --------------
// Constructor
joint_dist_param::joint_dist_param( )
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
// ----------- joint_dist_param::geometry --------------
// determines the geometry for 2d integration over outgoing cosine and energy
// returns true if the geometry makes sense
bool joint_dist_param::geometry( const current_data &this_data )
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
    Warning( "joint_dist_param::geometry", "all outgoing energies too high" );
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
    Warning( "joint_dist_param::geometry", "all outgoing energies too low" );
    geom_OK = false;
    no_upper_warning = false;
  }

  // set up common mu values
  lower_2d_hits.common_hits( upper_2d_hits );

  return geom_OK;
}
// ----------- joint_dist_param::common_mu_Eout --------------
// Interpolates data to common values of unit-base mu and outgoing energy
void joint_dist_param::common_mu_Eout( double UB_Eout, E_mu_P_data *Ein0_mu0_data,
  E_mu_P_data *Ein0_mu1_data, E_mu_P_data *Ein1_mu0_data, E_mu_P_data *Ein1_mu1_data )
{
  // set up the E_mu_P_data at given incident energy, mu and outgoing energy
  // temporary storage for interpolation in unit-base mu
  E_mu_P_data mu0_ENDF;  // for UB_Eout data at lower ENDF mu
  E_mu_P_data mu1_ENDF;  // for UB_Eout data at higher ENDF mu

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
    mu0_ENDF.mu_interp( lower_mu, mu1_ENDF, Ein0_mu0_data );
  }
  if( upper_mu == mu1_ENDF.UB_mu )
  {
    Ein0_mu1_data->copy( mu1_ENDF );
  }
  else
  {
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
    mu0_ENDF.mu_interp( lower_mu, mu1_ENDF, Ein1_mu0_data );
  }
  if( upper_mu == mu1_ENDF.UB_mu )
  {
    Ein1_mu1_data->copy( mu1_ENDF );
  }
  else
  {
    mu0_ENDF.mu_interp( upper_mu, mu1_ENDF, Ein1_mu1_data );
  }
}

// ********* class one_mu *********
// ----------- one_mu::copy --------------
// Make a copy for the ENDL_kludge routine
void one_mu::copy( const one_mu &to_copy )
{
  // copy the basic information
  set_mu( to_copy.get_mu( ) );
  interp_type = to_copy.interp_type;
  ubase_map.copy( to_copy.ubase_map );

  // copy the entries
  for( one_mu::const_iterator this_entry = to_copy.begin( );
         this_entry != to_copy.end( ); ++this_entry )
  {
    add_entry( this_entry->x, this_entry->y );
  }
}
// ----------- one_mu::set_E_mu_P_data --------------
// Sets up E_mu_P_data at unit-base outgoing energy UB_Eout
void one_mu::set_E_mu_P_data( double UB_Eout, double phys_mu,
    one_mu::const_iterator prev_data,
    one_mu::const_iterator next_data, E_mu_P_data *mid_data )
{
  mid_data->UB_mu = get_mu( );
  mid_data->phys_mu = phys_mu;
  mid_data->UB_Eout = UB_Eout;
  mid_data->phys_Eout = ubase_map.un_unit_base( UB_Eout );
  static double etol = Global.Value( "E_tol" );
  if( interp_type == HISTOGRAM )
  {
    if( UB_Eout > ( 1 - etol )*next_data->x )
    {
      mid_data->Prob = next_data->y;
    }
    else
    {
      mid_data->Prob = prev_data->y;
    }
  }
  else
  {
    mid_data->Prob = prev_data->linlin_interp( UB_Eout, *next_data );
  }
}

// ********* class one_joint_dist *********
// ----------- one_joint_dist::copy --------------
void one_joint_dist::copy( const one_joint_dist &to_copy )
// make a copy at the threshold, used by the ENDL_kludge routine
{
  for( one_joint_dist::const_iterator copy_ptr = to_copy.begin( );
       copy_ptr != to_copy.end( ); ++copy_ptr )
  {
    // make a new energy distribution
    one_joint_dist::iterator new_e_dist_ptr = insert( end( ), one_mu( ) );
    // copy into it
    new_e_dist_ptr->copy( *copy_ptr );
  }
}
// ----------- one_joint_dist::to_ENDF --------------
// Converts double-differential data from ENDL to ENDF format
void one_joint_dist::to_ENDF( angle_dist::const_iterator &angles )
{
  // first check the incident energy
  if( get_E_in( ) != angles->get_E_in( ) )
  {
    FatalError("one_joint_dist::check_data",
      pastenum("E_in for i=1: ", angles->get_E_in( ) ) +
      pastenum(" different from E_in for i=3: ", get_E_in( ) ) );
  }

  one_joint_dist::iterator joint_mu = begin( );
  one_joint_dist::iterator prev_joint_mu = begin( );
  dd_vector::const_iterator angle_mu = angles->begin( );

  // First, check to see whether we need to interpolate the table to intermediate mu
  for( ; ( joint_mu != end( ) ) && ( angle_mu != angles->end( ) ); ++angle_mu )
  {
    if( ( joint_mu->get_mu( ) < angle_mu->x ) ||
        ( ( joint_mu->get_mu( ) > angle_mu->x ) && ( joint_mu == prev_joint_mu ) ) )
    {
      FatalError("one_joint_dist::to_ENDF ",
        pastenum(" mu for i=1: ", angle_mu->x ) +
        pastenum(" different from for i=3: ", joint_mu->get_mu( ) ) );
    }
    else if( joint_mu->get_mu( ) > angle_mu->x )
    {
      // insert energy probability densities for an intermediate mu
      interpolate_mu( angle_mu->x, prev_joint_mu, joint_mu );
      ++prev_joint_mu;  // points to the inserted intermediate list
    }
    else
    {
      // they match; scale by the angular probability density
      prev_joint_mu = joint_mu;  // go to the next mu
      ++joint_mu;
    }
  }

  // now do the scaling
  joint_mu = begin( );
  angle_mu = angles->begin( );
  for( ; ( joint_mu != end( ) ) && ( angle_mu != angles->end( ) ); ++angle_mu )
  {
    if( ( joint_mu->get_mu( ) < angle_mu->x ) ||
        ( joint_mu->get_mu( ) > angle_mu->x ) )
    {
      FatalError("one_joint_dist::to_ENDF ",
        pastenum(" mu for i=1: ", angle_mu->x ) +
        pastenum(" different from for table: ", joint_mu->get_mu( ) ) );
    }
    else
    {
      // they match; scale by the angular probability density
      *joint_mu *= angle_mu->y;
      ++joint_mu;  // go to the next mu
    }
  }
}
// ----------- one_joint_dist::interpolate_mu --------------
// insert energy probability densities for an intermediate mu
void one_joint_dist::interpolate_mu( double mu,
  one_joint_dist::iterator &prev_joint_mu,
  one_joint_dist::iterator &joint_mu )
{
  one_joint_dist::iterator new_mu_ptr = insert( joint_mu, one_mu( ) );
  new_mu_ptr->interpolate( mu, *prev_joint_mu, *joint_mu );
  if( mu_interp.qualifier == UNITBASE )
  {
    double mu_diff = joint_mu->get_mu( ) - prev_joint_mu->get_mu( );
    if( mu_diff <= 0.0 )
    {
      FatalError( "one_joint_dist::interpolate_mu",
		  "mu values out of order" );
    }
    double alpha = ( mu - prev_joint_mu->get_mu( ) )/mu_diff;
    new_mu_ptr->ubase_map.interpolate( alpha, prev_joint_mu->ubase_map,
				       joint_mu->ubase_map );
  }
}
// ----------- one_joint_dist::to_unit_base --------------
// Maps the direction cosines to 0 <= mu <= 1
void one_joint_dist::to_unit_base( )
{
  one_joint_dist::iterator this_mu = begin( );
  one_joint_dist::iterator last_mu = end( );
  --last_mu;
  mu_ubase_map.Eout_min = this_mu->get_mu( );
  mu_ubase_map.Eout_max = last_mu->get_mu( );
  double scale_mu = mu_ubase_map.Eout_max - mu_ubase_map.Eout_min;

  for( ; this_mu != end( ); ++this_mu )
  {
    this_mu->set_mu( mu_ubase_map.to_unit_base( this_mu->get_mu( ) ) );
    *this_mu *= scale_mu;
  }
}

// ********* class joint_dist *********
// ----------- joint_dist::read_data --------------
// Reads double-differential energy-angle data
void joint_dist::read_data( data_parser &inFile, int num_Ein )
{
  interp_flag_F::read_3d_interpolation( inFile, &Ein_interp, &mu_interp,
                       &Eout_interp );

  joint_dist::iterator new_Ein_ptr;
  one_joint_dist::iterator new_mu_ptr;

  double E_out;
  double Prob;

  // read the data
  for( int Ein_count = 0; Ein_count < num_Ein; ++Ein_count )
  {
    // make a new link for this incident energy
    new_Ein_ptr = insert( end( ), one_joint_dist( ) );
    new_Ein_ptr->set_E_in( inFile.get_next_double( ) );
    new_Ein_ptr->mu_interp = mu_interp;
    new_Ein_ptr->Ein_interp = Ein_interp;
    new_Ein_ptr->Eout_interp = Eout_interp;

    // loop over cosine
    int num_mu = inFile.get_next_int( );
    for( int mu_count = 0; mu_count < num_mu; ++mu_count )
    {
      // make a new energy distribution
      new_mu_ptr = new_Ein_ptr->insert( new_Ein_ptr->end( ), one_mu( ) );
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
      if( mu_interp.qualifier == UNITBASE )
      {
        new_mu_ptr->unit_base( false, &new_mu_ptr->ubase_map );
      }
    }
  }
  // we may need to copy the first distribution at the threshold
  ENDL_kludge( );
  if( ENDL_data )
  {
    to_ENDF( );  // convert to ENDF format
  }
  if( Ein_interp.qualifier == UNITBASE )
  {
    to_unit_base( );
  }
}
// ----------- joint_dist::ENDL_kludge --------------
void joint_dist::ENDL_kludge( )
// ENDL has the convention that the distribution is assumed
// independent of energy at low incident energies.
{
  joint_dist::iterator first_dist = begin( );
  if( first_dist->get_E_in( ) > angle_data.threshold )
  {
    // prepend a copy at the threshold energy
    joint_dist::iterator second_dist = begin( );
    first_dist = insert( begin( ), one_joint_dist( ) );
    first_dist->copy( *second_dist );
    first_dist->set_E_in( angle_data.threshold );
  }
}
// ----------- joint_dist::to_ENDF --------------
// Converts double-differential data from ENDL to ENDF format
void joint_dist::to_ENDF( )
{
  joint_dist::iterator joint_mu = begin( );
  angle_dist::const_iterator angle_mu = angle_data.begin( );

  for( ; ( joint_mu != end( ) ) && ( angle_mu != angle_data.end( ) );
       ++joint_mu, ++angle_mu )
  {
    joint_mu->to_ENDF( angle_mu );
  }
}
// ----------- joint_dist::to_unit_base --------------
// Maps the direction cosines to 0 <= mu <= 1
void joint_dist::to_unit_base( )
{
  joint_dist::iterator joint_mu = begin( );
  for( ; joint_mu != end( ); ++joint_mu )
  {
    joint_mu->to_unit_base( );
  }
}
// ----------- joint_dist::get_T --------------
// Calculates the transfer matrix for this particle.
void joint_dist::get_T( const dd_vector& sigma, const dd_vector& multiple,
  const dd_vector& weight, T_matrix& transfer )
{
  if( ( Ein_interp.qualifier != UNITBASE ) || ( Ein_interp.flag != LINLIN ) )
  {
    FatalError( "joint_dist::get_T", "Ein_interp not implemented" );
  }
  if( ( mu_interp.qualifier != UNITBASE ) || ( mu_interp.flag != LINLIN ) )
  {
    FatalError( "joint_dist::get_T", "mu_interp not implemented" );
  }
  if( ( Eout_interp != LINLIN ) && ( Eout_interp != HISTOGRAM ) )
  {
    FatalError( "joint_dist::get_T", "Eout_interp not implemented" );
  }
  if( ENDL_data && ( angle_data.size( ) == 0 ) )
  {
    FatalError( "joint_dist::get_T",
      "ENDL table missing angular data" );
  }
  if( ( angle_data.size( ) > 0 ) &&
      ( Ein_interp.qualifier != angle_data.Ein_interp.qualifier ) )
  {
    FatalError( "joint_dist::get_T",
      "Incident energy interpolation not consistent with the angular data" );
  }

  bool done = get_Ein_range( sigma, multiple, weight, transfer.e_flux,
    transfer.in_groups );
  if( done )
  {
    transfer.zero_transfer( );
  }

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
    joint_dist_param Ein_param;
    // set up the data range for this bin
    Ein_param.setup_bin( Ein_bin, sigma, multiple, weight, transfer.e_flux,
                         transfer.in_groups );
    Ein_param.Eout_quad_method = transfer.Eout_quad_method;
    Ein_param.mu_quad_method = transfer.mu_quad_method;
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
  cout << "3d quadratures: " << quad_count << endl;
  cout << "joint_dist_F::Ein_F calls: " << Ein_F_count << endl;
  cout << "joint_dist_F::mu_F calls: " << mu_F_count << endl;
  cout << "joint_dist_F::Eout_F calls: " << Eout_F_count << endl;
  cout << "average joint_dist_F::Ein_F_count: " << 1.0*Ein_F_count/quad_count << endl;
  cout << "average joint_dist_F::mu_F_count: " << 1.0*mu_F_count/Ein_F_count << endl;
  cout << "average joint_dist_F::Eout_F_count: " << 1.0*Eout_F_count/mu_F_count << endl;
}
// ----------- joint_dist::setup_param ------------------
// Initializes the quadrature parameters
void joint_dist::setup_param( joint_dist_param *Ein_param )
{
  static double skip_tol = Global.Value( "abs_tol" );

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
      FatalError( "joint_dist::setup_param", "energies inconsistent" );
    }
  }
}
// ----------- joint_dist::get_Ein_range --------------
//  Gets the range of nontrivial incident energy bins; computes first_Ein and last_Ein
// returns true if the threshold is too high for the energy bins
bool joint_dist::get_Ein_range( const dd_vector& sigma, const dd_vector& mult,
    const dd_vector& weight,
    const Flux_List& e_flux, const Energy_groups& Ein_groups )
{
  double E_first;
  double E_last;
  joint_dist_param initial_param;
  bool done = initial_param.get_Ein_range( sigma, mult, weight, e_flux,
                                         Ein_groups, &E_first, &E_last );
  if( done ) return true;

  // check the range of incident energies for the probability data
  joint_dist::const_iterator Ein_data_ptr = begin( );
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
// -----------  joint_dist::set_Ein_range ------------------
// Sets the range of incident energies for this intergration
void joint_dist::set_Ein_range( joint_dist_param *Ein_param )
{
  Ein_param->set_Ein_range( );
  double this_E = Ein_param->this_Ein_ptr->get_E_in( );
  if( this_E > Ein_param->data_E_0 ) Ein_param->data_E_0 = this_E;
  this_E = Ein_param->next_Ein_ptr->get_E_in( );
  if( this_E < Ein_param->data_E_1 ) Ein_param->data_E_1 = this_E;

  if( Ein_param->data_E_1 < Ein_param->data_E_0 )
  {
    FatalError( "joint_dist_Ein_range", "check the incident energies" );
  }
}
// ----------- joint_dist::one_box ------------------
// Integrate over one E-E' box
void joint_dist::one_box( T_matrix& transfer, int Eout_count,
   joint_dist_param *Ein_param )
{
  // pointers to the hit_lists
  joint_dist_hits::iterator mu0_low_hit_ptr = Ein_param->lower_mu0_hits.begin( );
  joint_dist_hits::iterator next_mu0_low_ptr = mu0_low_hit_ptr;
  ++next_mu0_low_ptr;
  joint_dist_hits::iterator mu1_low_hit_ptr = Ein_param->lower_mu1_hits.begin( );
  joint_dist_hits::iterator next_mu1_low_ptr = mu1_low_hit_ptr;
  ++next_mu1_low_ptr;
  joint_dist_hits::iterator mu0_high_hit_ptr = Ein_param->upper_mu0_hits.begin( );
  joint_dist_hits::iterator next_mu0_high_ptr = mu0_high_hit_ptr;
  ++next_mu0_high_ptr;
  joint_dist_hits::iterator mu1_high_hit_ptr = Ein_param->upper_mu1_hits.begin( );
  joint_dist_hits::iterator next_mu1_high_ptr = mu1_high_hit_ptr;
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
    if( ( ( mu0_low_hit_ptr->hit_edge == ABOVE ) ||
          ( mu0_low_hit_ptr->hit_edge == TOP_OUT ) ) &&
	( ( mu1_low_hit_ptr->hit_edge == ABOVE ) ||
          ( mu1_low_hit_ptr->hit_edge == TOP_OUT ) ) )
    {
      // do nothing---we are above the E-E' box
    }
    else if( ( ( mu0_high_hit_ptr->hit_edge == BOTTOM_OUT ) ||
               (  mu0_high_hit_ptr->hit_edge == BELOW ) ) &&
	     ( ( mu1_high_hit_ptr->hit_edge == BOTTOM_OUT ) ||
               (  mu1_high_hit_ptr->hit_edge == BELOW ) ) )
    {
      // do nothing---we are below the E-E' box
    }
    else
    {
      update_T( transfer, Eout_count, lower_Ein, upper_Ein, Ein_param );
    }
    // update the energy range
    lower_Ein = upper_Ein;
    static double e_tol = Global.Value( "E_tol" );
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
// ----------- joint_dist::update_T --------------
// Increments the transfer matrix
void joint_dist::update_T( T_matrix &transfer, int Eout_count,
       double Ein0_orig, double Ein1_orig, joint_dist_param *Ein_param )
{
  static double tol = Global.Value( "quad_tol" );
  // differences of nearly-equal numbers can cause problems; when to skip an interval
  static double from_quad_tol = Global.Value( "abs_quad_tol" )/100;
  static double from_abs_tol = 1000*Global.Value( "abs_tol" );
  double skip_tol = ( from_quad_tol > from_abs_tol ) ? from_quad_tol : from_abs_tol;
  // a vector to store the integrals, one Legendre order
  coef_vector value( transfer.order, transfer.conserve );
  value.set_zero( );
  // parameters for the integration
  QuadParamBase *params = static_cast< QuadParamBase* >( Ein_param );

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
      quad_F::integrate( joint_dist_F::Ein_F, transfer.Ein_quad_method,
                         use_Ein0, use_Ein1, params, tol, &value );
    }
    // add this integral
    transfer( Ein_param->Ein_count, Eout_count ) += value;
    // increment the function counts
    Ein_param->Ein_F_count += Ein_param->func_count;
    ++Ein_param->quad_count;
  }
}
// ----------- joint_dist::next_Ein_pair --------------
// Go to the next pair of incident energies.  Returns "true" when finished.
bool joint_dist::next_Ein_pair( double E_in, joint_dist_param *Ein_param )
{
  bool done = Ein_param->update_bin_pointers( E_in );
  if( !done )
  {
    static double etol = Global.Value( "E_tol" );
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
// ----------- joint_dist::mu_data_ladder --------------
// Handles the ( cosine, Eout, probability ) data for one pair of incident energies
void joint_dist::mu_data_ladder( T_matrix& transfer, joint_dist_param *Ein_param )
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
// ----------- joint_dist::start_mu_data------------
// Initializes the pointers to the ( Eout, probability ) data for current cosines
void joint_dist::start_mu_data( joint_dist_param *Ein_param )
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
// ----------- joint_dist::next_mu_pairs --------------
// Go to the next pairs of direction cosines.  Returns "true" when finished.
bool joint_dist::next_mu_pairs( joint_dist_param *Ein_param )
{
  bool done = false;
  // Ein_param->lower_mu and Ein_param->upper_mu define the range of interpolation in unit-base mu
  Ein_param->lower_mu = Ein_param->upper_mu;  // the new lower limit for the integral over mu
  static double mu_tol = Global.Value( "E_tol" );
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
// ----------- joint_dist::Eout_data_ladder --------------
// Adds to the transfer matrix for all E_out bins for given Ein and mu ranges
void joint_dist::Eout_data_ladder( T_matrix& transfer, joint_dist_param *Ein_param )
{
  double dummy = 0.0;  // needed for Ein_param->lower_mu0_hits.hit_box
  start_Eout_data( Ein_param );
  // loop through the energy data
  for( ; ; )
  {
    // loop through the outgoing energy bins
    for( int Eout_count = 0; Eout_count < transfer.num_Eout_bins;
      ++Eout_count )
    {
      vector< double >::const_iterator Eout_ptr = transfer.out_groups.begin( )
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
    bool done = next_Eout( Ein_param );
    if( done ) break;
  }
}
// ----------- joint_dist::start_Eout_data --------------
// Starts one staircase of the Eout data
void joint_dist::start_Eout_data( joint_dist_param *Ein_param )
{
  // This code is OK for unit-base data
  Ein_param->lower_Eout = 0.0;

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

  // For 3-d interpolation set up the E_mu_P_data at given incident energies for
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
// ----------- joint_dist::next_Eout --------------
// go to the next sets of (E_out, probability) pairs
bool joint_dist::next_Eout( joint_dist_param *Ein_param )
{
  bool done = false;
  double this_Eout;

  // the new lower limit for the integral over Eout is the former upper limit
  Ein_param->lower_Eout = Ein_param->upper_Eout;
  Ein_param->Ein0_data.mu0_Eout0.copy( Ein_param->Ein0_data.mu0_Eout1 );
  Ein_param->Ein0_data.mu1_Eout0.copy( Ein_param->Ein0_data.mu1_Eout1 );
  Ein_param->Ein1_data.mu0_Eout0.copy( Ein_param->Ein1_data.mu0_Eout1 );
  Ein_param->Ein1_data.mu1_Eout0.copy( Ein_param->Ein1_data.mu1_Eout1 );

  static double e_tol = Global.Value( "E_tol" );
  if( Ein_param->Ein0_mu0_Eout1->x <= Ein_param->lower_Eout*(1 + e_tol ) )
  {
    Ein_param->Ein0_mu0_Eout0 = Ein_param->Ein0_mu0_Eout1;
    ++Ein_param->Ein0_mu0_Eout1;
    if( Ein_param->Ein0_mu0_Eout1 == Ein_param->this_Ein_this_mu->end( ) )
    {
      return true;
    }
  }
  Ein_param->upper_Eout = Ein_param->Ein0_mu0_Eout1->x;

  if( Ein_param->Ein0_mu1_Eout1->x <= Ein_param->lower_Eout*(1 + e_tol ) )
  {
    Ein_param->Ein0_mu1_Eout0 = Ein_param->Ein0_mu1_Eout1;
    ++Ein_param->Ein0_mu1_Eout1;
    if( Ein_param->Ein0_mu1_Eout1 == Ein_param->this_Ein_next_mu->end( ) )
    {
      return true;
    }
  }
  this_Eout = Ein_param->Ein0_mu1_Eout1->x;
  if( Ein_param->upper_Eout > this_Eout )
  {
    Ein_param->upper_Eout = this_Eout;
  }

  if( Ein_param->Ein1_mu0_Eout1->x <= Ein_param->lower_Eout*(1 + e_tol ) )
  {
    Ein_param->Ein1_mu0_Eout0 = Ein_param->Ein1_mu0_Eout1;
    ++Ein_param->Ein1_mu0_Eout1;
    if( Ein_param->Ein1_mu0_Eout1 == Ein_param->next_Ein_this_mu->end( ) )
    {
      return true;
    }
  }
  this_Eout = Ein_param->Ein1_mu0_Eout1->x;
  if( Ein_param->upper_Eout > this_Eout )
  {
    Ein_param->upper_Eout = this_Eout;
  }


  if( Ein_param->Ein1_mu1_Eout1->x <= Ein_param->lower_Eout*(1 + e_tol ) )
  {
    Ein_param->Ein1_mu1_Eout0 = Ein_param->Ein1_mu1_Eout1;
    ++Ein_param->Ein1_mu1_Eout1;
    if( Ein_param->Ein1_mu1_Eout1 == Ein_param->next_Ein_next_mu->end( ) )
    {
      return true;
    }
  }
  this_Eout = Ein_param->Ein1_mu1_Eout1->x;
  if( Ein_param->upper_Eout > this_Eout )
  {
    Ein_param->upper_Eout = this_Eout;
  }

  // reset the E_mu_P_data for the higher Eout
  Ein_param->common_mu_Eout( Ein_param->upper_Eout,
                  &Ein_param->Ein0_data.mu0_Eout1,
		  &Ein_param->Ein0_data.mu1_Eout1, &Ein_param->Ein1_data.mu0_Eout1,
		  &Ein_param->Ein1_data.mu1_Eout1 );

  // Reset the parameters for the intersection of the physical enegies with the quadrature box
  Ein_param->lower_mu0_hits.get_phys_Eout( );
  Ein_param->lower_mu1_hits.get_phys_Eout( );
  Ein_param->upper_mu0_hits.get_phys_Eout( );
  Ein_param->upper_mu1_hits.get_phys_Eout( );

  return done;
}

// **************** Functions to integrate *********************
// ---------------- joint_dist_F::Eout_F ------------------
// Function for the 1-d quadrature over outgoing energy
void joint_dist_F::Eout_F( double Eout, QuadParamBase *Eout_quad_param,
			   coef_vector *value )
{
  // the parameters are really joint_Eout_param
  joint_Eout_param *Eout_params = static_cast< joint_Eout_param* >( Eout_quad_param );
  Eout_params->func_count += 1;

  // for interpolation
  E_mu_P_data mid_data;
  if( Eout_params->Eout_interp == LINLIN )
  {
    Eout_params->Eout0_data.UB_Eout_interp( Eout, Eout_params->Eout1_data,
					     &mid_data);
  }
  else
  {
    Eout_params->Eout0_data.UB_Eout_histogram( Eout, Eout_params->Eout1_data,
					     &mid_data);
  }
  value->weight_1[ 0 ] = mid_data.Prob;
  value->weight_E[ 0 ] = mid_data.Prob*mid_data.phys_Eout;
}
// ---------------- joint_dist_F::mu_F ------------------
// Function for the 1-d quadrature over outgoing energy
void joint_dist_F::mu_F( double mu, QuadParamBase *mu_quad_param,
			   coef_vector *value )
{
  // the parameters are really joint_mu_param
  joint_mu_param *mu_params = static_cast< joint_mu_param* >( mu_quad_param );
  mu_params->func_count += 1;

  // interpolate data
  unit_base_map mid_ubase_map;
  mu_params->this_data.mu_interpolate( mu, &mu_params->Eout_params.Eout0_data,
					&mu_params->Eout_params.Eout1_data,
					&mid_ubase_map );

  // the range of integration
  double Eout0 = ( mu_params->Eout_params.Eout0_data.phys_Eout < *mu_params->Eout_bottom )?
    mid_ubase_map.to_unit_base( *mu_params->Eout_bottom ) :
    mu_params->Eout_params.Eout0_data.UB_Eout;

  double Eout1 = ( mu_params->Eout_params.Eout1_data.phys_Eout > *mu_params->Eout_top )?
    mid_ubase_map.to_unit_base( *mu_params->Eout_top ) :
    mu_params->Eout_params.Eout1_data.UB_Eout;

  // for integration over outgoing energy
  QuadParamBase *params = static_cast< QuadParamBase* >( &mu_params->Eout_params );
  // the mu dependence is scalar, but we might need the weight by ougtoing energy
  coef_vector Eout_integral( 0, BOTH );
  static double tol = Global.Value( "quad_tol" );
  // 2nd-order Gaussian quadrature is exact
  quad_F::integrate( joint_dist_F::Eout_F, GAUSS2, Eout0, Eout1,
		     params, tol, &Eout_integral );

  // the Legendre polynomials
  math_F::Legendre( mu_params->Eout_params.Eout0_data.phys_mu, value );
  // scale by the integrals over outgoing energy
  double Prob = Eout_integral.weight_1[ 0 ];
  *value *= Prob;
  // do the energy weighting if necessary
  if( ( Prob > 0.0 ) &&
      ( ( value->conserve == ENERGY ) || ( value->conserve == BOTH ) ) )
  {
    value->scale_E( Eout_integral.weight_E[ 0 ]/Prob );
  }
  mu_params->Eout_F_count += mu_params->Eout_params.func_count;
}
// ----------- joint_dist_F::Ein_F ------------------
// Function for the 3-d quadrature over incident energy, cosine, and outgoing energy
void joint_dist_F::Ein_F( double E_in, QuadParamBase *Ein_quad_param,
  coef_vector *value )
{
  // the parameters are really joint_dist_param
  joint_dist_param *Ein_params = static_cast< joint_dist_param* >( Ein_quad_param );
  Ein_params->func_count += 1;

  // Interpolate the data in incident energy
  Ein_params->Ein0_data.Ein_interpolate( E_in, Ein_params->Ein1_data,
					 &Ein_params->mu_params.this_data );

  // parameters for the integration over mu
  QuadParamBase *params = static_cast< QuadParamBase* >( &Ein_params->mu_params );
  static double tol = Global.Value( "quad_tol" );

  // We are integrating a quartic in mu, so 4-th order Gauss is exact
  Quadrature_Method mu_quad_method = GAUSS4;
  //  Quadrature_Method mu_quad_method = ADAPTIVE2;

  // to hold the value of the integral over a mu subinterval
  coef_vector temp_value( value->order, value->conserve );

  // Set the value to zero in case the geometry is bad
  value->set_zero( );

  // Determine how the data meets the outgoing physical energy range
  if( Ein_params->geometry( Ein_params->mu_params.this_data ) )
  {
    energy_hit_list::iterator low_hit_ptr = Ein_params->lower_2d_hits.begin( );
    energy_hit_list::iterator next_low_ptr = low_hit_ptr;
    ++next_low_ptr;
    energy_hit_list::iterator high_hit_ptr = Ein_params->upper_2d_hits.begin( );
    energy_hit_list::iterator next_high_ptr = high_hit_ptr;
    ++next_high_ptr;
    for( ; ( next_low_ptr != Ein_params->lower_2d_hits.end( ) ) &&
           ( next_high_ptr != Ein_params->upper_2d_hits.end( ) );
         low_hit_ptr = next_low_ptr, ++next_low_ptr,
           high_hit_ptr = next_high_ptr, ++next_high_ptr )
    {
      if( ( low_hit_ptr->hit_edge == ABOVE ) ||
          ( low_hit_ptr->hit_edge == TOP_OUT ) )
      {
        // do nothing---we are above the E-E' box
        continue;
      }
      else if( ( high_hit_ptr->hit_edge == BOTTOM_OUT ) ||
               ( high_hit_ptr->hit_edge == BELOW ) )
      {
        // do nothing---we are below the E-E' box
        continue;
      }
      else
      {
        quad_F::integrate( joint_dist_F::mu_F, mu_quad_method,
                         low_hit_ptr->E_in, next_low_ptr->E_in, 
                         params, tol, &temp_value );
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
}
