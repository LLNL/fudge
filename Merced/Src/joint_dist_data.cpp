/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 601 $
 * $Date: 2017-12-18 $
 * $Author: hedstrom $
 * $Id: joint_dist_data.hpp 601 2017-12-18Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// header for the classes for data used on joint energy-angle distributions

#include "joint_dist_data.hpp"
#include "messaging.hpp"
#include "global_params.hpp"

// ********* class Jdata::E_mu_P_data *+********
// ----------- Jdata::E_mu_P_data::UB_Eout_interp --------------
// Interpolates with respect to unit-base outgoing energy
bool Jdata::E_mu_P_data::UB_Eout_interp( double mid_UB_Eout,
   const Jdata::E_mu_P_data &next_data,
		       Jdata::E_mu_P_data *mid_data ) const
{
  // check the inputs
  bool E_OK = true;
  double E_diff = next_data.UB_Eout - UB_Eout;
  if( E_diff <= 0.0 )
  {
    Msg::DebugInfo( "Jdata::E_mu_P_data::UB_Eout_interp", "data out of order" );
    E_OK = false;
  }
  double alpha = ( mid_UB_Eout - UB_Eout )/E_diff;
  static double etol = Global.Value( "looser_tol" );

  if( alpha < 0.0 )
  {
    Msg::DebugInfo( "Jdata::E_mu_P_data::UB_Eout_interp", "extrapolation down" );
    alpha = 0.0;
  }
  else if( alpha > 1.0+etol )
  {
    Msg::DebugInfo( "Jdata::E_mu_P_data::UB_Eout_interp", "extrapolation up" );
    alpha = 1.0;
  }
  mid_data->UB_mu = UB_mu;
  mid_data->phys_mu = phys_mu;
  mid_data->UB_Eout = mid_UB_Eout;
  mid_data->phys_Eout = ( 1.0 - alpha )*phys_Eout +
    alpha*next_data.phys_Eout;
  mid_data->mu_Prob = ( 1.0 - alpha )*mu_Prob + alpha*next_data.mu_Prob;
  mid_data->Eout_Prob = ( 1.0 - alpha )*Eout_Prob + alpha*next_data.Eout_Prob;

  return( E_OK );
}
// ----------- Jdata::E_mu_P_data::UB_Eout_histogram --------------
// Interpolates histogram data with respect to unit-base outgoing energy
void Jdata::E_mu_P_data::UB_Eout_histogram( double mid_UB_Eout,
    const Jdata::E_mu_P_data &next_data,
		       Jdata::E_mu_P_data *mid_data ) const
{
  // check the inputs
  double E_diff = next_data.UB_Eout - UB_Eout;
  if( E_diff <= 0.0 )
  {
    Msg::FatalError( "Jdata::E_mu_P_data::UB_Eout_histogram",
		     "data out of order" );
  }
  double alpha = ( mid_UB_Eout - UB_Eout )/E_diff;
  static double etol = Global.Value( "looser_tol" );
  static bool no_warning = true;
  if( ( ( alpha < 0.0 ) || ( alpha > 1.0 ) ) &&
      ( E_diff > etol ) && no_warning )
  {
    Msg::Warning( "Jdata::E_mu_P_data::UB_Eout_histogram", "extrapolation" );
    no_warning = false;
  }
  mid_data->UB_mu = UB_mu;
  mid_data->phys_mu = phys_mu;
  mid_data->UB_Eout = mid_UB_Eout;
  mid_data->phys_Eout = ( 1.0 - alpha )*phys_Eout +
    alpha*next_data.phys_Eout;
  mid_data->mu_Prob = mu_Prob;
  mid_data->Eout_Prob = Eout_Prob;
}
// ----------- Jdata::E_mu_P_data::mu_interp --------------
// Interpolates with respect to direction cosine
bool Jdata::E_mu_P_data::mu_interp( double mid_mu, const Jdata::E_mu_P_data &next_data,
		       Jdata::E_mu_P_data *mid_data ) const
{
  // check the inputs
  bool mu_OK = true;
  double mu_diff = next_data.UB_mu - UB_mu;
  if( mu_diff <= 0.0 )
  {
    Msg::DebugInfo( "Jdata::E_mu_P_data::mu_interp", "data out of order" );
    mu_OK = false;
  }
  double alpha = ( mid_mu - UB_mu )/mu_diff;
  static const double mu_tol = Global.Value( "looser_tol" );
  static bool no_warning = true;
  if( ( ( alpha <= -mu_tol ) || ( alpha >= 1.0 + mu_tol ) ) && no_warning )
  {
    Msg::DebugInfo( "Jdata::E_mu_P_data::mu_interp", "extrapolation" );
    no_warning = false;
  }
  mid_data->UB_mu = mid_mu;
  mid_data->phys_mu = ( 1.0 - alpha )*phys_mu + alpha*next_data.phys_mu;
  mid_data->UB_Eout = UB_Eout;
  mid_data->phys_Eout = ( 1.0 - alpha )*phys_Eout +
    alpha*next_data.phys_Eout;
  mid_data->mu_Prob = ( 1.0 - alpha )*mu_Prob + alpha*next_data.mu_Prob;
  mid_data->Eout_Prob = ( 1.0 - alpha )*Eout_Prob + alpha*next_data.Eout_Prob;

  return( mu_OK && no_warning );
}
// ----------- Jdata::E_mu_P_data::Ein_interp --------------
// Interpolates with respect to unit-base outgoing energy
bool Jdata::E_mu_P_data::Ein_interp( double alpha, const Jdata::E_mu_P_data &next_data,
		       Jdata::E_mu_P_data *mid_data ) const
{
  // check the inputs
  if( next_data.UB_mu != UB_mu )
  {
    Msg::FatalError( "Jdata::E_mu_P_data::Ein_interp",
		"interpolation valid only for idential unit-base mu values" );
  }

  static bool no_warning = true;
  if( ( ( alpha < 0.0 ) || ( alpha > 1.0 ) ) && no_warning )
  {
    Msg::DebugInfo( "Jdata::E_mu_P_data::Ein_interp", "extrapolation" );
    no_warning = false;
  }
  //  mid_data->UB_mu = ( 1.0 - alpha )*UB_mu + alpha*next_data.UB_mu;
  mid_data->UB_mu = UB_mu;
  mid_data->phys_mu = ( 1.0 - alpha )*phys_mu + alpha*next_data.phys_mu;
  mid_data->UB_Eout = UB_Eout;
  mid_data->phys_Eout = ( 1.0 - alpha )*phys_Eout +
    alpha*next_data.phys_Eout;
  mid_data->mu_Prob = ( 1.0 - alpha )*mu_Prob + alpha*next_data.mu_Prob;
  mid_data->Eout_Prob = ( 1.0 - alpha )*Eout_Prob + alpha*next_data.Eout_Prob;

  return no_warning;
}
// ----------- Jdata::E_mu_P_data::copy --------------
// Copies the data
void Jdata::E_mu_P_data::copy( const Jdata::E_mu_P_data &to_copy )
{
  UB_mu = to_copy.UB_mu;
  phys_mu = to_copy.phys_mu;
  UB_Eout = to_copy.UB_Eout;
  phys_Eout = to_copy.phys_Eout;
  mu_Prob = to_copy.mu_Prob;
  Eout_Prob = to_copy.Eout_Prob;
}

// ********* class Jdata::current_data *********
// ----------- Jdata::current_data::Ein_interpolate --------------
// Interpolate in incident energy, returns mid_data
bool Jdata::current_data::Ein_interpolate( double mid_Ein,
                    const Jdata::current_data &next_data,
                    Jdata::current_data *mid_data ) const
{
  // check the inputs
  bool is_OK = true;
  double E_diff = next_data.get_E_in( ) - get_E_in( );
  if( E_diff <= 0.0 )
  {
    Msg::DebugInfo( "Jdata::current_data::Ein_interpolate", "data out of order" );
    return false;
  }
  double alpha = ( mid_Ein - get_E_in( ) )/E_diff;
  static const double e_tol = Global.Value( "looser_tol" );

  if( ( alpha < 0.0 ) || ( alpha > 1.0+e_tol ) )
  {
    Msg::DebugInfo( "Jdata::current_data::Ein_interpolate", "extrapolation" );
    return false;
  }
  mid_data->set_E_in( mid_Ein );

  is_OK = mu0_Eout0.Ein_interp( alpha, next_data.mu0_Eout0, &( mid_data->mu0_Eout0 ) );
  if( !is_OK ) return false;
  is_OK = mu0_Eout1.Ein_interp( alpha, next_data.mu0_Eout1, &( mid_data->mu0_Eout1 ) );
  if( !is_OK ) return false;
  is_OK = mu1_Eout0.Ein_interp( alpha, next_data.mu1_Eout0, &( mid_data->mu1_Eout0 ) );
  if( !is_OK ) return false;
  is_OK = mu1_Eout1.Ein_interp( alpha, next_data.mu1_Eout1, &( mid_data->mu1_Eout1 ) );
  if( !is_OK ) return false;

  // interpolate the unit-base map ranges of outgoing energy
  is_OK = mid_data->mu0_ubase_map.interpolate( alpha, mu0_ubase_map,
				       next_data.mu0_ubase_map );
  if( !is_OK ) return false;
  is_OK = mid_data->mu1_ubase_map.interpolate( alpha, mu1_ubase_map,
				       next_data.mu1_ubase_map );
  if( !is_OK ) return false;

  return true;
}
// ----------- Jdata::current_data::mu_interpolate --------------
// Interpolate in mu, returns Eout0_data, Eout1_data, and mid_ubase_map
bool Jdata::current_data::mu_interpolate( double mid_mu, Jdata::E_mu_P_data *Eout0_data,
                    Jdata::E_mu_P_data *Eout1_data,
				   Ddvec::unit_base_map *mid_ubase_map ) const
{
  bool is_OK = mu0_Eout0.mu_interp( mid_mu, mu1_Eout0, Eout0_data );
  if( !is_OK ) return false;
  is_OK = mu0_Eout1.mu_interp( mid_mu, mu1_Eout1, Eout1_data );
  if( !is_OK ) return false;

  // interpolate the unit-base map
  double alpha = ( mid_mu - mu0_Eout0.UB_mu )/( mu1_Eout0.UB_mu - mu0_Eout0.UB_mu );
  is_OK = mid_ubase_map->interpolate( alpha, mu0_ubase_map, mu1_ubase_map );
  if( !is_OK ) return false;

  return true;
}

