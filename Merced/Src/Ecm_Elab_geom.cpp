/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2010-11-09  (Tue., Sept. 15, 2009) $
 * $Author: hedstrom $
 * $Id: Ecm_Elab_geom.cpp 1 2010-11-09 hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
//! Implements the classes used for the geometry of the map from center of mass to lab frame

#include <cmath>

#include "Ecm_Elab_geom.hpp"
#include "math_util.hpp"
#include "messaging.hpp"
#include "global_params.hpp"

// ************* class Ecm_Elab_mu_param *****************
// ----------- Ecm_Elab_mu_param::setup --------------
// Sets up the data for this incident energy
void Ecm_Elab_mu_param::setup( double Ein, double Eoutcm, const Ecm_Elab_Ecm_param& Ecm_param )
{
  map = Ecm_param.map;
  E_in = Ein;
  Eout_cm = Eoutcm;

  // use the geometry to set the range of integration over center-of-mass cosine
  double E_trans = map->get_Etrans( E_in );
  double V_trans = sqrt( E_trans );
  double V_cm = sqrt( Eout_cm );
  double V_bottom = sqrt( Ecm_param.lab_Eout_min );
  double V_top = sqrt( Ecm_param.lab_Eout_max );
  if( V_trans < V_bottom )
  {
    if( ( V_trans + V_cm < V_bottom ) || ( V_cm > V_trans + V_top ) )
    {
      FatalError( "Ecm_Elab_mu_param::setup", "energy error for small Ein" );
    }
    if( V_trans + V_cm <= V_top )
    {
      mu_cm_max = 1.0;
    }
    else
    {
      mu_cm_max = map->get_mu_cm( Ein, Ecm_param.lab_Eout_max, Eout_cm );
    }
    if( V_cm <= V_trans + V_bottom )
    {
      mu_cm_min =  map->get_mu_cm( Ein, Ecm_param.lab_Eout_min, Eout_cm );
    }
    else if( V_cm <= V_trans + V_top )
    {
      mu_cm_min = -1.0;
    }
    else
    {
      mu_cm_min =  map->get_mu_cm( Ein, Ecm_param.lab_Eout_max, Eout_cm );
    }
  }
  else if( V_trans > V_top )
  {
    if( ( V_trans > V_top + V_cm ) || ( V_cm > V_trans + V_top ) )
    {
      FatalError( "Ecm_Elab_mu_param::setup", "energy error for large Ein" );
    }
    if( ( V_cm < V_trans + V_bottom ) && ( V_trans < V_cm + V_bottom ) )
    {
      mu_cm_min = map->get_mu_cm( Ein, Ecm_param.lab_Eout_min, Eout_cm );
    }
    else
    {
      mu_cm_min = -1.0;
    }
    mu_cm_max =  map->get_mu_cm( Ein, Ecm_param.lab_Eout_max, Eout_cm );
  }
  else
  {
    if( V_cm > V_trans + V_top )
    {
      FatalError( "Ecm_Elab_mu_param::setup", "energy error for intermediate Ein" );
    }
    if( ( V_cm < V_trans + V_bottom ) && ( V_trans < V_cm + V_bottom ) )
    {
      mu_cm_min = map->get_mu_cm( Ein, Ecm_param.lab_Eout_min, Eout_cm );
    }
    else
    {
      mu_cm_min = -1.0;
    }
    if( V_top < V_trans + V_cm )
    {
      mu_cm_max = map->get_mu_cm( Ein, Ecm_param.lab_Eout_max, Eout_cm );
    }
    else
    {
      mu_cm_max = 1.0;
    }
  }
}

// ************* class Ecm_Elab_Ecm_param *****************
// ----------- Ecm_Elab_Ecm_param::setup --------------
// Sets up the data for this incident energy
void Ecm_Elab_Ecm_param::setup( double Ein, double Eoutmin, double Eoutmax,
  double data_Ecmmin, double data_Ecmmax )
{
  E_in = Ein;
  lab_Eout_min = Eoutmin;  // laboratory energy range
  lab_Eout_max = Eoutmax;
  V_lab_max = sqrt( lab_Eout_max );
  V_lab_min = sqrt( lab_Eout_min );
  data_Ecm_min = data_Ecmmin;   //. center of mass energy range for the data
  data_Ecm_max = data_Ecmmax;
}
// ----------- Ecm_Elab_Ecm_param::V_lab_sectors --------------
// Identifies the regions of integration over Eout_lab and mu_cm
void Ecm_Elab_Ecm_param::V_lab_sectors( )
{
  if( V_cm_limits.size( ) > 0 )
  {
    V_cm_limits.erase( V_cm_limits.begin( ), V_cm_limits.end( ) );
  }
  Ecm_intersect one_intersection;
  one_intersection.gamma = map->gamma;
  one_intersection.set_energies( E_in, data_Ecm_min );
  min_V_cm = one_intersection.V_cm;
  one_intersection.set_energies( E_in, data_Ecm_max );
  max_V_cm = one_intersection.V_cm;

  Vcm_quadBox_Hit Vcm_quadBox_hit;
  if( V_lab_min > 0.0 )
  {
    Vcm_quadBox_hit.V_cm = abs( one_intersection.V_trans - V_lab_min );  // mu = 1 for V_lab_min
    Vcm_quadBox_hit.hit_corner = BOTTOM_FORWARD;
    V_cm_limits.push_back( Vcm_quadBox_hit );
    Vcm_quadBox_hit.V_cm = one_intersection.V_trans + V_lab_min;  // mu = -1 for V_lab_min
    Vcm_quadBox_hit.hit_corner = BOTTOM_BACKWARD;
    V_cm_limits.push_back( Vcm_quadBox_hit );
  }
  Vcm_quadBox_hit.V_cm = abs( one_intersection.V_trans - V_lab_max );  // mu = 1 for V_lab_max
  Vcm_quadBox_hit.hit_corner = TOP_FORWARD;
  V_cm_limits.push_back( Vcm_quadBox_hit );
  Vcm_quadBox_hit.V_cm = one_intersection.V_trans + V_lab_max;  // mu = -1 for V_lab_max
  Vcm_quadBox_hit.hit_corner = TOP_BACKWARD;
  V_cm_limits.push_back( Vcm_quadBox_hit );

  // add the values of V_cm given by the data
  Vcm_quadBox_hit.V_cm = min_V_cm;  // the lower V_cm
  Vcm_quadBox_hit.set_region( min_V_cm, one_intersection.V_trans, V_lab_min, V_lab_max );
  V_cm_limits.push_front( Vcm_quadBox_hit );

  Vcm_quadBox_hit.V_cm = max_V_cm;  // the higher V_cm
  Vcm_quadBox_hit.set_region( max_V_cm, one_intersection.V_trans, V_lab_min, V_lab_max );
  V_cm_limits.push_front( Vcm_quadBox_hit );

  V_cm_limits.sort( Vcm_quadBox_Hit_F::lessthan_F );
}
// ----------- Ecm_Elab_Ecm_param::Ecm_range --------------
// Computes the range of center-of-mass outgoing energies
// Returns the tolerance to use in the quadrature over center-of-mass energy and cosine
double Ecm_Elab_Ecm_param::Ecm_range( )
{
  double E_trans = map->get_Etrans( E_in );
  double V_trans = sqrt( E_trans );
  double V_bottom = sqrt( lab_Eout_min );
  double V_top = sqrt( lab_Eout_max );
  Ecm_min = get_Ecm( min_hit_corner, V_trans, V_bottom, V_top, data_Ecm_min );
  Ecm_max = get_Ecm( max_hit_corner, V_trans, V_bottom, V_top, data_Ecm_max );
  if( Ecm_min > Ecm_max )
  {
    Warning( "Ecm_Elab_Ecm_param::Ecm_range", "E_cm out of order" );
  }
  static double epsilon = Global.Value( "quad_tol" );
  if( E_trans > lab_Eout_max )
  {
    double energy_ratio = epsilon*E_trans/lab_Eout_max;
    return energy_ratio/( 1.0 + energy_ratio );
  }
  else
  {
    return epsilon;
  }
}
// ----------- Ecm_Elab_Ecm_param::get_Ecm --------------
// Computes the center-of-mass outgoing energy for a given Hit_Corner
double Ecm_Elab_Ecm_param::get_Ecm( Hit_Corner hit_corner, double V_trans,
  double V_bottom, double V_top, double data_Ecm )
{
  double V_cm = 0.0;
  switch( hit_corner )
  {
  case V_INSIDE:
    return data_Ecm;
  case BOTTOM_FORWARD:
    V_cm = abs( V_bottom - V_trans );
    break;
  case BOTTOM_BACKWARD:
    V_cm = V_bottom + V_trans;
    break;
  case TOP_FORWARD:
    V_cm = abs( V_top - V_trans );
    break;
  case TOP_BACKWARD:
    V_cm = V_top + V_trans;
    break;
  default:
    FatalError( "Ecm_Elab_Ecm_param::get_Ecm", "bad Hit_Corner" );
  }

  double E_cm = V_cm*V_cm;
  if( ( E_cm < data_Ecm_min ) || ( E_cm > data_Ecm_max ) )
  {
    Warning( "Ecm_Elab_Ecm_param::get_Ecm", "E_cm out of data interval" );
  }
  return E_cm;
}

// ************* class Ecm_Elab_Ein_param *****************
