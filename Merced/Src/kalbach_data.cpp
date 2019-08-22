/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2009-09-15 (Tue, 15 Sep 2009) $
 * $Author: hedstrom $
 * $Id: kalbach_data.cpp 1 2009-09-15 hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// implementation of the classes used to handle Kalbach data

#include <cmath>

#include "kalbach_data.hpp"
#include "messaging.hpp"
#include "global_params.hpp"

// ************* class kalbach_data *****************
// ----------- kalbach_data::linlin_interp --------------
// Does linear interpolation of this data with the next
void kalbach_data::linlin_interp( double e_in, const kalbach_data& left_data,
  const kalbach_data& right_data )
{
  if( ( e_in < left_data.E_in ) || ( e_in > right_data.E_in ) )
  {
    Warning( "kalbach_data::linlin_interp", "extrapolation");
  }
  E_in = e_in;
  Eout_interp = left_data.Eout_interp;
  double denom = right_data.E_in - left_data.E_in;
  double alpha;
  if( denom == 0.0 )
  {
    this_Ecm = left_data.this_Ecm;
    this_f0 = left_data.this_f0;
    this_r = left_data.this_r;
    next_Ecm = left_data.next_Ecm;
    next_f0 = left_data.next_f0;
    next_r = left_data.next_r;
    alpha = 0.0;
  }
  else
  {
    alpha = ( e_in - left_data.E_in )/denom;
    this_Ecm = ( 1.0 - alpha )*left_data.this_Ecm + alpha*right_data.this_Ecm;
    this_f0 = ( 1.0 - alpha )*left_data.this_f0 + alpha*right_data.this_f0;
    next_Ecm = ( 1.0 - alpha )*left_data.next_Ecm + alpha*right_data.next_Ecm;
    next_f0 = ( 1.0 - alpha )*left_data.next_f0 + alpha*right_data.next_f0;
    this_r = ( 1.0 - alpha )*left_data.this_r + alpha*right_data.this_r;
    next_r = ( 1.0 - alpha )*left_data.next_r + alpha*right_data.next_r;
  }
}
// ----------- kalbach_data::unit_base_interp --------------
// Does (unit-base) linear interpolation of this data with the next
void kalbach_data::unit_base_interp( double e_in, const kalbach_data& left_data,
  const kalbach_data& right_data )
{
  linlin_interp( e_in, left_data, right_data );

  double denom = right_data.E_in - left_data.E_in;
  double alpha;
  if( denom == 0.0 )
  {
    alpha = 0.0;
  }
  else
  {
    alpha = ( e_in - left_data.E_in )/denom;
  }
  unit_base.interpolate( alpha, left_data.unit_base, right_data.unit_base );
}
// ----------- kalbach_data::un_unit_base --------------
// Undoes the unit-base map for one outgoing energy
double kalbach_data::un_unit_base( double E_unit )
{
  return unit_base.un_unit_base( E_unit );
}
// ----------- kalbach_data::un_unit_base --------------
// Undoes the unit-base map; used on interpolated data
void kalbach_data::un_unit_base( )
{
  double scale = unit_base.Eout_max - unit_base.Eout_min;
  if( scale <= 0.0 )
  {
    FatalError( "kalbach_data::un_unit_base", "Bad unit base scale" );
  }
  this_Ecm = unit_base.un_unit_base( this_Ecm );
  next_Ecm = unit_base.un_unit_base( next_Ecm );
  this_f0 /= scale;
  next_f0 /= scale;
}
// ----------- kalbach_data::get_f0_r --------------
// Calculates r and the center-of-mass outgoing energy density
void kalbach_data::get_f0_r( double Eoutcm, double *Ecm_prob, double *r ) const
{
  if( ( Eoutcm < this_Ecm ) || ( Eoutcm > next_Ecm ) )
  {
    FatalError( "kalbach_data::get_f0_r", "Eoutcm outside its range" );
  }
  if( Eout_interp == HISTOGRAM )
  {
    *Ecm_prob = this_f0;
    *r = this_r;
  }
  else if( Eout_interp == LINLIN )
  {
    double alpha = ( Eoutcm - this_Ecm )/( next_Ecm - this_Ecm );
    *Ecm_prob = ( 1.0 - alpha )*this_f0 + alpha*next_f0;
    *r = ( 1.0 - alpha )*this_r + alpha*next_r;
  }
  else
  {
    FatalError( "kalbach_data::get_f0_r", "interpolation not implemented" );
  }
}

// ************* class nucleon *****************
// ----------- nucleon::operator= --------------
//! copies the data
nucleon& nucleon::operator=( const nucleon& to_copy )
{
  mass = to_copy.mass;
  Z = to_copy.Z;
  A = to_copy.A;
  Kalbach_I = to_copy.Kalbach_I;
  Kalbach_M = to_copy.Kalbach_M;
  Kalbach_m = to_copy.Kalbach_m;
  return *this;
}
// ----------- nucleon::set_params --------------
// Sets Z, A, I, M, m
void nucleon::set_params( int ZA )
{
  A = ZA % 1000;
  Z = ( ZA - A )/1000;
  switch( ZA )
  {
  case 0:  // gamma
    Kalbach_M = 0;
    Kalbach_m = 0.0;
    Kalbach_I = 0.0;
    break;
  case 1:  // neutron
    Kalbach_M = 1;
    Kalbach_m = 0.5;
    Kalbach_I = 0.0;
    break;
  case 1001:  // proton
    Kalbach_M = 1;
    Kalbach_m = 1.0;
    Kalbach_I = 0.0;
    break;
  case 1002:  // deuteron
    Kalbach_M = 1;
    Kalbach_m = 1.0;
    Kalbach_I = 2.22;
    break;
  case 1003:  // triton
    Kalbach_M = 0;
    Kalbach_m = 1.0;
    Kalbach_I = 8.48;
    break;
  case 2003:  // helium-3
    Kalbach_M = 0;
    Kalbach_m = 1.0;
    Kalbach_I = 7.72;
    break;
  case 2004:  // alpha
    Kalbach_M = 0;
    Kalbach_m = 2.0;
    Kalbach_I = 28.3;
    break;
  default:   // ZA > 2004
    Kalbach_M = 0;
    Kalbach_m = 0.0;
    Kalbach_I = 0.0;
    break;
  }
}
// ----------- nucleon::get_Sa --------------
// Computes the Kalbach S_a function, where *this is the compound nucleus
double nucleon::get_Sa( const nucleon &in_or_out )
{
  double S;
  // replace by neutron for incident gammas
  double use_A = ( in_or_out.A == 0 ) ? 1 : in_or_out.A;
  double resid_A = A - use_A;
  double resid_Z = Z - in_or_out.Z;
  double n_excess = A - 2*Z;  // how many more neutrons than protons
  double resid_n_excess = n_excess - ( use_A - 2*in_or_out.Z );

  S = 15.68*use_A -
    28.07*( (1.0*n_excess*n_excess)/A - (1.0*resid_n_excess*resid_n_excess)/resid_A ) -
    18.56*( pow( 1.0*A, 2.0/3 ) - pow( 1.0*resid_A, 2.0/3 ) ) +
    33.22*( (1.0*n_excess*n_excess)/pow( 1.0*A, 4.0/3 ) -
	    (1.0*resid_n_excess*resid_n_excess)/pow( 1.0*resid_A, 4.0/3 ) ) -
    0.717*( (1.0*Z*Z)/pow( 1.0*A, 1.0/3 ) - (1.0*resid_Z*resid_Z)/pow( 1.0*resid_A, 1.0/3 ) ) +
    1.211*( (1.0*Z*Z)/A - (1.0*resid_Z*resid_Z)/resid_A ) - in_or_out.Kalbach_I;
  return S;
}
// ----------- nucleon::check_data --------------
// Checks for proper initialization
bool nucleon::check_data( )
{
  return ( A >= 0 ) && ( mass >= 0.0 );
}

// ************* class Kalbach_a *****************
// -----------  Kalbach_a::setup_params ------------------
// Sets up the mass ratios in map
void Kalbach_a::setup_params( )
{
  // verify that we have particle data
  if( !projectile.check_data( ) )
  {
    FatalError( "Kalbach_a::setup_params", "insufficient projectile data" );
  }
  if( projectile.Z == 3 )
  {
    FatalError( "Kalbach_a::setup_params", 
      "Kalbach M value is undefined for incident triton and helium-3" );
  }
  if( !target.check_data( ) )
  {
    FatalError( "Kalbach_a::setup_params", "insufficient target data" );
  }
  int Z = projectile.Z + target.Z;
  int A = projectile.A + target.A;
  compound.set_params( 1000*Z + A );
  if( !compound.check_data( ) )
  {
    FatalError( "Kalbach_a::setup_params", "mass of compound not given" );
  }
  if( eject.get_ZA( ) == projectile.get_ZA( ) )
  {
    eject = projectile;
    residual = target;
  }
  else
  {
    Z = compound.Z - eject.Z;
    A = compound.A - eject.A;
    residual.set_params( 1000*Z + A );
    if( !residual.check_data( ) )
    {
      FatalError( "Kalbach_a::setup_params", "mass of residual not given" );
    }
  }
  // Do the masses make sense?
  double massSlop = 1.0e-6;  // the masses should agree by at least this much
  double refMass = ( projectile.mass == 0 )? target.mass * massSlop:
    projectile.mass / 40;
  if( abs( projectile.mass + target.mass - compound.mass ) > refMass )
  {
    FatalError( "Kalbach_a::setup_params", "Check the mass of the compound" );
  }
  refMass = ( eject.mass == 0 )? compound.mass * massSlop :
    eject.mass / 40;
  if( abs( eject.mass + residual.mass - compound.mass ) >  eject.mass/40 )
  {
    FatalError( "Kalbach_a::setup_params", "Check the mass of the residual" );
  }

  set_Sab( );
  map->setup_params( projectile.mass, target.mass, eject.mass, residual.mass );
}
// ----------- Kalbach_a::set_Sab --------------
// Computes the S_a and S_b functions
void Kalbach_a::set_Sab( )
{
  // There is special coding for photonuclear reactions
  if( projectile.A == 0 )
  {
    nucleon neutron;
    neutron.set_params( 1 );
    projectile_S = compound.get_Sa( neutron );
  }
  else
  {
    projectile_S = compound.get_Sa( projectile );
  }
  if( ( projectile.A > 0 ) && ( projectile.get_ZA( ) == eject.get_ZA( ) ) )
  {
    eject_S = projectile_S;
  }
  else
  {
    eject_S = compound.get_Sa( eject );
  }
}
// ----------- Kalbach_a::get_a --------------
// Computes the Kalbach a function
double Kalbach_a::get_a( double E_in, double E_out )
{
  double slope_a;
  // total center-of-mass kinetic energy of incident particles
  double total_cm_Ein = map->incident_cm_KE( E_in );
  // total outgoing kinetic energy for 2-body with unknown Q value
  double total_cm_Eout = map->total_Eout( E_out );

  double Ecm_proj = total_cm_Ein + projectile_S;  // available energy
  const double E_t1 = 130.0;  // MeV
  double R_1 = ( Ecm_proj < E_t1 ) ? Ecm_proj : E_t1;
  double Ecm_eject = total_cm_Eout + eject_S;
  double X_1 = ( Ecm_proj == 0 ) ? Ecm_eject : R_1*Ecm_eject/Ecm_proj;
  const double E_t3 = 41.0;   // MeV
  double R_3 = ( Ecm_proj < E_t3 ) ? Ecm_proj : E_t1;
  double X_3 = ( Ecm_proj == 0 ) ? Ecm_eject : R_3*Ecm_eject/Ecm_proj;
  slope_a = X_1*(C_1 + C_2*X_1*X_1) + 
    C_3*projectile.Kalbach_M*eject.Kalbach_m*X_3*X_3*X_3*X_3;

  // adjust slope_a parameter for photonuclear reactions as per
  // M. B. Chadwick, P. G. Young, S. Chiba, "Photonuclear
  // angular distribution systematics in the quasideuteron regime",
  // Journal of Nuclear Science and Technology 32, 1154 (1995).
  if( projectile.A == 0 )
  {
    static double m_neutron = Global.Value( "m_neutron" );
    slope_a *= sqrt( Ecm_proj/(2.0*m_neutron) ) * min( 4.0, max( 1.0, 9.3/sqrt( E_out ) ) );
  }
  return slope_a;
}
// ----------- Kalbach_a::copy_masses --------------
// Stores the masses
void Kalbach_a::copy_masses( const particleInfo &particle_info )
{
  projectile.mass = particle_info.mProj;
  target.mass = particle_info.mTarg;
  eject.mass = particle_info.mProd;
  residual.mass = particle_info.mRes;
}

