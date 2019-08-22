/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: mappings.cpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/

// implementation of the mappings class

#include <fstream>  // standard file stream package
#include <iostream>
#include <cmath>

#include "mappings.hpp"
#include "math_util.hpp"
#include "messaging.hpp"
#include "global_params.hpp"

using namespace std;

// ********* class particleInfo **********************
// ------------- particleInfo::check_data --------------------
// Check that we have masses
void particleInfo::check_data( ) const
{
  if( mProj < 0.0 )
  {
    FatalError( "particleInfo::check_data", "projectile mass not given" );
  }
  if( mTarg < 0.0 )
  {
    FatalError( "particleInfo::check_data", "target mass not given" );
  }
  if( mProd < 0.0 )
  {
    FatalError( "particleInfo::check_data", "product mass not given" );
  }
}
// ------------- particleInfo::copy --------------------
// Make a copy
void particleInfo::copy( const particleInfo &to_copy )
{
  mProj = to_copy.mProj;
  mTarg = to_copy.mTarg;
  mProd = to_copy.mProd;
  mRes = to_copy.mRes;
}

// ********* class map_cm_lab **********************
// ------------- map_cm_lab::setup_params --------------------
// Sets alpha, beta, and gamma
void map_cm_lab::setup_params( double mProj, double mTarg, double mEject, double mRes )
{
  if( mProj > 0.0 )
  {
    gamma = mProj * mEject /( ( mProj + mTarg ) * ( mProj + mTarg ) );
  }
  else  // incident photon
  {
    gamma =  mEject / mTarg;
  }
  beta_targ = mProj / ( mTarg + mProj );  // for cm kinetic energy of the target
  beta_proj = mTarg / ( mTarg + mProj );  // for cm kinetic energy of the projectile
  beta_eject = mRes / ( mEject + mRes );  // for cm kinetic energy of the ejectile
  alpha = beta_eject * mTarg / ( mProj + mTarg );
}
// ------------- map_cm_lab::setup_ratios --------------------
// Sets beta and gamma mass ratios
void map_cm_lab::setup_ratios( double mProj, double mTarg, double mEject )
{
  // alpha and beta_eject are not used, so we may ignore mRes
  setup_params( mProj, mTarg, mEject, 0.0 );
}
// ------------- map_cm_lab::get_Elab --------------------
// Gets the laboratory energy of outgoing particle
double map_cm_lab::get_Elab( double E_in, double E_cm, double mu_cm )
{
  double E0 = get_Etrans( E_in );
  double E_lab = E0 + E_cm + 2*mu_cm*sqrt( E0 * E_cm );
  return E_lab;
}
// ------------- map_cm_lab::get_E_mu_lab --------------------
// Gets the laboratory energy and cosine of outgoing particle
void map_cm_lab::get_E_mu_lab( double E_in, double E_cm, double mu_cm,
  double *Eout_lab, double *mu_lab )
{
  *Eout_lab = get_Elab( E_in, E_cm, mu_cm );
  if( *Eout_lab <= 0.0 )
  {
    *Eout_lab = 0.0;
    *mu_lab = 0.0;
  }
  else
  {
    // use trigonometry: Vout_lab*mu_lab = V0 + V_cm*mu_cm
    double V0 = sqrt( get_Etrans( E_in ) );
    double V_cm = sqrt( E_cm );
    static double etol = Global.Value( "E_tol" );
    *mu_lab = ( V0 + V_cm*mu_cm )/sqrt( *Eout_lab );
    if( *mu_lab > 1.0 )
    {
      if( *mu_lab > 1.0 + etol )
      {
        Warning( "map_cm_lab::get_E_mu_lab", "mu_lab > 1" );
      }
      *mu_lab = 1.0;
    }
    else if( *mu_lab < -1.0 )
    {
      if( *mu_lab < -1.0 - etol )
      {
        Warning( "map_cm_lab::get_E_mu_lab", "mu_lab < -1" );
      }
      *mu_lab = -1.0;
    }
  }
}
// ------------- map_cm_lab::get_mu_cm --------------------
// Gets the center-of-mass cosine for given incident energy, E_out_lab, and E_cm 
double map_cm_lab::get_mu_cm( double E_in, double E_out_lab, double E_cm )
{
  double E_trans = get_Etrans( E_in );
  if( E_cm*E_trans <= 0.0 )
  {
    FatalError( "map_cm_lab::get_mu_cm", "square root out of range" );
  }
  double mu = ( E_out_lab - E_cm - E_trans )/(2*sqrt( E_trans*E_cm ));
  if( mu > 1.0 ) mu = 1.0;
  if( mu < -1.0 ) mu = -1.0;
  return mu;
}

// ********* class two_body_map **********************
// ------------- two_body_map::set_map --------------------
void two_body_map::set_map( particleInfo &particle_info, double Q )
// set up the map from center-of-mass to laboratory coordinates
{
  particle_info.check_data( );
  Q_value = Q;
  Ecm_photon_residual = -1.0;  // this is reset if it is used

  static double etol = Global.Value( "E_tol" );

  // check for photon as residual
  if( particle_info.mRes == 0.0 )
  {
    double true_Q = ( particle_info.mProj + particle_info.mTarg ) -
      particle_info.mProd;
    if( abs( true_Q - Q_value ) > etol*abs( true_Q ) )
    {
      Warning( "two_body_map::set_map", "Q value inconsistent with photon residual" );
    }
    Q_value = true_Q;
    Ecm_photon_residual = Q_value * Q_value /
      ( 2.0 * ( particle_info.mProj + particle_info.mTarg ) );
  }
  else
  {
    // calculate the mass of the residual from the Q value
    double true_mRes = ( particle_info.mProj - particle_info.mProd ) +
      particle_info.mTarg - Q_value;
    if( ( particle_info.mRes > 0.0 ) && 
        ( abs( true_mRes - particle_info.mRes ) > etol*true_mRes ) )
    {
      Warning( "two_body_map::set_map", "mass of residual inconsistent with Q value" );
    }
    particle_info.mRes = true_mRes;
  }

  setup_params( particle_info.mProj, particle_info.mTarg, particle_info.mProd,
               particle_info.mRes );

  // parameters for min_Eout
  Newt_gamma = 2.0 * particle_info.mProd * beta_eject;
  mass_sum = particle_info.mTarg + particle_info.mProj;
  mass_product = 2.0 * particle_info.mProj * particle_info.mProd * particle_info.mProd;
}
// ------------- two_body_map::get_threshold --------------------
// Calculates the threshold energy
void two_body_map::get_threshold( )
{
  threshold = ( Q_value >= 0.0 ) ? 0.0: -beta_eject*Q_value/alpha;
}
// ------------- two_body_map::get_Ecm --------------------
// Gets the center-of-mass energy of outgoing particle
double two_body_map::get_Ecm( double E_in )
{
  double E_cm = alpha*E_in + beta_eject*Q_value;
  static double etol = Global.Value( "E_tol" );
  if( E_cm < -etol )
  {
    Warning( "two_body_map::get_Ecm",
           pastenum( "negative cm-energy ", E_cm ) );
  }
  // we could be just below threshold
  if ( E_cm < 0.0 ){
    E_cm = 0.0;
  }
  return E_cm;
}
// ------------- two_body_map::get_mu_cm --------------------
// Gets the value of mu_cm given the incident energy and the
// laboratory energy of outgoing particle.
double two_body_map::get_mu_cm( double E_in, double E_lab )
{
  double E0 = get_Etrans( E_in );
  double E_cm = get_Ecm( E_in );
  double denom = 2 * sqrt( E0 * E_cm );
  double mu_cm = ( denom == 0.0 ) ? 0.0 : ( E_lab - E0 - E_cm ) / denom;
  return mu_cm;
}
// ------------- two_body_map::two_body_get_E_mu_lab --------------------
// Gets the laboratory energy and cosine of outgoing particle for discrete 2-body reactions
void two_body_map::two_body_get_E_mu_lab( double E_in, double mu_cm, double *Eout_lab,
  double *mu_lab )
{
  double E_cm = get_Ecm( E_in );
  get_E_mu_lab( E_in, E_cm, mu_cm, Eout_lab, mu_lab );
}
// ------------- two_body_map::zero_Eout --------------------
// For endothermic reactions, returns the incident energy for zero outgoing lab energy
double two_body_map::zero_Eout( )
{
  double E_zero_in;

  if( ( Q_value == 0.0 ) && ( alpha == gamma ) )
  {
    E_zero_in = 0.0;
  }
  else if( ( ( Q_value > 0.0 ) && ( alpha > gamma ) ) ||
	   ( ( Q_value < 0.0 ) && ( alpha < gamma ) ) ||
	   ( ( Q_value < 0.0 ) && ( Ecm_photon_residual > 0.0 ) ) )
  {
    Warning( "two_body_map::zero_Eout", "no zero outgoing energy" );
    E_zero_in = threshold;
  }
  else if( Ecm_photon_residual < 0.0 ) // the residual is not a photon
  {
    E_zero_in = -beta_eject*Q_value/( alpha - gamma );
  }
  else  // Q > 0 and the residual is a photon
  {
    E_zero_in = Ecm_photon_residual / gamma;
  }

  return E_zero_in;
}
// ------------- two_body_map::min_Eout --------------------
// For exothermic reactions, returns the incident energy for minimal outgoing lab energy
double two_body_map::min_Eout( )
{
  if( ( Q_value <= 0.0 ) && ( Ecm_photon_residual < 0.0 ) )
  {
    Warning( "two_body_map::min_Eout", "only for exothermic reactions" );
    return threshold;
  }
  // For mu = -1 find E_in such that dV_0/dE_in = dvcm_out/dE_in
  // parameters Newt_gamma and Newt_delta such that
  //  pcm^2 = Newt_gamma * ( Q + Newt_delta * E_in )
  double E_min_in = mass_product * Q_value /
    ( beta_proj * ( Newt_gamma *  beta_proj * mass_sum * mass_sum - mass_product ) );
  return E_min_in;
}

// ********* class Newton_map_param **********************
// ------------- Newton_map_param::find_bottom --------------------
// For mu < 0 find the incident energy which minimizes Eout
dd_entry Newton_map_param::find_bottom( double mu )
{
  if( ( masses->Q_value == 0.0 ) || ( mu >= 0.0 ) )
  {
    FatalError( "two_body_map::find_bottom",
		"improper input" );
  }

  double alpha_gamma = masses->alpha + masses->gamma;
  double Delta = sqrt( alpha_gamma * alpha_gamma -
		       4 * mu * mu * masses->alpha * masses->gamma );
  double numerator;
  double denominator;

  if( masses->Ecm_photon_residual < 0.0 ) // the residual is not a photon
  {
    if( masses->Q_value > 0.0 )
    {
      double diff_sq = alpha_gamma - Delta;
      diff_sq *= diff_sq;
      numerator = masses->beta_eject * diff_sq * masses->Q_value;
      denominator = masses->alpha * ( 4 * mu * mu * masses->alpha * masses->gamma -
				    diff_sq );
    }
    else
    {
      double sum_sq = alpha_gamma + Delta;
      sum_sq *= sum_sq;
      numerator = masses->beta_eject * sum_sq * masses->Q_value;
      denominator = masses->alpha * ( 4 * mu * mu * masses->alpha * masses->gamma -
				      sum_sq );
    }
  }
  else
  {
    numerator = mu * mu * masses->Ecm_photon_residual;
    denominator = masses->gamma;
  }

  double Ein = numerator / denominator;

  // parameters for the function Newtonian_F::T_out_lab
  mu_cm = mu;
  void *params = static_cast< void *>( this );
  // don't go below the threshold
  if( Ein < masses->threshold )
  {
    Ein = masses->threshold;
  }
  double Eout = Newtonian_F::T_out_lab( Ein, params );

  dd_entry ans( Ein, Eout );
  return ans;
}
// ------------- Newton_map_param::find_lowest_bottom --------------------
// Find the incident energy which minimizes Eout for mu = -1
dd_entry Newton_map_param::find_lowest_bottom(  )
{
  dd_entry ans;

  // parameters for the function Newtonian_F::T_out_lab
  void *params = static_cast< void *>( this );
 
  if( ( masses->Q_value == 0.0 ) ||
      ( masses->alpha == masses->gamma ) )
  {
    ans.x = 0.0;
    ans.y = 0.0;
  }
  else if( ( ( masses->Q_value < 0.0 ) && ( masses->alpha > masses->gamma ) ) ||
	   ( ( masses->Q_value > 0.0 ) && ( masses->alpha < masses->gamma ) ) )
  {
    ans.x = masses->zero_Eout( );
    ans.y = 0.0;
  }
  else
  {
    ans.x = masses->min_Eout( );
    mu_cm = 0.0;
    ans.y = Newtonian_F::T_out_lab( ans.x, params );
  }
  return ans;
}
// ------------- Newton_map_param::find_hit --------------------
// Finds the incident energy for given T_lab and mu_cm.
// Returns the solutions.
double Newton_map_param::find_hit( double T_lab, double mucm, const dd_entry &pair_0,
		   const dd_entry &pair_1 )
{
  // set up the void* parameters for the root finder
  mu_cm = mucm;
  void *params = static_cast< void * >( this );

  // use copies of the initial pairs
  dd_entry use_pair_0 = pair_0;
  dd_entry use_pair_1 = pair_1;
  static double E_tol = Global.Value( "E_tol" );
  double this_Tin = math_F::zeroin( Newtonian_F::T_out_lab, T_lab,
				    use_pair_0, use_pair_1, params, E_tol );
  return this_Tin;
}

// ********* class Ecm_intersect **********************
// ------------- Ecm_intersect::set_energies --------------------
void Ecm_intersect::set_energies( double Ein, double Ecm_out )
{
  E_in = Ein;
  double E_trans = get_Etrans( );
  E_cm = Ecm_out;
  V_trans = sqrt( E_trans );
  V_cm = sqrt( E_cm );
}
// ------------- Ecm_intersect::set_hit_G --------------------
// This Ecm curve hits the E_lab=const. curve if and only if hit_G(E_lab) >= 0
void Ecm_intersect::set_hit_G( double E_lab )
{
  double E_trans = get_Etrans( );
  if( E_cm == 0.0 )
  {
    G = E_trans - E_lab;
    G *= -G;
  }
  else
  {
    G = 2*E_lab*( E_trans + E_cm ) - ( E_trans - E_cm )*( E_trans - E_cm ) - E_lab*E_lab;
  }
}
// ------------- Ecm_intersect::interpolate_Ein --------------------
// Sets up the parameters for an intermediate incident energy
void Ecm_intersect::interpolate_Ein( double Ein, double E_lab,
  const Ecm_intersect &low_data, const Ecm_intersect &high_data )
{
  gamma = low_data.gamma;
  double alpha_ratio = ( Ein - low_data.E_in )/( high_data.E_in - low_data.E_in );
  E_cm = ( 1.0 - alpha_ratio )*low_data.E_cm + alpha_ratio*high_data.E_cm;
  set_energies( Ein, E_cm );
  set_hit_G( E_lab );
}

// ********* class phase_space_map **********************
// ----------- phase_space_map::set_data --------------
// Sets the private members
void phase_space_map::set_data( int numParticles, double mEject, double totalMass, double Q )
{
  num_particles = numParticles;
  Q_value = Q;
  out_mass_ratio = ( totalMass - mEject )/totalMass;
  exponent = 1.5*num_particles - 4;
  // set the normalization factor
  switch( num_particles )
  {
  case 3:
    normalize = 4.0/M_PI;
    break;
  case 4:
    normalize = 105.0/32;
    break;
  case 5:
    normalize = 256/(14 * M_PI);
    break;
  default:
    if( num_particles < 3 )
    {
      FatalError( "phase_space_map::set_data",
		  "the model makes no sense for fewer than 3 outgoing particles" );
    }
    else
    {
      FatalError( "phase_space_map::set_data",
		"Implement more than 5 outgoing particles" );
    }
  }
}
// ----------- phase_space_map::get_threshold --------------
// Returns the reaction threshold for the phase-space model
double phase_space_map::get_threshold( )
{
  return ( - Q_value / beta_proj );
}
// ----------- phase_space_map::get_Ecm_max --------------
// Returns the maximum center-of-mass energy of the outgoing particle for the phase-space model
double phase_space_map::get_Ecm_max( double E_in )
{
  double availableE = beta_proj * E_in + Q_value;
  if( availableE < 0.0 ) availableE = 0.0;
  return out_mass_ratio * availableE;
}
// ----------- phase_space_map::get_prob --------------
// Returns the center-of-mass energy probability density
double phase_space_map::get_prob( double E_cm, double Ecm_max )
{
  if( Ecm_max <= 0.0 )
  {
    FatalError( "phase_space_map::get_prob", "Negative maximum energy" );
  }
  double norm = 0.0;
  switch( num_particles )
  {
  case 3:
    norm = normalize / ( Ecm_max * Ecm_max );
    break;
  case 4:
    norm = normalize / pow( Ecm_max, 3.5 );
    break;
  case 5:
    norm = normalize / pow( Ecm_max, 5 );
    break;
  }
  return norm * sqrt( E_cm ) * pow( Ecm_max-E_cm, exponent );
}


// ***************** functions ********************
// ------------- Newtonian_F::T_out_lab --------------------
// Returns the lab-frame kinetic energy of the outgoing particle
double Newtonian_F::T_out_lab( double T_lab_in, void *params )
{
  // The parameters are really a Newton_map_param *
  Newton_map_param *map = static_cast< Newton_map_param * >( params );
  // get the center-of-mass energy of the emitted particle
  double T_cm = map->masses->get_Ecm( T_lab_in );
  double T = map->masses->get_Elab( T_lab_in, T_cm, map->mu_cm );
  return T;
}
