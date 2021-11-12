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

// ********* class Maps::particleInfo **********************
// ------------- Maps::particleInfo::check_data --------------------
// Check that we have masses
void Maps::particleInfo::check_data( ) const
{
  if( mProj < 0.0 )
  {
    Msg::FatalError( "Maps::particleInfo::check_data",
		     "projectile mass not given" );
  }
  if( mTarg < 0.0 )
  {
    Msg::FatalError( "Maps::particleInfo::check_data",
		     "target mass not given" );
  }
  if( mProd < 0.0 )
  {
    Msg::FatalError( "Maps::particleInfo::check_data",
		     "product mass not given" );
  }
}
// ------------- Maps::particleInfo::copy --------------------
// Make a copy
void Maps::particleInfo::copy( const Maps::particleInfo &to_copy )
{
  mProj = to_copy.mProj;
  mTarg = to_copy.mTarg;
  mProd = to_copy.mProd;
  mRes = to_copy.mRes;
}
// ------------- Maps::particleInfo::photon_in_out --------------------
// Returns true if any particle is a photon
bool Maps::particleInfo::photon_in_out( )
{
  bool is_photon = false;
  if( mProd * mRes * mProj == 0.0 )
  {
    is_photon = true;
  }
  return is_photon;
}

// ********* class Maps::map_cm_lab **********************
// ------------- Maps::map_cm_lab::setup_params --------------------
// Sets alpha, beta, and gamma
void Maps::map_cm_lab::setup_params( double mProj, double mTarg, double mEject, double mRes )
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
// ------------- Maps::map_cm_lab::setup_ratios --------------------
// Sets beta and gamma mass ratios
void Maps::map_cm_lab::setup_ratios( double mProj, double mTarg, double mEject )
{
  // alpha and beta_eject are not used, so we may ignore mRes
  setup_params( mProj, mTarg, mEject, 0.0 );
}
// ------------- Maps::map_cm_lab::get_Elab --------------------
// Gets the laboratory energy of outgoing particle
double Maps::map_cm_lab::get_Elab( double E_in, double E_cm, double mu_cm )
{
  double E0 = get_Etrans( E_in );
  double E_lab = E0 + E_cm + 2*mu_cm*sqrt( E0 * E_cm );
  return E_lab;
}
// ------------- Maps::map_cm_lab::get_E_mu_lab --------------------
// Gets the laboratory energy and cosine of outgoing particle
void Maps::map_cm_lab::get_E_mu_lab( double E_in, double E_cm, double mu_cm,
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
    static double etol = Global.Value( "looser_tol" );
    *mu_lab = ( V0 + V_cm*mu_cm )/sqrt( *Eout_lab );
    if( *mu_lab > 1.0 )
    {
      if( *mu_lab > 1.0 + etol )
      {
        Msg::Warning( "Maps::map_cm_lab::get_E_mu_lab", "mu_lab > 1" );
      }
      *mu_lab = 1.0;
    }
    else if( *mu_lab < -1.0 )
    {
      if( *mu_lab < -1.0 - etol )
      {
        Msg::Warning( "Maps::map_cm_lab::get_E_mu_lab", "mu_lab < -1" );
      }
      *mu_lab = -1.0;
    }
  }
}
// ------------- Maps::map_cm_lab::get_mu_cm --------------------
// Gets the center-of-mass cosine for given incident energy, E_out_lab, and E_cm 
double Maps::map_cm_lab::get_mu_cm( double E_in, double E_out_lab, double E_cm )
{
  double E_trans = get_Etrans( E_in );
  if( E_cm*E_trans <= 0.0 )
  {
    Msg::FatalError( "Maps::map_cm_lab::get_mu_cm",
		     "square root out of range" );
  }
  double mu = ( E_out_lab - E_cm - E_trans )/(2*sqrt( E_trans*E_cm ));
  if( mu > 1.0 ) mu = 1.0;
  if( mu < -1.0 ) mu = -1.0;
  return mu;
}

// ********* class Maps::two_body_map **********************
// ------------- Maps::two_body_map::set_map --------------------
void Maps::two_body_map::set_map( Maps::particleInfo &particle_info, double Q )
// set up the map from center-of-mass to laboratory coordinates
{
  particle_info.check_data( );
  Q_value = Q;
  Ecm_photon_residual = -1.0;  // this is reset if it is used

  static double etol = Global.Value( "looser_tol" );

  // check for photon as residual
  if( particle_info.mRes == 0.0 )
  {
    double true_Q = ( particle_info.mProj + particle_info.mTarg ) -
      particle_info.mProd;
    if( std::abs( true_Q - Q_value ) > etol*std::abs( true_Q ) )
    {
      Msg::Warning( "Maps::two_body_map::set_map",
		    "Q value inconsistent with photon residual" );
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
        ( std::abs( true_mRes - particle_info.mRes ) > etol*true_mRes ) )
    {
      Msg::Warning( "Maps::two_body_map::set_map",
		    "mass of residual inconsistent with Q value" );
    }
    particle_info.mRes = true_mRes;
  }

  setup_params( particle_info.mProj, particle_info.mTarg, particle_info.mProd,
               particle_info.mRes );

}
// ------------- Maps::two_body_map::get_threshold --------------------
// Calculates the threshold energy
void Maps::two_body_map::get_threshold( )
{
  threshold = ( Q_value >= 0.0 ) ? 0.0: -beta_eject*Q_value/alpha;
}
// ------------- Maps::two_body_map::get_Ecm --------------------
// Gets the center-of-mass energy of outgoing particle
double Maps::two_body_map::get_Ecm( double E_in )
{
  double E_cm = alpha*E_in + beta_eject*Q_value;
  static double etol = Global.Value( "looser_tol" );
  if( E_cm < -etol )
  {
    Msg::Warning( "Maps::two_body_map::get_Ecm",
           Msg::pastenum( "negative cm-energy ", E_cm ) );
  }
  // we could be just below threshold
  if ( E_cm < 0.0 ){
    E_cm = 0.0;
  }
  return E_cm;
}
// ------------- Maps::two_body_map::get_mu_cm --------------------
// Gets the value of mu_cm given the incident energy and the
// laboratory energy of outgoing particle.
double Maps::two_body_map::get_mu_cm( double E_in, double E_lab )
{
  double E0 = get_Etrans( E_in );
  double E_cm = get_Ecm( E_in );
  double denom = 2 * sqrt( E0 * E_cm );
  double mu_cm = ( denom == 0.0 ) ? 0.0 : ( E_lab - E0 - E_cm ) / denom;
  return mu_cm;
}
// ------------- Maps::two_body_map::two_body_get_E_mu_lab --------------------
// Gets the laboratory energy and cosine of outgoing particle for discrete 2-body reactions
void Maps::two_body_map::two_body_get_E_mu_lab( double E_in, double mu_cm, double *Eout_lab,
  double *mu_lab )
{
  double E_cm = get_Ecm( E_in );
  get_E_mu_lab( E_in, E_cm, mu_cm, Eout_lab, mu_lab );
}
// ------------- Maps::two_body_map::zero_Eout --------------------
// Returns the incident energy for zero outgoing lab energy when this exists
double Maps::two_body_map::zero_Eout( )
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
    Msg::Warning( "Maps::two_body_map::zero_Eout",
		  "no zero outgoing energy" );
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
// ------------- Maps::two_body_map::min_Eout --------------------
// Returns the incident energy for minimal outgoing lab energy when this minimum is positive
double Maps::two_body_map::min_Eout( )
{
  double E_min_in = 0.0;

  // For mu = -1 find E_in such that dV_0/dE_in = dvcm_out/dE_in
  if( Q_value * ( alpha - gamma ) > 0.0 )
  {
    E_min_in = - gamma * beta_eject * Q_value / ( alpha * ( gamma - alpha ) );
  }
  else if( Q_value * ( alpha - gamma ) < 0.0 )
  {
    E_min_in = beta_eject * Q_value / ( gamma - alpha );
  }
  else //  if( Q_value * ( alpha - gamma ) == 0.0 )
  {
    E_min_in = threshold;
  }

  return E_min_in;
}

  
// ********* class Maps::Newton_map_param **********************
// ------------- Maps::Newton_map_param::find_bottom --------------------
// For mu < 0 find the incident energy which minimizes Eout
Ddvec::dd_entry Maps::Newton_map_param::find_bottom( double mu )
{
  if( ( masses->Q_value == 0.0 ) || ( mu >= 0.0 ) )
  {
    Msg::FatalError( "Maps::Newton_map_param::find_bottom",
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

  Ddvec::dd_entry ans( Ein, Eout );
  return ans;
}
// ------------- Maps::Newton_map_param::find_lowest_bottom --------------------
// Find the incident energy which minimizes Eout for mu = -1
Ddvec::dd_entry Maps::Newton_map_param::find_lowest_bottom(  )
{
  Ddvec::dd_entry ans;

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
// ------------- Maps::Newton_map_param::find_hit --------------------
// Finds the incident energy for given T_lab and mu_cm.
// Returns the solutions.
double Maps::Newton_map_param::find_hit( double T_lab, double mucm, const Ddvec::dd_entry &pair_0,
		   const Ddvec::dd_entry &pair_1 )
{
  // set up the void* parameters for the root finder
  mu_cm = mucm;
  void *params = static_cast< void * >( this );

  // use copies of the initial pairs
  Ddvec::dd_entry use_pair_0 = pair_0;
  Ddvec::dd_entry use_pair_1 = pair_1;
  static double E_tol = Global.Value( "looser_tol" );
  double this_Tin = math_F::zeroin( Newtonian_F::T_out_lab, T_lab,
				    use_pair_0, use_pair_1, params, E_tol );
  return this_Tin;
}

// ********* class Maps::two_step_map **********************
// ------------- Maps::two_step_map::set_map --------------------
void Maps::two_step_map::set_map( Maps::particleInfo &step1_particles, double Q,
			    Maps::particleInfo &step2_particles )
// set up the map from center-of-mass to laboratory coordinates
{
  static double etol = Global.Value( "looser_tol" );

  // information for step 1
  // Check the mass of the excited product
  double true_mProd = step1_particles.mProj + step1_particles.mTarg -
    step1_particles.mRes - Q;
  if( ( step1_particles.mProd > 0.0 ) &&
      ( std::abs( step1_particles.mProd - true_mProd ) > etol*true_mProd ) )
  {
    Msg::Warning( "Maps::two_step_map::set_map",
	     "mass of product from step 1 inconsistent" );
  }
  step1_particles.mProd = true_mProd;
  
  setup_params( step1_particles.mProj, step1_particles.mTarg, step1_particles.mProd,
		step1_particles.mRes );

  Q_value = Q;

  if( Q * ( gamma - alpha ) >= 0.0 )
  {
    Msg::FatalError( "Maps::two_step_map::set_map",
		"write code for this case" );
  }

  // set up step 2
  step2_particles.mProj = 0.0;
  step2_particles.mTarg = step1_particles.mProd;

  // Q for step 2
  double true_Q = step2_particles.mTarg -
    ( step2_particles.mProd +
      step2_particles.mRes );

  // set the mass ratios for step 2
  step_2_map.set_map( step2_particles, true_Q );

  // reset alpha and gamma for step 2 of the reaction
  // The energy of the final product in this frame depends only on Q for step 2.
  step_2_map.alpha = 0.0;
  
  // Use the energy of the residual from step 1 as E_trans for step 2,
  // scaled by the ratio of masses
  step_2_map.gamma = step2_particles.mProd / step2_particles.mTarg;

  //! the energy of the product in the frame of the second step
  E_cm2_prod = true_Q * step_2_map.beta_eject;
}
// ------------- Maps::two_step_map::two_step_get_E_lab --------------------
// Gets the laboratory energy  outgoing particle for 2-step reactions
double Maps::two_step_map::two_step_get_E_lab( double E_in, double mucm1, double mucm2 )
{
  Maps::two_step_map_param map_param;
  map_param.map = this;
  
  map_param.mucm_1 = mucm1;
  map_param.mucm_2 = mucm2;

  // Use this as a function parameter
  void *params = static_cast< void *>( &map_param );
  double Eout = Newtonian_F::two_step_T_out_lab( E_in, params );

  return Eout;
}
// ------------- Maps::two_step_map::two_step_get_mu_lab --------------------
// Gets the laboratory direction cosine of outgoing particle for a 2-step reactions
double Maps::two_step_map::two_step_get_mu_lab( double E_in, double mucm1,
		double mucm2, double Etrans2, double Eout_2, double w )
{
  double mu_lab;

  // We need mu_lab1 from step 1.
  double Eout_lab1;
  double mu_lab1;
  two_body_get_E_mu_lab( E_in, mucm1, &Eout_lab1, &mu_lab1 );
  
  // for the direction cosine after the second step
  // energy of product in frame of residual from step 1
  double V2_cm = std::sqrt( E_cm2_prod );
  double xi = V2_cm * mucm2;
  double eta = V2_cm * std::sqrt( 1 - mucm2*mucm2 ) * std::sin( w );
  double V2_trans = std::sqrt( Etrans2 );
  double V2_lab = std::sqrt( Eout_2 );
  mu_lab = ( mu_lab1 * ( V2_trans + xi ) -
	      eta * std::sqrt( 1 - mu_lab1*mu_lab1 ) ) / V2_lab;

  static double etol = Global.Value( "looser_tol" );
  if( mu_lab > 1.0 )
  {
    if( mu_lab > 1.0 + etol )
    {
      Msg::Warning( "Maps::two_step_map::two_step_get_mu_lab", "mu_lab > 1" );
    }
    mu_lab = 1.0;
  }
  else if( mu_lab < -1.0 )
  {
    if( mu_lab < -1.0 - etol )
    {
      Msg::Warning( "Maps::two_step_map::two_step_get_mu_lab", "mu_lab < -1" );
    }
    mu_lab = -1.0;
  }

  return mu_lab;
}
// ------------- Maps::two_step_map::get_mucm1 --------------------
// Returns the direction cosine for the first step for given incident energy and desired outgoing energy
double Maps::two_step_map::get_mucm1( double E_in, int mucm2, double Eout_lab )
{
  double V_desired = std::sqrt( Eout_lab );
  double V_lab1 = ( V_desired - mucm2 * std::sqrt( E_cm2_prod ) ) /
    std::sqrt( step_2_map.gamma );
  double E_lab1 = V_lab1 * V_lab1;
  double E_trans1 = get_Etrans( E_in );
  double E_cm1 = get_Ecm( E_in );
  
  double mucm1 = ( E_lab1 - E_trans1 - E_cm1 ) /
    ( 2 * std::sqrt( E_trans1 * E_cm1 ) );

  static double etol = Global.Value( "looser_tol" );
  
  if( mucm1 > 1.0 )
  {
    if( mucm1 > 1.0 + etol )
    {
      Msg::Warning( "Maps::two_step_map::get_mucm1", "mucm1 > 1" );
    }
    mucm1 = 1.0;
  }
  else if( mucm1 < -1.0 )
  {
    if( mucm1 < -1.0 - etol )
    {
      Msg::Warning( "Maps::two_step_map::get_mucm1", "mucm1 < -1" );
    }
    mucm1 = -1.0;
  }

  return mucm1;	   
}
// ------------- Maps::two_step_map::get_mucm2 --------------------
// Returns the direction cosine for the second step for given translational energy and desired outgoing energy
double Maps::two_step_map::get_mucm2( double Etrans2, double Eout_lab )
{
  double mucm2 = ( Eout_lab - Etrans2 - E_cm2_prod ) /
    ( 2 * std::sqrt( Etrans2 * E_cm2_prod ) );

  static double etol = Global.Value( "looser_tol" );
  
  if( mucm2 > 1.0 )
  {
    if( mucm2 > 1.0 + etol )
    {
      Msg::Warning( "Maps::two_step_map::get_mucm2", "mucm2 > 1" );
    }
    mucm2 = 1.0;
  }
  else if( mucm2 < -1.0 )
  {
    if( mucm2 < -1.0 - etol )
    {
      Msg::Warning( "Maps::two_step_map::get_mucm2", "mucm2 < -1" );
    }
    mucm2 = -1.0;
  }

  return mucm2;	   
}

// ********* class Maps::two_step_map_param **********************
// ------------- Maps::two_step_map_param::find_hit --------------------
// Finds the incident energy for given E_lab
double Maps::two_step_map_param::find_hit( double E_lab,
		    const Ddvec::dd_entry &pair_0, const Ddvec::dd_entry &pair_1 )
{
  // set up the void* parameters for the root finder
  void *params = static_cast< void * >( this );

  // use copies of the initial pairs
  Ddvec::dd_entry use_pair_0 = pair_0;
  Ddvec::dd_entry use_pair_1 = pair_1;
  static double E_tol = Global.Value( "looser_tol" );
  double this_Tin = math_F::zeroin( Newtonian_F::two_step_T_out_lab, E_lab,
				    use_pair_0, use_pair_1, params, E_tol );
  return this_Tin;
}
// ------------- Maps::two_step_map_param::find_lowest_bottom --------------------
// Find the incident energy which minimizes Eout for mucm_1 = mucm_2 = -1
Ddvec::dd_entry Maps::two_step_map_param::find_lowest_bottom(  )
{
  Ddvec::dd_entry ans;

  if( ( map->Q_value == 0.0 ) ||
      ( map->alpha == map->gamma ) )
  {
    Msg::FatalError( "Maps::two_step_map_param::find_lowest_bottom",
		"option 1 not implemented" );
    /*
    ans.x = 0.0;
    ans.y = 0.0;
    */
  }
  else if( ( ( map->Q_value < 0.0 ) && ( map->alpha > map->gamma ) ) ||
	   ( ( map->Q_value > 0.0 ) && ( map->alpha < map->gamma ) ) )
  {
    Msg::FatalError( "Maps::two_step_map_param::find_lowest_bottom",
		"option 2 not implemented" );
    /*
    ans.x = map->zero_Eout( );
    ans.y = 0.0;
    */
  }
  else
  {
    ans.x = map->min_Eout( );
    ans.y = map->two_step_get_E_lab( ans.x, -1.0, -1.0 );
  }
  return ans;
}

// ********* class Maps::Ecm_intersect **********************
// ------------- Maps::Ecm_intersect::set_energies --------------------
// Sets the energies
void Maps::Ecm_intersect::set_energies( double Ein, double Ecm_out )
{
  E_in = Ein;
  double E_trans = get_Etrans( );
  E_cm = Ecm_out;
  V_trans = sqrt( E_trans );
  V_cm = sqrt( E_cm );
}
// ------------- Maps::Ecm_intersect::set_hit_G --------------------
// This Ecm curve hits the E_lab=const. curve if and only if hit_G(E_lab) >= 0
void Maps::Ecm_intersect::set_hit_G( double E_lab )
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
// ------------- Maps::Ecm_intersect::interpolate_Ein --------------------
// Sets up the parameters for an intermediate incident energy
void Maps::Ecm_intersect::interpolate_Ein( double Ein, double E_lab,
  const Maps::Ecm_intersect &low_data, const Maps::Ecm_intersect &high_data )
{
  gamma = low_data.gamma;
  double alpha_ratio = ( Ein - low_data.E_in )/( high_data.E_in - low_data.E_in );
  E_cm = ( 1.0 - alpha_ratio )*low_data.E_cm + alpha_ratio*high_data.E_cm;
  set_energies( Ein, E_cm );
  set_hit_G( E_lab );
}

// ********* class Maps::phase_space_map **********************
// ----------- Maps::phase_space_map::set_data --------------
// Sets the private members
void Maps::phase_space_map::set_data( int numParticles, double mEject, double totalMass, double Q )
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
      Msg::FatalError( "Maps::phase_space_map::set_data",
		  "the model makes no sense for fewer than 3 outgoing particles" );
    }
    else
    {
      Msg::FatalError( "Maps::phase_space_map::set_data",
		"Implement more than 5 outgoing particles" );
    }
  }
}
// ----------- Maps::phase_space_map::get_threshold --------------
// Returns the reaction threshold for the phase-space model
double Maps::phase_space_map::get_threshold( )
{
  return ( - Q_value / beta_proj );
}
// ----------- Maps::phase_space_map::get_Ecm_max --------------
// Returns the maximum center-of-mass energy of the outgoing particle for the phase-space model
double Maps::phase_space_map::get_Ecm_max( double E_in )
{
  double availableE = beta_proj * E_in + Q_value;
  if( availableE < 0.0 ) availableE = 0.0;
  return out_mass_ratio * availableE;
}
// ----------- Maps::phase_space_map::get_prob --------------
// Returns the center-of-mass energy probability density
double Maps::phase_space_map::get_prob( double E_cm, double Ecm_max )
{
  if( Ecm_max <= 0.0 )
  {
    Msg::FatalError( "Maps::phase_space_map::get_prob",
		     "Negative maximum energy" );
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
  // The parameters are really a Maps::Newton_map_param *
  Maps::Newton_map_param *map = static_cast< Maps::Newton_map_param * >( params );
  // get the center-of-mass energy of the emitted particle
  double T_cm = map->masses->get_Ecm( T_lab_in );
  double T = map->masses->get_Elab( T_lab_in, T_cm, map->mu_cm );
  return T;
}
// ------------- Newtonian_F::two_step_T_out_lab --------------------
// For a 2-step reaction, returns the lab-frame kinetic energy of the outgoing particle
double Newtonian_F::two_step_T_out_lab( double T_lab_in, void *params )
{
  // The parameters are really a Maps::two_step_map_param *
  Maps::two_step_map_param *map_param = static_cast< Maps::two_step_map_param * >( params );
  // get the center-of-mass energy of the emitted particle for step 1
  double T_cm1 = map_param->map->get_Ecm( T_lab_in );
  double T1 = map_param->map->get_Elab( T_lab_in, T_cm1, map_param->mucm_1 );

  // for step 2
  double T2 = map_param->map->step_2_map.get_Elab( T1, map_param->map->E_cm2_prod,
						 map_param->mucm_2 );
  return T2;
}
