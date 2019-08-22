/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: mappings.cpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
  Copyright (c) 2017, Lawrence Livermore National Security, LLC.
  Produced at the Lawrence Livermore National Laboratory.
  Written by the LLNL Nuclear Data and Theory group
          (email: mattoon1@llnl.gov)
  LLNL-CODE-725546.
  All rights reserved.
  
  This file is part of the Merced package, used to generate nuclear reaction
  transfer matrices for deterministic radiation transport.
  
  
      Please also read this link - Our Notice and Modified BSD License
  
  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:
      * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the disclaimer below.
      * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the disclaimer (as noted below) in the
        documentation and/or other materials provided with the distribution.
      * Neither the name of LLNS/LLNL nor the names of its contributors may be used
        to endorse or promote products derived from this software without specific
        prior written permission.
  
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
  THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  
  
  Additional BSD Notice
  
  1. This notice is required to be provided under our contract with the U.S.
  Department of Energy (DOE). This work was produced at Lawrence Livermore
  National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
  
  2. Neither the United States Government nor Lawrence Livermore National Security,
  LLC nor any of their employees, makes any warranty, express or implied, or assumes
  any liability or responsibility for the accuracy, completeness, or usefulness of any
  information, apparatus, product, or process disclosed, or represents that its use
  would not infringe privately-owned rights.
  
  3. Also, reference herein to any specific commercial products, process, or services
  by trade name, trademark, manufacturer or otherwise does not necessarily constitute
  or imply its endorsement, recommendation, or favoring by the United States Government
  or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
  herein do not necessarily state or reflect those of the United States Government or
  Lawrence Livermore National Security, LLC, and shall not be used for advertising or
  product endorsement purposes.
  
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

  // calculate the mass of the residual from the Q value
  double true_mRes = ( particle_info.mProj - particle_info.mProd ) +
    particle_info.mTarg - Q_value;
  static double etol = Global.Value( "E_tol" );
  if( ( particle_info.mRes > 0.0 ) && 
      ( abs( true_mRes - particle_info.mRes ) > etol*true_mRes ) )
  {
    Warning( "two_body_map::set_map", "mass of residual inconsistent with Q value" );
  }
  particle_info.mRes = true_mRes;
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
  if( ( Q_value == 0.0 ) && ( alpha == gamma ) )
  {
    return 0.0;
  }
  if( ( ( Q_value > 0.0 ) && ( alpha > gamma ) ) ||
      ( ( Q_value < 0.0 ) && ( alpha < gamma ) ) )
  {
    Warning( "two_body_map::zero_Eout", "no zero outgoing energy" );
    return 0.0;
  }
  double E_zero_in = -beta_eject*Q_value/( alpha - gamma );
  return E_zero_in;
}
// ------------- two_body_map::min_Eout --------------------
// For exothermic reactions, returns the incident energy for minimal outgoing lab energy
double two_body_map::min_Eout( )
{
  if( Q_value <= 0.0 )
  {
    Warning( "two_body_map::min_Eout", "only for exothermic reactions" );
    return 0.0;
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
// On exit pair_0 and pair_1 are tighter bounds.
dd_entry Newton_map_param::find_bottom( double mu, dd_entry *pair_0,
  dd_entry *pair_1, double tol )
{
  if( ( masses->Q_value == 0.0 ) || ( mu >= 0.0 ) )
  {
    FatalError( "two_body_map::find_bottom",
		"improper input" );
  }
  double Ein;
  double Eout;
  // parameters for the function Newtonian_F::T_out_lab
  mu_cm = mu;
  void *params = static_cast< void *>( this );
  while( abs( pair_1->x - pair_0->x ) > tol * pair_1->x )
  {
    Ein = math_F::parabola_bottom( Newtonian_F::T_out_lab, pair_0,
					  pair_1, params );
    // don't go below the threshold
    if( Ein < masses->threshold )
    {
      Ein = masses->threshold;
    }
    Eout = Newtonian_F::T_out_lab( Ein, params );

    // set up the next iteration
    if( Ein <  pair_0->x )
    {
      pair_1->x = pair_0->x;
      pair_1->y = pair_0->y;
      pair_0->x = Ein;
      pair_0->y = Eout;
    }
    else if( Ein <  pair_1->x )
    {
      if( pair_0->y < pair_1->y )
      {
	pair_1->x = Ein;
	pair_1->y = Eout;
      }
      else
      {
	pair_0->x = Ein;
	pair_0->y = Eout;
      }
    }
    else // Ein >=  pair_1->x
    {
      pair_0->x = pair_1->x;
      pair_0->y = pair_1->y;
      pair_1->x = Ein;
      pair_1->y = Eout;
    }
  }
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
