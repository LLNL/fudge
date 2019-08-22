/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2015-05-20 -0800 (Wed, 20 May 2015) $
 * $Author: hedstrom $
 * $Id: relativistic.cpp 1 2015-05-20 03:06:56Z hedstrom $
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

#include <fstream>  // standard file stream package
#include <iostream>
#include <cmath>

#include "relativistic.hpp"
#include "math_util.hpp"
#include "messaging.hpp"
#include "global_params.hpp"

using namespace std;

// implementation of the classes used in relativistic mechanics

// ********* class relativistic_masses **********************
// ------------- relativistic_masses::setup_masses --------------------
// Saves the rest masses
void relativistic_masses::setup_masses( particleInfo *to_save, double file_Q )
{
  rest_masses = to_save;
  Q_value = file_Q;

  // calculate the mass of the residual from the Q value
  double true_mRes;
  if( rest_masses->mProd == rest_masses->mProj )
  {
    // elastic or inelastic collision
    true_mRes = rest_masses->mTarg - Q_value;
  }
  else
  {
    // knock-on reaction
    true_mRes = ( rest_masses->mTarg - Q_value ) + 
      ( rest_masses->mProj - rest_masses->mProd );
  }
  // check the input mass of residual
  static double etol = Global.Value( "E_tol" );
  if( ( rest_masses->mRes > 0.0 ) && ( abs( true_mRes - rest_masses->mRes ) > etol*true_mRes ) )
  {
    Warning( "relativistic_masses::set_map", "mass of residual inconsistent with Q value" );
  }
  to_save->mRes = true_mRes;

  // total mass of all particles
  M_total = rest_masses->mTarg + rest_masses->mRes + 
       rest_masses->mProj + rest_masses->mProd;
}
// ------------- relativistic_masses::get_threshold --------------------
// Calculates the threshold energy
void relativistic_masses::get_threshold( )
{
  if( Q_value >= 0.0 )
  {
    threshold = 0.0;
  }
  else
  {
    threshold = - M_total * Q_value / ( 2.0*rest_masses->mTarg );
  }
}
// ------------- relativistic_masses::Newton_min_Eout --------------------
// For exothermic reactions, returns the Newtonian incident energy for minimal outgoing lab energy
double relativistic_masses::Newton_min_Eout( )
{
  // For mu = -1 find E_in such that dV_0/dE_in = dvcm_out/dE_in
  // parameters Newt_gamma and Newt_delta such that
  //  pcm^2 = Newt_gamma * ( Q + Newt_delta * E_in )
  double mass_in = rest_masses->mTarg + rest_masses->mProj;
  double beta_eject = rest_masses->mRes / ( rest_masses->mProd + rest_masses->mRes );
  double Newt_gamma = 2.0 * rest_masses->mProd * beta_eject;
  double mass_product = 2.0 * rest_masses->mProj * rest_masses->mProd * rest_masses->mProd;
  double E_min_in = mass_product * Q_value /
    ( beta_eject * ( Newt_gamma * beta_eject * mass_in * mass_in - mass_product ) );
  // don't go below the threshold
  if( E_min_in < threshold )
  {
    E_min_in = threshold;
  }
  return E_min_in;
}
// ------------- relativistic_masses::Newton_zero_Eout --------------------
// Returns the incident energy for zero outgoing lab-frame energy, Newtonian
double relativistic_masses::Newton_zero_Eout( )
{
  double mass_in = rest_masses->mTarg + rest_masses->mProj;
  double gamma = rest_masses->mProj * rest_masses->mProd /( mass_in * mass_in );
  double beta_eject = rest_masses->mRes / ( rest_masses->mProd + rest_masses->mRes );
  double alpha = beta_eject * rest_masses->mTarg / mass_in;

  if( gamma == alpha )
  {
    return 0.0;
  }
  else
  {
    double E_zero_in = -beta_eject*Q_value/( alpha - gamma );
    return E_zero_in;
  }
}

// ************************ clsss relativistic_param ******************
// ------------- relativistic_param::p_from_T --------------------
// Gets the momentum from the kinetic energy
double relativistic_param::p_from_T( double T, double E0 )
{
  double p = sqrt( T*( 2*E0 + T ) );
  return p;
}
// ------------- relativistic_param::T_from_p --------------------
// Gets the momentum from the kinetic energy
double relativistic_param::T_from_p( double p, double E0 )
{
  double p_sq = p*p;
  double T = p_sq/( E0 + sqrt( E0*E0 + p_sq ) );
  return T;
}
// ------------- relativistic_param::set_boost --------------------
// Sets up the boost to the lab frame
void relativistic_param::set_boost( double T_in_lab )
{
  Tin_lab = T_in_lab;
  pin_lab = p_from_T( Tin_lab, masses->rest_masses->mProj );
  double initial_rest_energy = masses->rest_masses->mTarg +
    masses->rest_masses->mProj;
  Minkowski = initial_rest_energy * initial_rest_energy +
    2.0 * T_in_lab * masses->rest_masses->mTarg;
  Minkowski = sqrt( Minkowski );
  pin_cm = pin_lab*masses->rest_masses->mTarg / Minkowski;
  sinh_chi = pin_cm/masses->rest_masses->mTarg;
  cosh_chi = sqrt( 1.0 + sinh_chi*sinh_chi );
  get_p_cm_out( );
}
// ------------- relativistic_param::boost --------------------
// Boosts the ejected particle from the center-of-mass to the lab frame
void relativistic_param::boost( double mu_cm, double *T_lab, double *p_lab_parallel )
{
  double E_cm = masses->rest_masses->mProd + T_cm_out;
  double p_cm_parallel = p_cm_out*mu_cm;
  double p_cm_perpendicular = p_cm_out*sqrt( 1.0 - mu_cm*mu_cm );
  *p_lab_parallel = sinh_chi*E_cm + cosh_chi*p_cm_parallel;
  double p_lab_sq = ( *p_lab_parallel )*( *p_lab_parallel ) +
    p_cm_perpendicular * p_cm_perpendicular;
  *T_lab = T_from_p( sqrt( p_lab_sq ), masses->rest_masses->mProd );
}
// ------------- relativistic_param::get_p_cm_out --------------------
// Calculates the center-of-mass energy and momentum for discrete 2-body reactions
void relativistic_param::get_p_cm_out( )
{
  double temp = masses->M_total*masses->Q_value + 2*masses->rest_masses->mTarg* Tin_lab;
  if( temp < 0.0 ) temp = 0.0;
  p_cm_out = sqrt( temp*( temp + 4*masses->rest_masses->mRes*masses->rest_masses->mProd ) ) /
    ( 2.0 * Minkowski );
  T_cm_out = T_from_p( p_cm_out, masses->rest_masses->mProd );
}
// ------------- relativistic_param::get_T_lab_out --------------------
// Calculates the lab-frame kinetic energy for discrete 2-body reactions
double relativistic_param::get_T_lab_out( double mu_cm )
{
  // boost to lab frame
  double p_lab_parallel;
  double Tout_lab;
  boost( mu_cm, &Tout_lab, &p_lab_parallel );
  return Tout_lab;
}
// ------------- relativistic_param::get_mu_cm --------------------
// Gets the value of mu_cm given the laboratory energy of outgoing particle.
double relativistic_param::get_mu_cm( double T_lab )
{
  // This routine uses the fact that the equation for mu_cm is a quadratic
  double p_lab_sq = T_lab*( 2*masses->rest_masses->mProd + T_lab );
  double A = p_cm_out*sinh_chi;
  A *= A;  // square it
  double E_prod = masses->rest_masses->mProd + T_cm_out;
  double B = 2*E_prod*p_cm_out*sinh_chi*cosh_chi;
  double C = E_prod*sinh_chi;
  C = C*C + p_cm_out*p_cm_out - p_lab_sq;
  double root_1;
  double root_2;
  int num_roots = math_F::quadratic( A, B, C, &root_1, &root_2 );
  if( num_roots != 2 )
  {
    FatalError( "relativistic_param::get_mu_cm", "wrong number of roots" );
  }

  static double mu_tol = Global.Value( "E_tol" );
  if( root_2 > 1.0 )
  {
    if( root_2 > 1.0 + mu_tol )
    {
      Warning ( "relativistic_param::get_mu_cm", pastenum( "mu too big: ", root_2 ) );
    }
    root_2 = 1.0;
  }
  if( root_2 < -1.0 )
  {
    if( ( root_2 < -1.0 - mu_tol ) && ( Tin_lab - masses->threshold < mu_tol*T_lab ) )
    {
      Warning ( "relativistic_param::get_mu_cm", pastenum( "mu too small: ", root_2 ) );
    }
    root_2 = -1.0;
  }
  return root_2;
}
// ------------- relativistic_param::get_E_mu_lab --------------------
// Gets the laboratory energy and cosine of outgoing particle
void relativistic_param::get_E_mu_lab( double Mu_cm, double *Tout_lab, double *mu_lab )
{
  mu_cm = Mu_cm;
  double p_cm_perpendicular = p_cm_out*sqrt( 1.0 - mu_cm*mu_cm );
  double p_lab_parallel;
  boost( Mu_cm, Tout_lab, &p_lab_parallel );
  *mu_lab = p_lab_parallel/sqrt( p_lab_parallel*p_lab_parallel +
    p_cm_perpendicular*p_cm_perpendicular );
}
// ------------- relativistic_param::zero_Eout --------------------
// Returns the incident energy for zero outgoing lab-frame energy, relativistic
double relativistic_param::zero_Eout( )
{
  // set up the void* parameters for the root finder
  void *params = static_cast< void * >( this );

  // get the incident energy for the Newtonian approximation (It's low)
  double T_in = masses->Newton_zero_Eout( );
  if( T_in <= 0.0 )
  {
    return 0.0;
  }

  mu_cm = -1.0;
  double p_out = relativistic_F::p_out_lab( T_in, params );
  dd_entry low_Tin( T_in, p_out );
  // a high bound
  T_in *= 1.1;
  p_out = relativistic_F::p_out_lab( T_in, params );
  dd_entry high_Tin( T_in, p_out );

  // for the root finder
  static double E_tol = Global.Value( "E_tol" );
  double this_Tin = math_F::zeroin( relativistic_F::p_out_lab, 0.0,
				    low_Tin, high_Tin, params, E_tol );
  return this_Tin;
}
// ------------- relativistic_param::find_hit --------------------
// Finds the incident energy for given T_lab and mucm.
// Returns the solutions.
double relativistic_param::find_hit( double T_lab, double mucm, const dd_entry &pair_0,
		   const dd_entry &pair_1 )
{
  // set up the void* parameters for the root finder
  mu_cm = mucm;
  void *params = static_cast< void * >( this );

  // use copies of the initial pairs
  dd_entry use_pair_0 = pair_0;
  dd_entry use_pair_1 = pair_1;
  static double E_tol = Global.Value( "E_tol" );
  double this_Tin = math_F::zeroin( relativistic_F::T_out_lab, T_lab,
				    use_pair_0, use_pair_1, params, E_tol );
  return this_Tin;
}
// ------------- relativistic_param::find_bottom --------------------
// For mu < 0 find the incident energy which minimizes Eout
// On exit pair_0 and pair_1 are tighter bounds.
dd_entry relativistic_param::find_bottom( double mu, dd_entry *pair_0,
  dd_entry *pair_1, double tol )
{
  if( ( masses->Q_value == 0.0 ) || ( mu >= 0.0 ) )
  {
    FatalError( "relativistic_param::find_bottom",
		"improper input" );
  }
  double Ein;
  double Eout;
  // parameters for the function relativistic_F::T_out_lab
  mu_cm = mu;
  void *params = static_cast< void *>( this );
  while( abs( pair_1->x - pair_0->x ) > tol * pair_1->x )
  {
    Ein = math_F::parabola_bottom( relativistic_F::T_out_lab, pair_0,
					  pair_1, params );
    // don't go below the threshold
    if( Ein < masses->threshold )
    {
      Ein = masses->threshold;
    }
    Eout = relativistic_F::T_out_lab( Ein, params );

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
// ------------- relativistic_param::find_lowest_bottom --------------------
// Find the incident energy which minimizes Eout for mu = -1
dd_entry relativistic_param::find_lowest_bottom( )
{
  static double etol = Global.Value( "E_tol" );

  dd_entry ans;
  // pairs for the root finder, zeroin
  dd_entry pair_0; // ( Ein, Eout ) for Ein below the minimum
  dd_entry pair_1; // ( Ein, Eout ) for Ein above the minimum

  // parameters for the function relativistic_F::T_out_lab
  void *params = static_cast< void *>( this );

  if( masses->Q_value == 0.0 )
  {
    ans.x = 0.0;
    ans.y = 0.0;
  }
  else
  {
    double Root = zero_Eout( );
    if( Root > 0.0 )
    {
      ans.x = Root;
      ans.y = 0.0;
    }
    else
    {
      // use the root finder
      pair_0.x = 0.0;
      pair_0.y = relativistic_F::T_out_lab( 0.0, params );
      pair_1.x = masses->Newton_min_Eout( );
      mu_cm = -1.0;
      pair_1.y = relativistic_F::T_out_lab( pair_1.x, params );
      ans = find_bottom( -1.0, &pair_0, &pair_1, etol );
    }
  }
  return ans;
}

// ***************** functions ********************
// ------------- relativistic_F::T_out_lab --------------------
// Returns the lab-frame kinetic energy of the outgoing particle
double relativistic_F::T_out_lab( double T_lab_in, void *params )
{
  // The parameters are really a relativistic_param *
  relativistic_param *map = static_cast< relativistic_param * >( params );
  // Sets up the boost to the lab frame
  map->set_boost( T_lab_in );
  double T = map->get_T_lab_out( map->mu_cm );
  return T;
}

// ------------- relativistic_F::p_out_lab --------------------
// Returns the lab-frame  parallel momentum of the outgoing particle
double relativistic_F::p_out_lab( double T_lab_in, void *params )
{
  // The parameters are really a relativistic_param *
  relativistic_param *map = static_cast< relativistic_param * >( params );
  // Sets up the boost to the lab frame
  map->set_boost( T_lab_in );
  double p_parallel;  // parallel momentum in lab frame
  double T_lab;  // lab frame kinetic energy
  map->boost( map->mu_cm, &T_lab, &p_parallel );
  return p_parallel;
}
