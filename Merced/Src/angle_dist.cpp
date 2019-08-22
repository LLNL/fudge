/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: angle_dist.cpp 1 2006-02-02 03:06:56Z hedstrom $
 *
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

// implementation of the classes used to handle angular distributions

#include <string>
#include <cmath>
#ifdef _OPENMP
 #include <omp.h>
#endif

#include "angle_dist.hpp"
#include "messaging.hpp"
#include "global_params.hpp"

// ********* class two_body_Ein_param *********
// ---------------- two_body_Ein_param::set_rel_map ------------------
// Sets up the relativistic mapping from center-of-mass to lab frame
void two_body_Ein_param::set_rel_map( relativistic_masses *map )
{
  relativistic_map.masses = map;
  lower_hits.relativistic_map.masses = map;
  upper_hits.relativistic_map.masses = map;
}
// ---------------- two_body_Ein_param::set_Newton_map ------------------
// Sets up the Newtonian mapping from center-of-mass to lab frame
void two_body_Ein_param::set_Newton_map( two_body_map *map )
{
  Newton_map.masses = map;
  lower_hits.Newton_map.masses = map;
  upper_hits.Newton_map.masses = map;
}
// ---------------- two_body_Ein_param::reset_start ------------------
// Initializes the pointers to the angular probabilities for this E_in
void two_body_Ein_param::reset_start( )
{
  left_ptr = this_mucm_dist->begin( );
  next_left_ptr = left_ptr;
  ++next_left_ptr;
  right_ptr = next_mucm_dist->begin( );
  next_right_ptr = right_ptr;
  ++next_right_ptr;

  // get the range of mu_cm values
  double lower_mu = ( left_ptr->x > right_ptr->x )?
    left_ptr->x : right_ptr->x;
  double higher_mu = ( next_left_ptr->x < next_right_ptr->x )?
    next_left_ptr->x : next_right_ptr->x;
  if( higher_mu <= lower_mu )
  {
    FatalError( "two_body_Ein_param::reset_start", "Check the mu values." );
  }

  // Interpolate to the common mu_cm values
  common_low_mucm( lower_mu );
  common_high_mucm( higher_mu );

  // set up upper_hits
  // The first step in the main loop in angle_dist::mucm_ladder 
  // will transfer this information to lower_hits.
  set_upper_hits( lower_mu );
}
// ---------------- two_body_Ein_param::set_upper_hits ------------------
// Initializes upper_hits for the value of mucm
void two_body_Ein_param::set_upper_hits( double mucm )
{
  list< min_Eout_info >::const_iterator min_Eout_ptr;
  upper_hits.left_Ein_Eout.x = left_data_Ein;
  upper_hits.right_Ein_Eout.x = right_data_Ein;
  void *params;  // for the function parameters
  if( use_relativistic )
  {
    relativistic_map.mu_cm = mucm;
    params = static_cast< void * >( &relativistic_map );
    upper_hits.left_Ein_Eout.y = relativistic_F::T_out_lab( left_data_Ein, params );
    upper_hits.right_Ein_Eout.y = relativistic_F::T_out_lab( right_data_Ein, params );
    // for the minimal outgoing energy
    if( ( mucm >= 0.0 ) || ( relativistic_map.masses->Q_value == 0.0 ) )
    {
      upper_hits.flip.x = threshold;
      upper_hits.flip.y = threshold_out;
    }
    else
    {
      min_Eout_ptr = get_min_Eout_info( mucm );
      upper_hits.flip.x = min_Eout_ptr->Ein;
      upper_hits.flip.y = min_Eout_ptr->Eout;
    }
  }
  else
  {
    upper_hits.Newton_map.mu_cm = mucm;
    params = static_cast< void * >( &upper_hits.Newton_map );
    upper_hits.left_Ein_Eout.y = Newtonian_F::T_out_lab( left_data_Ein, params );
    upper_hits.right_Ein_Eout.y = Newtonian_F::T_out_lab( right_data_Ein, params );
    // for the minimal outgoing energy
    if( ( mucm >= 0.0 ) || ( Newton_map.masses->Q_value == 0.0 ) )
    {
      upper_hits.flip.x = threshold;
      upper_hits.flip.y = threshold_out;
    }
    else
    {
      min_Eout_ptr = get_min_Eout_info( mucm );
      upper_hits.flip.x = min_Eout_ptr->Ein;
      upper_hits.flip.y = min_Eout_ptr->Eout;
    }
  }
}
// ---------------- two_body_Ein_param::reset_hits ------------------
// Initializes lower_hits and upper_hits for the values of mucm
void two_body_Ein_param::reset_hits( )
{
  // move upper_hits data to lower_hits
  lower_hits.flip = upper_hits.flip;
  lower_hits.left_Ein_Eout = upper_hits.left_Ein_Eout;
  lower_hits.right_Ein_Eout = upper_hits.right_Ein_Eout;
  lower_hits.Newton_map.mu_cm = upper_hits.Newton_map.mu_cm;
  // set up upper_hits
  set_upper_hits( next_left_data.x );

  // set the range of outgoing energies
  if( upper_hits.flip.x <= upper_hits.left_Ein_Eout.x )
  {
    upper_data_max_Eout = upper_hits.right_Ein_Eout.y;
  }
  else if( upper_hits.flip.x >= upper_hits.right_Ein_Eout.x )
  {
    upper_data_max_Eout = upper_hits.left_Ein_Eout.y;
  }
  else
  {
    upper_data_max_Eout =
      ( upper_hits.left_Ein_Eout.y > upper_hits.right_Ein_Eout.y ) ?
      upper_hits.left_Ein_Eout.y : upper_hits.right_Ein_Eout.y;
  }

  if( lower_hits.flip.x <= lower_hits.left_Ein_Eout.x )
  {
    lower_data_min_Eout = lower_hits.left_Ein_Eout.y;
  }
  else if( lower_hits.flip.x >= lower_hits.right_Ein_Eout.x )
  {
    lower_data_min_Eout = lower_hits.right_Ein_Eout.y;
  }
  else
  {
    lower_data_min_Eout = lower_hits.flip.y;
  }
}
// ---------------- two_body_Ein_param::next_mucm ------------------
// Sets up the next interval of mu_cm values
bool two_body_Ein_param::next_mucm( )
{
  bool done = false;

  left_data = next_left_data;
  right_data = next_right_data;

  // update the pointers
  if( next_left_ptr->x <= left_data.x )
  {
    left_ptr = next_left_ptr;
    ++next_left_ptr;
    if( next_left_ptr == this_mucm_dist->end( ) )
    {
      return true;
    }
  }
  if( next_right_ptr->x <= right_data.x )
  {
    right_ptr = next_right_ptr;
    ++next_right_ptr;
    if( next_right_ptr == next_mucm_dist->end( ) )
    {
      return true;
    }
  }

  // set the common upper mu_cm value
  double upper_mu = ( next_left_ptr->x < next_right_ptr->x ) ?
    next_left_ptr->x : next_right_ptr->x;
  common_high_mucm( upper_mu );
  return done;
}
// ---------------- two_body_Ein_param::common_low_mucm ------------------
// Interpolates (mu_cm, probability) data to the lower common mu_cm value
void two_body_Ein_param::common_low_mucm( double lower_mu )
{
  if( left_ptr->x == lower_mu )
  {
    left_data = *left_ptr;
  }
  else
  {
    left_data.x = lower_mu;
    left_data.y = left_ptr->linlin_interp( lower_mu, *next_left_ptr );
  }

  if( right_ptr->x == lower_mu )
  {
    right_data = *right_ptr;
  }
  else
  {
    right_data.x = lower_mu;
    right_data.y = right_ptr->linlin_interp( lower_mu, *next_right_ptr );
  }
}
// ---------------- two_body_Ein_param::common_high_mucm ------------------
// Interpolates (mu_cm, probability) data to the higher common mu_cm value
void two_body_Ein_param::common_high_mucm( double higher_mu )
{
  if( next_left_ptr->x == higher_mu )
  {
    next_left_data = *next_left_ptr;
  }
  else
  {
    next_left_data.x = higher_mu;
    next_left_data.y = left_ptr->linlin_interp( higher_mu, *next_left_ptr );
  }

  if( next_right_ptr->x == higher_mu )
  {
    next_right_data = *next_right_ptr;
  }
  else
  {
    next_right_data.x = higher_mu;
    next_right_data.y = right_ptr->linlin_interp( higher_mu, *next_right_ptr );
  }
}
// ----------- two_body_Ein_param::get_min_Eout_info --------------
// Returns the min_Eout_info for given center-of-mass outgoing cosine
min_Eout_info_list::const_iterator two_body_Ein_param::get_min_Eout_info( double mu_cm ) const
{
  min_Eout_info_list::const_iterator min_Eout_info_ptr = min_Eout_list->begin( );
  while( min_Eout_info_ptr->mu < mu_cm )
  {
    ++min_Eout_info_ptr;
    if( min_Eout_info_ptr == min_Eout_list->end( ) )
    {
      FatalError( "two_body_Ein_param::get_min_Eout_info", "min_Eout_info not found" );
    }
  }
  return min_Eout_info_ptr;
}
// ********* class two_body_mucm_param *********

// ********* class angle_hit_list *********
// Calculates E_out from E_in, for testing the side of a box
// ----------- angle_hit_list::get_Eout --------------
double angle_hit_list::get_Eout( double E_in )
{
  double Eout;
  void *params;  // parameters for T_out_lab
  if( use_relativistic )
  {
    relativistic_map.mu_cm = eta;
    params = static_cast< void * >( &relativistic_map );
    Eout = relativistic_F::T_out_lab( E_in, params );
  }
  else
  {
    Newton_map.mu_cm = eta;
    params = static_cast< void * >( &Newton_map );
    Eout = Newtonian_F::T_out_lab( E_in, params );
  }

  return Eout;
}
// ----------- angle_hit_list::find_hit --------------
// Finds an intersection with the bottom or top of a box
double angle_hit_list::find_hit( double E_out, const dd_entry &pair_0,
				 const dd_entry &pair_1 )
{
  if( ( pair_0.y - E_out ) * ( pair_1.y - E_out ) > 0.0 )
  {
    return -1.0;
    // FatalError( "angle_hit_list::find_hit", "no root" );
  }

  double root;
  if( use_relativistic )
  {
    root = relativistic_map.find_hit( E_out, eta, pair_0, pair_1 );
  }
  else
  {
    root = Newton_map.find_hit( E_out, eta, pair_0, pair_1 );
  }
  return root;
}
// ----------- angle_hit_list::find_bottom_hits --------------
// Finds the intersections with the bottom of a box
void angle_hit_list::find_bottom_hits( double E_out,
  vector< Ein_Eta_Hit > *Ein_hits )
{
  // for new entries
  Ein_Eta_Hit Ein_mucm_hit;

  // where is the minimum?
  if( flip.x <= left_Ein_Eout.x )
  {
    // E_lab is increasing with E_in
    // append this entry
    Ein_mucm_hit.E_in = find_hit( E_out, left_Ein_Eout, right_Ein_Eout );
    Ein_mucm_hit.hit_edge = BOTTOM_IN;
    Ein_hits->push_back( Ein_mucm_hit );
  }
  else if( flip.x >= right_Ein_Eout.x )
  {
    // E_lab is decreasing with E_in
    // append this entry
    Ein_mucm_hit.E_in = find_hit( E_out, left_Ein_Eout, right_Ein_Eout );
    Ein_mucm_hit.hit_edge = BOTTOM_OUT;
    Ein_hits->push_back( Ein_mucm_hit );
  }
  else if( flip.y < E_out ) // omit a tangent contact
  {
    // We are here only for Q < 0
    // tolerance for short intervals
    static double etol = Global.Value( "E_tol" );
    double slop = etol*( right_Ein_Eout.x - left_Ein_Eout.x );
    // E_lab is at first decreasing with E_in
    // append this entry
    if( flip.x >= left_Ein_Eout.x + slop )
    {
      Ein_mucm_hit.E_in = find_hit( E_out, left_Ein_Eout, flip );
      Ein_mucm_hit.hit_edge = BOTTOM_OUT;
      Ein_hits->push_back( Ein_mucm_hit );
    }
    // E_lab is then increasing with E_in
    // append this entry
    if( flip.x <= right_Ein_Eout.x - slop )
    {
      Ein_mucm_hit.E_in = find_hit( E_out, flip, right_Ein_Eout );
      Ein_mucm_hit.hit_edge = BOTTOM_IN;
      Ein_hits->push_back( Ein_mucm_hit );
    }
  }
}
// ----------- angle_hit_list::find_top_hits --------------
// Finds the intersections with the top of a box
void angle_hit_list::find_top_hits( double E_out,
  vector< Ein_Eta_Hit > *Ein_hits )
{
  // treat it like the bottom of the box
  find_bottom_hits( E_out, Ein_hits );
  for( vector< Ein_Eta_Hit >::iterator this_hit = Ein_hits->begin( );
       this_hit != Ein_hits->end( ); ++this_hit )
  {
    if( this_hit->hit_edge == BOTTOM_OUT )
    {
      this_hit->hit_edge = TOP_IN;
    }
    else if( this_hit->hit_edge == BOTTOM_IN )
    {
      this_hit->hit_edge = TOP_OUT;
    }
  }
}

// ********* class angle_dist *********
// ----------- angle_dist::setup_map --------------
void angle_dist::setup_map( )
// set up the map from center-of-mass to laboratory coordinates
{
  // function parameter
  void *params;
  relativistic_mass.setup_masses( &particle_info, Q );

  if( use_relativistic )
  {
    relativistic_mass.setup_masses( &particle_info, Q );
    relativistic_mass.get_threshold( );
    threshold = relativistic_mass.threshold;

    // parameters for the relativistic_F functions
     relativistic_param relativistic_map;
    relativistic_map.masses = &relativistic_mass;
    params = static_cast< void * >( &relativistic_map );
    // we need to set the direction cosine
    relativistic_map.mu_cm = 0.0;
    threshold_out = relativistic_F::T_out_lab( threshold, params );
  }
  else
  {
    map.set_map( particle_info, Q );
    map.get_threshold( );
    threshold = map.threshold;

    // parameters for Newtonian_F::T_out_lab
    Newton_map_param Newton_map;
    Newton_map.masses = &map;
    params = static_cast< void * >( &Newton_map );
    threshold_out = Newtonian_F::T_out_lab( threshold, params );
  }
}
// ----------- angle_dist::init_min_Eout_list --------------
// Initializes min_Eout_list, used in the relativistic treatment of endothermic reactions
void angle_dist::init_min_Eout_list( )
{
  list< double > mu_list;
  list< double > add_list;
  angle_dist::const_iterator Ein_ptr = begin( );
  dd_vector::const_iterator mu_ptr;

  // Make the first list
  for( mu_ptr = Ein_ptr->begin( ); mu_ptr != Ein_ptr->end( ); ++mu_ptr )
  {
    if( mu_ptr->x >= 0 )
    {
      break;
    }
    mu_list.push_back( mu_ptr->x );
  }

  // Loop over the other incident energies
  for( ++Ein_ptr; Ein_ptr != end( ); ++Ein_ptr )
  {
    for( mu_ptr = Ein_ptr->begin( ); mu_ptr != Ein_ptr->end( ); ++mu_ptr )
    {
      if( mu_ptr->x >= 0 )
      {
        break;
      }
      add_list.push_back( mu_ptr->x );
    }
    mu_list.merge( add_list );
    mu_list.unique( );
  }
  // initialize min_Eout_list
  min_Eout_info min_Eout_item;
  list< double >::const_iterator muPtr;
  for( muPtr = mu_list.begin( ); muPtr != mu_list.end( ); ++muPtr )
  {
    min_Eout_list.push_back( min_Eout_item );
    min_Eout_info_list::iterator list_ptr = min_Eout_list.end( );
    --list_ptr;
    list_ptr->mu = *muPtr;
  }
}
// ----------- angle_dist::make_min_Eout_list --------------
// Computes min_Eout_list, used in the relativistic treatment of endothermic reactions
void angle_dist::make_min_Eout_list( )
{
  static double etol = Global.Value( "E_tol" );

  // this is not for elastic scattering
  if( Q == 0.0 )
  {
    return;
  }
  // First, set up the list
  init_min_Eout_list( );

  // pairs for the root finder, zeroin
  dd_entry pair_0( threshold, threshold_out ); // ( Ein, Eout ) for Ein below the minimum
  dd_entry pair_1; // ( Ein, Eout ) for Ein above the minimum
  dd_entry ans;

  // parameters for T_out_lab
  Newton_map_param Newton_map;
  relativistic_param relativistic_map;
  void *params;

  // use the threshold for the lower bound
  if( use_relativistic )
  {
    relativistic_map.masses = &relativistic_mass;
    ans = relativistic_map.find_lowest_bottom( );
    params = static_cast< void * >( &relativistic_map );
  }
  else // Newtonian
  {
    Newton_map.masses = &map;
    ans = Newton_map.find_lowest_bottom( );
    params = static_cast< void * >( &Newton_map );
  }

  min_Eout_info_list::iterator mu_ptr = min_Eout_list.begin( );
  // special for mu = -1
  if( mu_ptr->mu == -1.0 )
  {
    mu_ptr->Ein = ans.x;
    mu_ptr->Eout = ans.y;
    ++mu_ptr;
  }

  // do the rest of the cosines
  for( ; mu_ptr != min_Eout_list.end( ); ++mu_ptr )
  {
    pair_0.x = threshold;
    pair_0.y = threshold_out;
    pair_1.x = ans.x;
    if( use_relativistic )
    {
      relativistic_map.mu_cm = mu_ptr->mu;
      pair_1.y = relativistic_F::T_out_lab( pair_1.x, params );
      ans = relativistic_map.find_bottom( mu_ptr->mu, &pair_0, &pair_1, etol );
    }
    else
    {
      Newton_map.mu_cm = mu_ptr->mu;
      pair_1.y = Newtonian_F::T_out_lab( pair_1.x, params );
      ans = Newton_map.find_bottom( mu_ptr->mu, &pair_0, &pair_1, etol );
    }
    mu_ptr->Ein = ans.x;
    mu_ptr->Eout = ans.y;
  }
}
// ----------- angle_dist::set_threshold --------------
// Uses the mass difference to set the threshold
void angle_dist::set_threshold( )
{
  /*
  if( use_relativistic )
  {
    threshold = relativistic_mass.threshold;
  }
  else
  {
    threshold = map.threshold;
  }
  */

  // adjust the data if necessary
  angle_dist::iterator data_ptr = begin( );
  double first_Ein = data_ptr->get_E_in( );
  if( first_Ein < threshold )
  {
    data_ptr->set_E_in( threshold );
  }
}
// ----------- angle_dist::get_Ein_range --------------
//  Gets the range of nontrivial incident energy bins; computes first_Ein and last_Ein
// returns true if the threshold is too high for the energy bins
bool angle_dist::get_Ein_range( const dd_vector& sigma, const dd_vector& mult,
    const dd_vector& weight,
    const Flux_List& e_flux, const Energy_groups& Ein_groups )
{
  double E_first;
  double E_last;
  two_body_Ein_param initial_param;
  bool done = initial_param.get_Ein_range( sigma, mult, weight, e_flux,
					 Ein_groups, &E_first, &E_last );
  if( done ) return true;

  // check the range of incident energies for the probability data
  angle_dist::const_iterator mucm_ptr = begin( );
  double E_data = mucm_ptr->get_E_in( );
  if( E_data > E_first ) 
  {
    E_first = E_data;
  }
  first_Ein = Ein_groups.first_bin_ID( E_first );

  mucm_ptr = end( );
  --mucm_ptr;
  E_data = mucm_ptr->get_E_in( );
  if( E_data < E_last ) 
  {
    E_last = E_data;
  }
  last_Ein = Ein_groups.last_bin_ID( E_last );

  return false;
}
// ----------- angle_dist::initialize_param --------------
// Initializes the quadrature parameters
void angle_dist::initialize_param( Quadrature_Method mu_quad_method, 
  two_body_Ein_param *Ein_param )
{
  Ein_param->set_Newton_map( &map );
  Ein_param->set_rel_map( &relativistic_mass );
  Ein_param->use_relativistic = use_relativistic;
  Ein_param->lower_hits.use_relativistic = use_relativistic;
  Ein_param->upper_hits.use_relativistic = use_relativistic;
  Ein_param->threshold = threshold;
  Ein_param->threshold_out = threshold_out;
  Ein_param->mu_quad_method = mu_quad_method;
  Ein_param->min_Eout_list = &min_Eout_list;
}
// ----------- angle_dist::setup_data --------------
// Initializes the quadrature parameters
void angle_dist::setup_data( two_body_Ein_param *Ein_param )
{
  Ein_param->this_mucm_dist = begin( );
  Ein_param->next_mucm_dist = Ein_param->this_mucm_dist;
  ++Ein_param->next_mucm_dist;
  while( Ein_param->next_mucm_dist->get_E_in( ) <= *Ein_param->Ein_ptr )
  {
    Ein_param->this_mucm_dist = Ein_param->next_mucm_dist;
    ++Ein_param->next_mucm_dist;
  }

  double first_E = Ein_param->this_mucm_dist->get_E_in( );
  if( first_E > Ein_param->data_E_0 )
  {
    // because of the ENDL kludge, this should never happen
    Ein_param->data_E_0 = first_E;
  }
  Ein_param->left_data_Ein = Ein_param->this_mucm_dist->get_E_in( );
  Ein_param->right_data_Ein = Ein_param->next_mucm_dist->get_E_in( );
}
// ----------- angle_dist::set_Ein_range --------------
// Sets the range of incident energies for this intergration
void angle_dist::set_Ein_range( int Ein_bin, two_body_Ein_param &Ein_param )
{
  Ein_param.set_Ein_range( );
  double this_E = Ein_param.this_mucm_dist->get_E_in( );
  if( this_E > Ein_param.data_E_0 ) Ein_param.data_E_0 = this_E;
  this_E = Ein_param.next_mucm_dist->get_E_in( );
  if( this_E < Ein_param.data_E_1 ) Ein_param.data_E_1 = this_E;

  // the data may be below the threshold
  if( Ein_param.data_E_1 < threshold )
  {
    Ein_param.data_E_0 = threshold;
    Ein_param.data_E_1 = threshold;
    // Warning( "angle_dist::set_Ein_range", "ignoring 2 data below threshold" );
  }
  else if( Ein_param.data_E_0 < threshold )
  {
    Ein_param.data_E_0 = threshold;
    // Warning( "angle_dist::set_Ein_range", "ignoring data below threshold" );
  }

  if( Ein_param.data_E_1 < Ein_param.data_E_0 )
  {
    FatalError( "angle_dist::set_Ein_range", "check the I=1 incident energies" );
  }
  Ein_param.set_sigma_range( );
}
// ----------- angle_dist::next_ladder --------------
// go to the next interval
bool angle_dist::next_ladder( double E_in, two_body_Ein_param *Ein_param )
{
  bool done = Ein_param->update_bin_pointers( E_in );
  static double etol = Global.Value( "E_tol" );
  if( !done )
  {
    double E_tol = E_in * etol;
    if( E_in + E_tol >= Ein_param->next_mucm_dist->get_E_in( ) )
    {
      while( E_in + E_tol >= Ein_param->next_mucm_dist->get_E_in( ) )
      {
        // get the next angular data
        Ein_param->this_mucm_dist = Ein_param->next_mucm_dist;
        ++Ein_param->next_mucm_dist;
        if( Ein_param->next_mucm_dist == end ( ) )
        {
          return true;
        }
      }
      Ein_param->left_data_Ein = Ein_param->this_mucm_dist->get_E_in( );
      Ein_param->right_data_Ein = Ein_param->next_mucm_dist->get_E_in( );
    }
  }
  return done;
}
// ----------- angle_dist::read_data --------------
void angle_dist::read_data( data_parser &input_file, int num_Ein )
{
  // Read the interpolation rules
  interp_flag_F::read_2d_interpolation( input_file, &Ein_interp, &mu_interp );

  dd_vector new_angle_dist;  // angular distribution for one E_in
  angle_dist::iterator new_angle_ptr;

  // read the data
  for( int Ein_count = 0; Ein_count < num_Ein; ++Ein_count )
  {
    // insert a new angular distribution
    new_angle_ptr = insert( end( ), dd_vector( ) );
    // get the incident energy and the data pairs
    new_angle_ptr->set_E_in( input_file.get_next_double( ) );
    new_angle_ptr->interp_type = mu_interp;
    int num_mu = input_file.get_next_int( );
    new_angle_ptr->read_data( input_file, num_mu );
  }
  // we may need to copy the first distribution at the threshold
  ENDL_kludge( );
  //  print( );
}
// ----------- angle_dist::ENDL_kludge --------------
void angle_dist::ENDL_kludge( )
// ENDL has the convention that the distribution is assumed
// independent of energy at low incident energies.
{
  // Make sure that the threshold has been set
  if( threshold < 0.0 )
  {
    FatalError( "ENDL_kludge", "You need to set the threshold." );
  }
  angle_dist::iterator first_dist = begin( );
  if( first_dist->get_E_in( ) > threshold )
  {
    // prepend a copy at the threshold energy
    angle_dist::iterator second_dist = begin( );
    first_dist = insert( begin( ), dd_vector( ) );
    first_dist->copy( *second_dist );
    first_dist->set_E_in( threshold );
  }
}
// ----------- angle_dist::print --------------
void angle_dist::print( )
{
  for( angle_dist::iterator angle_ptr = begin( );
       angle_ptr != end( ); ++angle_ptr )
  {
    angle_ptr->print( );
  }
}
// ----------- angle_dist::get_T --------------
// Calculates the transfer matrix for this particle.
// sigma is the cross section.
void angle_dist::get_T( const dd_vector& sigma, dd_vector& weight, T_matrix& transfer )
{
  // set the flag for relativistic mechanics
  if( Global.Flag( "kinetics" ) == "newtonian" )
  {
    use_relativistic = false;
  }
  else
  {
    use_relativistic = true;
  }

  if( ( Ein_interp.flag != LINLIN ) || ( Ein_interp.qualifier != DIRECT ) )
  {
    FatalError( "angle_dist::get_T", "Incident energy interpolation not implemented" );
  }
  if( mu_interp != LINLIN )
  {
    FatalError( "angle_dist::get_T", "cosine interpolation not implemented" );
  }
  if( particle_info.mProd == 0.0 )
  {
    FatalError( "angle_dist::get_T", "gamma emission not implemented" );
  }

  // the multiplicity is one
  dd_vector multiple;
  multiple.make_flat( sigma, 1.0 );
  bool done = get_Ein_range( sigma, multiple, weight, transfer.e_flux,
    transfer.in_groups );
  if( done )
  {
    // *** Bail out, all angular data too high ***
    transfer.zero_transfer( );
  }

  // for center-of-mass data
  setup_map( );
  // use the threshold computed from the kinetics
  set_threshold( ); 
  dd_vector::const_iterator last_sigma = sigma.end( );
  --last_sigma;
  if( threshold >= last_sigma->x )
  {
    // *** Bail out, all angular data too high ***
    Warning( "angle_dist::get_T", "check the Q value" );
    transfer.zero_transfer( );
  }

  // for inelastic reactions find the minimal outgoing energies for backward emission
  make_min_Eout_list( );

  long int quad_count = 0;  // number of 2-d quadratures
  long int Ein_F_count= 0;  // number of calls to angle_dist_F::E_quad_F
  long int mu_F_count = 0;  // number of calls to angle_dist_F::mu_cm_quad_F

  // do the integrals incident bin by incident bin
#pragma omp parallel for schedule( dynamic, 1 ) default( none )	\
  shared( sigma, multiple, weight, transfer ) \
  reduction( +: quad_count ) reduction( +: Ein_F_count ) reduction( +: mu_F_count )
  for( int Ein_bin = first_Ein; Ein_bin < last_Ein; ++Ein_bin )
  {
    two_body_Ein_param Ein_param;
    //    cout << "bin count " << Ein_bin << endl;
    initialize_param( transfer.mu_quad_method, &Ein_param );

   // set up the data range for this bin
    Ein_param.setup_bin( Ein_bin, sigma, multiple, weight, transfer.e_flux,
			 transfer.in_groups );
    setup_data( &Ein_param );
    // work on this bin
    for( ; ; )
    {
      // get the incident energy interval common to all data
      set_Ein_range( Ein_bin, Ein_param );
      mucm_ladder( transfer, &Ein_param );
      // go to the next interval
      bool Done = next_ladder( Ein_param.data_E_1, &Ein_param );
      if( Done )
      {
        break;
      }
    }
    quad_count += Ein_param.quad_count;
    Ein_F_count += Ein_param.Ein_F_count;
    mu_F_count += Ein_param.mu_F_count;
  } // end of parallel loop

  // print the counts of function evaluations
  cout << "2d quadratures: " << quad_count << endl;
  cout << "angle_dist_F::E_quad_F calls: " << Ein_F_count << endl;
  cout << "angle_dist_F::mu_cm_quad_F calls: " << mu_F_count << endl;
  cout << "average angle_dist_F::E_quad_F calls: " << 1.0*Ein_F_count/quad_count << endl;
  cout << "average angle_dist_F::mu_cm_quad_F calls: " << 1.0*mu_F_count/Ein_F_count << endl;
}
// ----------- angle_dist::mucm_ladder --------------
// This routine uses the angular distributions this_mucm_dist and the
// next to calculate the contribution to the E_out boxes of the
// transfer matrix between incident energies Ein_param->data_E_0 and
// Ein_param->data_E_1.
void angle_dist::mucm_ladder( T_matrix& transfer,
  two_body_Ein_param *Ein_param )
{
  bool geom_OK;  // for checking the consistency of the geometry
  bool done = false;
  angle_hit_list test_hits;
  test_hits.Newton_map.masses = &map;
  test_hits.use_relativistic = use_relativistic;
  test_hits.relativistic_map.masses = &relativistic_mass;
  // loop through the angular data
  Ein_param->reset_start( );
  for( ; !done; done = Ein_param->next_mucm( ) )
  {
    // initialize lower_hits and upper_hits
    Ein_param->reset_hits( );
    // loop through the outgoing energies (column of transfer)
    for( int Eout_count = 0; Eout_count < transfer.num_Eout_bins;
      ++Eout_count )
    {
      vector< double >::const_iterator Eout_ptr = transfer.out_groups.begin( )
        + Eout_count;
      // the box may be below the data
      vector< double >::const_iterator next_Eout = Eout_ptr;
      ++next_Eout;
      if( ( Eout_count < transfer.num_Eout_bins - 1 ) &&
	  ( Ein_param->lower_data_min_Eout >= *next_Eout ) )
      {
	// go on to the next E-E' box
	continue;
      }
      // the box may be above the data
      if( ( Eout_count > 0 ) && ( Ein_param->upper_data_max_Eout <= *Eout_ptr ) )
      {
	// we are done with this pair of mucm values
	break;
      }
      // how does the lowest mucm = const curve meet this E-E' box?
      geom_OK = Ein_param->lower_hits.hit_box( Ein_param->left_data.x, Eout_ptr,
        Ein_param->data_E_0, Ein_param->data_E_1 );
      if( !geom_OK )
      {
	test_hits.eta = Ein_param->lower_hits.eta;
        test_hits.hit_box( Ein_param->left_data.x, Eout_ptr,
          Ein_param->data_E_0, Ein_param->data_E_1 );
        test_hits.print( );
        FatalError( "angle_dist::mucm_ladder", "Check the coding, 1" );
      }
      // how does the next mucm = const curve meet this E-E' box?
      geom_OK = Ein_param->upper_hits.hit_box( Ein_param->next_left_data.x, Eout_ptr,
        Ein_param->data_E_0, Ein_param->data_E_1 );
      if( !geom_OK )
      {
	test_hits.eta = Ein_param->upper_hits.eta;
        test_hits.hit_box( Ein_param->next_left_data.x, Eout_ptr,
			   Ein_param->data_E_0, Ein_param->data_E_1 );
        test_hits.print( );
        FatalError( "angle_dist::mucm_ladder", "Check the coding, 2" );
      }
      // integrate over this E-E' box
      one_Ebox( transfer, Eout_count, Ein_param );
    }
  }
}
// ----------- angle_dist::one_Ebox --------------
// Integrate over one E-E' box
void angle_dist::one_Ebox( T_matrix& transfer, int Eout_count,
  two_body_Ein_param *Ein_param )
{
  // the E' energy range
  Ein_param->Eout_min = transfer.out_groups[ Eout_count ];
  Ein_param->Eout_max = transfer.out_groups[ Eout_count + 1 ];

  // set up common incident energies
  Ein_param->lower_hits.common_hits( Ein_param->upper_hits );

  // integrate depending on how the curves mucm = const meet the box
  angle_hit_list::iterator low_hit_ptr = Ein_param->lower_hits.begin( );
  angle_hit_list::iterator next_low_ptr = low_hit_ptr;
  ++next_low_ptr;
  angle_hit_list::iterator high_hit_ptr = Ein_param->upper_hits.begin( );
  angle_hit_list::iterator next_high_ptr = high_hit_ptr;
  ++next_high_ptr;
  for( ; ( next_low_ptr != Ein_param->lower_hits.end( ) ) &&
	 ( next_high_ptr != Ein_param->upper_hits.end( ) );
       low_hit_ptr = next_low_ptr, ++next_low_ptr,
	 high_hit_ptr = next_high_ptr, ++next_high_ptr )
  {
    if( ( low_hit_ptr->hit_edge == ABOVE ) ||
	( low_hit_ptr->hit_edge == TOP_OUT ) )
    {
      // do nothing---we are above the E-E' box
      continue;
    }
    else if( ( low_hit_ptr->hit_edge == TOP_IN ) ||
             ( low_hit_ptr->hit_edge == BOTTOM_IN ) ||
	     ( low_hit_ptr->hit_edge == INSIDE ) )
    {
      // the lower mucm = const curve is inside the E-E' box
      Ein_param->use_Eout_min = false;
      // where is the upper curve?
      if( ( high_hit_ptr->hit_edge == ABOVE ) ||
          ( high_hit_ptr->hit_edge == TOP_OUT ) )
      {
	// integrate up to the top of the E-E' bin
        Ein_param->use_Eout_max = true;
      }
      else
      {
	// integrate up to the next mucm = const curve
        Ein_param->use_Eout_max = false;
      }
    }
    else
    {
      // the lower mucm = const curve is below the E-E' box;
      // integrate from Eout_min
      Ein_param->use_Eout_min = true;
      // where is the upper mucm = const curve?
      if( ( high_hit_ptr->hit_edge == BOTTOM_OUT ) ||
	  ( high_hit_ptr->hit_edge == BELOW ) )
      {
        // do nothing---we are below the E-E' box
	continue;
      }
      else if( ( high_hit_ptr->hit_edge == TOP_IN ) ||
               ( high_hit_ptr->hit_edge == BOTTOM_IN ) ||
	       ( high_hit_ptr->hit_edge == INSIDE ) )
      {
        // the upper mucm = const curve is inside the E-E' box
        Ein_param->use_Eout_max = false;
      }
      else
      {
        // the upper mucm = const curve is above the E-E' box
        Ein_param->use_Eout_max = true;
      }
    }
    // the range of integration in incident energy
    Ein_param->Ein_0 = low_hit_ptr->E_in;
    Ein_param->Ein_1 = next_low_ptr->E_in;
    update_T( transfer, Eout_count, Ein_param );
  }
}
// ----------- angle_dist::update_T --------------
void angle_dist::update_T( T_matrix &transfer, int Eout_count,
  two_body_Ein_param *Ein_param )
{
  static double tol = Global.Value( "quad_tol" );
  // a vector to store the integrals
  coef_vector value( transfer.order, transfer.conserve );
  value.set_zero( );

  // parameters for the integration
  QuadParamBase *params = static_cast< QuadParamBase* >( Ein_param );

  // loop over the cross section data
  Ein_param->this_sigma = Ein_param->first_ladder_sigma;
  Ein_param->next_sigma = Ein_param->this_sigma;
  ++Ein_param->next_sigma;
  // Ein_param->Ein_0 may be past Ein_param->next_sigma
  while( ( Ein_param->this_sigma != Ein_param->last_ladder_sigma ) &&
         ( Ein_param->next_sigma->x < Ein_param->Ein_0 ) )
  {
    Ein_param->this_sigma = Ein_param->next_sigma;
    ++Ein_param->next_sigma;
  }
  for( ; ( Ein_param->this_sigma != Ein_param->last_ladder_sigma ) &&
         ( Ein_param->this_sigma->x <  Ein_param->Ein_1 );
       Ein_param->this_sigma = Ein_param->next_sigma, ++Ein_param->next_sigma )
  {
    double left_E = ( Ein_param->this_sigma->x < Ein_param->Ein_0 ) ?
      Ein_param->Ein_0 : Ein_param->this_sigma->x;
    double right_E = ( Ein_param->next_sigma->x > Ein_param->Ein_1 ) ?
      Ein_param->Ein_1 : Ein_param->next_sigma->x;
    // evaluate the integral
    quad_F::integrate( angle_dist_F::E_quad_F, transfer.Ein_quad_method,
      left_E, right_E, params, tol, &value );

    // add this integral
    transfer( Ein_param->Ein_count, Eout_count ) += value;
    // increment the function counts
    Ein_param->Ein_F_count += Ein_param->func_count;
    ++Ein_param->quad_count;
  }
}
// ----------- angle_dist::isotropic --------------
// Checks whether an angular probability density is isotropic
bool angle_dist::isotropic( )
{
  bool iso = true;
  for( angle_dist::const_iterator angle_ptr = begin( );
       angle_ptr != end( ); ++angle_ptr )
  {
    if( !angle_ptr->isotropic( ) )
    {
      iso = false;
      break;
    }
  }
  return iso;
}

// **************** functions to integrate **********
// Function for the 1-d quadrature
// ---------------- angle_dist_F::mu_cm_quad_F ------------------
void angle_dist_F::mu_cm_quad_F( double mucm, QuadParamBase *mucm_quad_param,
   coef_vector *value )
{
  // the parameters are really two_body_mucm_param
  two_body_mucm_param *params = static_cast<two_body_mucm_param*>( mucm_quad_param );
  params->func_count += 1;
  //  if( params->func_count % 100 == 0 )
  //  {
  //    Info( "angle_dist_F::mu_cm_quad_F", pastenum( "got ",
  //       params->func_count ) + " evaluations");
  //  }
  // get Eout_lab and mu_lab
  double Eout_lab;
  double mu_lab;
  if( params->use_relativistic )
  {
    params->relativistic_map->get_E_mu_lab( mucm, &Eout_lab, &mu_lab );
  }
  else
  {
    double E_in = params->get_E_in( );
    params->Newton_map->masses->two_body_get_E_mu_lab( E_in, mucm, &Eout_lab,
       &mu_lab );
  }
  // the Legendre polynomials
  math_F::Legendre( mu_lab, value );

  // the probability density
  double Prob = params->value( mucm );
  *value *= Prob;

  // do the energy weighting if necessary
  if( ( value->conserve == ENERGY ) || ( value->conserve == BOTH ) )
  {
    value->scale_E( Eout_lab );
  }
}
// ----------------  angle_dist_F::E_quad_F ------------------
void  angle_dist_F::E_quad_F( double E_in, QuadParamBase *e_quad_param,
  coef_vector *value )
// Function for the 2-d quadrature
{
  // the parameters are really two_body_Ein_param *
  two_body_Ein_param *e_params = static_cast<two_body_Ein_param *>( e_quad_param );
  e_params->func_count += 1;
  //  if( e_params->func_count % 100 == 0 )
  //  {
  //    Info( "E_quad_F", pastenum( "got ", e_params->func_count ) + " evaluations");
  //  }

  // The value of E_quad_F is itself an integral over mucm.
  // *value comes in as 0.  

  // parameters for the integration over mucm
  two_body_mucm_param mucm_params;
  mucm_params.set_E_in( E_in );
  mucm_params.Newton_map = &e_params->Newton_map;
  mucm_params.relativistic_map = &e_params->relativistic_map;
  mucm_params.use_relativistic = e_params->use_relativistic;
  // interpolate the (mucm, probability) with respect to incident energy
  mucm_params.first.linlin_interp( E_in, e_params->left_data_Ein,
    e_params->left_data, e_params->right_data_Ein, e_params->right_data );
  mucm_params.second.linlin_interp( E_in, e_params->left_data_Ein,
    e_params->next_left_data, e_params->right_data_Ein, e_params->next_right_data );

  // the range of integration
  double mucm_0;
  double mucm_1;
  if( e_params->use_relativistic )
  {
    mucm_params.relativistic_map->set_boost( E_in );

    mucm_0 = ( e_params->use_Eout_min ) ?
      mucm_params.relativistic_map->get_mu_cm( e_params->Eout_min ) :
        mucm_params.first.x;
    mucm_1 = ( e_params->use_Eout_max ) ?
      mucm_params.relativistic_map->get_mu_cm( e_params->Eout_max ) :
      mucm_params.second.x;
  }
  else
  {
    mucm_0 = ( e_params->use_Eout_min ) ?
      mucm_params.Newton_map->masses->get_mu_cm( E_in, e_params->Eout_min ) :
        mucm_params.first.x;
    mucm_1 = ( e_params->use_Eout_max ) ?
      mucm_params.Newton_map->masses->get_mu_cm( E_in, e_params->Eout_max ) :
        mucm_params.second.x;
  }

  // evaluate the integral over mucm
  QuadParamBase *params = static_cast< QuadParamBase* >( &mucm_params );
  static double tol = Global.Value( "quad_tol" );
  quad_F::integrate( angle_dist_F::mu_cm_quad_F, e_params->mu_quad_method,
    mucm_0, mucm_1, params, tol, value );

  e_params->mu_F_count += mucm_params.func_count;
  // weight it by flux * cross section
  e_params->set_weight( E_in );
  *value *= e_params->current_weight;
  //  cout << "E_in: " << E_in << " mucm_0: " << mucm_0 << " mucm_1: " <<
  //    mucm_1 << endl;
  //  value->print( );
}
