/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: param_base.cpp 1 2006-02-02 03:06:56Z hedstrom $
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
// implementation of the classes used to handle uncorrelated energy-angle distributions

#include <cmath>

#include "param_base.hpp"
#include "messaging.hpp"
#include "global_params.hpp"

using namespace std;

// ************* class param_base *****************
// ----------- param_base::setup --------------
// Sets up the quadrature parameters
void param_base::setup( const dd_vector& sigma_, const dd_vector& mult_,
  const dd_vector& weight_,
  const Flux_List& e_flux_, const Energy_groups& Ein_groups )
{
  if( weight_.interp_type != HISTOGRAM )
  {
    FatalError( "param_base::setup",
		 "only histogram interpolation has been implemented for the model weight" );
  }
  // pointers to the cross section
  first_ladder_sigma = sigma_.begin( );
  sigma_end = sigma_.end( );
  if( first_ladder_sigma->y == 0.0 )
  {
    next_sigma = first_ladder_sigma;
    for( ++next_sigma; next_sigma != sigma_end;
	 first_ladder_sigma = next_sigma, ++next_sigma )
    {
      if( next_sigma->y != 0.0 )
      {
	break;
      }
    }
  }

  // pointers to the multiplicity
  this_mult = mult_.begin( );
  next_mult = this_mult;
  ++next_mult;
  mult_end = mult_.end( );

  // pointers to the model weight
  this_weight = weight_.begin( );
  next_weight = this_weight;
  ++next_weight;
  weight_end = weight_.end( );

  // pointers to the flux data
  flux_ptr = e_flux_.begin( );
  next_flux = flux_ptr;
  ++next_flux;
  flux_end = e_flux_.end( );
  order = flux_ptr->order;
  current_weight.initialize( order );

  // pointers to the incident energy boundaries
  Ein_ptr = Ein_groups.begin( );
  next_Ein = Ein_ptr;
  ++next_Ein;
  Ein_end = Ein_groups.end( );
  Ein_count = 0;

  // Set the first common incident energy
  common_E0( );
}
// ----------- param_base::setup_bin --------------
// Sets up the quadrature parameters for integration on a bin
void param_base::setup_bin( int Ein_bin, const dd_vector& sigma_,
  const dd_vector& mult_, const dd_vector& weight_,
  const Flux_List& e_flux_, const Energy_groups& Ein_groups )
{
  // the incident energy bin
  Ein_count = Ein_bin;
  Ein_ptr = Ein_groups.begin( ) + Ein_count;
  next_Ein = Ein_ptr;
  ++next_Ein;

  // pointers to the cross section
  first_ladder_sigma = sigma_.begin( );
  next_sigma = first_ladder_sigma;
  ++next_sigma;
  sigma_end = sigma_.end( );
  while( ( next_sigma->x <= *Ein_ptr ) && ( next_sigma != sigma_end ) )
  {
    first_ladder_sigma = next_sigma;
    ++next_sigma;
  }

  // pointers to the multiplicity
  this_mult = mult_.begin( );
  next_mult = this_mult;
  ++next_mult;
  mult_end = mult_.end( );
  while( ( next_mult->x <= *Ein_ptr ) && ( next_mult != mult_end ) )
  {
    this_mult = next_mult;
    ++next_mult;
  }

  // pointers to the model weight
  this_weight = weight_.begin( );
  next_weight = this_weight;
  ++next_weight;
  weight_end = weight_.end( );
  while( ( next_weight->x <= *Ein_ptr ) && ( next_weight != weight_end ) )
  {
    this_weight = next_weight;
    ++next_weight;
  }

  // pointers to the flux data
  flux_ptr = e_flux_.begin( );
  next_flux = flux_ptr;
  ++next_flux;
  flux_end = e_flux_.end( );
  while( ( next_flux->get_E_in( ) <= *Ein_ptr ) && ( next_flux != flux_end ) )
  {
    flux_ptr = next_flux;
    ++next_flux;
  }

  order = flux_ptr->order;
  current_weight.initialize( order );

  // Set the first common incident energy
  common_E0( );
}
// ----------- param_base::common_E0 --------------
// Sets the first common incident energy
void param_base::common_E0( )
{
  data_E_0 = ( first_ladder_sigma->x < *Ein_ptr ) ? *Ein_ptr :
    first_ladder_sigma->x;
  if( this_mult->x > data_E_0 ) data_E_0 = this_mult->x;
  if( this_weight->x > data_E_0 ) data_E_0 = this_weight->x;
  if( flux_ptr->get_E_in( ) > data_E_0 ) data_E_0 = flux_ptr->get_E_in( );
}
// ----------- param_base::set_Ein_range --------------
// Sets the range of interpolagtion over incident energy
void param_base::set_Ein_range( )
{
  data_E_0 = ( this_mult->x < *Ein_ptr ) ? *Ein_ptr :
    this_mult->x;
  if( this_weight->x > data_E_0 ) data_E_0 = this_weight->x;
  if( flux_ptr->get_E_in( ) > data_E_0 ) data_E_0 = flux_ptr->get_E_in( );

  // Ensure that energies are in order
  bool data_bad = update_pointers( data_E_0 );
  if( data_bad )
  {
    FatalError( "param_base::set_Ein_range", "incident energies out of order" );
  }

  data_E_1 = ( next_mult->x > *next_Ein ) ? *next_Ein :
    next_mult->x;
  if( next_weight->x < data_E_1 ) data_E_1 = next_weight->x;
  if( next_flux->get_E_in( ) < data_E_1 ) data_E_1 = next_flux->get_E_in( );
}
// ----------- param_base::set_sigma_range --------------
// Sets the range of pointers to cross sections for this set of data
void param_base::set_sigma_range( )
{
  dd_vector::const_iterator sigma_ptr = first_ladder_sigma;
  while( sigma_ptr->x <= data_E_0 )
  {
    ++sigma_ptr;
    if( sigma_ptr == sigma_end )
    {
      FatalError( "param_base::set_sigma_range", "we ran out of cross section data" );
    }
  }
  if( first_ladder_sigma != sigma_ptr )
  {
    first_ladder_sigma = sigma_ptr;
    --first_ladder_sigma;
  }

  while( sigma_ptr->x < data_E_1 )
  {
    ++sigma_ptr;
    if( sigma_ptr == sigma_end )
    {
      break;
    }
  }
  last_ladder_sigma = sigma_ptr;
}
// ----------- param_base::get_Ein_range --------------
// Gets the range of nontrivial incident energy bins
// returns true if the threshold is too high for the energy bins
bool param_base::get_Ein_range( const dd_vector& sigma, const dd_vector& mult,
  const dd_vector& weight,
  const Flux_List& e_flux, const Energy_groups& Ein_groups,
  double *first_Ein, double *last_Ein )
{
  bool done = false;
  // The cross section sometimes starts with several zeros.
  dd_vector::const_iterator this_sigma = sigma.begin( );
  dd_vector::const_iterator next_sigma = this_sigma;
  ++next_sigma;
  if( this_sigma->y == 0.0 )
  {
    for( ; next_sigma != sigma_end;
	 this_sigma = next_sigma, ++next_sigma )
    {
      if( next_sigma->y != 0.0 )
      {
	break;
      }
    }
  }
  double E_sigma = this_sigma->x;
  double E_mult = mult.begin( )->x;
  double E_weight = weight.begin( )->x;
  double E_flux = e_flux.begin( )->get_E_in( );

  double E_first = ( E_sigma < E_mult ) ? E_mult : E_sigma;
  if( E_weight > E_first ) E_first = E_weight;
  if( E_flux > E_first ) E_first = E_flux;
  Energy_groups::const_iterator last_bin = Ein_groups.end( );
  --last_bin;
  if( E_first >= *last_bin )
  {
    *first_Ein = 0.0;
    *last_Ein = 0.0;
    return true;
  }
  *first_Ein = E_first;

  dd_vector::const_iterator sigma_ptr = sigma.end( );
  --sigma_ptr;
  E_sigma = sigma_ptr->x;
  dd_vector::const_iterator mult_ptr = mult.end( );
  --mult_ptr;
  E_mult = mult_ptr->x;
  dd_vector::const_iterator weight_ptr = weight.end( );
  --weight_ptr;
  E_weight = weight_ptr->x;
  Flux_List::const_iterator Flux_ptr = e_flux.end( );
  --Flux_ptr;
  E_flux = Flux_ptr->get_E_in( );

  double E_last = ( E_sigma > E_mult ) ? E_mult : E_sigma;
  if( E_weight < E_last ) E_last = E_weight;
  if( E_flux < E_last ) E_last = E_flux;
  *last_Ein = E_last;
  return done;
}
// ----------- param_base::update_pointers --------------
// Sets the data pointers for a new incident energy interval.
// Returns "true" at the end of the data
bool param_base::update_pointers( double E_in )
{
  static double etol = Global.Value( "E_tol" );
  double E_tol = E_in * etol;
  //  double E_tol = 0.0;
    while( E_in + E_tol >= next_mult->x )
    {
      // increment multiplicity
      this_mult = next_mult;
      ++next_mult;
      if( next_mult == mult_end )
      {
        return true;
      }
    }
    while( E_in + E_tol >= next_weight->x )
    {
      // increment model weight
      this_weight = next_weight;
      ++next_weight;
      if( next_weight == weight_end )
      {
        return true;
      }
    }
    while( E_in + E_tol >= *next_Ein )
    {
      // go to the next incident energy bin
      ++Ein_count;
      Ein_ptr = next_Ein;
      ++next_Ein;
      if( next_Ein == Ein_end )
      {
        return true;
      }
    }
    while( E_in + E_tol >= next_flux->get_E_in( ) )
    {
      // increment flux
      flux_ptr = next_flux;
      ++next_flux;
      if( next_flux == flux_end )
      {
        return true;
      }
    }

  return false;
}
// ----------- param_base::update_bin_pointers --------------
// Sets the data pointers for a new incident energy interval in one bin.
// Returns "true" at the end of the data on this bin.
bool param_base::update_bin_pointers( double E_in )
{
  static double etol = Global.Value( "E_tol" );
  double E_tol = E_in * etol;
  //  double E_tol = 0.0;
    if( E_in + E_tol >= *next_Ein )
    {
      return true;
    }
    while( E_in + E_tol >= next_mult->x )
    {
      // increment multiplicity
      this_mult = next_mult;
      ++next_mult;
      if( next_mult == mult_end )
      {
        return true;
      }
    }
    while( E_in + E_tol >= next_weight->x )
    {
      // increment model weight
      this_weight = next_weight;
      ++next_weight;
      if( next_weight == weight_end )
      {
        return true;
      }
    }
    while( E_in + E_tol >= next_flux->get_E_in( ) )
    {
      // increment flux
      flux_ptr = next_flux;
      ++next_flux;
      if( next_flux == flux_end )
      {
        return true;
      }
    }

  return false;
}
// ----------- param_base::set_weight --------------
// Calculates the weight (cross section) * flux * multiplicity * model weight
void param_base::set_weight( double E_in )
{
  current_weight.linlin_interp( E_in, *flux_ptr, *next_flux );
  double scale_by = this_sigma->linlin_interp( E_in, *next_sigma ) *
     this_mult->linlin_interp( E_in, *next_mult ) * this_weight->y;
  current_weight *= scale_by;
}
// ----------- param_base::flux_weight --------------
// Calculates the weight for gammas: flux * multiplicity * model weight
void param_base::flux_weight( double E_in )
{
  current_weight.linlin_interp( E_in, *flux_ptr, *next_flux );
  double scale_by = this_mult->linlin_interp( E_in, *next_mult ) * this_weight->y;
  current_weight *= scale_by;
}
