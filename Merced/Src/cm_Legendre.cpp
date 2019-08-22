/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 381 $
 * $Date: 2013-06-14 (Fri, Jun 14, 2013) $
 * $Author: hedstrom $
 * $Id: cm_Legendre.cpp 1 2013-06-14 hedstrom $
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
// implementation of the classes used to handle Legendre expansions of energy 
// probability density given in the center-of-mass frame

#include <cmath>
#include <cstdlib>
#include <cfloat>
#ifdef _OPENMP
 #include <omp.h>
#endif

#include "cm_Legendre.hpp"
#include "messaging.hpp"
#include "global_params.hpp"
#include "math_util.hpp"

// ************* class cm_Legendre_Ein_param *****************
// ----------- cm_Legendre_Ein_param::set_data --------------
// Sets up the unit-base Legendre data for interpolation with respect to incident energy
void cm_Legendre_Ein_param::set_data( )
{
  cm_Legendre_vector::const_iterator next_left_data = left_Eout_data;
  ++next_left_data;
  cm_Legendre_vector::const_iterator next_right_data = right_Eout_data;
  ++next_right_data;

  Ein0_data.set_data( *left_Eout_data, *next_left_data, prev_data_Eout,
     next_data_Eout );
  Ein1_data.set_data( *right_Eout_data, *next_right_data, prev_data_Eout,
     next_data_Eout );
}
// ----------- cm_Legendre_Ein_param::ubase_interpolate --------------
// Do unit-base interpolation and return the physical Legendre coefficients
void cm_Legendre_Ein_param::ubase_interpolate( double Ein )
{
  mid_data.ubase_interpolate( Ein, Ein0_data, Ein1_data );
  // map to physical variables
  mid_data.un_unit_base( );
  Ecm_params.this_data = &mid_data;
}
// ----------- cm_Legendre_Ein_param::set_Ecm_param --------------
// Sets up the parameters for integration over center-of-mass outgoing energy
void cm_Legendre_Ein_param::set_Ecm_param( )
{
  Ecm_params.map = map;
  Ecm_params.mu_quad_method = mu_quad_method;  // method for integration over cosine
  Ecm_params.Eout_interp = Eout_interp;
  Ecm_params.setup( mid_data.get_E_in( ), Eout_min, Eout_max,
     mid_data.prev_data.get_E_out( ), mid_data.next_data.get_E_out( ) );
  Ecm_params.V_lab_sectors( );
}
// ----------- cm_Legendre_Ein_param::start_Eout_cm --------------
// Starts one staircase of the Eout_cm data
void cm_Legendre_Ein_param::start_Eout_cm( )
{
  Ein0_data.set_E_in( left_Ein_data->get_E_in( ) );
  Ein1_data.set_E_in( right_Ein_data->get_E_in( ) );

  left_Eout_data = left_Ein_data->begin( );
  next_left_Eout_data = left_Eout_data;
  ++next_left_Eout_data;
  right_Eout_data = right_Ein_data->begin( );
  next_right_Eout_data = right_Eout_data;
  ++next_right_Eout_data;

  if( Ein_interp.qualifier == UNITBASE )
  {
    setup_Ein_ubase( );
  }
  else if( Ein_interp.qualifier == CUMULATIVE_POINTS )
  {
    setup_Ein_cum_prob( );
  }
  else
  {
    FatalError( "cm_Legendre_Ein_param::start_Eout_cm",
		"Incident energy interpolation not implemented" );
  }

}
// ----------- cm_Legendre_Ein_param::setup_Ein_ubase --------------
// Sets up the data for unit-base interpolation in incident energy
void cm_Legendre_Ein_param::setup_Ein_ubase( )
{
  // The following coding is safe, because the unit-base outgoing energies are 0 <= E <= 1.
  prev_data_Eout = 0.0;
  double left_Eout = next_left_Eout_data->get_E_out( );
  double right_Eout = next_right_Eout_data->get_E_out( );
  static double etol = Global.Value( "E_tol" );
  if( left_Eout < right_Eout*(1 + etol ) )
  {
    next_data_Eout = left_Eout;
  }
  else
  {
    next_data_Eout = right_Eout;
  }

  // save the unit-base interpolation
  Ein0_data.ubase_map.copy( left_Ein_data->ubase_map );
  Ein1_data.ubase_map.copy( right_Ein_data->ubase_map );

  set_data( );
}
// ---------------- cm_Legendre_Ein_param::setup_Ein_cum_prob ------------------
// Sets up the data for cumulative points interpolation in incident energy
void cm_Legendre_Ein_param::setup_Ein_cum_prob( )
{
  // for cumulative points interpolation these are fixed unit-base values
  prev_data_Eout = 0.0;
  next_data_Eout = 1.0;

  // lower incident energy
  left_cum_prob = left_Ein_data->cum_prob.begin( );
  next_left_cum_prob = left_cum_prob;
  ++next_left_cum_prob;

  // skip zero probability intervals
  while( ( left_cum_prob->Prob == 0.0 ) &&
         ( left_cum_prob->slope == 0.0 ) )
  {
    left_cum_prob = next_left_cum_prob;
    ++next_left_cum_prob;
    left_Eout_data = next_left_Eout_data;
    ++next_left_Eout_data;
  }

  // higher incident energy
  right_cum_prob = right_Ein_data->cum_prob.begin( );
  next_right_cum_prob = right_cum_prob;
  ++next_right_cum_prob;

  // skip zero probability intervals
  while( ( right_cum_prob->Prob == 0.0 ) &&
         ( right_cum_prob->slope == 0.0 ) )
  {
    right_cum_prob = next_right_cum_prob;
    ++next_right_cum_prob;
    right_Eout_data = next_right_Eout_data;
    ++next_right_Eout_data;
  }

  // for the range of cumulative probabilities A
  //  double lower_A = 0.0;
  double higher_A = ( next_left_cum_prob->cum_prob < next_right_cum_prob->cum_prob )?
      next_left_cum_prob->cum_prob : next_right_cum_prob->cum_prob;

  // set up Ein0_data and Ein1_data
  setup_low_A( );
  setup_high_A( higher_A );
  Ein0_data.to_unit_base( );
  Ein1_data.to_unit_base( );
}
// ---------------- cm_Legendre_Ein_param::setup_low_A ------------------
// Sets (Eout, probability) data to the lower zero cumulative probability
void cm_Legendre_Ein_param::setup_low_A( )
{
  Ein0_data.prev_data.set_E_out( left_Eout_data->get_E_out( ) );
  Ein1_data.prev_data.set_E_out( right_Eout_data->get_E_out( ) );

  Ein0_data.prev_data.copy_coef( *left_Eout_data );
  Ein1_data.prev_data.copy_coef( *right_Eout_data );
}
// ---------------- cm_Legendre_Ein_param::setup_high_A ------------------
// Interpolates (Eout, probability) data to the higher common cumulative probability
void cm_Legendre_Ein_param::setup_high_A( double higher_A )
{
  double higher_Eout;

  if( next_left_cum_prob->cum_prob == higher_A )
  {
    Ein0_data.next_data.set_E_out( next_left_Eout_data->get_E_out( ) );
    Ein0_data.next_data.copy_coef( *next_left_Eout_data );
  }
  else
  {
    higher_Eout = left_cum_prob->get_cum_inv( higher_A );
    Ein0_data.next_data.set_E_out( higher_Eout );
    if( Ein0_data.Eout_interp == HISTOGRAM )
    {
      Ein0_data.next_data.copy_coef( *left_Eout_data );
    }
    else
    {
      // check space allocation
      Ein0_data.next_data.set_max_order( left_Eout_data->order, 
					 next_left_Eout_data->order );
      Ein0_data.next_data.linlin_interp( higher_Eout, *left_Eout_data, *next_left_Eout_data );
    }
  }

  if( next_right_cum_prob->cum_prob == higher_A )
  {
    Ein1_data.next_data.set_E_out( next_right_Eout_data->get_E_out( ) );
    Ein1_data.next_data.copy_coef( *next_right_Eout_data );
  }
  else
  {
    higher_Eout = right_cum_prob->get_cum_inv( higher_A );
    Ein1_data.next_data.set_E_out( higher_Eout );
    if( Ein1_data.Eout_interp == HISTOGRAM )
    {
      Ein1_data.next_data.copy_coef( *right_Eout_data );
    }
    else
    {
      // check space allocation
      Ein1_data.next_data.set_max_order( right_Eout_data->order, 
					 next_right_Eout_data->order );
      Ein1_data.next_data.linlin_interp( higher_Eout, *right_Eout_data,
                     *next_right_Eout_data );
    }
  }
}
// ----------- cm_Legendre_Ein_param::next_Ecm --------------
// Gets the next interval of unit-base outgoing energy
bool cm_Legendre_Ein_param::next_Ecm( )
{
  // which outgoing intervals do we increment?
  prev_data_Eout = next_data_Eout;
  double left_Eout = next_left_Eout_data->get_E_out( );
  static double etol = Global.Value( "E_tol" );
  if( left_Eout < prev_data_Eout*(1 + etol ) )
  {
    left_Eout_data = next_left_Eout_data;
    ++next_left_Eout_data;
    if( next_left_Eout_data == left_Ein_data->end( ) )
    {
      return( true );
    }
  }
  double right_Eout = next_right_Eout_data->get_E_out( );
  if( right_Eout < prev_data_Eout*(1 + etol ) )
  {
    right_Eout_data = next_right_Eout_data;
    ++next_right_Eout_data;
    if( next_right_Eout_data == right_Ein_data->end( ) )
    {
      return( true );
    }
  }

  // find the upper common outgoing energy
  left_Eout = next_left_Eout_data->get_E_out( );
  right_Eout = next_right_Eout_data->get_E_out( );
  if( left_Eout < right_Eout*(1 + etol ) )
  {
    next_data_Eout = left_Eout;
  }
  else
  {
    next_data_Eout = right_Eout;
  }

  set_data( );

  return( false );
}
// ----------- cm_Legendre_Ein_param::CP_next_Ecm --------------
// Gets the next interval of outgoing energy for cumulative points interpolation
bool cm_Legendre_Ein_param::CP_next_Ecm( )
{
  bool done = false;

  // ignore intervals with probability less than skip_tol
  static double skip_tol = Global.Value( "abs_tol" );

  // undo the unit-base map before copying data
  Ein0_data.un_unit_base( );
  Ein1_data.un_unit_base( );

  Ein0_data.prev_data.set_E_out( Ein0_data.next_data.get_E_out( ) );
  Ein0_data.prev_data.copy_coef( Ein0_data.next_data );
  Ein1_data.prev_data.set_E_out( Ein1_data.next_data.get_E_out( ) );
  Ein1_data.prev_data.copy_coef( Ein1_data.next_data );

  // update the pointers
  if( next_left_Eout_data->get_E_out( ) <= Ein0_data.prev_data.get_E_out( ) )
  {
    left_Eout_data = next_left_Eout_data;
    ++next_left_Eout_data;
    if( next_left_Eout_data == left_Ein_data->end( ) )
    {
      return true;
    }

    // Update the left cumulative probabilities
    left_cum_prob = next_left_cum_prob;
    ++next_left_cum_prob;

    bool do_skip = false;
    // skip intervals with essentially zero probability
    while( next_left_cum_prob->cum_prob - left_cum_prob->cum_prob <= skip_tol )
    {
      do_skip = true;
      left_cum_prob = next_left_cum_prob;
      ++next_left_cum_prob;
      left_Eout_data = next_left_Eout_data;
      ++next_left_Eout_data;
      if( next_left_Eout_data == left_Ein_data->end( ) )
      {
        return true;
      }
    }
    if( do_skip )
    {
      Ein0_data.prev_data.set_E_out( left_Eout_data->get_E_out( ) );
      Ein0_data.prev_data.copy_coef( *left_Eout_data );
    }
  }
  if( next_right_Eout_data->get_E_out( ) <= Ein1_data.prev_data.get_E_out( ) )
  {
    right_Eout_data = next_right_Eout_data;
    ++next_right_Eout_data;
    if( next_right_Eout_data == right_Ein_data->end( ) )
    {
      return true;
    }

    // Update the right cumulative probabilities
    right_cum_prob = next_right_cum_prob;
    ++next_right_cum_prob;

    bool do_skip = false;
    // skip intervals with essentially zero probability
    while( next_right_cum_prob->cum_prob - right_cum_prob->cum_prob <= skip_tol )
    {
      do_skip = true;
      right_cum_prob = next_right_cum_prob;
      ++next_right_cum_prob;
      right_Eout_data = next_right_Eout_data;
      ++next_right_Eout_data;
      if( next_right_Eout_data == right_Ein_data->end( ) )
      {
        return true;
      }
    }
    if( do_skip )
    {
      Ein1_data.prev_data.set_E_out( right_Eout_data->get_E_out( ) );
      Ein1_data.prev_data.copy_coef( *right_Eout_data );
    }
  }

  // Interpolate to the common higher cumulative probability
  double higher_A = ( next_left_cum_prob->cum_prob < next_right_cum_prob->cum_prob )?
      next_left_cum_prob->cum_prob : next_right_cum_prob->cum_prob;
  setup_high_A( higher_A );
  Ein0_data.to_unit_base( );
  Ein1_data.to_unit_base( );

  // We may have skipped an interval with zero probability
  return done;
}

// ************* class cm_Legendre_mu_param **********************
// ----------- cm_Legendre_mu_param::setup --------------
// Sets up the data for this incident energy
void cm_Legendre_mu_param::setup( double Ein, double Eoutcm,
   const cm_Legendre_Ecm_param& Ecm_param )
{
  Ecm_Elab_mu_param::setup( Ein, Eoutcm, Ecm_param );
  this_data.set_E_out( Eoutcm );

  if( ( Ecm_param.Eout_interp == HISTOGRAM ) ||
      ( Eoutcm <= Ecm_param.this_data->prev_data.get_E_out( ) ) )
  {
    this_data.set_max_order( Ecm_param.this_data->prev_data.order,
			   Ecm_param.this_data->prev_data.order );
    this_data.copy_coef( Ecm_param.this_data->prev_data );
  }
  else
  {
    this_data.set_max_order( Ecm_param.this_data->prev_data.order,
			   Ecm_param.this_data->next_data.order );
    this_data.linlin_interp( Eoutcm,
      Ecm_param.this_data->prev_data, Ecm_param.this_data->next_data );
  }
}

// ************* class cm_Legendre_vector *****************
// ----------- cm_Legendre_vector::read_coef --------------
// Reads the Legendre coefficients of the energy probability density
void cm_Legendre_vector::read_coef( data_parser& infile, int num_Eout )
{
  cm_Legendre_vector::iterator new_Eout_ptr;
  // read the data
  for( int Eout_count = 0; Eout_count < num_Eout; ++Eout_count )
  {
    // make a new set of Legendre coefficients
    new_Eout_ptr = insert( end( ), Legendre_coefs( ) );
    new_Eout_ptr->set_E_out( infile.get_next_double( ) );  // energy of outgoing particle
    int file_order = infile.get_next_int( ) - 1;  // Legendre order of input data
    new_Eout_ptr->initialize( file_order );
    for( int coef_count = 0; coef_count <= file_order; ++coef_count )
    {
      ( *new_Eout_ptr )[ coef_count ] = infile.get_next_double( );
    }
    new_Eout_ptr->truncate_zeros( );
  }
  // ensure that the norm is 1
  renorm( );

  if( Ein_interp.qualifier == UNITBASE )
  {
    to_unit_base( );
  }
  else if( Ein_interp.qualifier == CUMULATIVE_POINTS )
  {
    form_cum_prob( );
  }
}
// ----------- cm_Legendre_vector::read_order_zero --------------
// Reads the isotropic energy probability density
void cm_Legendre_vector::read_order_zero( data_parser& infile, int num_Eout )
{
  cm_Legendre_vector::iterator new_Eout_ptr;
  // read the data
  for( int Eout_count = 0; Eout_count < num_Eout; ++Eout_count )
  {
    // make a new set of Legendre coefficients
    new_Eout_ptr = insert( end( ), Legendre_coefs( ) );
    new_Eout_ptr->set_E_out( infile.get_next_double( ) );  // energy of outgoing particle
    new_Eout_ptr->initialize( 0 );  // only order 0
    ( *new_Eout_ptr )[ 0 ] = infile.get_next_double( );
  }

  // ensure that the norm is 1
  renorm( );

  if( Ein_interp.qualifier == UNITBASE )
  {
    to_unit_base( );
  }
  else if( Ein_interp.qualifier == CUMULATIVE_POINTS )
  {
    form_cum_prob( );
  }
}
// ----------- cm_Legendre_vector::form_cum_prob --------------
// Forms the list of cumulative probabilities
void cm_Legendre_vector::form_cum_prob( )
{
  // copy the data
  cum_prob.Eout_interp = Eout_interp;
  for( cm_Legendre_vector::const_iterator Eout_ptr = begin( );
       Eout_ptr != end( ); ++Eout_ptr )
  {
    cumulative_prob_list::iterator cum_prob_ptr = cum_prob.insert(
      cum_prob.end( ), cumulative_prob_entry( ) );
    cum_prob_ptr->E_out = Eout_ptr->get_E_out( );
    cum_prob_ptr->Prob = Eout_ptr->value( 0 );
  }
  // now form the slopes and cumulative probabilities
  if( Eout_interp == HISTOGRAM )
  {
    cum_prob.get_cum_prob_flat( );
  }
  else // lin-lin
  {
    cum_prob.get_cum_prob_linlin( );
  }
}

// ************* class cm_Legendre *****************
// ----------- cm_Legendre::read_data ------------------
// Reads the Legendre data
void cm_Legendre::read_data( data_parser& infile, int num_Ein )
{
  cm_Legendre::iterator new_Ein_ptr;
  for( int Ein_count = 0; Ein_count < num_Ein; ++Ein_count )
  {
    // make a new cm_Legendre_vector
    new_Ein_ptr = insert( end( ), cm_Legendre_vector( ) );
    new_Ein_ptr->set_E_in( infile.get_next_double( ) );  // energy of incident particle
    new_Ein_ptr->Ein_interp = Ein_interp;
    new_Ein_ptr->Eout_interp = Eout_interp;
    int num_Eout = infile.get_next_int( );  // how many outgoing energies
    new_Ein_ptr->read_coef( infile, num_Eout );
  }
}
// ----------- cm_Legendre::read_isotropic ------------------
// Reads isotropic data
void cm_Legendre::read_isotropic( data_parser& infile, int num_Ein )
{
  cm_Legendre::iterator new_Ein_ptr;
  for( int Ein_count = 0; Ein_count < num_Ein; ++Ein_count )
  {
    // make a new cm_Legendre_vector
    new_Ein_ptr = insert( end( ), cm_Legendre_vector( ) );
    new_Ein_ptr->set_E_in( infile.get_next_double( ) );  // energy of incident particle
    new_Ein_ptr->Ein_interp = Ein_interp;
    new_Ein_ptr->Eout_interp = Eout_interp;
    int num_Eout = infile.get_next_int( );  // how many outgoing energies
    new_Ein_ptr->read_order_zero( infile, num_Eout );
  }
}
// ----------- cm_Legendre::setup_map ------------------
// Sets up the map from center-of-mass to laboratory coordinates
void cm_Legendre::setup_map( )
{
  particles.check_data( );
  map.setup_ratios( particles.mProj, particles.mTarg, particles.mProd );
}
// ----------- cm_Legendre::get_Ein_range --------------
//  Gets the range of nontrivial incident energy bins; computes first_Ein and last_Ein
// returns true if the threshold is too high for the energy bins
bool cm_Legendre::get_Ein_range( const dd_vector& sigma, const dd_vector& mult,
    const dd_vector& weight,
    const Flux_List& e_flux, const Energy_groups& Ein_groups )
{
  double E_first;
  double E_last;
  cm_Legendre_Ein_param initial_param;
  bool done = initial_param.get_Ein_range( sigma, mult, weight, e_flux,
                                         Ein_groups, &E_first, &E_last );
  if( done ) return true;

  // check the range of incident energies for the probability data
  cm_Legendre::const_iterator Ein_data_ptr = begin( );
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
// ----------- cm_Legendre::get_T ------------------
// Computes the transfer matrix
void cm_Legendre::get_T( const dd_vector& sigma, const dd_vector& multiple, 
  const dd_vector& weight, T_matrix& transfer )
{
  bool interp_OK = ( ( Ein_interp.qualifier == UNITBASE ) &&
                     ( Ein_interp.flag == LINLIN ) ) ||
    ( ( Ein_interp.qualifier == CUMULATIVE_POINTS ) &&
      ( ( Ein_interp.flag == LINLIN ) ||
        ( Ein_interp.flag == HISTOGRAM ) ) );
  if( !interp_OK )
  {
    FatalError( "cm_Legendre::get_T", "Incident interpolation type not implemented" );
  }

  interp_OK = ( Eout_interp == LINLIN ) || ( Eout_interp == HISTOGRAM );
  if( !interp_OK )
  {
    FatalError( "cm_Legendre::get_T", "Outgoing interpolation type not implemented" );
  }
  if( particles.mProd == 0.0 )
  {
    FatalError( "cm_Legendre::get_T", "gamma emission not implemented" );
  }

  // print( );
  setup_map( );
  bool done = get_Ein_range( sigma, multiple, weight, transfer.e_flux,
    transfer.in_groups );
  if( done )
  {
    transfer.zero_transfer( );
  }

  long int quad_count = 0;  // number of 2-d quadratures
  long int Ein_F_count= 0;  // number of calls to cm_Legendre_F::Ein_F
  long int Ecm_F_count= 0;  // number of calls to cm_Legendre_F::Ecm_F
  long int mu_F_count = 0;  // number of calls to cm_Legendre_F::mu_cm_quad_F

  // now do the integrals bin by bin
#pragma omp parallel for schedule( dynamic, 1 ) default( none )	\
  shared( sigma, multiple, weight, transfer ) \
  reduction( +: quad_count ) reduction( +: Ein_F_count ) \
  reduction( +: Ecm_F_count ) reduction( +: mu_F_count )
  for( int Ein_bin = first_Ein; Ein_bin < last_Ein; ++Ein_bin )
  {
    cm_Legendre_Ein_param Ein_param;
    // set up the data range for this bin
    Ein_param.setup_bin( Ein_bin, sigma, multiple, weight, transfer.e_flux,
                         transfer.in_groups );
    Ein_param.Eout_quad_method = transfer.Eout_quad_method;
    Ein_param.mu_quad_method = transfer.mu_quad_method;
    setup_param( &Ein_param );
    for( ; ; )
    {
      set_Ein_range( &Ein_param );   // get the incident energy interval
      Eout_data_ladder( transfer, &Ein_param );  // loop over the outgoing energies
      bool Done = next_ladder( Ein_param.data_E_1, &Ein_param );   // go to the next interval
      if( Done )
      {
        break;
      }
    }
    quad_count += Ein_param.quad_count;
    Ein_F_count += Ein_param.Ein_F_count;
    Ecm_F_count += Ein_param.Ecm_F_count;
    mu_F_count += Ein_param.Ecm_params.mu_F_count;
  } // end of parallel loop

  // print the counts of function evaluations
  cout << "3d quadratures: " << quad_count << endl;
  cout << "cm_Legendre_F::Ein_F calls: " << Ein_F_count << endl;
  cout << "cm_Legendre_F::Ecm_F calls: " << Ecm_F_count << endl;
  cout << "cm_Legendre_F::mu_F calls: " << mu_F_count << endl;
  cout << "average cm_Legendre_F::Ein_F_count: " << 1.0*Ein_F_count/quad_count << endl;
  cout << "average cm_Legendre_F::Ecm_F_count: " << 1.0*Ecm_F_count/Ein_F_count << endl;
  cout << "average cm_Legendre_F::mu_F_count: " << 1.0*mu_F_count/Ecm_F_count << endl;
}
// ----------- cm_Legendre::setup_param ------------------
// Initializes the quadrature parameters; returns true if the threshold is too high
void cm_Legendre::setup_param( cm_Legendre_Ein_param *Ein_param )
{
  Ein_param->Ein_interp = Ein_interp;
  Ein_param->Eout_interp = Eout_interp;
  Ein_param->map = &map;

  Ein_param->Ein0_data.Ein_interp = Ein_interp;
  Ein_param->Ein0_data.Eout_interp = Eout_interp;
  Ein_param->Ein1_data.Ein_interp = Ein_interp;
  Ein_param->Ein1_data.Eout_interp = Eout_interp;

  Ein_param->left_Ein_data = begin( );
  Ein_param->right_Ein_data = Ein_param->left_Ein_data;
  ++Ein_param->right_Ein_data;
  while( Ein_param->right_Ein_data->get_E_in( ) <= Ein_param->data_E_0 )
  {
    Ein_param->left_Ein_data = Ein_param->right_Ein_data;
    ++Ein_param->right_Ein_data;
    if( Ein_param->right_Ein_data == end( ) )
    {
      FatalError( "cm_Legendre::setup_param",
             "incident energy ranges inconsistent" );
    }
  }
  double first_Ein = Ein_param->left_Ein_data->get_E_in( );
  if( first_Ein > Ein_param->data_E_0 )
  {
    Ein_param->data_E_0 = first_Ein;
  }
  // the Vcm_Vlab_hit_list objects need the gamma for the energy of translation of the center of mass
  Ein_param->lower_hits.G0_data.gamma = map.gamma;
  Ein_param->lower_hits.G1_data.gamma = map.gamma;
  Ein_param->upper_hits.G0_data.gamma = map.gamma;
  Ein_param->upper_hits.G1_data.gamma = map.gamma;
}
// -----------  cm_Legendre::set_Ein_range ------------------
// Sets the range of incident energies for this intergration
void cm_Legendre::set_Ein_range( cm_Legendre_Ein_param *Ein_param )
{
  Ein_param->set_Ein_range( );
  double this_E = Ein_param->left_Ein_data->get_E_in( );
  if( this_E > Ein_param->data_E_0 ) Ein_param->data_E_0 = this_E;
  this_E = Ein_param->right_Ein_data->get_E_in( );
  if( this_E < Ein_param->data_E_1 ) Ein_param->data_E_1 = this_E;

  if( Ein_param->data_E_1 < Ein_param->data_E_0 )
  {
    FatalError( "cm_Legendre::set_Ein_range", "check the incident energies" );
  }
}
// -----------  cm_Legendre::next_ladder ------------------
// Go to the next pair of incident energies.  Returns "true" when finished.
bool cm_Legendre::next_ladder( double E_in, cm_Legendre_Ein_param *Ein_param )
{
  bool done = Ein_param->update_bin_pointers( E_in );
  static double etol = Global.Value( "E_tol" );
  double E_tol = E_in * etol;
  if( !done )
  {
    if( E_in + E_tol >= Ein_param->left_Ein_data->get_E_in( ) )
    {
      while( E_in + E_tol >= Ein_param->right_Ein_data->get_E_in( ) )
      {
        // get the next E_in cm_Legendre data
        Ein_param->left_Ein_data = Ein_param->right_Ein_data;
        ++Ein_param->right_Ein_data;
        if( Ein_param->right_Ein_data == end( ) )
        {
          return true;
        }
      }
    }
  }
  return done;
}
// ----------- cm_Legendre::Eout_data_ladder --------------
// Loops through the center-of-mass data for a pair of incident energies
void cm_Legendre::Eout_data_ladder( T_matrix& transfer,
   cm_Legendre_Ein_param *Ein_param )
{
  Ein_param->start_Eout_cm( );
  bool done = false;
  double left_Ein = Ein_param->left_Ein_data->get_E_in( );
  double right_Ein = Ein_param->right_Ein_data->get_E_in( );
  // loop through the center-of-mass outgoing data
  while( !done )
  {
    double phys_Ecm = Ein_param->Ein0_data.ubase_map.un_unit_base( Ein_param->prev_data_Eout );
    Ein_param->lower_hits.G0_data.set_energies( left_Ein, phys_Ecm );

    phys_Ecm = Ein_param->Ein1_data.ubase_map.un_unit_base( Ein_param->prev_data_Eout );
    Ein_param->lower_hits.G1_data.set_energies( right_Ein, phys_Ecm );


    phys_Ecm = Ein_param->Ein0_data.ubase_map.un_unit_base( Ein_param->next_data_Eout );
    Ein_param->upper_hits.G0_data.set_energies( left_Ein, phys_Ecm );

    phys_Ecm = Ein_param->Ein1_data.ubase_map.un_unit_base( Ein_param->next_data_Eout );
    Ein_param->upper_hits.G1_data.set_energies( right_Ein, phys_Ecm );

    lab_Eout_ladder( transfer, Ein_param );
    if( Ein_interp.qualifier == CUMULATIVE_POINTS )
    {
      done = Ein_param->CP_next_Ecm( );
    }
    else
    {
      done = Ein_param->next_Ecm( );
    }
  }
}
// ----------- cm_Legendre::lab_Eout_ladder --------------
// Loops through the outgoing energy bins given cm_Eout data
void cm_Legendre::lab_Eout_ladder( T_matrix& transfer, cm_Legendre_Ein_param *Ein_param)
{
  //  bool check_geometry = true;
  bool check_geometry = false;
  bool geom_OK;  // for checking the consistency of the geometry
  bool upper_hits_set = false;
  bool lower_hits_set = false;
  Vcm_Vlab_hit_list test_hits;
  test_hits.G0_data.gamma = map.gamma;
  test_hits.G1_data.gamma = map.gamma;
  double left_Ein = Ein_param->left_Ein_data->get_E_in( );
  double right_Ein = Ein_param->right_Ein_data->get_E_in( );
  double dummy = 0.0;
  double phys_Ecm;
  int Eout_count = 0;
  vector< double >::const_iterator Eout_ptr = transfer.out_groups.begin( );
  vector< double >::const_iterator next_Eout = Eout_ptr;
  ++next_Eout;
  // Check for only forward emission
  if( Ein_param->upper_hits.G1_data.E_cm < Ein_param->upper_hits.G1_data.get_Etrans( ) )
  {
    geom_OK = Ein_param->upper_hits.hit_box( dummy, Eout_ptr, Ein_param->data_E_0, Ein_param->data_E_1 );
    upper_hits_set = true;
    if( check_geometry )
    {
      cout << "Forward with next_Eout: " << *next_Eout << endl;
      Ein_param->upper_hits.print( );
    }
    if( !geom_OK )
    {
      phys_Ecm = Ein_param->left_Ein_data->ubase_map.un_unit_base( Ein_param->next_data_Eout );
      test_hits.G0_data.set_energies( left_Ein, phys_Ecm );
      phys_Ecm = Ein_param->right_Ein_data->ubase_map.un_unit_base( Ein_param->next_data_Eout );
      test_hits.G1_data.set_energies( right_Ein, phys_Ecm );
      geom_OK = test_hits.hit_box( dummy, Eout_ptr, Ein_param->data_E_0, Ein_param->data_E_1 );
      test_hits.print( );
      FatalError( "cm_Legendre::lab_Eout_ladder", "Check the coding, 1" );
    }
    while( Ein_param->upper_hits.is_above( ) )
    {
      Eout_ptr = next_Eout;
      ++next_Eout;
      if( next_Eout == transfer.out_groups.end( ) )
      {
        return;
      }
      ++Eout_count;
      geom_OK = Ein_param->upper_hits.hit_box( dummy, Eout_ptr, Ein_param->data_E_0, Ein_param->data_E_1 );
      if( check_geometry )
      {
        cout << "next Forward with next_Eout: " << *next_Eout << endl;
        Ein_param->upper_hits.print( );
      }
      if( !geom_OK )
      {
        phys_Ecm = Ein_param->left_Ein_data->ubase_map.un_unit_base( Ein_param->next_data_Eout );
        test_hits.G0_data.set_energies( left_Ein, phys_Ecm );
        phys_Ecm = Ein_param->right_Ein_data->ubase_map.un_unit_base( Ein_param->next_data_Eout );
        test_hits.G1_data.set_energies( right_Ein, phys_Ecm );
        geom_OK = test_hits.hit_box( dummy, Eout_ptr, Ein_param->data_E_0, Ein_param->data_E_1 );
        test_hits.print( );
        FatalError( "cm_Legendre::lab_Eout_ladder", "Check the coding, 2" );
      }
    }
  }
  else if( Ein_param->lower_hits.G1_data.E_cm > Ein_param->lower_hits.G1_data.get_Etrans( ) )
  {
    // Check whether all emission is above the lab energy bin
    geom_OK = Ein_param->lower_hits.hit_box( dummy, Eout_ptr, Ein_param->data_E_0, Ein_param->data_E_1 );
    lower_hits_set = true;
    if( check_geometry )
    {
      cout << "Backward with E_out: " << *Eout_ptr << endl;
      Ein_param->lower_hits.print( );
    }
    if( !geom_OK )
    {
      phys_Ecm = Ein_param->left_Ein_data->ubase_map.un_unit_base( Ein_param->prev_data_Eout );
      test_hits.G0_data.set_energies( left_Ein, phys_Ecm );
      phys_Ecm = Ein_param->right_Ein_data->ubase_map.un_unit_base( Ein_param->prev_data_Eout );
      test_hits.G1_data.set_energies( right_Ein, phys_Ecm );
      geom_OK = test_hits.hit_box( dummy, Eout_ptr, Ein_param->data_E_0, Ein_param->data_E_1 );
      test_hits.print( );
      FatalError( "cm_Legendre::lab_Eout_ladder", "Check the coding, 3" );
    }
    while( Ein_param->lower_hits.is_below( ) )
    {
      Eout_ptr = next_Eout;
      ++next_Eout;
      if( next_Eout == transfer.out_groups.end( ) )
      {
        return;
      }
      ++Eout_count;
      geom_OK = Ein_param->lower_hits.hit_box( dummy, Eout_ptr, Ein_param->data_E_0, Ein_param->data_E_1 );
      if( check_geometry )
      {
        cout << "backward with E_out: " << *Eout_ptr << endl;
        Ein_param->lower_hits.print( );
      }
      if( !geom_OK )
      {
        phys_Ecm = Ein_param->left_Ein_data->ubase_map.un_unit_base( Ein_param->prev_data_Eout );
        test_hits.G0_data.set_energies( left_Ein, phys_Ecm );
        phys_Ecm = Ein_param->right_Ein_data->ubase_map.un_unit_base( Ein_param->prev_data_Eout );
        test_hits.G1_data.set_energies( right_Ein, phys_Ecm );
        geom_OK = test_hits.hit_box( dummy, Eout_ptr, Ein_param->data_E_0, Ein_param->data_E_1 );
        test_hits.print( );
        FatalError( "cm_Legendre::lab_Eout_ladder", "Check the coding, 4" );
      }
    }
  }
  // Now, compute integrals until the lab energy bin is above the E_cm data
  for( ; Eout_count < transfer.num_Eout_bins;
       ++Eout_count, Eout_ptr = next_Eout, ++next_Eout )
  {
    if( !upper_hits_set )
    {
      geom_OK = Ein_param->upper_hits.hit_box( dummy, Eout_ptr, Ein_param->data_E_0, Ein_param->data_E_1 );
      if( check_geometry )
      {
        cout << "upper_hits for Eout: " << *Eout_ptr << endl;
        Ein_param->upper_hits.print( );
      }
      if( !geom_OK )
      {
        phys_Ecm = Ein_param->left_Ein_data->ubase_map.un_unit_base( Ein_param->next_data_Eout );
        test_hits.G0_data.set_energies( left_Ein, phys_Ecm );
        phys_Ecm = Ein_param->right_Ein_data->ubase_map.un_unit_base( Ein_param->next_data_Eout );
        test_hits.G1_data.set_energies( right_Ein, phys_Ecm );
        geom_OK = test_hits.hit_box( dummy, Eout_ptr, Ein_param->data_E_0, Ein_param->data_E_1 );
        test_hits.print( );
        FatalError( "cm_Legendre::lab_Eout_ladder", "Check the coding, 5" );
      }
    }
    if( Ein_param->upper_hits.is_below( ) )
    {
      break;  // we are done
    }
    if( !lower_hits_set )
    {
      geom_OK = Ein_param->lower_hits.hit_box( dummy, Eout_ptr, Ein_param->data_E_0, Ein_param->data_E_1 );
      if( check_geometry )
      {
        cout << "lower_hits for Eout: " << *Eout_ptr << endl;
        Ein_param->lower_hits.print( );
      }
      if( !geom_OK )
      {
        phys_Ecm = Ein_param->left_Ein_data->ubase_map.un_unit_base( Ein_param->prev_data_Eout );
        test_hits.G0_data.set_energies( left_Ein, phys_Ecm );
        phys_Ecm = Ein_param->right_Ein_data->ubase_map.un_unit_base( Ein_param->prev_data_Eout );
        test_hits.G1_data.set_energies( right_Ein, phys_Ecm );
        geom_OK = test_hits.hit_box( dummy, Eout_ptr, Ein_param->data_E_0, Ein_param->data_E_1 );
        test_hits.print( );
        FatalError( "cm_Legendre::lab_Eout_ladder", "Check the coding, 6" );
      }
    }
    // integrate over this E-E' box
    one_Ebox( transfer, Eout_count, Ein_param );
    upper_hits_set = false;
    lower_hits_set = false;
  }
}
// -----------  cm_Legendre::one_Ebox ------------------
// Does the integration for one Eout_lab annulus between a pair of incident energies
void cm_Legendre::one_Ebox( T_matrix& transfer, int Eout_count,
  cm_Legendre_Ein_param *Ein_param )
{
  // the E' energy range
  Ein_param->Eout_min = transfer.out_groups[ Eout_count ];
  Ein_param->Eout_max = transfer.out_groups[ Eout_count + 1 ];
  //  cout << Ein_param->Eout_min << " < E_out < " << Ein_param->Eout_max << endl;
  // set up common incident energies
  Ein_param->lower_hits.common_hits( Ein_param->upper_hits );

  // integrate depending on how the arcs E_cm = const meet the box
  Vcm_Vlab_hit_list::iterator low_hit_ptr = Ein_param->lower_hits.begin( );
  Vcm_Vlab_hit_list::iterator next_low_ptr = low_hit_ptr;
  ++next_low_ptr;
  Vcm_Vlab_hit_list::iterator high_hit_ptr = Ein_param->upper_hits.begin( );
  Vcm_Vlab_hit_list::iterator next_high_ptr = high_hit_ptr;
  ++next_high_ptr;
  for( ; ( next_low_ptr != Ein_param->lower_hits.end( ) ) &&
         ( next_high_ptr != Ein_param->upper_hits.end( ) );
       low_hit_ptr = next_low_ptr, ++next_low_ptr,
         high_hit_ptr = next_high_ptr, ++next_high_ptr )
  {
    if( ( ( low_hit_ptr->hit_edge == ABOVE ) &&
          ( high_hit_ptr->hit_edge == ABOVE ) ) ||
        ( ( next_low_ptr->hit_edge == ABOVE ) &&
          ( next_high_ptr->hit_edge == ABOVE ) ) ||
        ( ( low_hit_ptr->hit_edge == ABOVE_FORWARD ) &&
          ( high_hit_ptr->hit_edge == ABOVE_FORWARD ) ) ||
        ( ( next_low_ptr->hit_edge == ABOVE_FORWARD ) &&
          ( next_high_ptr->hit_edge == ABOVE_FORWARD ) ) ||
        ( ( low_hit_ptr->hit_edge == BELOW ) &&
          ( high_hit_ptr->hit_edge == BELOW ) ) ||
        ( ( next_low_ptr->hit_edge == BELOW ) &&
          ( next_high_ptr->hit_edge == BELOW ) ) )
    {
      continue;
    }
    // the range of integration in incident energy
    Ein_param->Ein_0 = low_hit_ptr->E_in;
    Ein_param->Ein_1 = next_low_ptr->E_in;
    update_T( transfer, Eout_count, Ein_param );
  }
}
// -----------  cm_Legendre::update_T ------------------
// Adds to an element of transfer the integral between the intersections of 2 Eout_cm = const arcs with the Eout_lab box
void cm_Legendre::update_T( T_matrix &transfer, int Eout_count,
   cm_Legendre_Ein_param *Ein_param )
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
    double left_E = ( Ein_param->this_sigma->x < Ein_param->Ein_0 ) ? Ein_param->Ein_0 :
      Ein_param->this_sigma->x;
    double right_E = ( Ein_param->next_sigma->x > Ein_param->Ein_1 ) ? Ein_param->Ein_1 :
      Ein_param->next_sigma->x;
    // evaluate the integral
    quad_F::integrate( cm_Legendre_F::Ein_F, transfer.Ein_quad_method, left_E, right_E,
		       params, tol, &value );

    // add this integral
    transfer( Ein_param->Ein_count, Eout_count ) += value;
    // increment the function counts
    Ein_param->Ein_F_count += Ein_param->func_count;
    Ein_param->quad_count += Ein_param->Vcm_hit_count;
  }
}

// **************** Functions to integrate *********************
// --------------------  cm_Legendre_F::mu_F ------------------
// Function for the 1-d quadrature over cm cosine
void cm_Legendre_F::mu_F( double mu, QuadParamBase *mu_quad_param, coef_vector *value )
{
  // the parameters are really cm_Legendre_mu_param
  cm_Legendre_mu_param *mu_params = static_cast< cm_Legendre_mu_param* >( mu_quad_param );
  mu_params->func_count += 1;
  //   if( mu_params->func_count % 100 == 0 )
  //   {
  //     Info( "cm_Legendre_F::mu_F", pastenum( "got ", mu_params->func_count ) + " evaluations");
  //   }
  static bool negativeMu = false;
  double Eout_lab;
  double mu_lab;
  mu_params->map->get_E_mu_lab( mu_params->E_in, mu_params->Eout_cm, mu, &Eout_lab,
    &mu_lab );

  // the Legendre polynomials
  math_F::Legendre( mu_lab, value );
  // the energy-angle probability density
  double Prob = mu_params->this_data.sum_Legendre( mu );
  if(( Prob < 0.0 ) && !negativeMu )
  {
    Warning( "cm_Legendre_F::mu_F", pastenum( "probability negative for mu: ", mu ) +
	pastenum( "  Ecm: ", mu_params->Eout_cm) + pastenum( "  Ein: ", mu_params->E_in) );
    negativeMu = true;
  }
  *value *= Prob;

  // do the energy weighting if necessary
  if( ( value->conserve == ENERGY ) || ( value->conserve == BOTH ) )
  {
    value->scale_E( Eout_lab );
  }
}
// ------------------- cm_Legendre_F::Ecm_F ------------------
// Function for the 2-d quadrature over cm cosine and Eout_cm
void cm_Legendre_F::Ecm_F( double Eout_cm, QuadParamBase *Ecm_quad_param, coef_vector *value )
{
  // the parameters are really cm_Legendre_Ecm_param *
  cm_Legendre_Ecm_param *Ecm_param = static_cast<cm_Legendre_Ecm_param *>( Ecm_quad_param );
  Ecm_param->func_count += 1;

  //   if( Ecm_param->func_count % 100 == 0 )
  //   {
  //     Info( "cm_Legendre_F::Ecm_F", pastenum( "got ", Ecm_param->func_count ) + " evaluations");
  //   }

  // The value of cm_Legendre_F::Ecm_F is itself an integral over cm cosine.
  // *value comes in as 0.  

  // parameters for the integration over cm cosine
  cm_Legendre_mu_param mu_param;
  mu_param.setup( Ecm_param->E_in, Eout_cm, *Ecm_param );

  // evaluate the integral over eta
  QuadParamBase *params = static_cast< QuadParamBase* >( &mu_param );
  static double tol = Global.Value( "quad_tol" );
  quad_F::integrate( cm_Legendre_F::mu_F, Ecm_param->mu_quad_method, mu_param.mu_cm_min,
		     mu_param.mu_cm_max, params, tol, value );

  Ecm_param->mu_F_count += mu_param.func_count;
}
// ------------------- cm_Legendre_F::Ein_F ------------------
// Function for the 3-d quadrature over E_in, and Eout_cm and cm cosine
// The value of cm_Legendre_F::Ein_F is itself an integral over Eout_cm and cm cosine.
void cm_Legendre_F::Ein_F( double E_in, QuadParamBase *Ein_quad_param, coef_vector *value )
{
  value->set_zero( );  // initialize to zero
  // the parameters are really cm_Legendre_Ein_param *
  cm_Legendre_Ein_param *Ein_param = static_cast<cm_Legendre_Ein_param *>( Ein_quad_param );

  //   if( Ein_param->func_count % 100 == 0 )
  //   {
  //     Info( "cm_Legendre_F::Ein_F", pastenum( "got ", Ein_param->func_count ) + " evaluations");
  //   }
  Ein_param->Vcm_hit_count = 0;   // number of local calls to quad_F::integrate
  Ein_param->ubase_interpolate( E_in );

  // set up the parameters for the integration over Eout_cm and cm cosine
  Ein_param->set_Ecm_param( );
  QuadParamBase *params = static_cast< QuadParamBase* >( &Ein_param->Ecm_params );

  // Integrate over sectors of ( Eout_lab, mu_cm ) space
  coef_vector one_value( value->order, value->conserve );

  list< Vcm_quadBox_Hit >::const_iterator this_V_hit =
    Ein_param->Ecm_params.V_cm_limits.begin( );
  list< Vcm_quadBox_Hit >::const_iterator next_V_hit = this_V_hit;
  ++next_V_hit;
  for( ; next_V_hit != Ein_param->Ecm_params.V_cm_limits.end( );
         this_V_hit = next_V_hit, ++next_V_hit )
  {
    if( ( next_V_hit->V_cm <= Ein_param->Ecm_params.min_V_cm ) ||
        ( this_V_hit->hit_corner == V_BELOW ) )
    {
      continue;  // current V_cm values are below
    }
    else if( ( this_V_hit->V_cm >= Ein_param->Ecm_params.max_V_cm ) ||
             ( next_V_hit->hit_corner == V_ABOVE ) )
    {
      break;  // all remaining V_cm values are above
    }
    else
    {
      Ein_param->Vcm_hit_min = *this_V_hit;
      Ein_param->Vcm_hit_max = *next_V_hit;
    }
    //    cout << "integrate V_lab" << endl;
    //    Ein_param->Vcm_hit_min.print( );
    //    Ein_param->Vcm_hit_max.print( );
    Ein_param->Ecm_params.min_hit_corner = Ein_param->Vcm_hit_min.hit_corner;
    Ein_param->Ecm_params.max_hit_corner = Ein_param->Vcm_hit_max.hit_corner;
    double tol = Ein_param->Ecm_params.Ecm_range( );

    quad_F::integrate( cm_Legendre_F::Ecm_F, Ein_param->Eout_quad_method,
		       Ein_param->Ecm_params.Ecm_min,
		       Ein_param->Ecm_params.Ecm_max, params, tol, &one_value );
    *value += one_value;
    // we actually want to count the number of 3d integrals
    Ein_param->Vcm_hit_count += 1;
    Ein_param->func_count += 1;
    Ein_param->Ecm_F_count += Ein_param->Ecm_params.func_count;
  }
  // weight it by flux * cross section
  Ein_param->set_weight( E_in );
  *value *= Ein_param->current_weight;
  //  cout << "E_in: " << E_in << " eta_0: " << eta_0 << " eta_1: " <<
  //    eta_1 << endl;
  //  value->print( );
}
