/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: energy_function.cpp 1 2006-02-02 03:06:56Z hedstrom $
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
//! Implementation of the classes used to handle formulas for energy probability density.
//! All data is in laboratory coordinates.

#include "energy_function.hpp"
#include "messaging.hpp"
#include "global_params.hpp"
#include "math_util.hpp"

// ************* class E_function_param *****************
// ---------------- E_function_param::set_Ein_default --------------------
// Interpolate the parameters
void E_function_param::set_Ein_default( double E_in )
{
  Theta = this_Theta->linlin_interp( E_in, *next_Theta );
  multiplicity = this_mult->linlin_interp( E_in, *next_mult );
  if( ( Theta <= 0.0 ) || ( multiplicity < 0.0 ) )
  {
    FatalError( "E_function_param::set_Ein_default", "got a negative parameter" );
  }
  // the range of integration
  E_max = E_in - U;
  if( E_max > top_E_out )
  {
    E_max = top_E_out;
  }
  E_1 = ( use_Eout_max ) ? Eout_max : E_max;
  set_scales( );
  norm = get_norm( );
  Eout_0 = Eout_min;
  Eout_1 = E_1;
}

// ************* class U_Ein_hit_list *****************
// ----------- U_Ein_hit_list::find_bottom_hits --------------
// Finds the intersection with the bottom of a box
void U_Ein_hit_list::find_bottom_hits( double E_out,
  vector< Ein_Eta_Hit > *Ein_hits )
{
  // for new entries
  Ein_Eta_Hit Ein_eta_hit;
  double Ein = E_out + get_U( );

  if( Ein > 0.0 )
  {
    // append this entry
    Ein_eta_hit.E_in = Ein;
    Ein_eta_hit.hit_edge = BOTTOM_IN;
    Ein_hits->push_back( Ein_eta_hit );
  }
}
// ----------- U_Ein_hit_list::find_top_hits --------------
// Finds the intersection with the top of a box
void U_Ein_hit_list::find_top_hits( double E_out,
  vector< Ein_Eta_Hit > *Ein_hits )
{
  // for new entries
  Ein_Eta_Hit Ein_eta_hit;
  double Ein = E_out + get_U( );

  if( Ein > 0.0 )
  {
    // append this entry
    Ein_eta_hit.E_in = Ein;
    Ein_eta_hit.hit_edge = TOP_OUT;
    Ein_hits->push_back( Ein_eta_hit );
  }
}


// ************* class energy_function *****************
// ---------------- energy_function::setup_data_default --------------------
// Initializes the quadrature parameters
void energy_function::setup_data_default( const Energy_groups& Eout_groups,
  E_function_param *Ein_param )
{
  Ein_param->U = U;
  Energy_groups::const_iterator Eout_ptr = Eout_groups.end();
  --Eout_ptr;
  Ein_param->top_E_out = *Eout_ptr;
  Ein_param->this_Theta = begin( );
  Ein_param->next_Theta = Ein_param->this_Theta;
  ++Ein_param->next_Theta;
  while( Ein_param->next_Theta->x <= *Ein_param->Ein_ptr )
  {
    Ein_param->this_Theta = Ein_param->next_Theta;
    ++Ein_param->next_Theta;
  }
  Ein_param->Theta_end = end( );

  Ein_param->upper_hits.set_U( Ein_param->U );
  double first_E = Ein_param->this_Theta->x;
  if( first_E > Ein_param->data_E_0 )
  {
    Ein_param->data_E_0 = first_E;
  }
}
// ---------------- energy_function::set_Ein_range_default --------------------
// Sets the range of incident energies for this intergration
void energy_function::set_Ein_range_default( int Ein_bin, E_function_param *Ein_param )
{
  Ein_param->set_Ein_range( );
  double this_E = Ein_param->this_Theta->x;
  if( this_E > Ein_param->data_E_0 ) Ein_param->data_E_0 = this_E;
  this_E = Ein_param->next_Theta->x;
  if( this_E < Ein_param->data_E_1 ) Ein_param->data_E_1 = this_E;

  if( Ein_param->data_E_1 < Ein_param->data_E_0 )
  {
    FatalError( "energy_function::set_Ein_range_default", "check the Theta incident energies" );
  }
}
// ---------------- energy_function::next_ladder_default --------------------
// Default go to the next (incident energy, Theta).  Returns "true" when finished.
bool energy_function::next_ladder_default( double E_in, E_function_param *Ein_param )
{
  bool done = Ein_param->update_bin_pointers( E_in );
  static double etol = Global.Value( "E_tol" );
  double E_tol = E_in * etol;
  //    double E_tol = 0.0;
  if( !done )
  {
    if( E_in + E_tol >= Ein_param->next_Theta->x )
    {
      while( E_in + E_tol >= Ein_param->next_Theta->x )
      {
        // get the next (E_in, Theta) data
        Ein_param->this_Theta = Ein_param->next_Theta;
        ++Ein_param->next_Theta;
        if( Ein_param->next_Theta == end( ) )
        {
          return true;
        }
      }
    }
  }
  return done;
}
// ----------- energy_function::Eout_ladder --------------
// This routine uses the angular distributions this_eta_dist and the
// next to calculate the contribution to the E_out boxes of the
// transfer matrix between incident energies Ein_param->data_E_0 and
// Ein_param->data_E_1.
void energy_function::Eout_ladder( T_matrix& transfer, E_function_param *Ein_param )
{
  bool geom_OK;  // for checking the consistency of the geometry
  U_Ein_hit_list test_hits;
  // loop through the outgoing energies (column of transfer)
  for( int Eout_count = 0; Eout_count < transfer.num_Eout_bins;
    ++Eout_count )
  {
    vector< double >::const_iterator Eout_ptr = transfer.out_groups.begin( )
      + Eout_count;
    // how does the line Eout = Ein - U meet this E-E' box?
    double U = Ein_param->U;
    geom_OK = Ein_param->upper_hits.hit_box( U, Eout_ptr,
      Ein_param->data_E_0, Ein_param->data_E_1 );
    if( !geom_OK )
    {
      test_hits.set_U( U );
      test_hits.hit_box( U, Eout_ptr,
                         Ein_param->data_E_0, Ein_param->data_E_1 );
      test_hits.print( );
      FatalError( "energy_function::Eout_ladder", "Check the coding" );
    }
    if( ( Eout_count > 0 ) && ( Ein_param->upper_hits.is_below( ) ) )
    {
      // we are done with this incident energy bin
      break;
    }
    // integrate over this E-E' box
    one_Ebox( transfer, Eout_count, Ein_param );
  }
}
// ----------- energy_function::one_Ebox --------------
// Does the integration for one E-E' box
void energy_function::one_Ebox( T_matrix& transfer, int Eout_count,
   E_function_param *Ein_param )
{
  static double from_quad_tol = Global.Value( "abs_quad_tol" )/100;
  static double from_abs_tol = 1000*Global.Value( "abs_tol" );
  double skip_tol = ( from_quad_tol > from_abs_tol ) ? from_quad_tol : from_abs_tol;
  // the E' energy range
  Ein_param->Eout_min = transfer.out_groups[ Eout_count ];
  Ein_param->Eout_max = transfer.out_groups[ Eout_count + 1 ];

  // integrate depending on how the line Eout = Ein - U meets the box
  U_Ein_hit_list::iterator high_hit_ptr = Ein_param->upper_hits.begin( );
  U_Ein_hit_list::iterator next_high_ptr = high_hit_ptr;
  ++next_high_ptr;
  for( ; next_high_ptr != Ein_param->upper_hits.end( );
         high_hit_ptr = next_high_ptr, ++next_high_ptr )
  {
    // always integrate from Eout_min
    Ein_param->use_Eout_min = true;
    // where is the line Eout = Ein - U?
    if( high_hit_ptr->hit_edge == BELOW )
    {
      // do nothing---we are below the E-E' box
      continue;
    }
    else if( ( high_hit_ptr->hit_edge == BOTTOM_IN ) ||
             ( high_hit_ptr->hit_edge == INSIDE ) )
    {
      // the line Eout = Ein - U is inside the E-E' box
      Ein_param->use_Eout_max = false;
    }
    else
    {
      // the line Eout = Ein - U is above the E-E' box
      Ein_param->use_Eout_max = true;
    }
    // the range of integration in incident energy
    Ein_param->Ein_0 = high_hit_ptr->E_in;
    Ein_param->Ein_1 = next_high_ptr->E_in;
    if( Ein_param->Ein_1 - Ein_param->Ein_0 <= Ein_param->Ein_1 * skip_tol )
    {
      Warning( "energy_function::one_Ebox", "skipping a very short interval");
      continue;  // skip this interval
    }
    if( transfer.interpolate_Eout_integrals )
    {
      interp_update_T( transfer, Eout_count, Ein_param );
    }
    else
    {
      update_T( transfer, Eout_count, Ein_param );
    }
  }
}
// ----------- energy_function::update_T --------------
// Adds to an element of transfer the integral between over the E-E' box
void energy_function::update_T( T_matrix &transfer, int Eout_count,
   E_function_param *Ein_param )
{
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
    double tol = Ein_param->set_tol( left_E, right_E );

    quad_F::integrate( Energy_function_F::Ein_F, transfer.Ein_quad_method, left_E,
		       right_E, params, tol, &value );

    if( value.weight_1[ 0 ] < 0.0 )
    {
      Warning( "energy_function::update_T", pastenum( "negative integral ", left_E ) +
	       pastenum(" ", right_E ) );
      value.set_zero( );  // throw out these values
    }

    // add this integral
    transfer( Ein_param->Ein_count, Eout_count ) += value;
    // increment the function counts
    Ein_param->Ein_F_count += Ein_param->func_count;
    ++Ein_param->quad_count;
  }
}
// ----------- energy_function::interp_update_T --------------
// Adds to an element of transfer the interpolated integral over the E-E' box
void energy_function::interp_update_T( T_matrix &transfer, int Eout_count,
  E_function_param *Ein_param )
{
  static double etol = Global.Value( "E_tol" );
  // A list of integrals over E_out for interpolation over one energy bin
  Eout_integrals Eout_ints;
  Eout_ints.Eout_int_params = Ein_param;

  // initial integration over E_out
  Eout_ints.setup_Eout_ints( transfer.order, transfer.conserve,
    Ein_param->Ein_0, Ein_param->Ein_1 );
  Ein_param->this_Eout_int = Eout_ints.begin( );
  Ein_param->next_Eout_int = Ein_param->this_Eout_int;
  ++Ein_param->next_Eout_int;
  Ein_param->Eout_int_end = Eout_ints.end( );

  // a vector to store the integrals
  coef_vector value( transfer.order, transfer.conserve );

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
         ( Ein_param->this_sigma->x < Ein_param->Ein_1 ); )
  {
    double left_E = ( Ein_param->this_sigma->x < Ein_param->Ein_0 ) ?
      Ein_param->Ein_0 : Ein_param->this_sigma->x;
    double right_E = ( Ein_param->next_sigma->x > Ein_param->Ein_1 ) ?
      Ein_param->Ein_1 : Ein_param->next_sigma->x;
    if( left_E < Ein_param->this_Eout_int->E_in )
    {
      left_E = Ein_param->this_Eout_int->E_in;
    }
    if( right_E > Ein_param->next_Eout_int->E_in )
    {
      right_E = Ein_param->next_Eout_int->E_in;
    }
    // evaluate the integral exactly by 4th-order Gaussian quadrature
    double tol = 0.0;  // not used here
    quad_F::integrate( Energy_function_F::interp_Ein_F, GAUSS4, left_E,
		       right_E, params, tol, &value );

    // add this integral
    transfer( Ein_param->Ein_count, Eout_count ) += value;
    // increment the function counts
    Ein_param->Ein_F_count += Ein_param->func_count;
    ++Ein_param->quad_count;

    // the next subinterval
    double E_tol = right_E * etol;
    if( right_E > Ein_param->next_Eout_int->E_in - E_tol )
    {
      Ein_param->this_Eout_int = Ein_param->next_Eout_int;
      ++Ein_param->next_Eout_int;
      if( Ein_param->next_Eout_int == Ein_param->Eout_int_end )
      {
	break;
      }
    }
    if( right_E > Ein_param->next_sigma->x - E_tol )
    {
      Ein_param->this_sigma = Ein_param->next_sigma;
      ++Ein_param->next_sigma;
    }
  }
}

// ************** Probability density model ******************************
// ----------- Energy_function_F::Ein_F --------------
//! Integral function for the model
void Energy_function_F::Ein_F( double E_in, QuadParamBase *Ein_param,
   coef_vector *value )
{
  // the parameters are really E_function_param *
  E_function_param *e_params = static_cast<E_function_param *>( Ein_param );
  e_params->func_count += 1;
  e_params->set_Ein( E_in );  // interpolate the data

  e_params->get_integrals( e_params->Eout_min, e_params->E_1, *value );

  // weight it by flux * cross section * multiplicity * model_weight
  e_params->set_weight( E_in );
  *value *= e_params->current_weight;
}

// ----------- Energy_function_F::interp_Ein_F --------------
//! Interpolated integral function for the model
void Energy_function_F::interp_Ein_F( double E_in, QuadParamBase *Ein_param,
  coef_vector *value )
{
  // the parameters are really E_function_param *
  E_function_param *e_params = static_cast<E_function_param *>( Ein_param );
  e_params->func_count += 1;
  e_params->set_Ein( E_in );  // interpolate the data
  e_params->this_Eout_int->Interpolate( E_in,
    *(e_params->next_Eout_int), value );

  // weight it by flux * cross section * multiplicity * model_weight
  e_params->set_weight( E_in );
  *value *= e_params->current_weight;
}
