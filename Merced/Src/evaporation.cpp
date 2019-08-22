/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: evaporation.cpp 1 2006-02-02 03:06:56Z hedstrom $
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
//! Classes used for energy probability density given as a function

#include <cmath>
#ifdef _OPENMP
 #include <omp.h>
#endif

#include "evaporation.hpp"
#include "global_params.hpp"
#include "messaging.hpp"
#include "math_util.hpp"

#include "protos.h"

// ****************** class Evap_param **********************
// ----------- Evap_param::set_tol --------------
// Sets the tolerance for the quadrature over incident energy
double Evap_param::set_tol( double left_E, double right_E )
{
  static double tol = Global.Value( "quad_tol" );  // the default tolerance
  return tol;
}
// ---------------- Evap_param::get_integrals --------------------
// Gets the integrals over outgoing energy
void Evap_param::get_integrals( double Eout_0, double Eout_1, coef_vector &value )
{
  if( Eout_0 == 0.0 )
  {
    if( ( value.conserve == NUMBER ) || ( value.conserve == BOTH ) )
    {
      value.weight_1[ 0 ] = integral_prob( Eout_1 );
    }
    if( ( value.conserve == ENERGY ) || ( value.conserve == BOTH ) )
    {
      value.weight_E[ 0 ] = integral_Eprob( Eout_1 );
    }
  }
  else
  {
    static double tol = Global.Value( "quad_tol" );
    evap_params_1d params_1d;
    params_1d.Theta = Theta;
    QuadParamBase *void_params_1d = static_cast< QuadParamBase* >( &params_1d );
    quad_F::integrate( evaporation_F::Eout_F, Eout_quad_method, Eout_0, Eout_1,
		       void_params_1d, tol, &value );
    Eout_F_count += params_1d.func_count;
  }
  value *= 1.0/norm;
}
// ---------------- Evap_param::tol_get_integrals --------------------
// Gets the integrals over outgoing energy and returns the noise in the calculation
double Evap_param::tol_get_integrals( double Eout_0, double Eout_1, coef_vector &value )
{
  // the default is no noise---Watt and Madland-Nix are special
  get_integrals( Eout_0, Eout_1, value );
  return set_tol( Ein_0, Ein_1 );
}
// ---------------- Evap_param::integral_prob --------------------
// Integral of the probability density from 0 to E_1
double Evap_param::integral_prob( double E_1 )
{
  double y_1 = E_1/Theta;
  return prob_integral_scale*igam( 2.0, y_1 );
}
// ---------------- Evap_param::integral_Eprob --------------------
// Integral of energy times the probability density from 0 to E_1
double Evap_param::integral_Eprob( double E_1 )
{
  double y_1 = E_1/Theta;
  return Eprob_integral_scale*igam( 3.0, y_1 );
}
// ---------------- Evap_param::set_scales --------------------
//! Sets the scale factors for the integrals of probability and energy*probability
void Evap_param::set_scales( )
{
  // the scale factors used with the incomplete gamma functions
  prob_integral_scale = Theta*Theta;  // includes Gamma( 2 ) = 1
  Eprob_integral_scale = 2.0*prob_integral_scale*Theta;  // includes Gamma( 3 ) = 2
}
// ---------------- Evap_param::get_norm --------------------
// Integrate from 0 to E_max to get the norm
double Evap_param::get_norm( )
{
  return integral_prob( E_max );
}
  
// ****************** class evaporation **********************
// ----------- evaporation::get_Ein_range --------------
//  Gets the range of nontrivial incident energy bins; computes first_Ein and last_Ein
// returns true if the threshold is too high for the energy bins
bool evaporation::get_Ein_range( const dd_vector& sigma, const dd_vector& mult,
    const dd_vector& weight,
    const Flux_List& e_flux, const Energy_groups& Ein_groups )
{
  double E_first;
  double E_last;
  Evap_param initial_param;
  bool done = initial_param.get_Ein_range( sigma, mult, weight, e_flux,
                                         Ein_groups, &E_first, &E_last );
  if( done ) return true;

  // check the range of incident energies for the probability data
  evaporation::const_iterator Theta_ptr = begin( );
  double E_data = Theta_ptr->x;
  if( E_data > E_first )
  {
    E_first = E_data;
  }
  first_Ein = Ein_groups.first_bin_ID( E_first );

  Theta_ptr = end( );
  --Theta_ptr;
  E_data = Theta_ptr->x;
  if( E_data < E_last )
  {
    E_last = E_data;
  }
  last_Ein = Ein_groups.last_bin_ID( E_last );

  return false;
}
// ----------- evaporation::get_T --------------
// Calculates the transfer matrix for this particle.
// sigma is the cross section.
void evaporation::get_T( const dd_vector& sigma, const dd_vector& multiple,
  const dd_vector& weight, T_matrix& transfer )
{
  if( interp_type != LINLIN )
  {
    FatalError( "evaporation::get_T", "interp_type for Theta not implemented" );
  }
  transfer.getBinCrossSection( sigma );

  // synchronize the starting energies
  bool done = get_Ein_range( sigma, multiple, weight, transfer.e_flux,
    transfer.in_groups );
  if( done )
  {
    transfer.zero_transfer( );
  }

  long int quad_count = 0;  // number of quadratures
  long int Ein_F_count= 0;  // number of calls to Energy_function_F::Ein_F
  long int Eout_F_count= 0;  // number of calls to evaporation_F::Eout_F

  // do the integrals incident bin by incident bin
#pragma omp parallel for schedule( dynamic, 1 ) default( none )	\
  shared( sigma, multiple, weight, transfer ) \
  reduction( +: quad_count ) reduction( +: Ein_F_count ) reduction( +: Eout_F_count )
  // now do the integrals incident bin by incident bin
  for( int Ein_bin = first_Ein; Ein_bin < last_Ein; ++Ein_bin )
  {
    Evap_param Ein_param;
    Ein_param.Eout_quad_method = transfer.Eout_quad_method;
    // set up the data range for this bin
    Ein_param.setup_bin( Ein_bin, sigma, multiple, weight, transfer.e_flux,
                         transfer.in_groups );
    setup_data( transfer.out_groups, &Ein_param );
    // work on this bin
    for( ; ; )
    {
      set_Ein_range( Ein_bin, &Ein_param );   // get the incident energy interval
      Ein_param.set_sigma_range( );
      Eout_ladder( transfer, &Ein_param );
      bool done = next_ladder( Ein_param.data_E_1, &Ein_param );   // go to the next interval
      if( done )
      {
        break;
      }
    }
    quad_count += Ein_param.quad_count;
    Ein_F_count += Ein_param.Ein_F_count;
    Eout_F_count += Ein_param.Eout_F_count;
  } // end of parallel loop

  // print the counts of function evaluations
  cout << "2d quadratures: " << quad_count << endl;
  cout << "Energy_function_F::Ein_F calls: " << Ein_F_count << endl;
  cout << "evaporation_F::Eout_F calls: " << Eout_F_count << endl;
  cout << "average Energy_function_F::Ein_F calls: " << 1.0*Ein_F_count/quad_count << endl;
  cout << "average evaporation_F::Eout_F calls: " << 1.0*Eout_F_count/Ein_F_count << endl;
}

// **************** Function to integrate *********************
// --------------------  evaporation_F::Eout_F ------------------
// Function for the 1-d quadrature over outgoing energy
void evaporation_F::Eout_F( double E_out, QuadParamBase *E_out_param,
  coef_vector *value )
{   
  // the parameters are really evap_params_1d
  evap_params_1d *Eout_params = static_cast< evap_params_1d* >( E_out_param );
  Eout_params->func_count += 1;
  double Prob = E_out*exp( - E_out/Eout_params->Theta );
  if( ( value->conserve == NUMBER ) || ( value->conserve == BOTH ) )
  {
    value->weight_1[ 0 ] = Prob;
  }
  if( ( value->conserve == ENERGY ) || ( value->conserve == BOTH ) )
  {
    value->weight_E[ 0 ] = E_out*Prob;
  }
}
