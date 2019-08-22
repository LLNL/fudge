/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2011-01-18 11:06:56 -0800 (Tue, 18 Jan 2011) $
 * $Author: hedstrom $
 * $Id: Legendre2Body.cpp 1 2011-01-18 12:06:56 -0800 hedstrom $
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
// Implement the classes used for ENDF Legendre expansions for discrete 2-body reactions

#ifdef _OPENMP
 #include <omp.h>
#endif

#include "Legendre2Body.hpp"
#include "global_params.hpp"
#include "messaging.hpp"

using namespace std;

// ************* Legendre_param *************
// ---------------- Legendre_param::interpolate ------------------
// Interpolates between two incident energies
void Legendre_param::interpolate( double Ein,
   Legendre_angle_dist::const_iterator prev_coefs,
   Legendre_angle_dist::const_iterator next_coefs )
{
  if( Ein_interp == LINLIN )
  {
    coefs.linlin_interp( Ein, *prev_coefs, *next_coefs );
  }
  else
  {
    coefs.linlog_interp( Ein, *prev_coefs, *next_coefs );
  }
}

// ************* Legendre_angle_dist *************
// ---------------- Legendre_angle_dist::setup_map ------------------
void Legendre_angle_dist::setup_map( )
// set up the map from center-of-mass to laboratory coordinates
{
  // function parameter
  void *params;

  if( use_relativistic )
  {
    relativistic_mass.setup_masses( &particle_info, Q );
    relativistic_mass.get_threshold( );
    threshold = relativistic_mass.threshold;

    // parameters for the relativistic_F functions
    relativistic_param relativistic_map;
    relativistic_map.masses = &relativistic_mass;
    // we need to set the direction cosine
    relativistic_map.mu_cm = 0.0;
    params = static_cast< void * >( &relativistic_map );
    threshold_out = relativistic_F::T_out_lab( threshold, params );
    flip = relativistic_map.find_lowest_bottom( );
  }
  else
  {
    map.set_map( particle_info, Q );
    map.get_threshold( );
    threshold = map.threshold;

    // parameters for the Newtonian_F functions
    Newton_map_param Newton_param;
    Newton_param.masses = &map;
    params = static_cast< void * >( &Newton_param );
    threshold_out = Newtonian_F::T_out_lab( threshold, params );
    flip = Newton_param.find_lowest_bottom( );
  }
}
// ---------------- Legendre_angle_dist::setup_param_map ------------------
void Legendre_angle_dist::setup_param_map( Legendre2d_param *Ein_param )
// set up the Ein_param->map from center-of-mass to laboratory coordinates
{
  if( use_relativistic )
  {
    Ein_param->relativistic_map.masses = &relativistic_mass;
    Ein_param->lower_hits.relativistic_map.masses = &relativistic_mass;
    Ein_param->upper_hits.relativistic_map.masses = &relativistic_mass;
  }
  else
  {
    Ein_param->Newton_map.masses = &map;
    Ein_param->lower_hits.Newton_map.masses = &map;
    Ein_param->upper_hits.Newton_map.masses = &map;
  }
}
// ----------- Legendre_angle_dist::set_threshold --------------
// Uses the mass difference to set the threshold
void Legendre_angle_dist::set_threshold( )
{
  if( use_relativistic )
  {
    threshold = relativistic_mass.threshold;
  }
  else
  {
    threshold = map.threshold;
  }

  // adjust the data if necessary
  Legendre_angle_dist::iterator data_ptr = begin( );
  double first_Ein = data_ptr->get_E_in( );
  if( first_Ein < threshold )
  {
    data_ptr->set_E_in( threshold );
  }
}
// ----------- Legendre_angle_dist::initialize_param --------------
// Initializes the quadrature parameters
void Legendre_angle_dist::initialize_param( Quadrature_Method mu_quad_method, 
  Legendre2d_param *Ein_param )
{
  Ein_param->use_relativistic = use_relativistic;
  Ein_param->lower_hits.use_relativistic = use_relativistic;
  Ein_param->upper_hits.use_relativistic = use_relativistic;
  Ein_param->threshold = threshold;
  Ein_param->threshold_out = threshold_out;
  Ein_param->upper_hits.flip.x = threshold;  // for mu_cm = 1
  Ein_param->upper_hits.flip.y = threshold_out;
  Ein_param->lower_hits.flip.x = flip.x;  // for mu_cm = -1
  Ein_param->lower_hits.flip.y = flip.y;
  Ein_param->Ein_interp = Ein_interp;
  Ein_param->mu_quad_method = mu_quad_method;
}
// ----------- Legendre_angle_dist::read_data --------------
void Legendre_angle_dist::read_data( data_parser &input_file, int num_Ein )
{
  Ein_interp = interp_flag_F::read_1d_interpolation( input_file );

  Legendre_coefs new_angle_dist;  // angular distribution for one E_in
  Legendre_angle_dist::iterator new_angle_ptr;

  // read the data
  for( int Ein_count = 0; Ein_count < num_Ein; ++Ein_count )
  {
    // insert a new angular distribution
    new_angle_ptr = insert( end( ), Legendre_coefs( ) );
    // get the incident energy and the data pairs
    new_angle_ptr->set_E_in( input_file.get_next_double( ) );
    int Order = input_file.get_next_int( ) - 1;
    new_angle_ptr->initialize( Order );
    for( int coef_count = 0; coef_count <= Order; ++coef_count )
    {
      (*new_angle_ptr)[ coef_count ] = input_file.get_next_double( );
    }
  }
  new_angle_ptr->truncate_zeros( );
  //  print( );
}
// ----------- Legendre_angle_dist::get_T --------------
// Calculates the transfer matrix for this particle.
// sigma is the cross section.
void Legendre_angle_dist::get_T( const dd_vector& sigma, const dd_vector& weight,
  T_matrix& transfer )
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

  if( ( Ein_interp != LINLIN ) && ( Ein_interp != LINLOG ) )
  {
    FatalError( "Legendre_angle_dist::get_T",
      "Incident energy interpolation not implemented" );
  }
  if( particle_info.mProd == 0.0 )
  {
    FatalError( "Legendre_angle_dist::get_T", "gamma emission not implemented" );
  }

  // the multiplicity is one
  dd_vector multiple;
  multiple.make_flat( sigma, 1.0 );
  bool done = get_Ein_range( sigma, multiple, weight, transfer.e_flux,
    transfer.in_groups );
  if( done )
  {
    transfer.zero_transfer( );
  }

  // for center-of-mass data
  setup_map( );
  // use the threshold computed from the kinetics
  set_threshold( ); 
  // the computed threshold may be above 20 MeV
   dd_vector::const_iterator last_sigma = sigma.end( );
  --last_sigma;
  if( last_sigma->x < threshold )
  {
    Info( "Legendre_angle_dist::get_T", "computed threshold outside the data range" );
    transfer.zero_transfer( );
  }

  if( ( transfer.Ein_quad_method != ADAPTIVE2 ) &&
      ( transfer.Ein_quad_method != ADAPTIVE4 ) )
  {
    Warning( "Legendre_angle_dist::get_T",
	     "Gaussian quadrature is not advised for this data" );
  }

  int num_negative = 0; // number of negative sums of Legendre series
  long int quad_count = 0;  // number of 2-d quadratures
  long int Ein_F_count= 0;  // number of calls to Legendre2Body_F::E_quad_F
  long int mu_F_count = 0;  // number of calls to Legendre2Body_F::mu_cm_quad_F

  // now do the integrals incident bin by incident bin
#pragma omp parallel for schedule( dynamic, 1 ) default( none )            \
  shared( sigma, multiple, weight, transfer ) \
  reduction( +: quad_count ) reduction( +: Ein_F_count ) \
  reduction( +: mu_F_count ) reduction( +: num_negative )
  for( int Ein_bin = first_Ein; Ein_bin < last_Ein; ++Ein_bin )
  {
    Legendre2d_param Ein_param;
    initialize_param( transfer.mu_quad_method, &Ein_param );
    setup_param_map( &Ein_param );

    // set up the data range for this bin
    Ein_param.setup_bin( Ein_bin, sigma, multiple, weight, transfer.e_flux,
                         transfer.in_groups );
    setup_data( &Ein_param );
    for( ; ; )
    {
      set_Ein_range( &Ein_param );   // get the incident energy interval
      Eout_ladder( transfer, &Ein_param );
      bool done = next_ladder( Ein_param.data_E_1, &Ein_param );   // go to the next interval
      if( done )
      {
        break;
      
      }
    }
    num_negative += Ein_param.num_negative;
    quad_count += Ein_param.quad_count;
    Ein_F_count += Ein_param.Ein_F_count;
    mu_F_count += Ein_param.mu_F_count;

  }
  if( num_negative > 0 )
  {
    Info( "Legendre_angle_dist::get_T", pastenum( "got ", num_negative ) +
	  " negative Legendre sums." );
  }  // end of parallel loop

  // print the counts of function evaluations
  cout << "2d quadratures: " << quad_count << endl;
  cout << "Legendre2Body_F::E_quad_F calls: " << Ein_F_count << endl;
  cout << "Legendre2Body_F::mu_cm_quad_F calls: " << mu_F_count << endl;
  cout << "average Legendre2Body_F::E_quad_F calls: " << 1.0*Ein_F_count/quad_count << endl;
  cout << "average Legendre2Body_F::mu_cm_quad_F calls: " << 1.0*mu_F_count/Ein_F_count << endl;
}
// ----------- Legendre_angle_dist::get_Ein_range --------------
//  Gets the range of nontrivial incident energy bins; computes first_Ein and last_Ein
// returns true if the threshold is too high for the energy bins
bool Legendre_angle_dist::get_Ein_range( const dd_vector& sigma, const dd_vector& mult,
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
  Legendre_angle_dist::const_iterator Ein_data_ptr = begin( );
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
// ----------- Legendre_angle_dist::setup_data --------------
// Initializes the quadrature parameters
void Legendre_angle_dist::setup_data( Legendre2d_param *Ein_param )
{
  Ein_param->left_data = begin( );
  Ein_param->right_data = Ein_param->left_data;
  ++Ein_param->right_data;
  while( Ein_param->right_data->get_E_in( ) <= *Ein_param->Ein_ptr )
  {
    Ein_param->left_data = Ein_param->right_data;
    ++Ein_param->right_data;
  }

  double first_E = Ein_param->left_data->get_E_in( );
  if( first_E > Ein_param->data_E_0 )
  {
    Ein_param->data_E_0 = first_E;
  }
}
// ----------- Legendre_angle_dist::next_ladder --------------
bool Legendre_angle_dist::next_ladder( double E_in, Legendre2d_param *Ein_param )
{
  bool done = Ein_param->update_bin_pointers( E_in );
  static double etol = Global.Value( "E_tol" );
  if( !done )
  {
    double E_tol = E_in * etol;
    if( E_in + E_tol >= Ein_param->right_data->get_E_in( ) )
    {
      while( E_in + E_tol >= Ein_param->right_data->get_E_in( ) )
      {
        // get the next angular data
        Ein_param->left_data = Ein_param->right_data;
        ++Ein_param->right_data;
        if( Ein_param->right_data == end ( ) )
        {
          return true;
        }
      }
    }
  }
  return done;
}

// ----------- Legendre_angle_dist::set_Ein_range --------------
// Sets the range of incident energies for this intergration
void Legendre_angle_dist::set_Ein_range( Legendre2d_param *Ein_param )
{
  Ein_param->set_Ein_range( );
  double this_E = Ein_param->left_data->get_E_in( );
  if( this_E > Ein_param->data_E_0 ) Ein_param->data_E_0 = this_E;
  this_E = Ein_param->right_data->get_E_in( );
  if( this_E < Ein_param->data_E_1 ) Ein_param->data_E_1 = this_E;

  // the data may be below the threshold
 if( Ein_param->data_E_1 < threshold )
  {
    Ein_param->data_E_0 = threshold;
    Ein_param->data_E_1 = threshold;
    // Warning( "Legendre_angle_dist::set_Ein_range", "ignoring 2 data below threshold" );
  }
  else if( Ein_param->data_E_0 < threshold )
  {
    Ein_param->data_E_0 = threshold;
    // Warning( "Legendre_angle_dist::set_Ein_range", "ignoring data below threshold" );
  }

  if( Ein_param->data_E_1 < Ein_param->data_E_0 )
  {
    FatalError( "Legendre_angle_dist::set_Ein_range", "check the incident energies" );
  }
  Ein_param->set_sigma_range( );

  // set the energies for Ein_param->upper_hits and Ein_param->lower_hits
  double left_data_Ein = Ein_param->left_data->get_E_in( );
  Ein_param->upper_hits.left_Ein_Eout.x = left_data_Ein;
  Ein_param->lower_hits.left_Ein_Eout.x = left_data_Ein;
  double right_data_Ein = Ein_param->right_data->get_E_in( );
  Ein_param->upper_hits.right_Ein_Eout.x = right_data_Ein;
  Ein_param->lower_hits.right_Ein_Eout.x = right_data_Ein;
  void *params;  // for the function parameters
  if( use_relativistic )
  {
    Ein_param->relativistic_map.mu_cm = 1.0;
    params = static_cast< void * >( &Ein_param->relativistic_map );
    Ein_param->upper_hits.left_Ein_Eout.y = relativistic_F::T_out_lab( left_data_Ein, params );
    Ein_param->upper_hits.right_Ein_Eout.y = relativistic_F::T_out_lab( right_data_Ein, params );
    // for mu_cm = -1
    Ein_param->relativistic_map.mu_cm = -1.0;
    Ein_param->lower_hits.left_Ein_Eout.y = relativistic_F::T_out_lab( left_data_Ein, params );
    Ein_param->lower_hits.right_Ein_Eout.y = relativistic_F::T_out_lab( right_data_Ein, params );
   }
  else
  {
    Ein_param->Newton_map.mu_cm = 1.0;
    params = static_cast< void * >( &Ein_param->Newton_map );
    Ein_param->upper_hits.left_Ein_Eout.y = Newtonian_F::T_out_lab( left_data_Ein, params );
    Ein_param->upper_hits.right_Ein_Eout.y = Newtonian_F::T_out_lab( right_data_Ein, params );
    // for mu_cm = -1
    Ein_param->Newton_map.mu_cm = -1.0;
    Ein_param->lower_hits.left_Ein_Eout.y = Newtonian_F::T_out_lab( left_data_Ein, params );
    Ein_param->lower_hits.right_Ein_Eout.y = Newtonian_F::T_out_lab( right_data_Ein, params );
  }
}
// ----------- Legendre_angle_dist::Eout_ladder --------------
// This routine uses the this angular Legendre expansion and the
// next to calculate the contribution to the E_out boxes of the
// transfer matrix between incident energies Ein_param->data_E_0 and
// Ein_param->data_E_1.
void Legendre_angle_dist::Eout_ladder( T_matrix& transfer,
   Legendre2d_param *Ein_param )
{
  bool geom_OK;  // for checking the consistency of the geometry
  angle_hit_list test_hits;
  test_hits.Newton_map.masses = &map;
  test_hits.relativistic_map.masses = &relativistic_mass;

  // loop through the outgoing energies (column of transfer)
  for( int Eout_count = 0; Eout_count < transfer.num_Eout_bins;
    ++Eout_count )
  {
    vector< double >::const_iterator Eout_ptr = transfer.out_groups.begin( )
      + Eout_count;
    // how does the mu = -1 hyperbola meet this E-E' box?
    geom_OK = Ein_param->lower_hits.hit_box( -1.0, Eout_ptr,
      Ein_param->data_E_0, Ein_param->data_E_1 );
    if( !geom_OK )
    {
      test_hits.eta = -1.0;
      test_hits.hit_box( -1.0, Eout_ptr,
        Ein_param->data_E_0, Ein_param->data_E_1 );
      test_hits.print( );
      FatalError( "Legendre_angle_dist::Eout_ladder", "Check the coding, 1" );
    }
    if( ( Eout_count < transfer.num_Eout_bins - 1 ) &&
        ( Ein_param->lower_hits.is_above( ) ) )
    {
      // go on to the next E-E' box
      continue;
    }
    // how does the mu = 1 hyperbola meet this E-E' box?
    geom_OK = Ein_param->upper_hits.hit_box( 1.0, Eout_ptr,
      Ein_param->data_E_0, Ein_param->data_E_1 );
    if( !geom_OK )
    {
      test_hits.eta = 1.0;
      test_hits.hit_box( 1.0, Eout_ptr,
                         Ein_param->data_E_0, Ein_param->data_E_1 );
      test_hits.print( );
      FatalError( "Legendre_angle_dist::Eout_ladder", "Check the coding, 2" );
    }
    if( ( Eout_count > 0 ) && ( Ein_param->upper_hits.is_below( ) ) )
    {
      // we are done with this pair of eta values
      break;
    }
    // integrate over this E-E' box
    one_Ebox( transfer, Eout_count,Ein_param  );
  }
}
// ----------- Legendre_angle_dist::one_Ebox --------------
// Integrate over one E-E' box
void Legendre_angle_dist::one_Ebox( T_matrix& transfer, int Eout_count,
   Legendre2d_param *Ein_param )
{
  // the E' energy range
  Ein_param->Eout_min = transfer.out_groups[ Eout_count ];
  Ein_param->Eout_max = transfer.out_groups[ Eout_count + 1 ];

  // set up common incident energies
  Ein_param->lower_hits.common_hits( Ein_param->upper_hits );

  // integrate depending on how the hyperbolas mu_cm = const meet the box
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
      // the lower mu_cm = const hyperbola is inside the E-E' box
      Ein_param->use_Eout_min = false;
      // where is the upper hyperbola?
      if( ( high_hit_ptr->hit_edge == ABOVE ) ||
          ( high_hit_ptr->hit_edge == TOP_OUT ) )
      {
        // integrate up to the top of the E-E' bin
        Ein_param->use_Eout_max = true;
      }
      else
      {
        // integrate up to the next mu_cm = const hyperbola
        Ein_param->use_Eout_max = false;
      }
    }
    else
    {
      // the lower mu_cm = const hyperbola is below the E-E' box;
      // integrate from Eout_min
      Ein_param->use_Eout_min = true;
      // where is the upper mu_cm = const hyperbola?
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
        // the upper mu_cm = const hyperbola is inside the E-E' box
        Ein_param->use_Eout_max = false;
      }
      else
      {
        // the upper mu_cm = const hyperbola is above the E-E' box
        Ein_param->use_Eout_max = true;
      }
    }
    // the range of integration in incident energy
    Ein_param->Ein_0 = low_hit_ptr->E_in;
    Ein_param->Ein_1 = next_low_ptr->E_in;
    update_T( transfer, Eout_count, Ein_param );
  }
}
// ----------- Legendre_angle_dist::update_T --------------
void Legendre_angle_dist::update_T( T_matrix &transfer, int Eout_count,
   Legendre2d_param *Ein_param )
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
    quad_F::integrate( Legendre2Body_F::E_quad_F, transfer.Ein_quad_method,
       left_E, right_E, params, tol, &value );

    // Add this integral
    transfer( Ein_param->Ein_count, Eout_count ) += value;
    // increment the function counts
    Ein_param->Ein_F_count += Ein_param->func_count;
    ++Ein_param->quad_count;
  }
}

// **************** functions to integrate **********
// Function for the 1-d quadrature
// ---------------- Legendre2Body_F::mu_cm_quad_F ------------------
void Legendre2Body_F::mu_cm_quad_F( double mu_cm,
   QuadParamBase *mu_cm_quad_param, coef_vector *value )
{
  // the parameters are really Legendre_param
  Legendre_param *params = static_cast<Legendre_param*>( mu_cm_quad_param );
  params->func_count += 1;
  //  if( params->func_count % 100 == 0 )
  //  {
  //    Info( "mu_cm_quad_F", pastenum( "got ", params->func_count ) + " evaluations");
  //  }
  // get Eout_lab and mu_lab
  double Eout_lab;
  double mu_lab;
  double E_in = params->coefs.get_E_in( );
  if( params->use_relativistic )
  {
    params->relativistic_map->get_E_mu_lab( mu_cm, &Eout_lab, &mu_lab );
  }
  else
  {
    double E_in = params->coefs.get_E_in( );
    params->Newton_map->masses->two_body_get_E_mu_lab( E_in, mu_cm,
      &Eout_lab, &mu_lab );
  }

  // the Legendre polynomials
  math_F::Legendre( mu_lab, value );

  // the probability density
  double Prob = params->coefs.sum_Legendre( mu_cm );
  *value *= Prob;
  if( Prob < 0.0 )
  {
    ++(params->num_negative);
    if( ( !params->flag_set ) && ( params->num_negative == 1 ) )
    {
      string Ein_value = pastenum( "Negative Legendre sum for E_in:", E_in );
      string mu_value = pastenum( " and mu_cm:", mu_cm );
      Info( "Legendre_mu_cm_quad_F", Ein_value + mu_value );
    }
  
  }

  // do the energy weighting if necessary
  if( ( value->conserve == ENERGY ) || ( value->conserve == BOTH ) )
  {
    value->scale_E( Eout_lab );
  }
}
// ---------------- Legendre2Body_F::E_quad_F ------------------
void Legendre2Body_F::E_quad_F( double E_in, QuadParamBase *e_quad_param,
   coef_vector *value )
// Function for the 2-d quadrature
{
  // the parameters are really Legendre2d_param *
  Legendre2d_param *e_params = static_cast<Legendre2d_param *>( e_quad_param );
  e_params->func_count += 1;
  //  if( e_params->func_count % 100 == 0 )
  //  {
  //    Info( "Legendre2Body_F::E_quad_F", pastenum( "got ",
  //      e_params->func_count ) + " evaluations");
  //  }

  // The value of E_quad_F is itself an integral over mu_cm.
  // *value comes in as 0.  

  // parameters for the integration over mu_cm
  Legendre_param mu_cm_params;
  mu_cm_params.flag_set = ( e_params->num_negative > 0 );
  mu_cm_params.coefs.set_E_in( E_in );
  mu_cm_params.Newton_map = &e_params->Newton_map;
  mu_cm_params.relativistic_map = &e_params->relativistic_map;
  mu_cm_params.use_relativistic = e_params->use_relativistic;
  mu_cm_params.Ein_interp = e_params->Ein_interp;
  // interpolate the (mu_cm, probability) with respect to incident energy
  int max_order = ( e_params->left_data->order > e_params->right_data->order ) ?
    e_params->left_data->order : e_params->right_data->order;
  mu_cm_params.coefs.initialize( max_order );
  mu_cm_params.interpolate( E_in, e_params->left_data, e_params->right_data );

  // the range of integration
  double mu_cm_0;
  double mu_cm_1;
  if( e_params->use_relativistic )
  {
    mu_cm_params.relativistic_map->set_boost( E_in );

    mu_cm_0 = ( e_params->use_Eout_min ) ?
      mu_cm_params.relativistic_map->get_mu_cm( e_params->Eout_min ) : -1.0;
    mu_cm_1 = ( e_params->use_Eout_max ) ?
    mu_cm_params.relativistic_map->get_mu_cm( e_params->Eout_max ) : 1.0;
  }
  else
  {
    mu_cm_0 = ( e_params->use_Eout_min ) ?
      mu_cm_params.Newton_map->masses->get_mu_cm( E_in, e_params->Eout_min ) : -1.0;
    mu_cm_1 = ( e_params->use_Eout_max ) ?
      mu_cm_params.Newton_map->masses->get_mu_cm( E_in, e_params->Eout_max ) : 1.0;
  }

  // evaluate the integral over mu_cm
  QuadParamBase *params = static_cast< QuadParamBase* >( &mu_cm_params );
  static double tol = Global.Value( "quad_tol" );

  quad_F::integrate( Legendre2Body_F::mu_cm_quad_F, e_params->mu_quad_method,
         mu_cm_0, mu_cm_1, params, tol, value );

  e_params->num_negative += mu_cm_params.num_negative;
  e_params->mu_F_count += mu_cm_params.func_count;
  // weight it by flux * cross section
  e_params->set_weight( E_in );
  *value *= e_params->current_weight;
  //  cout << "E_in: " << E_in << " mu_cm_0: " << mu_cm_0 << " mu_cm_1: " <<
  //    mu_cm_1 << endl;
  //  value->print( );
}
