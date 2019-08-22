/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: Legendre_data.cpp 1 2006-02-02 03:06:56Z hedstrom $
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

// routines for the list of (incident energy, flux)

#include <cmath>
#include "Legendre_data.hpp"
#include "messaging.hpp"
#include "global_params.hpp"

// *************** class Legendre_base *************************
// ------------------ Legendre_base::clean_data ----------------
void Legendre_base::clean_data( )
{
  if( order >= 0 )
  {
    delete [] data;
  }
  order = -1;
}
// ------------------ Legendre_base::initialize ----------------
// Sets the incident energy and allocates space
void Legendre_base::initialize( int Order )
{
  if( order != Order )
  {
    if( order >= 0 )
    {
      clean_data( );
    }
    order = Order;
    data = new double[ order + 1 ];
    for( int L_count = 0; L_count <= order; ++L_count )
    {
      data[ L_count ] = 0.0;
    }
  }
}
// ------------------ Legendre_base::operator[ ] ----------------
// access routine
double& Legendre_base::operator[ ]( int N )
{
  if( ( N < 0 ) || ( N > order ) )
  {
    FatalError( "Legendre_base::operator[ ]",
		 "index out of range" );
  }
  return data[ N ];
}
// ------------------ Legendre_base::value ----------------
// access routine
double Legendre_base::value( int N ) const
{
  if( ( N < 0 ) || ( N > order ) )
  {
    FatalError( "Legendre_base::value",
		 "index out of range" );
  }
  return data[ N ];
}
// ------------------ Legendre_base::truncate_zeros ----------------
// Ignore zero high-order Legendre coefficients
void Legendre_base::truncate_zeros( )
{
  int N = order;
  for( N = order; N > 0; --N )
  {
    if( data[ N ] != 0.0 ) break;
  }
  order = N;
}
// ------------------ Legendre_base::operator*= ----------------
// Scales the vector
Legendre_base& Legendre_base::operator*=( double factor )
{
  for( int L_count = 0; L_count <= order; ++L_count )
  {
    data[ L_count ] *= factor;
  }
  return *this;
}
// --------------- Legendre_base::sum_Legendre ---------------
// Sums the Legenre series
double Legendre_base::sum_Legendre( double mu )
{
  if( order < 0 )
  {
    FatalError( "Legendre_data::sum_Legendre", "order not set" );
  }
  if( mu > 1.0 )
  {
    mu = 1.0;
  }
  if( mu < -1.0 )
  {
    mu = -1.0;
  }
  int ell = 0;
  double this_coef = value( ell );
  double sum = this_coef/2;
  if( order == 0 ) return sum;

  double prevP = 1.0;  // P_0
  double thisP = mu;   // P_1
  double nextP;        // P_{ell+1}
  ell = 1;
  this_coef = value( ell );
  sum += 1.5*this_coef*thisP;
  if( order == 1 ) return sum;

  for( ell = 1; ell < order; ++ell )
  {
    this_coef = value( ell + 1 );
    nextP = ( ( 2*ell + 1 )*mu*thisP - ell*prevP )/(ell + 1.0);
    sum += ( ell + 1.5 )*this_coef*nextP;
    prevP = thisP;
    thisP = nextP;
  }
  return sum;
}
// ------------------ Legendre_base::print ----------------
// For debugging
void Legendre_base::print( ) const
{
  cout << "E " << Energy << ":";
  for( int L_count = 0; L_count <= order; ++L_count )
  {
    cout << " " << value( L_count );
  }
  cout << endl;
}

// *************** class Legendre_coefs *************************
// ------------------ Legendre_coefs::initialize ----------------
// Allocates space
void Legendre_coefs::initialize( int order )
{
  Legendre_base::initialize( order );
}
// ------------------ Legendre_coefs::zero_data ----------------
// Sets all coefficients to zero
void Legendre_coefs::zero_data( )
{
  for( int L_count = 0; L_count <= order; ++L_count )
  {
    data[ L_count ] = 0.0;
  }
}
// ------------------ Legendre_coefs::copy_coef ----------------
// Copies the Legendre coefficients
void Legendre_coefs::copy_coef( const Legendre_coefs& to_copy )
{
  // reset the order if necessary
  if( order != to_copy.order )
  {
    if( order >= 0 )
    {
      clean_data( );
    }
    order = to_copy.order;
    data = new double[ order + 1 ];
  }
  for( int L_count = 0; L_count <= order; ++L_count )
  {
    data[ L_count ] = to_copy.data[ L_count ];
  }
}
// ------------------ Legendre_coefs::only_copy_coef ----------------
// Copies the Legendre coefficients; do not change the order
void Legendre_coefs::only_copy_coef( const Legendre_coefs& to_copy )
{
  // make sure that space is allocated
  if( order < 0 )
  {
    FatalError( "Legendre_coefs::only_copy_coef", "no space allocated" );
  }
  int use_order = ( order <= to_copy.order ) ? order : to_copy.order;
  for( int L_count = 0; L_count <= use_order; ++L_count )
  {
    data[ L_count ] = to_copy.data[ L_count ];
  }
}
// ------------------ Legendre_coefs::set_max_order ----------------
// Sets the order for interpolated data
void Legendre_coefs::set_max_order( int left_order, int right_order )
{
  int max_order = ( left_order >= right_order ) ? left_order : right_order;
  if( order != max_order )
  {
    if( order >= 0 )
    {
      clean_data( );
    }
    order = max_order;
    data = new double[ max_order + 1 ];
  }
}
// ------------------ Legendre_coefs::basic_linlin_interp ----------------
// Interpolates the flux with weight alpha
void Legendre_coefs::basic_linlin_interp( double alpha,
  const Legendre_coefs& prev_flux, const Legendre_coefs& next_flux )
{
  int max_order;
  int min_order;
  bool prev_longer;
  if( prev_flux.order > next_flux.order )
  {
    max_order = prev_flux.order;
    min_order = next_flux.order;
    prev_longer = true;
  }
  else
  {
    min_order = prev_flux.order;
    max_order = next_flux.order;
    prev_longer = false;
  }
  if( order < max_order )
  {
    max_order = order;
  }
  if( order < min_order )
  {
    min_order = order;
  }
  int L_count;
  for( L_count = 0; L_count <= min_order; ++L_count )
  {
    (*this)[ L_count ] = ( 1 - alpha ) * prev_flux.data[ L_count ] +
      alpha * next_flux.data[ L_count ];
  }
  if( prev_longer )
  {
    for( L_count = min_order+1; L_count <= max_order; ++L_count )
    {
      (*this)[ L_count ] = ( 1 - alpha ) * prev_flux.data[ L_count ];
    }
  }
  else
  {
    for( L_count = min_order+1; L_count <= max_order; ++L_count )
    {
      (*this)[ L_count ] = alpha * next_flux.data[ L_count ];
    }
  }
  for( L_count = max_order+1; L_count <= order; ++L_count )
  {
    (*this)[ L_count ] = 0.0;
  }
}
// ------------------ Legendre_coefs::linlin_interp ----------------
// Interpolates the flux at energy E_in
void Legendre_coefs::linlin_interp( double E_in, const Legendre_coefs& prev_flux,
  const Legendre_coefs& next_flux )
{
  if( order < 0 )
  {
    FatalError( "Legendre_coefs::linlin_interp", "no space for data allocated" );
  }
  set_E_in( E_in );
  // linearly interpolate
  double denom = next_flux.get_E_in( ) - prev_flux.get_E_in( );
  if( denom <= 0.0 )
  {
    FatalError( "Legendre_coefs::linlin_interp", "incident energies out of order" );
  }
  double alpha = ( E_in - prev_flux.get_E_in( ) ) / denom;
  basic_linlin_interp( alpha, prev_flux, next_flux );
}
// ------------------ Legendre_coefs::Ein_linlin_interp ----------------
// Interpolates Legendre-coefficient data at energy E_in
void Legendre_coefs::Ein_linlin_interp( double E_in, double left_Ein, 
    const Legendre_coefs& left_flux, double right_Ein,
    const Legendre_coefs& right_flux )
{
  set_E_in( left_flux.get_E_in( ) );
  // linearly interpolate
  double denom = right_Ein - left_Ein;
  if( denom <= 0.0 )
  {
    FatalError( "Legendre_coefs::Ein_linlin_interp", "incident energies out of order" );
  }
  double alpha = ( E_in - left_Ein ) / denom;
  basic_linlin_interp( alpha, left_flux, right_flux );
}
// ------------------ Legendre_coefs::unitbase_interp ----------------
// Interpolates unit-base Legendre-coefficient data
void Legendre_coefs::unitbase_interp( double E_in, double alpha, 
    const Legendre_coefs& left_flux,
    const Legendre_coefs& right_flux )
{
  set_E_in( left_flux.get_E_in( ) );
  basic_linlin_interp( alpha, left_flux, right_flux );
}
// ------------------ Legendre_coefs::linlog_interp ----------------
// Interpolates the flux linearly with respect to the logarithm of the energy
void Legendre_coefs::linlog_interp( double E_in, const Legendre_coefs& prev_flux,
  const Legendre_coefs& next_flux )
{
  if( order < 0 )
  {
    FatalError( "Legendre_coefs::linlog_interp", "no space for data allocated" );
  }
  set_E_in( E_in );
  // interpolate
  double denom = log( next_flux.get_E_in( ) / prev_flux.get_E_in( ) );
  if( denom <= 0.0 )
  {
    FatalError( "Legendre_coefs::linlog_interp", "incident energies out of order" );
  }
  int max_order;
  int min_order;
  bool prev_longer;
  if( prev_flux.order > next_flux.order )
  {
    max_order = prev_flux.order;
    min_order = next_flux.order;
    prev_longer = true;
  }
  else
  {
    min_order = prev_flux.order;
    max_order = next_flux.order;
    prev_longer = false;
  }
  if( order < max_order )
  {
    max_order = order;
  }
  if( order < min_order )
  {
    min_order = order;
  }
  int L_count;
  double alpha = log( E_in / prev_flux.get_E_in( ) ) / denom;
  for( L_count = 0; L_count <= min_order; ++L_count )
  {
    (*this)[ L_count ] = ( 1 - alpha ) * prev_flux.data[ L_count ] +
      alpha * next_flux.data[ L_count ];
  }
  if( prev_longer )
  {
    for( L_count = min_order+1; L_count <= max_order; ++L_count )
    {
      (*this)[ L_count ] = ( 1 - alpha ) * prev_flux.data[ L_count ];
    }
  }
  else
  {
    for( L_count = min_order+1; L_count <= max_order; ++L_count )
    {
      (*this)[ L_count ] = alpha * next_flux.data[ L_count ];
    }
  }
  for( L_count = max_order+1; L_count <= order; ++L_count )
  {
    (*this)[ L_count ] = 0.0;
  }
}

// *************** class Legendre_data_range *************************
// ----------- Legendre_data_range::new_Ein --------------
// Sets up a new incident energy
void Legendre_data_range::new_Ein( double Ein, const unit_base_map &ubasemap,
  Interp_Type Eoutinterp )
{
  E_in = Ein;
  Eout_interp = Eoutinterp;
  ubase_map.copy( ubasemap );
}
// ----------- Legendre_data_range::set_data --------------
// Sets up the data for a given range of outgloing energies
void Legendre_data_range::set_data(  const Legendre_coefs &prevdata,
  const Legendre_coefs &nextdata, double Eout_min, double Eout_max )
{
  if( Eout_min >= Eout_max )
  {
    FatalError( "Legendre_data_range::set_data", "improper energy range" );
  }
  static double E_tol = Global.Value( "E_tol" );
  if( Eout_min < ( 1.0 - E_tol )*prevdata.get_E_out( ) )
  {
    FatalError( "Legendre_data_range::set_data", "Eout_min too low" );
  }
  if( Eout_max > ( 1.0 + E_tol )*nextdata.get_E_out( ) )
  {
    FatalError( "Legendre_data_range::set_data", "Eout_max too high" );
  }

  if( Eout_interp == HISTOGRAM )
  {
    prev_data.set_E_out( Eout_min );
    prev_data.copy_coef( prevdata );
    next_data.set_E_out( Eout_max );
  }
  else
  {
    if( Eout_min < ( 1.0 + E_tol )*prevdata.get_E_out( ) )
    {
      prev_data.set_E_out( Eout_min );
      prev_data.copy_coef( prevdata );
    }
    else
    {
      prev_data.set_max_order( prevdata.order, nextdata.order );
      prev_data.linlin_interp( Eout_min, prevdata, nextdata );
    }
    if( Eout_max > ( 1.0 - E_tol )*nextdata.get_E_out( ) )
    {
      next_data.set_E_out( Eout_max );
      next_data.copy_coef( nextdata );
    }
    else
    {
      next_data.set_max_order( prevdata.order, nextdata.order );
      next_data.linlin_interp( Eout_max, prevdata, nextdata );
    }
  }
}
// ----------- Legendre_data_range::ubase_interpolate --------------
// Do unit-base interpolation between incident energies
void Legendre_data_range::ubase_interpolate( double E_in,
  const Legendre_data_range &left_data, const Legendre_data_range &right_data )
{
  Eout_interp = left_data.Eout_interp;
  double left_Ein = left_data.get_E_in( );
  double right_Ein = right_data.get_E_in( );
  double denom = right_Ein - left_Ein;
  if( denom <= 0.0 )
  {
    FatalError( "Legendre_data_range::ubase_interpolate",
                "incident energies out of order" );
  }
  double alpha = ( E_in - left_Ein )/denom;
  if( ( alpha < 0.0 ) || ( alpha > 1.0 ) )
  {
    FatalError( "Legendre_data_range::ubase_interpolate",
                "extrapolation" );
  }
  prev_data.set_max_order( left_data.prev_data.order, right_data.prev_data.order );
  prev_data.unitbase_interp( E_in, alpha, left_data.prev_data,
			     right_data.prev_data );
  if( Eout_interp == HISTOGRAM )
  {
    next_data.set_E_out( left_data.next_data.get_E_out( ) );
  }
  else
  {
    next_data.set_max_order( left_data.next_data.order, right_data.next_data.order );
    next_data.unitbase_interp( E_in, alpha, left_data.next_data,
			       right_data.next_data );
  }
  set_E_in( E_in );
  ubase_map.interpolate( alpha, left_data.ubase_map, right_data.ubase_map );
}
// ----------- Legendre_data_range::to_unit_base --------------
//  Maps from physical variables to unit-base
void Legendre_data_range::to_unit_base( )
{
  double Eout_min = prev_data.get_E_out( );
  double Eout_max = next_data.get_E_out( );
  double scale = Eout_max - Eout_min;
  if( scale <= 0.0 )
  {
    FatalError( "Legendre_data_range::to_unit_base",
                "bad energy range" );
  }
  ubase_map.Eout_min = Eout_min;
  ubase_map.Eout_max = Eout_max;

  prev_data.set_E_out( 0.0 );
  prev_data *= scale;
  next_data.set_E_out( 1.0 );
  next_data *= scale;
}
// ----------- Legendre_data_range::un_unit_base --------------
// Maps from unit-base to physical variables
void Legendre_data_range::un_unit_base( )
{
  double scale = ubase_map.Eout_max - ubase_map.Eout_min;
  if( scale <= 0.0 )
  {
    FatalError( "Legendre_data_range::un_unit_base",
                "bad unit-base map" );
  }
  double phys_Eout = ubase_map.un_unit_base( prev_data.get_E_out( ) );
  prev_data.set_E_out( phys_Eout );
  prev_data *= 1.0/scale;
  phys_Eout = ubase_map.un_unit_base( next_data.get_E_out( ) );
  next_data.set_E_out( phys_Eout );
  next_data *= 1.0/scale;
}
// ----------- Legendre_data_range::Eout_interpolate --------------
// Returns the Legendre coefficients for this outgoing energy
Legendre_coefs Legendre_data_range::Eout_interpolate( double E_out )
{ 
  Legendre_coefs interpolated_data;
  interpolated_data.set_max_order( prev_data.order, next_data.order );
  interpolated_data.linlin_interp( E_out, prev_data, next_data );
  return interpolated_data;
}

// *************** class Legendre_list_base *************************
// ----------- Legendre_list_base::to_unit_base --------------
// Maps the data to unit base
void Legendre_list_base::to_unit_base( )
{
  // find the energy range
  Legendre_list_base::iterator Eprob_ptr = begin( );
  double Eout_min = Eprob_ptr->get_E_out( );
  ubase_map.Eout_min = Eout_min;
  Eprob_ptr = end( );
  --Eprob_ptr;
  double Eout_max = Eprob_ptr->get_E_out( );
  ubase_map.Eout_max = Eout_max;
  double scale = Eout_max - Eout_min;
  if( scale <= 0.0 )
  {
    FatalError( "Legendre_list_base::to_unit_base",
                 "improper energy range" );
  }

  for( Eprob_ptr = begin( ); Eprob_ptr != end( ); ++Eprob_ptr )
  {
    Eprob_ptr->set_E_out( ( Eprob_ptr->get_E_out( ) - Eout_min ) / scale );
    *Eprob_ptr *= scale;
  }
}
// ----------- Legendre_list_base::get_norm --------------
// Finds the total probability
double Legendre_list_base::get_norm( ) const
{
  // coding based on dd_vector::get_norm
  Legendre_list_base::const_iterator entry = begin();
  
  double norm = 0.0;
  double last_E = entry->get_E_out( );
  double last_P = entry->value( 0 );
  double next_E;
  double next_P;
  
// loop through the vector
  for(++entry; entry != end(); ++entry)
  {
    // get the next energy and probability density
    next_E = entry->get_E_out( );
    next_P = entry->value( 0 );

    // accumulate the probability
    if( Eout_interp == LINLIN )
    {
      norm += 0.5*(next_P + last_P)*(next_E - last_E);
    }
    else // if( Eout_interp == HISTOGRAM )
    {
      norm += last_P*(next_E - last_E);
    }

    // save the values for the next interval
    last_E = next_E;
    last_P = next_P;
  }

  return norm;
}
// ----------- Legendre_list_base::renorm --------------
// Normalizes the total probability
void Legendre_list_base::renorm( )
{
  // coding based on dd_vector::renorm
  double norm = get_norm( );
  Legendre_list_base::iterator entry = begin( );

  // error check
  if(norm == 0.0)
  {
    Warning("Legendre_list_base::renorm",
	    pastenum( "zero norm in renorm routine for E_in: ", get_E_in( ) ) );
    // fudge something
    entry->data[ 0 ] = 1.0;
    norm = get_norm( );
  }
  static double etol = 1.0e+4 * Global.Value( "e_tol" );
  static bool norm_good = true;
  if( norm_good && ( abs( norm - 1.0 ) > etol ) )
  {
    Warning("Legendre_list_base::renorm",
            pastenum( "bad norm in renorm routine: ", norm ) );
    norm_good = false;
  }
  // Now, renormalize
  for( ; entry != end( ); ++entry )
  {
    *entry *= 1.0/norm;
  }
  //  cout << "E_in: " << get_E_in( ) << " norm: " << norm << endl;
}
// ----------- Legendre_list_base::print --------------
// For debugging
void Legendre_list_base::print( )
{
  cout << "E_in: " << E_in << endl;
  Legendre_list_base::const_iterator this_entry = begin( );
  for( ; this_entry != end( ); ++this_entry )
  {
    this_entry->print( );
  }
  cout << endl;
}

// *************** class Flux_List *************************
// ------------------ Flux_List::read_flux ----------------
// Constructs the list from the Python data
void Flux_List::read_flux( data_parser &infile, int num_Ein )
{
  // Read the interpolation
  interp = interp_flag_F::read_1d_interpolation( infile );

  order = Global.Value( "outputLegendreOrder" );
  if( order < 0 )
  {
    FatalError( "Flux_List::convert_flux", "Desired Legendre order not set" );
  }

  Flux_List::iterator next_ptr;  // point to the next data
  // read the flux data
  for( int Ein_count = 0; Ein_count < num_Ein; ++Ein_count )
  {
    // what is the next incident energy?
    double this_E = infile.get_next_double( );

    // append a link
    next_ptr = insert( end( ), Legendre_coefs( ) );
    next_ptr->initialize( order );
    next_ptr->set_E_in( this_E );

    int file_order = infile.get_next_int( ) - 1;  // Legendre order of input data
    int save_coefs = ( file_order > order ) ? order : file_order;
    int coef_count;
    double next_flux;
    for( coef_count = 0; coef_count <= save_coefs; ++coef_count )
    {
      next_flux = infile.get_next_double( );
      ( *next_ptr )[ coef_count ] = next_flux;
    }
    // we may need to fill in the higher-order Legendre coefficients
    if( file_order < order )
    {
      for( coef_count = save_coefs+1; coef_count <= order; ++coef_count )
      {
        ( *next_ptr )[ coef_count ] = next_flux;
      }
    }
    // we may need to discard high-order coefficients
    else if( file_order > order )
    {
      for( coef_count = save_coefs+1; coef_count <= file_order; ++coef_count )
      {
        next_flux = infile.get_next_double( );
      }
    }
  }
}
// ------------------ Flux_List::value ----------------
Legendre_coefs Flux_List::value( double E_in, Flux_List::const_iterator &ptr ) const
{
  Flux_List::const_iterator next_ptr = ptr;
  ++next_ptr;  // point to the next data
  if( next_ptr == end( ) )
  {
    FatalError( "Flux_List::value", pastenum( "energy ", E_in) +
		 " out of range" );
  }
  double denom = next_ptr->get_E_in( ) - ptr->get_E_in( );
  if( denom <= 0.0 )
  {
    FatalError( "Flux_List::value", "energies out of order" );
  }

  // linear interpolate
  double alpha = ( E_in - ptr->get_E_in( ) ) / denom;
  Legendre_coefs answer;
  answer.initialize( order );
  answer.set_E_in( E_in );
  for( int L_order = 0; L_order <= order; ++L_order )
  {
    answer[ L_order ] = ( 1 - alpha ) * ptr->data[ L_order ] +
      alpha * next_ptr->data[ L_order ];
  }
  return answer;
}

// *************** class weight_vector *************************
// ------------------ weight_vector::increment ----------------
// Adds the integrals over ( E_left, E_right )
void weight_vector::increment( Flux_List &e_flux, 
  Flux_List::const_iterator this_flux, double E_left, double E_right )
{
  Legendre_coefs left_value = e_flux.value( E_left, this_flux );
  Legendre_coefs right_value = e_flux.value( E_right, this_flux );
  // use the trapezoid rule
  for( int L_count = 0; L_count <= order; ++L_count )
  {
    data[ L_count ] += 0.5 * ( E_right - E_left ) *
      ( left_value[ L_count ] + right_value[ L_count ] );
  }
}
// ------------------ weight_vector::invert ----------------
// take the reciprocals
void weight_vector::invert( )
{
  for( int L_count = 0; L_count <= order; ++L_count )
  {
    if( data[ L_count ] <= 0.0 )
    {
      FatalError( "weight_vector::invert",
                   "negative flux weight" );
    }
    data[ L_count ] = 1.0 / data[ L_count ];
  }
}
