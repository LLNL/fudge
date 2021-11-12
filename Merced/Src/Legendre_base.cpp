/*
* ******** merced: calculate the transfer matrix *********
* $Revision: 1 $
* $Date: 2021-10-09 $
* $Author: hedstrom $
* $Id: Legendre_base.cpp 1 2021-10-09Z hedstrom $
*
* ******** merced: calculate the transfer matrix *********
*
* # <<BEGIN-copyright>>
* # <<END-copyright>>
*/

// Implementation for the base class for Legendre data

#include <iostream>

#include "Legendre_base.hpp"
#include "messaging.hpp"

// *************** class LgBase::Legendre_base **********************
// --------------- LgBase::Legendre_base::clean_data ----------------
void LgBase::Legendre_base::clean_data( )
{
  if( order >= 0 )
  {
    delete [] data;
  }
  order = -1;
}

// --------------- LgBase::Legendre_base::initialize ----------------
// Sets the incident energy and allocates space
void LgBase::Legendre_base::initialize( int Order )
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

// -------------- LgBase::Legendre_base::operator[ ] ----------------
// access routine
double& LgBase::Legendre_base::operator[ ]( int N )
{
  if( ( N < 0 ) || ( N > order ) )
  {
    Msg::FatalError( "LgBase::Legendre_base::operator[ ]",
		 "index out of range" );
  }
  return data[ N ];
}

// ------------------ LgBase::Legendre_base::value ----------------
// access routine
double LgBase::Legendre_base::value( int N ) const
{
  if( ( N < 0 ) || ( N > order ) )
  {
    Msg::FatalError( "LgBase::Legendre_base::value",
		 "index out of range" );
  }
  return data[ N ];
}

// -------------- LgBase::Legendre_base::truncate_zeros -------------
// Ignore zero high-order Legendre coefficients
void LgBase::Legendre_base::truncate_zeros( )
{
  int N = order;
  for( N = order; N > 0; --N )
  {
    if( data[ N ] != 0.0 ) break;
  }
  order = N;
}

// ------------- LgBase::Legendre_base::operator*= ----------------
// Scales the vector
LgBase::Legendre_base& LgBase::Legendre_base::operator*=( double factor )
{
  for( int L_count = 0; L_count <= order; ++L_count )
  {
    data[ L_count ] *= factor;
  }
  return *this;
}

// --------------- LgBase::Legendre_base::sum_Legendre ---------------
// Sums the Legenre series
double LgBase::Legendre_base::sum_Legendre( double mu )
{
  if( order < 0 )
  {
    Msg::FatalError( "LgBase::Legendre_data::sum_Legendre",
		     "order not set" );
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

// ------------------ LgBase::Legendre_base::print ----------------
// For debugging
void LgBase::Legendre_base::print( ) const
{
  std::cout << "E " << Energy << ":";
  for( int L_count = 0; L_count <= order; ++L_count )
  {
    std::cout << " " << value( L_count );
  }
  std::cout << std::endl;
}
