/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2009-09-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: coef_vector.cpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */

// implementation of coef_vector class

#include <cmath>

#include "coef_vector.hpp"
#include "messaging.hpp"
#include "global_params.hpp"

using namespace std;

// *************** class coef_vector *********************
// ------------------- coef_vector constructor ------------------------
coef_vector::coef_vector( int Order, Conserve cons )
{
  set_order( Order, cons );
}
// ------------------- coef_vector::set_order ------------------------
void coef_vector::set_order( int Order, Conserve cons )
{
  conserve = cons;
  order = Order;
  set_zero( );
}
// ------------------- coef_vector::set_zero ------------------------
void coef_vector::set_zero( )
{
  // Set the entries to zero
  if( ( conserve == NUMBER ) || ( conserve == BOTH ) ){
    for( int i = 0; i <= order; ++i )
    {
      weight_1[i] = 0.0;
    }
  }
  if( ( conserve == ENERGY ) || ( conserve == BOTH ) ){
    for( int i = 0; i <= order; ++i )
    {
      weight_E[i] = 0.0;
    }
  }
}
// ------------------- coef_vector::test_zero ------------------------
void coef_vector::test_zero( double *Norm_1, double *Norm_E, double length )
{
  // ensure that our rough estimate is nonzero
  if( *Norm_1 <= 0.0 )
  {
    *Norm_1 = length;
  }
  if( ( conserve == NUMBER ) || ( conserve == BOTH ) ){
    for( int i = 0; i <= order; ++i )
    {
      if( weight_1[i] == 0.0 )
      {
	weight_1[i] = *Norm_1;
      }
    }
  }
  if( *Norm_E <= 0.0 )
  {
    *Norm_E = length;
  }
  if( ( conserve == ENERGY ) || ( conserve == BOTH ) ){
    for( int i = 0; i <= order; ++i )
    {
      if( weight_E[i] == 0.0 )
      {
	weight_E[i] = *Norm_E;
      }
    }
  }
}
// ------------------- coef_vector::copy ------------------------
void coef_vector::copy( const coef_vector &to_copy )
{
  order = to_copy.order;
  conserve = to_copy.conserve;

  if( ( conserve == NUMBER ) || ( conserve == BOTH ) )
  {
    for( int i = 0; i <= order; ++i )
    {
      weight_1[i] = to_copy.weight_1[i];
    }
  }
  if( ( conserve == ENERGY ) || ( conserve == BOTH ) )
  {
    for( int i = 0; i <= order; ++i )
    {
      weight_E[i] = to_copy.weight_E[i];
    }
  }
}
// ------------------- coef_vector::operator= ------------------------
coef_vector& coef_vector::operator=( const coef_vector &to_copy )
{
  order = to_copy.order;
  conserve = to_copy.conserve;

  if( ( conserve == NUMBER ) || ( conserve == BOTH ) ){
    for( int i = 0; i <= order; ++i )
    {
      weight_1[i] = to_copy.weight_1[i];
    }
  }
  if( ( conserve == ENERGY ) || ( conserve == BOTH ) ){
    for( int i = 0; i <= order; ++i )
    {
      weight_E[i] = to_copy.weight_E[i];
    }
  }
  return *this;
}
// ------------------- coef_vector::operator+= ------------------------
coef_vector& coef_vector::operator+=( const coef_vector& to_add )
{
  if( to_add.order != order ){
    FatalError( "coef_vector::operator+=", "incompatible orders" );
  }
  else if( to_add.conserve != conserve ){
    FatalError( "coef_vector::operator+=", "incompatible conserve" );
  }

  if( ( conserve == NUMBER ) || ( conserve == BOTH ) ){
    for( int i = 0; i <= order; ++i )
    {
      weight_1[i] += to_add.weight_1[i];
    }
  }
  if( ( conserve == ENERGY ) || ( conserve == BOTH ) ){
    for( int i = 0; i <= order; ++i )
    {
      weight_E[i] += to_add.weight_E[i];
    }
  }
  return *this;
}
// ------------------- coef_vector::plus ------------------------
// Adds a scalar to the terms of order L_order
coef_vector& coef_vector::plus( const coef_vector& to_add, int L_order )
{
  if( to_add.order != 0 ){
    FatalError( "coef_vector::plus", "wrong order for summand" );
  }
  else if( to_add.conserve != conserve ){
    FatalError( "coef_vector::plus", "incompatible conserve" );
  }

  if( ( conserve == NUMBER ) || ( conserve == BOTH ) ){
    weight_1[ L_order ] += to_add.weight_1[ 0 ];
  }
  if( ( conserve == ENERGY ) || ( conserve == BOTH ) ){
    weight_E[ L_order ] += to_add.weight_E[ 0 ];
  }
  return *this;
}
// ------------------- coef_vector::operator*= ------------------------
coef_vector& coef_vector::operator*=( double factor )
{
  if( ( conserve == NUMBER ) || ( conserve == BOTH ) ){
    for( int i = 0; i <= order; ++i )
    {
      weight_1[i] *= factor;
    }
  }
  if( ( conserve == ENERGY ) || ( conserve == BOTH ) ){
    for( int i = 0; i <= order; ++i )
    {
      weight_E[i] *= factor;
    }
  }
  return *this;
}
// ------------------- coef_vector::scale_E ------------------------
// Scales the weight_E terms
void coef_vector::scale_E( double factor )
{
  for( int i = 0; i <= order; ++i )
  {
    weight_E[i] *= factor;
  }
}
// ------------------- coef_vector::operator*= ------------------------
coef_vector& coef_vector::operator*=( Legendre_base &factor )
{
  if( ( conserve == NUMBER ) || ( conserve == BOTH ) ){
    for( int i = 0; i <= order; ++i )
    {
      weight_1[i] *= factor[i];
    }
  }
  if( ( conserve == ENERGY ) || ( conserve == BOTH ) ){
    for( int i = 0; i <= order; ++i )
    {
      weight_E[i] *= factor[i];
    }
  }
  return *this;
}
// ------------------- coef_vector::max_norm ------------------------
// Calculates the max norms of weight_1 and weight_E
void coef_vector::max_norm( double *Norm_1, double *Norm_E )
{
  *Norm_1 = 0.0;
  if( ( conserve == NUMBER ) || ( conserve == BOTH ) ){
    for( int i = 0; i <= order; ++i )
    {
      weight_1[i] =  abs( weight_1[i] );
      if( weight_1[i] > *Norm_1 ) *Norm_1 = weight_1[i];
    }
  }
  *Norm_E = 0.0;
  if( ( conserve == ENERGY ) || ( conserve == BOTH ) ){
    for( int i = 0; i <= order; ++i )
    {
      weight_E[i] = abs( weight_E[i] );
      if( weight_E[i] > *Norm_E ) *Norm_E = weight_E[i];
    }
  }
}
// ------------------- coef_vector::print ------------------------
void coef_vector::print( )
{
  if( ( conserve == NUMBER ) || ( conserve == BOTH ) ){
    cout << " weight 1: ";
    for( int i = 0; i <= order; ++i )
    {
      cout << weight_1[i] << "  ";
    }
  }
  if( ( conserve == ENERGY ) || ( conserve == BOTH ) ){
    cout << " weight E: ";
    for( int i = 0; i <= order; ++i )
    {
      cout << weight_E[i] << "  ";
    }
  }
  cout << endl;
}
