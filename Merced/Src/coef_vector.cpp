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


// *************** class Coef::coef_vector *********************
// ------------------- Coef::coef_vector constructor ------------------------
Coef::coef_vector::coef_vector( int Order, Conserve cons )
{
  set_order( Order, cons );
}
// ------------------- Coef::coef_vector::set_order ------------------------
void Coef::coef_vector::set_order( int Order, Conserve cons )
{
  conserve = cons;
  order = Order;
  set_zero( );
}
// ------------------- Coef::coef_vector::set_zero ------------------------
void Coef::coef_vector::set_zero( )
{
  // Set the entries to zero
  if( ( conserve == Coef::NUMBER ) || ( conserve == Coef::BOTH ) ){
    for( int i = 0; i <= order; ++i )
    {
      weight_1[i] = 0.0;
    }
  }
  if( ( conserve == Coef::ENERGY ) || ( conserve == Coef::BOTH ) ){
    for( int i = 0; i <= order; ++i )
    {
      weight_E[i] = 0.0;
    }
  }
}
// ------------------- Coef::coef_vector::copy ------------------------
void Coef::coef_vector::copy( const Coef::coef_vector &to_copy )
{
  order = to_copy.order;
  conserve = to_copy.conserve;

  if( ( conserve == Coef::NUMBER ) || ( conserve == Coef::BOTH ) )
  {
    for( int i = 0; i <= order; ++i )
    {
      weight_1[i] = to_copy.weight_1[i];
    }
  }
  if( ( conserve == Coef::ENERGY ) || ( conserve == Coef::BOTH ) )
  {
    for( int i = 0; i <= order; ++i )
    {
      weight_E[i] = to_copy.weight_E[i];
    }
  }
}
// ------------------- Coef::coef_vector::operator= ------------------------
Coef::coef_vector& Coef::coef_vector::operator=( const Coef::coef_vector &to_copy )
{
  order = to_copy.order;
  conserve = to_copy.conserve;

  if( ( conserve == Coef::NUMBER ) || ( conserve == Coef::BOTH ) ){
    for( int i = 0; i <= order; ++i )
    {
      weight_1[i] = to_copy.weight_1[i];
    }
  }
  if( ( conserve == Coef::ENERGY ) || ( conserve == Coef::BOTH ) ){
    for( int i = 0; i <= order; ++i )
    {
      weight_E[i] = to_copy.weight_E[i];
    }
  }
  return *this;
}
// ------------------- Coef::coef_vector::operator+= ------------------------
Coef::coef_vector& Coef::coef_vector::operator+=( const Coef::coef_vector& to_add )
{
  if( to_add.order != order ){
    Msg::FatalError( "Coef::coef_vector::operator+=",
		     "incompatible orders" );
  }
  else if( to_add.conserve != conserve ){
    Msg::FatalError( "Coef::coef_vector::operator+=",
		     "incompatible conserve" );
  }

  if( ( conserve == Coef::NUMBER ) || ( conserve == Coef::BOTH ) ){
    for( int i = 0; i <= order; ++i )
    {
      weight_1[i] += to_add.weight_1[i];
    }
  }
  if( ( conserve == Coef::ENERGY ) || ( conserve == Coef::BOTH ) ){
    for( int i = 0; i <= order; ++i )
    {
      weight_E[i] += to_add.weight_E[i];
    }
  }
  return *this;
}
// ------------------- Coef::coef_vector::plus ------------------------
// Adds a scalar to the terms of order L_order
Coef::coef_vector& Coef::coef_vector::plus( const Coef::coef_vector& to_add, int L_order )
{
  if( to_add.order != 0 ){
    Msg::FatalError( "Coef::coef_vector::plus",
		     "wrong order for summand" );
  }
  else if( to_add.conserve != conserve ){
    Msg::FatalError( "Coef::coef_vector::plus",
		     "incompatible conserve" );
  }

  if( ( conserve == Coef::NUMBER ) || ( conserve == Coef::BOTH ) ){
    weight_1[ L_order ] += to_add.weight_1[ 0 ];
  }
  if( ( conserve == Coef::ENERGY ) || ( conserve == Coef::BOTH ) ){
    weight_E[ L_order ] += to_add.weight_E[ 0 ];
  }
  return *this;
}
// ------------------- Coef::coef_vector::operator*= ------------------------
Coef::coef_vector& Coef::coef_vector::operator*=( double factor )
{
  if( ( conserve == Coef::NUMBER ) || ( conserve == Coef::BOTH ) ){
    for( int i = 0; i <= order; ++i )
    {
      weight_1[i] *= factor;
    }
  }
  if( ( conserve == Coef::ENERGY ) || ( conserve == Coef::BOTH ) ){
    for( int i = 0; i <= order; ++i )
    {
      weight_E[i] *= factor;
    }
  }
  return *this;
}
// ----------------- Coef::coef_vector::scale_E ---------------------
// Scales the weight_E terms
void Coef::coef_vector::scale_E( double factor )
{
  for( int i = 0; i <= order; ++i )
  {
    weight_E[i] *= factor;
  }
}
// --------------- Coef::coef_vector::operator*= --------------------
Coef::coef_vector& Coef::coef_vector::operator*=( LgBase::Legendre_base &factor )
{
  if( ( conserve == Coef::NUMBER ) || ( conserve == Coef::BOTH ) ){
    for( int i = 0; i <= order; ++i )
    {
      weight_1[i] *= factor[i];
    }
  }
  if( ( conserve == Coef::ENERGY ) || ( conserve == Coef::BOTH ) ){
    for( int i = 0; i <= order; ++i )
    {
      weight_E[i] *= factor[i];
    }
  }
  return *this;
}
// ------------------- Coef::coef_vector::max_norm ------------------------
// Calculates the max norms of weight_1 and weight_E
void Coef::coef_vector::max_norm( double *Norm_1, double *Norm_E )
{
  *Norm_1 = 0.0;
  if( ( conserve == Coef::NUMBER ) || ( conserve == Coef::BOTH ) )
  {
    *Norm_1 = std::abs( weight_1[0] );
  }
  
  *Norm_E = 0.0;
  if( ( conserve == Coef::ENERGY ) || ( conserve == Coef::BOTH ) )
  {
    *Norm_E= std::abs( weight_E[0] );
  }
}
// ------------------- Coef::coef_vector::test_zero ------------------------
// Makes sure that no entry is zero
void Coef::coef_vector::test_zero( )
{
  // get the max norms
  double Norm_1;
  double Norm_E;
  max_norm( &Norm_1, &Norm_E );

  if( ( conserve == Coef::NUMBER ) || ( conserve == Coef::BOTH ) )
  {
    if( Norm_1 <= 0.0 )
    {
      weight_1[0] = 1.0;
    }
  }

  if( ( conserve == Coef::ENERGY ) || ( conserve == Coef::BOTH ) )
  {
    if( Norm_E <= 0.0 )
    {
      weight_E[0] = 1.0;
    }
  }
}
// ------------------- Coef::coef_vector::print ------------------------
void Coef::coef_vector::print( )
{
  if( ( conserve == Coef::NUMBER ) || ( conserve == Coef::BOTH ) ){
    std::cout << " weight 1: ";
    for( int i = 0; i <= order; ++i )
    {
      std::cout << weight_1[i] << "  ";
    }
  }
  if( ( conserve == Coef::ENERGY ) || ( conserve == Coef::BOTH ) ){
    std::cout << " weight E: ";
    for( int i = 0; i <= order; ++i )
    {
      std::cout << weight_E[i] << "  ";
    }
  }
  std::cout << std::endl;
}
