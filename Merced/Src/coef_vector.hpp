/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: coef_vector.hpp 1 2009-09-01 03:06:56Z hedstrom $
 * initial code
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/

// header for the coef_vector class
#ifndef COEF_VECTOR
#define COEF_VECTOR

#include "max_output_order.hpp"
#include "Legendre_data.hpp"

//! Do we conserve particle number, energy, or both?
enum Conserve{ NUMBER, ENERGY, BOTH, NOT_SET };

// --------------- coef_vector -------------------------
//! Class to hold Legendre coefficients j=0, 1, ..., order
class coef_vector
{
private:

public:
  //! Coefficients for conservation of number of particles
  double weight_1[ max_output_order+1 ];

  //! Coefficients for conservation of energy
  double weight_E[ max_output_order+1 ];

  //! actual Legendre order
  int order;

  //! flag for conservation of particle number or energy
  Conserve conserve;

  //! Constructor
  inline coef_vector( ): order( -1 ), conserve( NOT_SET ) {}

  //! Constructor
  coef_vector( int Order, Conserve cons );

  //! Destructor
  inline ~coef_vector() {}

  //! Sets the order and conservation flag
  //! \param Order the Legendre order of the output
  //! \param cons the conservation flag
  void set_order( int Order, Conserve cons );

  //! Sets the entries to zero
  void set_zero( );

  //! Ensure that our rough estimate is nonzero
  //! param Norm_1 the max norm of weight_1, set to length if initially 0
  //! param Norm_E the max norm of weight_E, set to length if initially 0
  //! param length default value for norms initially 0
  void test_zero( double *Norm_1, double *Norm_E, double length );

  //! Makes a copy
  //! param to_copy the values to copy
  void copy( const coef_vector &to_copy );

  //! Makes a copy
  //! param to_copy the values to copy
  coef_vector& operator=( const coef_vector& to_copy );

  //! Does vector addition
  //! param to_add the values to add
  coef_vector& operator+=( const coef_vector& to_add );

  //! Adds a scalar to the terms of order L_order
  //! param to_add the values to add (one order)
  //! param L_order the Legendre order to add to
  coef_vector& plus( const coef_vector& to_add, int L_order );

  //! Scales the vector
  //! param factor the scale factor
  coef_vector& operator*=( double factor );

  //! Scales the weight_E terms
  //! param factor the scale factor for weight_E
  void scale_E( double factor );

  //! Weights the vector 
  //! param 
  coef_vector& operator*=( Legendre_base &factor );

  //! Calculates the max norms of weight_1 and weight_E
  //! Replaces the entries by their absolute values
  //! param Norm_1 the max norm of weight_1
  //! param Norm_E the max norm of weight_E
  void max_norm( double *Norm_1, double *Norm_E );

  void print( );
};

#endif
