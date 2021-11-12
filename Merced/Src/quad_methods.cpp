/*
 * ******** merced: calculate the transfer matrix ********
 * $Revision: 601 $
 * $Date: 2018-05-04 $
 * $Author: hedstrom $
 * $Id: quad_methods.cpp 601  2018-05-04Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */

// Classes for Gaussian quadrature

#include <cmath>

#include "quad_methods.hpp"

// ******************** Qmeth::Quadrature_Rule *********************
// ------------------- Qmeth::Quadrature_Rule::operator= ---------------
// to copy
Qmeth::Quadrature_Rule& Qmeth::Quadrature_Rule::operator=( const Qmeth::Quadrature_Rule &to_copy )
{
  quad_method = to_copy.quad_method;
  adaptive = to_copy.adaptive;
  input_set = to_copy.input_set;
  return *this;
}

// ********** basic Gaussian quadrature routines **************
// ************* Qmeth::midpoint_rule *******************************
// The midpoint rule on the interval (A, B).
bool Qmeth::midpoint_rule( bool (*F)( double x, Qparam::QuadParamBase *params, Coef::coef_vector *Value ),
 		double A, double B, Qparam::QuadParamBase *params, Coef::coef_vector *value )
{
  value->set_zero( );
  double M = ( A + B ) / 2;  // midpoint
  bool is_OK = F( M, params, value );
  *value *= B - A;
  return is_OK;
}

// ************* Qmeth::Gauss_2 *******************************
// Second-order Gaussian quadrature on the interval (A, B).
bool Qmeth::Gauss_2( bool (*F)( double x, Qparam::QuadParamBase *params, Coef::coef_vector *Value ), 
		double A, double B, Qparam::QuadParamBase *params, Coef::coef_vector *value )
{
  value->set_zero( );

  // Gauss quadrature points
  const double x1 = 1.0/std::sqrt( 3.0 );
  double M = ( A + B ) / 2;  // midpoint
  double scale = B - M;
  double abcissae[ 2 ] = { M - scale * x1,
			   M + scale * x1 };
  // weights
  const double wt1 = 0.5;
  double weights[ 2 ] = { wt1, wt1 };

  // to hold the function values
  Coef::coef_vector wert( value->order, value->conserve );
  for( int x_count = 0; x_count < 2; ++x_count )
  {
    bool is_OK = F( abcissae[ x_count ], params, &wert );
    if( !is_OK )
    {
      return false;
    }
    wert *= weights[ x_count ];
    *value += wert;
  }
  *value *= B - A;
  return true;
}

// ************* Qmeth::Gauss_3 *******************************
// Third-order Gaussian quadrature on the interval (A, B).
bool Qmeth::Gauss_3( bool (*F)( double x, Qparam::QuadParamBase *params, Coef::coef_vector *Value ), 
		double A, double B, Qparam::QuadParamBase *params, Coef::coef_vector *value )
{
  value->set_zero( );

  // Gauss quadrature points
  const double x1 =std::sqrt( 0.6 );
  double M = ( A + B ) / 2;  // midpoint
  double scale = B - M;
  double abcissae[ 3 ] = { M - scale * x1,
			   M,
			   M + scale * x1 };
  // weights
  const double wt0 = 4.0/9.0;
  const double wt1 = 5.0/18.0;
  double weights[ 3 ] = { wt1, wt0, wt1 };

  // to hold the function values
  Coef::coef_vector wert( value->order, value->conserve );
  for( int x_count = 0; x_count < 3; ++x_count )
  {
    bool is_OK = F( abcissae[ x_count ], params, &wert );
    if( !is_OK )
    {
      return false;
    }
    wert *= weights[ x_count ];
    *value += wert;
  }
  *value *= B - A;
  return true;
}

// ************* Qmeth::Gauss_4 *******************************
// Fourth-order Gaussian quadrature on the interval (A, B).
bool Qmeth::Gauss_4( bool (*F)( double x, Qparam::QuadParamBase *params, Coef::coef_vector *Value ), 
		double A, double B, Qparam::QuadParamBase *params, Coef::coef_vector *value )
{
  value->set_zero( );

  // Gauss quadrature points
  const double x1 = 0.339981043584856;
  const double x2 = 0.861136311594053;
  double M = ( A + B ) / 2;  // midpoint
  double scale = B - M;
  double abcissae[ 4 ] = { M - scale * x2,
			   M - scale * x1,
			   M + scale * x1,
			   M + scale * x2 };
  // weights
  const double wt1 = 0.5*0.652145154862546;
  const double wt2 = 0.5*0.347854845137454;
  double weights[ 4 ] = { wt2, wt1, wt1, wt2 };

  // to hold the function values
  Coef::coef_vector wert( value->order, value->conserve );
  for( int x_count = 0; x_count < 4; ++x_count )
  {
    bool is_OK = F( abcissae[ x_count ], params, &wert );
    if( !is_OK )
    {
      return false;
    }
    wert *= weights[ x_count ];
    *value += wert;
  }
  *value *= B - A;
  return true;
}

// ************* Qmeth::Gauss_6 *******************************
// Sixth-order Gaussian quadrature on the interval (A, B).
bool Qmeth::Gauss_6( bool (*F)( double x, Qparam::QuadParamBase *params, Coef::coef_vector *Value ), 
		double A, double B, Qparam::QuadParamBase *params, Coef::coef_vector *value )
{
  value->set_zero( );

  // Gauss quadrature points from "Gaussian quadrature formulas" by Stroud and Secrest
  const double x1 = 0.238619186083197;
  const double x2 = 0.661209386466265;
  const double x3 = 0.932469514203152;
  double M = ( A + B ) / 2;  // midpoint
  double scale = B - M;
  double abcissae[ 6 ] = { M - scale * x3,
			   M - scale * x2,
			   M - scale * x1,
			   M + scale * x1,
			   M + scale * x2,
			   M + scale * x3 };
  // weights
  const double wt1 = 0.5*0.467913934572691;
  const double wt2 = 0.5*0.360761573048139;
  const double wt3 = 0.5*0.171324492379170;
  double weights[ 6 ] = { wt3, wt2, wt1, wt1, wt2, wt3 };

  // to hold the function values
  Coef::coef_vector wert( value->order, value->conserve );
  for( int x_count = 0; x_count < 6; ++x_count )
  {
    bool is_OK = F( abcissae[ x_count ], params, &wert );
    if( !is_OK )
    {
      return false;
    }
    wert *= weights[ x_count ];
    *value += wert;
  }
  *value *= B - A;
  return true;
}

// ************* Qmeth::Gauss_10 *******************************
// Tenth-order Gaussian quadrature on the interval (A, B).
bool Qmeth::Gauss_10( bool (*F)( double x, Qparam::QuadParamBase *params, Coef::coef_vector *Value ), 
		double A, double B, Qparam::QuadParamBase *params, Coef::coef_vector *value )
{
  value->set_zero( );

  // Gauss quadrature points from "Gaussian quadrature formulas" by Stroud and Secrest
  const double x1 = 0.148874338981631;
  const double x2 = 0.433395394129247;
  const double x3 = 0.679409568299024;
  const double x4 = 0.865063366688985;
  const double x5 = 0.973906528517172;

  double M = ( A + B ) / 2;  // midpoint
  double scale = B - M;
  double abcissae[ 10 ] = { M - scale * x5,
			   M - scale * x4,
			   M - scale * x3,
			   M - scale * x2,
			   M - scale * x1,
			   M + scale * x1,
			   M + scale * x2,
			   M + scale * x3,
			   M + scale * x4,
			   M + scale * x5 };
  // weights
  const double wt1 = 0.5*0.295524224714753;
  const double wt2 = 0.5*0.269266719309996;
  const double wt3 = 0.5*0.219086362515982;
  const double wt4 = 0.5*0.149451349150581;
  const double wt5 = 0.5*0.0666713443086881;
  double weights[ 10 ] = { wt5, wt4, wt3, wt2, wt1, wt1, wt2, wt3, wt4, wt5 };

  // to hold the function values
  Coef::coef_vector wert( value->order, value->conserve );
  for( int x_count = 0; x_count < 10; ++x_count )
  {
    bool is_OK = F( abcissae[ x_count ], params, &wert );
    if( !is_OK )
    {
      return false;
    }
    wert *= weights[ x_count ];
    *value += wert;
  }
  *value *= B - A;
  return true;
}

// ************* Qmeth::weight_L1 ******************************
// First-order Gaussian quadrature on the interval $(A, B)$ with sqrt{x} singularity
bool Qmeth::weight_L1( bool (*F)( double x, Qparam::QuadParamBase *params, Coef::coef_vector *Value ), 
		double A, double B, Qparam::QuadParamBase *params, Coef::coef_vector *value )
{
  value->set_zero( );

  // The evaluation point is at
  // M = 0.6 * B * Z_4 / Z_2 with
  // zeta = sqrt( A / B ) and
  // Z_n = \sum_{j=0}^n zeta^j
  double zeta = std::sqrt( A / B );
  double Z_2 = 1 + zeta*( 1 + zeta );
  double Z_4 = 1 + zeta*( 1 + zeta*Z_2 );
  double M = 0.6 * B * Z_4 / Z_2;
  bool is_OK = F( M, params, value );
  // the weight is 2*( B^{3/2} - A^{3/2} ) / 3
  // but also account for the weight \sqrt{ M }
  *value *= 2 * B * std::sqrt(B) * ( 1 - zeta ) * Z_2 / ( 3 * std::sqrt( M ) );
  return is_OK;
}
