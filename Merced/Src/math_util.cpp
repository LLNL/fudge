/*
 * ******** merced: calculate the transfer matrix ********
 * $Revision: 355 $
 * $Date: 2013-03-15 19:06:56 -0800 (Fri, 15 Mar 2013) $
 * $Author: hedstrom $
 * $Id: math_util.cpp 355  2013-03-15 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */

// implementation of quadrature routine

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include "math_util.hpp"
#include "messaging.hpp"
#include "global_params.hpp"

using namespace std;

// ************* math_F::quadratic *******************************
// Solve the quadratic A*alpha^2 + B*alpha + C = 0; returns the number of real roots
int math_F::quadratic( double A, double B, double C, double *alpha_1,
  double *alpha_2 )
{
  int num_roots;
  if( A == 0.0 )
  {
    if( B == 0.0 )
    {
      if( C == 0.0 )
      {
        num_roots = 3;
        Warning( "quadratic", "reduces to 0 = 0" );
	*alpha_1 = 0.0;
	*alpha_2 = 1.0;
      }
      else
      {
        num_roots = 0;
        Warning( "quadratic", "reduces to 0 = 1" );
      }
    }
    else
    {
      num_roots = 1;
      //      Warning( "quadratic", "only linear terms" );
      *alpha_1 = -C/B;
    }
    return num_roots;
  }

  double root_1;
  double root_2;
  double discriminant = B*B - 4*A*C;
  static double abs_tol = Global.Value( "abs_tol" );

  if( discriminant < -B*B*abs_tol )
  {
    num_roots = 0;
    root_1 = 0.0;
    root_2 = 0.0;
  }
  else if( discriminant < B*B*abs_tol )
  {
    num_roots = 1;
    *alpha_1 = -B/(2*A);
    *alpha_2 = 0.0;
    return num_roots;
  }
  else
  {
    num_roots = 2;
    if( B < 0.0 )
    {
      root_1 = ( -B + sqrt( discriminant ) )/( 2.0*A );
      root_2 = C/( A*root_1 );
    }
    else
    {
      root_1 = ( -B - sqrt( discriminant ) )/( 2.0*A );
      root_2 = C/( A*root_1 );
    }
  }
  // return the roots in increasing order
  if( root_1 < root_2 )
  {
    *alpha_1 = root_1;
    *alpha_2 = root_2;
  }
  else
  {
    *alpha_1 = root_2;
    *alpha_2 = root_1;
  }
  return num_roots;
}

// ************* quad_F::integrate *******************************
// Integrates F over the interval (A, B).
void quad_F::integrate( void (*F)( double x, QuadParamBase *params, coef_vector *Value ), 
		  Quadrature_Method use_quad, double A, double B,
		  QuadParamBase *params, double tol, coef_vector *value )
{
  // we may integrate over mu, so B may be zero
  double abs_AB = 0.2*( abs(A) + abs( B ) );
  value->set_zero( );
  params->func_count = 0;
  // We may have a trivial integral
  if( A >= B ){
    return;
  }

  static double from_quad_tol = Global.Value( "abs_quad_tol" )/100;
  static double from_abs_tol = 1000*Global.Value( "abs_tol" );
  static double skip_tol = ( from_quad_tol > from_abs_tol ) ? from_quad_tol :
      from_abs_tol;
  // data for adaptive quadrature
  quad_F::adapt_quad_info quad_info;  // information for adaptive quadrature
  quad_F::one_interval top_interval( value->order, value->conserve );
  switch( use_quad )
  {
  case ADAPTIVE2:
  case ADAPTIVE4:  // not yet implemented
    // really short intergvals can cause problems
    if( B - A <= abs_AB * skip_tol/100 )
    {
      //    quad_F::midpoint_rule( F, A, B, params, value );
      break;
    }
    else if( B - A <=  abs_AB * skip_tol )
    {
      quad_F::Gauss_2( F, A, B, params, value );
      break;
    }

    quad_info.set_order( value->order, value->conserve );
    top_interval.initialize_top( F, A, B, params, &quad_info );
    quad_info.set_tolerance( tol );
    quad_F::adapt_quad2( F, A, B, &top_interval, value );
    break;
    //  case ADAPTIVE4:  not yet implemented
    //    params->func_count = 0;
    //    quad_F::adapt_quad4( F, A, B, params, value );
    //    break;
  case ADAPT_HALF:
    // really short intergvals can cause problems
    if( B - A <=  abs_AB * skip_tol )
    {
      quad_F::Gauss_half( F, A, B, params, value );
      break;
    }

    quad_info.set_order( value->order, value->conserve );
    top_interval.initialize_top_half( F, A, B, params, &quad_info );
    quad_info.set_tolerance( tol );
    quad_F::adapt_quad_half( F, A, B, &top_interval, value );
    break;
  case GAUSS2:
    quad_F::Gauss_2( F, A, B, params, value );
    break;
  case GAUSS4:
    quad_F::Gauss_4( F, A, B, params, value );
    break;
  case GAUSS6:
    quad_F::Gauss_6( F, A, B, params, value );
    break;
  case GAUSS10:
    quad_F::Gauss_10( F, A, B, params, value );
    break;
  default:
    FatalError("quad_F::integrate", "bad quadrature method");
  }
  if( quad_info.warning_set )
  {
    Warning("quad_F::integrate", "The integral may not be accurate");
  }
}

// ************* quad_F::midpoint_rule *******************************
// The midpoint rule on the interval (A, B).
void quad_F::midpoint_rule( void (*F)( double x, QuadParamBase *params, coef_vector *Value ), 		double A, double B, QuadParamBase *params, coef_vector *value )
{
  value->set_zero( );
  double M = ( A + B ) / 2;  // midpoint
  F( M, params, value );
  *value *= B - A;
}

// ************* quad_F::Gauss_2 *******************************
// Second-order Gaussian quadrature on the interval (A, B).
void quad_F::Gauss_2( void (*F)( double x, QuadParamBase *params, coef_vector *Value ), 
		double A, double B, QuadParamBase *params, coef_vector *value )
{
  value->set_zero( );
  // really short intergvals can cause problems
  static double from_quad_tol = Global.Value( "abs_quad_tol" )/100;
  static double from_abs_tol = 1000*Global.Value( "abs_tol" );
  static double skip_tol = ( from_quad_tol > from_abs_tol ) ? from_quad_tol : from_abs_tol;
  if( B - A <= abs( B )* skip_tol/100 )
  {
    //    quad_F::midpoint_rule( F, A, B, params, value );
    return;
  }

  // Gauss quadrature points
  const double x1 = 1.0/sqrt( 3.0 );
  double M = ( A + B ) / 2;  // midpoint
  double scale = B - M;
  double abcissae[ 2 ] = { M - scale * x1,
			   M + scale * x1 };
  // weights
  const double wt1 = 0.5;
  double weights[ 2 ] = { wt1, wt1 };

  // to hold the function values
  coef_vector wert( value->order, value->conserve );
  for( int x_count = 0; x_count < 2; ++x_count )
  {
    F( abcissae[ x_count ], params, &wert );
    wert *= weights[ x_count ];
    *value += wert;
  }
  *value *= B - A;
}

// ************* quad_F::Gauss_4 *******************************
// Fourth-order Gaussian quadrature on the interval (A, B).
void quad_F::Gauss_4( void (*F)( double x, QuadParamBase *params, coef_vector *Value ), 
		double A, double B, QuadParamBase *params, coef_vector *value )
{
  value->set_zero( );
  // really short intergvals can cause problems
  static double from_quad_tol = Global.Value( "abs_quad_tol" )/100;
  static double from_abs_tol = 1000*Global.Value( "abs_tol" );
  static double skip_tol = ( from_quad_tol > from_abs_tol ) ? from_quad_tol : from_abs_tol;
  if( B - A <= abs( B )* skip_tol/100 )
  {
    quad_F::midpoint_rule( F, A, B, params, value );
    return;
  }

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
  coef_vector wert( value->order, value->conserve );
  for( int x_count = 0; x_count < 4; ++x_count )
  {
    F( abcissae[ x_count ], params, &wert );
    wert *= weights[ x_count ];
    *value += wert;
  }
  *value *= B - A;
}

// ************* quad_F::Gauss_6 *******************************
// Sixth-order Gaussian quadrature on the interval (A, B).
void quad_F::Gauss_6( void (*F)( double x, QuadParamBase *params, coef_vector *Value ), 
		double A, double B, QuadParamBase *params, coef_vector *value )
{
  value->set_zero( );
  // really short intergvals can cause problems
  static double from_quad_tol = Global.Value( "abs_quad_tol" )/100;
  static double from_abs_tol = 1000*Global.Value( "abs_tol" );
  static double skip_tol = ( from_quad_tol > from_abs_tol ) ? from_quad_tol : from_abs_tol;
  if( B - A <= abs( B )* skip_tol/100 )
  {
    quad_F::midpoint_rule( F, A, B, params, value );
    return;
  }

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
  coef_vector wert( value->order, value->conserve );
  for( int x_count = 0; x_count < 6; ++x_count )
  {
    F( abcissae[ x_count ], params, &wert );
    wert *= weights[ x_count ];
    *value += wert;
  }
  *value *= B - A;
}

// ************* quad_F::Gauss_10 *******************************
// Tenth-order Gaussian quadrature on the interval (A, B).
void quad_F::Gauss_10( void (*F)( double x, QuadParamBase *params, coef_vector *Value ), 
		double A, double B, QuadParamBase *params, coef_vector *value )
{
  value->set_zero( );
  // really short intergvals can cause problems
  static double from_quad_tol = Global.Value( "abs_quad_tol" )/100;
  static double from_abs_tol = 1000*Global.Value( "abs_tol" );
  static double skip_tol = ( from_quad_tol > from_abs_tol ) ? from_quad_tol : from_abs_tol;
  if( B - A <= abs( B )* skip_tol/100 )
  {
    quad_F::midpoint_rule( F, A, B, params, value );
    return;
  }

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
  coef_vector wert( value->order, value->conserve );
  for( int x_count = 0; x_count < 10; ++x_count )
  {
    F( abcissae[ x_count ], params, &wert );
    wert *= weights[ x_count ];
    *value += wert;
  }
  *value *= B - A;
}

// ************* quad_F::Gauss_half *******************************
// 1-point Gaussian quadrature on (A, B) with singularity $\sqrt{1 - x}$
void quad_F::Gauss_half(
    void (*F)( double x, QuadParamBase *params, coef_vector *Value ),
 		double A, double B, QuadParamBase *params, coef_vector *value )
{
  value->set_zero( );
  // The evaluation point is at
  // M = 1 - 0.6*(1 - A)*(1 - beta^5)/(1 - beta^3)
  // with beta = sqrt{(1 - B)/(1 - A) }.
  double beta = sqrt( (1 - B)/(1 - A) );
  double denom = 1 + beta*( 1 + beta );
  double numerator = 1 + beta*( 1 + beta*denom );
  double M = 1 - 0.6*( 1 - A )*numerator/denom;
  F( M, params, value );
  *value *= B - A;
}

// ************* quad_F::adapt_quad2 *******************************
// Adaptive Gaussian quadrature on the interval (A, B).
void quad_F::adapt_quad2( void (*F)( double x, QuadParamBase *params, coef_vector *Value ), 
		    double A, double B, quad_F::one_interval *top_interval,
                    coef_vector *value )
{
  // set up the data for the subintervals
  double m = 0.5*( A + B );

  quad_F::one_interval left_interval( value->order, value->conserve );
  left_interval.initialize_sub( F, A, m, top_interval->params,
     top_interval->quad_info );
  quad_F::one_interval right_interval( value->order, value->conserve );
  right_interval.initialize_sub( F, m, B, top_interval->params,
     top_interval->quad_info );

  // test for trivial subdivision
  if( ( m <= A ) || ( B <= m ) )
  {
    top_interval->quad_info->warning_set = true;
    *value += left_interval.integral;  // update the running integral
    *value += right_interval.integral;  // update the running integral
    return;
  }

  // test for too many subdivisions
  static bool warned = false;
  top_interval->quad_info->interval_count += 2;
  static int max_count = Global.Value( "max_divisions" );
  if( top_interval->quad_info->interval_count > max_count )
  {
    if( !warned )
    {
      Warning("quad_F::adapt_quad2", "Subdivision limit exceeded");
      warned = true;
    }
    top_interval->quad_info->warning_set = true;
    *value += left_interval.integral;  // update the running integral
    *value += right_interval.integral;  // update the running integral
    return;
  }

  // test the accuracy
  if( top_interval->test_OK( left_interval.integral, right_interval.integral, 4 ) )
  {
    *value += left_interval.integral;  // update the running integral
    *value += right_interval.integral;  // update the running integral
  }
  else
  {
    // iterate
    quad_F::adapt_quad2( F, A, m, &left_interval, value );
    quad_F::adapt_quad2( F, m, B, &right_interval, value );
  }
}

// ************* quad_F::adapt_quad_half *******************************
// Adaptive Gaussian quadrature on the interval (A, B)
// with $sqrt{1 - x}$ singularity.
void quad_F::adapt_quad_half( void (*F)( double x, QuadParamBase *params, coef_vector *Value ), 
		    double A, double B, quad_F::one_interval *top_interval,
                    coef_vector *value )
{
  // set up the data for the subintervals
  double m = 0.5*( A + B );

  quad_F::one_interval left_interval( value->order, value->conserve );
  left_interval.initialize_sub_half( F, A, m, top_interval->params,
     top_interval->quad_info );
  quad_F::one_interval right_interval( value->order, value->conserve );
  right_interval.initialize_sub_half( F, m, B, top_interval->params,
     top_interval->quad_info );

  // test for trivial subdivision
  if( ( m <= A ) || ( B <= m ) )
  {
    top_interval->quad_info->warning_set = true;
    *value += left_interval.integral;  // update the running integral
    *value += right_interval.integral;  // update the running integral
    return;
  }

  // test for too many subdivisions
  top_interval->quad_info->interval_count += 2;
  static int max_count = Global.Value( "max_divisions" );
  if( top_interval->quad_info->interval_count > max_count )
  {
    if( !top_interval->quad_info->warning_set )
    {
      Warning("quad_F::adapt_quad_half", "Subdivision limit exceeded");
      top_interval->quad_info->warning_set = true;
    }
    *value += left_interval.integral;  // update the running integral
    *value += right_interval.integral;  // update the running integral
    return;
  }

  // test the accuracy
  if( top_interval->test_OK( left_interval.integral, right_interval.integral, 2 ) )
  {
    *value += left_interval.integral;  // update the running integral
    *value += right_interval.integral;  // update the running integral
  }
  else
  {
    // iterate
    quad_F::adapt_quad_half( F, A, m, &left_interval, value );
    quad_F::adapt_quad_half( F, m, B, &right_interval, value );
  }
}

// ************* quad_F::adapt_quad_info *******************************
// ------------------- quad_F::adapt_quad_info::set_order ------------------
// Sets the Legendre order and the conservation flag
void quad_F::adapt_quad_info::set_order( int Order, Conserve cons )
{
  order = Order;
  conserve = cons;
  quad_tolerance.order = order;
  quad_tolerance.conserve = conserve;
  rough_est.order = order;
  rough_est.conserve = conserve;
}
// ------------------- quad_F::adapt_quad_info::set_tolerance ------------------
// Sets the tolerances for the different Legendre orders
void quad_F::adapt_quad_info::set_tolerance( double tol )
{
  if( ( Norm_1 <= 0.0 ) || ( Norm_E <= 0.0 ) )
  {
    FatalError( "quad_F::adapt_quad_info::set_tolerance",
		"Set the norms of the approximate integrals" );
  }
  double EPS = 2*DBL_EPSILON;  // 2*(machine accuracy)
  double tol_ = ( tol < EPS ) ? EPS : tol;

  if( ( conserve == NUMBER ) || ( conserve == BOTH ) )
  {
    for( int i = 0; i <= order; ++i )
    {
      if( rough_est.weight_1[ i ] < Norm_1*tol_ )
      {
        quad_tolerance.weight_1[ i ] = 1.0;
      }
      else
      {
        quad_tolerance.weight_1[ i ] = Norm_1*tol_/rough_est.weight_1[ i ];
      }
    }
  }

  if( ( conserve == ENERGY ) || ( conserve == BOTH ) )
  {
    for( int i = 0; i <= order; ++i )
    {
      if( rough_est.weight_E[ i ] < Norm_E*tol_ )
      {
        quad_tolerance.weight_E[ i ] = 1.0;
      }
      else
      {
        quad_tolerance.weight_E[ i ] = Norm_E*tol_/rough_est.weight_E[ i ];
      }
    }
  }
}
// ************ one_interval ******************
// ------------------- quad_F::one_interval::one_interval ------------------
// constructor
quad_F::one_interval::one_interval( int Order, Conserve cons )
{
  set_order( Order, cons );
}
// ------------------- quad_F::one_interval::set_order ------------------
// constructor
void quad_F::one_interval::set_order( int Order, Conserve cons )
{
  integral.set_order( Order, cons );
}
// ------------------- quad_F::one_interval::initialize_top ------------------
// Get a rough approximate integral on the initial full interval
void quad_F::one_interval::initialize_top( void (*F)( double x, QuadParamBase *params, coef_vector *Value ), 
				       double A, double B, QuadParamBase *Params,
				       quad_F::adapt_quad_info *Quad_info )
{
  params = Params;
  quad_info = Quad_info;
  // get an initial estimate of the integral
  quad_F::Gauss_2( F, A, B, params, &integral );
  // save a copy as our rough estimate
  quad_info->rough_est.copy( integral );

  // get the max norms
  quad_info->rough_est.max_norm( &quad_info->Norm_1, &quad_info->Norm_E );
  // ensure nonzero comparisons 
  quad_info->rough_est.test_zero( &quad_info->Norm_1, &quad_info->Norm_E, B - A );
}
// ------------------- quad_F::one_interval::initialize_top_half ----------------
// Get a rough approximate integral on the initial full interval 
// for the case of $sqrt{1 - x}$ singularity
void quad_F::one_interval::initialize_top_half(
     void (*F)( double x, QuadParamBase *params, coef_vector *Value ), 
				       double A, double B, QuadParamBase *Params,
				       quad_F::adapt_quad_info *Quad_info )
{
  params = Params;
  quad_info = Quad_info;
  // get an initial estimate of the integral
  quad_F::Gauss_half( F, A, B, params, &integral );
  // save a copy as our rough estimate
  quad_info->rough_est.copy( integral );

  // get the max norms
  quad_info->rough_est.max_norm( &quad_info->Norm_1, &quad_info->Norm_E );
  // ensure nonzero comparisons 
  quad_info->rough_est.test_zero( &quad_info->Norm_1, &quad_info->Norm_E, B - A );
}
// ------------------- quad_F::one_interval::initialize_sub ------------------
// Initialize the integral on a subinterval
void quad_F::one_interval::initialize_sub( void (*F)( double x, QuadParamBase *params, coef_vector *Value ), 
				       double A, double B, QuadParamBase *Params,
				       quad_F::adapt_quad_info *Quad_info )
{
  params = Params;
  quad_info = Quad_info;
  // get an initial estimate of the integral
  quad_F::Gauss_2( F, A, B, params, &integral );
}
// ------------------- quad_F::one_interval::initialize_sub_half ----------------
// Initialize the integral on a subinterval
// with $sqrt{1 - x}$ singularity.
void quad_F::one_interval::initialize_sub_half(
     void (*F)( double x, QuadParamBase *params, coef_vector *Value ), 
				       double A, double B, QuadParamBase *Params,
				       quad_F::adapt_quad_info *Quad_info )
{
  params = Params;
  quad_info = Quad_info;
  // get an initial estimate of the integral
  quad_F::Gauss_half( F, A, B, params, &integral );
}
// ------------------- quad_F::one_interval::test_OK ------------------
// Tests whether two estimates are sufficiently close
bool quad_F::one_interval::test_OK( const coef_vector& left_half,
				    const coef_vector& right_half, int quad_order )
{
  int order = integral.order;
  Conserve cons = integral.conserve;
  coef_vector extrapolate( order, cons ); 
  // Do a Richardson extrapolation
  //  For quad_F::Gauss_4 use
  //  extrapolate = ( 64*(left_half + right_half) - course ) / 63;
  //  For quad_F::Gauss_2 use
  //  extrapolate = ( 16*(left_half + right_half) - course ) / 15;
  double scale = 1.0/( quad_order * quad_order );
  extrapolate.copy( integral );
  extrapolate *= -scale;
  extrapolate += left_half;
  extrapolate += right_half;
  extrapolate *= 1.0/( 1 - scale );

  static double abs_quad_tol = Global.Value( "abs_quad_tol" );

  if( ( cons == NUMBER ) || ( cons == BOTH ) ){
    for( int i = 0; i <= order; ++i )
    {
      if( abs( extrapolate.weight_1[ i ] - left_half.weight_1[ i ] - right_half.weight_1[ i ] ) >
	    quad_info->quad_tolerance.weight_1[ i ]*quad_info->rough_est.weight_1[ i ] &&
	  ( abs( extrapolate.weight_1[ i ] ) > quad_info->Norm_1*abs_quad_tol ) )
      {
	return false;
      }
    }
  }
  if( ( cons == ENERGY ) || ( cons == BOTH ) ){
    for( int i = 0; i <= order; ++i )
    {
      if( abs( extrapolate.weight_E[ i ] - left_half.weight_E[ i ] - right_half.weight_E[ i ] ) >
	  quad_info->quad_tolerance.weight_E[ i ]*quad_info->rough_est.weight_E[ i ] &&
	  ( abs( extrapolate.weight_E[ i ] ) > quad_info->Norm_E*abs_quad_tol ) )
      {
	return false;
      }
    }
  }
  return true;
}

// ************* Legendre polynomials *******************************
// ---------------- math_F::Legendre ------------------
void math_F::Legendre( double mu, coef_vector *value )
{
  // Use the iteration formula to compute the Legendre functions
  //   (n+1) P_{n+1}(mu) - (2n+1) mu P_n(mu) + n P_{n-1}(mu) = 0
  // See Courant and Hilbert, Methods of Mathematical Physics,
  // vol. 1, p. 86.

  // set the zero-order coefficient
  if( ( value->conserve == NUMBER ) || ( value->conserve == BOTH ) ){
    value->weight_1[ 0 ] = 1.0;
  }
  if( ( value->conserve == ENERGY ) || ( value->conserve == BOTH ) ){
    value->weight_E[ 0 ] = 1.0;
  }
  if( value->order == 0 ){
    return;
  }

  // set the first-order coefficient
  if( ( value->conserve == NUMBER ) || ( value->conserve == BOTH ) ){
    value->weight_1[ 1 ] = mu;
  }
  if( ( value->conserve == ENERGY ) || ( value->conserve == BOTH ) ){
    value->weight_E[ 1 ] = mu;
  }
  if( value->order == 1 ){
    return;
  }

  double P_prev = 1.0;  // P_0(mu)
  double P_this = mu;   // P_1(mu)
  double P_next;
  for(int ell = 1; ell < value->order; ++ell)
  {
    P_next = ((2*ell + 1)*mu*P_this - ell*P_prev)/(ell + 1.0);
    if( ( value->conserve == NUMBER ) || ( value->conserve == BOTH ) ){
      value->weight_1[ ell + 1 ] = P_next;
    }
    if( ( value->conserve == ENERGY ) || ( value->conserve == BOTH ) ){
      value->weight_E[ ell + 1 ] = P_next;
    }
    // shift the values
    P_prev = P_this;
    P_this = P_next;
  }
}

// ------------------ math_F::Gamma_up --------------------------------
// Increments the incomplete Gamma function \int_A^\infty t^{kappa-1} e^{-t}
double  math_F::Gamma_up( double kappa, double A, double Gamma_kappa )
{
  // Returns Gamma( kappa + 1, A ); this is a stable iteration
  double Gamma = pow( A, kappa )*exp( -A ) + kappa*Gamma_kappa;
  return Gamma;
}

// ------------------  math_F::gamma_down --------------------------------
// Decrements the incomplete Gamma function \int_0^A t^{kappa-1} e^{-t}
double  math_F::gamma_down( double kappa, double A, double gamma_kappa )
{
  // Returns gamma( kappa - 1, A ); this is a stable iteration
  if( kappa <= 1.0 )
  {
    FatalError( "gamma_down", pastenum( "improper value of kappa: ", kappa) );
  }
  double gamma = ( pow( A, kappa-1 )*exp( -A ) + gamma_kappa )/(kappa - 1);
  return gamma;
}

// ------------------ math_F::zeroin --------------------------------
// Find a root of func(x, params) = target between BB and CC.
// Return the root with accuracy of tol + 4*EPS*|root|
double math_F::zeroin(double (*func)(double, void*),
              double target,
              const dd_entry& BB,
              const dd_entry& CC,
	      void *params,
              double tol)

/* This routine is a translation from the Brent zeroin
 * from netlib.  The original is available by e-mail:
 *   e-mail: netlib@ornl.gov
 *   subject: send zeroin from go
*/
{
  //Stash the original 1d_links
  dd_entry B = BB;
  dd_entry C = CC;

  const int MAX_ITER = 100;  // the maximum number of tries

  double rel_tol;     // bound on relative error
  double diff_cb;   // (c - b)/2
  double new_width; // |diff_cb|

  double num;    // fraction for the secant rule
  double denom;

  int num_poor = 0;  // how many poor improvements

  // use a reasonable tolerance
  static double EPS = 4.0 * DBL_EPSILON;
  double use_tol = (tol < EPS) ? EPS : tol;

  // the current interval
  double old_width = abs(C.x - B.x);

  // shift by the target
  B.y -= target;
  C.y -= target;

  // worse error
  double worse = (abs(B.y) < abs(C.y)) ? abs(C.y) : abs(B.y);

  // the oldest estimate
  dd_entry A = C;

  for(int count = 0; count < MAX_ITER; ++count)
  {
    // make B be the best estimate so far
    if(abs(C.y) < abs(B.y))
    {
      A = B;
      B = C;
      C = A;
    }

    // how close are we?
    diff_cb = 0.5*(C.x - B.x);
    new_width = abs(diff_cb);
    rel_tol = use_tol*abs(B.x) + EPS;
    if(new_width <= rel_tol)
    {
      if(B.y*C.y > 0.0)
      {
        SevereError("math_F::zeroin", "No root found" );
      }
      else if(abs(B.y) > worse)
      {
        SevereError("math_F::zeroin", "This looks like a pole.");
      }
      else
      {
	if( ( B.x < BB.x ) || ( B.x > CC.x ) )
	{
	  Warning( "math_F::zeroin", "root outside original interval" );
	}
        return B.x;
      }
    }

    // set up the secant iteration
    num = (B.x - A.x)*B.y;
    denom = A.y - B.y;

    // arrange so that num >= 0
    if(num < 0.0)
    {
      num *= -1;
      denom *= -1;
    }

    // save the best so far
    A = B;

    // have we had too many poor ones?
    ++num_poor;
    if((num_poor >= 4) && (8.0*new_width < old_width))
    {
      num_poor = 0;
      old_width = new_width;
    }

    // which type of iteration?
    if(num_poor >= 4)
    {
      B.x = 0.5*(C.x + B.x);
    }
    else if(num < abs(denom)*rel_tol)
    {
      // too small a change
      B.x += (diff_cb > 0.0) ? rel_tol : -rel_tol;
    }
    else if(num < denom*diff_cb)
    {
      // secant rule if x between B and (C + B)/2
      B.x += num/denom;
    }
    else
    {
      // bisection
      B.x = 0.5*(C.x + B.x);
    }

    // the new function value
    B.y = func(B.x, params) - target;

    // did we hit it?
    if(B.y == 0.0)
    {
      return B.x;
    }

    // which old point do we keep?
    if(B.y*C.y > 0.0)
    {
      C = A;
    }
    //    cout << B.x << "  " << B.y << endl;
  }

  // if we got here, there were too many iterations
  Warning("math_F::zeroin","Too many iterations in zeroin");
  return B.x;  //However, we MUST return something or exit()?
}

// ------------------ math_F::parabola_bottom --------------------------------
// Fit func(x, params) by a parabola at AA, CC, and their midpoint BB.
// Returns the minimum MM of the parabola.
// If MM < BB, set CC = BB, otherwise set AA = BB.
double math_F::parabola_bottom(double (*func)(double, void*),
              dd_entry *AA,
              dd_entry *CC,
	      void *params )
{
  static double etol = Global.Value( "E_tol" ); // to check a flat bottom

  double Ein_mid = 0.5*( AA->x + CC->x );
  dd_entry BB( Ein_mid, func( Ein_mid, params ) );
  // fit with c + b*x + a*x*x with x = Ein - Ein_mid
  double h = 0.5*( CC->x - AA->x );
  if( h <= 0.0 )
  {
    FatalError( "math_F::parabola_bottom", 
		"incident energies out of order" );
  }
  double a = ( AA->y - BB.y ) - ( BB.y - CC->y );
  double b = CC->y - AA->y;
  if( a <= 0.0 )
  {
    // check for a flat bottom
    if( -a < etol*abs( BB.y ) )
    {
      if( b > 0.0 )  // small positive slope
      {
	*CC = *AA;  // copy
        return AA->x;
      }
      else if( b < 0.0 )  // small negative slope
      {
	*AA = *CC;  // copy
        return CC->x;
      }
      else  // really flat
      {
	*AA = BB;  // copy
	*CC = BB;  // copy
        return BB.x;
      }
    }
    else
    {
      FatalError( "math_F::parabola_bottom", 
		"parabola not convex up" );
    }
  }
  a /= ( 2*h*h );
  b /= ( 2*h );
  // bottom of the parabola at
  double bottom = -b/(2*a);
  // reset the bounds
  if( b < 0.0 )
  {
    *CC = BB;  // do a copy
  }
  else
  {
    *AA = BB;  // do a copy
  }
  return BB.x + bottom;
 }
