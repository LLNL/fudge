/*
 * ******** merced: calculate the transfer matrix ********
 * $Revision: 601 $
 * $Date: 2017-12-08 $
 * $Author: hedstrom $
 * $Id: adapt_quad.cpp 601  2017-12-08Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */

// Implements the classes for adaptive quadrature

#include "adapt_quad.hpp"

// ************ quad_F::one_interval ******************
// ------------------- quad_F::one_interval::one_interval ------------------
// constructor
quad_F::one_interval::one_interval( int Order, Coef::Conserve cons )
{
  set_order( Order, cons );
}
// ------------------- quad_F::one_interval::set_order ------------------
// constructor
void quad_F::one_interval::set_order( int Order, Coef::Conserve cons )
{
  integral.set_order( Order, cons );
}
// ------------------- quad_F::one_interval::test_OK ------------------
// Tests whether two estimates are sufficiently close
bool quad_F::one_interval::test_OK( const Coef::coef_vector& left_half,
				    const Coef::coef_vector& right_half ) const
{
  int order = integral.order;
  Coef::Conserve cons = integral.conserve;
  Coef::coef_vector extrapolate( order, cons ); 

  // Do a Richardson extrapolation
  double scale = 1.0/ quad_info->Richardson;
  extrapolate.copy( integral );
  extrapolate *= -scale;
  extrapolate += left_half;
  extrapolate += right_half;
  extrapolate *= 1.0/( 1 - scale );

  static double quad_tol = Global.Value( "quad_tol" );
  static double quad_tolfloor = Global.Value( "quad_tol_floor" );
  double Norm_1;
  double Norm_E;
  quad_info->rough_est.max_norm( &Norm_1, &Norm_E );
  
  if( ( cons == Coef::NUMBER ) || ( cons == Coef::BOTH ) ){
    for( int i = 0; i <= order; ++i )
    {
      if( ( std::abs( extrapolate.weight_1[ i ] - left_half.weight_1[ i ] - right_half.weight_1[ i ] ) >
	    quad_tol*quad_info->quad_weight[ i ]*quad_info->rough_est.weight_1[ i ] ) &&
	  ( std::abs( extrapolate.weight_1[ i ] ) > quad_tolfloor * Norm_1 ) )
      {
	return false;
      }
    }
  }
  if( ( cons == Coef::ENERGY ) || ( cons == Coef::BOTH ) ){
    for( int i = 0; i <= order; ++i )
    {
      if( ( std::abs( extrapolate.weight_E[ i ] - left_half.weight_E[ i ] - right_half.weight_E[ i ] ) >
	  quad_tol*quad_info->quad_weight[ i ]*quad_info->rough_est.weight_E[ i ] ) &&
	  ( std::abs( extrapolate.weight_E[ i ] ) > quad_tolfloor * Norm_E ) )
      {
	return false;
      }
    }
  }
  return true;
}

// ************ quad_F::quad_list ******************
// ------------------- quad_F::quad_list::set_up ------------------
// Evaluate the integral
void quad_F::quad_list::set_up( bool (*F)( double x, Qparam::QuadParamBase *params,
					       Coef::coef_vector *Value ),
                    double A, double B, Qparam::QuadParamBase *Params )
{
  // save the information we need
  FF = F;
  AA = A;
  BB = B;
  params = Params;
  
}
// ------------------- quad_F::quad_list::adapt_quad ------------------
// Evaluate the integral
bool quad_F::quad_list::adapt_quad( double tol, Coef::coef_vector *Value )
{
  quad_info.set_order( Value->order, Value->conserve );
  
  // initialize
  bool is_OK = initialize( );
  if( !is_OK )
  {
    return false;
  }
  quad_info.depth = 0;
  
  // set up the rough estimate
  is_OK = get_rough( tol, Value );
  if( !is_OK )
  {
    return false;
  }

  // process the subintervals
  while( !empty( ) )
  {
    is_OK = one_level( Value );
  }
  if( quad_info.warning_set )
  {
    Msg::DebugInfo("quad_list::adapt_quad",
		   "The integral may not be accurate");
    return false;
  }

  return is_OK;
}
// ------------------- quad_F::quad_list::initialize ------------------
// set up the quadrature
bool quad_F::quad_list::initialize( )
{
  // start with a new list
  if( size( ) > 0 )
  {
    erase( begin( ), end( ) );
  }

  // set up the weights for the error tolerance
  quad_info.set_weights( );
  
  // to fill in the first link
  quad_F::quad_list::iterator top_interval = insert( end( ), one_interval( ) );
  // allocate space for the vectors after insertion into the list
  top_interval->set_order( quad_info.order, quad_info.conserve );
  top_interval->a = AA;
  top_interval->b = BB;
  top_interval->params = params;
  top_interval->quad_info = &quad_info;
  
  bool is_OK = Quad_Method( FF, AA, BB, params, &(top_interval->integral) );
  return is_OK;
}
// ------------------- quad_F::quad_list::get_rough ------------------
// set up the rough estimate for comparison
bool quad_F::quad_list::get_rough( double tol, Coef::coef_vector *Value )
{
  double m = 0.5 * ( AA + BB );
  Coef::coef_vector left_interval( quad_info.order, quad_info.conserve );
  Coef::coef_vector right_interval( quad_info.order, quad_info.conserve );

  bool is_OKL = Quad_Method( FF, AA, m, params, &left_interval );
  bool is_OKR = Quad_Method( FF, m, BB, params, &right_interval );
  quad_F::quad_list::iterator top_interval = begin( );
  if( !is_OKL || !is_OKR )
  {
    *Value += top_interval->integral;
    return false;
  }

  ++quad_info.depth;  // add a layer
  
  int order = quad_info.order;
  if( ( quad_info.conserve == Coef::NUMBER ) || ( quad_info.conserve == Coef::BOTH ) )
  {
    for( int i = 0; i <= order; ++i )
    {
      quad_info.rough_est.weight_1[ i ] = 0.5 * std::abs( ( top_interval->integral.weight_1[ i ] +
	      left_interval.weight_1[ i ] + right_interval.weight_1[ i ] ) );
    }
  }
  if( ( quad_info.conserve == Coef::ENERGY ) || ( quad_info.conserve == Coef::BOTH ) )
  {
    for( int i = 0; i <= order; ++i )
    {
      quad_info.rough_est.weight_E[ i ] = 0.5 * std::abs( ( top_interval->integral.weight_E[ i ] +
	     left_interval.weight_E[ i ] + right_interval.weight_E[ i ] ) );
    }
  }
  
  // ensure that rough_est is not zero
  quad_info.rough_est.test_zero( );

  // This may be good enough
  test_link( top_interval, left_interval, right_interval, Value );
  return true;
}
// ------------------- quad_F::quad_list::test_link ------------------
// test a link for accuracy
void quad_F::quad_list::test_link( quad_F::quad_list::iterator link,
      const Coef::coef_vector &left_half, const Coef::coef_vector &right_half,
      Coef::coef_vector *Value )
{
  // Set up for possible abort
  static const int max_divisions = Global.Value( "max_divisions" );
  
  // short intervals cause trouble with floating-point arithmetic
  static double too_short = Global.Value( "short_interval" );
  double abs_ab = std::abs( link->a ) + std::abs( link->b );
  
  bool link_OK = link->test_OK( left_half, right_half );
  if( link_OK || ( link->b - link->a < too_short * abs_ab ) )
  {
    *Value += left_half;  // update the global integral
    *Value += right_half;
    erase( link );  // we are done with this interval
  }
  else
  {
    // split this interval
    ++quad_info.subdivisions;
    if( quad_info.subdivisions >= max_divisions )
    {
      quad_info.warning_set = true;
      *Value += left_half;  // update the global integral
      *Value += right_half;
      erase( link );  // we are done with this interval
    }
    else
    {
      double m = 0.5 * ( link->a + link->b );
      quad_list::iterator left_ptr = insert( link, one_interval( ) );
      // Allocate space for integral after inserting the link.
      left_ptr->set_order( quad_info.order, quad_info.conserve );
      left_ptr->a = link->a;
      left_ptr->b = m;
      left_ptr->params = params;
      left_ptr->quad_info = &quad_info;
      left_ptr->integral.copy( left_half );

      // update the right half-interval
      link->a = m;
      link->integral.copy( right_half );
    }
  }
}
// ------------------- quad_F::quad_list::one_level ------------------
// Process one level of subdivision
bool quad_F::quad_list::one_level( Coef::coef_vector *Value )
{
  // for subdividing an interval
  // We do one level of subdivision at a time, in case something goes wrong.
  Coef::coef_vector left_interval( quad_info.order, quad_info.conserve );
  Coef::coef_vector right_interval( quad_info.order, quad_info.conserve );

  static int flag_depth = Global.Value( "flag_adapt_quad_depth" );
  ++quad_info.depth;  // add a layer
  bool flagged = false;
  if( quad_info.depth == flag_depth )
  {
    if( !flagged )
    {
      flagged = true;
      Msg::DebugInfo( "quad_F::quad_list::one_level",
	       "adapt_quad at flagged level" );
    }
  }

  // go through the current list
  quad_list::iterator link = begin( );
  quad_list::iterator next_link = link;
  ++next_link;
  for( ; link != end( ); link = next_link, ++next_link )
  {
    double m = 0.5 * ( link->a + link->b );

    bool is_OKL = Quad_Method( FF, link->a, m, params, &left_interval );
    bool is_OKR = Quad_Method( FF, m, link->b, params, &right_interval );
    if( !is_OKL || !is_OKR )
    {
      *Value += link->integral;
      erase( link );  // we are done with this interval
      return false;
    }
 
    // This may be good enough
    test_link( link, left_interval, right_interval, Value );
  }
  return true;
}

// ************* quad_F::integrate *******************************
// Integrates F over the interval (A, B).
bool quad_F::integrate( bool (*F)( double x, Qparam::QuadParamBase *params, Coef::coef_vector *Value ), 
		  Qmeth::Quadrature_Rule use_quad, double A, double B,
		  Qparam::QuadParamBase *params, double tol, Coef::coef_vector *value )
{
  value->set_zero( );
  params->func_count = 0;
  // We may have a trivial integral
  if( A >= B ){
    return false;
  }

  // really short intergvals can cause problems
  // we may integrate over mu, so B may be zero
  static double skip_tol = 10*Global.Value( "tight_tol" );
  double abs_AB = 0.5*( std::abs(A) + std::abs( B ) );
  if( B - A <= abs_AB * skip_tol )
  {
    Qmeth::midpoint_rule( F, A, B, params, value );
    return false;
  }

  bool is_OK = true;
  // which integration option
  if( use_quad.adaptive )
  {
    if( use_quad.quad_method == Qmeth::GAUSS1 )
    {
      quad_F::adapt_Gauss_1 this_method;
      this_method.quad_info.Richardson = 4.0;
      this_method.set_up( F, A, B, params );
      is_OK = this_method.adapt_quad( tol, value );
    }
    else if( use_quad.quad_method == Qmeth::GAUSS2 )
    {
      quad_F::adapt_Gauss_2 this_method;
      this_method.quad_info.Richardson = 16.0;
      this_method.set_up( F, A, B, params );
      is_OK = this_method.adapt_quad( tol, value );
    }
    else if( use_quad.quad_method == Qmeth::GAUSS3 )
    {
      quad_F::adapt_Gauss_3 this_method;
      this_method.quad_info.Richardson = 36.0;
      this_method.set_up( F, A, B, params );
      is_OK = this_method.adapt_quad( tol, value );
    }
    else if( use_quad.quad_method == Qmeth::GAUSS4 )
    {
      quad_F::adapt_Gauss_4 this_method;
      this_method.quad_info.Richardson = 64.0;
      this_method.set_up( F, A, B, params );
      is_OK = this_method.adapt_quad( tol, value );
    }
    else if( use_quad.quad_method == Qmeth::GAUSS6 )
    {
      Msg::FatalError( "quad_F::integrate",
		  "adaptive Gauss6 quadrature not implemented" );
    }
    else if( use_quad.quad_method == Qmeth::GAUSS10 )
    {
      Msg::FatalError( "quad_F::integrate",
		  "adaptive Gauss10 quadrature not implemented" );
    }
    else if( use_quad.quad_method == Qmeth::WEIGHT_L1 )
    {
      quad_F::adapt_Gauss_wt_L1 this_method;
      this_method.quad_info.Richardson = 4.0;
      this_method.set_up( F, A, B, params );
      is_OK = this_method.adapt_quad( tol, value );
    }
    else
    {
      Msg::FatalError( "quad_F::integrate",
		  "adaptive quadrature method not implemented" );
    }
  }
  else
  {
  switch( use_quad.quad_method )
  {
    case Qmeth::GAUSS1:
      is_OK = Qmeth::midpoint_rule( F, A, B, params, value );
      break;
    case Qmeth::GAUSS2:
      is_OK = Qmeth::Gauss_2( F, A, B, params, value );
      break;
    case Qmeth::GAUSS3:
      is_OK = Qmeth::Gauss_3( F, A, B, params, value );
      break;
    case Qmeth::GAUSS4:
      is_OK = Qmeth::Gauss_4( F, A, B, params, value );
      break;
    case Qmeth::GAUSS6:
      is_OK = Qmeth::Gauss_6( F, A, B, params, value );
      break;
    case Qmeth::GAUSS10:
      is_OK = Qmeth::Gauss_10( F, A, B, params, value );
      break;
  case Qmeth::WEIGHT_L1:
      is_OK = Qmeth::weight_L1( F, A, B, params, value );
      break;
    }
  }

  return is_OK;
}
