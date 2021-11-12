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
        Msg::Warning( "quadratic", "reduces to 0 = 0" );
	*alpha_1 = 0.0;
	*alpha_2 = 1.0;
      }
      else
      {
        num_roots = 0;
        Msg::Warning( "quadratic", "reduces to 0 = 1" );
      }
    }
    else
    {
      num_roots = 1;
      //      Msg::Warning( "quadratic", "only linear terms" );
      *alpha_1 = -C/B;
    }
    return num_roots;
  }

  double root_1;
  double root_2;
  double discriminant = B*B - 4*A*C;
  static double abs_tol = Global.Value( "tight_tol" );

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


// ************* Legendre polynomials *******************************
// ---------------- math_F::Legendre ------------------
void math_F::Legendre( double mu, Coef::coef_vector *value )
{
  // Use the iteration formula to compute the Legendre functions
  //   (n+1) P_{n+1}(mu) - (2n+1) mu P_n(mu) + n P_{n-1}(mu) = 0
  // See Courant and Hilbert, Methods of Mathematical Physics,
  // vol. 1, p. 86.

  // set the zero-order coefficient
  if( ( value->conserve == Coef::NUMBER ) || ( value->conserve == Coef::BOTH ) ){
    value->weight_1[ 0 ] = 1.0;
  }
  if( ( value->conserve == Coef::ENERGY ) || ( value->conserve == Coef::BOTH ) ){
    value->weight_E[ 0 ] = 1.0;
  }
  if( value->order == 0 ){
    return;
  }

  // set the first-order coefficient
  if( ( value->conserve == Coef::NUMBER ) || ( value->conserve == Coef::BOTH ) ){
    value->weight_1[ 1 ] = mu;
  }
  if( ( value->conserve == Coef::ENERGY ) || ( value->conserve == Coef::BOTH ) ){
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
    if( ( value->conserve == Coef::NUMBER ) || ( value->conserve == Coef::BOTH ) ){
      value->weight_1[ ell + 1 ] = P_next;
    }
    if( ( value->conserve == Coef::ENERGY ) || ( value->conserve == Coef::BOTH ) ){
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
    Msg::FatalError( "gamma_down",
		     Msg::pastenum( "improper value of kappa: ", kappa) );
  }
  double gamma = ( pow( A, kappa-1 )*exp( -A ) + gamma_kappa )/(kappa - 1);
  return gamma;
}

// ------------------ math_F::zeroin --------------------------------
// Find a root of func(x, params) = target between BB and CC.
// Return the root with accuracy of tol + 4*EPS*|root|
double math_F::zeroin(double (*func)(double, void*),
              double target,
              const Ddvec::dd_entry& BB,
              const Ddvec::dd_entry& CC,
	      void *params,
              double tol)

/* This routine is a translation from the Brent zeroin
 * from netlib.  The original is available by e-mail:
 *   e-mail: netlib@ornl.gov
 *   subject: send zeroin from go
*/
{
  if( BB.x >= CC.x )
  {
    Msg::SevereError( "math_F::zeroin", "data out of order" );
  }
  
  //Stash the original 1d_links
  Ddvec::dd_entry B = BB;
  Ddvec::dd_entry C = CC;

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
  double old_width = std::abs(C.x - B.x);

  // shift by the target
  B.y -= target;
  C.y -= target;

  // worse error
  double worse = (std::abs(B.y) < std::abs(C.y)) ? std::abs(C.y) : std::abs(B.y);

  // the oldest estimate
  Ddvec::dd_entry A = C;

  for(int count = 0; count < MAX_ITER; ++count)
  {
    // make B be the best estimate so far
    if(std::abs(C.y) < std::abs(B.y))
    {
      A = B;
      B = C;
      C = A;
    }

    // how close are we?
    diff_cb = 0.5*(C.x - B.x);
    new_width = std::abs(diff_cb);
    rel_tol = use_tol*std::abs(B.x) + EPS;
    if(new_width <= rel_tol)
    {
      if(B.y*C.y > 0.0)
      {
        Msg::SevereError("math_F::zeroin", "No root found" );
      }
      else if(std::abs(B.y) > worse)
      {
        Msg::SevereError("math_F::zeroin", "This looks like a pole.");
      }
      else
      {
	if( ( B.x < BB.x ) || ( B.x > CC.x ) )
	{
	  Msg::Warning( "math_F::zeroin", "root outside original interval" );
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
    else if(num < std::abs(denom)*rel_tol)
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
  Msg::Warning("math_F::zeroin","Too many iterations in zeroin");
  return B.x;  //However, we MUST return something or exit()?
}

// ------------------ math_F::parabola_bottom --------------------------------
// Fit func(x, params) by a parabola at AA, CC, and their midpoint BB.
// Returns the minimum MM of the parabola.
double math_F::parabola_bottom(double (*func)(double, void*),
             const Ddvec::dd_entry &AA,
             const Ddvec::dd_entry &CC,
	      void *params )
{
  double Ein_mid = 0.5*( AA.x + CC.x );
  Ddvec::dd_entry BB( Ein_mid, func( Ein_mid, params ) );
  // fit with c + b*x + a*x*x with x = Ein - Ein_mid
  double h = 0.5*( CC.x - AA.x );
  if( h <= 0.0 )
  {
    Msg::FatalError( "math_F::parabola_bottom", 
		"incident energies out of order" );
  }
  double a = ( AA.y - BB.y ) - ( BB.y - CC.y );
  double b = CC.y - AA.y;
  if( a <= 0.0 )
  {
    return Ein_mid;
  }
  a /= ( 2*h*h );
  b /= ( 2*h );
  // bottom of the parabola at
  double bottom = -b/(2*a);

  return BB.x + bottom;
 }
