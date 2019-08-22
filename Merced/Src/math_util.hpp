/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: math_util.hpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */

// header for my math routines
#ifndef MATH_UTIL
#define MATH_UTIL

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <list>

#include "param_base.hpp"
#include "coef_vector.hpp"
#include "quadrature.hpp"

using namespace std;

namespace quad_F
{
  //! We use Walter Gander's adaptive quadrature routine.
  //! See W. Gander and W. Gautschi, "Adaptive quadrature--revisited",
  //! BIT 40 (2000), 84-101.
  //! The outines are modified to integrate coef_vector-valued functions and 
  //! to use Gauss_2 in place of Simpson's rule.

  // ************* integrate *******************************
  //! Integrates F over the interval $(A, B)$
  //! \param F the vector function to integrate: $\int_A^B F$
  //! \param use_quad which quadrature method to use
  //! \param A lower limit of integration
  //! \param B upper limit of integration
  //! \param params parameters for the function, F
  //! \param tol relative tolerance for adaptive quadrature
  //! \param value vector of approximate integrals: $\int_A^B F$
  void integrate( void (*F)( double x, QuadParamBase *params, coef_vector *Value ), 
		  Quadrature_Method use_quad, double A, double B,
                  QuadParamBase *params,
		  double tol, coef_vector *value );

  // ************* midpoint_rule *******************************
  //! The midpoint rule on the interval $(A, B)$.  used on really short intervals
  //! \param F the vector function to integrate: $\int_A^B F$
  //! \param A lower limit of integration
  //! \param B upper limit of integration
  //! \param params parameters for the function, F
  //! \param value vector of approximate integrals: $\int_A^B F$
  void midpoint_rule( void (*F)( double x, QuadParamBase *params, coef_vector *Value ), 
    double A, double B, QuadParamBase *params, coef_vector *value );


  // ************* Gauss_2 *******************************
  //! Second-order Gaussian quadrature on the interval $(A, B)$.
  //! This is used because it is exact for polynomials of degree 3.
  //! \param F the vector function to integrate: $\int_A^B F$
  //! \param A lower limit of integration
  //! \param B upper limit of integration
  //! \param params parameters for the function, F
  //! \param value vector of approximate integrals: $\int_A^B F$
  void Gauss_2( void (*F)( double x, QuadParamBase *params, coef_vector *Value ), 
		double A, double B, QuadParamBase *params, coef_vector *value );

  // ************* Gauss_4 *******************************
  //! Fourth-order Gaussian quadrature on the interval $(A, B)$.
  //! This is used because it is exact for polynomials of degree 7.
  //! \param F the vector function to integrate: $\int_A^B F$
  //! \param A lower limit of integration
  //! \param B upper limit of integration
  //! \param params parameters for the function, F
  //! \param value vector of approximate integrals: $\int_A^B F$
  void Gauss_4( void (*F)( double x, QuadParamBase *params, coef_vector *Value ), 
		double A, double B, QuadParamBase *params, coef_vector *value );

  // ************* Gauss_6 *******************************
  //! Sixth-order Gaussian quadrature on the interval $(A, B)$.
  //! This is used because it is exact for polynomials of degree 11.
  //! \param F the vector function to integrate: $\int_A^B F$
  //! \param A lower limit of integration
  //! \param B upper limit of integration
  //! \param params parameters for the function, F
  //! \param value vector of approximate integrals: $\int_A^B F$
  void Gauss_6( void (*F)( double x, QuadParamBase *params, coef_vector *Value ), 
		double A, double B, QuadParamBase *params, coef_vector *value );

  // ************* Gauss_10 *******************************
  //! Tenth-order Gaussian quadrature on the interval $(A, B)$.
  //! This is used because it is exact for polynomials of degree 19.
  //! \param F the vector function to integrate: $\int_A^B F$
  //! \param A lower limit of integration
  //! \param B upper limit of integration
  //! \param params parameters for the function, F
  //! \param value vector of approximate integrals: $\int_A^B F$
  void Gauss_10( void (*F)( double x, QuadParamBase *params, coef_vector *Value ), 
		double A, double B, QuadParamBase *params, coef_vector *value );

  // ************* Gauss_half *******************************
  //! Second-order Gaussian quadrature on the interval $(A, B)$
  //! with $sqrt{1 - x}$ singularity.
  //! This is used because it is exact for polynomials of degree 3.
  //! \param F the vector function to integrate: $\int_A^B F$
  //! \param A lower limit of integration
  //! \param B upper limit of integration
  //! \param params parameters for the function, F
  //! \param value vector of approximate integrals: $\int_A^B F$
  void Gauss_half( void (*F)( double x, QuadParamBase *params, coef_vector *Value ), 
		double A, double B, QuadParamBase *params, coef_vector *value );

  // ************* adapt_quad_info *******************************
  //! Class holds information for adaptive quadrature
  class adapt_quad_info
  {
  public:
    coef_vector quad_tolerance;  // tolerances for the quadrature
    coef_vector rough_est;  // a rough approximation to the full integral
    double Norm_1;   // max norm of rough_est.weight_1
    double Norm_E;   // max norm of rough_est.weight_E
    int order;  // the Legendre order of the integral
    int interval_count;  // the number of intervals used
    Conserve conserve;  // do we conserve energy or particle number?
    bool warning_set;  // warn of inaccurate adaptive quadrature

    //! Constructor
    adapt_quad_info( ): Norm_1( 0.0 ), Norm_E( 0.0 ), order( -1 ),
			interval_count( 1 ), conserve( NOT_SET ),
                        warning_set( false ) {}

    //! Destructor
    ~adapt_quad_info( ) {}

    //! Sets the Legendre order and the conservation flag
    //! \param Order the Legendre order of the integrand
    //! \param cons conserve energy or particle number
    void set_order( int Order, Conserve cons );

    //! Sets the tolerances for the different Legendre orders
    //! \param tol the relative error tolerance for adaptive quadrature
    void set_tolerance( double tol );
  };

  // ************* one_interval *******************************
  //! Information for one interval in the quad_F::adapt_quad2 routine
  class one_interval
  {
  private:

  public:
    coef_vector integral;  // the approximate integral in this subinterval
    QuadParamBase *params;  // the function parameters for this integration
    quad_F::adapt_quad_info *quad_info;  // information for adaptive quadrature

    //! Constructor
    one_interval( );

    //! Constructor
    //! \param Order the Legendre order of the integrand
    //! \param cons do we conserve energy or particle number?
    one_interval( int Order, Conserve cons );

    //! Destructor
    inline ~one_interval( ) {}

    //! Allocate space for the vectors
    //! \param Order the Legendre order of the integrand
    //! \param cons do we conserve energy or particle number?
    void set_order( int Order, Conserve cons );

    //! Get a rough approximate integral on the initial full interval
    //! \param F the vector function to integrate: $\int_A^B F$
    //! \param A lower limit of integration
    //! \param B upper limit of integration
    //! \param Params the function parameters for F
    //! \param Quad_info tolerance information for the adaptive quadrature
    void initialize_top( void (*F)( double x, QuadParamBase *params, coef_vector *Value ), 
		    double A, double B, QuadParamBase *Params,
		    quad_F::adapt_quad_info *Quad_info );

    //! Initialize data for a subinterval
    //! \param F the vector function to integrate: $\int_A^B F$
    //! \param A lower limit of integration
    //! \param B upper limit of integration
    //! \param Params the function parameters for F
    //! \param Quad_info tolerance information for the adaptive quadrature
    void initialize_sub( void (*F)( double x, QuadParamBase *params, coef_vector *Value ), 
		    double A, double B, QuadParamBase *Params,
		    quad_F::adapt_quad_info *Quad_info );

    //! Initialize the integral on a subinterval with $sqrt{1 - x}$ singularity.
    //! \param F the vector function to integrate: $\int_A^B F$
    //! \param A lower limit of integration
    //! \param B upper limit of integration
    //! \param Params the function parameters for F
    //! \param Quad_info tolerance information for the adaptive quadrature
    void initialize_top_half(
             void (*F)( double x, QuadParamBase *params, coef_vector *Value ), 
		    double A, double B, QuadParamBase *Params,
		    quad_F::adapt_quad_info *Quad_info );

    //! Initialize data for a subinterval with $sqrt{1 - x}$ singularity.
    //! \param F the vector function to integrate: $\int_A^B F$
    //! \param A lower limit of integration
    //! \param B upper limit of integration
    //! \param Params the function parameters for F
    //! \param Quad_info tolerance information for the adaptive quadrature
    void initialize_sub_half(
               void (*F)( double x, QuadParamBase *params, coef_vector *Value ), 
		    double A, double B, QuadParamBase *Params,
		    quad_F::adapt_quad_info *Quad_info );

    //! Tests whether two estimates are sufficiently close
    //! \param left_half approximate integral over left half subinterval
    //! \param right_half approximate integral over right half subinterval
    //! \param quad_order the order of accuracy of the quadrature method
    bool test_OK( const coef_vector& left_half, const coef_vector& right_half,
          int quad_order );

  };

  // ************* adapt_quad2 *******************************
  //! Adaptive 2nd-order Gaussian quadrature on the interval $(A, B)$.
  //! \param F the vector function to integrate: $\int_A^B F$
  //! \param A lower limit of integration
  //! \param B upper limit of integration
  //! \param top_interval information for adaptive quadrature
  //! \param value vector of approximate integrals: $\int_A^B F$
  void adapt_quad2( void (*F)( double x, QuadParamBase *params, coef_vector *Value ), 
		    double A, double B, quad_F::one_interval *top_interval,
                    coef_vector *value );

  // ************* adapt_quad4 *******************************
  //! Adaptive 4th-order Gaussian quadrature on the interval $(A, B)$.
  //! \param F the vector function to integrate: $\int_A^B F$
  //! \param A lower limit of integration
  //! \param B upper limit of integration
  //! \paramtop_interval information for adaptive quadrature
  //! \param value vector of approximate integrals: $\int_A^B F$
  void adapt_quad4( void (*F)( double x, QuadParamBase *params, coef_vector *Value ), 
		    double A, double B, quad_F::one_interval *top_interval,
                    coef_vector *value );

  // ************* adapt_quad_half *******************************
  //! Adaptive 1st-order Gaussian quadrature on the interval $(A, B)$
  //! with $sqrt{1 - x}$ singularity.
  //! \param F the vector function to integrate: $\int_A^B F$
  //! \param A lower limit of integration
  //! \param B upper limit of integration
  //! \param top_interval information for adaptive quadrature
  //! \param value vector of approximate integrals: $\int_A^B F$
  void adapt_quad_half(
           void (*F)( double x, QuadParamBase *params, coef_vector *Value ), 
		    double A, double B, quad_F::one_interval *top_interval,
                    coef_vector *value );

}

namespace math_F
{
  // ---------------- Legendre ------------------
  //! Legendre polynomials
  //! Computes a vector of Legendre polynomials P_\ell( \mu )
  //! \param mu the independent variable, \mu
  //! \param value a vector of values, P_\ell( \mu )
  void Legendre( double mu, coef_vector *value );

  // ************ quadratic ******************
  //! Solves the quadratic A*alpha^2 + B*alpha + C = 0; returns the number of real roots
  //! \param A coefficient of alpha^2
  //! \param B coefficient of alpha
  //! \param C constant term
  //! \param alpha_1 the smaller real root, if there is any
  //! \param alpha_2 the larger real root, if there are 2
  int quadratic( double A, double B, double C, double *alpha_1, double *alpha_2 );

  // ------------------ Gamma_up --------------------------------
  //! Increments the incomplete Gamma function
  //!   $\Gamma( \kappa, A ) = \int_A^\infty t^{kappa-1} e^{-t}\, dt$
  //! \param kappa parameter $\kappa$
  //! \param A lower limit of the integral
  //! \param Gamma_kappa value of $\Gamma( \kappa, A )$
  //! Returns $\gamma( \kappa + 1, A )$
  double Gamma_up( double kappa, double A, double Gamma_kappa );

  // ------------------ gamma_down --------------------------------
  //! Decrements the incomplete Gamma function
  //!   $\gamma( \kappa, A ) = \int_0^A t^{kappa-1} e^{-t}\, dt$
  //! \param kappa parameter $\kappa$
  //! \param A upper limit of the integral
  //! \param Gamma_kappa value of $\gamma( \kappa, A )$
  //! Returns $\gamma( \kappa - 1, A )$
  double gamma_down( double kappa, double A, double gamma_kappa );

  // ------------------ zeroin --------------------------------
  //! Find a root of func(x, params) = target between BB and CC.
  //! \param func the function to take the root of
  //! \param target the desired equality
  //! \param BB the minimum value to consider
  //! \param CC the maximum value to consider
  //! \param params the function parameters
  //! \param tol the parameter used to set the tolerance of the seek
  //! \retval the root with accuracy of tol + 4*tol*|root|
  double zeroin(double (*func)(double, void*),
              double target,
              const dd_entry& BB,
              const dd_entry& CC,
	      void *params,
              double tol);

  // ------------------ parabola_bottom --------------------------------
  //! Fit func(x, params) by a parabola at AA, CC, and their midpoint BB.
  //! Returns the minimum MM of the parabola.
  //! If MM < BB, set CC = BB, otherwise set AA = BB.
  //! \param func the function to approximate
  //! \param AA the left-hand endpoint
  //! \param CC the right-hand endpoint
  //! \param params the function parameters
  double parabola_bottom(double (*func)(double, void*),
              dd_entry *AA,
              dd_entry *CC,
	      void *params );

}

#endif
