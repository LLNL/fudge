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

#include "dd_vector.hpp"
#include "coef_vector.hpp"

namespace math_F
{
  // ---------------- Legendre ------------------
  //! Legendre polynomials
  //! Computes a vector of Legendre polynomials P_\ell( \mu )
  //! \param mu the independent variable, \mu
  //! \param value a vector of values, P_\ell( \mu )
  void Legendre( double mu, Coef::coef_vector *value );

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
              const Ddvec::dd_entry& BB,
              const Ddvec::dd_entry& CC,
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
             const Ddvec::dd_entry &AA,
             const Ddvec::dd_entry &CC,
	      void *params );

}

#endif
