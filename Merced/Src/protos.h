/*
 * ******** getTransferMatrix: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: protos.h 1 2009-09-01 03:06:56Z hedstrom $
 * ******** getTransferMatrix: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
  Copyright 2022, Lawrence Livermore National Security, LLC.
  See the top-level COPYRIGHT file for details.
  
  SPDX-License-Identifier: BSD-3-Clause
 * # <<END-copyright>>
 */
/* header for the C code */
#ifndef CROUTINES
#define CROUTINES

#include <string>

namespace Proto
{
/* 
 * exponential integral $\int_1^\infty dt \, \exp{-xt}/t^n$
 * \param n the exponent in the denominator
 * \param x parameter in the exponent
*/
 double expn ( int n, double x );

/* 
 * incomplete gamma function $\int_0^x dt \, t^{a - 1} \exp{-t}$
 * \param a for the exponent of t
 * \param x upper limit of integration
*/
 double igam ( double a, double x );

/* 
 * incomplete gamma function $\int_x^\infty dt \, t^{a - 1} \exp{-t}$
 * \param a for the exponent of t
 * \param x lower limit of integration
*/
 double igamc ( double a, double x );
  
/*
 * flags a math error in a C routine
 * \param name the name of the routine in which the error occurs
 * \param code the type of error
 */
 int mtherr ( const std::string &name, int code );

} // end of namespace Proto

#endif
