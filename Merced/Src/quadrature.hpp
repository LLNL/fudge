/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2014-03-21 19:06:56 -0800 (Fri, 21 Mar 2014) $
 * $Author: hedstrom $
 * $Id: quadrature.hpp 1 2008-07-01 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// Define the quadrature methods

#ifndef QUADRATURE_METHODS
#define QUADRATURE_METHODS

// ------------------------ Quadrature_Method ---------------
//! Specifies the quadrature rule.  The options are
//! ADAPTIVE2: adaptive quadrature based on 2nd-order Gaussian quadrature
//! ADAPTIVE4: adaptive quadrature based on 4th-order Gaussian quadrature
//! GAUSS2: 2nd-order Gaussian quadrature
//! GAUSS4: 4th-order Gaussian quadrature
//! GAUSS6: 6th-order Gaussian quadrature
//! GAUSS10: 10th-order Gaussian quadrature
//! ADAPT_HALF: 1st-order Gaussian quadrature with $sqrt{1 - x}$ singularity
//! EXACT: For Legendre double-differential data in the lab frame do exact
//!   integration.  This option is no longer implemented---too slow.
enum Quadrature_Method{ ADAPTIVE2, ADAPTIVE4, ADAPT_HALF, GAUSS2, GAUSS4,
  GAUSS6, GAUSS10, EXACT };

#endif
