/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2015-08-04 19:06:56 -0800 (Tue, 04 Aug 2015) $
 * $Author: hedstrom $
 * $Id: max_output_order.hpp 1 2015-08-04 03:06:56Z hedstrom $
 * initial code
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/

// header for the coef_vector class
#ifndef MAX_LEGENDRE
#define MAX_LEGENDRE

namespace Order
{
  //! Sets the maximum Legendre order of the output transfer matrix
  //! The code runs a lot faster if this is hard-wired.
  static const int max_output_order = 15;
}

#endif
