/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2008-04-16 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: x_vector.hpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 * 
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// define the classes used for the data used in Compton and coherent scattering

#ifndef X_VECTOR_CLASS
#define X_VECTOR_CLASS

#include "dd_vector.hpp"

namespace Xvec
{
// ----------------------- class x_vector -------------------
//! Special Ddvec::dd_vector used for Compton and coherent scattering data
class x_vector : public Ddvec::dd_vector
{
private:

public:
  //! Default constructor
  x_vector( ) {}

  //! Default destructor
  ~x_vector( ) {}

  //! Scales x from 1/(wave length) to energy
  //! \param x_to_energy the scale factor
  void scale_x_to_energy( double x_to_energy );
};
} // end of namespace Xvec

namespace x_vector_F
{
  // *************** function get_mu_from_x ************
  //! function to get mu from x and E
  //! \param x ( E_in/(c*h) )*\sqrt{ ( 1 - \mu )/2}
  //! \aram E_in energy of incident gamma
  double get_mu_from_x( double x, double E_in );
}

#endif
