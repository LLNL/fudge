/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2008-04-16 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: x_vector.cpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// implement the classes used for the data used in Compton and coherent scattering

#include "x_vector.hpp"
#include "messaging.hpp"
#include "global_params.hpp"

// **************** class x_vector ******************
// ----------------------- x_vector::scale_x_to_energy -------------------
// Scales x from 1/(wave length) to energy
void x_vector::scale_x_to_energy( double x_to_energy )
{
  for( x_vector::iterator ptr = begin( ); ptr != end( ); ++ptr )
  {
    ptr->x *= x_to_energy;
  }
}

// ****************** x_vector_F::get_mu_from_x ***********************
// function to get mu from x and E
double x_vector_F::get_mu_from_x( double x, double E_in )
{
  if( E_in <= 0.0 )
  {
    FatalError( "x_vector_F::get_mu_from_x",
      pastenum( "improper incident energy: ", E_in ) );
  }
  double y = x/E_in;
  double mu = 1.0 - 2*y*y;
  if( mu < -1.0 )
  {
    Warning( "x_vector_F::get_mu_from_x",
      pastenum( "changing mu to -1: ", mu ) );
    mu = -1.0;
  }
  return mu;
}
