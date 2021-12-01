/*
 * ******** merced: calculate the transfer matrix ********
 * $Revision: 601 $
 * $Date: 2017-12-05 $
 * $Author: hedstrom $
 * $Id: adapt_quad_info.cpp 601  2017-12-05Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */

#include "adapt_quad_info.hpp"
#include "global_params.hpp"

// ************* AQinfo::adapt_quad_info *******************************
// ------------------- AQinfo::adapt_quad_info::~adapt_quad_info ----------
// Destructor
AQinfo::adapt_quad_info::~adapt_quad_info( )
{
  if( order >= 0 )
  {
    delete [] quad_weight;
  }
}
// ------------------- AQinfo::adapt_quad_info::set_order ------------------
// Sets the Legendre order and the conservation flag
void AQinfo::adapt_quad_info::set_order( int Order, Coef::Conserve cons )
{
  order = Order;
  conserve = cons;
  rough_est.order = order;
  rough_est.conserve = conserve;
}
// ------------------- AQinfo::adapt_quad_info::set_weights ------------------
// Sets the weights of the tolerances for the different Legendre orders
void AQinfo::adapt_quad_info::set_weights( )
{
  quad_weight = new double[ order + 1 ];
  static double wt = Global.Value( "quad_weight_increase" );
  double term = 1.0;
  for( int i = 0; i <= order; ++i )
  {
    quad_weight[ i ] = term;
    term *= wt;
  }
}
