/*
 * ******** merced: calculate the transfer matrix ********
 * $Revision: 601 $
 * $Date: 2017-12-05 $
 * $Author: hedstrom $
 * $Id: adapt_quad_info.hpp 601  2017-12-05Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */

// Class for information used in adaptive Gaussian quadrature
#ifndef ADAPT_QUAD_INFO
#define ADAPT_QUAD_INFO

#include "coef_vector.hpp"

namespace AQinfo
{
  // ************* adapt_quad_info *******************************
  //! Class holds information for adaptive quadrature
  class adapt_quad_info
  {
  public:
    Coef::coef_vector rough_est;  // a rough approximation to the full integral
    double *quad_weight;  // the weights for the Legendre orders
    
    //! parameter for Richardson extrapolation
    //! For a method of order n, use Richardson = n^2
    double Richardson;

    int order;  // the Legendre order of the integral
    int depth;  // the number of levels of subdivision
    int subdivisions;  // the number of subdivisions
    Coef::Conserve conserve;  // do we conserve energy or particle number?
    bool warning_set;  // warn of inaccurate adaptive quadrature

    //! Constructor
    adapt_quad_info( ): order( -1 ),
			depth( 0 ), subdivisions( 0 ), conserve( Coef::NOT_SET ),
                      warning_set( false ) {}

    //! Destructor
    ~adapt_quad_info( );

    //! Sets the Legendre order and the conservation flag
    //! \param Order the Legendre order of the integrand
    //! \param cons conserve energy or particle number
    void set_order( int Order, Coef::Conserve cons );

    //! Sets the weights of the tolerances for the different Legendre orders
    void set_weights( );

    
  };
}  // end of namespace AQinfo

#endif
