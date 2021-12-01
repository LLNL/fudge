/*
 * ******** merced: calculate the transfer matrix ********
 * $Revision: 601 $
 * $Date: 2017-12-05 $
 * $Author: hedstrom $
 * $Id: adapt_quad.hpp 601  2017-12-05Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */

// Defines the general class for adaptive quadrature

#ifndef ADAPT_QUAD_CLASS
#define ADAPT_QUAD_CLASS

#include <cmath>

#include "adapt_quad_info.hpp"
#include "param_base.hpp"
#include "coef_vector.hpp"
#include "quad_methods.hpp"
#include "global_params.hpp"
#include "messaging.hpp"

//! We use Walter Gander's adaptive quadrature routine.
//! See W. Gander and W. Gautschi, "Adaptive quadrature--revisited",
//! BIT 40 (2000), 84-101.
//! The outines are modified to integrate Coef::coef_vector-valued functions and 
//! to let the user choose the quadrature method in place of Simpson's rule.

namespace quad_F
{
  // ************* one_interval *******************************
  //! Information for one interval in the quad_F::adapt_quad routine
  class one_interval
  {
  private:

  public:
    double a;  // left end point
    double b;  // right end point
    
    Coef::coef_vector integral;  // the approximate integral in this subinterval
    Qparam::QuadParamBase *params;  // the function parameters for this integration
    AQinfo::adapt_quad_info *quad_info;  // information for adaptive quadrature

    //! Constructor
    one_interval( ) {}

    //! Constructor
    //! \param Order the Legendre order of the integrand
    //! \param cons do we conserve energy or particle number?
    one_interval( int Order, Coef::Conserve cons );

    //! Destructor
    inline ~one_interval( ) {}

    //! Allocate space for the vectors
    //! \param Order the Legendre order of the integrand
    //! \param cons do we conserve energy or particle number?
    void set_order( int Order, Coef::Conserve cons );

    //! Tests whether two estimates are sufficiently close
    //! \param left_half approximate integral over left half subinterval
    //! \param right_half approximate integral over right half subinterval
    bool test_OK( const Coef::coef_vector& left_half, const Coef::coef_vector& right_half )
      const;

  };

  // ************* quad_list *******************************
  // The linked list used to hold the quadrature routine adapt_quad
  class quad_list : public std::list< one_interval >
  {
  private:
    bool (*FF)( double x, Qparam::QuadParamBase *params, Coef::coef_vector *value );  // the function to integrate
    double AA;   // left end point
    double BB;   // right end point
    Qparam::QuadParamBase *params;  // the function paramters, including conservation and Legendre order

    //! set up the quadrature
    //! Returns false if the quadrature fails
    bool initialize( );

    //! set up the rough estimate for comparison
    //! Returns false if the quadrature fails
    //! tol, the quadrature tolerance
    //! \param Value, returned if the rough estimate is OK
    bool get_rough( double tol, Coef::coef_vector *Value );

    //! test a link for accuracy
    //! \param link, the link to test
    //! \param left_half, the integral on the left half interval
    //! \param right_half, the integral on the right half interval
    //! \param Value, returned if link is OK
    void test_link( quad_F::quad_list::iterator link,
        const Coef::coef_vector &left_half, const Coef::coef_vector &right_half,
        Coef::coef_vector *Value );
    
    //! Process one level of subdivision
    //! Returns false if the quadrature fails
    //! \param Value, returned if link is OK
    bool one_level( Coef::coef_vector *Value );
    
  public:
    AQinfo::adapt_quad_info quad_info;

    quad_list( ) {}

    virtual ~quad_list() {}

    // ************** virtual routines ******************************
    //! The quadrature method on the interval $(A, B)$.
    //! Returns false if the quadrature fails
    //! \param F the vector function to integrate: $\int_A^B F$
    //! \param A lower limit of integration
    //! \param B upper limit of integration
    //! \param params parameters for the function, F
    //! \param value vector of approximate integrals: $\int_A^B F$
    virtual bool Quad_Method( bool (*F)( double x, Qparam::QuadParamBase *params, Coef::coef_vector *Value ), 
       double A, double B, Qparam::QuadParamBase *params,
			      Coef::coef_vector *value ) = 0;

    // ************** end of virtual routines ***********************
    
    //! Set up the adaptive integration
    //! \param F the vector function to integrate: $\int_A^B F$
    //! \param A lower limit of integration
    //! \param B upper limit of integration
    //! \param Params, the quadrature parameters
    void set_up( bool (*F)( double x, Qparam::QuadParamBase *params,
			    Coef::coef_vector *Value ), 
		 double A, double B, Qparam::QuadParamBase *Params );

    //! Evaluate the integral
    //! Returns false if the quadrature fails
    //! tol, the quadrature tolerance
    //! \param Value, the computed integral
    bool adapt_quad( double tol, Coef::coef_vector *Value );
  };

  // ------------- adapt_Gauss_1 ----------------------
  // Adaptive 1st-order Gaussian quadrature on the interval $(A, B)$.
  class adapt_Gauss_1 : public quad_F::quad_list
  {
  private:
  public:

    adapt_Gauss_1( ) {}

    ~adapt_Gauss_1( ) { }

    bool Quad_Method( bool (*F)( double x, Qparam::QuadParamBase *params,
				 Coef::coef_vector *Value ), 
		 double A, double B, Qparam::QuadParamBase *params,
		      Coef::coef_vector *value )
    {
      return Qmeth::midpoint_rule( F, A, B, params, value );
    }
  };

  // ------------- adapt_Gauss_2 ----------------------
  // Adaptive 2nd-order Gaussian quadrature on the interval $(A, B)$.
  class adapt_Gauss_2 : public quad_F::quad_list
  {
  private:
  public:

    adapt_Gauss_2( ) {}

    ~adapt_Gauss_2( ) { }

    bool Quad_Method( bool (*F)( double x, Qparam::QuadParamBase *params,
				 Coef::coef_vector *Value ), 
		 double A, double B, Qparam::QuadParamBase *params,
		      Coef::coef_vector *value )
    {
      return Qmeth::Gauss_2( F, A, B, params, value );
    }
  };

  // ------------- adapt_Gauss_3 ----------------------
  // Adaptive 3rd-order Gaussian quadrature on the interval $(A, B)$.
  class adapt_Gauss_3 : public quad_F::quad_list
  {
  private:
  public:

    adapt_Gauss_3( ) {}

    ~adapt_Gauss_3( ) { }

    bool Quad_Method( bool (*F)( double x, Qparam::QuadParamBase *params,
				 Coef::coef_vector *Value ), 
		 double A, double B, Qparam::QuadParamBase *params,
		      Coef::coef_vector *value )
    {
      return Qmeth::Gauss_3( F, A, B, params, value );
    }
  };

  // ------------- adapt_Gauss_4 ----------------------
  //! Adaptive 4th-order Gaussian quadrature on the interval $(A, B)$.
  class adapt_Gauss_4 : public quad_F::quad_list
  {
  private:
  public:

    adapt_Gauss_4( ) {}

    ~adapt_Gauss_4( ) { }

    bool Quad_Method( bool (*F)( double x, Qparam::QuadParamBase *params,
				 Coef::coef_vector *Value ), 
		 double A, double B, Qparam::QuadParamBase *params,
		      Coef::coef_vector *value )
    {
      return Qmeth::Gauss_4( F, A, B, params, value );
    }
  };

  // ------------- adapt_Gauss_wt_L1 ----------------------
  //! Adaptive 1st-order Gaussian quadrature on the interval $(A, B)$.
  class adapt_Gauss_wt_L1 : public quad_F::quad_list
  {
  private:
  public:

    adapt_Gauss_wt_L1( ) {}

    ~adapt_Gauss_wt_L1( ) { }

    bool Quad_Method( bool (*F)( double x, Qparam::QuadParamBase *params,
				 Coef::coef_vector *Value ), 
		 double A, double B, Qparam::QuadParamBase *params,
		      Coef::coef_vector *value )
    {
      return Qmeth::weight_L1( F, A, B, params, value );
    }
  };

  // ************* integrate *******************************
  //! Integrates F over the interval (A, B).
  //! Returns false if the quadrature fails
  //! \param F the vector function to integrate: $\int_A^B F$
  //! \param use_quad, the quadrature rule to use
  //! \param A lower limit of integration
  //! \param B upper limit of integration
  //! \param params, the quadrature parameters
  //! \param tol, tolerance for adaptive quadrature
  //! \param value vector of approximate integrals: $\int_A^B F$
  bool integrate( bool (*F)( double x, Qparam::QuadParamBase *params,
			     Coef::coef_vector *Value ),
                  Qmeth::Quadrature_Rule use_quad, double A, double B,
                  Qparam::QuadParamBase *params, double tol, Coef::coef_vector *value );

}  // end of namespace quad_F

#endif


