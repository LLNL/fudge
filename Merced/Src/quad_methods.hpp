/*
 * ******** merced: calculate the transfer matrix ********
 * $Revision: 601 $
 * $Date: 2017-12-12 $
 * $Author: hedstrom $
 * $Id: quad_methods.hpp 601  2017-12-12Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */

// Classes for Gaussian quadrature

#ifndef QUAD_METHODS
#define QUAD_METHODS

#include "param_base.hpp"
#include "coef_vector.hpp"

namespace Qmeth
{
// ------------------------ Quadrature_Method ---------------
//! Specifies the quadrature rule.  The options are
//! GAUSS1: the midpoint rule
//! GAUSS2: 2nd-order Gaussian quadrature
//! GAUSS3: 3nd-order Gaussian quadrature
//! GAUSS4: 4th-order Gaussian quadrature
//! GAUSS6: 6th-order Gaussian quadrature
//! GAUSS10: 10th-order Gaussian quadrature
//! WEIGHT_L1: 1st-order Gaussian quadrature with $sqrt{x}$ singularity
enum Quadrature_Method{ GAUSS1, GAUSS2, GAUSS3, GAUSS4, GAUSS6, GAUSS10,
			WEIGHT_L1 };


//! Class to specify the quadrature rule
// ------------------------ Quadrature_Rule ---------------
class Quadrature_Rule
{
private:
public:
  Quadrature_Method quad_method;

  //! Do we use adaptive quadrature?
  bool adaptive;

  //! Is this rule set at input?
  bool input_set;

  Quadrature_Rule( ): adaptive( true ), input_set( false ) {}
  ~Quadrature_Rule( ) { }

  //! to copy
  //! \param to_copy, the data to copy
  Quadrature_Rule& operator=( const Quadrature_Rule &to_copy );
};

  // ********** basic Gaussian quadrature routines **************
  // ------------- midpoint_rule -------------------------------
  //! The midpoint rule on the interval $(A, B)$.
  //! Returns true if there were no problems
  //! \param F the vector function to integrate: $\int_A^B F$
  //! \param A lower limit of integration
  //! \param B upper limit of integration
  //! \param params parameters for the function, F
  //! \param value vector of approximate integrals: $\int_A^B F$
  bool midpoint_rule( bool (*F)( double x, Qparam::QuadParamBase *params,
				 Coef::coef_vector *Value ), 
    double A, double B, Qparam::QuadParamBase *params,
		      Coef::coef_vector *value );


  // ------------- Gauss_2 -------------------------------
  //! Second-order Gaussian quadrature on the interval $(A, B)$.
  //! Returns true if there were no problems
  //! This is used because it is exact for polynomials of degree 3.
  //! \param F the vector function to integrate: $\int_A^B F$
  //! \param A lower limit of integration
  //! \param B upper limit of integration
  //! \param params parameters for the function, F
  //! \param value vector of approximate integrals: $\int_A^B F$
  bool Gauss_2( bool (*F)( double x, Qparam::QuadParamBase *params,
			   Coef::coef_vector *Value ), 
		double A, double B, Qparam::QuadParamBase *params,
		Coef::coef_vector *value );

  // ------------- Gauss_3 -------------------------------
  //! Third-order Gaussian quadrature on the interval $(A, B)$.
  //! Returns true if there were no problems
  //! This is used because it is exact for polynomials of degree 3.
  //! \param F the vector function to integrate: $\int_A^B F$
  //! \param A lower limit of integration
  //! \param B upper limit of integration
  //! \param params parameters for the function, F
  //! \param value vector of approximate integrals: $\int_A^B F$
  bool Gauss_3( bool (*F)( double x, Qparam::QuadParamBase *params,
			   Coef::coef_vector *Value ), 
		double A, double B, Qparam::QuadParamBase *params,
		Coef::coef_vector *value );

  // ------------- Gauss_4 -------------------------------
  //! Fourth-order Gaussian quadrature on the interval $(A, B)$.
  //! Returns true if there were no problems
  //! This is used because it is exact for polynomials of degree 7.
  //! \param F the vector function to integrate: $\int_A^B F$
  //! \param A lower limit of integration
  //! \param B upper limit of integration
  //! \param params parameters for the function, F
  //! \param value vector of approximate integrals: $\int_A^B F$
  bool Gauss_4( bool (*F)( double x, Qparam::QuadParamBase *params,
			   Coef::coef_vector *Value ), 
		double A, double B, Qparam::QuadParamBase *params,
		Coef::coef_vector *value );

  // ------------- Gauss_6 -------------------------------
  //! Sixth-order Gaussian quadrature on the interval $(A, B)$.
  //! Returns true if there were no problems
  //! This is used because it is exact for polynomials of degree 11.
  //! \param F the vector function to integrate: $\int_A^B F$
  //! \param A lower limit of integration
  //! \param B upper limit of integration
  //! \param params parameters for the function, F
  //! \param value vector of approximate integrals: $\int_A^B F$
  bool Gauss_6( bool (*F)( double x, Qparam::QuadParamBase *params,
			   Coef::coef_vector *Value ), 
		double A, double B, Qparam::QuadParamBase *params,
		Coef::coef_vector *value );

  // ------------- Gauss_10 -------------------------------
  //! Tenth-order Gaussian quadrature on the interval $(A, B)$.
  //! Returns true if there were no problems
  //! This is used because it is exact for polynomials of degree 19.
  //! \param F the vector function to integrate: $\int_A^B F$
  //! \param A lower limit of integration
  //! \param B upper limit of integration
  //! \param params parameters for the function, F
  //! \param value vector of approximate integrals: $\int_A^B F$
  bool Gauss_10( bool (*F)( double x, Qparam::QuadParamBase *params,
			    Coef::coef_vector *Value ), 
		double A, double B, Qparam::QuadParamBase *params,
		 Coef::coef_vector *value );

  // ------------- weight_L1 -------------------------------
  //! First-order Gaussian quadrature on the interval $(A, B)$
  //! with $sqrt{x}$ singularity.
  //! Returns true if there were no problems
  //! \param F the vector function to integrate: $\int_A^B F$
  //! \param A lower limit of integration
  //! \param B upper limit of integration
  //! \param params parameters for the function, F
  //! \param value vector of approximate integrals: $\int_A^B F$
  bool weight_L1( bool (*F)( double x, Qparam::QuadParamBase *params,
			     Coef::coef_vector *Value ), 
		double A, double B, Qparam::QuadParamBase *params,
		  Coef::coef_vector *value );

} // end of namespace Qmeth

#endif
