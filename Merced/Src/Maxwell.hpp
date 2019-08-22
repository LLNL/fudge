/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: Maxwell.hpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */
//! Classes used for Maxwell energy probability density

#ifndef MAXWELL_DEF
#define MAXWELL_DEF

#include "evaporation.hpp"

class Maxwell_param;  // forward reference

//! Class for the maxwell model
// ---------------- class maxwell ------------------
class maxwell : public evaporation
{
private:

public:
  maxwell( ) {}
 
  ~maxwell( ) {}

  //! Calculates the transfer matrix for this particle
  //! \param sigma the cross section data
  //! \param multiple the multiplicity of the outgoing particle
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param transfer the computed transfer matrix
  void get_T( const dd_vector& sigma, const dd_vector& multiple,
    const dd_vector& weight, T_matrix& transfer );
};

//! Class for parameters for the Maxwell model
// ---------------- class Maxwell_param ------------------
class Maxwell_param : public Evap_param
{
private:
  //! Integral of the probability density from 0 to E_1
  //! \param E_1 upper limit of integration in outgoing energy
  double integral_prob( double E_1 );

  //! Integral of energy times the probability density from 0 to E_1
  //! \param E_1 upper limit of integration in outgoing energy
  double integral_Eprob( double E_1 );

public:
  inline Maxwell_param() {}
  inline ~Maxwell_param() {}

  // *** Implement the virtual functions ***
  //! Gets the integrals over outgoing energy
  //! \param Eout_0 lower limit of the integral over outgoing energy
  //! \param Eout_1 upper limit of the integral over outgoing energy
  //! \param value the computed integral
  void get_integrals( double Eout_0, double Eout_1, coef_vector &value );

  //! Gets the integrals over outgoing energy and returns the noise in the calculation
  //! \param Eout_0 lower limit of the integral over outgoing energy
  //! \param Eout_1 upper limit of the integral over outgoing energy
  //! \param value the computed integral
  inline double tol_get_integrals( double Eout_0, double Eout_1, coef_vector &value )
  {
    return Evap_param::tol_get_integrals( Eout_0, Eout_1, value );
  }

  //! Integrate from 0 to E_max to get the norm
  double get_norm( );

  //! Sets the scale factors for the integrals of probability and energy*probability
  void set_scales( );

  //! Sets the tolerance for the quadrature over incident energy
  //! \param left_E lower limit of the integral over incident energy
  //! \param right_E upper limit of the integral over incident energy
  inline double set_tol( double left_E, double right_E )
  {
    return Evap_param::set_tol( left_E, right_E );
  }
  // *** End of the virtual functions ***
};

// **************** Function to integrate *********************
namespace Maxwell_F
{
  // --------------------  Eout_F ------------------
  // Function for the 1-d quadrature over outgoing energy
  //! \param E_out energy of outgoing particle
  //! \param E_out_param parameters for this function
  //! \param value the computed contribution to the transfer matrix
  void Eout_F( double E_out, QuadParamBase *E_out_param, coef_vector *value );
}

#endif
