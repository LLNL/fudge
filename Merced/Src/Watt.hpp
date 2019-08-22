/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: Watt.hpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */
//! Classes used for Watt energy probability density

#ifndef WATT_DEF
#define WATT_DEF

#include "energy_function.hpp"

class Watt_param;  // forward reference

//! Class for the watt model
// ---------------- class watt ------------------
class watt : public energy_function
{
private:

public:
  dd_vector b_data;

  watt( ) {}
  ~watt( ) {}

  //!  Gets the range of nontrivial incident energy bins; computes first_Ein and last_Ein
  //! returns true if the threshold is too high for the energy bins
  //! \param sigma the cross section data
  //! \param mult the multiplicity of the outgoing particle
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param e_flux approximate flux used to weight the transfer matrix
  //! \param Ein_groups the boundaries of the incident energy groups
  bool get_Ein_range( const dd_vector& sigma, const dd_vector& mult,
    const dd_vector& weight,
    const Flux_List& e_flux, const Energy_groups& Ein_groups );

  //! Calculates the transfer matrix for this particle
  //! \param sigma the cross section data
  //! \param multiple the multiplicity of the outgoing particle
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param transfer the computed transfer matrix
  void get_T( const dd_vector& sigma, const dd_vector& multiple,
    const dd_vector& weight, T_matrix& transfer );

  // *** Implement the virtual functions ***
  //! Initializes the quadrature parameters
  //! \param Eout_groups the boundaries of the outgoing energy groups
  //! \param Ein_param parameters for integration over incident energy
  void setup_data( const Energy_groups& Eout_groups,
		   E_function_param *Ein_param );

  //! \param Ein_bin identifies the incident energy bin
  //! \param Ein_param parameters for integration over incident energy
  void set_Ein_range( int Ein_bin, E_function_param *Ein_param );

  //! Go to the next pair of incident energies.  Returns "true" when finished.
  //! \param E_in the next incident energy
  //! \param Ein_param parameters for integration over incident energy
  bool next_ladder( double E_in, E_function_param *Ein_param );
  // *** End of the virtual functions ***
};

//! Class for parameters for the Watt model
// ---------------- class Watt_param ------------------
class Watt_param : public E_function_param
{
private:
  //! Indefinite integral of the probability density
  //! \param E_0 lower limit of the integral over outgoing energy
  //! \param E_1 upper limit of the integral over outgoing energy
  double integral_prob( double E_0, double E_1 );

  //! Indefinite integral of energy times the probability density
  //! \param E_0 lower limit of the integral over outgoing energy
  //! \param E_1 upper limit of the integral over outgoing energy
  double integral_Eprob( double E_0, double E_1 );

  //! A robust integration of \int_{\alpha_0}^{\alpha_1} w \exp{-w^2} \, dw
  //! \param alpha_0 lower limit of the integral
  //! \param alpha_1 upper limit of the integral
  double w_integral( double alpha_0, double alpha_1 );

  //! A robust integration of \int_{\alpha_0}^{\alpha_1} \exp{-w^2} \, dw
  //! \param alpha_0 lower limit of the integral
  //! \param alpha_1 upper limit of the integral
  double one_integral( double alpha_0, double alpha_1 );

  // Do the integration by adaptive quadrature
  //! \param E_0 lower limit of the integral over outgoing energy
  //! \param E_1 upper limit of the integral over outgoing energy
  double check_Watt( double E_0, double E_1 );

public:
  //! pointers to the b data
  dd_vector::const_iterator this_b;
  dd_vector::const_iterator next_b;
  dd_vector::const_iterator b_end;
  double b;  // the interpolated value

  inline Watt_param() {}
  inline ~Watt_param() {}

  // *** Implement the virtual functions ***
  //! Interpolates the data to the given incident energy
  //! \param E_in the next incident energy
  void set_Ein( double E_in );

  //! Gets the integrals over outgoing energy
  //! \param E_0 lower limit of the integral over outgoing energy
  //! \param E_1 upper limit of the integral over outgoing energy
  //! \param value the computed integral
  void get_integrals( double Eout_0, double Eout_1, coef_vector &value );

  // Integrate from 0 to E_max to get the norm
  double get_norm( );

  //! Sets the scale factors for the integrals of probability and energy*probability
  void set_scales( );

  //! Gets the integrals over outgoing energy and returns the noise in the calculation
  //! \param E_0 lower limit of the integral over outgoing energy
  //! \param E_1 upper limit of the integral over outgoing energy
  //! \param value the computed integral
  double tol_get_integrals( double Eout_0, double Eout_1, coef_vector &value );

  //! Sets the tolerance for the quadrature over incident energy
  //! \param left_E lower limit of the integral over incident energy
  //! \param right_E upper limit of the integral over incident energy
  double set_tol( double left_E, double right_E );
  // *** End of the virtual functions ***
};

//! Class for parameters to check the Watt model
// ---------------- class check_Watt_params ------------------
class check_Watt_params : public QuadParamBase
{
public:
  double a;  // the Watt parameter
  double b;  // the Watt parameter

  inline check_Watt_params() {}
  inline ~check_Watt_params() {}

};

namespace Watt_F
{
  // **************** Function to integrate *********************
  // Function for the 1-d quadrature over outgoing energy
  void check_Watt_F( double E_out, QuadParamBase *E_out_param, coef_vector *value );
}

#endif
