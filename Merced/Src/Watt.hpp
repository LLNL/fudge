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

namespace Wtt
{
class Watt_param;  // forward reference

//! Class for the watt model
// ---------------- class watt ------------------
class watt : public Efunc::energy_function
{
private:

public:
  Ddvec::dd_vector b_data;

  watt( ) {}
  ~watt( ) {}

  //!  Gets the range of nontrivial incident energy bins; computes first_Ein and last_Ein
  //! returns true if the threshold is too high for the energy bins
  //! \param sigma the cross section data
  //! \param mult the multiplicity of the outgoing particle
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param e_flux approximate flux used to weight the transfer matrix
  //! \param Ein_groups the boundaries of the incident energy groups
  bool get_Ein_range( const Ddvec::dd_vector& sigma, const Ddvec::dd_vector& mult,
    const Ddvec::dd_vector& weight,
    const Lgdata::Flux_List& e_flux, const Egp::Energy_groups& Ein_groups );

  //! Calculates the transfer matrix for this particle
  //! \param sigma the cross section data
  //! \param multiple the multiplicity of the outgoing particle
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param transfer the computed transfer matrix
  void get_T( const Ddvec::dd_vector& sigma, const Ddvec::dd_vector& multiple,
    const Ddvec::dd_vector& weight, Trf::T_matrix& transfer );

  // *** Implement the virtual functions ***
  //! Initializes the quadrature parameters
  //! \param Eout_groups the boundaries of the outgoing energy groups
  //! \param Ein_param parameters for integration over incident energy
  void setup_data( const Egp::Energy_groups& Eout_groups,
		   Efunc::E_function_param *Ein_param );

  //! \param Ein_bin identifies the incident energy bin
  //! \param Ein_param parameters for integration over incident energy
  void set_Ein_range( int Ein_bin, Efunc::E_function_param *Ein_param );

  //! Go to the next pair of incident energies.  Returns "true" when finished.
  //! \param E_in the next incident energy
  //! \param Ein_param parameters for integration over incident energy
  bool next_ladder( double E_in, Efunc::E_function_param *Ein_param );
  // *** End of the virtual functions ***
};

//! Class for parameters for the Watt model
// ---------------- class Watt_param ------------------
class Watt_param : public Efunc::E_function_param
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
  Ddvec::dd_vector::const_iterator this_b;
  Ddvec::dd_vector::const_iterator next_b;
  Ddvec::dd_vector::const_iterator b_end;
  double b;  // the interpolated value

  inline Watt_param() {}
  inline ~Watt_param() {}

  // *** Implement the virtual functions ***
  //! Interpolates the data to the given incident energy
  //! Returns true if the interpolation is OK
  //! \param E_in the next incident energy
  bool set_Ein( double E_in );

  //! Gets the integrals over outgoing energy
  //! Returns true if the interpolation is OK
  //! \param E_0 lower limit of the integral over outgoing energy
  //! \param E_1 upper limit of the integral over outgoing energy
  //! \param value the computed integral
  bool get_integrals( double Eout_0, double Eout_1, Coef::coef_vector &value );

  // Integrate from 0 to E_max to get the norm
  double get_norm( );

  //! Sets the scale factors for the integrals of probability and energy*probability
  void set_scales( );

  //! Sets the tolerance for the quadrature over incident energy
  //! \param left_E lower limit of the integral over incident energy
  //! \param right_E upper limit of the integral over incident energy
  double set_tol( double left_E, double right_E );
  // *** End of the virtual functions ***
};

//! Class for parameters to check the Watt model
// ---------------- class check_Watt_params ------------------
class check_Watt_params : public Qparam::QuadParamBase
{
public:
  double a;  // the Watt parameter
  double b;  // the Watt parameter

  inline check_Watt_params() {}
  inline ~check_Watt_params() {}

};

} // end of namespace Wtt

namespace Watt_F
{
  // **************** Function to integrate *********************
  // Function for the 1-d quadrature over outgoing energy
  //! Returns true if the interpolation is OK
  //! \param E_out, the outgoing energy
  //! \param E_out_param, the quadrature parameters
  //! \param value, the computed integrand
  bool check_Watt_F( double E_out, Qparam::QuadParamBase *E_out_param,
		     Coef::coef_vector *value );
}

#endif
