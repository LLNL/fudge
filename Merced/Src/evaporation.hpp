/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: evaporation.hpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */
//! Classes used for the evaporation model energy probability density

#ifndef ENERGY_MODELS_DEF
#define ENERGY_MODELS_DEF

#include "energy_function.hpp"

class Evap_param;  // forward reference

namespace Evap
{
//! Class for parameters for integration over outgoing energy
// ---------------- class evap_params_1d ------------------
class evap_params_1d : public Qparam::QuadParamBase
{
public:
  double Theta;  // the model parameter

  inline evap_params_1d() {}
  inline ~evap_params_1d() {}
};

//! Class for the evaporation model
// ---------------- class evaporation ------------------
class evaporation : public Efunc::energy_function
{
protected:

public:

  evaporation( ) {}
  ~evaporation( ) {}

  //!  Gets the range of nontrivial incident energy bins; computes first_Ein and last_Ein
  //! returns true if the threshold is too high for the energy bins
  //! \param sigma the cross section data
  //! \param mult the multiplicity of the outgoing particle
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param e_flux approximate flux used to weight the transfer matrix
  //! \param Ein_groups the boundaries of the incident energy groups
  bool get_Ein_range( const dd_vector& sigma, const dd_vector& mult,
    const dd_vector& weight,
    const Lgdata::Flux_List& e_flux, const Egp::Energy_groups& Ein_groups );

  //! Calculates the transfer matrix for this particle
  //! \param sigma the cross section data
  //! \param multiple the multiplicity of the outgoing particle
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param transfer the computed transfer matrix
  void get_T( const dd_vector& sigma, const dd_vector& multiple,
    const dd_vector& weight, Trf::T_matrix& transfer );

// ************** implement virtual routines ******************************
  //! Initializer for the quadrature parameters
  //! \param Eout_groups the boundaries of the outgoing energy groups
  //! \param Ein_param parameters for integration over incident energy
  inline void setup_data( const Egp::Energy_groups& Eout_groups,
      Efunc::E_function_param *Ein_param )
  {
    setup_data_default( Eout_groups, Ein_param );
  }

  //! Sets the range of incident energies for this intergration
  //! \param Ein_bin identifies the incident energy bin
  //! \param Ein_param parameters for integration over incident energy
  inline void set_Ein_range( int Ein_bin, Efunc::E_function_param *Ein_param )
  {
    set_Ein_range_default( Ein_bin, Ein_param );
  }

  //! Go to the next pair of incident energies.  Returns "true" when finished.
  //! \param E_in the next incident energy
  //! \param Ein_param parameters for integration over incident energy
  inline bool next_ladder( double E_in, Efunc::E_function_param *Ein_param )
  {
    return next_ladder_default( E_in, Ein_param );
  }

// ************** end of virtual routines ******************************
};
} // end of namespace Evap

//! Class for parameters for the evaporation model
// ---------------- class Evap_param ------------------
class Evap_param : public Efunc::E_function_param
{
private:
  //! Integral of the probability density from 0 to E_1
  //! \param E_1 upper limit of the integral over outgoing energy
  double integral_prob( double E_1 );

  //! Integral of energy times the probability density from 0 to E_1
  //! \param E_1 upper limit of the integral over outgoing energy
  double integral_Eprob( double E_1 );

public:
  inline Evap_param() {}
  inline ~Evap_param() {}

  // *** Implement the virtual functions ***
  //! Interpolate the parameters
  //! Returns true if the interpolation is OK
  //! \param E_in the energy of the incident particle
  inline bool set_Ein( double E_in )
  {
    return set_Ein_default( E_in );
  }

  //! Gets the integrals over outgoing energy
  //! Returns true if the interpolation is OK
  //! \param Eout_0 lower limit of the integral over outgoing energy
  //! \param Eout_1 upper limit of the integral over outgoing energy
  //! \param value the computed integral
  bool get_integrals( double Eout_0, double Eout_1, Coef::coef_vector &value );

  //! Integrate from 0 to E_max to get the norm
  double get_norm( );

  //! Sets the scale factors for the integrals of probability and energy*probability
  void set_scales( );

  //! Sets the tolerance for the quadrature over incident energy
  //! \param left_E lower limit of the integral over incident energy
  //! \param right_E upper limit of the integral over incident energy
  double set_tol( double left_E, double right_E );
  // *** End of the virtual functions ***
};

namespace evaporation_F
{
  // **************** Function to integrate *********************
  // --------------------  Eout_F ------------------
  // Function for the 1-d quadrature over outgoing energy
  //! Returns true if the interpolation is OK
  //! \param E_out energy of outgoing particle
  //! \param E_out_param parameters for this function
  //! \param value the computed contribution to the transfer matrix
  bool Eout_F( double E_out, Qparam::QuadParamBase *E_out_param,
	       Coef::coef_vector *value );
}

#endif
