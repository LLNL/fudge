/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: general_evap.hpp 1 2010-06-18 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */
//! Classes used for outgoing an energy distribution independent of the initial energy

#ifndef GENERAL_EVAP_DEF
#define GENERAL_EVAP_DEF

#include "dd_vector.hpp"
#include "transfer.hpp"
#include "param_base.hpp"

namespace Gevap
{
class general_evap_param;  // forward reference

//! Class for a single energy spectrum
// ---------------- class general_evap ------------------
class general_evap : public Ddvec::dd_vector
{
private:
  //! The smallest incident energy for cross section, multiplicity,
  //! model weight, flux weight, and energy groups
  double E_first;

  int first_Ein;  // index of the left-hand end of the first significant energy bin
  int last_Ein;  // index of the right-hand end of the last significant energy bin

  //! A list of integrals over E_out, independent of incident energy
  std::list< Coef::coef_vector > E_out_ints;

  //!  Gets the range of nontrivial incident energy bins; computes first_Ein and last_Ein
  //! returns true if the threshold is too high for the energy bins
  //! \param sigma the cross section data
  //! \param multiple the multiplicity of the outgoing particle
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param e_flux approximate flux used to weight the transfer matrix
  //! \param Ein_groups the boundaries of the incident energy groups
  bool get_Ein_range( const Ddvec::dd_vector& sigma, const Ddvec::dd_vector& multiple,
    const Ddvec::dd_vector& weight,
    const Lgdata::Flux_List& e_flux, const Egp::Energy_groups& Ein_groups );

  //! Sets up the integrals over the E_out bins
  //! \param order the Legendre order of the output transfer matrix
  //! \param conserve the flag for conservation of energy or particle number
  //! \param Eout_groups the boundaries of the outgoing energy groups
  void setup_Eout_ints( int order, Coef::Conserve conserve,
    const vector< double >& Eout_groups );

  //! The probability integals for outgoing energy should add to 1
  void check_sum( );

  //! Go to the next pair of incident energies.  Returns "true" when finished.
  //! \param E_in the energy of the incident particle
  //! \param Ein_param parameters for intgration over incident energy
  bool next_ladder( double E_in, Gevap::general_evap_param *Ein_param );

  //! Adds to the transfer matrix for all E_out bins for a pair of incident energies.
  //! \param transfer the computed transfer matrix
  //! \param Ein_param parameters for intgration over incident energy
  void Eout_ladder( Trf::T_matrix& transfer, Gevap::general_evap_param *Ein_param );

  //! Adds to an element of transfer the integral over the E-E' box
  //! \param transfer the computed transfer matrix
  //! \param Eout_count current row number of the computed transfer matrix
  //! \param Ein_param parameters for intgration over incident energy
  void update_T( Trf::T_matrix &transfer, int Eout_count,
		 Gevap::general_evap_param *Ein_param );

  //! The code currently handles only theta = 1
  bool theta_OK( );

public:

  //! The temperature parameter
  Ddvec::dd_vector theta;

  general_evap() {}
  ~general_evap() {}

  //! Calculates the transfer matrix for this particle
  void get_T( const Ddvec::dd_vector& sigma, const Ddvec::dd_vector& multiple, 
	      const Ddvec::dd_vector& weight, Trf::T_matrix& transfer );
};

//! Class for parameters for a single energy spectrum
// ---------------- class general_evap_param ------------------
class general_evap_param : public Pbase::param_base
{
private:

public:
  //! which outgoing energy bin
  std::list< Coef::coef_vector >::iterator current_Eout_int;

  long int quad_count;  // number of 2-d quadratures
  long int Ein_F_count;  // number of calls to general_evap_F::E_quad_F

  general_evap_param(): quad_count( 0 ), Ein_F_count( 0 ) {}
  ~general_evap_param() {}

  //! Gets the integrals over this E_out bin
  //! \param value the computed integrals of probability density over outgoing energy
  void get_integrals( Coef::coef_vector *value );
};
  
} // end of namespace Gevap

namespace general_evap_F
{
  // ************** general_evap_F::Ein_F ******************************
  //! Integral function for the model
  bool Ein_F( double E_in, Qparam::QuadParamBase *e_quad_param,
    Coef::coef_vector *value );
}

#endif
