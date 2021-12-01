/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: energy_function.hpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */
//! Base classes used for energy probability density given as a function

#ifndef ENERGY_FUNCTION_DEF
#define ENERGY_FUNCTION_DEF

#include "dd_vector.hpp"
#include "transfer.hpp"
#include "box_geom.hpp"
#include "param_base.hpp"

namespace Efunc
{
class E_function_param;  // forward declaration

//! Class for the list of intersections of the line Eout = U - Ein with an integration box
// ---------------- class U_Ein_hit_list ------------------
class U_Ein_hit_list : public Box::hit_list
{
private:
  // Implement the virtual functions
  //! Calculates E_out from E_in, for testing the side of a box
  //! \param E_in energy of the incident particle
  inline double get_Eout( double E_in ) { return E_in - eta; }

  //! Finds the intersections with the bottom of a box
  //! \param E_out bottom of the outgoing energy bin
  //! \param Ein_hits the computed intersections with E_out
  void find_bottom_hits( double E_out, std::vector< Box::Ein_Eta_Hit > *Ein_hits );

  //! Finds the intersections with the top of a box
  //! \param E_out top of the outgoing energy bin
  //! \param Ein_hits the computed intersections with E_out
  void find_top_hits( double E_out, std::vector< Box::Ein_Eta_Hit > *Ein_hits );

public:

  inline U_Ein_hit_list( ) {}
  inline ~U_Ein_hit_list( ) {}

  //! Use the Box::hit_list eta as our U value
  //! \param U the parameter used in many energy function models
  inline void set_U( double U ) { eta = U; }
  inline double get_U( ) { return eta; }
};

//! Class for a general energy probability function
// ---------------- class energy_function ------------------
class energy_function : public Ddvec::dd_vector
{
protected:
  //! The smallest incident energy for cross section, multiplicity,
  //! model weight, flux weight, and energy groups
  double E_first;

  int first_Ein;  // index of the left-hand end of the first significant energy bin
  int last_Ein;  // index of the right-hand end of the last significant energy bin

  //! Do we use interpolation?
  bool use_Eout_ints;

  //! Adds to the transfer matrix for all E_out bins for a pair of incident energies.
  //! \param transfer the computed transfer matrix
  void Eout_ladder( Trf::T_matrix& transfer, Efunc::E_function_param *Ein_param );

  //! Does the integration for one E-E' box
  //! \param transfer the computed transfer matrix
  //! \param Eout_count the current row of the transfer matrix
  void one_Ebox( Trf::T_matrix& transfer, int Eout_count,
		 Efunc::E_function_param *Ein_param );

  //! Adds to an element of transfer the integral over the E-E' box
  //! \param transfer the computed transfer matrix
  //! \param Eout_count the current row of the transfer matrix
  void update_T( Trf::T_matrix &transfer, int Eout_count,
		 Efunc::E_function_param *Ein_param );

public:
  //! For incident energy E the outgoing energy range is 0 <= E' <= E-U.
  double U;
  double threshold;

  inline energy_function( ) {}
  inline virtual ~energy_function( ) {}

  //! Default initializer for the quadrature parameters
  //! \param Eout_groups the boundaries of the outgoing energy groups
  //! \param Ein_param parameters for integration over incident energy
  void setup_data_default( const Egp::Energy_groups& Eout_groups,
    Efunc::E_function_param *Ein_param );

  //! Default set the range of incident energies for this intergration
  //! \param Ein_bin identifies the incident energy bin
  //! \param Ein_param parameters for integration over incident energy
  void set_Ein_range_default( int Ein_bin, Efunc::E_function_param *Ein_param );

  //! Default go to the next pair of incident energies.  Returns "true" when finished.
  //! \param E_in the next incident energy
  //! \param Ein_param parameters for integration over incident energy
  bool next_ladder_default( double E_in, Efunc::E_function_param *Ein_param );

// ************** virtual routines ******************************
  //! Initializer for the quadrature parameters
  //! \param Eout_groups the boundaries of the outgoing energy groups
  //! \param Ein_param parameters for integration over incident energy
  virtual void setup_data( const Egp::Energy_groups& Eout_groups,
    Efunc::E_function_param *Ein_param ) = 0;

  //! Sets the range of incident energies for this intergration
  //! \param Ein_bin identifies the incident energy bin
  //! \param Ein_param parameters for integration over incident energy
  virtual void set_Ein_range( int Ein_bin, Efunc::E_function_param *Ein_param ) = 0;

  //! Go to the next pair of incident energies.  Returns "true" when finished.
  //! \param E_in the next incident energy
  //! \param Ein_param parameters for integration over incident energy
  virtual bool next_ladder( double E_in, Efunc::E_function_param *Ein_param ) = 0;
// ************** end of virtual routines ******************************
};

//! Class for parameters for the energy spectrum models in 2d
// ---------------- class E_function_param ------------------
class E_function_param : public Pbase::param_base
{
public:
  double U;
  double top_E_out;  // top of the highest outgoing energy bin
  double E_max;  // top range of integration
  double E_1;    // actual upper limit of integration
  double multiplicity;
  double prob_integral_scale;  // scale factor for integral of probability density
  double Eprob_integral_scale;  // scale factor for integral of energy*probability density
  double norm;   // integral from 0 to top_E_out
  Coef::Conserve conserve_flag;  // conserving number or energy (used only for quadrature over E_out)
  Qmeth::Quadrature_Rule Eout_quad_rule;  // rule for integrating outgoing energy

  long int quad_count;  // number of 2-d quadratures
  long int Ein_F_count;  // number of calls to Energy_function_F::Ein_F
  long int Eout_F_count;  // number of calls to the routine for outgoing energy

  // where the upper data limit hits the quadrature box
  Efunc::U_Ein_hit_list upper_hits;

  //! pointers to the Theta data
  //! For the Watt spectrum Theta = a, and for Madland-Nix Theta = TM.
  Ddvec::dd_vector::const_iterator this_Theta;
  Ddvec::dd_vector::const_iterator next_Theta;
  Ddvec::dd_vector::const_iterator Theta_end;
  double Theta;  // the interpolated value

  E_function_param(): quad_count( 0 ), Ein_F_count( 0 ), Eout_F_count( 0 ) {}
  virtual ~E_function_param() {}

  //! Interpolate the parameters
  //! Returns true if the interpolation is OK
  //! \param E_in energy of the incident particle
  bool set_Ein_default( double E_in );


// ************** virtual routines ******************************
  //! Interpolate the parameters
  //! Returns true if the interpolation is OK
  //! \param E_in energy of the incident particle
  virtual bool set_Ein( double E_in ) = 0;

  //! Gets the integrals over outgoing energy
  //! returns true if the interpolation is OK
  //! \param Eout_0 lower limit of the outgoing energy range
  //! \param Eout_1 upper limit of the outgoing energy range
  //! \param value computed value of the indegrals
  virtual bool get_integrals( double Eout_0, double Eout_1, Coef::coef_vector &value ) = 0;

  //! Integrate from 0 to E_max to get the norm
  virtual double get_norm( ) = 0;

  //! Sets the scale factors for the integrals of probability and energy*probability
  virtual void set_scales( ) = 0;

  //! Sets the tolerance for the quadrature over incident energy
  //! \param left_E lower limit of the incident energy range
  //! \param right_E upper limit of the incident energy range
  virtual double set_tol( double left_E, double right_E ) = 0;
// ************** end of virtual routines ***********************
};

} // end of namespace Efunc

namespace Energy_function_F
{
  // ************** Probability density model ******************************
  //! Integral function for the model
  bool Ein_F( double E_in, Qparam::QuadParamBase *e_quad_param,
	      Coef::coef_vector *value );

}

#endif
