/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2008-07-01 19:06:56 -0800 (Tue, 01 Jul 2008) $
 * $Author: hedstrom $
 * $Id: coherent.hpp 1 2008-07-01 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// define the classes used for coherent scattering

#ifndef COHERENT_CLASS
#define COHERENT_CLASS

#include <iostream>
#include "data_parser.hpp"
#include "transfer.hpp"
#include "param_base.hpp"
#include "box_geom.hpp"
#include "x_vector.hpp"

namespace Coh
{
  class coherent_Ein_param;  // forward declaration

//! Class to hold pointers to real and imaginary anomalous scattering factors
// ---------------- class anomalous_ptrs ------------------
class anomalous_ptrs
{
public:

  // pointers to real anomalous factor data
  Ddvec::dd_vector::const_iterator Re_anomalous_ptr;
  Ddvec::dd_vector::const_iterator next_Re_anomalous;
  Ddvec::dd_vector::const_iterator start_Re_anomalous;
  Ddvec::dd_vector::const_iterator end_Re_anomalous;
  Ddvec::dd_vector::const_iterator first_ladder_Re_anomalous;
  Ddvec::dd_vector::const_iterator last_ladder_Re_anomalous;

  // pointers to imaginary anomalous factor data
  Ddvec::dd_vector::const_iterator Im_anomalous_ptr;
  Ddvec::dd_vector::const_iterator next_Im_anomalous;
  Ddvec::dd_vector::const_iterator start_Im_anomalous;
  Ddvec::dd_vector::const_iterator end_Im_anomalous;
  Ddvec::dd_vector::const_iterator first_ladder_Im_anomalous;
  Ddvec::dd_vector::const_iterator last_ladder_Im_anomalous;

  anomalous_ptrs( ) {}
  ~anomalous_ptrs( ) {}

  //! Gets the real anomalous factor at this energy
  //! \param E_in the incident energy
  //! \param interp_OK, true if the interpolation is OK
  double get_Re_anomalous( double E_in, bool *interp_OK );

  //! Gets the imaginary anomalous factor at this energy
  //! \param E_in the incident energy
  //! \param interp_OK, true if the interpolation is OK
  double get_Im_anomalous( double E_in, bool *interp_OK );

  //! Initializes the pointers to real anomalous factor data
  //! \param Re_anomalous the real anomalous factor data
  void init_Re_anomalous( Ddvec::dd_vector &Re_anomalous );

  //! Initializes the pointers to imaginary anomalous factor data
  //! \param Im_anomalous the imaginary anomalous factor data
  void init_Im_anomalous( Ddvec::dd_vector &Im_anomalous );

  //! Initializes the pointers to real anomalous factor data
  //! \param Re_anomalous the real anomalous factor data, count down
  void init_Re_down( Ddvec::dd_vector &Re_anomalous );

  //! Initializes the pointers to imaginary anomalous factor data
  //! \param Im_anomalous the imaginary anomalous factor data, count down
  void init_Im_down( Ddvec::dd_vector &Im_anomalous );

  //! Sets the range of data for this quadrature
  //! \param E_0 lower incident energy
  //! \param E_1 higher incident energy
  void set_range( double E_0, double E_1 );

  //! Increments the pointers to real anomalous factor data
  //! \param E_in the next incident energy, count down
  void get_next_Re_down( double E_in );

  //! Increments the pointers to imaginary anomalous factor data
  //! \param  E_in the next incident energy, count down
  void get_next_Im_down( double E_in );

  //! Increments the pointers to anomalous factor data, count up
  //! \param right_E, right-hand end of most recent subinterval
  //! \param Ein_1, right-hand end of original interval
  bool next_range( double right_E, double Ein_1 );
};

//! Class for parameters for the 1-d quadrature of the coherent form factor over mu
//! class for (x_0, FF_0) and (x_1, FF_1) scattering factors
// ---------------- class coherent_mu_param ------------------
class coherent_mu_param : public Ddvec::dd_pair, public Pbase::param_base
{
public:
  Coh::anomalous_ptrs anomalous;

  coherent_mu_param( ) {}
  inline ~coherent_mu_param( ) {}

  //! gets the angular probability density
  //! \param mu the direction cosine
  //! \param interp_OK, true if the interpolation is OK
  double get_sigma( double mu, bool *interp_OK );
};

//! Class for the list of intersections of an x=const curve with an integration box
// ---------------- class coherent_hit_list ------------------
class coherent_hit_list : public Box::hit_list
{
private:
  //! The virtual functions are not used for this class
  //! \param E_in energy of incident gamma
  double get_Eout( double E_in );

  //! \param E_out energy of outgoing gamma
  //! \param Ein_hits incident energies which produce this value of E_out
  void find_bottom_hits( double E_out, std::vector< Box::Ein_Eta_Hit > *Ein_hits );

  //! \param E_out energy of outgoing gamma
  //! \param Ein_hits incident energies which produce this value of E_out
  void find_top_hits( double E_out, std::vector< Box::Ein_Eta_Hit > *Ein_hits );

public:
  coherent_hit_list( ) {}
  ~coherent_hit_list( ) {}

  //! Finds intersections of the curve x = const with the E'-mu box.
  //! \param x ( E_in/(c*h) )*\sqrt{ ( 1 - \mu )/2}
  //! \param E_in_left the lower incident gamma energy
  //! \param E_in_right the higher incident gamma energy
  void hit_box( double x, double E_in_left, double E_in_right );
};

//! Class for coherent scattering data
//--------------- class coherent ----------------
class coherent : public Xvec::x_vector
{
private:
  //! minimum common energy for anomalous data
  double anomalous_min_Ein;

  //! Extrapolates the scattering factor to 0 and anomalous to max_Ein
  //! \param max_Ein, the top of the highest incident energy bin
  void extrapolate_data( double max_Ein );

  //! Adds the result of one integration
  //! \param transfer: the entries in the transfer matrix get updated
  //! \param Ein_count the row of transfer to update
  //! \param Eout_count the column of transfer to update
  //! \param Ein_param data for quadrature over incident energy
  void update_T( Trf::T_matrix &transfer, int Ein_count, int Eout_count,
    Coh::coherent_Ein_param *Ein_param );

  //! Calculates the cross section
  //! \param mu_quad_rule the rule of quadrature over direction cosine
  //! \param sigma the computed cross section data
  void get_xsec( Qmeth::Quadrature_Rule mu_quad_rule, Ddvec::dd_vector& sigma );

public:
  //! scattering factor in the file
  Xvec::x_vector file_data;

  //! real anomalous factor
  Ddvec::dd_vector realAnomalous;

  //! imaginary anomalous factor
  Ddvec::dd_vector imaginaryAnomalous;

  coherent( ) {}

  ~coherent( ) {}

  //! Calculates the cross section and transfer matrix
  //! \param transfer: the entries in the transfer matrix get updated
  //! \param xsec the cross section data
  void get_T( Trf::T_matrix& transfer, Ddvec::dd_vector& xsec );

  //! Climbs up the outgoing energy bins
  //! \param transfer: the entries in the transfer matrix get updated
  //! \param Ein_count the row of transfer to update
  //! \param Ein_param data for quadrature over incident energy
  void Eout_ladder( Trf::T_matrix& transfer, int Ein_count,
     Coh::coherent_Ein_param *Ein_param );

  //! Sets the range of incident energies for this intergration
  //! \param Ein_param data for quadrature over incident energy
  void set_Ein_range( Coh::coherent_Ein_param *Ein_param );

  //! Integrates over one x-E box; loop over the (E_in, mu) region
  //! \param transfer: the entries in the transfer matrix get updated
  //! \param Ein_count the row of transfer to update
  //! \param Eout_count the column of transfer to update
  //! \param Ein_param data for quadrature over incident energy
  void one_Ebox( Trf::T_matrix& transfer, int Ein_count, int Eout_count,
    Coh::coherent_Ein_param *Ein_param );

};


//! Class for parameters for the 2-d quadrature of the coherent form factor
// ---------------- class coherent_Ein_param ------------------
class coherent_Ein_param : public Pbase::param_base
{
private:

public:
  // parameters for quadrature over mu
  Coh::coherent_mu_param  mu_param;

  // where we are in the x-scattering function data
  Xvec::x_vector::const_iterator x_ptr;
  Xvec::x_vector::const_iterator next_x;

  // the current integration interval
  double left_E;
  double right_E;

  // pointers to the outgoing energy boundaries
  std::vector< double >::const_iterator Eout_ptr;
  std::vector< double >::const_iterator next_Eout;

  bool use_mu_minus1;  // Is the lower limit mu = -1?
  bool use_mu1;        // Is the upper limit mu = 1?
  Qmeth::Quadrature_Rule mu_quad_rule; // quadrature rule for outgoing cosine

  // where the lower and upper x values hit the quadrature box
  Coh::coherent_hit_list lower_hits;
  Coh::coherent_hit_list upper_hits;
 
  // number of 2-d quadratures
  long int quad_count;
  // number of calls to coherent_F::Ein_F or coherent_F::sigma_F
  long int Ein_F_count;
  long int mu_F_count;  // number of calls to coherent_F::mu_F

  coherent_Ein_param( ): quad_count( 0 ), Ein_F_count( 0 ), mu_F_count( 0 ) {}
  ~coherent_Ein_param( ) {}

  //! Finds the lower mu, given E_in (uses the larger x-value)
  //! \param E_in energy of incident gamma
  inline double bottom_mu( double E_in )
  { return x_vector_F::get_mu_from_x( next_x->x, E_in ); }

  //! Finds the upper mu, given E_in (uses the smaller x-value)
  //! \param E_in energy of incident gamma
  inline double top_mu( double E_in )
  { return x_vector_F::get_mu_from_x( x_ptr->x, E_in ); }

  //! Sets up the loop over cross section and anomalous data
  //! Returns false if the data misses the energy bin.
  bool start_sigma( );

  //! Gets the current range of incident energies
  void get_range( );

  //! Increments the data for the next range of incident energies
  bool next_range( );
};

} // end of namespace Coh

namespace coherent_F
{
  // ----------------- functions to integrate --------------------
  //! Function for computing the cross section
  //! Returns false if there are problems
  //! \param mu the center-of-mass direction cosine of the outgoing particle
  //! \param FF_param the function parameters
  //! \param value the value of the integrand
  bool sigma_F( double mu, Qparam::QuadParamBase *FF_param,
		Coef::coef_vector *value );

  //! Function for the 1-d quadrature over mu
  //! Returns false if there are problems
  //! \param mu the center-of-mass direction cosine of the outgoing particle
  //! \param FF_param the function parameters
  //! \param value the value of the integrand
  bool mu_F( double mu, Qparam::QuadParamBase *FF_param,
	     Coef::coef_vector *value );

  //! Function for computing the cross section with singularity sqrt( flip_mu )
  //! Returns false if there are problems
  //! \param flip_mu = 1 - mu
  //! \param FF_param the function parameters
  //! \param value the value of the integrand
  bool flip_sigma_F( double flip_mu, Qparam::QuadParamBase *FF_param,
		     Coef::coef_vector *value );

  //! Function for the 1-d quadrature over mu with singularity sqrt( flip_mu )
  //! Returns false if there are problems
  //! \param flip_mu = 1 - mu
  //! \param FF_param the function parameters
  //! \param value the value of the integrand
  bool flip_mu_F( double flip_mu, Qparam::QuadParamBase *FF_param,
		  Coef::coef_vector *value );

  //! Function for the 2-d quadrature over E_in
  //! Returns false if there are problems
  //! \param E_in the energy of the incident gamma
  //! \param FF_param the function parameters
  //! \param value the value of the integrand
  bool Ein_F( double E_in, Qparam::QuadParamBase *coherent_param,
	      Coef::coef_vector *value );
}

#endif
