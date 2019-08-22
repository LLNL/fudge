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

using namespace std;

class coherent_Ein_param;  // forward declaration

//! Class to hold pointers to real and imaginary anomalous scattering factors
// ---------------- class anomalous_ptrs ------------------
class anomaolous_ptrs
{
public:

  // pointers to real anomalous factor data
  dd_vector::const_iterator Re_anomalous_ptr;
  dd_vector::const_iterator next_Re_anomalous;
  dd_vector::const_iterator start_Re_anomalous;
  dd_vector::const_iterator end_Re_anomalous;
  dd_vector::const_iterator first_ladder_Re_anomalous;
  dd_vector::const_iterator last_ladder_Re_anomalous;

  // pointers to imaginary anomalous factor data
  dd_vector::const_iterator Im_anomalous_ptr;
  dd_vector::const_iterator next_Im_anomalous;
  dd_vector::const_iterator start_Im_anomalous;
  dd_vector::const_iterator end_Im_anomalous;
  dd_vector::const_iterator first_ladder_Im_anomalous;
  dd_vector::const_iterator last_ladder_Im_anomalous;

  anomaolous_ptrs( ) {}
  ~anomaolous_ptrs( ) {}

  //! Gets the real anomalous factor at this energy
  //! \param E_in the incident energy
  double get_Re_anomalous( double E_in );

  //! Gets the imaginary anomalous factor at this energy
  //! \param E_in the incident energy
  double get_Im_anomalous( double E_in );

  //! Initializes the pointers to real anomalous factor data
  //! \param Re_anomalous the real anomalous factor data
  void init_Re_anomalous( dd_vector &Re_anomalous );

  //! Initializes the pointers to imaginary anomalous factor data
  //! \param Im_anomalous the imaginary anomalous factor data
  void init_Im_anomalous( dd_vector &Im_anomalous );

  //! Initializes the pointers to real anomalous factor data
  //! \param Re_anomalous the real anomalous factor data, count down
  void init_Re_down( dd_vector &Re_anomalous );

  //! Initializes the pointers to imaginary anomalous factor data
  //! \param Im_anomalous the imaginary anomalous factor data, count down
  void init_Im_down( dd_vector &Im_anomalous );

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
class coherent_mu_param : public dd_pair, public param_base
{
public:
  anomaolous_ptrs anomalous;

  coherent_mu_param( ) {}
  inline ~coherent_mu_param( ) {}

  //! gets the angular probability density
  //! \param mu the direction cosine
  double get_sigma( double mu );
};
//! Class for the list of intersections of an x=const curve with an integration box
// ---------------- class coherent_hit_list ------------------
class coherent_hit_list : public hit_list
{
private:
  //! The virtual functions are not used for this class
  //! \param E_in energy of incident gamma
  double get_Eout( double E_in );

  //! \param E_out energy of outgoing gamma
  //! \param Ein_hits incident energies which produce this value of E_out
  void find_bottom_hits( double E_out, vector< Ein_Eta_Hit > *Ein_hits );

  //! \param E_out energy of outgoing gamma
  //! \param Ein_hits incident energies which produce this value of E_out
  void find_top_hits( double E_out, vector< Ein_Eta_Hit > *Ein_hits );

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
class coherent : public x_vector
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
  void update_T( T_matrix &transfer, int Ein_count, int Eout_count,
    coherent_Ein_param *Ein_param );

  //! Calculates the cross section
  //! \param mu_quad_method the method of quadrature over direction cosine
  //! \param sigma the computed cross section data
  void get_xsec( Quadrature_Method mu_quad_method, dd_vector& sigma );

public:
  //! scattering factor in the file
  x_vector file_data;

  //! real anomalous factor
  dd_vector realAnomalous;

  //! imaginary anomalous factor
  dd_vector imaginaryAnomalous;

  coherent( ) {}

  ~coherent( ) {}

  //! Calculates the cross section and transfer matrix
  //! \param transfer: the entries in the transfer matrix get updated
  //! \param xsec the cross section data
  void get_T( T_matrix& transfer, dd_vector& xsec );

  //! Climbs up the outgoing energy bins
  //! \param transfer: the entries in the transfer matrix get updated
  //! \param Ein_count the row of transfer to update
  //! \param Ein_param data for quadrature over incident energy
  void Eout_ladder( T_matrix& transfer, int Ein_count,
     coherent_Ein_param *Ein_param );

  //! Sets the range of incident energies for this intergration
  //! \param Ein_param data for quadrature over incident energy
  void set_Ein_range( coherent_Ein_param *Ein_param );

  //! Integrates over one x-E box; loop over the (E_in, mu) region
  //! \param transfer: the entries in the transfer matrix get updated
  //! \param Ein_count the row of transfer to update
  //! \param Eout_count the column of transfer to update
  //! \param Ein_param data for quadrature over incident energy
  void one_Ebox( T_matrix& transfer, int Ein_count, int Eout_count,
    coherent_Ein_param *Ein_param );

};


//! Class for parameters for the 2-d quadrature of the coherent form factor
// ---------------- class coherent_Ein_param ------------------
class coherent_Ein_param : public param_base
{
private:

public:
  // parameters for quadrature over mu
  coherent_mu_param  mu_param;

  // where we are in the x-scattering function data
  x_vector::const_iterator x_ptr;
  x_vector::const_iterator next_x;

  // the current integration interval
  double left_E;
  double right_E;

  // pointers to the outgoing energy boundaries
  vector< double >::const_iterator Eout_ptr;
  vector< double >::const_iterator next_Eout;

  bool use_mu_minus1;  // Is the lower limit mu = -1?
  bool use_mu1;        // Is the upper limit mu = 1?
  Quadrature_Method mu_quad_method; // quadrature method for outgoing cosine

  // where the lower and upper x values hit the quadrature box
  coherent_hit_list lower_hits;
  coherent_hit_list upper_hits;
 
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
  void start_sigma( );

  //! Gets the current range of incident energies
  void get_range( );

  //! Increments the data for the next range of incident energies
  bool next_range( );
};

namespace coherent_F
{
  // ----------------- functions to integrate --------------------
  //! Function for computing the cross section
  //! \param mu the center-of-mass direction cosine of the outgoing particle
  //! \param FF_param the function parameters
  //! \param value the value of the integrand
  void sigma_F( double mu, QuadParamBase *FF_param, coef_vector *value );

  //! Function for the 1-d quadrature over mu
  //! \param mu the center-of-mass direction cosine of the outgoing particle
  //! \param FF_param the function parameters
  //! \param value the value of the integrand
  void mu_F( double mu, QuadParamBase *FF_param, coef_vector *value );

  //! Function for the 2-d quadrature over E_in
  //! \param E_in the energy of the incident gamma
  //! \param FF_param the function parameters
  //! \param value the value of the integrand
  void Ein_F( double E_in, QuadParamBase *coherent_param, coef_vector *value );
}

#endif
