/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2009-09-15  (Tue., Sept. 15, 2009) $
 * $Author: hedstrom $
 * $Id: kalbach.hpp 1 2009-09-15 hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
//! Defines the class used for the Kalbach model for energy-angle probability density

#ifndef KALBACH_CLASS
#define KALBACH_CLASS

#include "kalbach_data.hpp"
#include "Ecm_Elab_geom.hpp"
#include "transfer.hpp"

namespace Kbach
{
  class kalbach_Ein_param;  // forward declaration

//! Class for parameters for the 2-d quadrature over cm cosine and Eout_cm
// ---------------- class kalbach_Ecm_param ------------------
class kalbach_Ecm_param : public Egeom::Ecm_Elab_Ecm_param
{
public:
  // the data entries for this incident energy
  Kdata::kalbach_data Ein_data;  // the current Kalbach data
  Kdata::Kalbach_a *kalbach_a;
  Qmeth::Quadrature_Rule mu_quad_rule;  // quadrature rule for outgoing cosine
  long int mu_F_count;  // number of calls to Kalbach_F::mu_F

  inline kalbach_Ecm_param( ): mu_F_count ( 0 ) {}
  inline ~kalbach_Ecm_param( ) {}

  //! Sets up the data for this incident energy
  //! \param E_in energy of the incident particle
  //! \param Ein0_data Kalbach data at a lower incident energy
  //! \param Ein1_data Kalbach data at a higher incident energy
  //! \param Eoutmin minimum outgoing energy in the lab frame
  //! \param Eoutmax maximum outgoing energy in the lab frame
  //! \param kalbach_A data for the Kalbach a parameter
  void setup( double E_in, const Kdata::kalbach_data &Ein0_data,
              const Kdata::kalbach_data &Ein1_data,
              double Eoutmin, double Eoutmax, Kdata::Kalbach_a *kalbach_A );

};

//! Class for parameters for the 1-d quadrature over cm cosine
// ---------------- class kalbach_mu_param ------------------
class kalbach_mu_param : public Egeom::Ecm_Elab_mu_param
{
public:
  double r;
  double a;
  // Is the incident particle a gamma?
  bool gamma_in;

  inline kalbach_mu_param( ) {}
  inline ~kalbach_mu_param( ) {}

  // Sets up the data for this incident energy
  //! \param E_in energy of the incident particle
  //! \param Eoutcm center-of-mass energy of the outgoing particle
  //! \param Ecm_param the computed paramters for quadrature over outgoing energy
  void setup( double Ein, double Eoutcm, const Kbach::kalbach_Ecm_param& Ecm_param );
};

//! class for Kalbach data at one incident energy
// ---------------- class Kalbach_one_Ein ------------------
class Kalbach_one_Ein : public Ddvec::dd_vector
{
private:

public:
  //! the values of (E_out, r)
  Ddvec::dd_vector Eout_r;

  Ddvec::unit_base_map unit_base;  // unit_base for this incident energy

  Kalbach_one_Ein( ) {}
  ~Kalbach_one_Ein( ){}

  //! Reads the Kalbach energy probability density for one incident energy
  //! \param infile input file
  //! \param num_Eout number of outgoing energies for this data
  void read_probability( Dpar::data_parser &input_file, int num_Eout );

  //! Reads the Kalbach r parameters for one incident energy
  //! \param infile input file
  //! \param num_Eout number of outgoing energies for this data
  void read_r( Dpar::data_parser &input_file, int num_Eout );

   // Truncates histogram data at the maximum energy
  //! \param EMax the maximum incident energy for the transfer matrix
  void chop_histogram( double EMax );
};

//! class for Kalbach data
// ---------------- class Kalbach ------------------
class Kalbach : public std::list< Kbach::Kalbach_one_Ein >
{
private:
  Maps::map_cm_lab map;

  //! The smallest incident energy for cross section, multiplicity,
  //! model weight, flux weight, and energy groups
  double E_first;

  int first_Ein; // index of the left-hand end of the first significant energy bin
  int last_Ein;  // index of the right-hand end of the last significant energy bin

  //! Checks the data for consistency
  void check_data( );

  //! Gets the range of nontrivial incident energy bins; computes first_Ein and last_Ein
  //! returns true if the threshold is too high for the energy bins
  //! \param sigma_ the cross section data
  //! \param mult- the outgoing particle multiplicity data
  //! \param weight- the weighting to apply to the transfer matrix entries
  //! \param e_flux_ the initial approximation to apply to the particle flux
  //! \param Ein_groups the boundaries of the incident energy groups
  bool get_Ein_range( const Ddvec::dd_vector& sigma_, const Ddvec::dd_vector& mult_,
		      const Ddvec::dd_vector& weight_,
		      const Lgdata::Flux_List& e_flux_,
		      const Egp::Energy_groups& Ein_groups );

  //! Adds to the transfer matrix for a pair of incident energies.
  //! \param transfer the entries in the transfer matrix get updated
  //! \param Ein_param the quadrature parameters
  void cm_Eout_ladder( Trf::T_matrix& transfer
		       , Kbach::kalbach_Ein_param *Ein_param );

  //! Loops over the E_out bins for a given cm_Eout bin
  //! \param transfer: the entries in the transfer matrix get updated
  //! \param Ein_param the quadrature parameters
  void lab_Eout_ladder( Trf::T_matrix& transfer,
			Kbach::kalbach_Ein_param &Ein_param );

  //! Initializes the quadrature parameters
  //! \param sigma the cross section data
  //! \param mult the outgoing particle multiplicity data
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param e_flux the initial approximation to apply to the particle flux
  //! \param Ein_groups the boundaries of the incident energy groups
  void setup_param( const Ddvec::dd_vector& sigma, const Ddvec::dd_vector& mult,
    const Ddvec::dd_vector& weight, const Lgdata::Flux_List& e_flux,
    const Egp::Energy_groups& Ein_groups );

  //! Initializes the quadrature parameters
  //! \param Ein_param the quadrature parameters
  void setup_data( Kbach::kalbach_Ein_param *Ein_param );

  //! Sets the range of incident energies for this intergration
  //! \param Ein_bin the incident energy bin
  //! \param Ein_param the quadrature parameters
  void set_Ein_range( int Ein_bin, Kbach::kalbach_Ein_param &Ein_param );

  //! Does the integration for one E-E' box between 2 eta = const hyperbolas
  //! \param transfer: the entries in the transfer matrix get updated
  //! \param Eout_count identifies the matrix entry to update
  //! \param Ein_param the quadrature parameters
  void one_Ebox( Trf::T_matrix& transfer, int Eout_count,
		 Kbach::kalbach_Ein_param &Ein_param );

  //! Go to the next pair of incident energies.  Returns "true" when finished.
  //! \param  E_in the next incident energy
  //! \param Ein_param the quadrature parameters
  bool next_ladder( double E_in, Kbach::kalbach_Ein_param *Ein_param );

  //! Adds to an element of transfer the integral between the intersections of 2 Eout_cm = const arcs with the Eout_lab box
  //! \param transfer: the entries in the transfer matrix get updated
  //! \param Eout_count identifies the matrix entry to update
  //! \param Ein_param the quadrature parameters
  void update_T( Trf::T_matrix &transfer, int Eout_count,
		 Kbach::kalbach_Ein_param &Ein_param );

  //! Starts one staircase of the Eout_cm histogram
  void start_Eout_cm( );

public:
  Kdata::Kalbach_a kalbach_a;
  Terp::two_d_interp Ein_interp;
  Terp::Interp_Type Eout_interp;
  Terp::two_d_interp Ein_r_interp;
  Terp::Interp_Type r_interp;

  Kalbach( );
  ~Kalbach( ) {}

  //! Reads the Kalbach energy probability density
  //! \param infile input file
  //! \param num_Ein number of incident energies for this reaction
  void read_probability( Dpar::data_parser &input_file, int num_Ein );

  //! Reads the Kalbach r parameters
  //! \param infile input file
  //! \param num_Ein number of incident energies for this reaction
  void read_r( Dpar::data_parser &input_file, int num_Ein );

  //! Calculates the transfer matrix for this particle
  //! \param sigma the cross section data
  //! \param multiple the outgoing particle multiplicity data
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param transfer the transfer matrix
  void get_T( const Ddvec::dd_vector& sigma, const Ddvec::dd_vector& multiple, 
	      const Ddvec::dd_vector& weight, Trf::T_matrix& transfer );
};

//! Class for parameters for the 3-d quadrature over Ein, cm cosine, and Eout_cm
// ---------------- class kalbach_Ein_param ------------------
class kalbach_Ein_param : public Egeom::Ecm_Elab_Ein_param
{
private:
  //! Extrapolated data used with linlin interpolation of incident energy
  Kbach::Kalbach_one_Ein low_linlin;  // (Eout, probability) at lower incident energy
  Kbach::Kalbach_one_Ein high_linlin;  // (Eout, probability) at higher incident energy

  //! Interpolates (Eout_cm, probability) data to the lower common Eout_cm value
  //! \param lower_Eout the lower end of the Eout_cm range
  void common_low_Eoutcm( double lower_Eout );

  //! Interpolates (Eout_cm, probability) data to the higher common Eout_cm value
  //! \param higher_Eout the upper end of the Eout_cm range
  void common_high_Eoutcm( double higher_Eout );

  //! Sets up the data for unit-base interpolation in incident energy
  void setup_Ein_ubase( );

  //! Sets up the data for linlin interpolation in incident energy
  void setup_Ein_linlin( );

public:
  // parameters for integration over center-of-mass energy
  Kbach::kalbach_Ecm_param Ecm_params;

  // (Eout_cm, probability) data interpolated to common Eout_cm values
  Kdata::kalbach_data Ein0_data;  // Kalbach data for lower incident energy
  Kdata::kalbach_data Ein1_data;  // Kalbach data for higher incident energy

  Kdata::Kalbach_a *kalbach_a;
  Qmeth::Quadrature_Rule Eout_quad_rule;  // quadrature rule for outgoing energy
  Qmeth::Quadrature_Rule mu_quad_rule;  // quadrature rule for outgoing cosine

  // where we are in the list
  Kbach::Kalbach::iterator this_Ein;
  Kbach::Kalbach::iterator next_Ein;

  // pointers to the data entries
  Kbach::Kalbach_one_Ein::const_iterator left_ptr;  //(lower Eoutcm, probability), low Ein
  Kbach::Kalbach_one_Ein::const_iterator next_left_ptr; //(higher Eoutcm, probability), low Ein
  Kbach::Kalbach_one_Ein::const_iterator last_left_ptr; // end of (Eoutcm, probability),low Ein
  Ddvec::dd_vector::const_iterator left_r_ptr; //(lower Eoutcm, r), low Ein
  Ddvec::dd_vector::const_iterator next_left_r_ptr; //(higher Eoutcm, r), low Ein
  Ddvec::dd_vector::const_iterator last_left_r_ptr; //end of (Eoutcm, r), low Ein
  Kbach::Kalbach_one_Ein::const_iterator right_ptr;  //(lower Eoutcm, probability), high Ein
  Kbach::Kalbach_one_Ein::const_iterator next_right_ptr; //(higher Eoutcm, probability), high Ein
  Kbach::Kalbach_one_Ein::const_iterator last_right_ptr; //end of (Eoutcm, probability), high Ein
  Ddvec::dd_vector::const_iterator right_r_ptr; //(lower Eoutcm, r), high Ein
  Ddvec::dd_vector::const_iterator next_right_r_ptr; //(higher Eoutcm, r), high Ein
  Ddvec::dd_vector::const_iterator last_right_r_ptr; //end of (Eoutcm, r), high Ein

  Vhit::Vcm_Vlab_hit_list upper_hits;
  Vhit::Vcm_Vlab_hit_list lower_hits;

  long int quad_count;  // number of 3-d quadratures
  long int Ein_F_count;  // number of calls to Kalbach_F::Ein_F
  long int Ecm_F_count;  // number of calls to Kalbach_F::Ecm_F
  long int mu_F_count;  // number of calls to Kalbach_F::mu_F
  int Vcm_hit_count;     // number of calls to quad_F::integrate in Kalbach_F::Ein_F

  kalbach_Ein_param( ): quad_count(0), Ein_F_count(0), Ecm_F_count( 0 ), \
			mu_F_count(0) {}
  ~kalbach_Ein_param( ) {}

  //! Sets up the parameters for integration over center-of-mass outgoing energy
  //! \param E_in energy of the incident particle
  void set_Ecm_param( double E_in );

  //! Initializes the pointers to the angular probabilities for this E_in
  void start_Eout_cm( );

  //! Sets up the next interval of Eout_cm values
  //! Returns true when finished
  bool next_Eoutcm( );
};
} // end of namespace Kbach

// ************* functions to integrate ******************
namespace Kalbach_F
{
  // ---------------- Kalbach_F::mu_F ------------------
  //! Function for the 1-d quadrature over cm cosine
  //! Returns true if the interpolation is OK.
  //! \param mu the center-of-mass direction cosine of the outgoing particle
  //! \param mu_quad_param the function parameters
  //! \param value the value of the integrand, a set of Legendre coefficients
  bool mu_F( double mu, Qparam::QuadParamBase *mu_quad_param,
	     Coef::coef_vector *value );

  // ---------------- Kalbach_F::Ecm_F ------------------
  //! Function for the 2-d quadrature over cm cosine and outgoing energy
  //! Function for the 1-d quadrature over cm cosine
  //! \param Eout_cm the center-of-mass energy of the outgoing particle
  //! \param Ecm_quad_param the function parameters
  //! \param value the value of the integrand, a set of Legendre coefficients
  bool Ecm_F( double Eout_cm, Qparam::QuadParamBase *Ecm_quad_param,
	      Coef::coef_vector *value );

  // ---------------- Kalbach_F::Ein_F ------------------
  //! Function for the 3-d quadrature over incident energy, cm cosine, and outgoing energy
  //! Function for the 1-d quadrature over cm cosine
  //! \param E_in the energy of the incident particle
  //! \param quad_param the function parameters
  //! \param value the value of the integrand, a set of Legendre coefficients
  bool Ein_F( double E_in, Qparam::QuadParamBase *quad_param,
	      Coef::coef_vector *value );
}

#endif
