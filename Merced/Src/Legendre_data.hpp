/*
* ******** merced: calculate the transfer matrix *********
* $Revision: 1 $
* $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
* $Author: hedstrom $
* $Id: Legendre_data.hpp 1 2006-02-02 03:06:56Z hedstrom $
*
* ******** merced: calculate the transfer matrix *********
*
* # <<BEGIN-copyright>>
* # <<END-copyright>>
*/

// header for the list of (incident energy, flux)
#ifndef E_FLUX
#define E_FLUX

#include <iostream>
#include <list>
#include <vector>

#include "Legendre_base.hpp"
#include "dd_vector.hpp"
#include "quad_param.hpp"
#include "coef_vector.hpp"
#include "data_parser.hpp"

namespace Lgdata
{

  // --------------- class Legendre_coefs -------------------------
//! Class to hold Legendre coefficients of flux j=0, 1, ..., order
  class Legendre_coefs : public LgBase::Legendre_base
{
private:
//! Returns true if the interpolation is OK
//! \param alpha the weight for next_flux
//! \param prev_flux Legendre coefficients at a lower energy
//! \param next_flux Legendre coefficients at a higher energy
void basic_linlin_interp( double alpha, const Legendre_coefs& prev_flux,
const Legendre_coefs& next_flux );

public:
//! Constructor
inline Legendre_coefs( )
{}

//! Constructor
//! \param Order the Legendre order
inline Legendre_coefs( int Order )
{ initialize( Order ); }

//! Destructor
inline ~Legendre_coefs( )
{}

//! Allocates space
//! \param Order the Legendre order
void initialize( int order );

//! Sets all coefficients to zero
void zero_data( );

//! Copies the Legendre coefficients
//! Resets the order, if they are inconsistent
//! \param to_copy the Legendre data to copy
void copy_coef( const Legendre_coefs& to_copy );

//! Only copies the Legendre coefficients, order not reset
//! \param to_copy the Legendre data to copy
void only_copy_coef( const Legendre_coefs& to_copy );

///! Sets the order for interpolated data
//! \param left_order the Legendre order for one set of coefficients
//! \param right_order the Legendre order for another set of coefficients
void set_max_order( int left_order, int right_order );

//! Interpolates the flux at energy E_in
//! Returns true if the interpolation is OK
//! \param E_in intermediate energy
//! \param prev_flux Legendre coefficients at a lower energy
//! \param next_flux Legendre coefficients at a higher energy
bool linlin_interp( double E_in, const Legendre_coefs& prev_flux,
const Legendre_coefs& next_flux );

//! Interpolates Legendre-coefficient data at energy E_in
//! It is required that left_flux and right_flux be at the same outgoing energy
//! Returns true if the interpolation is OK
//! \param E_in intermediate incident energy
//! \param left_Ein lower incident energy
//! \param left_flux Legendre coefficients at incident energy left_Ein
//! \param right_Ein higher incident energy
//! \param right_flux Legendre coefficients at incident energy right_Ein
bool Ein_linlin_interp( double E_in, double left_Ein, 
const Legendre_coefs& left_flux, double right_Ein,
const Legendre_coefs& right_flux );

//! Interpolates unit-base Legendre-coefficient data
//! It is required that left_flux and right_flux be at the same outgoing energy
//! Returns true if the interpolation is OK
//! \param E_in intermediate incident energy
//! \param alpha the proportionality factor
//! \param left_flux Legendre coefficients at incident energy left_Ein
//! \param right_flux Legendre coefficients at incident energy right_Ein
bool unitbase_interp( double E_in, double alpha, 
const Legendre_coefs& left_flux,
const Legendre_coefs& right_flux );

//! Interpolates the flux linearly with respect to the logarithm of the energy
//! Returns true if the interpolation is OK
//! \param E_in intermediate energy
//! \param left_flux Legendre coefficients at a lower energy
//! \param right_flux Legendre coefficients at a higher energy
bool linlog_interp( double E_in, 
const Legendre_coefs& left_flux,
const Legendre_coefs& right_flux );

    //! Addition operator
    //! \param to_add the Legendre moments to add
    Legendre_coefs& operator+=( const Coef::coef_vector &to_add );

};

  //! Class for parameters for the quadrature over mu
  // ---------------- class Lgdata::mu_param ------------------
  class mu_param: public Qparam::QuadParamBase
  {
  public:
    Ddvec::dd_vector::const_iterator left_data;
    Ddvec::dd_vector::const_iterator right_data;

    inline mu_param() {}
    inline ~mu_param() {}

    //! Evaluate by linear-linear interpolation
    //! \param mu an intermediate direction cosine
    //! \param is_OK, true if the intertpolation is OK
    double value( double mu, bool *is_OK );
  };

// ------------ class Lgdata::Legendre_data_range ------------------
//! Class to hold Legendre coefficient data for a range on outgoing energies
class Legendre_data_range
{
private:
double E_in; // the incident energy

public:
Legendre_coefs prev_data;  // Legendre data for lower outgoing energy
Legendre_coefs next_data;  // Legendre data for upper outgoing energy
Ddvec::unit_base_map ubase_map;  // unit-base map for this data

Terp::two_d_interp Ein_interp;
Terp::Interp_Type Eout_interp;

inline Legendre_data_range( ) {}
inline ~Legendre_data_range( ) {}

//! Sets the energy of the incident particle
//! \param Ein the energy of the incident particle
inline void set_E_in( double Ein )
{
E_in = Ein;
}

//! Gets the energy of the incident particle
inline double get_E_in( ) const
{
return E_in;
}

//! Sets up a new incident energy
//! \param Ein the energy of the incident particle
//! \param ubasemap the unit-base-map for this data
//! \param Eoutinterp rule for interpolation is outgoing energy
void new_Ein( double Ein, const Ddvec::unit_base_map &ubasemap, Terp::Interp_Type Eoutinterp );

//! Sets up prev_data and next_data for a given range of outgoing energies
//! \param prevdata Legendre coefficients at a lower energy
//! \param nextdata Legendre coefficients at a higher energy
//! \param Eout_min bottom of the desired energy range
//! \param Eout_max top of the desired energy range
void set_data( const Legendre_coefs &prevdata, const Legendre_coefs &nextdata,
double Eout_min, double Eout_max );

//! Does unit-base interpolation between incident energies
//! Returns true if the interpolation is OK
//! It is required that left_data and right_data be at the same outgoing energy
//! \param E_in an intermediate incident energy
//! \param left_data Legendre coefficients at a lower incident energy
//! \param right_data Legendre coefficients at a higher incident energy
bool ubase_interpolate( double E_in,
const Legendre_data_range &left_data, const Legendre_data_range &right_data );

//! Maps from physical variables to unit-base
bool to_unit_base( );

//! Used by cumulative points interpolation for intervals of length zero
//! \param dA, the cumulative probability for this interval
void short_to_unit_base( double dA );

//! Maps from unit-base to physical variables
//! Returns true if the interpolation is OK
bool un_unit_base( );


};

// --------------- class Legendre_list_base -------------------------
//! Class to hold an array Legendre coefficient data
class Legendre_list_base : public std::list< Legendre_coefs >
{
private:
  //! The energy of the incident particle
  double E_in;

  //! Finds the total probability
  double get_norm( ) const;

public:
  int order;
  Ddvec::unit_base_map ubase_map;
  Terp::two_d_interp Ein_interp; // interpolation rule for incident energy
  Terp::Interp_Type Eout_interp; // interpolation rule for outgoing energy

  //! Constructor
  inline Legendre_list_base( ): order( -1 ) {}

  //! Destructor
  inline ~Legendre_list_base( ) {}

  //! Returns the energy for the Legendre data
  inline double get_E_in( ) const { return E_in; }

  //! Sets the energy for the Legendre data
  //! \param Ein the energy for the Legendre data
  inline void set_E_in( double Ein ) { E_in = Ein; }

  //! Normalizes the total probability
  //! Returns true if the norm is positive
  //! \param truncated, true if this is data for interpolation with truncation
  bool renorm( bool truncated );

  //! Maps the data to unit base
  void to_unit_base( );

  // For debugging
  void print( );

};

// --------------- class Flux_List -------------------------
//! Class to hold Legendre coefficients of flux j=0, 1, ..., order
class Flux_List : public std::list< Legendre_coefs >
{
private:

public:
  int order;

  Terp::Interp_Type interp;

  //! Constructor
  inline Flux_List( ): order( -1 ) {}

  //! Destructor
  inline ~Flux_List( ) {}

  //! Constructs the list from the python data
  //! \param infile the input file
  //! \param num_Ein the number of incident energies
  void read_flux( Dpar::data_parser &infile, int num_Ein );

  //! Evaluates the flux at energy E_in with search starting at ptr
  //! \param E_in incident energy requested
  //! \param ptr pointer to the most recent successful search
  Legendre_coefs value( double E_in, Flux_List::const_iterator &ptr ) const;
};

// --------------- class weight_vector -------------------------
//! Class to hold Legendre flux weights for one energy bin
  class weight_vector : public LgBase::Legendre_base
{
private:

public:
  //! Constructor
  inline weight_vector( )
  {}

  //! Destructor
  inline ~weight_vector( )
  {}

  //! Adds the integrals over ( E_left, E_right )
  //! \param e_flux the pairs ( incident energy, Legendre coefficients ) to integrate
  //! \param this_flux pointer to the most recent successful search
  //! \param E_left lower limit of integration
  //! \param E_right upper limit of integration
  void increment( Flux_List &e_flux, 
    Flux_List::const_iterator this_flux, double E_left, double E_right );

  //! Takes the reciprocals
  void invert( );
};
// --------------- class weight_list -------------------------
//! Class to hold all of the Legendre flux weights
class weight_list : public std::list<  weight_vector >
{
private:

public:
  //! Constructor
  inline weight_list( )
  { }

  //! Destructor
  inline ~weight_list( )
  { }

};

} // End of namespace Lgdata

namespace to_Legendre_F
{
  // ************* functions to integrate ***************
  // ---------------- mu_F ------------------
  //! Function for the quadrature over mu
  //! \param mu the direction cosine
  //! \param param the function parameters
  //! \param value the value of the integrand, a set of Legendre coefficients
  bool mu_F( double mu_in, Qparam::QuadParamBase *param, Coef::coef_vector *value );

  // ------------------ to_Legendre_F::from_table --------------
  //! Evaluates the Legendre moments
  //! Returns the number of function evaluations used in the quadrature
  //! \param mu_table pairs ( direction cosine, probability density )
  //! \param coefs, the corresponding Legendre coefficients
  int from_table( const Ddvec::dd_vector& mu_table,
		  Lgdata::Legendre_coefs *coefs );

}  // end of namespace to_Legendre_F

#endif
