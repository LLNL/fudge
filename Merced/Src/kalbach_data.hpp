/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2009-09-15  (Tue., Sept. 15, 2009) $
 * $Author: hedstrom $
 * $Id: kalbach_data.hpp 1 2009-09-15 hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
//! Defines the class used for data for the Kalbach model

#ifndef KALBACH_DATA_CLASS
#define KALBACH_DATA_CLASS

#include "dd_vector.hpp"
#include "mappings.hpp"

namespace Kdata
{
//! class to identify a nucleon
// ---------------- class nucleon --------------------
class nucleon
{
private:

public:
  double mass;
  double Kalbach_I;  // particle separation energy
  double Kalbach_m;
  int Kalbach_M;
  int Z;
  int A;

  nucleon( ): mass( -1.0 ), A( -1 ) {}
  ~nucleon( ) {}

  //! copies the data
  nucleon& operator=( const Kdata::nucleon& to_copy );

  //! Sets Z, A, I, M, m
  //! \param ZA the ENDL ZA number: A + 1000*Z
  void set_params( int ZA );

  //! Returns ZA
  inline int get_ZA( )
    {return A + 1000*Z;}

  //! Computes the Kalbach S_a function
  //! \param proj the identifiers for the projectile
  double get_Sa( const Kdata::nucleon &proj );

  //! Checks for proper initialization
  bool check_data( );
};

//! class for the Kalbach a coefficient
// ---------------- class Kalbach_a ------------------
class Kalbach_a
{
private:
  // Coefficients of the Kalbach polynomial
  double C_1;
  double C_2;
  double C_3;

  double projectile_S;
  double eject_S;

public:
  Kdata::nucleon target;
  Kdata::nucleon projectile;
  Kdata::nucleon eject;
  Kdata::nucleon compound;
  Kdata::nucleon residual;
  Maps::map_cm_lab *map;

  Kalbach_a( ): C_1( 0.04 ), C_2( 1.8e-6 ), C_3( 6.7e-7 ) {}
  ~Kalbach_a( ) {}

  //! Initializes the mass ratios in map
  void setup_params( );

  //! Computes the S_a and S_b functions
  void set_Sab( );

  //! Computes the Kalbach a function
  //! \param E_in the energy of the incident particle
  //! \param E_out the energy of the outgoing particle
  double get_a( double E_in, double E_out );

  //! Stores the masses
  //! \param particle_info the identities of the particles involved in the reaction
  void copy_masses( const Maps::particleInfo &particle_info );
};

//! Class for the current data entries for one incident energy
// ---------------- class kalbach_data -----------------------
class kalbach_data
{
public:
  double E_in;
  // the data entries for this incident energy
  double this_Ecm; // center-of-mass energy
  double this_f0;  // center-of-mass energy probability density
  double this_r;   // r-value
  double next_Ecm;
  double next_f0;
  double next_r;

  Terp::two_d_interp Ein_interp; // interpolation with respect to incident energy
  Terp::Interp_Type Eout_interp; // interpolation with respect to outgoing energy

  inline kalbach_data( ): Eout_interp( Terp::HISTOGRAM ) {}
  inline ~kalbach_data( ) {}

  Ddvec::unit_base_map unit_base;  // unit_base for this incident energy

  //! Does linear interpolation of this data with the next
  //! Returns true if the interpolation is OK.
  //! \param E_in the energy of the incident particle
  //! \param left_data Kalbach data at a lower incident energy
  //! \param right_data Kalbach data at a higher incident energy
  bool linlin_interp( double E_in, const Kdata::kalbach_data& left_data,
    const Kdata::kalbach_data& right_data );

  //! Does linear interpolation of this unit_base data with the next
  //! Returns true if the interpolation is OK.
  //! \param E_in the energy of the incident particle
  //! \param left_data Kalbach data at a lower incident energy
  //! \param right_data Kalbach data at a higher incident energy
  bool unit_base_interp( double E_in, const Kdata::kalbach_data& left_data,
    const Kdata::kalbach_data& right_data );

  //! Undoes the unit-base map for one outgoing energy
  //! \param E_unit the unit-base energy of the outgoing particle
  double un_unit_base( double E_unit );

  //! Undoes the unit-base map; used on interpolated data
  void un_unit_base( );

  //! Calculates r and the center-of-mass outgoing probability density
  //! \param Eoutcm center-of-mass outgoing energy
  //! \param computed probability density for this energy
  //! \param computed Kalbach r parameter for this energy
  void get_f0_r( double Eoutcm, double *Ecm_prob, double *r ) const;

};
} // end of namespace Kdata

#endif
