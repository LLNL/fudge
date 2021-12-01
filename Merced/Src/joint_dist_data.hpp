/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 601 $
 * $Date: 2017-12-18 $
 * $Author: hedstrom $
 * $Id: joint_dist_data.hpp 601 2017-12-18Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// header for the classes for data used on joint energy-angle distributions

#ifndef JOINT_DIST_DATA
#define JOINT_DIST_DATA

#include "dd_vector.hpp"

namespace Jdata
{
// ----------- class E_mu_P_data  -----------------
//! Class for one set of energy-mu-probability data
class E_mu_P_data
{
public:
  double UB_mu;   // the unit-base direction cosine of the outgoing particle
  double phys_mu;   // the physical direction cosine of the outgoing particle
  double UB_Eout; // the unit-base energy of the outgoing particle
  double phys_Eout; // the physical energy of the outgoing particle
  double mu_Prob;  // the cosine probability density of the outgoing particle
  double Eout_Prob;  // the energy probability density of the outgoing particle

  //! Default constructor
  inline E_mu_P_data( )
  {}

  //! Default destructor
  inline ~E_mu_P_data( )
  {}

  //! Interpolates with respect to unit-base outgoing energy
  //! Returns true if the interpolation is OK
  //! \param mid_UB_Eout an intermediate unit-base outgoing energy
  //! \param next_data Jdata::E_mu_P_data at a higher unit-base outgoing energy
  //! \param mid_data Jdata::E_mu_P_data at unit-base outgoing energy mid_UB_Eout
  bool UB_Eout_interp( double mid_UB_Eout, const Jdata::E_mu_P_data &next_data,
		       Jdata::E_mu_P_data *mid_data ) const;

  //! Interpolates histogram data with respect to unit-base outgoing energy
  //! \param mid_UB_Eout an intermediate unit-base outgoing energy
  //! \param next_data Jdata::E_mu_P_data at a higher unit-base outgoing energy
  //! \param mid_data Jdata::E_mu_P_data at unit-base outgoing energy mid_UB_Eout
  void UB_Eout_histogram( double mid_UB_Eout, const Jdata::E_mu_P_data &next_data,
		       Jdata::E_mu_P_data *mid_data ) const;

  //! Interpolates with respect to direction cosine
  //! Returns true if the interpolation is OK
  //! \param mid_mu an intermediate direction cosine
  //! \param next_data Jdata::E_mu_P_data at a higher direction cosine
  //! \param mid_data Jdata::E_mu_P_data at direction cosine
  bool mu_interp( double mid_mu, const Jdata::E_mu_P_data &next_data,
		  Jdata::E_mu_P_data *mid_data ) const;

  //! Interpolates with respect to incident energy
  //! Returns true if the interpolation is OK
  //! \param alpha the weight for next_data
  //! \param next_data Jdata::E_mu_P_data at a higher incident energy
  //! \param mid_data Jdata::E_mu_P_data at an an intermediate incident energy
  bool Ein_interp( double alpha, const Jdata::E_mu_P_data &next_data,
		   Jdata::E_mu_P_data *mid_data ) const;

  //! Copies the data
  //! \param to_copy Jdata::E_mu_P_data to copy
  void copy( const Jdata::E_mu_P_data &to_copy );
};

// ----------- class current_data -----------------
//! class for Jdata::E_mu_P_data at current values of Ein, mu, and Eout
class current_data
{
private:
  double E_in;  // energy of the incident particle

public:
  //! ( mu, outgoing energy, probability density ) data
  Jdata::E_mu_P_data mu0_Eout0;  // lower mu, lower Eout
  Jdata::E_mu_P_data mu0_Eout1;  // lower mu, higher Eout
  Jdata::E_mu_P_data mu1_Eout0;  // higher mu, lower Eout
  Jdata::E_mu_P_data mu1_Eout1;  // higher mu, higher Eout

  //! range of outgoing energies for the unit-base maps
  Ddvec::unit_base_map mu0_ubase_map;  // lower mu
  Ddvec::unit_base_map mu1_ubase_map;  // higher mu

  //! Constructor
  inline current_data( )
  {}

  //! Default destructor
  inline ~current_data( )
  {}

  //! Sets the incident energy
  //! \param Ein the energy of the incident particle
  inline void set_E_in( double Ein )
  { E_in = Ein; }

  //! Returns the energy of the incident particle
  inline double get_E_in( ) const
  { return E_in; }

  //! Interpolate in incident energy, returns mid_data
  //! returns true if the interpolation is OK
  //! \param mid_Ein an intermediate incident energy
  //! \param next_data (E_out, probability) at a higher incident energy
  //! \param mid_data (E_out, probability) at incident energy mid_Ein
  bool Ein_interpolate( double mid_Ein, const current_data &next_data,
		    current_data *mid_data ) const;

  //! Interpolate in mu; returns Eout0_data, Eout1_data, and mid_ubase_map
  //! Returns true if the interpolation is OK
  //! \param mid_mu an intermediate direction cosine
  //! \param Eout0_data (E_out, probability) at mid_mu and lower outgoing energy
  //! \param Eout1_data (E_out, probability) at mid_mu and higher outgoing energy
  //! \param mid_ubase_map the outgoing energy range at mid_mu
  bool mu_interpolate( double mid_mu, Jdata::E_mu_P_data *Eout0_data,
		       Jdata::E_mu_P_data *Eout1_data,
		       Ddvec::unit_base_map *mid_ubase_map ) const;

};

} // end of namespace Jdata

#endif
