/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: energy_dist_base.hpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
//! Defines the base class used for energy distributions and uncorrelated energy-angle distributions.

#ifndef ENERGY_BASE_CLASS
#define ENERGY_BASE_CLASS

#include "dd_vector.hpp"
#include "box_geom.hpp"
#include "transfer.hpp"
#include "cumulative_points.hpp"
#include "param_base.hpp"

namespace Ebase
{
//! Class for one energy distribution
//--------------- class Eprob_vector ----------------
class Eprob_vector : public Ddvec::dd_vector
{
public:
  // cumulative probabilities for cumulative points interpolation
  Cum::cumulative_prob_list cum_prob;

  Ddvec::unit_base_map ubase_map;

  inline Eprob_vector( ) {}

  inline ~Eprob_vector( ) {}

  //! Transforms one energy distribuiton to unit base
  //! \param L_order the Legendre order of this data
  void unit_base( int L_order );

  //! Forms the list of cumulative probabilities
  void form_cum_prob( );
};

//! Class for energy distributions
//--------------- class energy_dist_base ----------------
class energy_dist_base : public std::list< Ebase::Eprob_vector >
{
private:

public:
  Terp::two_d_interp Ein_interp;  // interpolation between incident enrgies
  Terp::Interp_Type Eout_interp;

  inline energy_dist_base( ) {}

  inline ~energy_dist_base( ) {}

  //! Transforms the outgoing energy distributions to the interval [0, 1]
  //! \param L_order the Legendre order of this data
  void unit_base( int L_order );

  // Prints the lists for debugging
  void print( );
};

} // end of namespace Ebase

#endif
