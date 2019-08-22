/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2011-03-07 (Mon, Mar 7, 2011) $
 * $Author: hedstrom $
 * $Id: cumulative_points.hpp 1 2011-03-07 hedstrom $
 *
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// classes used to handle interpolation by cumulative points

#ifndef CUMULATIVE_POINTS_DEF
#define CUMULATIVE_POINTS_DEF

#include <list>

#include "dd_vector.hpp"  // for Interp_Type

using namespace std;

//! Class for cumulative probability at one outgoing energy
//--------------- class cumulative_prob_entry ----------------
class cumulative_prob_entry
{
public:
  // ! the outgoing energy
  double E_out;

  //! the probability density
  double Prob;

  //! The derivative of the zero-order coefficient with respect to outgoing energy
  double slope;

  //! The cumulative probability up to this entry
  double cum_prob;

  inline cumulative_prob_entry( ) {}

  inline ~cumulative_prob_entry( ) {}

  //! Gets the probability density at incident energy E_in
  //! \param E the current outgoing energy
  double get_prob( double E );

  //! Gets the cumulative probability at incident energy E_in
  //! \param E the current outgoing energy
  double get_cum_prob( double E );

  //! Gets the energy corresponding to cumulative probability A
  //! \param A the current cumulative probability
  double get_cum_inv( double A ) const;
};

//! Class for cumulative probabilities at one incident energy
//--------------- class cumulative_prob_list ----------------
class cumulative_prob_list : public list< cumulative_prob_entry >
{
 private:
  //! The energy of the incident particle
  double E_in;

 public:
  Interp_Type Eout_interp;

  cumulative_prob_list( ) {}

  ~cumulative_prob_list( ) {}

  //! Returns the energy of the incident particle
  inline double get_E_in( ) const { return E_in; }

  //! Sets the energy of the incident particle
  //! \param E_in energy of incident particle
  inline void set_E_in( double Ein ) { E_in = Ein; }

  //! Computes the cumulative probabilities and slopes for histogram data
  void get_cum_prob_flat( );

  //! Computes the cumulative probabilities and slopes for lin-lin data
  void get_cum_prob_linlin( );

};


#endif
