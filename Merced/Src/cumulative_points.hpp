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

namespace Cum
{

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

  //! Gets the probability density at outgoing energy E_in
  //! \param E the current outgoing energy
  double get_prob( double E ) const;

  //! Gets the cumulative probability at outgoing energy E_in
  //! \param E the current outgoing energy
  double get_cum_prob( double E ) const;

  //! Gets the energy corresponding to cumulative probability A
  //! \param A the current cumulative probability
  double get_cum_inv( double A ) const;
  
  //! Sets the slope for this interval
  //! \param next_E, the outgoing energy for the next interval
  //! \param next_A, the cumulative probability for the next interval
  void set_slope( double next_E, double next_A );

  //! Interpolates between two cumulative_prob_entrys
  //! \param alpha, the weight
  //! \param this_A, the base cumulative probability for this entry
  //! \param prev_CP, cumulative_prob_entry for smaller mu
  //! \param next_CP, cumulative_prob_entry for larger mu
  void CP_interpolate( double alpha, double this_A,
		       const cumulative_prob_entry &prev_CP,
		       const cumulative_prob_entry &next_CP );

  //! Prints the entry
  void print( );
};

//! Class for cumulative probabilities at one incident energy
//--------------- class cumulative_prob_list ----------------
class cumulative_prob_list : public std::list< cumulative_prob_entry >
{
 private:
  //! Ususally the energy of the incident particle
  double tag;

 public:
  Terp::Interp_Type Eout_interp;

  cumulative_prob_list( ) {}

  ~cumulative_prob_list( ) {}

  //! Returns the tag as energy of the incident particle
  inline double get_E_in( ) const { return tag; }

  //! Sets the tag as energy of the incident particle
  //! \param E_in energy of incident particle
  inline void set_E_in( double Ein ) { tag = Ein; }

  //! Returns the tag as direction cosine
  inline double get_mu( ) const { return tag; }

  //! Sets the tag as direction cosine
  //! \param mu energy of incident particle
  inline void set_mu( double mu ) { tag = mu; }

  //! Computes the cumulative probabilities and slopes for histogram data
  void get_cum_prob_flat( );

  //! Computes the cumulative probabilities and slopes for lin-lin data
  void get_cum_prob_linlin( );

  //! Sets the slopes for all of the entries
  void set_slopes( );
};

}  // end of namespace Cum

#endif
