/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2011-01-28 (Fri, Jan 28, 2011) $
 * $Author: hedstrom $
 * $Id: cumulative_points.cpp 1 2011-01-28 hedstrom $
 *
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// classes used to handle interpolation by cumulative points

#include <cmath>
#include <cstdlib>
#include <cfloat>
#ifdef _OPENMP
 #include <omp.h>
#endif

#include "cumulative_points.hpp"
#include "messaging.hpp"

// ************* class cumulative_prob_entry *****************
// ----------- cumulative_prob_entry::get_prob --------------
// Gets the probability density at incident energy E_in
double cumulative_prob_entry::get_prob( double E )
{
  double this_prob = Prob + slope*( E - E_out );
  return this_prob;
}
// ----------- cumulative_prob_entry::get_cum_prob --------------
// Gets the cumulative probability at incident energy E_in
double cumulative_prob_entry::get_cum_prob( double E )
{
  double dE = E - E_out;
  double this_cum = cum_prob + dE*( Prob + 0.5*slope*dE );
  return this_cum;
}
// ----------- cumulative_prob_entry::get_cum_inv --------------
// Gets the energy corresponding to cumulative probability A
double cumulative_prob_entry::get_cum_inv( double A ) const
{
  double dA = A - cum_prob;
  double this_E = ( slope == 0.0 ) ?
    E_out + dA/ Prob :
    E_out + 2.0*dA/( Prob + sqrt( Prob*Prob + 2.0*slope*dA ) );
  return this_E;
}

// ************* class cumulative_prob_list *****************
// ----------- cumulative_prob_list::get_cum_prob_flat --------------
// Computes the cumulative probabilities and slopes for histogram data
void cumulative_prob_list::get_cum_prob_flat( )
{
  if( Eout_interp != HISTOGRAM )
  {
    FatalError( "cumulative_prob_list::get_cum_prob_flat",
		"wrong interpolation type" );
  }

  // pointers to the data
  static bool zero_found = false;
  cumulative_prob_list::iterator this_link = begin( );
  cumulative_prob_list::iterator next_link = this_link;
  ++next_link;

  // set up the list of cumulative probabilities
  double sum = 0.0;
  this_link->cum_prob = 0.0;
  this_link->slope = 0.0;

  for( ; next_link != end( ); this_link = next_link, ++next_link )
  {
    // norm for histogram data
    double to_add = this_link->Prob * ( next_link->E_out - this_link->E_out );
    if( ( to_add <= 0.0 ) && ( !zero_found ) )
    {
      Warning( "cumulative_prob_list::get_cum_prob_flat",
	       "data have an interval with zero probability" );
      zero_found = true;
    }
    sum += to_add;
    next_link->cum_prob = sum;
    next_link->slope = 0.0;
  }
  // ensure that cum_prob <= 1
  for( this_link = begin( ); this_link != end( ); ++this_link )
  {
    if( this_link->cum_prob > 1.0 )
    {
      this_link->cum_prob = 1.0;
    }
  }
  this_link = end( );
  --this_link;
  this_link->cum_prob = 1.0;
}
// ----------- cumulative_prob_list::get_cum_prob_linlin --------------
// Computes the cumulative probabilities and slopes for lin-lin data
void cumulative_prob_list::get_cum_prob_linlin( )
{
  if( Eout_interp != LINLIN )
  {
    FatalError( "cumulative_prob_list::get_cum_prob_linlin",
		"wrong interpolation type" );
  }

  // pointers to the data
  static bool zero_found = false;
  cumulative_prob_list::iterator this_link = begin( );
  cumulative_prob_list::iterator next_link = this_link;
  ++next_link;
  double E0 = this_link->E_out;
  double P0 = this_link->Prob;

  // set up the list of cumulative probabilities
  double sum = 0.0;
  this_link->cum_prob = 0.0;

  for( ; next_link != end( ); this_link = next_link, ++next_link )
  {
    double E1 = next_link->E_out;
    double P1 = next_link->Prob;
    double dE = E1 - E0;
    // norm for linlin data
    double to_add = 0.5*dE*( P1 + P0 );
    if( ( to_add <= 0.0 ) && ( !zero_found ) )
    {
      Warning( "cumulative_prob_list::get_cum_prob_linlin",
	       "data have an interval with zero probability" );
      zero_found = true;
    }
    sum += to_add;
    next_link->cum_prob = sum;
    if( dE == 0.0 )
    {
      this_link->slope = 0.0;
    }
    else
    {
      this_link->slope = ( P1 - P0 )/dE;
    }
    E0 = E1;
    P0 = P1;
  }
  // ensure that cum_prob <= 1
  for( this_link = begin( ); this_link != end( ); ++this_link )
  {
    if( this_link->cum_prob > 1.0 )
    {
      this_link->cum_prob = 1.0;
    }
  }
  this_link = end( );
  --this_link;
  this_link->cum_prob = 1.0;
  this_link->slope = 0.0;
}

