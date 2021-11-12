/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: energy_dist.cpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// implementation of the base classes used to handle energy distributions

#include <cmath>

#include "energy_dist_base.hpp"
#include "messaging.hpp"
#include "global_params.hpp"


// ************* class Ebase::Eprob_vector *****************
// ----------- Ebase::Eprob_vector::unit_base --------------
// Transform one energy distribuiton to unit base
void Ebase::Eprob_vector::unit_base( int L_order )
{
  bool Renorm = ( L_order == 0 );
  Ddvec::dd_vector::unit_base( Renorm, &ubase_map );
}
// ----------- Ebase::Eprob_vector::form_cum_prob --------------
// Forms the list of cumulative probabilities
void Ebase::Eprob_vector::form_cum_prob( )
{
  // copy the data
  cum_prob.Eout_interp = interp_type;
  for( Ebase::Eprob_vector::const_iterator Eout_ptr = begin( );
       Eout_ptr != end( ); ++Eout_ptr )
  {
    Cum::cumulative_prob_list::iterator cum_prob_ptr = cum_prob.insert(
      cum_prob.end( ), Cum::cumulative_prob_entry( ) );
    cum_prob_ptr->E_out = Eout_ptr->x;
    cum_prob_ptr->Prob = Eout_ptr->y;
  }
  // now form the slopes and cumulative probabilities
  if( interp_type == Terp::HISTOGRAM )
  {
    cum_prob.get_cum_prob_flat( );
  }
  else // lin-lin
  {
    cum_prob.get_cum_prob_linlin( );
  }
}

// ************* class Ebase::energy_dist_base *****************
// ----------- Ebase::energy_dist_base::print --------------
void Ebase::energy_dist_base::print( )
{
  for( Ebase::energy_dist_base::iterator energy_ptr = begin( );
       energy_ptr != end( ); ++energy_ptr )
  {
    energy_ptr->print( );
  }
}
// ----------- Ebase::energy_dist_base::unit_base --------------
// Transform the outgoing energy distribution to the interval [0, 1]
void Ebase::energy_dist_base::unit_base( int L_order )
{
  // do the mapping for each incident energy
  for( Ebase::energy_dist_base::iterator Ein_ptr = begin( );
       Ein_ptr != end( ); ++Ein_ptr )
  {
    Ein_ptr->unit_base( L_order );
  }
}
