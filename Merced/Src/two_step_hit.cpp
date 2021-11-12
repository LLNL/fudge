/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2018-11-19  (Mon., Nov. 19, 2018) $
 * $Author: hedstrom $
 * $Id: two_step_hit.cpp 1 2018-11-19 hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
//! Implements the classes used for the geometry of the 2-step reaction

#include "two_step_hit.hpp"
#include "messaging.hpp"
#include "global_params.hpp"

// ************* class Twohit::two_step_hit_list *****************
// ----------- Twohit::two_step_hit_list::set_mucm_1 --------------
// Sets the direction cosine for step 1
void Twohit::two_step_hit_list::set_mucm_1( double mucm1 )
{
  if( use_relativistic )
  {
    relTwoStepParam.mu_cm_1 = mucm1;
  }
  else
  {
    twoStepMap.mucm_1 = mucm1;
  }
}
// ----------- Twohit::two_step_hit_list::set_mucm_2 --------------
// Sets the direction cosine for step 2
void Twohit::two_step_hit_list::set_mucm_2( double mucm2 )
{
  if( use_relativistic )
  {
    relTwoStepParam.mu_cm_2 = mucm2;
  }
  else
  {
    twoStepMap.mucm_2 = mucm2;
  }
}
// ----------- Twohit::two_step_hit_list::get_Eout --------------
// Calculates E_out from E_in, for testing the side of a box
double Twohit::two_step_hit_list::get_Eout( double E_in )
{
  double Eout;
  void *params;  // parameters for T_out_lab
  if( use_relativistic )
  {
    Eout = relTwoStepParam.relTwoStepMap.two_step_get_E_lab( E_in,
	       relTwoStepParam.mu_cm_1, relTwoStepParam.mu_cm_2 );
  }
  else
  {
    params = static_cast< void * >( &twoStepMap );
    Eout = Newtonian_F::two_step_T_out_lab( E_in, params );
  }

  return Eout;
}
// ----------- Twohit::two_step_hit_list::find_hit --------------
// Finds an intersection with the bottom or top of a box
double Twohit::two_step_hit_list::find_hit( double E_out, const Ddvec::dd_entry &pair_0,
				 const Ddvec::dd_entry &pair_1 )
{
  if( ( pair_0.y - E_out ) * ( pair_1.y - E_out ) > 0.0 )
  {
    return -1.0;
  }

  double root;
  // The computer arithmetic may do weird things
  if( ( ( pair_1.x > flip.x ) && ( pair_1.y <= pair_0.y ) ) ||
      ( ( pair_0.x < flip.x ) && ( pair_1.y >= pair_0.y ) ) )
  {
    root = flip.x;
  }
  else if( use_relativistic )
  {
    root = relTwoStepParam.find_hit( E_out, pair_0, pair_1 );
  }
  else
  {
    root = twoStepMap.find_hit( E_out, pair_0, pair_1 );
  }
  return root;
}
// ----------- Twohit::two_step_hit_list::find_bottom_hits --------------
// Finds the intersections with the bottom of a box
void Twohit::two_step_hit_list::find_bottom_hits( double E_out_lab,
					  std::vector< Box::Ein_Eta_Hit > *Ein_hits )
{
  // for new entries
  Box::Ein_Eta_Hit Ein_mucm_hit;

  // where is the minimum?
  if( flip.x <= left_Ein_Eout.x )
  {
    // E_lab is increasing with E_in
    // append this entry
    Ein_mucm_hit.E_in = find_hit( E_out_lab, left_Ein_Eout, right_Ein_Eout );
    Ein_mucm_hit.hit_edge = Box::BOTTOM_IN;
    Ein_hits->push_back( Ein_mucm_hit );
  }
  else if( flip.x >= right_Ein_Eout.x )
  {
    // E_lab is decreasing with E_in
    // append this entry
    Ein_mucm_hit.E_in = find_hit( E_out_lab, left_Ein_Eout, right_Ein_Eout );
    Ein_mucm_hit.hit_edge = Box::BOTTOM_OUT;
    Ein_hits->push_back( Ein_mucm_hit );
  }
  else if( flip.y < E_out_lab ) // omit a tangent contact
  {
    // We are here only for Q < 0
    // tolerance for short intervals
    static double etol = Global.Value( "looser_tol" );
    double slop = etol*( right_Ein_Eout.x - left_Ein_Eout.x );
    // E_lab is at first decreasing with E_in
    // append this entry
    if( flip.x >= left_Ein_Eout.x + slop )
    {
      Ein_mucm_hit.E_in = find_hit( E_out_lab, left_Ein_Eout, flip );
      Ein_mucm_hit.hit_edge = Box::BOTTOM_OUT;
      Ein_hits->push_back( Ein_mucm_hit );
    }
    // E_lab is then increasing with E_in
    // append this entry
    if( flip.x <= right_Ein_Eout.x - slop )
    {
      Ein_mucm_hit.E_in = find_hit( E_out_lab, flip, right_Ein_Eout );
      Ein_mucm_hit.hit_edge = Box::BOTTOM_IN;
      Ein_hits->push_back( Ein_mucm_hit );
    }
  }
}
// ----------- Twohit::two_step_hit_list::find_top_hits --------------
// Finds the intersections with the top of a box
void Twohit::two_step_hit_list::find_top_hits( double E_out_lab,
  std::vector< Box::Ein_Eta_Hit > *Ein_hits )
{
  // treat it like the bottom of the box
  find_bottom_hits( E_out_lab, Ein_hits );
  for( std::vector< Box::Ein_Eta_Hit >::iterator this_hit = Ein_hits->begin( );
       this_hit != Ein_hits->end( ); ++this_hit )
  {
    if( this_hit->hit_edge == Box::BOTTOM_OUT )
    {
      this_hit->hit_edge = Box::TOP_IN;
    }
    else if( this_hit->hit_edge == Box::BOTTOM_IN )
    {
      this_hit->hit_edge = Box::TOP_OUT;
    }
  }
}

// ************* class Twohit::two_step_hit_list_array *****************
// ----------- Twohit::two_step_hit_list_array::get_array --------------
// allocates space
void Twohit::two_step_hit_list_array::get_array( int Laengd )
{
  laengd =  Laengd;
  hit_lists = new Twohit::two_step_hit_list[ Laengd ];
}
// ----------- Twohit::two_step_hit_list_array::common_hits --------------
//! Ensures that these hit_lists are defined at the same incident energies.  
void Twohit::two_step_hit_list_array::common_hits( )
{
  // First, get all of the incident energies
  std::list< double > Ein_list;
  for( int list_count = 0; list_count < laengd; ++list_count )
  {
    for( Box::hit_list_base::iterator hit_ptr = hit_lists[ list_count ].begin( );
         hit_ptr != hit_lists[ list_count ].end( ); ++hit_ptr )
    {
      Ein_list.push_back( hit_ptr->E_in );
    }
  }

  Ein_list.sort( );
  Ein_list.unique( );

  // Convert the energy list to a vector
  std::vector< double > Ein_vector;
  for( std::list< double >::const_iterator list_ptr = Ein_list.begin( );
       list_ptr != Ein_list.end( ); ++list_ptr )
  {
    Ein_vector.push_back( *list_ptr );
  }

  // Fill in each hit_list
  for( int list_count = 0; list_count < laengd; ++list_count )
  {
    hit_lists[ list_count ].fill_in( Ein_vector );
  }
}
