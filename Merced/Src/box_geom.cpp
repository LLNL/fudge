/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: angle_dist.cpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */

#include "messaging.hpp"
#include "box_geom.hpp"
#include "global_params.hpp"

// ********* class Box::Ein_Eta_Hit *********
// ---------------- Box::Ein_Eta_Hit copy constructor ------------------
Box::Ein_Eta_Hit::Ein_Eta_Hit( const Box::Ein_Eta_Hit& ein_eta_hit )
{
  E_in = ein_eta_hit.E_in;
  hit_edge = ein_eta_hit.hit_edge;
}
// ---------------- Box::Ein_Eta_Hit::print ------------------
// Prints a hit
void Box::Ein_Eta_Hit::print( ) const
{
  std::cout << " E_in: " << E_in << " " << edge_name( );
}

// ---------------- Box::Ein_Eta_Hit::edge_name ------------------
// Prints a hit
std::string Box::Ein_Eta_Hit::edge_name( ) const {

  std::string name;

  switch( hit_edge )
  {
  case Box::INSIDE:
    name = "INSIDE";
    break;
  case Box::BOTTOM_IN:
    name = "BOTTOM_IN";
    break;
  case Box::BOTTOM_OUT:
    name = "BOTTOM_OUT";
    break;
  case Box::TOP_IN:
    name = "TOP_IN";
    break;
  case Box::TOP_OUT:
    name = "TOP_OUT";
    break;
  case Box::BELOW:
    name = "BELOW";
    break;
  case Box::ABOVE:
    name = "ABOVE";
    break;
  case Box::ABOVE_FORWARD:
    name = "ABOVE_FORWARD";
    break;
  case Box::BOTTOM_CORNER_IN:
    name = "BOTTOM_CORNER_IN";
    break;
  case Box::BOTTOM_CORNER_OUT:
    name = "BOTTOM_CORNER_OUT";
    break;
  case Box::TOP_CORNER_IN:
    name = "TOP_CORNER_IN";
    break;
  case Box::TOP_CORNER_OUT:
    name = "TOP_CORNER_OUT";
    break;
  case Box::MISSED:
    name = "MISSED";
    break;
  }

  return name;
}

// ********* class Box::hit_list_base *********
// ---------------- Box::hit_list_base::hit_box ------------------
// Finds intersections of the curve eta = const with the E-E' box
bool Box::hit_list_base::hit_box( double Eta, std::vector< double >::const_iterator Eout_ptr,
  double E_in_left, double E_in_right )
{
  double E_measure = ( (  E_in_right - E_in_left ) < E_in_right ) ?
    (  E_in_right - E_in_left ) : E_in_right;
  static double etol = Global.Value( "looser_tol" );
  box_tol = etol * E_measure;

  // for the top of the E-E' box
  std::vector< double >::const_iterator next_Eout = Eout_ptr;
  ++next_Eout;

  // start with a new list
  if( size( ) > 0 )
  {
    erase( begin( ), end( ) );
  }

  eta = Eta;

  // test the bottom of the E-E' box
  test_bottom( Eout_ptr, E_in_left, E_in_right );

  // test the top of the E-E' box
  test_top( next_Eout, E_in_left, E_in_right );

  if( empty( ) )
  {
    // we missed but above, below, or through the middle?
    test_sides( Eout_ptr, E_in_left, E_in_right );
  }
  else
  {
    // test the left side of the eta-E_in region
    test_left( E_in_left, Eout_ptr );

    // test the right side of the E-E' box
    test_right( E_in_right, Eout_ptr );
  }

  // test for consistency
  bool is_OK = consistency( );
  //  if( !is_OK ) print( );
  return is_OK;
}
// ---------------- Box::hit_list_base::test_bottom ------------------
// Finds intersections of the curve eta = const with the bottom of the E-E' box
void Box::hit_list_base::test_bottom( std::vector< double >::const_iterator Eout_ptr,
    double E_in_left, double E_in_right )
{
  // the intersections
  std::vector< Box::Ein_Eta_Hit > Ein_hits;

  // find the intersections
  find_bottom_hits( *Eout_ptr, &Ein_hits );

  // are the hits inside the range of E_in?
  std::vector< Box::Ein_Eta_Hit >::iterator hit_ptr = Ein_hits.begin( );
  for( ; hit_ptr != Ein_hits.end( ); ++hit_ptr )
  {
    // interior intersections
    if( ( hit_ptr->E_in > E_in_left + box_tol ) &&
	( hit_ptr->E_in < E_in_right - box_tol ) )
    {
      stick_in( *hit_ptr );
    }
  }
}
// ---------------- Box::hit_list_base::test_top ------------------
// Finds intersections of the curve eta = const with the top of the E-E' box
void Box::hit_list_base::test_top( std::vector< double >::const_iterator Eout_ptr,
    double E_in_left, double E_in_right )
{
  // the intersections
  std::vector< Box::Ein_Eta_Hit > Ein_hits;

  // find the intersections
  find_top_hits( *Eout_ptr, &Ein_hits );

  // are the hits inside the range of E_in?
  std::vector< Box::Ein_Eta_Hit >::iterator hit_ptr = Ein_hits.begin( );
  for( ; hit_ptr != Ein_hits.end( ); ++hit_ptr )
  {
    // interior intersections
    if( ( hit_ptr->E_in > E_in_left + box_tol ) &&
	( hit_ptr->E_in < E_in_right - box_tol ) )
    {
      stick_in( *hit_ptr );
    }
  }
}
// ---------------- Box::hit_list_base::common_hits ------------------
// Ensures that these hits and upper_hits are defined at the same incident energies.
void Box::hit_list_base::common_hits( Box::hit_list_base &upper_hits )
{
  // set up Ein_hits, the vector of incident energies for which
  // the quadrature limits eta = const meet the edges of the energy
  // quadrature box.
  std::vector< double > Ein_lower;
  for( Box::hit_list_base::iterator hit_ptr = begin( );
       hit_ptr != end( ); ++hit_ptr )
  {
    Ein_lower.push_back( hit_ptr->E_in );
  }
  upper_hits.fill_in( Ein_lower );

  std::vector< double > Ein_upper;
  for( Box::hit_list_base::iterator hit_ptr = upper_hits.begin( );
       hit_ptr != upper_hits.end( ); ++hit_ptr )
  {
    Ein_upper.push_back( hit_ptr->E_in );
  }
  fill_in( Ein_upper );
}
// ---------------- Box::hit_list_base::stick_in ------------------
void Box::hit_list_base::stick_in( const Box::Ein_Eta_Hit& ein_eta_hit )
{
  if( empty( ) )
  {
    push_back( ein_eta_hit );
  }
  else
  {
    // find the next bigger incident energy
    Box::hit_list_base::iterator next_hit = begin( );
    if( next_hit->E_in > ein_eta_hit.E_in )
    {
      push_front( ein_eta_hit );
    }
    else
    {
      Box::hit_list_base::iterator prev_hit = next_hit;
      for( ++next_hit; next_hit != end( ); prev_hit = next_hit, ++next_hit )
      {
        if( next_hit->E_in > ein_eta_hit.E_in )
        {
          break;
        }
      }
      if( prev_hit->E_in < ein_eta_hit.E_in )
      {
        insert( next_hit, ein_eta_hit );
      }
    }
  }
}
// ----------- Box::hit_list_base::is_below --------------
// Checks whether this eta = const hyperbola is below the E-E' box
bool Box::hit_list_base::is_below( )
{
  for( Box::hit_list_base::iterator this_link = begin( );
       this_link != end( ); ++this_link )
  {
    if( this_link->hit_edge != Box::BELOW )
    {
      return false;
    }
  }
  return true;
}
// ----------- Box::hit_list_base::is_above --------------
// Checks whether this eta = const hyperbola is above the E-E' box
bool Box::hit_list_base::is_above( )
{
  Box::hit_list_base::iterator this_link = begin( );
  bool is_above = false;
  if( this_link->hit_edge == Box::ABOVE )
  {
    is_above = true;
    for( ++this_link; this_link != end( ); ++this_link )
    {
      if( this_link->hit_edge != Box::ABOVE )
      {
        is_above = false;
	break;
      }
    }
  }
  else if( this_link->hit_edge == Box::ABOVE_FORWARD )
  {
    is_above = true;
    for( ++this_link; this_link != end( ); ++this_link )
    {
      if( this_link->hit_edge != Box::ABOVE_FORWARD )
      {
        is_above = false;
	break;
      }
    }
  }

  return ( is_above );
}
// ----------- Box::hit_list_base::fill_in --------------
// fills in the list with energies from Ein_hits
void Box::hit_list_base::fill_in( const std::vector< double >& Ein_hits )
{
  std::vector< double >::const_iterator this_Ein = Ein_hits.begin( );
  Box::hit_list_base::iterator this_hit = begin( );
  Box::hit_list_base::iterator next_hit = this_hit;
  ++next_hit;
  // for adding new links
  Box::Ein_Eta_Hit new_link;

  // we may assume that our list and Ein_hits start together
  // loop through the rest of Ein_hits
  for( ++this_Ein; ( this_Ein != Ein_hits.end( ) ) &&
         ( next_hit != end( ) );  )
  {
    // do we have an intermediate energy?
    if( *this_Ein < next_hit->E_in )
    {
      // set up a new link
      new_link.E_in = *this_Ein;
      // the eta = const hyperbola at this E_in value must be
      // either below, inside, or above the E-E' quadrature box
      if( ( this_hit->hit_edge == Box::ABOVE ) ||
          ( next_hit->hit_edge == Box::ABOVE ) )
      {
        new_link.hit_edge = Box::ABOVE;
      }
      else if( ( this_hit->hit_edge == Box::ABOVE_FORWARD ) ||
          ( next_hit->hit_edge == Box::ABOVE_FORWARD ) )
      {
        new_link.hit_edge = Box::ABOVE_FORWARD;
      }
      else if( ( this_hit->hit_edge == Box::BELOW ) ||
               ( this_hit->hit_edge == Box::BOTTOM_OUT ) ||
               ( next_hit->hit_edge == Box::BELOW ) )
      {
        new_link.hit_edge = Box::BELOW;
      }
      else
      {
        new_link.hit_edge = Box::INSIDE;
      }
      insert( next_hit, new_link );
      ++this_Ein;
    }
    else if( *this_Ein > next_hit->E_in )
    {
      // we already have this incident energy; go on to the next
      this_hit = next_hit;
      ++next_hit;
    }
    else
    {
      // go on to the next E_in and the next link
      ++this_Ein;
      this_hit = next_hit;
      ++next_hit;
    }
  }
}
// ---------------- Box::hit_list_base::test_left_default ------------------
// Where do we hit the left-hand side of the box?
void Box::hit_list_base::test_left_default( double E_in, std::vector< double >::const_iterator Eout_ptr )
{
  // for a new entry
  Box::Ein_Eta_Hit Ein_eta_hit;
  Ein_eta_hit.E_in = E_in;

  // we have found intersections with the top or bottom
  hit_list::iterator first_hit = begin( );
  if( first_hit->hit_edge == Box::TOP_IN )
  {
    Ein_eta_hit.hit_edge = Box::ABOVE;
  }
  else if( first_hit->hit_edge == Box::BOTTOM_IN )
  {
    Ein_eta_hit.hit_edge = Box::BELOW;
  }
  else if( ( first_hit->hit_edge == Box::BOTTOM_OUT ) ||
           ( first_hit->hit_edge == Box::TOP_OUT ) )
  {
    Ein_eta_hit.hit_edge = Box::INSIDE;
  }
  push_front( Ein_eta_hit );
}
// ---------------- Box::hit_list_base::test_right_default ------------------
// Where do we hit the right-hand side of the box?
void Box::hit_list_base::test_right_default( double E_in, std::vector< double >::const_iterator Eout_ptr )
{
  // for a new entry
  Box::Ein_Eta_Hit Ein_eta_hit;
  Ein_eta_hit.E_in = E_in;

  // we already have at least one entry (at the left edge)
  hit_list::iterator last_hit = end( );
  --last_hit;

  if( ( last_hit->hit_edge == Box::TOP_OUT ) ||
      ( last_hit->hit_edge == Box::ABOVE ) )
  {
    Ein_eta_hit.hit_edge = Box::ABOVE;
  }
  else if( ( last_hit->hit_edge == Box::BOTTOM_OUT ) ||
      ( last_hit->hit_edge == Box::BELOW ) )
  {
    Ein_eta_hit.hit_edge = Box::BELOW;
  }
  else
  {
    Ein_eta_hit.hit_edge = Box::INSIDE;
  }

  // append this entry
  push_back( Ein_eta_hit );
}
// ---------------- Box::hit_list_base::print ------------------
// Prints the list
void Box::hit_list_base::print( )
{
  for( Box::hit_list_base::const_iterator ptr = begin( ); ptr != end( ); ++ptr )
  {
    ptr->print( );
  }
  std::cout << std::endl;
}
// ----------- Box::hit_list_base::consistency_default --------------
// Checks for inconsistencies caused by peculiarities of real arithmetic
bool Box::hit_list_base::consistency_default( )
{
  Box::hit_list_base::iterator this_link = begin( );
  Box::hit_list_base::iterator next_link = this_link;
  ++next_link;

  // loop through the intersections
  for( ;
       next_link != end( );
       this_link = next_link, ++next_link )
  {
    if( ( this_link->hit_edge == Box::BELOW ) || ( this_link->hit_edge == Box::BOTTOM_OUT ) )
    {
      if( ( next_link->hit_edge != Box::BELOW ) && ( next_link->hit_edge != Box::BOTTOM_IN ) )
      {
	Msg::Warning( "Box::hit_list_base::consistency_default",
		      "Check coding, 1" );
        return false;
      }
    }
    else if( this_link->hit_edge == Box::INSIDE )
    {
      if( ( next_link->hit_edge !=  Box::TOP_OUT ) &&
          ( next_link->hit_edge != Box::INSIDE ) &&
	  ( next_link->hit_edge != Box::BOTTOM_OUT ) )
      {
        Msg::Warning( "Box::hit_list_base::consistency_default",
		      "Check coding, 2" );
        return false;
      }
    }
    else if( ( this_link->hit_edge == Box::ABOVE ) || ( this_link->hit_edge == Box::TOP_OUT ) )
    {
      if( ( next_link->hit_edge !=  Box::ABOVE ) && ( next_link->hit_edge != Box::TOP_IN ) )
      {
	Msg::Warning( "Box::hit_list_base::consistency_default",
		      "Check coding, 3" );
        return false;
      }
    }
    else if( this_link->hit_edge == Box::BOTTOM_IN )
    {
      if( ( next_link->hit_edge != Box::TOP_OUT ) &&
	  ( next_link->hit_edge != Box::INSIDE ) &&
	  ( next_link->hit_edge != Box::BOTTOM_OUT ) )
      {
        Msg::Warning( "Box::hit_list_base::consistency_default",
		      "Check coding, 4" );
        return false;
      }
    }
    else if( this_link->hit_edge == Box::TOP_IN )
    {
      if( ( next_link->hit_edge != Box::BOTTOM_OUT ) &&
	  ( next_link->hit_edge != Box::INSIDE ) &&
	  ( next_link->hit_edge != Box::TOP_OUT ) )
      {
        Msg::Warning( "Box::hit_list_base::consistency_default",
		      "Check coding, 5" );
        return false;
      }
    }
  }
  return true;   // no inconsistencies
}

// ********* class Box::hit_list *********
// ---------------- Box::hit_list::test_sides ------------------
// No top or bottom intersections, so are we above, below, or through the middle?
void Box::hit_list::test_sides( std::vector< double >::const_iterator Eout_ptr,
  double E_in_left, double E_in_right )
{
  // for new entries
  Box::Ein_Eta_Hit Ein_eta_hit;
  std::vector< double >::const_iterator next_Eout = Eout_ptr;
  ++next_Eout;
  // test in the middle, because the curve may have hit a corner
  double Eout = get_Eout( 0.5*( E_in_left + E_in_right ) );
  if( Eout > *next_Eout )
  {
    Ein_eta_hit.E_in = E_in_left;
    Ein_eta_hit.hit_edge = Box::ABOVE;
    push_back( Ein_eta_hit );
    Ein_eta_hit.E_in = E_in_right;
    push_back( Ein_eta_hit );
  }
  else if( Eout < *Eout_ptr )
  {
    Ein_eta_hit.E_in = E_in_left;
    Ein_eta_hit.hit_edge = Box::BELOW;
    push_back( Ein_eta_hit );
    Ein_eta_hit.E_in = E_in_right;
    push_back( Ein_eta_hit );
  }
  else
  {
    Ein_eta_hit.E_in = E_in_left;
    Ein_eta_hit.hit_edge = Box::INSIDE;
    push_back( Ein_eta_hit );
    Ein_eta_hit.E_in = E_in_right;
    push_back( Ein_eta_hit );
  }
}

// ************* class Box::energy_hit_list *****************
// ----------- Box::energy_hit_list::get_Eout --------------
// Calculates E_out from E_in, for testing the side of a box
double Box::energy_hit_list::get_Eout( double e_in )
{
  // Ignore this bool
  bool is_OK;
  double ans = E_Eout.value( e_in, &is_OK );
  return ans;
}
// ----------- Box::energy_hit_list::find_bottom_hits --------------
// Finds the intersections with the bottom of a box
void Box::energy_hit_list::find_bottom_hits( double E_out,
  std::vector< Box::Ein_Eta_Hit > *Ein_hits )
{
  // for new entries
  Box::Ein_Eta_Hit Ein_eta_hit;
  int flag;
  double Ein = E_Eout.find_Ein( E_out, &flag );

  // how many intersections?
  if( flag == 1 )
  {
    // E_lab is increasing with E_in
    // append this entry
    Ein_eta_hit.E_in = Ein;
    Ein_eta_hit.hit_edge = ( E_Eout.second.y > E_Eout.first.y ) ?
      Box::BOTTOM_IN : Box::BOTTOM_OUT;
    Ein_hits->push_back( Ein_eta_hit );
  }
  else if( flag == 2 )
  {
    // This case should be covered by the test_side routine
  }
}
// ----------- Box::energy_hit_list::find_top_hits --------------
// Finds the intersections with the top of a box
void Box::energy_hit_list::find_top_hits( double E_out,
  std::vector< Box::Ein_Eta_Hit > *Ein_hits )
{
  // for new entries
  Box::Ein_Eta_Hit Ein_eta_hit;
  int flag;
  double Ein = E_Eout.find_Ein( E_out, &flag );

  // how many intersections?
  if( flag == 1 )
  {
    // E_lab is increasing with E_in
    // append this entry
    Ein_eta_hit.E_in = Ein;
    Ein_eta_hit.hit_edge = ( E_Eout.second.y > E_Eout.first.y ) ?
      Box::TOP_OUT : Box::TOP_IN;
    Ein_hits->push_back( Ein_eta_hit );
  }
  else if( flag == 2 )
  {
    // This case should be covered by the test_side routine
  }
}
