/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2014-08-05 19:06:56 -0800 (Tue, 05 Aug 2014) $
 * $Author: hedstrom $
 * $Id: Energy_groups.hpp 1 2014-08-05 03:06:56Z hedstrom $
 *
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */
// implementation of the class for energy groups

#include "Energy_groups.hpp"
#include "global_params.hpp"
#include "messaging.hpp"

//****************** class Energy_groups  *****************
//--------------- Energy_groups::read_bd -----------------
// Reads the energy group boundaries
void Energy_groups::read_bd( data_parser &input_file, int num_bd )
{
  for( int bd_count = 0; bd_count < num_bd; ++bd_count )
  {
    double bound = input_file.get_next_double( );
    push_back( bound );
  }
}
//--------------- Energy_groups::first_bin_ID -----------------
// Returns the index of the left-hand end of the energy bin containing this_E
int Energy_groups::first_bin_ID( double this_E ) const
{
  int bin_ID = 0;
  Energy_groups::const_iterator this_bin = begin( );
  for( ++this_bin; this_bin != end( ); ++this_bin, ++bin_ID )
  {
    if( *this_bin > this_E )
    {
      break;
    }
  }
  if( this_bin == end( ) )
  {
    FatalError( "Energy_groups::first_bin_ID", "energy too high" );
  }
  return bin_ID;
}
//--------------- Energy_groups::last_bin_ID -----------------
// Returns the index of the right-hand end of the energy bin containing this_E
int Energy_groups::last_bin_ID( double this_E ) const
{
  int bin_ID = size( ) - 1;
  Energy_groups::const_iterator this_bin = end( );
  --this_bin;
  Energy_groups::const_iterator lower_bin = this_bin;
  --lower_bin;

  static double skip_tol = Global.Value( "abs_tol" );
  for( ; this_bin != begin( ); this_bin = lower_bin, --lower_bin, --bin_ID )
  {
    if( *lower_bin < this_E * ( 1.0 - skip_tol ) )
    {
      break;
    }
  }
  if( this_bin == begin( ) )
  {
    FatalError( "Energy_groups::first_bin_ID", "energy too high" );
  }
  return bin_ID;
}
