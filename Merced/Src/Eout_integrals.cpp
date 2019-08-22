/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2009-08-11 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: Eout_integrals.cpp 1  2009-08-11 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */
// implementation for the array of integrals over Eout/mu

#include <cmath>

#include "Eout_integrals.hpp"
#include "messaging.hpp"
#include "global_params.hpp"

// *************** class Eout_link ***************************
// --------------- Eout_link::Interpolate ---------------------------
// Linear interpolation
void Eout_link::Interpolate( double Ein, const Eout_link& next_link,
		    coef_vector *interp ) const
{
  if( next_link.E_in == E_in )
  {
    FatalError( "Eout_link::Interpolate", "division by zero" );
  }
  if( ( interp->conserve != conserve ) || ( interp->order != order ) )
  {
    FatalError( "Eout_link::Interpolate", "inconsistent vector type" );
  }

  double alpha = ( Ein - E_in )/( next_link.E_in - E_in );
  int i;
  if( ( conserve == NUMBER ) || ( conserve == BOTH ) )
  {
    for( i = 0; i <= order; ++i )
    {
      interp->weight_1[ i ] = ( 1.0 - alpha ) * weight_1[ i ] +
	alpha * next_link.weight_1[ i ];
    }
  }
  if( ( conserve == ENERGY ) || ( conserve == BOTH ) )
  {
    for( i = 0; i <= order; ++i )
    {
      interp->weight_E[ i ] = ( 1.0 - alpha ) * weight_E[ i ] +
	alpha * next_link.weight_E[ i ];
    }
  }
}

// *************** class Eout_integrals ***************************
// --------------- Eout_integrals::setup_Eout_ints ---------------------------
// Sets up the integrals over E_out for interpolation
void Eout_integrals::setup_Eout_ints( int order, Conserve conserve, double Ein_0,
  double Ein_1 )
{
  if( !empty( ) )
  {
    erase( begin( ), end( ) );
  }
  static double E_tol = Global.Value( "E_tol" );
  static int max_Eout_ints_size = Global.Value( "max_Eout_ints_size" );
  Eout_integrals::iterator where = end( );  // where to insert the first entry
  new_Eout_int( where, order, conserve, Ein_0 );   // first incident energy
  new_Eout_int( where, order, conserve, Ein_1 );   // last incident energy
  // insert a link and test the noise in the function evaluation
  --where;
  double tol = tol_new_Eout_int( where, order, conserve, 0.5*( Ein_0 + Ein_1 ) );
  // now thicken as in fete
  Eout_integrals::iterator prev_link = begin( );
  where = prev_link;
  ++where;
  Eout_integrals::iterator next_link = where;
  ++next_link;
  // loop through the list
  for(; next_link != end( ); )
  {
    // keep checking until the interval passes the test
    for(;;)
    {
      if( interp_OK( *where, prev_link, next_link, tol ) )
      {
        prev_link = where;
        where = next_link;
        ++next_link;
        break;
      }
      else
      {
        // don't subdivide really small intervals
        if( where->E_in < prev_link->E_in + E_tol )
        {
          prev_link = where;
          where = next_link;
          ++next_link;
          break;
        }
        // safety check
        if( static_cast<int>( size( ) ) == max_Eout_ints_size )
        {
          Warning( "Eout_integrals::setup_Eout_ints", 
                   pastenum("got ", max_Eout_ints_size ) +
                   " in the Eout_ints list" );
          prev_link = where;
          where = next_link;
          ++next_link;
          break;
        }
        // insert new links at the midpoint energies
        new_Eout_int( where, order, conserve, 0.5*( prev_link->E_in + where->E_in ) );
        new_Eout_int( next_link, order, conserve, 0.5*( next_link->E_in + where->E_in ) );
        next_link = where;
        --where;
      }
    }
  }
}
// --------------- Eout_integrals::new_Eout_int ---------------------------
// Adds a new link to the Eout_ints list
void Eout_integrals::new_Eout_int( Eout_integrals::iterator where, int order, 
  Conserve conserve, double E_in )
{
  Eout_integrals::iterator new_entry = insert( where, Eout_link( ) );  // insert before where
  new_entry->set_order( order, conserve );
  Eout_int_params->set_Ein( E_in );
  new_entry->E_in = E_in;
  Eout_int_params->get_integrals( Eout_int_params->Eout_0, Eout_int_params->Eout_1,
    *new_entry );
}
// --------------- Eout_integrals::tol_new_Eout_int ---------------------------
// Adds a new link to the Eout_ints list and returns the noise in the calculation
double Eout_integrals::tol_new_Eout_int( Eout_integrals::iterator where, int order, 
  Conserve conserve, double E_in )
{
  Eout_integrals::iterator new_entry = insert( where, Eout_link( ) );  // insert before where
  new_entry->set_order( order, conserve );
  Eout_int_params->set_Ein( E_in );
  new_entry->E_in = E_in;
  double noise = Eout_int_params->tol_get_integrals( Eout_int_params->Eout_0,
    Eout_int_params->Eout_1, *new_entry );
  return noise;
}
// --------------- Eout_integrals::interp_OK ---------------------------
// Check the accuracy of linear interpolation
bool Eout_integrals::interp_OK( const Eout_link& interp_link,
		  Eout_integrals::iterator prev_link,
		  Eout_integrals::iterator next_link, double tol )
{
  bool isOK = true;
  double alpha = ( interp_link.E_in - prev_link->E_in )/
    ( next_link->E_in - prev_link->E_in );
  int i;
  if( ( interp_link.conserve == NUMBER ) || ( interp_link.conserve == BOTH ) )
  {
    for( i = 0; i <= interp_link.order; ++i )
    {
      if( abs( ( 1.0 - alpha ) * prev_link->weight_1[ i ] +
	    alpha * next_link->weight_1[ i ] - interp_link.weight_1[ i ] ) >
	  tol*abs( interp_link.weight_1[ i ] ) )
      {
	isOK = false;
	break;
      }
    }
  }
  if( ( interp_link.conserve == ENERGY ) || ( interp_link.conserve == BOTH ) )
  {
    for( i = 0; i <= interp_link.order; ++i )
    {
      if( abs( ( 1.0 - alpha ) * prev_link->weight_E[ i ] +
	    alpha * next_link->weight_E[ i ] - interp_link.weight_E[ i ] ) >
	  tol*abs( interp_link.weight_E[ i ] ) )
      {
	isOK = false;
	break;
      }
    }
  }
  return isOK;
}
