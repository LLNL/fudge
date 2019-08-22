/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2016-02-12 11:40:00 -0800 (Fri, 12 Feb 2016) $
 * $Author: hedstrom $
 * $Id: isotropic.hpp 1 2016-02-12 11:40:00Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// implement the classes for isotropic distributions in the lab frame.

#ifdef _OPENMP
 #include <omp.h>
#endif

#include "isotropic.hpp"
#include "messaging.hpp"
#include "global_params.hpp"

// *********************** class isotropic *******************
// ------------------ isotropic::read_data -----------------
void isotropic::read_data( data_parser& input_file, int num_Ein )
{
  data_order = 0;
  // this is a kludge to use old code
  ENDL_data = new energy_dist[ data_order + 1 ];
  energy_dist *new_moment_ptr = &ENDL_data[ 0 ];
  new_moment_ptr->read_data( input_file, num_Ein );

  // we need to set the outgoing interpolation
  Eprob_vector *Ein_ptr;
  for( int Ein_count = 0; Ein_count < num_Ein; ++Ein_count )
  {
    Ein_ptr = &new_moment_ptr->EProb_data[ Ein_count ];
    Ein_ptr->interp_type = Eout_interp;
    // ensure proper normalization
    Ein_ptr->renorm( );
  }

  // the energy_dist read_data doesn't do unit-base interpolation
  if( Ein_interp.qualifier == UNITBASE )
  {
    new_moment_ptr->unit_base( 0 );
  }
}
// ------------------ isotropic::get_T -----------------
// Calculates the transfer matrix for this particle.
// sigma is the cross section.
void isotropic::get_T( const dd_vector& sigma,
  const dd_vector& mult, const dd_vector& weight, T_matrix& transfer )
{ 
  bool interp_OK = ( ( Ein_interp.qualifier == UNITBASE ) &&
		     ( ( Ein_interp.flag == LINLIN ) ||
     		       ( Ein_interp.flag == LINLOG ) ) ) ||
                   ( ( Ein_interp.qualifier == CUMULATIVE_POINTS ) &&
		     ( ( Ein_interp.flag == LINLIN ) ||
		       ( Ein_interp.flag == LINLOG ) ) ) ||
    ( ( Ein_interp.qualifier == DIRECT ) &&
      ( ( Ein_interp.flag == LINLIN ) ||
        ( Ein_interp.flag == HISTOGRAM ) ) );

  if( !interp_OK )
  {
    FatalError( "isotropic::get_T",
      "Incident energy interpolation not implemented" );
  }
  interp_OK = ( Eout_interp == LINLIN ) || ( Eout_interp == HISTOGRAM );
  if( !interp_OK )
  { 
    FatalError( "isotropic::get_T",
      "Outgoing energy interpolation not implemented" );
  }

  // convert to ENDF format
  zero_order( );

  ENDF_data.Ein_interp = Ein_interp;
  ENDF_data.Eout_interp = Eout_interp;
  ENDF_data.get_T( sigma, mult, weight, transfer );
}
