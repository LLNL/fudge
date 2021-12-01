/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2011-04-28 19:06:56 -0800 (Fri, 28 Feb 2011) $
 * $Author: hedstrom $
 * $Id: energy_dist.cpp 1 2011-04-28 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// implementation of the classes used to handle energy distributions

#include <cmath>

#include "energy_dist.hpp"
#include "global_params.hpp"
#include "messaging.hpp"

// ----------- Edist::energy_dist::destructor --------------
Edist::energy_dist::~energy_dist( )
{
  if( number_Ein > 0 )
  {
    delete [] EProb_data;
  }
}
// ----------- Edist::energy_dist::read_data --------------
// Reads the ENDL data for one Legendre order
void Edist::energy_dist::read_data( Dpar::data_parser& input_file, int num_Ein )
{
  number_Ein = num_Ein;
  EProb_data = new Ebase::Eprob_vector[ num_Ein ];
  Ebase::Eprob_vector *new_energy_ptr;
  // loop over incident energy
  for( int Ein_count = 0; Ein_count < num_Ein; ++Ein_count )
  {
    // make a new energy distribution
    new_energy_ptr = &( EProb_data[ Ein_count ] );
    new_energy_ptr->set_E_in( input_file.get_next_double( ) );
    new_energy_ptr->interp_type = Eout_interp;
    // read the (energy, probability density) pairs
    int num_Eout = input_file.get_next_int( );
    for( int Eout_count = 0; Eout_count < num_Eout; ++Eout_count )
    {
      double E_out = input_file.get_next_double( );
      double Prob = input_file.get_next_double( );
      new_energy_ptr->add_entry( E_out, Prob );
    }
  }
}
// ----------- Edist::energy_dist::unit_base --------------
// Maps the data to unit base
void Edist::energy_dist::unit_base( int L_order )
{
  for( int Ein_count = 0; Ein_count < number_Ein; ++Ein_count )
  {
    EProb_data[ Ein_count ].unit_base( L_order );
  }
}

// ************* class Edist::energy_moments *****************
// ----------- Edist::energy_moments::destructor --------------
Edist::energy_moments::~energy_moments( )
{
  if( data_order >= 0 )
  {
    delete [] ENDL_data;
  }
}
// ----------- Edist::energy_moments::read_data --------------
void Edist::energy_moments::read_data( Dpar::data_parser& input_file, int num_moments )
{
  // Read the interpolation rules
  interp_flag_F::read_2d_interpolation( input_file, &Ein_interp, &Eout_interp );
  data_order = num_moments - 1;
  ENDL_data = new Edist::energy_dist[ data_order + 1 ];
  Edist::energy_dist *new_moment_ptr;

  // read the data
  for( int Legendre_order = 0; Legendre_order <= data_order;
    ++Legendre_order )
  {
    // the next Legendre order
    new_moment_ptr = &ENDL_data[ Legendre_order ];
    // is the Legendre order consistent?
    int file_order = input_file.get_next_int( );
    if( file_order != Legendre_order )
    {
      Msg::FatalError( "Edist::energy_moments::read_data",
		  Msg::pastenum( "improper Legendre order: ", file_order ) );
    }
    new_moment_ptr->Ein_interp = Ein_interp;
    new_moment_ptr->Eout_interp = Eout_interp;
    // how many incident energies
    int num_Ein = input_file.get_next_int( );
    new_moment_ptr->read_data( input_file, num_Ein );
    if( Ein_interp.qualifier == Terp::UNITBASE )
    {
      new_moment_ptr->unit_base( Legendre_order );
    }
  }
}
// ----------- Edist::energy_moments::zero_order --------------
// Converts isotropic ENDL data to ENDF format
void Edist::energy_moments::zero_order( )
{
  // copy the interpolation rules
  ENDF_data.Ein_interp.flag = Ein_interp.flag;
  ENDF_data.Ein_interp.qualifier = Ein_interp.qualifier;
  ENDF_data.Eout_interp = Eout_interp;

  Edist::energy_dist *ENDL_ptr = &ENDL_data[ 0 ];
  int num_Ein = ENDL_ptr->number_Ein;
  Ebase::Eprob_vector *ENDL_Ein_ptr;
  Ebase::Eprob_vector::const_iterator ENDL_Eout_ptr;
  StdLg::standard_Legendre::iterator ENDF_Ein_ptr;
  StdLg::standard_Legendre_vector::iterator ENDF_Eout_ptr;
  // iterate through the incident energies
  for( int Ein_count = 0; Ein_count < num_Ein; ++Ein_count )
  {
    ENDL_Ein_ptr = &( ENDL_ptr->EProb_data[ Ein_count ] );
    // make a new StdLg::standard_Legendre_vector
    ENDF_Ein_ptr = ENDF_data.insert( ENDF_data.end( ), StdLg::standard_Legendre_vector( ) );
    ENDF_Ein_ptr->set_E_in( ENDL_Ein_ptr->get_E_in( ) );  // energy of incident particle
    ENDF_Ein_ptr->Eout_interp = Eout_interp;
    if( Ein_interp.qualifier == Terp::UNITBASE )
    {
      ENDF_Ein_ptr->ubase_map.copy( ENDL_Ein_ptr->ubase_map );
    }

    // go through the pairs ( E_out, probability_density )
    for( ENDL_Eout_ptr = ENDL_Ein_ptr->begin( );
	 ENDL_Eout_ptr != ENDL_Ein_ptr->end( ); ++ENDL_Eout_ptr )
    {
      // make a new set of Legendre coefficients
      ENDF_Eout_ptr = ENDF_Ein_ptr->insert( ENDF_Ein_ptr->end( ), Lgdata::Legendre_coefs( ) );
      ENDF_Eout_ptr->initialize( 0 );
      ENDF_Eout_ptr->set_E_out( ENDL_Eout_ptr->x );  // energy of outgoing particle
      ( *ENDF_Eout_ptr )[ 0 ] = ENDL_Eout_ptr->y;
    }
    if( Ein_interp.qualifier == Terp::CUMULATIVE_POINTS )
    {
      ENDF_Ein_ptr->form_cum_prob( );
    }
  }
}
// ----------- Edist::energy_moments::to_ENDF --------------
// Converts ENDL data to ENDF format
void Edist::energy_moments::to_ENDF( )
{
  // Do we need to interpolate with respect to incident energy?
  check_Ein( );
  // copy the interpolation rules
  ENDF_data.Ein_interp.flag = Ein_interp.flag;
  ENDF_data.Ein_interp.qualifier = Ein_interp.qualifier;
  ENDF_data.Eout_interp = Eout_interp;

  // iterate over incident energy
  int num_Ein = ENDL_data[ 0 ].number_Ein;
  for( int energy_count = 0; energy_count < num_Ein; ++energy_count )
  {
    one_Ein_to_ENDF( energy_count );
  }
}
// ----------- Edist::energy_moments::check_Ein --------------
// Checks to see that the incident energies are consistent for all Legendre orders
void Edist::energy_moments::check_Ein( ) const
{
  // first, check the amount of data
  int num_Ein = ENDL_data[ 0 ].number_Ein;
  for( int L_order = 1; L_order <= data_order; ++L_order )
  {
    if( ENDL_data[ L_order ].number_Ein != num_Ein )
    {
      Msg::FatalError( "Edist::energy_moments::check_Ein",
		"Implement differing numbers of incident energies" );
    }
  }
  // now check the incident energies
  static double tol = Global.Value( "tight_tol" );
  for( int Ein_count = 0; Ein_count < num_Ein; ++Ein_count )
  {
    double this_Ein = ENDL_data[ 0 ].EProb_data[ Ein_count ].get_E_in( );
    for( int L_order = 1; L_order <= data_order; ++L_order )
    {
      if( std::abs( ENDL_data[ L_order ].EProb_data[ Ein_count ].get_E_in( ) - this_Ein )>
	  tol*this_Ein )
      {
        Msg::FatalError( "Edist::energy_moments::check_Ein",
		"Implement interpolation with respect to incident energy" );
      }
    }
  }
}
// ----------- Edist::energy_moments::one_Ein_to_ENDF --------------
//! Converts ENDL data to ENDF format for one incident energy
void Edist::energy_moments::one_Ein_to_ENDF( int Ein_count )
{
  int L_order;
  // make a new StdLg::standard_Legendre_vector for this incident energy
  StdLg::standard_Legendre::iterator ENDF_Ein_ptr;
  StdLg::standard_Legendre_vector::iterator ENDF_Eout_ptr;
  // where we are in the ENDL data
  Ebase::Eprob_vector::const_iterator *ENDL_Eout_ptr =
    new Ebase::Eprob_vector::const_iterator[ data_order + 1 ];
  Ebase::Eprob_vector::const_iterator *next_ENDL_ptr =
    new Ebase::Eprob_vector::const_iterator[ data_order + 1 ];
  for( L_order = 0; L_order <= data_order; ++L_order )
  {
    ENDL_Eout_ptr[ L_order ] = ENDL_data[ L_order ].EProb_data[ Ein_count ].begin( );
    next_ENDL_ptr[ L_order ] = ENDL_Eout_ptr[ L_order ];
    ++next_ENDL_ptr[ L_order ];
  }
  ENDF_Ein_ptr = ENDF_data.insert( ENDF_data.end( ), StdLg::standard_Legendre_vector( ) );
  ENDF_Ein_ptr->set_E_in( ENDL_data[ 0 ].EProb_data[ Ein_count ].get_E_in( ) );
  ENDF_Ein_ptr->Eout_interp = Eout_interp;
  if( Ein_interp.qualifier == Terp::UNITBASE )
  {
    ENDF_Ein_ptr->ubase_map.copy( ENDL_data[ 0 ].EProb_data[ Ein_count ].ubase_map );
  }
  // go through the ENDL data
  bool done = false;
  double next_Eout = ENDL_Eout_ptr[ 0 ]->x;  // this has to be zero
  for( ; ; )
  {
    // make a new set of Legendre coefficients
    ENDF_Eout_ptr = ENDF_Ein_ptr->insert( ENDF_Ein_ptr->end( ), Lgdata::Legendre_coefs( ) );
    ENDF_Eout_ptr->initialize( output_order );
    ENDF_Eout_ptr->set_E_out( next_Eout );  // energy of outgoing particle
    for( L_order = 0; L_order <= data_order; ++L_order )
    {
      if( Ein_interp.qualifier == Terp::UNITBASE )
      {
	bool is_OK;
        ( *ENDF_Eout_ptr )[ L_order ] = 
  	  ENDL_Eout_ptr[ L_order ]->linlin_interp( next_Eout,
			   *next_ENDL_ptr[ L_order ], &is_OK );
	if( !is_OK )
	{
	  Msg::FatalError( "Edist::energy_moments::one_Ein_to_ENDF",
			   "bad interpolation" );
	}
      }
      else
      {
	( *ENDF_Eout_ptr )[ L_order ] = ENDL_Eout_ptr[ L_order ]->y;
      }
    }
    if( done )
    {
      break;
    }
    // get the next outgoing energy
    static double tol = Global.Value( "tight_tol" );
    double huge_E = 1.0e20;  // a dummy large value
    next_Eout = huge_E;
    for( L_order = 0; L_order <= data_order; ++L_order )
    {
      if( next_ENDL_ptr[ L_order ]->x < ( 1.0 - tol )*next_Eout )
      {
	next_Eout = next_ENDL_ptr[ L_order ]->x;
      }
    }
    // update the pointers
    if( next_Eout > 1.0 - tol )
    {
      done = true;
    }
    else
    {
      for( L_order = 0; L_order <= data_order; ++L_order )
      {
        if( next_ENDL_ptr[ L_order ]->x < ( 1.0 + tol )*next_Eout )
        {
	  ENDL_Eout_ptr[ L_order ] = next_ENDL_ptr[ L_order ];
	  ++next_ENDL_ptr[ L_order ];
	}
      }
    }
  }
  delete [] ENDL_Eout_ptr;
  delete [] next_ENDL_ptr;
}
// ----------- Edist::energy_moments::get_T --------------
// Calculates the transfer matrix for this particle.
// sigma is the cross section.
void Edist::energy_moments::get_T( const Ddvec::dd_vector& sigma,
  const Ddvec::dd_vector& mult, const Ddvec::dd_vector& weight,
				   Trf::T_matrix& transfer )
{
  bool interp_OK = ( ( Ein_interp.qualifier == Terp::UNITBASE ) && 
                     ( Ein_interp.flag == Terp::LINLIN ) ) ||
    ( ( Ein_interp.qualifier == Terp::DIRECT ) &&
      ( Ein_interp.flag == Terp::HISTOGRAM ) );
  if( !interp_OK )
  {
    Msg::FatalError( "Edist::energy_moments::get_T",
      "Incident energy interpolation not implemented" );
  }
  interp_OK = ( Eout_interp == Terp::LINLIN ) || ( Eout_interp == Terp::HISTOGRAM );
  if( !interp_OK )
  {
    Msg::FatalError( "Edist::energy_moments::get_T",
      "Outgoing energy interpolation not implemented" );
  }
  output_order = transfer.order;
  if( data_order > output_order )
  {
    data_order = output_order;  // truncate the data
  }
  // convert to ENDF format
  if( output_order == 0 )
  {
    zero_order( );
  }
  else
  {
    to_ENDF( );
  }
  ENDF_data.Ein_interp = Ein_interp;
  ENDF_data.Eout_interp = Eout_interp;
  ENDF_data.get_T( sigma, mult, weight, transfer );
}
