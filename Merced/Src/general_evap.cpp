/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: general_evap.cpp 1 2010-06-18 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 * 
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */
//! Implementation of classes used for outgoing an energy distribution independent of the initial energy

#ifdef _OPENMP
 #include <omp.h>
#endif

#include "general_evap.hpp"
#include "adapt_quad.hpp"
#include "messaging.hpp"
#include "global_params.hpp"

// **************** class Gevap::general_evap_param ********************
// ---------------- Gevap::general_evap_param::get_integrals ------------------
// Gets the integrals over this E_out bin
void Gevap::general_evap_param::get_integrals( Coef::coef_vector *value )
{
  *value = *current_Eout_int;
}

// **************** class Gevap::general_evap ********************
// ----------- Gevap::general_evap::get_Ein_range --------------
//  Gets the range of nontrivial incident energy bins; computes first_Ein and last_Ein
// returns true if the threshold is too high for the energy bins
bool Gevap::general_evap::get_Ein_range( const Ddvec::dd_vector& sigma, const Ddvec::dd_vector& multiple,
    const Ddvec::dd_vector& weight,
    const Lgdata::Flux_List& e_flux, const Egp::Energy_groups& Ein_groups )
{
  double E_last;

  Gevap::general_evap_param initial_param;
  bool done = initial_param.get_Ein_range( sigma, multiple, weight, e_flux,
                                         Ein_groups, &E_first, &E_last );
  if( done ) return true;

  first_Ein = Ein_groups.first_bin_ID( E_first );
  last_Ein = Ein_groups.last_bin_ID( E_last );
  return false;
}
// ---------------- general_evap::setup_Eout_ints ------------------
// Sets up the integrals over the E_out bins
void Gevap::general_evap::setup_Eout_ints( int order, Coef::Conserve conserve,
  const vector< double >& Eout_groups )
{
  // ***** CODING ASSUMPTIONS *****
  // It is assumed that the probability density is given as a histogram.
  // It is also assumed that the probability density is independent of incident energy.
  // *****************************

  Coef::coef_vector new_link;
  E_out_ints.push_back( new_link );
  std::list< Coef::coef_vector >::iterator new_entry = E_out_ints.end( );
  --new_entry;   // points to the new Coef::coef_vector
  new_entry->set_order( order, conserve );  // identify space for the integrals
  vector< double >::const_iterator Eout_0 = Eout_groups.begin( );
  vector< double >::const_iterator Eout_1 = Eout_0;
  ++Eout_1;
  Ddvec::dd_vector::const_iterator data_0 = begin( );
  Ddvec::dd_vector::const_iterator data_1 = data_0;
  ++data_1;

  // synchronize the data
  while( data_1->x <= *Eout_0 )
  {
    data_0 = data_1;
    ++data_1;
    if( data_1 == end( ) )
    {
      Msg::FatalError( "Gevap::general_evap::setup_Eout_ints",
		       "All data below the lowest bin" );
    }
  }
  while( *Eout_1 <= data_0->x )
  {
    Eout_0 = Eout_1;
    ++Eout_1;
    if( Eout_1 == Eout_groups.end( ) )
    {
      Msg::FatalError( "Gevap::general_evap::setup_Eout_ints",
		       "All data above the highest bin" );
    }
  }

  for( ; ; )
  {
    double E_0 = ( data_0->x < *Eout_0 ) ? *Eout_0 : data_0->x;
    double E_1 = ( data_1->x > *Eout_1 ) ? *Eout_1 : data_1->x;
    double probability = data_0->y;
    double dE = E_1 - E_0;
    if( ( conserve == Coef::NUMBER ) || ( conserve == Coef::BOTH ) )
    {
      new_entry->weight_1[ 0 ] += probability * dE;
    }
    if( ( conserve == Coef::ENERGY ) || ( conserve == Coef::BOTH ) )
    {
      new_entry->weight_E[ 0 ] += probability * dE * ( E_0 + E_1 ) / 2;
    }
    if( data_1->x < *Eout_1 )
    {
      data_0 = data_1;  // go to the next data interval
      ++data_1;
      if( data_1 == end( ) )
      {
	break;  // no more data
      }
    }
    else if( data_1->x > *Eout_1 )
    {
      Eout_0 = Eout_1;  // go to the next energy bin
      ++Eout_1;
      if( Eout_1 == Eout_groups.end( ) )
      {
	break;  // no more energy groups
      }
      E_out_ints.push_back( new_link );
      new_entry = E_out_ints.end( );
      --new_entry;   // points to the new Coef::coef_vector
      new_entry->set_order( order, conserve );  // identify space for the integrals
    }
    else if( data_1->x > *Eout_1 )
    {
      data_0 = data_1;  // go to the next data interval
      ++data_1;
      if( data_1 == end( ) )
      {
	break;  // no more data
      }
      Eout_0 = Eout_1;  // go to the next energy bin and the next data interval
      ++Eout_1;
      if( Eout_1 == Eout_groups.end( ) )
      {
	break;  // no more energy groups
      }
      E_out_ints.push_back( new_link );
      new_entry = E_out_ints.end( );
      --new_entry;   // points to the new Coef::coef_vector
      new_entry->set_order( order, conserve );  // identify space for the integrals
    }
  }
  check_sum( );
}
// ---------------- Gevap::general_evap::check_sum ------------------
// The probability integals for outgoing energy should add to 1
void Gevap::general_evap::check_sum( )
{
  double sum = 0.0;
  static double quad_tol = Global.Value( "quad_tol" );
  for( std::list< Coef::coef_vector >::const_iterator one_term = E_out_ints.begin( );
       one_term != E_out_ints.end( ); ++one_term )
  {
    sum += one_term->weight_1[ 0 ];
  }
  if( std::abs( sum - 1.0 ) > quad_tol )
  {
    Msg::Warning( "Gevap::general_evap::check_sum",
           Msg::pastenum( "The norm should be 1; it is ", sum ) );
  }
}
// ---------------- Gevap::general_evap::next_ladder ------------------
// Sets up the integrals over the E_out bins
// Go to the next pair of incident energies.  Returns "true" when finished.
bool Gevap::general_evap::next_ladder( double E_in, Gevap::general_evap_param *Ein_param )
{
  bool done = Ein_param->update_bin_pointers( E_in );
  return done;
}
// ---------------- Gevap::general_evap::Eout_ladder ------------------
// Sets up the integrals over the E_out bins
// Adds to the transfer matrix for all E_out bins for a pair of incident energies.
void Gevap::general_evap::Eout_ladder( Trf::T_matrix& transfer,
				       Gevap::general_evap_param *Ein_param )
{
  Ein_param->current_Eout_int = E_out_ints.begin( );
  vector< double >::const_iterator Eout_ptr = transfer.out_groups.begin( );
  unsigned int Eout_count = 0;
  for( ; Eout_count < E_out_ints.size( );
       ++Ein_param->current_Eout_int, ++Eout_ptr, ++Eout_count )
  {
    // integrate over this E-E' box
    update_T( transfer, Eout_count, Ein_param );
  }
}
// ---------------- Gevap::general_evap::update_T ------------------
// Sets up the integrals over the E_out bins
// Adds to an element of transfer the integral over the E-E' box
void Gevap::general_evap::update_T( Trf::T_matrix &transfer, int Eout_count,
     Gevap::general_evap_param *Ein_param )
{
  // a vector to store the integrals
  Coef::coef_vector value( transfer.order, transfer.conserve );

  // parameters for the integration
  Qparam::QuadParamBase *params = static_cast< Qparam::QuadParamBase* >( Ein_param );

  // loop over the cross section data
  Ein_param->this_sigma = Ein_param->first_ladder_sigma;
  Ein_param->next_sigma = Ein_param->this_sigma;
  ++Ein_param->next_sigma;
  // Ein_param->data_E_0 may be past Ein_param->next_sigma
  while( ( Ein_param->this_sigma != Ein_param->last_ladder_sigma ) &&
         ( Ein_param->next_sigma->x < Ein_param->data_E_0 ) )
  {
    Ein_param->this_sigma = Ein_param->next_sigma;
    ++Ein_param->next_sigma;
  }
  for( ; ( Ein_param->this_sigma != Ein_param->last_ladder_sigma ) &&
         ( Ein_param->this_sigma->x <  Ein_param->data_E_1 );
       Ein_param->this_sigma = Ein_param->next_sigma, ++Ein_param->next_sigma )
  {
    double left_E = ( Ein_param->this_sigma->x < Ein_param->data_E_0 ) ? Ein_param->data_E_0 :
      Ein_param->this_sigma->x;
    double right_E = ( Ein_param->next_sigma->x > Ein_param->data_E_1 ) ? Ein_param->data_E_1 :
      Ein_param->next_sigma->x;
    static double tol = Global.Value( "quad_tol" );  // the default tolerance
    // evaluate the integral
    quad_F::integrate( general_evap_F::Ein_F, transfer.Ein_quad_rule, left_E,
		       right_E, params, tol, &value );

    if( value.weight_1[ 0 ] < 0.0 )
    {
      Msg::Warning( "energy_function::update_T",
		    Msg::pastenum( "negative integral ", left_E ) +
               Msg::pastenum(" ", right_E ) );
      value.set_zero( );  // throw out these values
    }

    // add this integral
    transfer( Ein_param->Ein_count, Eout_count ) += value;
    // increment the function counts
    Ein_param->Ein_F_count += Ein_param->func_count;
    Ein_param->quad_count += 1;
  }
}
// ---------------- Gevap::general_evap::theta_OK ------------------
// The code currently handles only theta = constant
bool Gevap::general_evap::theta_OK( )
{
  bool is_OK = true;
  static double abs_tol = Global.Value( "tight_tol" );
  Ddvec::dd_vector::const_iterator this_theta = theta.begin( );
  double first_theta = this_theta->y;

  for( ++this_theta; this_theta != theta.end( ); ++this_theta )
  {
    if( std::abs( this_theta->y - first_theta ) > abs_tol*first_theta )
    {
      is_OK = false;
      break;
    }
  }
  return is_OK;
}
// ---------------- Gevap::general_evap::get_T ------------------
// Sets up the integrals over the E_out bins
// Calculates the transfer matrix for this particle
void Gevap::general_evap::get_T( const Ddvec::dd_vector& sigma, const Ddvec::dd_vector& multiple, 
    const Ddvec::dd_vector& weight, Trf::T_matrix& transfer )
{
  if( interp_type != Terp::LINLIN )
  {
    Msg::FatalError( "Gevap::general_evap::get_T",
		     "interp_type not implemented" );
  }
  if( !theta_OK( ) )
  {
    Msg::FatalError( "Gevap::general_evap::get_T",
		     "only theta = constant is implemented" );
  }

  // Scale the x-data by 1/Theta and renorm
  scale_E( theta.begin( )->y );
  renorm( false );

  bool done = get_Ein_range( sigma, multiple, weight, transfer.e_flux,
    transfer.in_groups );
  if( done )
  {
    transfer.zero_transfer( );
  }

  // Set up the integrals over outgoing energy
  setup_Eout_ints( transfer.order, transfer.conserve, transfer.out_groups );
  transfer.threshold = sigma.begin( )->x;

  long int quad_count = 0;  // number of 2-d quadratures
  long int Ein_F_count = 0;  // number of calls to general_evap_F::E_quad_F

  // do the integrals incident bin by incident bin
#pragma omp parallel for schedule( dynamic, 1 ) default( none )	\
  shared( sigma, multiple, weight, transfer ) \
  reduction( +: quad_count ) reduction( +: Ein_F_count )
  // now do the integrals over the incident energy bins
  for( int Ein_bin = first_Ein; Ein_bin < last_Ein; ++Ein_bin )
  {
    Gevap::general_evap_param Ein_param;
    // set up the data range for this bin
    Ein_param.setup_bin( Ein_bin, sigma, multiple, weight, transfer.e_flux,
                         transfer.in_groups );
    // work on this bin
    for( ; ; )
    {
      Ein_param.set_Ein_range( );   // get the incident energy interval
      Ein_param.set_sigma_range( );
      Eout_ladder( transfer, &Ein_param );
      bool Done = next_ladder( Ein_param.data_E_1, &Ein_param );   // go to the next interval
      if( Done )
      {
        break;
      }
    }
    quad_count += Ein_param.quad_count;
    Ein_F_count += Ein_param.Ein_F_count;
  } // end of parallel loop

  // print the counts of function evaluations
  std::cout << "1d quadratures: " << quad_count << std::endl;
  std::cout << "general_evap_F::Ein_F calls: " << Ein_F_count << std::endl;
  std::cout << "average general_evap_F::Ein_F calls: " << 1.0*Ein_F_count/quad_count << std::endl;
 }

// ************** general_evap_F::Ein_F ******************************
// Integral function for the model
bool general_evap_F::Ein_F( double E_in, Qparam::QuadParamBase *Ein_param,
   Coef::coef_vector *value )
{
  // the parameters are really Gevap::general_evap_param *
  Gevap::general_evap_param *e_params = static_cast<Gevap::general_evap_param *>( Ein_param );
  e_params->func_count += 1;

  e_params->get_integrals( value );

  // weight it by flux * cross section * multiplicity * model_weight
  e_params->set_weight( E_in );
  *value *= e_params->current_weight;

  return true;
}
