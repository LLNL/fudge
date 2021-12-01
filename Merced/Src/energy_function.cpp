/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: energy_function.cpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */
//! Implementation of the classes used to handle formulas for energy probability density.
//! All data is in laboratory coordinates.

#include "energy_function.hpp"
#include "adapt_quad.hpp"
#include "messaging.hpp"
#include "global_params.hpp"

// ************* class Efunc::E_function_param *****************
// ---------------- Efunc::E_function_param::set_Ein_default --------------------
// Interpolate the parameters
bool Efunc::E_function_param::set_Ein_default( double E_in )
{
  bool Theta_OK;
  Theta = this_Theta->linlin_interp( E_in, *next_Theta, &Theta_OK );

  bool mult_OK;
  multiplicity = this_mult->linlin_interp( E_in, *next_mult, &mult_OK );
  if( ( !Theta_OK ) || ( !mult_OK ) )
  {
    Msg::DebugInfo( "Efunc::E_function_param::set_Ein_default",
		     "got a negative parameter" );
    return false;
  }
  // the range of integration
  E_max = E_in - U;
  if( E_max > top_E_out )
  {
    E_max = top_E_out;
  }
  E_1 = ( use_Eout_max ) ? Eout_max : E_max;
  set_scales( );
  norm = get_norm( );

  return true;
}

// ************* class Efunc::U_Ein_hit_list *****************
// ----------- Efunc::U_Ein_hit_list::find_bottom_hits --------------
// Finds the intersection with the bottom of a box
void Efunc::U_Ein_hit_list::find_bottom_hits( double E_out,
  std::vector< Box::Ein_Eta_Hit > *Ein_hits )
{
  // for new entries
  Box::Ein_Eta_Hit Ein_eta_hit;
  double Ein = E_out + get_U( );

  if( Ein > 0.0 )
  {
    // append this entry
    Ein_eta_hit.E_in = Ein;
    Ein_eta_hit.hit_edge = Box::BOTTOM_IN;
    Ein_hits->push_back( Ein_eta_hit );
  }
}
// ----------- Efunc::U_Ein_hit_list::find_top_hits --------------
// Finds the intersection with the top of a box
void Efunc::U_Ein_hit_list::find_top_hits( double E_out,
  std::vector< Box::Ein_Eta_Hit > *Ein_hits )
{
  // for new entries
  Box::Ein_Eta_Hit Ein_eta_hit;
  double Ein = E_out + get_U( );

  if( Ein > 0.0 )
  {
    // append this entry
    Ein_eta_hit.E_in = Ein;
    Ein_eta_hit.hit_edge = Box::TOP_OUT;
    Ein_hits->push_back( Ein_eta_hit );
  }
}


// ************* class Efunc::energy_function *****************
// ---------------- Efunc::energy_function::setup_data_default --------------------
// Initializes the quadrature parameters
void Efunc::energy_function::setup_data_default( const Egp::Energy_groups& Eout_groups,
  Efunc::E_function_param *Ein_param )
{
  static double skip_tol = Global.Value( "tight_tol" );

  Ein_param->U = U;
  Ein_param->top_E_out = Eout_groups.get_top_E( );
  Ein_param->this_Theta = begin( );
  Ein_param->next_Theta = Ein_param->this_Theta;
  ++Ein_param->next_Theta;
  while( Ein_param->next_Theta->x < (*Ein_param->Ein_ptr) *
	 ( 1.0 + skip_tol ) )
  {
    Ein_param->this_Theta = Ein_param->next_Theta;
    ++Ein_param->next_Theta;
  }
  Ein_param->Theta_end = end( );

  Ein_param->upper_hits.set_U( Ein_param->U );
  double first_E = Ein_param->this_Theta->x;
  if( first_E > Ein_param->data_E_0 )
  {
    Ein_param->data_E_0 = first_E;
    bool data_bad = Ein_param->update_pointers( first_E );
    if( data_bad )
    {
      Msg::FatalError( "Efunc::energy_function::setup_data_default",
		       "energies inconsistent" );
    }
  }
}
// ---------------- Efunc::energy_function::set_Ein_range_default --------------------
// Sets the range of incident energies for this intergration
void Efunc::energy_function::set_Ein_range_default( int Ein_bin,
			     Efunc::E_function_param *Ein_param )
{
  Ein_param->set_Ein_range( );
  double this_E = Ein_param->this_Theta->x;
  if( this_E > Ein_param->data_E_0 ) Ein_param->data_E_0 = this_E;
  this_E = Ein_param->next_Theta->x;
  if( this_E < Ein_param->data_E_1 ) Ein_param->data_E_1 = this_E;

  if( Ein_param->data_E_1 < Ein_param->data_E_0 )
  {
    Msg::FatalError( "Efunc::energy_function::set_Ein_range_default",
		     "check the Theta incident energies" );
  }
}
// ---------------- Efunc::energy_function::next_ladder_default --------------------
// Default go to the next (incident energy, Theta).  Returns "true" when finished.
bool Efunc::energy_function::next_ladder_default( double E_in,
				      Efunc::E_function_param *Ein_param )
{
  bool done = Ein_param->update_bin_pointers( E_in );
  static double etol = Global.Value( "tight_tol" );
  double E_tol = E_in * etol;
  //    double E_tol = 0.0;
  if( !done )
  {
    if( E_in + E_tol >= Ein_param->next_Theta->x )
    {
      while( E_in + E_tol >= Ein_param->next_Theta->x )
      {
        // get the next (E_in, Theta) data
        Ein_param->this_Theta = Ein_param->next_Theta;
        ++Ein_param->next_Theta;
        if( Ein_param->next_Theta == end( ) )
        {
          return true;
        }
      }
    }
  }
  return done;
}
// ----------- Efunc::energy_function::Eout_ladder --------------
// This routine uses the angular distributions this_eta_dist and the
// next to calculate the contribution to the E_out boxes of the
// transfer matrix between incident energies Ein_param->data_E_0 and
// Ein_param->data_E_1.
void Efunc::energy_function::Eout_ladder( Trf::T_matrix& transfer,
					  Efunc::E_function_param *Ein_param )
{
  bool geom_OK;  // for checking the consistency of the geometry
  Efunc::U_Ein_hit_list test_hits;
  // loop through the outgoing energies (column of transfer)
  for( int Eout_count = 0; Eout_count < transfer.num_Eout_bins;
    ++Eout_count )
  {
    std::vector< double >::const_iterator Eout_ptr = transfer.out_groups.begin( )
      + Eout_count;
    // how does the line Eout = Ein - U meet this E-E' box?
    double U = Ein_param->U;
    geom_OK = Ein_param->upper_hits.hit_box( U, Eout_ptr,
      Ein_param->data_E_0, Ein_param->data_E_1 );
    if( !geom_OK )
    {
      test_hits.set_U( U );
      test_hits.hit_box( U, Eout_ptr,
                         Ein_param->data_E_0, Ein_param->data_E_1 );
      test_hits.print( );
      Msg::FatalError( "Efunc::energy_function::Eout_ladder",
		       "Check the coding" );
    }
    if( ( Eout_count > 0 ) && ( Ein_param->upper_hits.is_below( ) ) )
    {
      // we are done with this incident energy bin
      break;
    }
    // integrate over this E-E' box
    one_Ebox( transfer, Eout_count, Ein_param );
  }
}
// ----------- Efunc::energy_function::one_Ebox --------------
// Does the integration for one E-E' box
void Efunc::energy_function::one_Ebox( Trf::T_matrix& transfer, int Eout_count,
   Efunc::E_function_param *Ein_param )
{
  static double skip_tol = Global.Value( "tight_tol" );
  // the E' energy range
  Ein_param->Eout_min = transfer.out_groups[ Eout_count ];
  Ein_param->Eout_max = transfer.out_groups[ Eout_count + 1 ];

  // integrate depending on how the line Eout = Ein - U meets the box
  Efunc::U_Ein_hit_list::iterator high_hit_ptr = Ein_param->upper_hits.begin( );
  Efunc::U_Ein_hit_list::iterator next_high_ptr = high_hit_ptr;
  ++next_high_ptr;
  for( ; next_high_ptr != Ein_param->upper_hits.end( );
         high_hit_ptr = next_high_ptr, ++next_high_ptr )
  {
    // always integrate from Eout_min
    Ein_param->use_Eout_min = true;
    // where is the line Eout = Ein - U?
    if( high_hit_ptr->hit_edge == Box::BELOW )
    {
      // do nothing---we are below the E-E' box
      continue;
    }
    else if( ( high_hit_ptr->hit_edge == Box::BOTTOM_IN ) ||
             ( high_hit_ptr->hit_edge == Box::INSIDE ) )
    {
      // the line Eout = Ein - U is inside the E-E' box
      Ein_param->use_Eout_max = false;
    }
    else
    {
      // the line Eout = Ein - U is above the E-E' box
      Ein_param->use_Eout_max = true;
    }
    // the range of integration in incident energy
    Ein_param->Ein_0 = high_hit_ptr->E_in;
    Ein_param->Ein_1 = next_high_ptr->E_in;
    if( Ein_param->Ein_1 - Ein_param->Ein_0 <= Ein_param->Ein_1 * skip_tol )
    {
      Msg::Warning( "Efunc::energy_function::one_Ebox",
		    "skipping a very short interval");
      continue;  // skip this interval
    }

    update_T( transfer, Eout_count, Ein_param );
  }
}
// ----------- Efunc::energy_function::update_T --------------
// Adds to an element of transfer the integral between over the E-E' box
void Efunc::energy_function::update_T( Trf::T_matrix &transfer, int Eout_count,
   Efunc::E_function_param *Ein_param )
{
  // a vector to store the integrals
  Coef::coef_vector value( transfer.order, transfer.conserve );
  value.set_zero( );

  // parameters for the integration
  Qparam::QuadParamBase *params = static_cast< Qparam::QuadParamBase* >( Ein_param );

  // loop over the cross section data
  Ein_param->this_sigma = Ein_param->first_ladder_sigma;
  Ein_param->next_sigma = Ein_param->this_sigma;
  ++Ein_param->next_sigma;
  // Ein_param->Ein_0 may be past Ein_param->next_sigma
  while( ( Ein_param->this_sigma != Ein_param->last_ladder_sigma ) &&
         ( Ein_param->next_sigma->x < Ein_param->Ein_0 ) )
  {
    Ein_param->this_sigma = Ein_param->next_sigma;
    ++Ein_param->next_sigma;
  }
  for( ; ( Ein_param->this_sigma != Ein_param->last_ladder_sigma ) &&
         ( Ein_param->this_sigma->x <  Ein_param->Ein_1 );
       Ein_param->this_sigma = Ein_param->next_sigma, ++Ein_param->next_sigma )
  {
    double left_E = ( Ein_param->this_sigma->x < Ein_param->Ein_0 ) ? Ein_param->Ein_0 :
      Ein_param->this_sigma->x;
    double right_E = ( Ein_param->next_sigma->x > Ein_param->Ein_1 ) ? Ein_param->Ein_1 :
      Ein_param->next_sigma->x;
    double tol = Ein_param->set_tol( left_E, right_E );

    quad_F::integrate( Energy_function_F::Ein_F, transfer.Ein_quad_rule, left_E,
		       right_E, params, tol, &value );

    if( value.weight_1[ 0 ] < 0.0 )
    {
      Msg::Warning( "Efunc::energy_function::update_T",
		    Msg::pastenum( "negative integral ", left_E ) +
	       Msg::pastenum(" ", right_E ) );
      value.set_zero( );  // throw out these values
    }

    // add this integral
    transfer( Ein_param->Ein_count, Eout_count ) += value;
    // increment the function counts
    Ein_param->Ein_F_count += Ein_param->func_count;
    ++Ein_param->quad_count;
  }
}

// ************** Probability density model ******************************
// ----------- Energy_function_F::Ein_F --------------
//! Integral function for the model
bool Energy_function_F::Ein_F( double E_in, Qparam::QuadParamBase *Ein_param,
   Coef::coef_vector *value )
{
  // the parameters are really Efunc::E_function_param *
  Efunc::E_function_param *e_params = static_cast<Efunc::E_function_param *>( Ein_param );
  e_params->func_count += 1;
  bool Ein_OK = e_params->set_Ein( E_in );  // interpolate the data

  bool integral_OK = e_params->get_integrals( e_params->Eout_min, e_params->E_1, *value );

  if( !Ein_OK || !integral_OK )
  {
    return false;
  }
  // weight it by flux * cross section * multiplicity * model_weight
  e_params->set_weight( E_in );
  *value *= e_params->current_weight;

  return true;
}
