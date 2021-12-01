/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: Watt.cpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 * 
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */
// Implement classes used for Watt energy probability density

#include <cmath>
#include <cfloat>
#ifdef _OPENMP
 #include <omp.h>
#endif

#include "Watt.hpp"
#include "global_params.hpp"
#include "messaging.hpp"
#include "adapt_quad.hpp"

#include "protos.h"

// ****************** class Wtt::Watt_param **********************
// ---------------- Wtt::Watt_param::set_Ein --------------------
bool Wtt::Watt_param::set_Ein( double E_in )
// Interpolates the data to the given incident energy
{
  // We need to set b first, because set_Ein_default uses it in calculating the norm.
  bool b_OK;
  b = this_b->linlin_interp( E_in, *next_b, &b_OK );
  bool a_OK = set_Ein_default( E_in );

  return( a_OK && b_OK );
}
// ---------------- Wtt::Watt_param::get_integrals --------------------
// Gets the integrals over outgoing energy
bool Wtt::Watt_param::get_integrals( double Eout_0, double Eout_1,
				     Coef::coef_vector &value )
{
  double numerator;

  if( ( value.conserve == Coef::NUMBER ) || ( value.conserve == Coef::BOTH ) )
  {
    // integrate probability density
    conserve_flag = Coef::NUMBER;
    numerator = integral_prob( Eout_0, Eout_1 );
    value.weight_1[ 0 ] = numerator/norm;
  }
  if( ( value.conserve == Coef::ENERGY ) || ( value.conserve == Coef::BOTH ) )
  {
    // integrate E*probability density
    //    conserve_flag = Coef::ENERGY;
    //    numerator = integral_Eprob( Eout_0, Eout_1 );
    //    value.weight_E[ 0 ] = numerator/norm;
    Msg::FatalError( "Watt_param::get_integrals",
		     "value.conserve == Coef::ENERGY not implemented" );
  }
  return true;
}
// ---------------- Wtt::Watt_param::set_scales --------------------
//! Sets the scale factors for the integrals of probability and energy*probability
void Wtt::Watt_param::set_scales( )
{
  // the scale factors used with the incomplete gamma functions
  prob_integral_scale = 1.0; // not used
  Eprob_integral_scale = 1.0; // not used
}
// ---------------- Wtt::Watt_param::get_norm --------------------
// Integrate from 0 to E_max to get the norm
double Wtt::Watt_param::get_norm( )
{
  return integral_prob( 0.0, E_max );
}
// ---------------- Wtt::Watt_param::integral_prob --------------------
// Indefinite integral of the probability density
double Wtt::Watt_param::integral_prob( double E_0, double E_1 )
{
  // The factor a*exp{beta^2} cancels, so ignore it.
  // For compatibility with other models use Theta = a.
  double beta = sqrt( b*Theta )/2.0;
  double alpha0_minus = sqrt( E_0/Theta ) - beta;
  double alpha0_plus = alpha0_minus + 2*beta;
  double alpha1_minus = sqrt( E_1/Theta ) - beta;

  double alpha1_plus = alpha1_minus + 2*beta;
  double integral;
  double sinh_term;

  // what we do depends on the sizes of the alphas
  if( alpha0_minus >= 0.0 )
  {
    integral = beta * ( one_integral( alpha0_minus, alpha1_minus ) +
			      one_integral( alpha0_plus, alpha1_plus ) );
    if( alpha1_minus <= alpha0_plus )
    {
      sinh_term = w_integral( alpha0_minus, alpha1_minus ) -
                  w_integral( alpha0_plus, alpha1_plus );
      integral += sinh_term;
    }
    else
    {
      integral += w_integral( alpha0_minus, alpha0_plus ) - w_integral( alpha1_minus, alpha1_plus );
    }
  }
  else if( alpha1_minus >= 0.0 )
  {
    integral = beta * ( one_integral( 0.0, -alpha0_minus ) + 
                              one_integral( 0.0, alpha1_minus ) +
			      one_integral( alpha0_plus, alpha1_plus ) );
    if( alpha1_minus <= alpha0_plus )
    {
      sinh_term = w_integral( -alpha0_minus, alpha1_minus ) -
                  w_integral( alpha0_plus, alpha1_plus );
      integral += sinh_term;
    }
    else
    {
      sinh_term = w_integral( -alpha0_minus, alpha0_plus ) - w_integral( alpha1_minus, alpha1_plus );
      integral += sinh_term;
    }
  }
  else
  {
    integral = beta * ( one_integral( -alpha1_minus, -alpha0_minus ) +
			      one_integral( alpha0_plus, alpha1_plus ) );
    sinh_term = w_integral( -alpha1_minus, -alpha0_minus ) + w_integral( alpha0_plus, alpha1_plus );
    integral -= sinh_term;
  }

  // do a check
  //  double zz = check_Watt( E_0, E_1 );
  return integral;
}
// ---------------- Wtt::Watt_param::integral_Eprob --------------------
// Indefinite integral of energy times the probability density
double Wtt::Watt_param::integral_Eprob( double E_0, double E_1 )
{
  return 0.0;
}
// ---------------- Wtt::Watt_param::w_integral -----------------------
// A robust integration of \int_{\alpha_0}^{\alpha_1} w \exp{-w^2} \, dw
double Wtt::Watt_param::w_integral( double alpha_0, double alpha_1 )
{
  double integral;
  if( alpha_0 > 1.0 )
  {
    integral = Proto::igamc( 1.0, alpha_0*alpha_0 ) - Proto::igamc( 1.0, alpha_1*alpha_1 );
  }
  else
  {
    integral = Proto::igam( 1.0, alpha_1*alpha_1 ) - Proto::igam( 1.0, alpha_0*alpha_0 );
  }
  // scale igam by Tgamma( 1 ) = 1
  return 0.5 * integral;
}
// ---------------- Wtt::Watt_param::one_integral -----------------------
// A robust integration of \int_{\alpha_0}^{\alpha_1} \exp{-w^2} \, dw
double Wtt::Watt_param::one_integral( double alpha_0, double alpha_1 )
{
  double integral;
  if( alpha_0 > 1.0 )
  {
    integral = Proto::igamc( 0.5, alpha_0*alpha_0 ) - Proto::igamc( 0.5, alpha_1*alpha_1 );
  }
  else
  {
    integral = Proto::igam( 0.5, alpha_1*alpha_1 ) - Proto::igam( 0.5, alpha_0*alpha_0 );
  }
  // scale igam by Tgamma( 1/2 )
  return 0.5 * tgamma( 0.5 ) * integral;
}
// ----------- Wtt::Watt_param::set_tol --------------
// Sets the tolerance for the quadrature over incident energy
double Wtt::Watt_param::set_tol( double left_E, double right_E )
{
  static double tol = Global.Value( "quad_tol" ); // the default tolerance
  return tol;
}
// ---------------- Wtt::Watt_param::check_Watt --------------------
// Do the integration by adaptive quadrature
double Wtt::Watt_param::check_Watt( double E_0, double E_1 )
{
  static double tol = Global.Value( "quad_tol" );
  // a vector to store the integrals
  Coef::coef_vector value( 0, Coef::NUMBER );  // Legendre order 0, conserve number
  Wtt::check_Watt_params check_Watt_param;
  check_Watt_param.a = Theta;
  check_Watt_param.b = b;
  // parameters for the integration
  Qparam::QuadParamBase *params =
    static_cast< Qparam::QuadParamBase* >( &check_Watt_param );

  // the linked list for adaptive quadrature
  quad_F::integrate( Watt_F::check_Watt_F, Eout_quad_rule, E_0, E_1,
		     params, tol, &value );

  Eout_F_count += check_Watt_param.func_count;
  return value.weight_1[ 0 ]/( Theta*exp(b*Theta/4) );
}

// ****************** class Wtt::watt **********************
// ---------------- Wtt::watt::setup_data --------------------
// Initializes the quadrature parameters
void Wtt::watt::setup_data( const Egp::Energy_groups& Eout_groups,
   Efunc::E_function_param *ein_param )
{
  static double skip_tol = Global.Value( "tight_tol" );

  // the parameters are really Wtt::Watt_param *
  Wtt::Watt_param *Ein_param = static_cast<Watt_param *>( ein_param );

  setup_data_default( Eout_groups, Ein_param );
  Ein_param->this_b = b_data.begin( );
  Ein_param->next_b = Ein_param->this_b;
  ++Ein_param->next_b;

  while( ( Ein_param->next_b->x < E_first * ( 1.0 + skip_tol ) ) ||
	 ( Ein_param->next_b->x < (*Ein_param->Ein_ptr) *
           ( 1.0 + skip_tol ) ) )
  {
    Ein_param->this_b = Ein_param->next_b;
    ++Ein_param->next_b;
  }
  Ein_param->b_end = b_data.end( );

  // Set the maximum energy of fission neutrons "large"
  //  upper_hits.set_U( -Ein_param->top_E_out );
  //  Ein_param->U = -Ein_param->top_E_out;
}
// ----------- Wtt::watt::get_Ein_range --------------
//  Gets the range of nontrivial incident energy bins; computes first_Ein and last_Ein
// returns true if the threshold is too high for the energy bins
bool Wtt::watt::get_Ein_range( const Ddvec::dd_vector& sigma,
			       const Ddvec::dd_vector& mult,
    const Ddvec::dd_vector& weight,
    const Lgdata::Flux_List& e_flux, const Egp::Energy_groups& Ein_groups )
{
  double E_last;

  Wtt::Watt_param initial_param;
  bool done = initial_param.get_Ein_range( sigma, mult, weight, e_flux,
                                         Ein_groups, &E_first, &E_last );
  if( done ) return true;

  // check the range of incident energies for the probability data
  Wtt::watt::const_iterator Theta_ptr = begin( );
  double E_data = Theta_ptr->x;
  if( E_data > E_first )
  {
    E_first = E_data;
  }
  Ddvec::dd_vector::const_iterator b_ptr = b_data.begin( );
  E_data = b_ptr->x;
  if( E_data > E_first )
  {
    E_first = E_data;
  }
  first_Ein = Ein_groups.first_bin_ID( E_first );

  Theta_ptr = end( );
  --Theta_ptr;
  E_data = Theta_ptr->x;
  if( E_data < E_last )
  {
    E_last = E_data;
  }
  b_ptr = b_data.end( );
  --b_ptr;
  E_data = b_ptr->x;
  if( E_data < E_last )
  {
    E_last = E_data;
  }
  last_Ein = Ein_groups.last_bin_ID( E_last );

  return false;
}
// -----------  Wtt::watt::set_Ein_range ------------------
// Sets the range of incident energies for this intergration
void Wtt::watt::set_Ein_range( int Ein_bin, Efunc::E_function_param *ein_param )
{
  // the parameters are really Wtt::Watt_param *
  Wtt::Watt_param *Ein_param = static_cast<Watt_param *>( ein_param );

  set_Ein_range_default( Ein_bin, Ein_param );
  double this_E = Ein_param->this_b->x;
  if( this_E > Ein_param->data_E_0 ) Ein_param->data_E_0 = this_E;
  this_E = Ein_param->next_b->x;
  if( this_E < Ein_param->data_E_1 ) Ein_param->data_E_1 = this_E;

  if( Ein_param->data_E_1 < Ein_param->data_E_0 )
  {
    Msg::FatalError( "Wtt::watt::set_Ein_range",
		     "check the Watt incident energies" );
  }
}
// ----------- Wtt::watt::get_T --------------
// Calculates the transfer matrix for this particle.
// sigma is the cross section.
void Wtt::watt::get_T( const Ddvec::dd_vector& sigma, const Ddvec::dd_vector& multiple,
  const Ddvec::dd_vector& weight, Trf::T_matrix& transfer )
{
  if( interp_type != Terp::LINLIN )
  {
    Msg::FatalError( "evaporation::get_T",
		     "interp_type for Theta not implemented" );
  }
  transfer.getBinCrossSection( sigma );
  transfer.threshold = sigma.begin( )->x;

  // synchronize the starting energies
  bool done = get_Ein_range( sigma, multiple, weight, transfer.e_flux,
    transfer.in_groups );
  if( done )
  {
    transfer.zero_transfer( );
  }

  long int quad_count = 0;  // number of quadratures
  long int Ein_F_count= 0;  // number of calls to Energy_function_F::Ein_F

  // do the integrals incident bin by incident bin
#pragma omp parallel for schedule( dynamic, 1 ) default( none )	\
  shared( sigma, multiple, weight, transfer ) \
  reduction( +: quad_count ) reduction( +: Ein_F_count )
  for( int Ein_bin = first_Ein; Ein_bin < last_Ein; ++Ein_bin )
  {
    Wtt::Watt_param Ein_param;
    Ein_param.Eout_quad_rule = transfer.Eout_quad_rule;
    // set up the data range for this bin
    Ein_param.setup_bin( Ein_bin, sigma, multiple, weight, transfer.e_flux,
                         transfer.in_groups );
    setup_data( transfer.out_groups, &Ein_param );
    // work on this bin
    for( ; ; )
    {
      set_Ein_range( Ein_bin, &Ein_param );   // get the incident energy interval
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
  std::cout << "2d quadratures: " << quad_count << std::endl;
  std::cout << "Energy_function_F::Ein_F calls: " << Ein_F_count << std::endl;
  std::cout << "average Energy_function_F::Ein_F calls: " << 1.0*Ein_F_count/quad_count << std::endl;
}
// -----------  Wtt::watt::next_ladder ------------------
// Go to the next pair of incident energies.  Returns "true" when finished.
bool Wtt::watt::next_ladder( double E_in, Efunc::E_function_param *ein_param )
{
  // the parameters are really Wtt::Watt_param *
  Wtt::Watt_param *Ein_param = static_cast<Wtt::Watt_param *>( ein_param );

  bool done = next_ladder_default( E_in, Ein_param );
  static double etol = Global.Value( "tight_tol" );
  double E_tol = E_in * etol;
  if( !done )
  {
    if( E_in + E_tol >= Ein_param->this_b->x )
    {
      while( E_in + E_tol >= Ein_param->next_b->x )
      {
        // get the next E_in Watt data
        Ein_param->this_b = Ein_param->next_b;
        ++Ein_param->next_b;
        if( Ein_param->next_b == Ein_param->b_end )
        {
          return true;
        }
      }
    }
  }
  return done;
}

// **************** Function to integrate *********************
// --------------------  Watt_F::check_Watt_F ------------------
// Function for the 1-d quadrature over outgoing energy
bool Watt_F::check_Watt_F( double E_out, Qparam::QuadParamBase *E_out_param,
			   Coef::coef_vector *value )
{   
  // the parameters are really Wtt::check_Watt_params
  Wtt::check_Watt_params *Eout_params =
    static_cast< Wtt::check_Watt_params* >( E_out_param );
  Eout_params->func_count += 1;

  value->weight_1[ 0 ] = exp( -E_out/Eout_params->a ) *
    sinh( sqrt( Eout_params->b * E_out ) );

  return true;
}
