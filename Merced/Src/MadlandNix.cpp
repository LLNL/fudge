/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: MadlandNix.cpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */
// Implement classes used for Madland-Nix energy probability density

#include <cmath>
#include <cfloat>
#ifdef _OPENMP
 #include <omp.h>
#endif

#include "MadlandNix.hpp"
#include "global_params.hpp"
#include "messaging.hpp"
#include "math_util.hpp"

#include "protos.h"

// ****************** class MadlandNix_param **********************
// ---------------- MadlandNix_param::set_Ein --------------------
void MadlandNix_param::set_Ein( double E_in )
{
  // Don't use set_Ein_default, because we need scale_g_light and scale_g_heavy to compute the norm
  Theta = this_Theta->linlin_interp( E_in, *next_Theta );
  multiplicity = this_mult->linlin_interp( E_in, *next_mult );
  if( ( Theta <= 0.0 ) || ( multiplicity < 0.0 ) )
  {
    FatalError( "MadlandNix_param::set_Ein", "got a negative parameter" );
  }

  E_1 = ( use_Eout_max ) ? Eout_max : E_max;

  alpha = sqrt( Theta );
  set_scales( );
  norm = get_norm( );
  Eout_0 = Eout_min;
  Eout_1 = E_1;
}
// ---------------- MadlandNix_param::get_integrals --------------------
// Gets the integrals over outgoing energy
void MadlandNix_param::get_integrals( double Eout_0, double Eout_1, coef_vector &value )
{
  double numerator;

  if( ( value.conserve == NUMBER ) || ( value.conserve == BOTH ) )
  {
    // integrate probability density
    conserve_flag = NUMBER;
    numerator = integral_prob( Eout_0, Eout_1 );
    value.weight_1[ 0 ] = numerator/norm;
  }
  if( ( value.conserve == ENERGY ) || ( value.conserve == BOTH ) )
  {
    // integrate E*probability density
    conserve_flag = ENERGY;
    numerator = integral_Eprob( Eout_0, Eout_1 );
    value.weight_E[ 0 ] = numerator/norm;
  }
}
// ---------------- MadlandNix_param::tol_get_integrals --------------------
// Gets the integrals over outgoing energy and returns the noise in the calculation
double MadlandNix_param::tol_get_integrals( double Eout_0, double Eout_1,
  coef_vector &value )
{
  double int_g_light;
  double noise = get_noise( &int_g_light );
  // for the heavier fragment
  beta = beta_heavy;
  double int_g_heavy = integral_g( );
  value.weight_1[0] = 0.5*( scale_g_light*int_g_light + scale_g_heavy*int_g_heavy );
  return noise;
}
// ---------------- MadlandNix_param::set_scales --------------------
//! Sets the scale factors for the integrals of probability and energy*probability
void MadlandNix_param::set_scales( )
{
  // the scale factors
  prob_integral_scale = 1.0;  // not used in this model
  Eprob_integral_scale = 1.0;  // not used in this model
  scale_g_light = 1.0/( 3*alpha*beta_light );
  scale_g_heavy = 1.0/( 3*alpha*beta_heavy );
}
// ---------------- MadlandNix_param::get_norm --------------------
// Interpolate the norm
double MadlandNix_param::get_norm( )
{
  return integral_prob( 0.0, E_max );
}
// ---------------- MadlandNix_param::integral_prob --------------------
// Integral of the probability density
double MadlandNix_param::integral_prob( double E_0, double E_1 )
{
  if( E_0 >= E_1 )
  {
    Info( "MadlandNix_param::integral_prob", "energies out of order" );
    return 0.0;
  }
  sqrt_E0 = sqrt( E_0 );
  sqrt_E1 = sqrt( E_1 );
  // for the lighter fragment
  beta = beta_light;
  double int_g_light = integral_g( );

  // for the heavier fragment
  beta = beta_heavy;
  double int_g_heavy = integral_g( );
  return 0.5*( scale_g_light*int_g_light + scale_g_heavy*int_g_heavy );
}
// ---------------- MadlandNix_param::integral_g --------------------
// Integral of the g function
double MadlandNix_param::integral_g( )
{
  // the parameters for this fragment
  A = ( sqrt_E0 + beta )*( sqrt_E0 + beta )/Theta;
  B = ( sqrt_E1 + beta )*( sqrt_E1 + beta )/Theta;
  Aprime = ( sqrt_E0 - beta )*( sqrt_E0 - beta )/Theta;
  Bprime = ( sqrt_E1 - beta )*( sqrt_E1 - beta )/Theta;

  double int_g = integral_u2( );
  // which incomplete gamma function to use
  use_tail = false;
  if( ( sqrt_E0 > beta ) && ( Aprime > 1.0 ) ) use_tail = true;
  if( ( sqrt_E1 < beta ) && ( Bprime > 1.0 ) ) use_tail = true;

  if( sqrt_E1 > beta )
  {
    int_g -= integral_u1_above( Bprime, sqrt_E1 );
  }
  else if( sqrt_E1 < beta )
  {
    int_g += integral_u1_below( Bprime, sqrt_E1 );
  }
  if( sqrt_E0 > beta )
  {
    int_g += integral_u1_above( Aprime, sqrt_E0 );
  }
  else if( sqrt_E0 < beta )
  {
    int_g -= integral_u1_below( Aprime, sqrt_E0 );
  }
  if( int_g < 0.0 ) int_g = 0.0;
  return int_g;
}
// ---------------- MadlandNix_param::integral_u2 --------------------
// Contribution of u_2 to integral_prob
double MadlandNix_param::integral_u2( )
{
  double pow15_A = A*sqrt( A );
  double pow15_B = B*sqrt( B );
  double E1_term = ( 0.4*Theta*B*pow15_B - 0.5*alpha*beta*B*B )*expn( 1, B ) -
    ( 0.4*Theta*A*pow15_A - 0.5*alpha*beta*A*A )*expn( 1, A );
  double gamma15_term;  // gamma(1.5, B) - gamma(1.5, A)
  double gamma20_term;  // gamma(2.0, B) - gamma(2.0, A)
  double gamma25_term;  // gamma(2.5, B) - gamma(2.5, A)
  if( A <= 1.0 )
  {
    // use the incomplete gamma function \int_0^A
    double gamma25 = tgamma( 2.5 );
    double igam25_A = gamma25*igam( 2.5, A );
    double igam25_B = gamma25*igam( 2.5, B );
    double igam15_A = math_F::gamma_down( 2.5, A, igam25_A );
    double igam15_B = math_F::gamma_down( 2.5, B, igam25_B );
    gamma15_term = ( sqrt_E1 - beta )*( sqrt_E1 + beta )*igam15_B -
      ( sqrt_E0 - beta )*( sqrt_E0 + beta )*igam15_A;
    gamma20_term = igam( 2.0, B ) - igam( 2.0, A );
    //   gamma25_term = igam( 2.5, B ) - igam( 2.5, A );
    gamma25_term = igam25_B - igam25_A;
  }
  else
  {
    // use the incomplete gamma function \int_A^\infty
    double gamma15 = tgamma( 1.5 );
    double igamc15_A = gamma15*igamc( 1.5, A );
    double igamc15_B = gamma15*igamc( 1.5, B );
    double igamc25_A = math_F::Gamma_up( 1.5, A, igamc15_A );
    double igamc25_B = math_F::Gamma_up( 1.5, B, igamc15_B );
    gamma15_term = gamma15*( sqrt_E1 - sqrt_E0 )*( sqrt_E1 + sqrt_E0 ) +
      ( sqrt_E0 - beta )*( sqrt_E0 + beta )*igamc15_A -
      ( sqrt_E1 - beta )*( sqrt_E1 + beta )*igamc15_B;
    gamma20_term = igamc( 2.0, A ) - igamc( 2.0, B );
    //    gamma25_term = igamc( 2.5, A ) - igamc( 2.5, B );
    gamma25_term = igamc25_A - igamc25_B;
  }
  double svar = E1_term + gamma15_term + 1.5*alpha*beta*gamma20_term -
    0.6*Theta*gamma25_term;
  return svar;
}
// ---------------- MadlandNix_param::integral_u1_above --------------------
// Contribution of u_1 to integral_prob for E > Ef
double MadlandNix_param::integral_u1_above( double Bprim, double sqrt_En )
{
  double pow15_B = Bprim*sqrt( Bprim );
  // This is the integral of u_1 from Ef to E
  double E1_term = ( 0.4*Theta*pow15_B*Bprim + 0.5*alpha*beta*Bprim*Bprim )*expn( 1, Bprim );
  double gamma15_term;
  double gamma20_term;
  double gamma25_term;
  if( use_tail )
  {
    // use the incomplete gamma function \int_A^\infty
    double gamma15 = tgamma( 1.5 );
    double igamc15_B = gamma15*igamc( 1.5, Bprim );
    gamma15_term = gamma15*sqrt_En*sqrt_En -
      ( sqrt_En - beta )*( sqrt_En + beta )*igamc15_B;
    gamma20_term = 1.5*alpha*beta*igamc( 2.0, Bprim );
    //    gamma25_term = 0.6*Theta*tgamma( 2.5 )*igamc( 2.5, Bprim );
    double igamc25_B = math_F::Gamma_up( 1.5, Bprim, igamc15_B );
    gamma25_term = 0.6*Theta*igamc25_B;
  }
  else
  {
    // use the incomplete gamma function \int_0^A
    double gamma25 = tgamma( 2.5 );
    double igam25_B = gamma25*igam( 2.5, Bprim );
    double igam15_B = math_F::gamma_down( 2.5, Bprim, igam25_B );
    gamma15_term = ( sqrt_En - beta )*( sqrt_En + beta )*igam15_B;
    gamma20_term = -1.5*alpha*beta*igam( 2.0, Bprim );
    //    gamma25_term = -0.6*Theta*tgamma( 2.5 )*igam( 2.5, Bprim );
    gamma25_term = -0.6*Theta*igam25_B;
  }
  double svar = E1_term + gamma15_term + gamma20_term + gamma25_term;
  return svar;
}
// ---------------- MadlandNix_param::integral_u1_below --------------------
// Contribution of u_1 to integral_prob for E < Ef
double MadlandNix_param::integral_u1_below( double Aprim, double sqrt_En )
{
  double pow15_A = Aprim*sqrt( Aprim );
  // This is the integral of u_1 from E to Ef
  double E1_term = ( -0.4*Theta*pow15_A*Aprim + 0.5*alpha*beta*Aprim*Aprim )*expn( 1, Aprim );
  double gamma15_term;
  double gamma20_term;
  double gamma25_term;
  if( use_tail )
  {
    // use the incomplete gamma function \int_A^\infty
    double gamma15 = tgamma( 1.5 );
    double igamc15_A = gamma15*igamc( 1.5, Aprim );
    gamma15_term = -gamma15*sqrt_En*sqrt_En -
      ( beta + sqrt_En )*( beta - sqrt_En )*igamc15_A;
    gamma20_term = 1.5*alpha*beta*igamc( 2.0, Aprim );
    //    gamma25_term = -0.6*Theta*tgamma( 2.5 )*igamc( 2.5, Aprim );
    double igamc25_A = math_F::Gamma_up( 1.5, Aprim, igamc15_A );
    gamma25_term = -0.6*Theta*igamc25_A;
  }
  else
  {
    // use the incomplete gamma function \int_0^A
    double gamma25 = tgamma( 2.5 );
    double igam25_A = gamma25*igam( 2.5, Aprim );
    double igam15_A = math_F::gamma_down( 2.5, Aprim, igam25_A );
    gamma15_term = ( beta + sqrt_En )*( beta - sqrt_En )*igam15_A;
    gamma20_term = -1.5*alpha*beta*igam( 2.0, Aprim );
    //    gamma25_term = 0.6*Theta*tgamma( 2.5 )*igam( 2.5, Aprim );
    gamma25_term = 0.6*Theta*( 1.5*igam15_A - pow15_A*exp( -Aprim) );
 }
  double svar = E1_term + gamma15_term + gamma20_term + gamma25_term;
  return svar;
}

// ---------------- MadlandNix_param::integral_Eprob --------------------
// Integral of energy times the probability density
double MadlandNix_param::integral_Eprob( double E_0, double E_1 )
{
  if( E_0 >= E_1 )
  {
    Info( "MadlandNix_param::integral_Eprob", "energies out of order" );
    return 0.0;
  }
  return 0.0;
}
// ---------------- MadlandNix_param::integral_Eu2 --------------------
// Contribution of u_2 to integral_Eprob
double MadlandNix_param::integral_Eu2( )
{
  return 0.0;
}
// ---------------- MadlandNix_param::integral_Eu1_above --------------------
// Contribution of u_1 to integral_Eprob for E > Ef
double MadlandNix_param::integral_Eu1_above( double Bprim, double sqrt_En )
{
  return 0.0;
}
// ---------------- MadlandNix_param::integral_Eu1_below --------------------
// Contribution of u_1 to integral_Eprob for E < Ef
double MadlandNix_param::integral_Eu1_below( double Aprim, double sqrt_En )
{
  return 0.0;
}
// ---------------- MadlandNix_param::set_tol --------------------
// Sets the tolerance for the quadrature over incident energy
double MadlandNix_param::set_tol( double left_E, double right_E )
{
  // This routine is used when transfer.interpolate_Eout_integrals == false
  // check the numerics for a "typical" incident energy
  // The incident energy is left_E < E_in < right_E
  set_Ein( 0.5*( left_E + right_E ) );
  double int_g_light;
  double noise = get_noise( &int_g_light );
  return noise;
}
// ---------------- MadlandNix_param::get_noise --------------------
// Calculates one integral and returns the noise in the calculation
double MadlandNix_param::get_noise( double *int_g_light )
{
  // The outgoing energy is Eout_min < E_out < E_1
  // the parameters for the lighter fragment
  sqrt_E0 = sqrt( param_base::Eout_min );
  sqrt_E1 = sqrt( E_1 );
  beta = beta_light;
  A = ( sqrt_E0 + beta )*( sqrt_E0 + beta )/Theta;
  B = ( sqrt_E1 + beta )*( sqrt_E1 + beta )/Theta;
  Aprime = ( sqrt_E0 - beta )*( sqrt_E0 - beta )/Theta;
  Bprime = ( sqrt_E1 - beta )*( sqrt_E1 - beta )/Theta;

  double int_g = integral_u2( );
  double biggest = abs( int_g );
  // which incomplete gamma function to use
  use_tail = false;
  if( ( sqrt_E0 > beta ) && ( Aprime > 1.0 ) ) use_tail = true;
  if( ( sqrt_E1 < beta ) && ( Bprime > 1.0 ) ) use_tail = true;

  if( sqrt_E1 > beta )
  {
    int_g -= integral_u1_above( Bprime, sqrt_E1 );
    if( abs( int_g ) > biggest ) biggest = abs( int_g );
  }
  else if( sqrt_E1 < beta )
  {
    int_g += integral_u1_below( Bprime, sqrt_E1 );
    if( abs( int_g ) > biggest ) biggest = abs( int_g );
  }
  if( sqrt_E0 > beta )
  {
    int_g += integral_u1_above( Aprime, sqrt_E0 );
    if( abs( int_g ) > biggest ) biggest = abs( int_g );
  }
  else if( sqrt_E0 < beta )
  {
    int_g -= integral_u1_below( Aprime, sqrt_E0 );
    if( abs( int_g ) > biggest ) biggest = abs( int_g );
  }
  if( int_g < 0.0 ) int_g = 0.0;
  *int_g_light = int_g;
  double delta = abs( int_g )/biggest;
  double noise = 2*DBL_EPSILON;
  static double quad_tol = Global.Value( "quad_tol" );
  if( delta <= noise ) return 1.0;  // really noisy integral
  else if( delta >= noise/quad_tol ) return quad_tol;  // normal
  else return noise/delta;
}

// ****************** class MadlandNix **********************
// ----------- MadlandNix::get_Ein_range --------------
//  Gets the range of nontrivial incident energy bins; computes first_Ein and last_Ein
// returns true if the threshold is too high for the energy bins
bool MadlandNix::get_Ein_range( const dd_vector& sigma, const dd_vector& mult,
    const dd_vector& weight,
    const Flux_List& e_flux, const Energy_groups& Ein_groups )
{
  double E_last;

  MadlandNix_param initial_param;
  bool done = initial_param.get_Ein_range( sigma, mult, weight, e_flux,
                                         Ein_groups, &E_first, &E_last );
  if( done ) return true;

  // check the range of incident energies for the probability data
  MadlandNix::const_iterator TM_ptr = begin( );
  double E_data = TM_ptr->x;
  if( E_data > E_first )
  {
    E_first = E_data;
  }
  first_Ein = Ein_groups.first_bin_ID( E_first );

  TM_ptr = end( );
  --TM_ptr;
  E_data = TM_ptr->x;
  if( E_data < E_last )
  {
    E_last = E_data;
  }
  last_Ein = Ein_groups.last_bin_ID( E_last );

  return false;
}
// ---------------- MadlandNix::setup_data --------------------
// Initializes the quadrature parameters
void MadlandNix::setup_data( const Energy_groups& Eout_groups,
			  E_function_param *ein_param )
{
  // the parameters are really MadlandNix_param *
  MadlandNix_param *Ein_param = static_cast<MadlandNix_param *>( ein_param );

  setup_data_default( Eout_groups, Ein_param );

  if( maxEout <= 0.0 ) // parameter not set
  {
    maxEout = 0.5*Ein_param->top_E_out;
  }
  // the range of integration
  Ein_param->E_max = maxEout;
  if( Ein_param->E_max > Ein_param->top_E_out )
  {
    Ein_param->E_max = Ein_param->top_E_out;
  }

  // kludge the model
  Ein_param->U = -Ein_param->E_max;
  Ein_param->beta_light = sqrt( Efl );
  Ein_param->beta_heavy = sqrt( Efh );
}
// ----------- MadlandNix::get_T --------------
// Calculates the transfer matrix for this particle.
// sigma is the cross section.
void MadlandNix::get_T( const dd_vector& sigma, const dd_vector& multiple,
  const dd_vector& weight, T_matrix& transfer )
{
  if( interp_type != LINLIN )
  {
    FatalError( "evaporation::get_T", "interp_type for Theta not implemented" );
  }
  transfer.getBinCrossSection( sigma );

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
    MadlandNix_param Ein_param;
    Ein_param.Eout_quad_method = transfer.Eout_quad_method;
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
      bool done = next_ladder( Ein_param.data_E_1, &Ein_param );   // go to the next interval
      if( done )
      {
        break;
      }
    }
    quad_count += Ein_param.quad_count;
    Ein_F_count += Ein_param.Ein_F_count;
  } // end of parallel loop

  // print the counts of function evaluations
  cout << "2d quadratures: " << quad_count << endl;
  cout << "Energy_function_F::Ein_F calls: " << Ein_F_count << endl;
  cout << "average Energy_function_F::Ein_F calls: " << 1.0*Ein_F_count/quad_count << endl;
}

