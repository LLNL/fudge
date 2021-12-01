/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: Maxwell.cpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */
// Implement classes used for Maxwell energy probability density

#include <cmath>
#ifdef _OPENMP
 #include <omp.h>
#endif

#include "Maxwell.hpp"
#include "adapt_quad.hpp"
#include "global_params.hpp"
#include "messaging.hpp"

#include "protos.h"

// ****************** class Max::Maxwell_param **********************
// ---------------- Max::Maxwell_param::get_integrals --------------------
// Gets the integrals over outgoing energy
bool Max::Maxwell_param::get_integrals( double Eout_0, double Eout_1,
					Coef::coef_vector &value )
{
  bool is_OK = true;
  
  if( Eout_0 == 0.0 )
  {
    if( ( value.conserve == Coef::NUMBER ) || ( value.conserve == Coef::BOTH ) )
    {
      value.weight_1[ 0 ] = integral_prob( Eout_1 );
    }
    if( ( value.conserve == Coef::ENERGY ) || ( value.conserve == Coef::BOTH ) )
    {
      value.weight_E[ 0 ] = integral_Eprob( Eout_1 );
    }
  }
  else
  {
    static double tol = Global.Value( "quad_tol" );
    Evap::evap_params_1d params_1d;
    params_1d.Theta = Theta;
    Qparam::QuadParamBase *void_params_1d =
      static_cast< Qparam::QuadParamBase* >( &params_1d );

    is_OK = quad_F::integrate( Maxwell_F::Eout_F, Eout_quad_rule, Eout_0, Eout_1,
		       void_params_1d, tol, &value );
    Eout_F_count += params_1d.func_count;
  }
  value *= 1.0/norm;

  return is_OK;
}
// ---------------- Max::Maxwell_param::set_scales --------------------
//! Sets the scale factors for the integrals of probability and energy*probability
void Max::Maxwell_param::set_scales( )
{
  // the scale factors used with the incomplete gamma functions
  prob_integral_scale = tgamma( 1.5 )*Theta*sqrt( Theta );
  Eprob_integral_scale = 1.5*prob_integral_scale*Theta;  // includes Gamma( 2.5 )
}
// ---------------- Max::Maxwell_param::get_norm --------------------
// Integrate from 0 to E_max to get the norm
double Max::Maxwell_param::get_norm( )
{
  return integral_prob( E_max );
}
// ---------------- Max::Maxwell_param::integral_prob --------------------
// Integral of the probability density from 0 to E_1
double Max::Maxwell_param::integral_prob( double E_1 )
{
  // This is Theta times the integral of y^{1/2}*exp( -y ) with y = E_out/Theta
  double y_1 = E_1/Theta;
  return prob_integral_scale*Proto::igam( 1.5, y_1 );
}
// ---------------- Max::Maxwell_param::integral_Eprob --------------------
// Integral of energy times the probability density from 0 to E_1
double Max::Maxwell_param::integral_Eprob( double E_1 )
{
  // This is Theta times the integral of y^{3/2}*exp( -y ) with y = E_out/Theta
  double y_1 = E_1/Theta;
  return Eprob_integral_scale*Proto::igam( 2.5, y_1 );
}

// ****************** class Max::maxwell **********************
// ----------- Max::maxwell::get_T --------------
// Calculates the transfer matrix for this particle.
// sigma is the cross section.
void Max::maxwell::get_T( const dd_vector& sigma, const dd_vector& multiple,
  const dd_vector& weight, Trf::T_matrix& transfer )
{
  if( interp_type != Terp::LINLIN )
  {
    Msg::FatalError( "Max::maxwell::get_T",
		     "interp_type for Theta not implemented" );
  }
  transfer.getBinCrossSection( sigma );
  threshold =sigma.begin( )->x;
  
  // synchronize the starting energies
  bool done = get_Ein_range( sigma, multiple, weight, transfer.e_flux,
    transfer.in_groups );
  if( done )
  {
    transfer.zero_transfer( );
  }

  long int quad_count = 0;  // number of quadratures
  long int Ein_F_count= 0;  // number of calls to Energy_function_F::Ein_F
  long int Eout_F_count= 0;  // number of calls to evaporation_F::Eout_F

  int num_processors;
#ifdef _OPENMP
  num_processors = omp_get_num_procs( );
#else
  num_processors = 1;
#endif
  std::cout << "number of processors: " << num_processors << std::endl;

  // do the integrals incident bin by incident bin
#pragma omp parallel for schedule( dynamic, 1 ) default( none )	\
  shared( sigma, multiple, weight, transfer ) \
  reduction( +: quad_count ) reduction( +: Ein_F_count ) reduction( +: Eout_F_count )
  // now do the integrals incident bin by incident bin
  for( int Ein_bin = first_Ein; Ein_bin < last_Ein; ++Ein_bin )
  {
    Max::Maxwell_param Ein_param;
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
      bool done = next_ladder( Ein_param.data_E_1, &Ein_param );   // go to the next interval
      if( done )
      {
        break;
      }
    }
    quad_count += Ein_param.quad_count;
    Ein_F_count += Ein_param.Ein_F_count;
    Eout_F_count += Ein_param.Eout_F_count;
  }
  // print the counts of function evaluations
  std::cout << "2d quadratures: " << quad_count << std::endl;
  std::cout << "Energy_function_F::Ein_F calls: " << Ein_F_count << std::endl;
  std::cout << "evaporation_F::Eout_F calls: " << Eout_F_count << std::endl;
  std::cout << "average Energy_function_F::Ein_F calls: " << 1.0*Ein_F_count/quad_count << std::endl;
  std::cout << "average evaporation_F::Eout_F calls: " << 1.0*Eout_F_count/Ein_F_count << std::endl;
}

// **************** Function to integrate *********************
// --------------------  Maxwell_F::Eout_F ------------------
// Function for the 1-d quadrature over outgoing energy
bool Maxwell_F::Eout_F( double E_out, Qparam::QuadParamBase *E_out_param,
			Coef::coef_vector *value )
{   
  // the parameters are really Evap::evap_params_1d
  Evap::evap_params_1d *Eout_params = static_cast< Evap::evap_params_1d* >( E_out_param );
  Eout_params->func_count += 1;
  double Prob = sqrt( E_out )*exp( - E_out/Eout_params->Theta );
  if( ( value->conserve == Coef::NUMBER ) || ( value->conserve == Coef::BOTH ) )
  {
    value->weight_1[ 0 ] = Prob;
  }
  if( ( value->conserve == Coef::ENERGY ) || ( value->conserve == Coef::BOTH ) )
  {
    value->weight_E[ 0 ] = E_out*Prob;
  }

  return true;
}
