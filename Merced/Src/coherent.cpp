/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2008-04-16 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: coherent.cpp 1 2008-07-01 03:06:56Z hedstrom $
 * initial code
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// implementation of the classes used for coherent scattering

#include <cmath>
#ifdef _OPENMP
 #include <omp.h>
#endif

#include "coherent.hpp"
#include "adapt_quad.hpp"
#include "math_util.hpp"
#include "messaging.hpp"
#include "global_params.hpp"

// **************** class Coh::anomalous_ptrs ******************
// ----------------  Coh::anomalous_ptrs::get_Re_anomalous ------------------------
// Gets the real anomalous factor at this energy
double Coh::anomalous_ptrs::get_Re_anomalous( double E_in, bool *interp_OK )
{
  double answer;
  if( next_Re_anomalous == end_Re_anomalous )
  {
    answer = 0.0;
  }
  else
  {
    answer = Re_anomalous_ptr->linlin_interp( E_in, *next_Re_anomalous, interp_OK );
  }
  return answer;
}
// ----------------  Coh::anomalous_ptrs::get_Im_anomalous ------------------------
// Gets the imaginary anomalous factor at this energy
double Coh::anomalous_ptrs::get_Im_anomalous( double E_in, bool *interp_OK )
{
  double answer;
  if( next_Im_anomalous == end_Im_anomalous )
  {
    answer = 0.0;
  }
  else
  {
    answer = Im_anomalous_ptr->linlin_interp( E_in, *next_Im_anomalous, interp_OK );
  }
  return answer;
}
// ----------------  Coh::anomalous_ptrs::init_Re_anomalous ------------------------
// Initializes the pointers to real anomalous factor data
void Coh::anomalous_ptrs::init_Re_anomalous( Ddvec::dd_vector &Re_anomalous )
{
  start_Re_anomalous = Re_anomalous.begin( );
  first_ladder_Re_anomalous = Re_anomalous.begin( );
  next_Re_anomalous = first_ladder_Re_anomalous;
  Re_anomalous_ptr = first_ladder_Re_anomalous;
  ++next_Re_anomalous;
  end_Re_anomalous = Re_anomalous.end( );
}
// ----------------  Coh::anomalous_ptrs::init_Im_anomalous ------------------------
// Initializes the pointers to imaginary anomalous factor data
void Coh::anomalous_ptrs::init_Im_anomalous( Ddvec::dd_vector &Im_anomalous )
{
  start_Im_anomalous = Im_anomalous.begin( );
  first_ladder_Im_anomalous = Im_anomalous.begin( );
  Im_anomalous_ptr = first_ladder_Im_anomalous;
  next_Im_anomalous = first_ladder_Im_anomalous;
  ++next_Im_anomalous;
  end_Im_anomalous = Im_anomalous.end( );
}
// ----------------  Coh::anomalous_ptrs::init_Re_down ------------------------
// Initializes the pointers to real anomalous factor data, count down
void Coh::anomalous_ptrs::init_Re_down( Ddvec::dd_vector &Re_anomalous )
{
  start_Re_anomalous = Re_anomalous.begin( );
  end_Re_anomalous = Re_anomalous.end( );
  next_Re_anomalous = end_Re_anomalous;
  --next_Re_anomalous;
  Re_anomalous_ptr = next_Re_anomalous;
  --Re_anomalous_ptr;
}
// ----------------  Coh::anomalous_ptrs::init_Im_down ------------------------
// Initializes the pointers to imaginary anomalous factor data, count down
void Coh::anomalous_ptrs::init_Im_down( Ddvec::dd_vector &Im_anomalous )
{
  start_Im_anomalous = Im_anomalous.begin( );
  end_Im_anomalous = Im_anomalous.end( );
  next_Im_anomalous = end_Im_anomalous;
  --next_Im_anomalous;
  Im_anomalous_ptr = next_Im_anomalous;
  --Im_anomalous_ptr;
}
// ----------------  Coh::anomalous_ptrs::set_range ------------------------
// Sets the range of data for this quadrature
void Coh::anomalous_ptrs::set_range( double E_0, double E_1 )
{
  // first do real anomalous
  Ddvec::dd_vector::const_iterator data_ptr = first_ladder_Re_anomalous;
  ++data_ptr;
  while( ( data_ptr != end_Re_anomalous ) &&
	 ( data_ptr->x <= E_0 ) )
  {
    ++data_ptr;
  }
  next_Re_anomalous = data_ptr;
  first_ladder_Re_anomalous = data_ptr;
  --first_ladder_Re_anomalous;
  Re_anomalous_ptr = first_ladder_Re_anomalous;

  while( ( data_ptr != end_Re_anomalous ) &&
	 ( data_ptr->x < E_1 ) )
  {
    ++data_ptr;
  }
  last_ladder_Re_anomalous = data_ptr;

  // now imaginary anomalous
  data_ptr = first_ladder_Im_anomalous;
  ++data_ptr;
  while( ( data_ptr != end_Im_anomalous ) &&
	 ( data_ptr->x <= E_0 ) )
  {
    ++data_ptr;
  }
  next_Im_anomalous = data_ptr;
  first_ladder_Im_anomalous = data_ptr;
  --first_ladder_Im_anomalous;
  Im_anomalous_ptr = first_ladder_Im_anomalous;

  while( ( data_ptr != end_Im_anomalous ) &&
	 ( data_ptr->x < E_1 ) )
  {
    ++data_ptr;
  }
  last_ladder_Im_anomalous = data_ptr;

}
// ----------------  Coh::anomalous_ptrs::get_next_Re_down ------------------------
// Increments the pointers to real anomalous factor data, count down
void Coh::anomalous_ptrs::get_next_Re_down( double E_in )
{
  while( ( Re_anomalous_ptr != start_Re_anomalous ) &&
	 ( E_in < Re_anomalous_ptr->x ) )
  {
    next_Re_anomalous = Re_anomalous_ptr;
    --Re_anomalous_ptr;
  }
}
// ----------------  Coh::anomalous_ptrs::get_next_Im_down ------------------------
// Increments the pointers to imaginary anomalous factor data, count down
void Coh::anomalous_ptrs::get_next_Im_down( double E_in )
{
  while( ( Im_anomalous_ptr != start_Im_anomalous ) &&
	 ( E_in < Im_anomalous_ptr->x ) )
  {
    next_Im_anomalous = Im_anomalous_ptr;
    --Im_anomalous_ptr;
  }
}
// ----------------  Coh::anomalous_ptrs::next_range ------------------------
// Increments the pointers to anomalous factor data, count up
bool Coh::anomalous_ptrs::next_range( double right_E, double Ein_1 )
{
  bool done = false;

  // increment real anomalous
  if( next_Re_anomalous->x <= right_E )
  {
    Re_anomalous_ptr = next_Re_anomalous;
    ++next_Re_anomalous;
    if( ( Re_anomalous_ptr == last_ladder_Re_anomalous ) ||
	( Re_anomalous_ptr->x >= Ein_1 ) )
    {
      return true;
    }
  }

  // increment imaginary anomalous
  if( next_Im_anomalous->x <= right_E )
  {
    Im_anomalous_ptr = next_Im_anomalous;
    ++next_Im_anomalous;
    if( ( Im_anomalous_ptr == last_ladder_Im_anomalous ) ||
	( Im_anomalous_ptr->x >= Ein_1 ) )
    {
      return true;
    }
  }

  return done;
}

// **************** class Coh::coherent_mu_param ******************
// ----------------  Coh::coherent_mu_param::get_sigma ------------------------
// Gets the angular probability density
double Coh::coherent_mu_param::get_sigma( double mu, bool *interp_OK )
{
  double E_in = get_E_in();
  double x = E_in * sqrt( 0.5 * ( 1.0 - mu ) );
  bool FF_OK;
  double FF = value( x, &FF_OK );   // the coherent scattering factor

  bool Re_OK;
  bool Im_OK;
  
  // get the coherent scattering cross section
  static double Thompson = Global.Value( "Thompson" );  // Thompson cross section
  double Re_anomalous = anomalous.get_Re_anomalous( E_in, &Re_OK );
  double Im_anomalous = anomalous.get_Im_anomalous( E_in, &Im_OK );
  double Re_sq = ( FF + Re_anomalous ) * ( FF + Re_anomalous );
  double Im_sq = Im_anomalous * Im_anomalous;
  double xs = 3.0 * Thompson * ( Re_sq + Im_sq ) *( 1.0 + mu*mu ) / 8;

  *interp_OK = FF_OK & Re_OK & Im_OK;
  return xs;
}

// **************** class Coh::coherent_Ein_param ******************
// ------------------ Coh::coherent_Ein_param::start_sigma --------------
// Sets up the loop over cross section and anomalous data
// Returns false if Ein_1 is below the data
bool Coh::coherent_Ein_param::start_sigma( )
{
  this_sigma = first_ladder_sigma;
  next_sigma = this_sigma;
  ++next_sigma;

  // Ein_1 may be below this_sigma
  if( this_sigma->x >= Ein_1 )
  {
    return false;
  }
  
  // Ein_0 may be past next_sigma
  while( ( this_sigma != last_ladder_sigma ) &&
         ( next_sigma->x < Ein_0 ) )
  {
    this_sigma = next_sigma;
    ++next_sigma;
  }

  // set up the pointers to anomalous data
  mu_param.anomalous.set_range( Ein_0, Ein_1 );

  return true;
}
// ------------------ Coh::coherent_Ein_param::get_range --------------
// Gets the current range of incident energies
void Coh::coherent_Ein_param::get_range( )
{
  left_E = ( this_sigma->x < Ein_0 ) ? Ein_0 : this_sigma->x;
  if( mu_param.anomalous.Re_anomalous_ptr->x > left_E )
  {
    left_E = mu_param.anomalous.Re_anomalous_ptr->x;
  }
  if( mu_param.anomalous.Im_anomalous_ptr->x > left_E )
  {
    left_E = mu_param.anomalous.Im_anomalous_ptr->x;
  }

  right_E = ( next_sigma->x > Ein_1 ) ? Ein_1 : next_sigma->x;
  if( mu_param.anomalous.next_Re_anomalous->x < right_E )
  {
    right_E = mu_param.anomalous.next_Re_anomalous->x;
  }
  if( mu_param.anomalous.next_Im_anomalous->x < right_E )
  {
    right_E = mu_param.anomalous.next_Im_anomalous->x;
  }

  if( right_E < left_E )
  {
    Msg::FatalError( "Coh::coherent_Ein_param::get_range",
		"Energies out of order" );
  }
}
// ------------------ Coh::coherent_Ein_param::next_range --------------
// Increments the data for the next range of incident energies
bool Coh::coherent_Ein_param::next_range( )
{
  // we may be done
  if( right_E >= Ein_1 )
  {
    return true;
  }

  // increment cross section data
  if( next_sigma->x <= right_E )
  {
    this_sigma = next_sigma;
    ++next_sigma;
  }
  // increment anomalous data
  bool done = mu_param.anomalous.next_range( right_E, Ein_1 );

  return done;
}

// **************** class Coh::coherent_hit_list ******************
// ------------------ Coh::coherent_hit_list::get_Eout --------------
// virtual function not used
double Coh::coherent_hit_list::get_Eout( double E_in )
{
  return 0.0;
}
// ------------------ Coh::coherent_hit_list::find_bottom_hits --------------
// virtual function not used
void Coh::coherent_hit_list::find_bottom_hits( double E_out,
					  std::vector< Box::Ein_Eta_Hit > *Ein_hits )
{
}
// ------------------ Coh::coherent_hit_list::find_top_hits --------------
// virtual function not used
void Coh::coherent_hit_list::find_top_hits( double E_out,
				       std::vector< Box::Ein_Eta_Hit > *Ein_hits )
{
}
// ------------------ Coh::coherent_hit_list::hit_box --------------
// Finds intersections of the curve x = const with the E'-mu box.
void Coh::coherent_hit_list::hit_box( double x, double E_in_left, double E_in_right )
{
  // start with a new list
  if( size( ) > 0 )
  {
    erase( begin( ), end( ) );
  }
  eta = x;
  Box::Ein_Eta_Hit new_hit;  // for insertion into the list
  if( x <= E_in_left )
  {
    new_hit.E_in = E_in_left;
    new_hit.hit_edge = Box::INSIDE;  // the curve enters on the left side
    push_back( new_hit );
    new_hit.E_in = E_in_right;  // the curve xits on the right side
    push_back( new_hit );
  }
  else if( x < E_in_right )
  {
    new_hit.E_in = E_in_left;
    new_hit.hit_edge = Box::BELOW;  // the curve enters the bottom
    push_back( new_hit );
    new_hit.E_in = x;
    new_hit.hit_edge = Box::BOTTOM_IN;  // the curve enters the bottom
    push_back( new_hit );
    new_hit.E_in = E_in_right;  // the curve xits on the right side
    new_hit.hit_edge = Box::INSIDE;  // the curve enters on the left side
    push_back( new_hit );
  }
  else // the curve misses the box
  {
    new_hit.E_in = E_in_left;
    new_hit.hit_edge = Box::BELOW;
    push_back( new_hit );
    new_hit.E_in = E_in_right;
    push_back( new_hit );
  }
}

// **************** class Coh::coherent ******************
// ---------------- Coh::coherent::get_T ------------------
// Calculates the transfer matrix fror this particle.
void Coh::coherent::get_T( Trf::T_matrix& transfer, Ddvec::dd_vector& sigma )
{
  // :::::::: kludge ::::::::::::::::
  // scale x from 1/(wave length) to energy
  double x_to_energy = Global.Value( "x_to_energy" );
  file_data.scale_x_to_energy( x_to_energy );
  // :::::::: end kludge ::::::::::::::::

  // ensure that the data starts at x=0 and extrapolate anomalous
  double last_Ein = transfer.in_groups.get_top_E( );
  extrapolate_data( last_Ein );

  // don't scatter to below the lowest energy group
  *transfer.out_groups.begin( ) = 0.0;
  // set up multiplicity and cross section
  Ddvec::dd_vector multiple;
  multiple.make_flat( transfer.in_groups, 1.0 );
  Ddvec::dd_vector weight;
  weight.make_flat( transfer.in_groups, 1.0 );

  // compute the cross section
  get_xsec( transfer.mu_quad_rule, sigma );
  transfer.threshold = sigma.begin( )->x;
  
  // number of 2-d quadratures
  long int quad_count = 0;
  // number of calls to coherent_F::Ein_F or coherent_F::sigma_F
  long int Ein_F_count = 0;
  long int mu_F_count = 0;  // number of calls to coherent_F::mu_F

  // loop over the incoming energy groups
#pragma omp parallel for schedule( dynamic, 1 ) default( none ) \
  shared( sigma, multiple, weight, transfer )   \
  reduction( +: quad_count ) reduction( +: Ein_F_count ) \
  reduction( +: mu_F_count )
  for( int Ein_count = 0; Ein_count < transfer.num_Ein_bins;
     ++Ein_count )
  {
    Coh::coherent_Ein_param Ein_param;
    Ein_param.mu_quad_rule = transfer.mu_quad_rule;
      Ein_param.mu_param.anomalous.init_Re_anomalous( realAnomalous );
      Ein_param.mu_param.anomalous.init_Im_anomalous( imaginaryAnomalous );

    // loop over the x-values
    for( unsigned int x_count = 0; x_count < this->size( ); ++x_count )
    {
      // set the initial pointers
      Ein_param.setup( sigma, multiple, weight, transfer.e_flux, transfer.in_groups );
      Ein_param.this_sigma = sigma.begin( );
      Ein_param.next_sigma = sigma.begin( );
        ++Ein_param.next_sigma;
      std::vector< double >::const_iterator Ein_ptr = transfer.in_groups.begin( ) +
         Ein_count;
      Ein_param.Ein_ptr = Ein_ptr;
      Ein_param.next_Ein = Ein_ptr;
      ++Ein_param.next_Ein;
      Ein_param.x_ptr = begin( ) + x_count;
      Ein_param.next_x = begin( ) + x_count + 1;

      if( *Ein_param.next_Ein <= Ein_param.x_ptr->x )
      {
        // do nothing---the energy is below the data
      }
      else
      {
        Eout_ladder( transfer, Ein_count, &Ein_param );
      }
    }
    quad_count += Ein_param.quad_count;
    Ein_F_count += Ein_param.Ein_F_count;
    mu_F_count += Ein_param.mu_F_count;
  }  // end of parallel loop

  // print the counts of function evaluations
  std::cout << "2d quadratures: " << quad_count << std::endl;
  std::cout << "coherent_F::Ein_F calls: " << Ein_F_count << std::endl;
  std::cout << "coherent_F::mu_F calls: " << mu_F_count << std::endl;
  std::cout << "average coherent_F::Ein_F calls: " << 1.0*Ein_F_count/quad_count << std::endl;
  std::cout << "average coherent_F::mu_F calls: " << 1.0*mu_F_count/Ein_F_count << std::endl;
}
// ---------------- Coh::coherent::extrapolate_data ------------------
// Extrapolates the scattering factor to 0 and anomalous to max_Ein
void Coh::coherent::extrapolate_data( double max_Ein )
{
  // find the lowest energy for which we have both anomalous data
  anomalous_min_Ein = realAnomalous.begin( )->x;
  if( imaginaryAnomalous.begin( )->x > anomalous_min_Ein )
  {
    anomalous_min_Ein = imaginaryAnomalous.begin( )->x;
  }
  // ensure that the data starts at x=0 and includes anomalous_min_Ein
  prepend( file_data, anomalous_min_Ein );

  // extrapolate the anomalous factors to high energy as zero
  Ddvec::dd_vector::const_iterator anomalousPtr = realAnomalous.end( );
  --anomalousPtr;
  if( anomalousPtr->x < max_Ein )
  {
    realAnomalous.add_entry( max_Ein, 0.0 );
  }
  anomalousPtr = imaginaryAnomalous.end( );
  --anomalousPtr;
  if( anomalousPtr->x < max_Ein )
  {
    imaginaryAnomalous.add_entry( max_Ein, 0.0 );
  }
}
// ---------------- Coh::coherent::Eout_ladder ------------------
// Climbs through the outgoing energy groups
void Coh::coherent::Eout_ladder( Trf::T_matrix& transfer, int Ein_count,
  Coh::coherent_Ein_param *Ein_param )
{
  // set up the geometry
  Ein_param->upper_hits.hit_box( Ein_param->x_ptr->x, *Ein_param->Ein_ptr,
		      *Ein_param->next_Ein );
  Ein_param->lower_hits.hit_box( Ein_param->next_x->x, *Ein_param->Ein_ptr,
		      *Ein_param->next_Ein );
  Ein_param->lower_hits.common_hits( Ein_param->upper_hits );
  // loop through the outgoing energies (column of transfer)
  for( int Eout_count = 0; Eout_count < transfer.num_Eout_bins;
    ++Eout_count )
  {
    Ein_param->Eout_ptr = transfer.out_groups.begin( )
      + Eout_count;
    Ein_param->next_Eout = Ein_param->Eout_ptr;
    ++Ein_param->next_Eout;
    // sets the range of integration for this incident energy bin
    set_Ein_range( Ein_param );
    Ein_param->set_sigma_range( );
    if( *Ein_param->next_Eout <= *Ein_param->Ein_ptr )
    {
      continue;  // outgoing energy too low---go to the next
    }
    else if( *Ein_param->Eout_ptr >= *Ein_param->next_Ein )
    {
      break;  // outgoing energy too high
    }
    else
    {
      one_Ebox( transfer, Ein_count, Eout_count, Ein_param );
    }
  }
}
// ----------- Coh::coherent::one_Ebox --------------
// Integrate over one x-E box; loop over the (E_in, mu) coherent region
void Coh::coherent::one_Ebox( Trf::T_matrix& transfer, int Ein_count, int Eout_count,
    Coh::coherent_Ein_param *Ein_param  )
{
  // integrate depending on how the curves x = const meet the box
  Coh::coherent_hit_list::iterator low_hit_ptr = Ein_param->lower_hits.begin( );
  Coh::coherent_hit_list::iterator next_low_ptr = low_hit_ptr;
  ++next_low_ptr;
  Coh::coherent_hit_list::iterator high_hit_ptr = Ein_param->upper_hits.begin( );
  Coh::coherent_hit_list::iterator next_high_ptr = high_hit_ptr;
  ++next_high_ptr;
  for( ; ( next_low_ptr != Ein_param->lower_hits.end( ) ) &&
         ( next_high_ptr != Ein_param->upper_hits.end( ) );
       low_hit_ptr = next_low_ptr, ++next_low_ptr,
         high_hit_ptr = next_high_ptr, ++next_high_ptr )
  {
    if( low_hit_ptr->E_in >= Ein_param->data_E_1 )
    {
      // the incident energy is too high
      break;
    }
    if( next_low_ptr->E_in <= Ein_param->data_E_0 )
    {
      // do nothing---the incident energy is too low
      continue;
    }
    Ein_param->Ein_0 = ( low_hit_ptr->E_in >= Ein_param->data_E_0 ) ?
      low_hit_ptr->E_in : Ein_param->data_E_0;
    Ein_param->Ein_1 = ( next_low_ptr->E_in <= Ein_param->data_E_1 ) ?
      next_low_ptr->E_in : Ein_param->data_E_1;
    // the lower limit of mu integration could be mu = -1
    Ein_param->use_mu_minus1 =
      ( ( next_low_ptr->hit_edge == Box::BOTTOM_IN ) ||
	( next_low_ptr->hit_edge == Box::BELOW ) );
    // the upper limit of mu integration could be mu = 1
    Ein_param->use_mu1 = ( Ein_param->x_ptr->x == 0.0 );
    update_T( transfer, Ein_count, Eout_count, Ein_param );
  }
}
// ---------------- Coh::coherent::set_Ein_range ------------------
// Sets the range of incident energies for this intergration
void Coh::coherent::set_Ein_range( Coh::coherent_Ein_param *Ein_param )
{
  // sets the range based on the flux data
  Ein_param->set_Ein_range( );
  if( *Ein_param->Ein_ptr > Ein_param->data_E_0 )
  {
    Ein_param->data_E_0 = *Ein_param->Ein_ptr;
  }
  if( *Ein_param->Eout_ptr > Ein_param->data_E_0 )
  {
    Ein_param->data_E_0 = *Ein_param->Eout_ptr;
  }
  if( Ein_param->x_ptr->x > Ein_param->data_E_0 )
  {
    Ein_param->data_E_0 = Ein_param->x_ptr->x;
  }
  if( *Ein_param->next_Ein < Ein_param->data_E_1 )
  {
    Ein_param->data_E_1 = *Ein_param->next_Ein;
  }
  if( *Ein_param->next_Eout < Ein_param->data_E_1 )
  {
    Ein_param->data_E_1 = *Ein_param->next_Eout;
  }
}
// ---------------- Coh::coherent::update_T ------------------
// Adds the result of one integration
void Coh::coherent::update_T( Trf::T_matrix &transfer, int Ein_count, int Eout_count,
    Coh::coherent_Ein_param *Ein_param )
{
  static double tol = Global.Value( "quad_tol" );
  // a vector to store the integrals
  Coef::coef_vector value( transfer.order, transfer.conserve );
  value.set_zero( );

  // parameters for the integration
  Qparam::QuadParamBase *params = static_cast< Qparam::QuadParamBase* >( Ein_param );

  // Do the data overlap this energy bin?
  bool bin_OK = Ein_param->start_sigma( );
  if( !bin_OK )
  {
    return;
  }

  // loop over the cross section and anomalous data
  bool done = false;
  while( !done )
  {
    Ein_param->get_range( );

    // evaluate the integral
    quad_F::integrate( coherent_F::Ein_F, transfer.Ein_quad_rule, Ein_param->left_E,
		       Ein_param->right_E, params, tol, &value );
    // add this integral
    transfer( Ein_count, Eout_count ) += value;
    // increment the function counts
    Ein_param->Ein_F_count += Ein_param->func_count;
    ++Ein_param->quad_count;

    // the next subinterval
    done = Ein_param->next_range( );
  }
}
// ---------------- Coh::coherent::get_xsec ------------------
// Calculates the cross section
void Coh::coherent::get_xsec( Qmeth::Quadrature_Rule mu_quad_rule, Ddvec::dd_vector& sigma )
{
  // the cutoff for the use of weight sqrt(x)
  static double cutoff = 1.0 - Global.Value( "sqrt_wt_cutoff" );
  Qmeth::Quadrature_Rule wt_mu_rule;
  wt_mu_rule.quad_method = Qmeth::WEIGHT_L1;
  wt_mu_rule.adaptive = mu_quad_rule.adaptive;
  wt_mu_rule.input_set = false;
  
  static double tol = Global.Value( "quad_tol" );

  // storage for quadrature
  Coef::coef_vector value( 0, Coef::NUMBER );
  // parameters for the integration
  Coh::coherent_mu_param mu_params;
  mu_params.anomalous.init_Re_down( realAnomalous );
  mu_params.anomalous.init_Im_down( imaginaryAnomalous );
  Qparam::QuadParamBase *params = static_cast< Qparam::QuadParamBase* >( &mu_params );

  // As incident energies use the x-values with mu = -1
  // But start at the lowest anomalous data energy
  static double e_tol = Global.Value( "tight_tol" );
  Coh::coherent::const_iterator this_entry = begin( );
  while( this_entry->x < anomalous_min_Ein*( 1.0 - e_tol ) )
  {
    ++this_entry;
  }
  for( ; this_entry != end( ); ++this_entry )
  {
    sigma.add_entry( this_entry->x, 0.0 );
  }
  long int sigma_quad_count = 0;  // number of quadratures for cross section
  long int sigma_F_count = 0;     // number of calls to coherent_F::sigma_F

  // loop downward through the incident energies
  Ddvec::dd_vector::iterator this_sigma = sigma.end( );
  --this_sigma;

  Coh::coherent::const_iterator this_x_entry;
  Coh::coherent::const_iterator next_x_entry;
  Coh::coherent::const_iterator last_x_entry = end( );
  --last_x_entry;

  while( this_sigma->x > 0.0 )
  {
    mu_params.this_sigma = this_sigma;
    // loop downward through the x-entries
    next_x_entry = last_x_entry;
    this_x_entry = next_x_entry;
    --this_x_entry;
    for( ; next_x_entry != begin( ); --next_x_entry, --this_x_entry )
    {
      // parameters for the integration over mu
      double E_in = this_sigma->x;
      mu_params.set_E_in( E_in );
      mu_params.set_pair( *this_x_entry, *next_x_entry );
      mu_params.anomalous.get_next_Re_down( E_in );
      mu_params.anomalous.get_next_Im_down( E_in );

      // mu decreases as x increases
      double mu_0 = x_vector_F::get_mu_from_x( next_x_entry->x, this_sigma->x );
      double mu_1 = x_vector_F::get_mu_from_x( this_x_entry->x, this_sigma->x );
      if( mu_0 >= cutoff )
      {
	quad_F::integrate( coherent_F::flip_sigma_F, wt_mu_rule,
			   1 - mu_1, 1 - mu_0,
  			 params, tol, &value );
      }
      else if( mu_1 <= cutoff )
      {
        quad_F::integrate( coherent_F::sigma_F, mu_quad_rule,
			 mu_0, mu_1, params, tol, &value );
      }
      else
      {
        quad_F::integrate( coherent_F::sigma_F, mu_quad_rule,
			 mu_0, cutoff, params, tol, &value );
	this_sigma->y += value.weight_1[0];
        // increment the function counts
	sigma_F_count += mu_params.func_count;
	++sigma_quad_count;
	  
        quad_F::integrate( coherent_F::flip_sigma_F, wt_mu_rule,
			   1 - mu_1, 1 - cutoff,
  			 params, tol, &value );
      }
      
      this_sigma->y += value.weight_1[0];
      // increment the function counts
      sigma_F_count += mu_params.func_count;
      ++sigma_quad_count;
    }
    if( this_sigma == sigma.begin( ) )
    {
      break;
    }
    --this_sigma;
    --last_x_entry;
  }
  // special for x = 0
  this_sigma = sigma.begin( );
  this_x_entry = begin( );
  static double Thompson = Global.Value( "Thompson" );  // Thompson cross section
  if( ( this_sigma->x == 0.0 ) && ( this_x_entry->x == 0.0 ) )
  {
    double FF = this_x_entry->y;
    this_sigma->y = Thompson * FF * FF;
  }
  // print the counts of function evaluations
  std::cout << "quadratures for cross section: " << sigma_quad_count << std::endl;
  std::cout << "coherent_F::sigma_F calls: " << sigma_F_count << std::endl;
  std::cout << "average coherent_F::sigma_F calls: " << 1.0*sigma_F_count/sigma_quad_count << std::endl;

  // use lin-lin interpolation
  sigma.interp_type = Terp::LINLIN;
}

// **************** functions to integrate **********
// ----------------  coherent_F::sigma_F ------------------------
// Function for integrating the cross section
bool coherent_F::sigma_F( double mu, Qparam::QuadParamBase *FF_param,
			  Coef::coef_vector *value )
{
  // the parameters are really Coh::coherent_mu_param
  Coh::coherent_mu_param *params = static_cast<Coh::coherent_mu_param*>( FF_param );
  params->func_count += 1;

  // get the coherent scattering cross section
  bool interp_OK;
  double xs = params->get_sigma( mu, &interp_OK );

  value->weight_1[0] = xs;
  return interp_OK;
}
// ----------------  coherent_F::flip_sigma_F ------------------------
// Function for integrating the cross section with singularity sqrt( flip_mu )
bool coherent_F::flip_sigma_F( double flip_mu, Qparam::QuadParamBase *FF_param,
			       Coef::coef_vector *value )
{
  return coherent_F::sigma_F( 1 - flip_mu, FF_param, value );
}
// ----------------  coherent_F::mu_F ------------------------
// Function for the 1-d quadrature over mu
bool coherent_F::mu_F( double mu, Qparam::QuadParamBase *FF_param,
		       Coef::coef_vector *value )
{
  // the parameters are really Coh::coherent_mu_param
  Coh::coherent_mu_param *params = static_cast<Coh::coherent_mu_param*>( FF_param );
  params->func_count += 1;

  // get the coherent scattering cross section
  bool interp_OK;
  double xs = params->get_sigma( mu, &interp_OK );

  // the Legendre polynomial
  math_F::Legendre( mu, value );
  *value *= xs;
  // do the energy weighting if necessary
  if( ( value->conserve == Coef::ENERGY ) || ( value->conserve == Coef::BOTH ) )
  {
    // E_out =  params->get_E_in
    value->scale_E( params->get_E_in() );
  }

  return interp_OK;
}
// ----------------  coherent_F::flip_mu_F ------------------------
// Function for the 1-d quadrature over mu
bool coherent_F::flip_mu_F( double flip_mu, Qparam::QuadParamBase *FF_param,
			    Coef::coef_vector *value )
{
  return coherent_F::mu_F( 1 - flip_mu, FF_param, value );
}
// ----------------  coherent_F::Ein_F ------------------------
// Function for the 2-d quadrature over E_in
bool coherent_F::Ein_F( double E_in, Qparam::QuadParamBase *coherent_param,
   Coef::coef_vector *value )
{
  // the parameters are really Coh::coherent_Ein_param *
  Coh::coherent_Ein_param *e_params =
    static_cast<Coh::coherent_Ein_param *>( coherent_param );
  e_params->func_count += 1;

  // The value of coherent_Ein_F is itself an integral over mu.
  // *value comes in as 0.  

  // parameters for the integration over mu
  e_params->mu_param.func_count = 0;
  e_params->mu_param.set_E_in( E_in );
  e_params->mu_param.set_pair( *e_params->x_ptr, *e_params->next_x );
  // the range of integration
  double mu_0 = ( e_params->use_mu_minus1 ) ? -1.0 : e_params->bottom_mu( E_in );
  double mu_1 = ( e_params->use_mu1 ) ? 1.0 : e_params->top_mu( E_in );

  // the cutoff for the use of weight sqrt(x)
  static double cutoff = 1.0 - Global.Value( "sqrt_wt_cutoff" );
  Qmeth::Quadrature_Rule wt_mu_rule;
  wt_mu_rule.quad_method = Qmeth::WEIGHT_L1;
  wt_mu_rule.adaptive = true;
  wt_mu_rule.input_set = false;
  
  Coef::coef_vector temp_value( value->order, value->conserve );
  temp_value.set_zero( );

  // evaluate the integral over mu
  Qparam::QuadParamBase *params =
    static_cast< Qparam::QuadParamBase* >( &e_params->mu_param );
  static double tol = Global.Value( "quad_tol" );
  bool is_OK;

  if( mu_0 >= cutoff )
  {
    is_OK = quad_F::integrate( coherent_F::flip_mu_F, wt_mu_rule,
			   1 - mu_1, 1 - mu_0,
  			 params, tol, value );
  }
  else if( mu_1 <= cutoff )
  {
    is_OK = quad_F::integrate( coherent_F::mu_F, e_params->mu_quad_rule,
			 mu_0, mu_1, params, tol, value );
  }
  else
  {
    bool is_OK1 = quad_F::integrate( coherent_F::mu_F, e_params->mu_quad_rule,
			 mu_0, cutoff, params, tol, &temp_value );
    // increment the function count
    e_params->mu_F_count += e_params->mu_param.func_count;
	  
    bool is_OK2 = quad_F::integrate( coherent_F::flip_mu_F, wt_mu_rule,
			   1 - mu_1, 1 - cutoff,
  			 params, tol, value );
    *value += temp_value;
    is_OK = is_OK1 && is_OK2;
  }
    
  e_params->mu_F_count += e_params->mu_param.func_count;
  // weight it by the flux
  e_params->flux_weight( E_in );
  *value *= e_params->current_weight;

  return is_OK;
}
