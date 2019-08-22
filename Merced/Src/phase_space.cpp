/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2010-10-01 (Mon, Oct 1, 2010) $
 * $Author: hedstrom $
 * $Id: phase_space.cpp 1 2010-10-01 hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// implementation of the classes used to handle phase-space energy model


#include <cmath>
#ifdef _OPENMP
 #include <omp.h>
#endif

#include "phase_space.hpp"
#include "math_util.hpp"
#include "messaging.hpp"
#include "global_params.hpp"

// ************* class phase_space_Ein_param *****************
// ----------- phase_space_Ein_param::set_Ecm_param --------------
// Sets up the parameters for integration over center-of-mass outgoing energy
void phase_space_Ein_param::set_Ecm_param( double E_in )
{
  Ecm_params.setup( E_in, Eout_min, Eout_max );
  Ecm_params.mu_quad_method = mu_quad_method;
}
// ----------- phase_space_Ein_param::setup_map --------------
// Copies the data for mapping to the lab frame
void phase_space_Ein_param::setup_map( phase_space_map *PSmap )
{
  Ecm_params.mu_quad_method = mu_quad_method;
  map = PSmap;
  Ecm_params.map = map;
  // use the same mass parameters for the phase-space functions
  Ecm_params.PSmap = map;
  // the Vcm_Vlab_hit_list objects need the gamma for the energy of translation of the center of mass
  lower_hits.G0_data.gamma = PSmap->gamma;
  lower_hits.G1_data.gamma = PSmap->gamma;
  upper_hits.G0_data.gamma = PSmap->gamma;
  upper_hits.G1_data.gamma = PSmap->gamma;
}

// ************* class phase_space_Ecm_param *****************
// ----------- phase_space_Ecm_param::setup --------------
// Sets up the data for this incident energy
void phase_space_Ecm_param::setup( double E_in, double Eoutmin, double Eoutmax )
{
  double max_Ecm_out = PSmap->get_Ecm_max( E_in );
  Ecm_Elab_Ecm_param::setup( E_in, Eoutmin, Eoutmax, 0.0, max_Ecm_out );
  V_lab_sectors( );
}

// ************ class phase_space **********************
// ----------- phase_space::check_input ------------------
// Ensures that the input data is complete
void phase_space::check_input( )
{
  if( mProj < 0.0 )
  {
    FatalError( "phase_space::check_input", "Projectile mass not set" );
  }
  if( mTarg < 0.0 )
  {
    FatalError( "phase_space::check_input", "Target mass not set" );
  }
  if( mEject < 0.0 )
  {
    FatalError( "phase_space::check_input", "Ejected particle mass not set" );
  }
  if( totalMass < 0.0 )
  {
    FatalError( "phase_space::check_input", "Total particle mass not set" );
  }
  if( Q_value < -1.0e10 )
  {
    FatalError( "phase_space::check_input", "Reaction Q value not set" );
  }
  if( numParticles < 0 )
  {
    FatalError( "phase_space::check_input", "Number of emitted particles not set" );
  }
}
// ----------- phase_space::get_Ein_range --------------
//  Gets the range of nontrivial incident energy bins; computes first_Ein and last_Ein
// returns true if the threshold is too high for the energy bins
bool phase_space::get_Ein_range( const dd_vector& sigma, const dd_vector& multiple,
    const dd_vector& weight,
    const Flux_List& e_flux, const Energy_groups& Ein_groups )
{
  double E_first;
  double E_last;
  phase_space_Ein_param initial_param;
  bool done = initial_param.get_Ein_range( sigma, multiple, weight, e_flux,
                                         Ein_groups, &E_first, &E_last );
  if( done ) return true;

  // don't try to work below the threshold of the reaction
  double threshold = PSmap.get_threshold( );
  if( threshold > E_first ) E_first = threshold;

  first_Ein = Ein_groups.first_bin_ID( E_first );
  last_Ein = Ein_groups.last_bin_ID( E_last );

  return false;
}
// ----------- phase_space::get_T ------------------
// Calculates the transfer matrix for this particle
void phase_space::get_T( const dd_vector& sigma, const dd_vector& multiple, 
  const dd_vector& weight, T_matrix& transfer )
{
  if( mEject == 0.0 )
  {
    FatalError( "phase_space::get_T", "Photon emission not implemented" );
  }
  check_input( );

  // for center-of-mass data
  setup_map( );
  bool done = get_Ein_range( sigma, multiple, weight, transfer.e_flux,
    transfer.in_groups );
  if( done )
  {
    Info( "phase_space::get_T", "No data in energy range of transfer matrix" );
    transfer.zero_transfer( );
  }

  long int quad_count = 0;  // quadratures
  long int Ein_F_count = 0;  // number of phase_space_F::Ein_F calls
  long int Ecm_F_count = 0;  // number of phase_space_F::Ecm_F calls
  long int mu_F_count = 0;   // number of phase_space_F::mu_F calls

  // now do the integrals incident bin by incident bin
#pragma omp parallel for schedule( dynamic, 1 ) default( none )	\
  shared( sigma, multiple, weight, transfer ) \
  reduction( +: quad_count ) reduction( +: Ein_F_count ) \
  reduction( +: Ecm_F_count ) reduction( +: mu_F_count )
  for( int Ein_bin = first_Ein; Ein_bin < last_Ein; ++Ein_bin )
  {
    phase_space_Ein_param Ein_param;
    Ein_param.setup_map( &PSmap );
    Ein_param.Eout_quad_method = transfer.Eout_quad_method;  // quadrature method for outgoing energy
    Ein_param.mu_quad_method = transfer.mu_quad_method;  // quadrature method for outgoing cosine
    // set up the data range for this bin
    Ein_param.setup_bin( Ein_bin, sigma, multiple, weight, transfer.e_flux,
                         transfer.in_groups );
    for( ; ; )
    {
      Ein_param.set_Ein_range( );   // get the incident energy interval
      Ein_param.lower_hits.G0_data.set_energies( Ein_param.data_E_0, 0.0 );
      double data_Ecm_max = PSmap.get_Ecm_max( Ein_param.data_E_0 );
      Ein_param.upper_hits.G0_data.set_energies( Ein_param.data_E_0, data_Ecm_max );
      Ein_param.lower_hits.G1_data.set_energies( Ein_param.data_E_1, 0.0 );
      data_Ecm_max = PSmap.get_Ecm_max( Ein_param.data_E_1 );
      Ein_param.upper_hits.G1_data.set_energies( Ein_param.data_E_1, data_Ecm_max );
      Eout_ladder( transfer, &Ein_param );  // loop over the outgoing lab energy bins
      bool Done = Ein_param.update_bin_pointers( Ein_param.data_E_1 );   // go to the next interval
      if( Done )
      {
        break;
      }
    }
    quad_count += Ein_param.quad_count;
    Ein_F_count += Ein_param.Ein_F_count;
    Ecm_F_count += Ein_param.Ecm_F_count;
    mu_F_count += Ein_param.mu_F_count;
  } // end of parallel loop

  // print the counts of function evaluations
  cout << "3d quadratures: " << quad_count << endl;
  cout << "phase_space_F::Ein_F calls: " << Ein_F_count << endl;
  cout << "phase_space_F::Ecm_F calls: " << Ecm_F_count << endl;
  cout << "phase_space_F::mu_F calls: " << mu_F_count << endl;
  cout << "average phase_space_F::Ein_F calls: " << 1.0*Ein_F_count/quad_count << endl;
  cout << "average phase_space_F::Ecm_F calls: " << 1.0*Ecm_F_count/Ein_F_count << endl;
  cout << "average phase_space_F::mu_F calls: " << 1.0*mu_F_count/Ecm_F_count << endl;
}
// ----------- phase_space::setup_map --------------
// Sets up the map from center-of-mass to laboratory coordinates
void phase_space::setup_map( )
{
  PSmap.setup_ratios( mProj, mTarg, mEject );
  PSmap.set_data( numParticles, mEject, totalMass, Q_value );
}
// ----------- phase_space::Eout_ladder ------------------
// Loops over the outgoing lab energy bins for one pair of outgoing cm energies
void phase_space::Eout_ladder( T_matrix& transfer, phase_space_Ein_param *Ein_param )
{
  //  bool check_geometry = true;
  bool check_geometry = false;
  bool geom_OK;  // for checking the consistency of the geometry
  bool upper_hits_set = false;
  Vcm_Vlab_hit_list test_hits;
  double dummy = 0.0;
  double Ecm;
  int Eout_count = 0;
  vector< double >::const_iterator Eout_ptr = transfer.out_groups.begin( );
  vector< double >::const_iterator next_Eout = Eout_ptr;
  ++next_Eout;

  // Ein_param->lower_hits is for emission at 0 center-of-mass energy
  geom_OK = Ein_param->lower_hits.hit_box( dummy, Eout_ptr, Ein_param->data_E_0, Ein_param->data_E_1 );
  if( check_geometry )
  {
    cout << "lower_hits for Eout: " << *Eout_ptr << endl;
    Ein_param->lower_hits.print( );
  }

  // Check for only forward emission
  if( Ein_param->upper_hits.G1_data.E_cm < Ein_param->upper_hits.G1_data.get_Etrans( ) )
  {
    geom_OK = Ein_param->upper_hits.hit_box( dummy, Eout_ptr, Ein_param->data_E_0, Ein_param->data_E_1 );
    upper_hits_set = true;
    if( check_geometry )
    {
      cout << "Forward with next_Eout: " << *next_Eout << endl;
      Ein_param->upper_hits.print( );
    }
    if( !geom_OK )
    {
      Ecm = PSmap.get_Ecm_max( Ein_param->data_E_0 );
      test_hits.G0_data.set_energies( Ein_param->data_E_0, Ecm );
      Ecm = PSmap.get_Ecm_max( Ein_param->data_E_1 );
      test_hits.G1_data.set_energies( Ein_param->data_E_1, Ecm );
      geom_OK = test_hits.hit_box( dummy, Eout_ptr, Ein_param->data_E_0, Ein_param->data_E_1 );
      test_hits.print( );
      FatalError( "phase_space::Eout_ladder", "Check the coding, 1" );
    }
    while( Ein_param->upper_hits.is_above( ) )
    {
      Eout_ptr = next_Eout;
      ++next_Eout;
      if( next_Eout == transfer.out_groups.end( ) )
      {
        return;
      }
      ++Eout_count;
      geom_OK = Ein_param->upper_hits.hit_box( dummy, Eout_ptr, Ein_param->data_E_0, Ein_param->data_E_1 );
      if( check_geometry )
      {
        cout << "next Forward with next_Eout: " << *next_Eout << endl;
        Ein_param->upper_hits.print( );
      }
      if( !geom_OK )
      {
        Ecm = PSmap.get_Ecm_max( Ein_param->data_E_0 );
        test_hits.G0_data.set_energies( Ein_param->data_E_0, Ecm );
        Ecm = PSmap.get_Ecm_max( Ein_param->data_E_1 );
        test_hits.G1_data.set_energies( Ein_param->data_E_1, Ecm );
        geom_OK = test_hits.hit_box( dummy, Eout_ptr, Ein_param->data_E_0, Ein_param->data_E_1 );
        test_hits.print( );
        FatalError( "phase_space::Eout_ladder", "Check the coding, 2" );
      }
    }
  }
  // Now, compute integrals until the lab energy bin is above the E_cm data
  for( ; Eout_count < transfer.num_Eout_bins;
       ++Eout_count, Eout_ptr = next_Eout, ++next_Eout )
  {
    if( !upper_hits_set )
    {
      geom_OK = Ein_param->upper_hits.hit_box( dummy, Eout_ptr, Ein_param->data_E_0, Ein_param->data_E_1 );
      if( check_geometry )
      {
        cout << "upper_hits for Eout: " << *Eout_ptr << endl;
        Ein_param->upper_hits.print( );
      }
      if( !geom_OK )
      {
        Ecm = PSmap.get_Ecm_max( Ein_param->data_E_0 );
        test_hits.G0_data.set_energies( Ein_param->data_E_0, Ecm );
        Ecm = PSmap.get_Ecm_max( Ein_param->data_E_1 );
        test_hits.G1_data.set_energies( Ein_param->data_E_1, Ecm );
        geom_OK = test_hits.hit_box( dummy, Eout_ptr, Ein_param->data_E_0, Ein_param->data_E_1 );
        test_hits.print( );
        FatalError( "phase_space::Eout_ladder", "Check the coding, 3" );
      }
    }
    if( Ein_param->upper_hits.is_below( ) )
    {
      break;  // we are done
    }
    // integrate over this E-E' box
    one_Ebox( transfer, Eout_count, Ein_param );
    upper_hits_set = false;
  }
}
// -----------  phase_space::one_Ebox ------------------
// Does the integration for one Eout_lab annulus between a pair of incident energies
void phase_space::one_Ebox( T_matrix& transfer, int Eout_count,
   phase_space_Ein_param *Ein_param )
{
  // the E' energy range
  Ein_param->Eout_min = transfer.out_groups[ Eout_count ];
  Ein_param->Eout_max = transfer.out_groups[ Eout_count + 1 ];
  //  cout << Ein_param->Eout_min << " < E_out < " << Ein_param->Eout_max << endl;
  // set up common incident energies
  Ein_param->lower_hits.common_hits( Ein_param->upper_hits );
  //  cout << "in one_Ebox upper_hits for Eout: " << Ein_param->Eout_min << endl;
  //  Ein_param->upper_hits.print( );
  //  cout << "in one_Ebox lower_hits for Eout: " << Ein_param->Eout_min << endl;
  //  Ein_param->lower_hits.print( );

  // integrate depending on how the arcs E_cm = const meet the box
  Vcm_Vlab_hit_list::iterator low_hit_ptr = Ein_param->lower_hits.begin( );
  Vcm_Vlab_hit_list::iterator next_low_ptr = low_hit_ptr;
  ++next_low_ptr;
  Vcm_Vlab_hit_list::iterator high_hit_ptr = Ein_param->upper_hits.begin( );
  Vcm_Vlab_hit_list::iterator next_high_ptr = high_hit_ptr;
  ++next_high_ptr;
  for( ; ( next_low_ptr != Ein_param->lower_hits.end( ) ) &&
         ( next_high_ptr != Ein_param->upper_hits.end( ) );
       low_hit_ptr = next_low_ptr, ++next_low_ptr,
         high_hit_ptr = next_high_ptr, ++next_high_ptr )
  {
    if( ( ( low_hit_ptr->hit_edge == ABOVE ) &&
          ( high_hit_ptr->hit_edge == ABOVE ) ) ||
        ( ( next_low_ptr->hit_edge == ABOVE ) &&
          ( next_high_ptr->hit_edge == ABOVE ) ) ||
        ( ( low_hit_ptr->hit_edge == ABOVE_FORWARD ) &&
          ( high_hit_ptr->hit_edge == ABOVE_FORWARD ) ) ||
        ( ( next_low_ptr->hit_edge == ABOVE_FORWARD ) &&
          ( next_high_ptr->hit_edge == ABOVE_FORWARD ) ) ||
        ( ( low_hit_ptr->hit_edge == BELOW ) &&
          ( high_hit_ptr->hit_edge == BELOW ) ) ||
        ( ( next_low_ptr->hit_edge == BELOW ) &&
          ( next_high_ptr->hit_edge == BELOW ) ) )
    {
      continue;
    }
    // the range of integration in incident energy
    Ein_param->Ein_0 = low_hit_ptr->E_in;
    Ein_param->Ein_1 = next_low_ptr->E_in;
    update_T( transfer, Eout_count, Ein_param );
  }
}
// -----------  phase_space::update_T ------------------
// Adds to an element of transfer the integral between the intersections of 2 Eout_cm = const arcs with the Eout_lab box
  void phase_space::update_T( T_matrix &transfer, int Eout_count,
   phase_space_Ein_param *Ein_param )
{
  static double tol = Global.Value( "quad_tol" );
  // a vector to store the integrals
  coef_vector value( transfer.order, transfer.conserve );
  value.set_zero( );

  // parameters for the integration
  QuadParamBase *params = static_cast< QuadParamBase* >( Ein_param );

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
    // evaluate the integral
    quad_F::integrate( phase_space_F::Ein_F, transfer.Ein_quad_method, left_E, right_E,
		       params, tol, &value );

    // add this integral
    transfer( Ein_param->Ein_count, Eout_count ) += value;
    // increment the function counts
    Ein_param->Ein_F_count += Ein_param->func_count;
    Ein_param->quad_count += Ein_param->Vcm_hit_count;
  }
}
// -----------  phase_space::copy_masses ------------------
// Stores the masses
void phase_space::copy_masses( const particleInfo &particle_info )
{
  mProj = particle_info.mProj;
  mTarg = particle_info.mTarg;
  mEject = particle_info.mProd;
}

// **************** Functions to integrate *********************
// --------------------  phase_space_F::mu_F ------------------
// Function for the 1-d quadrature over cm cosine
void phase_space_F::mu_F( double mu, QuadParamBase *mu_quad_param, coef_vector *value )
{
  // the parameters are really Ecm_Elab_mu_param
  Ecm_Elab_mu_param *mu_params = static_cast< Ecm_Elab_mu_param* >( mu_quad_param );
  mu_params->func_count += 1;

  double Eout_lab;
  double mu_lab;
  mu_params->map->get_E_mu_lab( mu_params->E_in, mu_params->Eout_cm, mu, &Eout_lab,
    &mu_lab );

  // the Legendre polynomials
  math_F::Legendre( mu_lab, value );
  *value *= mu_params->Ecm_prob;

  // do the energy weighting if necessary
  if( ( value->conserve == ENERGY ) || ( value->conserve == BOTH ) )
  {
    value->scale_E( Eout_lab );
  }
}
// ------------------- phase_space_F::Ecm_F ------------------
// Function for the 2-d quadrature over cm cosine and Eout_cm
void phase_space_F::Ecm_F( double Eout_cm, QuadParamBase *Ecm_quad_param, coef_vector *value )
{
  // the parameters are really phase_space_Ecm_param *
  phase_space_Ecm_param *Ecm_param = static_cast<phase_space_Ecm_param *>( Ecm_quad_param );
  Ecm_param->func_count += 1;
  if( Ecm_param->func_count % 500 == 0 )
  {
    Info( "phase_space_F::Ecm_F", pastenum( "got ", Ecm_param->func_count ) + " evaluations");
    cout << "lab_Eout_max: " << Ecm_param->lab_Eout_max <<
      " E_in: " << Ecm_param->E_in <<
      " Ecm_max: " << Ecm_param->Ecm_max << endl;
  }

  // The value of phase_space_Ecm_F is itself an integral over cm cosine.
  // *value comes in as 0.  

  // parameters for the integration over cm cosine
  Ecm_Elab_mu_param mu_param;
  mu_param.setup( Ecm_param->E_in, Eout_cm, *Ecm_param );
  mu_param.Ecm_prob = Ecm_param->get_Ecm_prob( Eout_cm, Ecm_param->data_Ecm_max );

  // evaluate the integral over eta
  QuadParamBase *params = static_cast< QuadParamBase* >( &mu_param );
  static double tol = Global.Value( "quad_tol" );
  quad_F::integrate( phase_space_F::mu_F, Ecm_param->mu_quad_method, mu_param.mu_cm_min,
		     mu_param.mu_cm_max, params, tol, value );
  Ecm_param->mu_F_count += mu_param.func_count;
}
// ------------------- phase_space_F::Ein_F ------------------
// Function for the 3-d quadrature over E_in, and Eout_cm and cm cosine
// The value of phase_space_Ein_F is itself an integral over Eout_cm and cm cosine.
void phase_space_F::Ein_F( double E_in, QuadParamBase *Ein_quad_param, coef_vector *value )
{
  //  bool check_sectors = true;
  bool check_sectors = false;
  value->set_zero( );  // initialize to zero
  // the parameters are really phase_space_Ein_param *
  phase_space_Ein_param *Ein_param = static_cast<phase_space_Ein_param *>( Ein_quad_param );
  Ein_param->Vcm_hit_count = 0;   // number of local calls to quad_F::integrate

  // set up parameters for the integration over Eout_cm and cm cosine
  Ein_param->set_Ecm_param( E_in );
  QuadParamBase *params = static_cast< QuadParamBase* >( &Ein_param->Ecm_params );

  // Integrate over sectors of ( Eout_cm, mu_cm ) space
  coef_vector one_value( value->order, value->conserve );

  list< Vcm_quadBox_Hit >::const_iterator this_V_hit = 
    Ein_param->Ecm_params.V_cm_limits.begin( );
  list< Vcm_quadBox_Hit >::const_iterator next_V_hit = this_V_hit;
  ++next_V_hit;
  //  this_V_hit->print();
  for( ; next_V_hit != Ein_param->Ecm_params.V_cm_limits.end( );
         this_V_hit = next_V_hit, ++next_V_hit )
  {
    if( ( next_V_hit->V_cm <= Ein_param->Ecm_params.min_V_cm ) ||
        ( this_V_hit->hit_corner == V_BELOW ) )
    {
      continue;  // current V_cm values are below
    }
    else if( ( this_V_hit->V_cm >= Ein_param->Ecm_params.max_V_cm ) ||
             ( next_V_hit->hit_corner == V_ABOVE ) )
    {
      break;  // all remaining V_cm values are above
    }
    else
    {
      Ein_param->Vcm_hit_min = *this_V_hit;
      Ein_param->Vcm_hit_max = *next_V_hit;
    }
    if( check_sectors )
    {
      cout << "integrate V_lab" << endl;
      Ein_param->Vcm_hit_min.print( );
      Ein_param->Vcm_hit_max.print( );
    }
    Ein_param->Ecm_params.min_hit_corner = Ein_param->Vcm_hit_min.hit_corner;
    Ein_param->Ecm_params.max_hit_corner = Ein_param->Vcm_hit_max.hit_corner;
    Ein_param->Ecm_params.mu_F_count = 0;
    double tol = Ein_param->Ecm_params.Ecm_range( );

    quad_F::integrate( phase_space_F::Ecm_F, Ein_param->Eout_quad_method,
                       Ein_param->Ecm_params.Ecm_min,
		       Ein_param->Ecm_params.Ecm_max, params, tol, &one_value );
    *value += one_value;
    // we actually want to count the number of 3d integrals
    Ein_param->Vcm_hit_count += 1;
    Ein_param->func_count += 1;
    Ein_param->Ecm_F_count += Ein_param->Ecm_params.func_count;
    Ein_param->mu_F_count += Ein_param->Ecm_params.mu_F_count;
    if( Ein_param->func_count % 500 == 0 )
    {
      Info( "phase_space_F::Ein_F", pastenum( "got ", Ein_param->func_count ) +
        " evaluations");
    }
  }
  // weight it by flux * cross section
  Ein_param->set_weight( E_in );
  *value *= Ein_param->current_weight;
  //  cout << "E_in: " << E_in << " eta_0: " << eta_0 << " eta_1: " <<
  //    eta_1 << endl;
  //  value->print( );
}
