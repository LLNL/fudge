/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2013-12-10 11:06:56 -0800 (Tue, 10 Dec 2013) $
 * $Author: hedstrom $
 * $Id: gamma2body.cpp 1 2013-12-10 11:06:56 -0800 hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */
// implement the classes used for Legendre expansions for gammas from capture reactions

#ifdef _OPENMP
 #include <omp.h>
#endif

#include "gamma2body.hpp"
#include "global_params.hpp"
#include "messaging.hpp"

using namespace std;

// ************* red_blue_shift_map *************
// Gets the turn-over for backward emission
// ---------------- red_blue_shift_map::get_turnover ------------------
void red_blue_shift_map::get_turnover( dd_entry below, dd_entry above )
{
  static double E_tol = Global.Value( "E_tol" );
  // set up the void* parameters for the root finder
  void *params = static_cast<void *>( this );
  // get the derivative at the limits of the range
  below.y = Capture_Gamma_F::d_red_shift( below.x, params );
  above.y = Capture_Gamma_F::d_red_shift( above.x, params );
  if( below.y >= 0.0 )  // may happen for an endothermic reaction
  {
    if( above.y <= 0.0 )  // really weird
    {
      FatalError( "red_blue_shift_map::get_turnover",
         "check Capture_Gamma_F::d_red_shift" );
    }
    else
    {
      turnover_T = below.x;  // backward T_lab_out always increasing
    }
  }
  else if( above.y <= 0.0 )  // not expected
  {
    turnover_T = above.x;  // backward T_lab_out always decreasing
  }
  else
  {
    turnover_T = math_F::zeroin( Capture_Gamma_F::d_red_shift, 0.0, below,
      above, params, E_tol );
  }
  // the gamma energy at the minimum
  min_E_back = Capture_Gamma_F::red_shift( turnover_T, params );
}

// ************* capture_gamma_Ein_param *************
// ---------------- capture_gamma_Ein_param::set_map ------------------
// Sets up the mapping from center-of-mass to lab frame
void capture_gamma_Ein_param::set_map( red_blue_shift_map &map_ )
{
  map.rest_masses = map_.rest_masses;
  map.turnover_T = map_.turnover_T;
  map.min_E_back = map_.min_E_back;
  backward_hits.map = &map;
  forward_hits.map = &map;
}

// ************* capture_gamma_mu_param *************
// ---------------- capture_gamma_mu_param::interpolate ------------------
// Interpolates between two incident energies
void capture_gamma_mu_param::interpolate( double Ein,
   capture_gamma::const_iterator prev_coefs,
   capture_gamma::const_iterator next_coefs )
{
  if( Ein_interp == LINLIN )
  {
    coefs.linlin_interp( Ein, *prev_coefs, *next_coefs );
  }
  else
  {
    coefs.linlog_interp( Ein, *prev_coefs, *next_coefs );
  }
}

// ********* class gamma_hit_list *********
// ----------- gamma_hit_list::find_decreasing_hit --------------
// Finds the incident energy producing a backward gamma with energy lab_E_out.
// Use this routine when the lab-frame gamma energy decreases in incident energy.
// \param lab_E_out target energy for the outgoing backward gamma
// \param low_Ein pair ( lower incident energy, lab-frame gamma energy )
// \param high_Ein pair ( higher incident energy, lab-frame gamma energy )
double gamma_hit_list::find_decreasing_hit( double lab_E_out,
  const dd_entry &low_Ein, const dd_entry &high_Ein )
{
  // for the root finder
  static double E_tol = Global.Value( "E_tol" );
  // set up the void* parameters for the root finder
  void *params = static_cast<void *>( map );

  // This routine is used only for backward gammas, mu_cm = -1
  double this_Ein = math_F::zeroin( Capture_Gamma_F::red_shift, lab_E_out,
      low_Ein, high_Ein, params, E_tol );
  return this_Ein;
}
// ----------- gamma_hit_list::find_increasing_hit --------------
// Finds the incident energy producing a gamma with energy lab_E_out
// Use this routine when the lab-frame gamma energy increases in incident energy.
// \param lab_E_out target energy for the outgoing backward gamma
// \param low_Ein pair ( lower incident energy, lab-frame gamma energy )
// \param high_Ein pair ( higher incident energy, lab-frame gamma energy )
double gamma_hit_list::find_increasing_hit( double lab_E_out,
  const dd_entry &low_Ein, const dd_entry &high_Ein )
{
  // for the root finder
  static double E_tol = Global.Value( "E_tol" );
  // set up the void* parameters for the root finder
  void *params = static_cast<void *>( map );

  double this_Ein;
  if( eta < 0.0 )
  {
    this_Ein = math_F::zeroin( Capture_Gamma_F::red_shift, lab_E_out,
        low_Ein, high_Ein, params, E_tol );
  }
  else
  {
    this_Ein = math_F::zeroin( Capture_Gamma_F::blue_shift, lab_E_out,
        low_Ein, high_Ein, params, E_tol );
  }
  return this_Ein;
}
// ------------------ gamma_hit_list::test_lab_Eout --------------
// Finds the relation between the outgoing gamma energy and the value lab_Eout
// \param lab_Eout the desired lab-frame energy of the outgoing gamma
// \param Ein_hits the output pairs ( incident energy, relation to lab_Eout )
void gamma_hit_list::test_lab_Eout( double lab_Eout,
  vector< Ein_Eta_Hit > *Ein_hits )
{
  // the intersections
  Ein_Eta_Hit Ein_hit;

  if( eta < 0.0 )
  {
    // backward emission, actually, eta = mu_cm = -1
    if( right_Ein_Egamma.x <= map->turnover_T )
    {
      // The red-shifted gamma energy is decreasing
      if( ( left_Ein_Egamma.y > lab_Eout ) && ( right_Ein_Egamma.y < lab_Eout ) )
      {
        Ein_hit.E_in = find_decreasing_hit( lab_Eout, left_Ein_Egamma, right_Ein_Egamma );
	Ein_hit.hit_edge = BOTTOM_OUT;
	Ein_hits->push_back( Ein_hit );
      }
    }
    else if( left_Ein_Egamma.x > map->turnover_T )
    {
      // The red-shifted gamma energy is increasing
      if( ( left_Ein_Egamma.y < lab_Eout ) && ( right_Ein_Egamma.y > lab_Eout ) )
      {
        Ein_hit.E_in = find_increasing_hit( lab_Eout, left_Ein_Egamma, right_Ein_Egamma );
	Ein_hit.hit_edge = BOTTOM_IN;
	Ein_hits->push_back( Ein_hit );
      }
    }
    else
    {
      // The red-shifted gamma energy has an interior minimum
      // for the root finders find_decreasing_hit and find_increasing_hit
      dd_entry mid_Ein;
      mid_Ein.x = map->turnover_T;
      mid_Ein.y = map->min_E_back;
      // on the decreasing interval
      if( ( left_Ein_Egamma.y > lab_Eout ) && ( mid_Ein.y < lab_Eout ) )
      {
        Ein_hit.E_in = find_decreasing_hit( lab_Eout, left_Ein_Egamma, mid_Ein );
	Ein_hit.hit_edge = BOTTOM_OUT;
	Ein_hits->push_back( Ein_hit );
      }
      // on the increasing interval
      if( ( mid_Ein.y < lab_Eout ) && ( right_Ein_Egamma.y > lab_Eout ) )
      {
        Ein_hit.E_in = find_increasing_hit( lab_Eout, mid_Ein, right_Ein_Egamma );
	Ein_hit.hit_edge = BOTTOM_IN;
	Ein_hits->push_back( Ein_hit );
      }
    }
  }
  else
  {
    // forward emission, eta = mu_cm = 1
    if( ( left_Ein_Egamma.y < lab_Eout ) && ( right_Ein_Egamma.y > lab_Eout ) )
    {
      Ein_hit.E_in = find_increasing_hit( lab_Eout, left_Ein_Egamma, right_Ein_Egamma );
      Ein_hit.hit_edge = BOTTOM_IN;
      Ein_hits->push_back( Ein_hit );
    }
  }
}
// ------------------ gamma_hit_list::set_E_range --------------
//! Sets left_Ein_Egamma and right_Ein_Egamma, the pairs ( incident energy, lab gamma energy )
//! \param mu_cm direction cosine of gamma in the center-of-mass frame
//! \param E_in_left lower end of the incident energy range
//! \param E_in_right upper end of the incident energy range
void gamma_hit_list::set_E_range( double mu_cm, double E_in_left, double E_in_right )
{
  eta = mu_cm;
  // set up the void* parameters for the functions
  void *params = static_cast<void *>( map );

  left_Ein_Egamma.x = E_in_left;
  right_Ein_Egamma.x = E_in_right;

  if( mu_cm < 0.0 )
  {
    // backward emission, actually, mu_cm = -1
    left_Ein_Egamma.y = Capture_Gamma_F::red_shift( left_Ein_Egamma.x, params );
    right_Ein_Egamma.y = Capture_Gamma_F::red_shift( right_Ein_Egamma.x, params );
  }
  else
  {
    // forward emission, mu_cm = 1
    left_Ein_Egamma.y = Capture_Gamma_F::blue_shift( left_Ein_Egamma.x, params );
    right_Ein_Egamma.y = Capture_Gamma_F::blue_shift( right_Ein_Egamma.x, params );
  }
}
// ------------------ gamma_hit_list::find_bottom_hits --------------
//! Finds the intersections with the bottom of a box
//! \param lab_E_out the lower desired outgoing energy
//! \param Ein_hits output: relations between the gamma energy and lab_E_out
void gamma_hit_list::find_bottom_hits( double lab_E_out,
  vector< Ein_Eta_Hit > *Ein_hits )
{
  // Find the intersections
  test_lab_Eout( lab_E_out, Ein_hits );
}
// ------------------ gamma_hit_list::find_top_hits --------------
//! Finds the intersections with the top of a box
//! \param lab_E_out the higher desired outgoing energy
//! \param Ein_hits output: relations between the gamma energy and lab_E_out
void gamma_hit_list::find_top_hits( double lab_E_out,
  vector< Ein_Eta_Hit > *Ein_hits )
{
  // Find the intersections
  test_lab_Eout( lab_E_out, Ein_hits );
  // these are intersections with the top, not the bottom of the energy bin
  for( vector< Ein_Eta_Hit >::iterator this_hit = Ein_hits->begin( );
       this_hit != Ein_hits->end( ); ++this_hit )
  {
    if( this_hit->hit_edge == BOTTOM_OUT )
    {
      this_hit->hit_edge = TOP_IN;
    }
    else if( this_hit->hit_edge == BOTTOM_IN )
    {
      this_hit->hit_edge = TOP_OUT;
    }
  }
}
// ------------------ gamma_hit_list::test_sides --------------
//! The code appends the previously calculated left_hit and right_hit
//! \param Eout_ptr desired lower energy of the outgoing particle
//! \param E_in_left lower incident energy
//! \param E_in_right higher incident energy
void gamma_hit_list::test_sides( vector< double >::const_iterator Eout_ptr,
   double E_in_left, double E_in_right )
{
  Ein_Eta_Hit this_hit;  // to add to the list
  vector< double >::const_iterator next_Eout = Eout_ptr;
  ++next_Eout;
  this_hit.E_in = left_Ein_Egamma.x;

  // are we below or above the box?
  if( ( left_Ein_Egamma.y < *Eout_ptr ) || ( right_Ein_Egamma.y < *Eout_ptr ) )
  {
    this_hit.hit_edge = BELOW;
  }
  else if( ( left_Ein_Egamma.y > *next_Eout ) || ( right_Ein_Egamma.y > *next_Eout ) )
  {
    this_hit.hit_edge = ABOVE;
  }
  else
  {
    this_hit.hit_edge = INSIDE;
  }
  push_back( this_hit );
  this_hit.E_in = right_Ein_Egamma.x;
  push_back( this_hit );
}

// ************* capture_gamma *************
// ----------- capture_gamma::get_Ein_range --------------
//  Gets the range of nontrivial incident energy bins; computes first_Ein and last_Ein
// returns true if the threshold is too high for the energy bins
bool capture_gamma::get_Ein_range( const dd_vector& sigma, const dd_vector& mult,
    const dd_vector& weight,
    const Flux_List& e_flux, const Energy_groups& Ein_groups )
{
  double E_last;

  capture_gamma_Ein_param initial_param;
  bool done = initial_param.get_Ein_range( sigma, mult, weight, e_flux,
                                         Ein_groups, &E_first, &E_last );
  if( done ) return true;

  // check the range of incident energies for the probability data
  capture_gamma::const_iterator data_ptr = begin( );
  double E_data = data_ptr->get_E_in( );
  if( E_data > E_first )
  {
    E_first = E_data;
  }
  first_Ein = Ein_groups.first_bin_ID( E_first );

  data_ptr = end( );
  --data_ptr;
  E_data = data_ptr->get_E_in( );
  if( E_data < E_last )
  {
    E_last = E_data;
  }
  last_Ein = Ein_groups.last_bin_ID( E_last );

  return false;
}
// ---------------- capture_gamma::min_E_back ------------------
void capture_gamma::get_min_E_back( const T_matrix& transfer )
// Finds the minimum energy of backward gammas
{
    dd_entry below;
    dd_entry above;
    below.x = map.threshold;
    vector< double >::const_iterator Ein_ptr = transfer.in_groups.end( );
    --Ein_ptr;
    above.x = *Ein_ptr;
    map.get_turnover( below, above );
}
// ----------- capture_gamma::read_data --------------
void capture_gamma::read_data( data_parser &input_file, int num_Ein )
{
  Ein_interp = interp_flag_F::read_1d_interpolation( input_file );

  Legendre_coefs new_angle_dist;  // angular distribution for one E_in
  capture_gamma::iterator new_angle_ptr;

  // read the data
  for( int Ein_count = 0; Ein_count < num_Ein; ++Ein_count )
  {
    // insert a new angular distribution
    new_angle_ptr = insert( end( ), Legendre_coefs( ) );
    // get the incident energy and the data pairs
    new_angle_ptr->set_E_in( input_file.get_next_double( ) );
    int Order = input_file.get_next_int( ) - 1;
    new_angle_ptr->initialize( Order );
    for( int coef_count = 0; coef_count <= Order; ++coef_count )
    {
      (*new_angle_ptr)[ coef_count ] = input_file.get_next_double( );
    }
  }
  new_angle_ptr->truncate_zeros( );
  //  print( );
}
// ----------- capture_gamma::get_T --------------
// Calculates the transfer matrix for this particle.
// sigma is the cross section.
void capture_gamma::get_T( const dd_vector& sigma, const dd_vector& weight,
  T_matrix& transfer )
{
  if( ( Ein_interp != LINLIN ) && ( Ein_interp != LINLOG ) )
  {
    FatalError( "capture_gamma::get_T",
      "Incident energy interpolation not implemented" );
  }
  // the multiplicity is one
  dd_vector multiple;
  multiple.make_flat( sigma, 1.0 );
  bool done = get_Ein_range( sigma, multiple, weight, transfer.e_flux,
    transfer.in_groups );
  if( done )
  {
    transfer.zero_transfer( );
  }
  // copy the masses for the map to the lab frame
  map.setup_params( &particle_info );
  // the computed threshold may be above 20 MeV
  double our_threshold = map.threshold;
  dd_vector::const_iterator last_sigma = sigma.end( );
  --last_sigma;
  if( last_sigma->x < our_threshold )
  {
    Info( "capture_gamma::get_T", "computed threshold outside the data range" );
    transfer.zero_transfer( );
  }
  // Find the minimum energy of backward gammas
  get_min_E_back( transfer );

  long int quad_count = 0;  // number of 2-d quadratures
  long int Ein_F_count= 0;  // number of calls to Capture_Gamma_F::Ein_F
  long int mu_F_count = 0;  // number of calls to Capture_Gamma_F::mu_F
  int num_negative = 0;  // number of negative sums of Legendre series

  // do the integrals incident bin by incident bin
#pragma omp parallel for schedule( dynamic, 1 ) default( none ) \
  shared( sigma, multiple, weight, transfer ) \
  reduction( +: quad_count ) reduction( +: Ein_F_count ) reduction( +: mu_F_count ) \
  reduction( +: num_negative )
  // now do the integrals incident bin by incident bin
  for( int Ein_bin = first_Ein; Ein_bin < last_Ein; ++Ein_bin )
  {
    capture_gamma_Ein_param Ein_param;
    Ein_param.set_map( map );
    Ein_param.Ein_interp = Ein_interp;
    Ein_param.mu_quad_method = transfer.mu_quad_method;
    // set up the data range for this bin
    Ein_param.setup_bin( Ein_bin, sigma, multiple, weight, transfer.e_flux,
                         transfer.in_groups );
    setup_data( &Ein_param );
    // work on this bin
    for( ; ; )
    {
      // get the incident energy interval common to all data
      set_Ein_range( &Ein_param );
      Eout_ladder( transfer, &Ein_param );
      // go to the next interval
      bool Done = next_ladder( Ein_param.data_E_1, &Ein_param );
      if( Done )
      {
        break;
      }
    }
    quad_count += Ein_param.quad_count;
    Ein_F_count += Ein_param.Ein_F_count;
    mu_F_count += Ein_param.mu_F_count;
    num_negative += Ein_param.num_negative;
  }  // End of parallel loop

  // print the counts of function evaluations
  cout << "2d quadratures: " << quad_count << endl;
  cout << "Capture_Gamma_F::Ein_F calls: " << Ein_F_count << endl;
  cout << "Capture_Gamma_F::mu_F calls: " << mu_F_count << endl;
  cout << "average Capture_Gamma_F::Ein_F calls: " << 1.0*Ein_F_count/quad_count << endl;
  cout << "average Capture_Gamma_F::mu_F calls: " << 1.0*mu_F_count/Ein_F_count << endl;

  if( num_negative > 0 )
  {
    Info( "capture_gamma::get_T", pastenum( "got ", num_negative ) +
          " negative Legendre sums." );
  }
}
// ----------- capture_gamma::setup_data --------------
// Initializes the quadrature parameters; returns true if the threshold is too high
void capture_gamma::setup_data( capture_gamma_Ein_param *Ein_param )
{
  static double skip_tol = Global.Value( "abs_tol" );

  // for center-of-mass data
  Ein_param->set_map( map );

  Ein_param->left_data = begin( );
  Ein_param->right_data = Ein_param->left_data;
  ++Ein_param->right_data;

  while( ( Ein_param->right_data->get_E_in( ) < E_first * ( 1.0 + skip_tol ) ) ||
	 ( Ein_param->right_data->get_E_in( ) < (*Ein_param->Ein_ptr) *
           ( 1.0 + skip_tol ) ) )
  {
    Ein_param->left_data = Ein_param->right_data;
    ++Ein_param->right_data;
  }

  double first_E = Ein_param->left_data->get_E_in( );
  if( first_E > Ein_param->data_E_0 )
  {
    // this should never happen
    Ein_param->data_E_0 = first_E;
    bool data_bad = Ein_param->update_pointers( first_E );
    if( data_bad )
    {
      FatalError( "capture_gamma::setup_data", "energies inconsistent" );
    }
  }
}
// ----------- capture_gamma::next_ladder --------------
bool capture_gamma::next_ladder( double E_in, capture_gamma_Ein_param *Ein_param )
{
  bool done = Ein_param->update_bin_pointers( E_in );
  static double etol = Global.Value( "E_tol" );
  if( !done )
  {
    double E_tol = E_in * etol;
    if( E_in + E_tol >= Ein_param->right_data->get_E_in( ) )
    {
      while( E_in + E_tol >= Ein_param->right_data->get_E_in( ) )
      {
        // get the next angular data
        Ein_param->left_data = Ein_param->right_data;
        ++Ein_param->right_data;
        if( Ein_param->right_data == end ( ) )
        {
          return true;
        }
      }
    }
  }
  return done;
}
// ----------- capture_gamma::set_Ein_range --------------
// Sets the range of incident energies for this intergration
void capture_gamma::set_Ein_range( capture_gamma_Ein_param *Ein_param )
{
  Ein_param->set_Ein_range( );
  double this_E = Ein_param->left_data->get_E_in( );
  if( this_E > Ein_param->data_E_0 ) Ein_param->data_E_0 = this_E;
  this_E = Ein_param->right_data->get_E_in( );
  if( this_E < Ein_param->data_E_1 ) Ein_param->data_E_1 = this_E;

  // the data may be below the threshold
  double our_threshold = map.threshold;
  if( Ein_param->data_E_1 < our_threshold )
  {
    Ein_param->data_E_0 = our_threshold;
    Ein_param->data_E_1 = our_threshold;
    // Warning( "capture_gamma::set_Ein_range", "ignoring 2 data below threshold" );
  }
  else if( Ein_param->data_E_0 < our_threshold )
  {
    Ein_param->data_E_0 = our_threshold;
    // Warning( "capture_gamma::set_Ein_range", "ignoring data below threshold" );
  }

  if( Ein_param->data_E_1 < Ein_param->data_E_0 )
  {
    FatalError( "capture_gamma::set_Ein_range", "check the incident energies" );
  }
  Ein_param->set_sigma_range( );
}
// ----------- capture_gamma::Eout_ladder --------------
// This routine uses the this angular Legendre expansion and the
// next to calculate the contribution to the E_out boxes of the
// transfer matrix between incident energies Ein_param->data_E_0 and
// Ein_param->data_E_1.
void capture_gamma::Eout_ladder( T_matrix& transfer,
  capture_gamma_Ein_param *Ein_param )
{
  bool geom_OK;  // for checking the consistency of the geometry
  gamma_hit_list test_hits;
  test_hits.map = &map;

    // loop through the outgoing energies (column of transfer)
    for( int Eout_count = 0; Eout_count < transfer.num_Eout_bins;
      ++Eout_count )
    {
      vector< double >::const_iterator Eout_ptr = transfer.out_groups.begin( )
        + Eout_count;
      // how does the mu = -1 curve meet this E-E' box?
      Ein_param->backward_hits.set_E_range( -1.0, Ein_param->data_E_0, Ein_param->data_E_1 );
      geom_OK = Ein_param->backward_hits.hit_box( -1.0, Eout_ptr,
        Ein_param->data_E_0, Ein_param->data_E_1 );
      if( !geom_OK )
      {
        test_hits.set_E_range( -1.0, Ein_param->data_E_0, Ein_param->data_E_1 );
        test_hits.hit_box( -1.0, Eout_ptr,
          Ein_param->data_E_0, Ein_param->data_E_1 );
        test_hits.print( );
        FatalError( "capture_gamma::Eout_ladder", "Check the coding, 1" );
      }
      if( ( Eout_count < transfer.num_Eout_bins - 1 ) &&
          ( Ein_param->backward_hits.is_above( ) ) )
      {
        // go on to the next E-E' box
        continue;
      }
      // how does the mu = 1 curve meet this E-E' box?
      Ein_param->forward_hits.set_E_range( 1.0, Ein_param->data_E_0, Ein_param->data_E_1 );
      geom_OK = Ein_param->forward_hits.hit_box( 1.0, Eout_ptr,
        Ein_param->data_E_0, Ein_param->data_E_1 );
      if( !geom_OK )
      {
        test_hits.set_E_range( 1.0, Ein_param->data_E_0, Ein_param->data_E_1 );
        test_hits.hit_box( 1.0, Eout_ptr,
                           Ein_param->data_E_0, Ein_param->data_E_1 );
        test_hits.print( );
        FatalError( "capture_gamma::Eout_ladder", "Check the coding, 2" );
      }
      if( ( Eout_count > 0 ) && ( Ein_param->forward_hits.is_below( ) ) )
      {
        // we are done with this pair of eta values
        break;
      }
      // integrate over this E-E' box
      one_Ebox( transfer, Eout_count, Ein_param );
    }
}
// ----------- capture_gamma::one_Ebox --------------
// Integrate over one E-E' box
void capture_gamma::one_Ebox( T_matrix& transfer, int Eout_count,
   capture_gamma_Ein_param *Ein_param )
{
  // the E' energy range
  Ein_param->Eout_min = transfer.out_groups[ Eout_count ];
  Ein_param->Eout_max = transfer.out_groups[ Eout_count + 1 ];

  // set up common incident energies
  Ein_param->backward_hits.common_hits( Ein_param->forward_hits );

  // integrate depending on how the curves mu_cm = const meet the box
  gamma_hit_list::iterator low_hit_ptr = Ein_param->backward_hits.begin( );
  gamma_hit_list::iterator next_low_ptr = low_hit_ptr;
  ++next_low_ptr;
  gamma_hit_list::iterator high_hit_ptr = Ein_param->forward_hits.begin( );
  gamma_hit_list::iterator next_high_ptr = high_hit_ptr;
  ++next_high_ptr;
  for( ; ( next_low_ptr != Ein_param->backward_hits.end( ) ) &&
         ( next_high_ptr != Ein_param->forward_hits.end( ) );
       low_hit_ptr = next_low_ptr, ++next_low_ptr,
         high_hit_ptr = next_high_ptr, ++next_high_ptr )
  {
    if( ( low_hit_ptr->hit_edge == ABOVE ) ||
        ( low_hit_ptr->hit_edge == TOP_OUT ) )
    {
      // do nothing---we are above the E-E' box
      continue;
    }
    else if( ( low_hit_ptr->hit_edge == TOP_IN ) ||
             ( low_hit_ptr->hit_edge == BOTTOM_IN ) ||
             ( low_hit_ptr->hit_edge == INSIDE ) )
    {
      // the backward mu_cm = const curve is inside the E-E' box
      Ein_param->use_Eout_min = false;
      // where is the forward curve?
      if( ( high_hit_ptr->hit_edge == ABOVE ) ||
          ( high_hit_ptr->hit_edge == TOP_OUT ) )
      {
        // integrate up to the top of the E-E' bin
        Ein_param->use_Eout_max = true;
      }
      else
      {
        // integrate up to the next mu_cm = const curve
        Ein_param->use_Eout_max = false;
      }
    }
    else
    {
      // the backward mu_cm = const curve is below the E-E' box;
      // integrate from Eout_min
      Ein_param->use_Eout_min = true;
      // where is the forward mu_cm = const curve?
      if( ( high_hit_ptr->hit_edge == BOTTOM_OUT ) ||
          ( high_hit_ptr->hit_edge == BELOW ) )
      {
        // do nothing---we are below the E-E' box
        continue;
      }
      else if( ( high_hit_ptr->hit_edge == TOP_IN ) ||
               ( high_hit_ptr->hit_edge == BOTTOM_IN ) ||
               ( high_hit_ptr->hit_edge == INSIDE ) )
      {
        // the forward mu_cm = const curve is inside the E-E' box
        Ein_param->use_Eout_max = false;
      }
      else
      {
        // the forward mu_cm = const curve is above the E-E' box
        Ein_param->use_Eout_max = true;
      }
    }
    // the range of integration in incident energy
    Ein_param->Ein_0 = low_hit_ptr->E_in;
    Ein_param->Ein_1 = next_low_ptr->E_in;
    update_T( transfer, Eout_count, Ein_param );
  }
}
// ----------- capture_gamma::update_T --------------
void capture_gamma::update_T( T_matrix &transfer, int Eout_count,
   capture_gamma_Ein_param *Ein_param )
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
    double left_E = ( Ein_param->this_sigma->x < Ein_param->Ein_0 ) ?
      Ein_param->Ein_0 : Ein_param->this_sigma->x;
    double right_E = ( Ein_param->next_sigma->x > Ein_param->Ein_1 ) ?
      Ein_param->Ein_1 : Ein_param->next_sigma->x;
    // evaluate the integral
    quad_F::integrate( Capture_Gamma_F::Ein_F, transfer.Ein_quad_method, left_E,
		       right_E, params, tol, &value );

    // add this integral
    transfer( Ein_param->Ein_count, Eout_count ) += value;
    // increment the function counts
    Ein_param->Ein_F_count += Ein_param->func_count;
    ++Ein_param->quad_count;
  }
}


// ************* functions *************
// ------------------ Capture_Gamma_F::red_shift -----------------------
// Returns the energy in lab frame of the backward emitted gamma.
double Capture_Gamma_F::red_shift( double T_in_lab, void *params )
{
  // the parameters are really red_blue_shift_map *
  red_blue_shift_map *shift_params = static_cast<red_blue_shift_map *>( params );
  // Sets up the boost to the lab frame
  shift_params->set_boost( T_in_lab );
  // red shift
  double factor = shift_params->cosh_chi - shift_params->sinh_chi;
  return factor*shift_params->T_cm_out;
}

// ------------------ Capture_Gamma_F::d_red_shift -----------------------
// Returns the derivative of the energy of the backward emitted gamma.
double Capture_Gamma_F::d_red_shift( double T_in_lab, void *params )
{
  // the parameters are really red_blue_shift_map *
  red_blue_shift_map *shift_params = static_cast<red_blue_shift_map *>( params );
  // Sets up the boost to the lab frame
  shift_params->set_boost( T_in_lab );
  // \partial Egamma_cm/ \partial T_in_lab
  double dEgamma_cm = 0.5*( 1.0 + 
      shift_params->rest_masses->mRes*shift_params->rest_masses->mRes /
			    ( shift_params->Minkowski*shift_params->Minkowski ) )*
    shift_params->rest_masses->mTarg/shift_params->Minkowski;

  // \partial chi/ \partial T_in_lab
  double dchi = ( ( T_in_lab + shift_params->rest_masses->mProj )/shift_params->sinh_chi -
		  shift_params->sinh_chi* shift_params->rest_masses->mTarg )/
    ( shift_params->cosh_chi*shift_params->Minkowski*shift_params->Minkowski );

  // red shift
  double factor = shift_params->cosh_chi - shift_params->sinh_chi;

  return factor*( dEgamma_cm - dchi*shift_params->T_cm_out );
}

// ------------------ Capture_Gamma_F::blue_shift -----------------------
// Returns the energy in lab frame of the forward emitted gamma.
double Capture_Gamma_F::blue_shift( double T_in_lab, void *params )
{
  // the parameters are really red_blue_shift_map *
  red_blue_shift_map *shift_params = static_cast<red_blue_shift_map *>( params );
  // Sets up the boost to the lab frame
  shift_params->set_boost( T_in_lab );
  // blue shift
  double factor = shift_params->cosh_chi + shift_params->sinh_chi;
  return factor*shift_params->T_cm_out;
}

// Function for the 1-d quadrature
// ---------------- Capture_Gamma_F::mu_F ------------------
void Capture_Gamma_F::mu_F( double mu_cm, QuadParamBase *mu_cm_quad_param,
  coef_vector *value )
{
  // the parameters are really capture_gamma_mu_param
  capture_gamma_mu_param *params = static_cast<capture_gamma_mu_param*>( mu_cm_quad_param );
  params->func_count += 1;
  //  if( params->func_count % 100 == 0 )
  //  {
  //    Info( "Capture_Gamma_F::mu_F", pastenum( "got ",
  //       params->func_count ) + " evaluations");
  //  }
  // get Eout_lab and mu_lab
  double Eout_lab;
  double mu_lab;
  params->map->get_E_mu_lab( mu_cm, &Eout_lab, &mu_lab );

  // the Legendre polynomials
  math_F::Legendre( mu_lab, value );

  // the probability density
  double Prob = params->coefs.sum_Legendre( mu_cm );
  *value *= Prob;
  if( Prob < 0.0 )
  {
    ++(params->num_negative);
    if( ( !params->flag_set ) && ( params->num_negative == 1 ) )
    {
      double E_in = params->coefs.get_E_in( );
      string Ein_value = pastenum( "Negative Legendre sum for E_in:", E_in );
      string mu_value = pastenum( " and mu_cm:", mu_cm );
      Info( "Capture_Gamma_F::mu_F", Ein_value + mu_value );
    }
  
  }
  // do the energy weighting if necessary
  if( ( value->conserve == ENERGY ) || ( value->conserve == BOTH ) )
  {
    value->scale_E( Eout_lab );
  }
}

// ---------------- Capture_Gamma_F::Ein_F ------------------
void Capture_Gamma_F::Ein_F( double E_in, QuadParamBase *e_quad_param,
  coef_vector *value )
// Function for the 2-d quadrature
{
  // the parameters are really capture_gamma_Ein_param *
  capture_gamma_Ein_param *e_params = static_cast<capture_gamma_Ein_param *>( e_quad_param );
  e_params->func_count += 1;
  //  if( e_params->func_count % 100 == 0 )
  //  {
  //    Info( "Capture_Gamma_F::Ein_F", pastenum( "got ", e_params->func_count ) + " evaluations");
  //  }

  // The value of Capture_Gamma_F::Ein_F is itself an integral over mu_cm.
  // *value comes in as 0.  

  // parameters for the integration over mu_cm
  capture_gamma_mu_param mu_cm_params;
  mu_cm_params.flag_set = ( e_params->num_negative > 0 );
  mu_cm_params.coefs.set_E_in( E_in );
  mu_cm_params.map = &e_params->map;
  mu_cm_params.map->set_boost( E_in );
  mu_cm_params.Ein_interp = e_params->Ein_interp;
  // interpolate the (mu_cm, probability) with respect to incident energy
  int max_order = ( e_params->left_data->order > e_params->right_data->order ) ?
    e_params->left_data->order : e_params->right_data->order;
  mu_cm_params.coefs.initialize( max_order );
  mu_cm_params.interpolate( E_in, e_params->left_data, e_params->right_data );

  // the range of integration
  double mu_cm_0 = ( e_params->use_Eout_min ) ?
    mu_cm_params.map->get_mu_cm( e_params->Eout_min ) : -1.0;
  double mu_cm_1 = ( e_params->use_Eout_max ) ?
    mu_cm_params.map->get_mu_cm( e_params->Eout_max ) : 1.0;

  // evaluate the integral over mu_cm
  QuadParamBase *params = static_cast< QuadParamBase* >( &mu_cm_params );
  static double tol = Global.Value( "quad_tol" );
  quad_F::integrate( Capture_Gamma_F::mu_F, e_params->mu_quad_method,
		      mu_cm_0, mu_cm_1, params, tol, value );
  e_params->num_negative += mu_cm_params.num_negative;
  e_params->mu_F_count +=  mu_cm_params.func_count;
  // weight it by flux * cross section
  e_params->set_weight( E_in );
  *value *= e_params->current_weight;
  //  cout << "E_in: " << E_in << " mu_cm_0: " << mu_cm_0 << " mu_cm_1: " <<
  //    mu_cm_1 << endl;
  //  value->print( );
}
