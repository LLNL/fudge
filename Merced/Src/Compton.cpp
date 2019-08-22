/*
 * ******** merced: calculate the transfer matrix ********* *
 * $Revision: 1 $
 * $Date: 2008-04-16 19:06:56 -0800 (Wed, 01Feb 2006) $
 * $Author: hedstrom $
 * $Id: Compton.cpp 1 2008-04-16 19:06:56 -0800Z hedstrom $
 * initial code
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// implementation of the classes used for Compton scattering

#include <cmath>
#ifdef _OPENMP
 #include <omp.h>
#endif

#include "Compton.hpp"
#include "math_util.hpp"
#include "messaging.hpp"
#include "global_params.hpp"

// **************** class Compton_Ein_param ******************

// **************** class Eout_curve ******************
// ------------------ Eout_curve::setup --------------
// Calculates the incident energy for backscatter, mu = -1
void Eout_curve::setup( double Eout )
{
  E_out = Eout;
  static double m_e = Global.Value( "m_electron" );
  if( E_out < m_e/2 )
  {
    backscatter = true;
    Ein_back = E_out / ( 1.0 - ( 2*E_out / m_e ) );
  }
  else
  {
    backscatter = false;
  }
}
// ------------------ Eout_curve::hit_x --------------
// Finds (Ein, mu) where this Eout curve hits the curve x = constant.
dd_entry Eout_curve::hit_x( double x ) const
{
  static double m_e = Global.Value( "m_electron" );
  dd_entry intersection;
  if( x == 0.0 )
  {
    intersection.y = 1.0;
    intersection.x = E_out;
  }
  else
  {
    // Solve the quadratic for xi = sqrt(1-mu)
    //    1/E_in = xi/(x * sqrt(2) ),
    //    1/E_in = 1/E_out - xi^2/m_e
    double c = m_e/E_out;
    double b = m_e/(x * sqrt(2.0) );
    double xi = ( 2.0 * c )/( b + sqrt( b*b + 4.0*c ) );
    double mu = 1.0 - xi*xi;
    if( mu <= -1.0)
    {
      // the curves don't intersect inside the region -1 <= mu <= 1
      intersection.y = -1.0;
      intersection.x = x;
    }
    else
    {
      intersection.y = mu;
      intersection.x = x*sqrt( 2.0 ) / xi;
    }
  }
  return intersection;
}
// ------------------ Eout_curve::hit_Ein --------------
// Finds mu where this Eout curve hits the line Ein = constant.
double Eout_curve::hit_Ein( double E_in ) const
{
  static double m_e = Global.Value( "m_electron" );
  if( ( E_in == 0.0 ) || ( E_out == 0.0 ) )
  {
    FatalError( "Eout_curve::hit_Ein", "zero energy" );
  }
  else if( E_in < E_out )
  {
    FatalError( "Eout_curve::hit_Ein", "outgoing energy too big" );
  }
  // Solve the equation for mu:
  //  1/E_out - 1/E_in = (1 - mu)/m_e.
  double mu = 1.0 - m_e*( 1.0/E_out - 1.0/E_in );
  return mu;
}

// **************** class Boundary_corner ******************
// ---------------- Boundary_corner constructor ------------------
Boundary_corner::Boundary_corner( )
{
}
// ------------------ Boundary_corner copy constructor --------------
Boundary_corner::Boundary_corner( const Boundary_corner& boundary_corner )
{
  E_in = boundary_corner.E_in;
  boundary_ID = boundary_corner.boundary_ID;
}

// **************** class corner_list ******************
// ---------------- corner_list::new_corner ------------------
// Append a corner to the list
void corner_list::new_corner( double x, const Eout_curve &Eoutcurve, Boundary_ID boundary_ID )
{
  Boundary_corner boundary_corner;  // for new entry
  boundary_corner.boundary_ID = boundary_ID;
  dd_entry corner_E_mu = Eoutcurve.hit_x( x );  // ( Ein, mu ) for a corner
  boundary_corner.E_in = corner_E_mu.x;
  push_back( boundary_corner );
}
// ---------------- corner_list::new_corner ------------------
// Append a corner to the list
void corner_list::new_corner( double E, Boundary_ID boundary_ID )
{
  Boundary_corner boundary_corner;  // for new entry
  boundary_corner.boundary_ID = boundary_ID;
  boundary_corner.E_in = E;
  push_back( boundary_corner );
}

// ---------------- corner_list::copy_last ------------------
// Copy the last entry from copy_from but change the boundary ID
void corner_list::copy_last( const corner_list &copy_from, Boundary_ID boundary_ID )
{
  Boundary_corner boundary_corner;  // for new entry
  corner_list::const_iterator to_copy = copy_from.end( );
  --to_copy;
  boundary_corner.E_in = to_copy->E_in;
  boundary_corner.boundary_ID = boundary_ID;
  push_back( boundary_corner );
}

// **************** class Compton_region ******************
// ---------------- Compton_region::setup_Eout ------------------
// constructs upper_Eout and lower_Eout
void Compton_region::setup_Eout( vector< double >::const_iterator Eout_ptr )
{
  // reset lower_Eout
  lower_Eout.setup( *Eout_ptr );

  // reset upper_Eout
  vector< double >::const_iterator next_Eout = Eout_ptr;
  ++next_Eout;
  upper_Eout.setup( *next_Eout );
}
// ---------------- Compton_region::setup_region ------------------
// Identifies the top and bottom of the quadrature region
// The incident energy for the x-data ranges over (x_0, x_1).
// The outgoing energy ranges over (lower_Eout.E_out, upper_Eout.E_out),
// but E_out = E_in for mu = 1.
void Compton_region::setup_region( double x_0, double x_1 )
{
  x0 = x_0;
  x1 = x_1;
  // reset bottom_corners
  if( bottom_corners.size( ) > 0 )
  {
    bottom_corners.erase( bottom_corners.begin( ), bottom_corners.end( ) );
  }
  // reset top_corners
  if( top_corners.size( ) > 0 )
  {
    top_corners.erase( top_corners.begin( ), top_corners.end( ) );
  }

  Boundary_corner boundary_corner;  // for new entry
  // test for void intersection
  if( upper_Eout.backscatter && x0 >= upper_Eout.Ein_back )
  {
    FatalError( "Compton_region::setup_Eout", "intersection void" );
    //    return;
  }
  // special coding for x_0 == 0
  if( x_0 == 0 )
  {
    top_corners.new_corner( lower_Eout.E_out, USE_MU1 );
    top_corners.new_corner( upper_Eout.E_out, USE_EOUT );
    if( lower_Eout.E_out == 0.0 )
    {
      bottom_corners.new_corner( 0.0, USE_MU_MINUS1 );
      if( ( upper_Eout.backscatter && x1 < upper_Eout.Ein_back )  ||
          !upper_Eout.backscatter )
      {
        bottom_corners.new_corner( x1, USE_X );
        bottom_corners.new_corner( x1, upper_Eout, USE_X );
      }
      else
      {
        bottom_corners.new_corner( upper_Eout.Ein_back, USE_X );
      }
    }
    else
    {
      bottom_corners.new_corner( lower_Eout.E_out, USE_EOUT );
      if( ( lower_Eout.backscatter && x1 < lower_Eout.Ein_back )  ||
          !lower_Eout.backscatter )
      {
        bottom_corners.new_corner( x1, lower_Eout, USE_X );
        bottom_corners.new_corner( x1, upper_Eout, USE_X );
      }
      else if( ( upper_Eout.backscatter && x1 < upper_Eout.Ein_back )  ||
                !upper_Eout.backscatter )
      {
        bottom_corners.new_corner( lower_Eout.Ein_back, USE_MU_MINUS1 );
        bottom_corners.new_corner( x1, USE_X );
        bottom_corners.new_corner( x1, upper_Eout, USE_X );
      }
      else
      {
        bottom_corners.new_corner( lower_Eout.Ein_back, USE_MU_MINUS1 );
        bottom_corners.new_corner( upper_Eout.Ein_back, USE_MU_MINUS1 );
      }
    }
    top_corners.copy_last( bottom_corners, USE_X );
  }
  // the most common case is for -1 < mu < 1
  else if( ( lower_Eout.backscatter && x1 < lower_Eout.Ein_back )  ||
           !lower_Eout.backscatter )
  {
    bottom_corners.new_corner( x0, lower_Eout, USE_EOUT );
    top_corners.copy_last( bottom_corners, USE_X );
    top_corners.new_corner( x0, upper_Eout, USE_EOUT );
    bottom_corners.new_corner( x1, lower_Eout, USE_X );
    bottom_corners.new_corner( x1, upper_Eout, USE_X );
    top_corners.copy_last( bottom_corners, USE_X );
  }
  else if( ( upper_Eout.backscatter && x1 < upper_Eout.Ein_back )  ||
      !upper_Eout.backscatter )
  {
    if( lower_Eout.backscatter && x_0 >= lower_Eout.Ein_back )
    {
      bottom_corners.new_corner( x0, USE_MU_MINUS1 );
      bottom_corners.new_corner( x1, USE_X );
      bottom_corners.new_corner( x1, upper_Eout, USE_X );
      top_corners.new_corner( x0, USE_X );
      top_corners.new_corner( x0, upper_Eout, USE_EOUT );
      top_corners.copy_last( bottom_corners, USE_X );
    }
    else //  0.0 < x_0 < lower_Eout.Ein_back
    {
      bottom_corners.new_corner( x0, lower_Eout, USE_EOUT );
      top_corners.copy_last( bottom_corners, USE_X );
      bottom_corners.new_corner( lower_Eout.Ein_back, USE_MU_MINUS1 );
      bottom_corners.new_corner( x1, USE_X );
      bottom_corners.new_corner( x1, upper_Eout, USE_X );
      top_corners.new_corner( x0, upper_Eout, USE_EOUT );
      top_corners.copy_last( bottom_corners, USE_X );
    }
  }
  else  //  upper_Eout.back_scatter && (x1 >= upper_Eout.Ein_back)
  {
    // We know that lower_Eout.backscatter is true.
    if( x0 >= lower_Eout.Ein_back )
    {
      bottom_corners.new_corner( x0, USE_MU_MINUS1 );
      bottom_corners.new_corner( upper_Eout.Ein_back, USE_MU_MINUS1 );
      top_corners.new_corner( x0, USE_X );
      top_corners.new_corner( x0, upper_Eout, USE_EOUT );
      top_corners.new_corner( upper_Eout.Ein_back, USE_MU_MINUS1 );
    }
    else
    {
      bottom_corners.new_corner( x0, lower_Eout, USE_EOUT );
      top_corners.copy_last( bottom_corners, USE_X );
      bottom_corners.new_corner( lower_Eout.Ein_back, USE_MU_MINUS1 );
      bottom_corners.new_corner( upper_Eout.Ein_back, USE_MU_MINUS1 );
      top_corners.new_corner( x0, upper_Eout, USE_EOUT );
      top_corners.new_corner( upper_Eout.Ein_back, USE_MU_MINUS1 );
    }
  }
  left_corner = bottom_corners.begin( );
  right_corner = bottom_corners.end( );
  --right_corner;
}

// ********* class Compton *********
// ---------------- Compton::get_T ------------------
// Calculates the transfer matrix for this particle.
void Compton::get_T( T_matrix& transfer, dd_vector& sigma )
{
  // ensure that the data starts at x=0
  prepend( file_data, 0.0 );

   // :::::::: kludge ::::::::::::::::
  // scale x from 1/(wave length) to energy
  double x_to_energy = Global.Value( "x_to_energy" );
  scale_x_to_energy( x_to_energy );
  // :::::::: end kludge ::::::::::::::::

  // set up multiplicity and cross section
  dd_vector multiple;
  multiple.make_flat( transfer.in_groups, 1.0 );
  dd_vector weight;
  weight.make_flat( transfer.in_groups, 1.0 );

  // compute the cross section
  get_xsec( transfer.mu_quad_method, sigma );

  vector< double >::const_iterator last_Ein = transfer.in_groups.end( );
  --last_Ein;

  long int quad_count = 0;
  // number of calls to Compton_F::Ein_F
  long int Ein_F_count = 0;
  long int mu_F_count = 0;  // number of calls to Compton_F::mu_F

  // loop over the outgoing energy groups
#pragma omp parallel for schedule( dynamic, 1 ) default( none ) \
  shared( sigma, multiple, weight, transfer, last_Ein )		\
  reduction( +: quad_count ) reduction( +: Ein_F_count ) \
  reduction( +: mu_F_count )
  for( int Eout_count = 0; Eout_count < transfer.num_Eout_bins; ++Eout_count )
  {
    vector< double >::const_iterator Eout_ptr = transfer.out_groups.begin( ) +
      Eout_count;

    Compton_Ein_param Ein_param;
    Ein_param.mu_quad_method = transfer.mu_quad_method;
    Ein_param.quad_box.Eout_count = Eout_count;
    Ein_param.quad_box.setup_Eout( Eout_ptr );
    Ein_param.Eout_ptr = &Ein_param.quad_box.lower_Eout;

    // reset the initial pointers
    Ein_param.setup( sigma, multiple, weight, transfer.e_flux, transfer.in_groups );
    Ein_param.this_sigma = sigma.begin( );
    Ein_param.next_sigma = sigma.begin( );
    ++Ein_param.next_sigma;

    // loop over the x-values
    Ein_param.x_ptr = begin( );
    Ein_param.next_x = Ein_param.x_ptr;
    ++Ein_param.next_x;
    for( ; Ein_param.next_x != end( ); 
         Ein_param.x_ptr = Ein_param.next_x, ++Ein_param.next_x )
    {
      Ein_param.next_Eout = &Ein_param.quad_box.upper_Eout;
      if( Ein_param.x_ptr->x >= *last_Ein )
      {
        // The remaining data doesn't contribute to the transfer matrix
        break;
      }
      if( Ein_param.quad_box.upper_Eout.backscatter &&
	  ( ( Ein_param.quad_box.upper_Eout.Ein_back <= *(transfer.in_groups.begin( ) ) ) ||
	    ( Ein_param.quad_box.upper_Eout.Ein_back <= Ein_param.x_ptr->x ) ) )
      {
        continue; // the outgoing energy is too low
      }
      else
      {
        Ein_param.quad_box.setup_region( Ein_param.x_ptr->x, Ein_param.next_x->x );
        Ein_param.first_ladder_sigma = sigma.begin( );  // initialize again
        Ein_ladder( transfer, &Ein_param );
      }
    }
    quad_count += Ein_param.quad_count;
    Ein_F_count += Ein_param.Ein_F_count;
    mu_F_count += Ein_param.mu_F_count;
  } // end of parallel loop

  // print the counts of function evaluations
  cout << "2d quadratures: " << quad_count << endl;
  cout << "Compton_F::Ein_F calls: " << Ein_F_count << endl;
  cout << "Compton_F::mu_F calls: " << mu_F_count << endl;
  cout << "average Compton_F::Ein_F calls: " << 1.0*Ein_F_count/quad_count << endl;
  cout << "average Compton_F::mu_F calls: " << 1.0*mu_F_count/Ein_F_count << endl;
}
// ---------------- Compton::Ein_ladder ------------------
// Climbs through the incident energy groups
void Compton::Ein_ladder( T_matrix& transfer, Compton_Ein_param *Ein_param )
{
  // loop through the incident energies (row of transfer)
  for( Ein_param->quad_box.Ein_count = 0;
    Ein_param->quad_box.Ein_count < transfer.num_Ein_bins;
    ++Ein_param->quad_box.Ein_count )
  {
    Ein_param->Ein_ptr = transfer.in_groups.begin( )
      + Ein_param->quad_box.Ein_count;
    Ein_param->next_Ein = Ein_param->Ein_ptr;
    ++Ein_param->next_Ein;
    if( Ein_param->quad_box.upper_Eout.backscatter &&
	( Ein_param->quad_box.upper_Eout.Ein_back <= Ein_param->quad_box.left_corner->E_in ) )
    {
      continue;  // incident energy too low---go to the next
    }
    else if( *(Ein_param->next_Ein) <= Ein_param->quad_box.left_corner->E_in )
    {
      continue;  // incident energy too low---go to the next
    }
    else if( *(Ein_param->Ein_ptr) >= Ein_param->quad_box.right_corner->E_in )
    {
      break;  // incident energy too high
    }
    else
    {
      one_Ebox( transfer, Ein_param );
    }
  }
}
// ----------- Compton::one_Ebox --------------
// Integrate over one x-E box; loop over the (E_in, mu) Compton region
void Compton::one_Ebox( T_matrix& transfer, Compton_Ein_param *Ein_param )
{
  // sets the range of integration
  set_Ein_range( Ein_param );
  Ein_param->set_sigma_range( );
  // point to the low end of the Compton_region
  Ein_param->bottom_ptr = Ein_param->quad_box.bottom_corners.begin( );
  Ein_param->next_bottom_ptr = Ein_param->bottom_ptr;
  ++Ein_param->next_bottom_ptr;
  Ein_param->top_ptr = Ein_param->quad_box.top_corners.begin( );
  Ein_param->next_top_ptr = Ein_param->top_ptr;
  ++Ein_param->next_top_ptr;

  // loop over the corners of the (E_in, mu) Compton region
  while( ( Ein_param->top_ptr->E_in < Ein_param->data_E_1 ) &&
	 ( Ein_param->next_top_ptr != Ein_param->quad_box.top_corners.end( ) ) &&
         ( Ein_param->bottom_ptr->E_in < Ein_param->data_E_1 ) &&
	 ( Ein_param->next_bottom_ptr != Ein_param->quad_box.bottom_corners.end( ) ) )
  {
    if( Ein_param->next_top_ptr->E_in <= Ein_param->data_E_0 )
    {
      // top_ptr too low
      Ein_param->top_ptr = Ein_param->next_top_ptr;
      ++Ein_param->next_top_ptr;
    }
    else if( Ein_param->next_bottom_ptr->E_in <= Ein_param->data_E_0 )
    {
      // bottom_ptr too low
      Ein_param->bottom_ptr = Ein_param->next_bottom_ptr;
      ++Ein_param->next_bottom_ptr;
    }
    else
    {
      if( ( Ein_param->bottom_ptr->E_in <= Ein_param->data_E_0 ) &&
	  ( Ein_param->top_ptr->E_in <= Ein_param->data_E_0 ) )
      {
	Ein_param->Ein_0 = Ein_param->data_E_0;
	Ein_param->Ein_1 = 
	  ( Ein_param->next_bottom_ptr->E_in <= Ein_param->next_top_ptr->E_in ) ?
	  Ein_param->next_bottom_ptr->E_in :
	  Ein_param->next_top_ptr->E_in;
	if( Ein_param->Ein_1 > Ein_param->data_E_1 )
	{
	  Ein_param->Ein_1 = Ein_param->data_E_1;
	}
      }
      else if( Ein_param->bottom_ptr->E_in < Ein_param->top_ptr->E_in )
      {
	Ein_param->Ein_0 = Ein_param->top_ptr->E_in;
	Ein_param->Ein_1 = 
	  ( Ein_param->next_bottom_ptr->E_in <= Ein_param->next_top_ptr->E_in ) ?
	  Ein_param->next_bottom_ptr->E_in :
	  Ein_param->next_top_ptr->E_in;
	if( Ein_param->Ein_1 > Ein_param->data_E_1 )
	{
	  Ein_param->Ein_1 = Ein_param->data_E_1;
	}
      }
      else if( Ein_param->bottom_ptr->E_in > Ein_param->top_ptr->E_in )
      {
	Ein_param->Ein_0 = Ein_param->bottom_ptr->E_in;
	Ein_param->Ein_1 = 
	  ( Ein_param->next_top_ptr->E_in <= Ein_param->next_bottom_ptr->E_in ) ?
	  Ein_param->next_top_ptr->E_in :
	  Ein_param->next_bottom_ptr->E_in;
	if( Ein_param->Ein_1 > Ein_param->data_E_1 )
	{
	  Ein_param->Ein_1 = Ein_param->data_E_1;
	}
      }
      else // Ein_param->bottom_ptr->E_in == Ein_param->top_ptr->E_in
      {
	Ein_param->Ein_0 = Ein_param->top_ptr->E_in;
	Ein_param->Ein_1 = 
	  ( Ein_param->next_top_ptr->E_in <= Ein_param->next_bottom_ptr->E_in ) ?
	  Ein_param->next_top_ptr->E_in :
	  Ein_param->next_bottom_ptr->E_in;
	if( Ein_param->Ein_1 > Ein_param->data_E_1 )
	{
	  Ein_param->Ein_1 = Ein_param->data_E_1;
	}
      }
      update_T( transfer, Ein_param );
      // update the pointers to the integration region
      if( Ein_param->next_top_ptr->E_in < Ein_param->next_bottom_ptr->E_in )
      {
        Ein_param->top_ptr = Ein_param->next_top_ptr;
        ++Ein_param->next_top_ptr;
      }
      else if( Ein_param->next_top_ptr->E_in > Ein_param->next_bottom_ptr->E_in )
      {
        Ein_param->bottom_ptr = Ein_param->next_bottom_ptr;
        ++Ein_param->next_bottom_ptr;
      }
      else
      {
        Ein_param->top_ptr = Ein_param->next_top_ptr;
        ++Ein_param->next_top_ptr;
        Ein_param->bottom_ptr = Ein_param->next_bottom_ptr;
        ++Ein_param->next_bottom_ptr;
      }
    }
  }
}
// ---------------- Compton::set_Ein_range ------------------
// Sets the range of incident energies for this intergration
void Compton::set_Ein_range( Compton_Ein_param *Ein_param )
{
  // sets the range based on the flux data
  Ein_param->set_Ein_range( );
  if( *(Ein_param->Ein_ptr) > Ein_param->data_E_0 )
  {
    Ein_param->data_E_0 = *(Ein_param->Ein_ptr);
  }
  if( Ein_param->quad_box.left_corner->E_in > Ein_param->data_E_0 )
  {
    Ein_param->data_E_0 = Ein_param->quad_box.left_corner->E_in;
  }
  if( *(Ein_param->next_Ein) < Ein_param->data_E_1 )
  {
    Ein_param->data_E_1 = *(Ein_param->next_Ein);
  }
  if( Ein_param->quad_box.right_corner->E_in < Ein_param->data_E_1 )
  {
    Ein_param->data_E_1 = Ein_param->quad_box.right_corner->E_in;
  }
  Ein_param->Ein_0 = Ein_param->data_E_0;
  Ein_param->Ein_1 = Ein_param->data_E_1;
}
// ---------------- Compton::update_T ------------------
// Adds the result of one integration
void Compton::update_T( T_matrix &transfer, Compton_Ein_param *Ein_param )
{
  // Don't bother if the interval is very small.
  static double E_tol = Global.Value( "E_tol" );
  if( Ein_param->Ein_1 - Ein_param->Ein_0 < E_tol*Ein_param->Ein_1 )
  {
    return;
  }
  static double tol = Global.Value( "quad_tol" );

  // a vector to store the integrals
  coef_vector value( transfer.order, transfer.conserve );
  value.set_zero( );

  // parameters for the integration
  QuadParamBase *params = static_cast< QuadParamBase* >( Ein_param );

  double left_E = Ein_param->Ein_0;
  double right_E = Ein_param->Ein_1;

  // evaluate the integral
  quad_F::integrate( Compton_F::Ein_F, transfer.Ein_quad_method, left_E, right_E,
		       params, tol, &value );
  // add this integral
  transfer( Ein_param->quad_box.Ein_count, Ein_param->quad_box.Eout_count ) += value;

  // increment the function counts
  Ein_param->Ein_F_count += Ein_param->func_count;
  ++Ein_param->quad_count;
}
// ---------------- Compton::get_xsec ------------------
// Calculates the cross section
void Compton::get_xsec( Quadrature_Method mu_quad_method, dd_vector& sigma )
{
  static double tol = Global.Value( "quad_tol" );
  static double Thompson = Global.Value( "Thompson" );  // Thompson cross section
  // storage for quadrature
  coef_vector value( 0, NUMBER );
  // parameters for the integration
  Compton_mu_param mu_params;
  QuadParamBase *params = static_cast< QuadParamBase* >( &mu_params );
  // As incident energies use the x-values with mu = -1
  for( Compton::const_iterator this_entry = begin( ); this_entry != end( );
       ++this_entry )
  {
    sigma.add_entry( this_entry->x, 0.0 );
  }

  long int quad_count = 0;
  // number of calls to Compton_F::sigma_F
  long int Ein_F_count = 0;

  // loop downward through the incident energies
  dd_vector::iterator next_sigma = sigma.end( );
  dd_vector::iterator this_sigma = next_sigma;
  --this_sigma;
  Compton::const_iterator this_x_entry;
  Compton::const_iterator next_x_entry;
  Compton::const_iterator last_x_entry = end( );
  --last_x_entry;
  for( ; next_sigma != sigma.begin( ); --last_x_entry, --next_sigma, --this_sigma )
  {
    mu_params.this_sigma = this_sigma;
    mu_params.next_sigma = next_sigma;
    // loop downward through the x-entries
    next_x_entry = last_x_entry;
    this_x_entry = next_x_entry;
    --this_x_entry;
    for( ; next_x_entry != begin( ); --next_x_entry, --this_x_entry )
    {
      // parameters for the integration over mu
      mu_params.set_E_in( this_sigma->x );
      mu_params.set_pair( *this_x_entry, *next_x_entry );
      // mu decreases as x increases
      double mu_0 = x_vector_F::get_mu_from_x( next_x_entry->x, this_sigma->x );
      double mu_1 = x_vector_F::get_mu_from_x( this_x_entry->x, this_sigma->x );
      quad_F::integrate( Compton_F::sigma_F, mu_quad_method, mu_0, mu_1,
			 params, tol, &value );
      this_sigma->y += value.weight_1[0];
      // increment the function counts
      Ein_F_count += mu_params.func_count;
      ++quad_count;
    }
  }
  this_sigma = sigma.begin( );
  this_x_entry = begin( );
  // if the first x and sigma energies are zero, use the Thompson cross section
  if( ( this_sigma->x == 0.0 ) && ( this_x_entry->x == 0.0 ) )
  {
    this_sigma->y = Thompson * this_x_entry->y;
  }
  // print the counts of function evaluations
  cout << "quadratures for cross section: " << quad_count << endl;
  cout << "Compton_F::sigma_F calls: " << Ein_F_count << endl;
  cout << "average Compton_F::sigma_F calls: " << 1.0*Ein_F_count/quad_count << endl;
  Ein_F_count = 0;
  quad_count = 0;

  // use lin-lin interpolation
  sigma.interp_type = LINLIN;
}

// **************** functions to integrate **********
// ----------------  Compton_F::sigma_F ------------------------
// Function for computing the cross section
void Compton_F::sigma_F( double mu, QuadParamBase *FF_param, coef_vector *value )
{
  static double m_electron = Global.Value( "m_electron" );
  static double Thompson = Global.Value( "Thompson" );
  // the parameters are really Compton_mu_param
  Compton_mu_param *params = static_cast<Compton_mu_param*>( FF_param );
  params->func_count += 1;
  double x = params->get_E_in() * sqrt( 0.5 * ( 1.0 - mu ) );
  double SF = params->value( x );   // the Compton scattering factor

  // get the Compton scattering cross section
  double kappa = params->get_E_in() / m_electron;
  double denominator = 1.0 + kappa * ( 1.0 - mu );
  double xs = 3.0 * Thompson * SF / ( 8 * denominator * denominator ) *
    ( 1.0 + mu*mu + kappa*kappa*( 1.0 - mu )*( 1.0 - mu ) / denominator );

  value->weight_1[0] = xs;
}
// Function for the 1-d quadrature over mu
// ---------------- Compton_F::mu_F ------------------
void Compton_F::mu_F( double mu, QuadParamBase *SF_param, coef_vector *value )
{
  static double m_electron = Global.Value( "m_electron" );
  static double Thompson = Global.Value( "Thompson" );
  // the parameters are really Compton_mu_param
  Compton_mu_param *params = static_cast<Compton_mu_param*>( SF_param );
  params->func_count += 1;
  double x = params->E_in * sqrt( 0.5 * ( 1.0 - mu ) );
  double SF = params->value( x );   // the Compton scattering factor

  // get the Compton scattering cross section
  double kappa = params->E_in / m_electron;
  double denominator = 1.0 + kappa * ( 1.0 - mu );
  double xs = 3.0 * Thompson * SF / ( 8 * denominator * denominator ) *
    ( 1.0 + mu*mu + kappa*kappa*( 1.0 - mu )*( 1.0 - mu ) / denominator );

  // the Legendre polynomial
  math_F::Legendre( mu, value );
  *value *= xs;
  // do the energy weighting if necessary
  if( ( value->conserve == ENERGY ) || ( value->conserve == BOTH ) )
  {
    double E_out =  params->E_in / denominator;
    value->scale_E( E_out );
  }
}
// Function for the 2-d quadrature over E_in
// ---------------- Compton_F::Ein_F ------------------
void Compton_F::Ein_F( double E_in, QuadParamBase *compton_param, coef_vector *value )
{
  // the parameters are really Compton_Ein_param
  Compton_Ein_param *params = static_cast<Compton_Ein_param*>( compton_param );
  params->func_count += 1;

  // parameters for the integration over mu
  Compton_mu_param SF_params;
  SF_params.set_pair( *params->x_ptr, *params->next_x );
  SF_params.E_in = E_in;

  // the range of integration
  double mu_0 = 0.0;
  double x_over_E;
  switch( params->bottom_ptr->boundary_ID )
  {
  case USE_MU_MINUS1:
    mu_0 = -1.0;
    break;
  case USE_X:
    x_over_E = params->next_x->x / E_in;
    mu_0 = 1.0 - 2*x_over_E*x_over_E;
    if( mu_0 < -1.0 ) mu_0 = -1.0;
    break;
  case USE_EOUT:
    mu_0 = params->Eout_ptr->hit_Ein( E_in );
    if( mu_0 < -1.0 ) mu_0 = -1.0;
    break;
  default:
    FatalError( "Compton_F::mu_F", "bad bottom pointer" );
  }

  double mu_1 = 0.0;
  switch( params->top_ptr->boundary_ID )
  {
  case USE_MU1:
    mu_1 = 1.0;
    break;
  case USE_X:
    x_over_E = params->x_ptr->x / E_in;
    mu_1 = 1.0 - 2*x_over_E*x_over_E;
    if( mu_1 < -1.0 ) mu_1 = -1.0;
    break;
  case USE_EOUT:
    mu_1 = params->next_Eout->hit_Ein( E_in );
    if( mu_1 < -1.0 ) mu_1 = -1.0;
    break;
  default:
    FatalError( "Compton_F::mu_F", "bad top pointer" );
  }

  // evaluate the integral over mu
  QuadParamBase *mu_params = static_cast< QuadParamBase* >( &SF_params );
  static double tol = Global.Value( "quad_tol" );
  quad_F::integrate( Compton_F::mu_F, params->mu_quad_method, mu_0, mu_1, 
       	       mu_params, tol, value );

  // weight it by flux
  params->flux_weight( E_in );
  *value *= params->current_weight;
  params->mu_F_count += SF_params.func_count;
}
