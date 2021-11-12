/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2021-09-29 $
 * $Author: hedstrom $
 * $Id: lab_doubleDiff.cpp 1 2021-09-21Z hedstrom $
 *
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// Implementation of the classes used to handle pointwise energy-angle
// probability density given in the lab frame.

#include "lab_doubleDiff.hpp"
#include "global_params.hpp"

// *************class labDD::lab_joint_dist *********************
// ------------- labDD::lab_joint_dist::to_Legendre --------------
// converts to Legendre data
void labDD::lab_joint_dist::to_Legendre( )
{
  int mu_F_count = 0;
  int from_table_count = 0;

  // the desired Legendre order
  int order = Global.Value( "outputLegendreOrder" );
  
  // interpolation rules
  Legendre_data.Ein_interp = Ein_interp;
  Legendre_data.Eout_interp = Eout_interp.flag;
  Legendre_data.order = order;

  // iterate over incident energy
  for( labDD::lab_joint_dist::const_iterator table_ptr = begin( );
       table_ptr != end( ); ++table_ptr )
  {
    // make a new StdLg::standard_Legendre_vector
    StdLg::standard_Legendre::iterator Legendre_ptr =
      Legendre_data.insert( Legendre_data.end( ),
			     StdLg::standard_Legendre_vector( ) );
    Legendre_ptr->set_E_in( table_ptr->get_E_in( ) );  // energy of incident particle
    Legendre_ptr->Ein_interp = Ein_interp;
    Legendre_ptr->Eout_interp = Eout_interp.flag;
    Legendre_ptr->order = order;

    Legendre_ptr->ubase_map.copy( table_ptr->Eout_ubase_map );

    // iterate over outgoing energy
    for( DDbase::one_joint_dist::const_iterator table_Eout = table_ptr->begin( );
	 table_Eout != table_ptr->end( ); ++table_Eout )
    {
      // insert new Legendre coefficients
      Lgdata::Legendre_list_base::iterator Legendre_Eout =
        Legendre_ptr->insert( Legendre_ptr->end( ),
			     Lgdata::Legendre_coefs( ) );

      Legendre_Eout->initialize( order );
      
      mu_F_count += to_Legendre_F::from_table( *table_Eout,
					       &(*Legendre_Eout) );
      ++from_table_count;
    }

  }
  // print the counts
  std::cout << "quadratures over mu: " << from_table_count << std::endl;
  std::cout << "to_Legendre_F::mu_F calls: " << mu_F_count << std::endl;
  std::cout << "average to_Legendre_F::mu_F calls: " <<
    mu_F_count / from_table_count << std::endl;
}

// ------------- labDD::lab_joint_dist::get_T --------------
// alculates the transfer matrix for this particle.
void labDD::lab_joint_dist::get_T( const Ddvec::dd_vector& sigma,
				  const Ddvec::dd_vector& mult,
                const Ddvec::dd_vector& weight,
				  Trf::T_matrix& transfer )
{
  // Convert the angular table to Legendre coefficients
  to_Legendre( );

  // Get the transfer matrix
  Legendre_data.get_T( sigma, mult, weight, transfer );
}
