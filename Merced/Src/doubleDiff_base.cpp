/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2021-10-05 $
 * $Author: hedstrom $
 * $Id: doubleDiff_base.cpp 1 2021-10-05Z hedstrom $
 *
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// Implementation of the basic classes for pointwise energy-angle
// probability density

#include <cmath>

#include "doubleDiff_base.hpp"
#include "global_params.hpp"
#include "messaging.hpp"


// ************* class DDbase::one_Eout ****************************

// ************* class DDbase::one_joint_dist **********************
// ----------- DDbase::one_joint_dist::check_norm --------------
// Checks the norm of the data
void DDbase::one_joint_dist::check_norm( )
{
  DDbase::one_joint_dist::iterator this_Eout = begin( );
  double sum = 0.0;
  // use the trapezoid rule---it's lin-lin interpolation in outgoing energy
  double prev_norm = this_Eout->get_norm( );
  double next_norm;
  double prev_E = this_Eout->get_Eout( );
  double next_E;

  for( ; this_Eout != end( ); ++this_Eout, prev_norm = next_norm,
      prev_E = next_E )
  {
    next_norm = this_Eout->get_norm( );
    next_E = this_Eout->get_Eout( );
    sum += 0.5 * ( prev_norm + next_norm ) * ( next_E - prev_E );
  }

  static bool warned = false;
  static double norm_tol = Global.Value( "norm_tol" );

  if( !warned && ( std::abs( sum - 1.0 ) > norm_tol ) )
  {
    Msg::Warning("DDbase::one_joint_dist::check_norm",
            Msg::pastenum( "bad norm in data: ", sum ) );
    warned = true;
  }

  double scale_Eout = 1.0 / sum;
  for( this_Eout = begin( ); this_Eout != end( ); ++this_Eout )
  {
    *this_Eout *= scale_Eout;
  }
}
// ----------- DDbase::one_joint_dist::to_unit_base --------------
// Maps the outgoing energies to [0, 1]
void DDbase::one_joint_dist::to_unit_base( )
{
  DDbase::one_joint_dist::iterator this_Eout = begin( );
  DDbase::one_joint_dist::iterator last_Eout = end( );
  --last_Eout;
  Eout_ubase_map.Eout_min = this_Eout->get_Eout( );
  Eout_ubase_map.Eout_max = last_Eout->get_Eout( );
  double scale_Eout = Eout_ubase_map.Eout_max - Eout_ubase_map.Eout_min;

  for( ; this_Eout != end( ); ++this_Eout )
  {
    bool interp_OK;
    this_Eout->set_Eout( Eout_ubase_map.to_unit_base( this_Eout->get_Eout( ),
						      &interp_OK ) );
    *this_Eout *= scale_Eout;
  }
}

// ************* class DDbase::joint_dist_base *********************
// ----------- DDbase::joint_dist_base::read_data --------------
// Reads double-differential energy-angle data
void DDbase::joint_dist_base::read_data( Dpar::data_parser &inFile, int num_Ein )
{
  interp_flag_F::read_3d_interpolation_GND( inFile, &Ein_interp, &Eout_interp,
                       &mu_interp );

  if( Eout_interp.qualifier == Terp::UNITBASE )
  {
    // The cosine ought to range over -1 <= mu <= 1
    Msg::FatalError( "DDbase::joint_dist_base::read_data",
       "unit-base interpolation of outgoing energy not implemented" );
    //new_Eout_ptr->unit_base( false, &new_Eout_ptr->ubase_map );
  }

  DDbase::joint_dist_base::iterator new_Ein_ptr;
  DDbase::one_joint_dist::iterator new_Eout_ptr;

  double mu;
  double Prob;

  // read the data
  for( int Ein_count = 0; Ein_count < num_Ein; ++Ein_count )
  {
    // make a new link for this incident energy
    new_Ein_ptr = insert( end( ), DDbase::one_joint_dist( ) );
    new_Ein_ptr->set_E_in( inFile.get_next_double( ) );

    // loop over outgoing energy
    int num_Eout = inFile.get_next_int( );
    for( int Eout_count = 0; Eout_count < num_Eout; ++Eout_count )
    {
      // make a new energy distribution
      new_Eout_ptr = new_Ein_ptr->insert( new_Ein_ptr->end( ), DDbase::one_Eout( ) );
      new_Eout_ptr->set_Eout( inFile.get_next_double( ) );
      new_Eout_ptr->interp_type = mu_interp;
      // read the (cosine, probability density) pairs
      int num_mu = inFile.get_next_int( );
      for( int mu_count = 0; mu_count < num_mu; ++mu_count )
      {
        mu = inFile.get_next_double( );
        Prob = inFile.get_next_double( );
        new_Eout_ptr->add_entry( mu, Prob );
      }
    }
  }
  check_norm( );

  if( Ein_interp.qualifier == Terp::UNITBASE )
  {
    to_unit_base( );
  }
}

// ----------- DDbase::joint_dist_base::check_norm --------------
// Checks the norm of the data
void DDbase::joint_dist_base::check_norm( )
{
  DDbase::joint_dist_base::iterator joint_ptr = begin( );
  for( ; joint_ptr != end( ); ++joint_ptr )
  {
    joint_ptr->check_norm( );
  }
}

// ----------- DDbase::joint_dist_base::to_unit_base --------------
// Maps the outgoing energies to 0 <= mu <= 1
void DDbase::joint_dist_base::to_unit_base( )
{
  DDbase::joint_dist_base::iterator joint_ptr = begin( );
  for( ; joint_ptr != end( ); ++joint_ptr )
  {
    joint_ptr->to_unit_base( );
  }
}
