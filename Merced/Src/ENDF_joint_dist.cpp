/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: ENDF_joint_dist.cpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
  Copyright (c) 2017, Lawrence Livermore National Security, LLC.
  Produced at the Lawrence Livermore National Laboratory.
  Written by the LLNL Nuclear Data and Theory group
          (email: mattoon1@llnl.gov)
  LLNL-CODE-725546.
  All rights reserved.
  
  This file is part of the Merced package, used to generate nuclear reaction
  transfer matrices for deterministic radiation transport.
  
  
      Please also read this link - Our Notice and Modified BSD License
  
  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:
      * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the disclaimer below.
      * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the disclaimer (as noted below) in the
        documentation and/or other materials provided with the distribution.
      * Neither the name of LLNS/LLNL nor the names of its contributors may be used
        to endorse or promote products derived from this software without specific
        prior written permission.
  
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
  THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  
  
  Additional BSD Notice
  
  1. This notice is required to be provided under our contract with the U.S.
  Department of Energy (DOE). This work was produced at Lawrence Livermore
  National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
  
  2. Neither the United States Government nor Lawrence Livermore National Security,
  LLC nor any of their employees, makes any warranty, express or implied, or assumes
  any liability or responsibility for the accuracy, completeness, or usefulness of any
  information, apparatus, product, or process disclosed, or represents that its use
  would not infringe privately-owned rights.
  
  3. Also, reference herein to any specific commercial products, process, or services
  by trade name, trademark, manufacturer or otherwise does not necessarily constitute
  or imply its endorsement, recommendation, or favoring by the United States Government
  or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
  herein do not necessarily state or reflect those of the United States Government or
  Lawrence Livermore National Security, LLC, and shall not be used for advertising or
  product endorsement purposes.
  
 * # <<END-copyright>>
*/
// implementation of the classes used to handle ENDF_joint energy-angle distributions

#include <iostream>

#include "ENDF_joint_dist.hpp"
#include "messaging.hpp"

using namespace std;

// ********* class ENDF_joint_dist *********
// ----------- ENDF_joint_dist::read_data --------------
void ENDF_joint_dist::read_data( data_parser &inFile, int num_Ein )
{
  interp_flag_F::read_3d_interpolation( inFile, &Ein_interp, &mu_interp,
		       &Eout_interp );

  ENDF_joint_dist::iterator new_joint_ptr;
  one_joint_dist::iterator new_e_dist_ptr;

  double E_out;
  double Prob;

  // read the data
  for( int Ein_count = 0; Ein_count < num_Ein; ++Ein_count )
  {
    // make a new link for this incident energy
    new_joint_ptr = insert( end( ), one_joint_dist( ) );
    new_joint_ptr->set_E_in( inFile.get_next_double( ) );
    new_joint_ptr->mu_interp = mu_interp;
    new_joint_ptr->Eout_interp = Eout_interp;
    new_joint_ptr->mu_quad_param.mu_interp = mu_interp;
    new_joint_ptr->mu_quad_param.Eout_interp = Eout_interp;

    // loop over cosine
    int num_mu = inFile.get_next_int( );
    for( int mu_count = 0; mu_count < num_mu; ++mu_count )
    {
      // make a new energy distribution
      new_e_dist_ptr = new_joint_ptr->insert( new_joint_ptr->end( ), one_mu( ) );
      new_e_dist_ptr->set_mu( inFile.get_next_double( ) );
      new_e_dist_ptr->interp_type = Eout_interp;
      // read the (energy, probability density) pairs
      int num_Eout = inFile.get_next_int( );
      for( int Eout_count = 0; Eout_count < num_Eout; ++Eout_count )
      {
        E_out = inFile.get_next_double( );
        Prob = inFile.get_next_double( );
        new_e_dist_ptr->add_entry( E_out, Prob );
      }
    }
  }
  // we need to set up the ENDL I=1 angular probability densities
  get_angular_data( );
}
// ----------- ENDF_joint_dist::get_angular_data --------------
void ENDF_joint_dist::get_angular_data( )
// set up the i=1 angular probability densities
{
  // fill angle_data
  angle_data.Ein_interp = Ein_interp;
  angle_data.mu_interp = mu_interp;
  for( ENDF_joint_dist::iterator joint_mu = begin( ); joint_mu != end( ); ++joint_mu )
  {
    angle_dist::iterator angle_dist_ptr = angle_data.insert( angle_data.end( ), dd_vector( ) );
    angle_dist_ptr->set_E_in( joint_mu->get_E_in( ) );
    angle_dist_ptr->interp_type = LINLIN;
    for( one_joint_dist::iterator e_mu_dist_ptr = joint_mu->begin( );
         e_mu_dist_ptr != joint_mu->end( ); ++e_mu_dist_ptr )
    {
      double mu = e_mu_dist_ptr->get_mu( );
      double norm = e_mu_dist_ptr->get_norm( );
      angle_dist_ptr->add_entry( mu, norm );
      if( ( mu_interp == UNITBASE ) && ( norm > 0 ) )
      {
	e_mu_dist_ptr->unit_base( true, &e_mu_dist_ptr->ubase_map );
      }
    }
    angle_dist_ptr->renorm( );
  }
}
