/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: energy_dist.cpp 1 2006-02-02 03:06:56Z hedstrom $
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
// implementation of the base classes used to handle energy distributions

#include <cmath>

#include "energy_dist_base.hpp"
#include "messaging.hpp"
#include "global_params.hpp"

using namespace std;

// ************* class Eprob_vector *****************
// ----------- Eprob_vector::unit_base --------------
// Transform one energy distribuiton to unit base
void Eprob_vector::unit_base( int L_order )
{
  bool Renorm = ( L_order == 0 );
  dd_vector::unit_base( Renorm, &ubase_map );
}
// ----------- Eprob_vector::form_cum_prob --------------
// Forms the list of cumulative probabilities
void Eprob_vector::form_cum_prob( )
{
  // copy the data
  cum_prob.Eout_interp = interp_type;
  for( Eprob_vector::const_iterator Eout_ptr = begin( );
       Eout_ptr != end( ); ++Eout_ptr )
  {
    cumulative_prob_list::iterator cum_prob_ptr = cum_prob.insert(
      cum_prob.end( ), cumulative_prob_entry( ) );
    cum_prob_ptr->E_out = Eout_ptr->x;
    cum_prob_ptr->Prob = Eout_ptr->y;
  }
  // now form the slopes and cumulative probabilities
  if( interp_type == HISTOGRAM )
  {
    cum_prob.get_cum_prob_flat( );
  }
  else // lin-lin
  {
    cum_prob.get_cum_prob_linlin( );
  }
}

// ************* class energy_dist_base *****************
// ----------- energy_dist_base::print --------------
void energy_dist_base::print( )
{
  for( energy_dist_base::iterator energy_ptr = begin( );
       energy_ptr != end( ); ++energy_ptr )
  {
    energy_ptr->print( );
  }
}
// ----------- energy_dist_base::unit_base --------------
// Transform the outgoing energy distribution to the interval [0, 1]
void energy_dist_base::unit_base( int L_order )
{
  // do the mapping for each incident energy
  for( energy_dist_base::iterator Ein_ptr = begin( );
       Ein_ptr != end( ); ++Ein_ptr )
  {
    Ein_ptr->unit_base( L_order );
  }
}
