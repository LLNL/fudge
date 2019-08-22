/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2016-02-12 11:40:00 -0800 (Fri, 12 Feb 2016) $
 * $Author: hedstrom $
 * $Id: isotropic.hpp 1 2016-02-12 11:40:00Z hedstrom $
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
// implement the classes for isotropic distributions in the lab frame.

#ifdef _OPENMP
 #include <omp.h>
#endif

#include "isotropic.hpp"
#include "messaging.hpp"
#include "global_params.hpp"

// *********************** class isotropic *******************
// ------------------ isotropic::read_data -----------------
void isotropic::read_data( data_parser& input_file, int num_Ein )
{
  data_order = 0;
  // this is a kludge to use old code
  ENDL_data = new energy_dist[ data_order + 1 ];
  energy_dist *new_moment_ptr = &ENDL_data[ 0 ];
  new_moment_ptr->read_data( input_file, num_Ein );

  // we need to set the outgoing interpolation
  Eprob_vector *Ein_ptr;
  for( int Ein_count = 0; Ein_count < num_Ein; ++Ein_count )
  {
    Ein_ptr = &new_moment_ptr->EProb_data[ Ein_count ];
    Ein_ptr->interp_type = Eout_interp;
    // ensure proper normalization
    Ein_ptr->renorm( );
  }

  // the energy_dist read_data doesn't do unit-base interpolation
  if( Ein_interp.qualifier == UNITBASE )
  {
    new_moment_ptr->unit_base( 0 );
  }
}
// ------------------ isotropic::get_T -----------------
// Calculates the transfer matrix for this particle.
// sigma is the cross section.
void isotropic::get_T( const dd_vector& sigma,
  const dd_vector& mult, const dd_vector& weight, T_matrix& transfer )
{ 
  bool interp_OK = ( ( Ein_interp.qualifier == UNITBASE ) &&
		     ( ( Ein_interp.flag == LINLIN ) ||
     		       ( Ein_interp.flag == LINLOG ) ) ) ||
                   ( ( Ein_interp.qualifier == CUMULATIVE_POINTS ) &&
		     ( ( Ein_interp.flag == LINLIN ) ||
		       ( Ein_interp.flag == LINLOG ) ) ) ||
    ( ( Ein_interp.qualifier == DIRECT ) &&
      ( ( Ein_interp.flag == LINLIN ) ||
        ( Ein_interp.flag == HISTOGRAM ) ) );

  if( !interp_OK )
  {
    FatalError( "isotropic::get_T",
      "Incident energy interpolation not implemented" );
  }
  interp_OK = ( Eout_interp == LINLIN ) || ( Eout_interp == HISTOGRAM );
  if( !interp_OK )
  { 
    FatalError( "isotropic::get_T",
      "Outgoing energy interpolation not implemented" );
  }

  // convert to ENDF format
  zero_order( );

  ENDF_data.Ein_interp = Ein_interp;
  ENDF_data.Eout_interp = Eout_interp;
  ENDF_data.get_T( sigma, mult, weight, transfer );
}
