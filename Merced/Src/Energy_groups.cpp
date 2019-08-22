/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2014-08-05 19:06:56 -0800 (Tue, 05 Aug 2014) $
 * $Author: hedstrom $
 * $Id: Energy_groups.hpp 1 2014-08-05 03:06:56Z hedstrom $
 *
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
// implementation of the class for energy groups

#include "Energy_groups.hpp"
#include "messaging.hpp"

//****************** class Energy_groups  *****************
//--------------- Energy_groups::read_bd -----------------
// Reads the energy group boundaries
void Energy_groups::read_bd( data_parser &input_file, int num_bd )
{
  for( int bd_count = 0; bd_count < num_bd; ++bd_count )
  {
    double bound = input_file.get_next_double( );
    push_back( bound );
  }
}
//--------------- Energy_groups::first_bin_ID -----------------
// Returns the index of the left-hand end of the energy bin containing this_E
int Energy_groups::first_bin_ID( double this_E ) const
{
  int bin_ID = 0;
  Energy_groups::const_iterator this_bin = begin( );
  for( ++this_bin; this_bin != end( ); ++this_bin, ++bin_ID )
  {
    if( *this_bin > this_E )
    {
      break;
    }
  }
  if( this_bin == end( ) )
  {
    FatalError( "Energy_groups::first_bin_ID", "energy too high" );
  }
  return bin_ID;
}
//--------------- Energy_groups::last_bin_ID -----------------
// Returns the index of the right-hand end of the energy bin containing this_E
int Energy_groups::last_bin_ID( double this_E ) const
{
  int bin_ID = size( ) - 1;
  Energy_groups::const_iterator this_bin = end( );
  --this_bin;
  Energy_groups::const_iterator lower_bin = this_bin;
  --lower_bin;
  for( ; this_bin != begin( ); this_bin = lower_bin, --lower_bin, --bin_ID )
  {
    if( *lower_bin <= this_E )
    {
      break;
    }
  }
  if( this_bin == begin( ) )
  {
    FatalError( "Energy_groups::first_bin_ID", "energy too high" );
  }
  return bin_ID;
}
