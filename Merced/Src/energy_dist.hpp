/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2011-04-29 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: energy_dist.hpp 1 2011-04-29 03:06:56Z hedstrom $
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
// define the classes used for ENDL i=4 energy distributions

#ifndef ENERGY_DIST_CLASS
#define ENERGY_DIST_CLASS

#include <iostream>
#include "data_parser.hpp"
#include "standard_Legendre.hpp"
#include "energy_dist_base.hpp"

using namespace std;

//! Class for energy distributions, one Legendre order
//--------------- class energy_dist ----------------
class energy_dist
{
private:

public:
  //! the ENDL data
  Eprob_vector *EProb_data;

  //! the number of incident energies
  int number_Ein;
  two_d_interp Ein_interp;  // interpolation between incident enrgies
  Interp_Type Eout_interp;  // interpolation between outgoing enrgies

  inline energy_dist( ): number_Ein( 0 ) {}

  ~energy_dist( );

  //! Reads the ENDL data
  //! \param infile input file
  //! \param num_Ein number of incident energies for this reaction
  void read_data( data_parser& infile, int num_Ein );

  //! Maps the data to unit base
  //! \param L_order the Legendre order of the current data
  void unit_base( int L_order );
};

//! list: one link for each Legendre order
//--------------- class energy_moments ----------------
class energy_moments
{
 private:
  int output_order;  // the Legendre order of the output

  //! Converts ENDL data to ENDF format
  void to_ENDF( );

  //! Converts ENDL data to ENDF format for one incident energy
  //! \param Ein_count identifies the current incident energy bin
  void one_Ein_to_ENDF( int Ein_count );

  //! Checks to see that the incident energies are consistent for all Legendre orders
  void check_Ein( ) const;

 protected:
  //! Use the coding for this data in ENDF format
  standard_Legendre ENDF_data;

  //! The original ENDL data
  energy_dist *ENDL_data;

  int data_order;  // the Legendre order of the data
 
  //! Converts isotropic ENDL data to ENDF format
  void zero_order( );

 public:
  two_d_interp Ein_interp;  // interpolation between incident enrgies
  Interp_Type Eout_interp;  // interpolation between outgoing enrgies

  inline energy_moments( ): output_order( -1 ), data_order( -1 ) {}
  ~energy_moments( );

  //! Reads the ENDL data
  //! \param infile input file
  //! \param num_moments number of Legendre moments for this reaction
  void read_data( data_parser& input_file, int num_moments );

   // Calculates the transfer matrix for this particle.
  //! \param sigma the cross section data
  //! \param mult the outgoing particle multiplicity data
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param transfer the transfer matrix
  void get_T( const dd_vector& sigma, const dd_vector& mult,
    const dd_vector& weight, T_matrix& transfer );

  // Prints the lists for debugging
  void print( );
};

#endif
