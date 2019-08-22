/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: transfer.hpp 1 2006-02-02 03:06:56Z hedstrom $
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
// the class for transfer matrices
#ifndef T_MATRIX_DEF
#define T_MATRIX_DEF

#include <vector>
#include <string>
#include <list>
#include <fstream>

#include "coef_vector.hpp"
#include "Legendre_data.hpp"
#include "data_parser.hpp"
#include "quadrature.hpp"
#include "Energy_groups.hpp"

using namespace std;

//! Class for the transfer matrix
//-----------class T_matrix ----------
class T_matrix
{
public:
  coef_vector *data;
  coef_vector *row_checks;  // used for verification
  coef_vector *row_sums;  // used for verification
  Flux_List e_flux;  // the Legendre coefficients of the incident flux
  int order;   // the Legendre order
  int num_Ein_bins;
  int num_Eout_bins;
  Conserve conserve;
  Quadrature_Method Ein_quad_method;  // quadrature method for incident energy
  Quadrature_Method Eout_quad_method;  // quadrature method for outgoing energy
  Quadrature_Method mu_quad_method;  // quadrature method for outgoing cosine

  bool interpolate_Eout_integrals;  // do we interpolate the integrals over Eout?

  //! The energy groups
  Energy_groups in_groups;
  Energy_groups out_groups;

  //! the weights 1 / \int_(energy group) flux by Legendre order
  weight_list flux_weight;

  //! Average cross section over energy bins, for checking the contribution of an interval
  vector< double > averageCrossSections;

  //! The maximum average cross section
  double maximumCrossSection;
  ofstream *output_file;  // pointer to the output file

  //! Constructor
  T_matrix( );

  //! Destructor
  ~T_matrix( );

  //! Allocates space
  void allocate( );

  //! Operator to grab an entry:
  //! \param Ein_count: incident energy group
  //! \param Eout_count: exit energy group
  coef_vector& operator()( int Ein_count, int Eout_count );

  //! Operator to grab an entry:
  //! \param Ein_count: incident energy group
  //! \param Eout_count: exit energy group
  //! \param ell: the Legendre order of this term
  coef_vector& operator()( int Ein_count, int Eout_count, int ell );

  //! Calculates the weights 1 / \int_(energy group) flux
  void get_flux_weight( );

  //! Applies the weights 1 / \int_(energy group) flux
  void use_weight( );

  //! Scales the weight_E terms by the average E_out, to check gamma output
  void scale_E( );

  // Prints the sums of the rows of the zero-order term.
  // They should agree with the flux-weighted average cross sections
  void check_ell0( );

  // Scales the row check sums by the weights 
  void scale_row_check( );

  //! Computes the net cross section for each energy bin
  //! \param sigma a vector of pairs ( incident energy, cross section )
  void getBinCrossSection( const dd_vector& sigma );

  //! Prints the matrix to the output file
  void write_transfer( );

  //! Prints zeros to the output file for a reaction with high threshold
  void zero_transfer( );
};

#endif
