/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: coef_vector.hpp 1 2009-09-01 03:06:56Z hedstrom $
 * initial code
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

// header for the coef_vector class
#ifndef COEF_VECTOR
#define COEF_VECTOR

#include "max_output_order.hpp"
#include "Legendre_data.hpp"

//! Do we conserve particle number, energy, or both?
enum Conserve{ NUMBER, ENERGY, BOTH, NOT_SET };

// --------------- coef_vector -------------------------
//! Class to hold Legendre coefficients j=0, 1, ..., order
class coef_vector
{
private:

public:
  //! Coefficients for conservation of number of particles
  double weight_1[ max_output_order+1 ];

  //! Coefficients for conservation of energy
  double weight_E[ max_output_order+1 ];

  //! actual Legendre order
  int order;

  //! flag for conservation of particle number or energy
  Conserve conserve;

  //! Constructor
  inline coef_vector( ): order( -1 ), conserve( NOT_SET ) {}

  //! Constructor
  coef_vector( int Order, Conserve cons );

  //! Destructor
  inline ~coef_vector() {}

  //! Sets the order and conservation flag
  //! \param Order the Legendre order of the output
  //! \param cons the conservation flag
  void set_order( int Order, Conserve cons );

  //! Sets the entries to zero
  void set_zero( );

  //! Ensure that our rough estimate is nonzero
  //! param Norm_1 the max norm of weight_1, set to length if initially 0
  //! param Norm_E the max norm of weight_E, set to length if initially 0
  //! param length default value for norms initially 0
  void test_zero( double *Norm_1, double *Norm_E, double length );

  //! Makes a copy
  //! param to_copy the values to copy
  void copy( const coef_vector &to_copy );

  //! Makes a copy
  //! param to_copy the values to copy
  coef_vector& operator=( const coef_vector& to_copy );

  //! Does vector addition
  //! param to_add the values to add
  coef_vector& operator+=( const coef_vector& to_add );

  //! Adds a scalar to the terms of order L_order
  //! param to_add the values to add (one order)
  //! param L_order the Legendre order to add to
  coef_vector& plus( const coef_vector& to_add, int L_order );

  //! Scales the vector
  //! param factor the scale factor
  coef_vector& operator*=( double factor );

  //! Scales the weight_E terms
  //! param factor the scale factor for weight_E
  void scale_E( double factor );

  //! Weights the vector 
  //! param 
  coef_vector& operator*=( Legendre_base &factor );

  //! Calculates the max norms of weight_1 and weight_E
  //! Replaces the entries by their absolute values
  //! param Norm_1 the max norm of weight_1
  //! param Norm_E the max norm of weight_E
  void max_norm( double *Norm_1, double *Norm_E );

  void print( );
};

#endif
