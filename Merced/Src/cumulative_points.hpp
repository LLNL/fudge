/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2011-03-07 (Mon, Mar 7, 2011) $
 * $Author: hedstrom $
 * $Id: cumulative_points.hpp 1 2011-03-07 hedstrom $
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
// classes used to handle interpolation by cumulative points

#ifndef CUMULATIVE_POINTS_DEF
#define CUMULATIVE_POINTS_DEF

#include <list>

#include "dd_vector.hpp"  // for Interp_Type

using namespace std;

//! Class for cumulative probability at one outgoing energy
//--------------- class cumulative_prob_entry ----------------
class cumulative_prob_entry
{
public:
  // ! the outgoing energy
  double E_out;

  //! the probability density
  double Prob;

  //! The derivative of the zero-order coefficient with respect to outgoing energy
  double slope;

  //! The cumulative probability up to this entry
  double cum_prob;

  inline cumulative_prob_entry( ) {}

  inline ~cumulative_prob_entry( ) {}

  //! Gets the probability density at incident energy E_in
  //! \param E the current outgoing energy
  double get_prob( double E );

  //! Gets the cumulative probability at incident energy E_in
  //! \param E the current outgoing energy
  double get_cum_prob( double E );

  //! Gets the energy corresponding to cumulative probability A
  //! \param A the current cumulative probability
  double get_cum_inv( double A ) const;
};

//! Class for cumulative probabilities at one incident energy
//--------------- class cumulative_prob_list ----------------
class cumulative_prob_list : public list< cumulative_prob_entry >
{
 private:
  //! The energy of the incident particle
  double E_in;

 public:
  Interp_Type Eout_interp;

  cumulative_prob_list( ) {}

  ~cumulative_prob_list( ) {}

  //! Returns the energy of the incident particle
  inline double get_E_in( ) const { return E_in; }

  //! Sets the energy of the incident particle
  //! \param E_in energy of incident particle
  inline void set_E_in( double Ein ) { E_in = Ein; }

  //! Computes the cumulative probabilities and slopes for histogram data
  void get_cum_prob_flat( );

  //! Computes the cumulative probabilities and slopes for lin-lin data
  void get_cum_prob_linlin( );

};


#endif
