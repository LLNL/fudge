/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: general_evap.hpp 1 2010-06-18 03:06:56Z hedstrom $
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
//! Classes used for outgoing an energy distribution independent of the initial energy

#ifndef GENERAL_EVAP_DEF
#define GENERAL_EVAP_DEF

#include "dd_vector.hpp"
#include "transfer.hpp"
#include "param_base.hpp"

class general_evap_param;  // forward reference

//! Class for a single energy spectrum
// ---------------- class general_evap ------------------
class general_evap : public dd_vector
{
private:
  int first_Ein;  // index of the left-hand end of the first significant energy bin
  int last_Ein;  // index of the right-hand end of the last significant energy bin

  //! A list of integrals over E_out, independent of incident energy
  list< coef_vector > E_out_ints;

  //!  Gets the range of nontrivial incident energy bins; computes first_Ein and last_Ein
  //! returns true if the threshold is too high for the energy bins
  //! \param sigma the cross section data
  //! \param multiple the multiplicity of the outgoing particle
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param e_flux approximate flux used to weight the transfer matrix
  //! \param Ein_groups the boundaries of the incident energy groups
  bool get_Ein_range( const dd_vector& sigma, const dd_vector& multiple,
    const dd_vector& weight,
    const Flux_List& e_flux, const Energy_groups& Ein_groups );

  //! Sets up the integrals over the E_out bins
  //! \param order the Legendre order of the output transfer matrix
  //! \param conserve the flag for conservation of energy or particle number
  //! \param Eout_groups the boundaries of the outgoing energy groups
  void setup_Eout_ints( int order, Conserve conserve,
    const vector< double >& Eout_groups );

  //! The probability integals for outgoing energy should add to 1
  void check_sum( );

  //! Go to the next pair of incident energies.  Returns "true" when finished.
  //! \param E_in the energy of the incident particle
  //! \param Ein_param parameters for intgration over incident energy
  bool next_ladder( double E_in, general_evap_param *Ein_param );

  //! Adds to the transfer matrix for all E_out bins for a pair of incident energies.
  //! \param transfer the computed transfer matrix
  //! \param Ein_param parameters for intgration over incident energy
  void Eout_ladder( T_matrix& transfer, general_evap_param *Ein_param );

  //! Adds to an element of transfer the integral over the E-E' box
  //! \param transfer the computed transfer matrix
  //! \param Eout_count current row number of the computed transfer matrix
  //! \param Ein_param parameters for intgration over incident energy
  void update_T( T_matrix &transfer, int Eout_count, general_evap_param *Ein_param );

  //! The code currently handles only theta = 1
  bool theta_OK( );

public:

  //! The temperature parameter
  dd_vector theta;

  general_evap() {}
  ~general_evap() {}

  //! Calculates the transfer matrix for this particle
  void get_T( const dd_vector& sigma, const dd_vector& multiple, 
	      const dd_vector& weight, T_matrix& transfer );
};

//! Class for parameters for a single energy spectrum
// ---------------- class general_evap_param ------------------
class general_evap_param : public param_base
{
private:

public:
  //! which outgoing energy bin
  list< coef_vector >::iterator current_Eout_int;

  long int quad_count;  // number of 2-d quadratures
  long int Ein_F_count;  // number of calls to general_evap_F::E_quad_F

  general_evap_param(): quad_count( 0 ), Ein_F_count( 0 ) {}
  ~general_evap_param() {}

  //! Gets the integrals over this E_out bin
  //! \param value the computed integrals of probability density over outgoing energy
  void get_integrals( coef_vector *value );
};

namespace general_evap_F
{
  // ************** general_evap_F::Ein_F ******************************
  //! Integral function for the model
  void Ein_F( double E_in, QuadParamBase *e_quad_param,
    coef_vector *value );
}

#endif
