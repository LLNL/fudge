/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: Maxwell.hpp 1 2006-02-02 03:06:56Z hedstrom $
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
//! Classes used for Maxwell energy probability density

#ifndef MAXWELL_DEF
#define MAXWELL_DEF

#include "evaporation.hpp"

class Maxwell_param;  // forward reference

//! Class for the maxwell model
// ---------------- class maxwell ------------------
class maxwell : public evaporation
{
private:

public:
  maxwell( ) {}
 
  ~maxwell( ) {}

  //! Calculates the transfer matrix for this particle
  //! \param sigma the cross section data
  //! \param multiple the multiplicity of the outgoing particle
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param transfer the computed transfer matrix
  void get_T( const dd_vector& sigma, const dd_vector& multiple,
    const dd_vector& weight, T_matrix& transfer );
};

//! Class for parameters for the Maxwell model
// ---------------- class Maxwell_param ------------------
class Maxwell_param : public Evap_param
{
private:
  //! Integral of the probability density from 0 to E_1
  //! \param E_1 upper limit of integration in outgoing energy
  double integral_prob( double E_1 );

  //! Integral of energy times the probability density from 0 to E_1
  //! \param E_1 upper limit of integration in outgoing energy
  double integral_Eprob( double E_1 );

public:
  inline Maxwell_param() {}
  inline ~Maxwell_param() {}

  // *** Implement the virtual functions ***
  //! Gets the integrals over outgoing energy
  //! \param Eout_0 lower limit of the integral over outgoing energy
  //! \param Eout_1 upper limit of the integral over outgoing energy
  //! \param value the computed integral
  void get_integrals( double Eout_0, double Eout_1, coef_vector &value );

  //! Gets the integrals over outgoing energy and returns the noise in the calculation
  //! \param Eout_0 lower limit of the integral over outgoing energy
  //! \param Eout_1 upper limit of the integral over outgoing energy
  //! \param value the computed integral
  inline double tol_get_integrals( double Eout_0, double Eout_1, coef_vector &value )
  {
    return Evap_param::tol_get_integrals( Eout_0, Eout_1, value );
  }

  //! Integrate from 0 to E_max to get the norm
  double get_norm( );

  //! Sets the scale factors for the integrals of probability and energy*probability
  void set_scales( );

  //! Sets the tolerance for the quadrature over incident energy
  //! \param left_E lower limit of the integral over incident energy
  //! \param right_E upper limit of the integral over incident energy
  inline double set_tol( double left_E, double right_E )
  {
    return Evap_param::set_tol( left_E, right_E );
  }
  // *** End of the virtual functions ***
};

// **************** Function to integrate *********************
namespace Maxwell_F
{
  // --------------------  Eout_F ------------------
  // Function for the 1-d quadrature over outgoing energy
  //! \param E_out energy of outgoing particle
  //! \param E_out_param parameters for this function
  //! \param value the computed contribution to the transfer matrix
  void Eout_F( double E_out, QuadParamBase *E_out_param, coef_vector *value );
}

#endif
