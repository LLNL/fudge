/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: MadlandNix.hpp 1 2006-02-02 03:06:56Z hedstrom $
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
//! Classes used for the Madland-Nix model energy probability density

#ifndef MADLAND_NIX_DEF
#define MADLAND_NIX_DEF

#include "energy_function.hpp"

class MadlandNix_param;

//! Class for the Madland-Nix model
// ---------------- class MadlandNix ------------------
class MadlandNix : public energy_function
{
private:
public:
  double Efl;  // average kinetic energy of light fission fragments
  double Efh;  // average kinetic energy of heavy fission fragments

  double maxEout;  // maximum outgoing energy, independent of incident energy

  MadlandNix( ): Efl( 0.0 ), Efh( 0.0 ), maxEout( 0.0 ) {}
  ~MadlandNix( ) {}

  //!  Gets the range of nontrivial incident energy bins; computes first_Ein and last_Ein
  //! returns true if the threshold is too high for the energy bins
  //! \param sigma the cross section data
  //! \param mult the multiplicity of the outgoing particle
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param e_flux approximate flux used to weight the transfer matrix
  //! \param Ein_groups the boundaries of the incident energy groups
  bool get_Ein_range( const dd_vector& sigma, const dd_vector& mult,
    const dd_vector& weight,
    const Flux_List& e_flux, const Energy_groups& Ein_groups );

  //! Calculates the transfer matrix for this particle
  //! \param sigma the cross section data
  //! \param multiple the multiplicity of the outgoing particle
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param transfer the computed transfer matrix
  void get_T( const dd_vector& sigma, const dd_vector& multiple,
    const dd_vector& weight, T_matrix& transfer );

  // *** Implement the virtual functions ***
  // Initializes the quadrature parameters
  //! \param Eout_groups the boundaries of the outgoing energy groups
  //! \param Ein_param parameters for integration over incident energy
  void setup_data( const Energy_groups& Eout_groups,
			  E_function_param *Ein_param );

  //! Sets the range of incident energies for this intergration
  //! \param Ein_bin identifies the incident energy bin
  //! \param Ein_param parameters for integration over incident energy
  inline void set_Ein_range( int Ein_bin, E_function_param *Ein_param )
  {
    set_Ein_range_default( Ein_bin, Ein_param );
  }

  //! Go to the next pair of incident energies.  Returns "true" when finished.
  //! \param E_in the next incident energy
  //! \param Ein_param parameters for integration over incident energy
  inline bool next_ladder( double E_in, E_function_param *Ein_param )
  {
    return next_ladder_default( E_in, Ein_param );
  }
  // *** End of the virtual functions ***
};

//! Class for parameters for the Madland-Nix model
// ---------------- class MadlandNix_param ------------------
class MadlandNix_param : public E_function_param
{
private:
  //! Parameter alpha = sqrt( TM )
  double alpha;
  //! Parameter beta = sqrt( Ef )
  double beta;

  //! Integration parameters
  double sqrt_E0;  // sqrt( E_0 )
  double sqrt_E1;  // sqrt( E_1 )
  double A;  // ( sqrt_E0 + beta )^2/Theta;
  double B;  // ( sqrt_E1 + beta )^2/Theta;
  double Aprime;  // ( sqrt_E0 - beta )^2/Theta;
  double Bprime;  // ( sqrt_E1 - beta )^2/Theta;

  //! Flags which incomplete gamma function to use
  bool use_tail;

  //! The scale factor for g for the lighter fission fragment
  double scale_g_light;
  //! The scale factor for g for the heavier fission fragment
  double scale_g_heavy;

  //! Interpolate the integral of the probability density
  //! \param E_0 lower limit of the integral over outgoing energy
  //! \param E_1 upper limit of the integral over outgoing energy
  double integral_prob( double E_0, double E_1 );

  //! Interpolate the integral of energy times the probability density
  //! \param E_0 lower limit of the integral over outgoing energy
  //! \param E_1 upper limit of the integral over outgoing energy
  double integral_Eprob( double E_0, double E_1 );

  //! Integral of the g function
  double integral_g( );

  //! Contribution of u_2 to integral_prob
  double integral_u2( );

  //! Contribution of u_2 to integral_Eprob
  double integral_Eu2( );

  //! Contribution of u_1 to integral_prob if E > Ef
  //! See the documentation
  //! \param Bprim the upper limit B' of the integral
  //! \param sqrt_En the square root of the outgoing energy
  double integral_u1_above( double Bprim, double sqrt_En );

  //! Contribution of u_1 to integral_Eprob if E > Ef
  //! See the documentation
  //! \param Bprim the upper limit B' of the integral
  //! \param sqrt_En the square root of the outgoing energy
  double integral_Eu1_above( double Bprim, double sqrt_En );

  //! Contribution of u_1 to integral_prob if E < Ef
  //! See the documentation
  //! \param Aprim the upper limit A' of the integral
  //! \param sqrt_En the square root of the outgoing energy
  double integral_u1_below( double Aprim, double sqrt_En );

  //! Contribution of u_1 to integral_Eprob if E < Ef
  //! See the documentation
  //! \param Aprim the upper limit A' of the integral
  //! \param sqrt_En the square root of the outgoing energy
  double integral_Eu1_below( double Aprim, double sqrt_En );

  //! Calculates the contribution of the light fission fragment
  //! to the integral and returns the noise in the calculation.
  //! See the documentation
  //! \param int_g_light contribution of the light fission fragment
  double get_noise( double *int_g_light );

public:
  //! This model uses Theta as Tm, the maximum excitation level
  //! Square root of average kinetic energy of the lighter fission fragment
  double beta_light;  // sqrt( Efl )
  //! Square root of average kinetic energy of the heavier fission fragment
  double beta_heavy;  // sqrt( Efh )

  inline MadlandNix_param( ) {}
  inline ~MadlandNix_param() {}

  // *** Implement the virtual functions ***
  //! Interpolate the parameters
  //! \param E_in energy of the incident particle
  void set_Ein( double E_in );

  //! Gets the integrals over outgoing energy
  //! \param Eout_0 lower limit of the integral over outgoing energy
  //! \param Eout_1 upper limit of the integral over outgoing energy
  //! \param value the computed integral
  void get_integrals( double Eout_0, double Eout_1, coef_vector &value );

  //! Gets the integrals over outgoing energy and returns the noise in the calculation
  //! \param Eout_0 lower limit of the integral over outgoing energy
  //! \param Eout_1 upper limit of the integral over outgoing energy
  //! \param value the computed integral
  double tol_get_integrals( double Eout_0, double Eout_1, coef_vector &value );

  //! Integrate from 0 to E_max to get the norm
  double get_norm( );

  //! Sets the scale factors for the integrals of probability and energy*probability
  void set_scales( );

  //! Sets the tolerance for the quadrature over incident energy
  //! \param left_E lower limit of the integral over incident energy
  //! \param right_E upper limit of the integral over incident energy
  double set_tol( double left_E, double right_E );
  // *** End of the virtual functions ***
};

#endif
