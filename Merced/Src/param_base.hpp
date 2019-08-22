/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id:param_base.hpp 1 2006-02-02 03:06:56Z hedstrom $
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

//header for param_base, the base class for the quadrature parameters

#ifndef PARAM_BASE_CLASS
#define PARAM_BASE_CLASS

#include "dd_vector.hpp"
#include "Legendre_data.hpp"
#include "Eout_integrals.hpp"
#include "Energy_groups.hpp"

// ----------- class QuadParamBase ------------------
//! Base lass for the quadrature parameters
class QuadParamBase
{
public:
  int func_count;      // the number of function calls

  //! Default constructor
  inline QuadParamBase(): func_count( 0 )
   {}

  //! Default destructor
  inline ~QuadParamBase() {}

};
// ----------- class param_base ------------------
//! Base lass for the parameters for integration over incident energy
class param_base: public QuadParamBase
{
public:
  double data_E_0;     // the lower incident energy for the eta ladder
  double data_E_1;     // the upper incident energy for the eta ladder
  double Ein_0;        // the lower incident energy for quadrature
  double Ein_1;        // the upper incident energy for quadrature
  double Eout_min;     // the bottom of the E-E' quadrature box
  double Eout_max;     // the top of the E-E' quadrature box
  bool use_Eout_min;   // if true, get minimum eta from Eout_min
                       // otherwise get minimum eta from the lower data point
  bool use_Eout_max;   // if true, get maximum eta from Eout_max
                       // otherwise get maximum eta from the upper data point
  
  // pointers to the cross section
  dd_vector::const_iterator this_sigma;
  dd_vector::const_iterator next_sigma;
  dd_vector::const_iterator first_ladder_sigma;  // first sigma for this eta ladder
  dd_vector::const_iterator last_ladder_sigma;   // last sigma for this eta ladder
  dd_vector::const_iterator sigma_end;           // final sigma

  // pointers to the multiplicity
  dd_vector::const_iterator this_mult;
  dd_vector::const_iterator next_mult;
  dd_vector::const_iterator mult_end;

  // pointers to the model weights
  dd_vector::const_iterator this_weight;
  dd_vector::const_iterator next_weight;
  dd_vector::const_iterator weight_end;

  // pointers to the flux data
  Flux_List::const_iterator flux_ptr;
  Flux_List::const_iterator next_flux;
  Flux_List::const_iterator flux_end;

  // pointers to the incident energy boundaries
  Energy_groups::const_iterator Ein_ptr;
  Energy_groups::const_iterator next_Ein;
  Energy_groups::const_iterator Ein_end;
  int Ein_count;

  // pointers to the integrals over E_out
  Eout_integrals::const_iterator this_Eout_int;
  Eout_integrals::const_iterator next_Eout_int;
  Eout_integrals::const_iterator Eout_int_end;

  //! Holds the weight: (cross section) * flux * multiplicity * model weight
  Legendre_coefs current_weight;

  int order;             // The Legendre order of the problem

  //! Default constructor
  inline param_base()
   {}

  //! Default destructor
  inline ~param_base()
    {}

  //!  Gets the range of nontrivial incident energy bins; computes first_Ein and last_Ein
  //! returns true if the threshold is too high for the energy bins
  //! \param sigma_ the cross section data
  //! \param mult_ the outgoing particle multiplicity data
  //! \param weight_ the weighting to apply to the transfer matrix entries
  //! \param e_flux_ the initial approximation to apply to the particle flux
  //! \param Ein_groups the boundaried of the incident energy groups
  //! \param first_Ein the first common incident energy bin
  //! \param last_Ein the last common incident energy bin
  bool get_Ein_range( const dd_vector& sigma_, const dd_vector& mult_,
    const dd_vector& weight_,
    const Flux_List& e_flux_, const Energy_groups& Ein_groups,
    double *first_Ein, double *last_Ein );

  //!  Sets up the initial quadrature parameters
  //! \param sigma_ the cross section data
  //! \param mult_ the outgoing particle multiplicity data
  //! \param weight_ the weighting to apply to the transfer matrix entries
  //! \param e_flux_ the initial approximation to apply to the particle flux
  //! \param Ein_groups the boundaries of the incident energy groups
   void setup( const dd_vector& sigma_, const dd_vector& mult_,
    const dd_vector& weight_,
    const Flux_List& e_flux_, const Energy_groups& Ein_groups );

  //! Sets the first common incident energy
  void common_E0( );

  //!  Sets up the initial quadrature parameters for parallel computing by bin
  //! \param Ein_bin_ the number of the incident energy bin
  //! \param sigma_ the cross section data
  //! \param mult_ the outgoing particle multiplicity data
  //! \param weight_ the weighting to apply to the transfer matrix entries
  //! \param e_flux_ the initial approximation to apply to the particle flux
  //! \param Ein_groups the boundaries of the incident energy groups
  void setup_bin( int Ein_bin_, const dd_vector& sigma_, const dd_vector& mult_,
    const dd_vector& weight_,
    const Flux_List& e_flux_, const Energy_groups& Ein_groups );

  //! Sets the range of integration over incident energy
  void set_Ein_range( );

  //! Sets the data pointers for a new incident energy interval.
  //! Returns "true" if we are finished with the data.
  //! \param E_in the next incident energy
  bool update_pointers( double E_in );

  //! Sets the data pointers for a new incident energy interval in one bin.
  //! Returns "true" if we are finished with the data on this bin.
  //! \param E_in the next incident energy
  bool update_bin_pointers( double E_in );

  //! Sets the range of pointers to cross sections for this set of data
  void set_sigma_range( );

  //! Calculates the weight: (cross section) * flux * multiplicity * model weight
  //! \param E_in the incident energy
  void set_weight( double E_in );

  //! Calculates the weight for gammas: flux * multiplicity * model weight
  //! \param E_in the incident energy
  void flux_weight( double E_in );

};

#endif
