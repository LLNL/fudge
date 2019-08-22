/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2011-01-18 11:06:56 -0800 (Tue, 18 Jan 2011) $
 * $Author: hedstrom $
 * $Id: Legendre2Body.hpp 1 2011-01-18 11:06:56 -0800 hedstrom $
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
// define the classes used for ENDF Legendre expansions for discrete 2-body reactions

#ifndef LEGENDRE_2BODY_DEF
#define LEGENDRE_2BODY_DEF

// we use the angle_hit_list class
#include "angle_dist.hpp"

class Legendre2d_param;  // forward declaration

//! Class for angular distributions as Legendre expansions
// ---------------- class Legendre_angle_dist ------------------
class Legendre_angle_dist : public list< Legendre_coefs >
{
private:
  double threshold;
  double threshold_out;  // outgoing energy at the threshold
  dd_entry flip;  // Ein for minimal Eout

  two_body_map map;
  relativistic_masses relativistic_mass;

  int first_Ein;  // index of the left-hand end of the first significant energy bin
  int last_Ein;  // index of the right-hand end of the last significant energy bin

  //! Uses the mass difference to set the threshold
  void set_threshold( );

  //! Initializes the quadrature parameters
  //! \param mu_quad_method method for integration over direction cosine
  //! \param Ein_param parameters for quadrature over incident energy
  void initialize_param( Quadrature_Method mu_quad_method, 
     Legendre2d_param *Ein_param );

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

  //! Sets up the map from center-of-mass to laboratory coordinates
  void setup_map( );

  //! Sets up the Ein_param->map from center-of-mass to laboratory coordinates
  //! \param Ein_param parameters for quadrature over incident energy
  void setup_param_map( Legendre2d_param *Ein_param );

  //! Sets the range of incident energies for this intergration
  //! \param Ein_param parameters for quadrature over incident energy
  void set_Ein_range( Legendre2d_param *Ein_param );

  //! Initializes the quadrature parameters
  //! \param Ein_param parameters for quadrature over incident energy
  void setup_data( Legendre2d_param *Ein_param );

  //! Adds to the transfer matrix for all E_out bins for a pair of incident energies.
  //! \param transfer the transfer matrix
  //! \param Ein_param parameters for quadrature over incident energy
  void Eout_ladder( T_matrix& transfer, Legendre2d_param *Ein_param );

  //! Does the integration for one E-E' box between 2 eta = const hyperbolas
  //! \param transfer the transfer matrix
  //! \param Eout_count count of the current outgoing energy bin
  //! \param Ein_param parameters for quadrature over incident energy
  void one_Ebox( T_matrix& transfer, int Eout_count, Legendre2d_param *Ein_param );

  //! Go to the next pair of incident energies.  Returns "true" when finished.
  //! \param E_in the next incident energy
  //! \param Ein_param parameters for quadrature over incident energy
  bool next_ladder( double E_in, Legendre2d_param *Ein_param );

  //! Adds to an element of transfer the integral between the intersections of 2 eta = const hyperbolas with the E-E' box
  //! \param transfer the transfer matrix
  //! \param Eout_count identify the row of the transfer matrix
  //! \param Ein_param parameters for quadrature over incident energy
  void update_T( T_matrix &transfer, int Eout_count, Legendre2d_param *Ein_param );

public:
  bool use_relativistic;
  Interp_Type Ein_interp;  // interpolation between incident energies

  Legendre_angle_dist( ): use_relativistic( false ), Ein_interp( NOTSET ) {}
  ~Legendre_angle_dist( ) {}

  particleInfo particle_info;

  double Q;  // the reaction Q value

  //! Reads the data from Python
  //! \param input_file input file
  //! \param num_Ein number of incident energies for this reaction
  void read_data( data_parser &input_file, int num_Ein );

  //! Calculates the transfer matrix for this particle
  //! \param sigma the cross section data
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param transfer the transfer matrix
  void get_T( const dd_vector& sigma, const dd_vector& weight, T_matrix& transfer );

  // for debugging
  void print( );

};

//! Class for parameters for the 2-d quadrature
// ---------------- class Legendre2d_param ------------------
class Legendre2d_param : public param_base
{
private:
public:
  double threshold;
  double threshold_out;  // outgoing energy at the threshold

  Newton_map_param Newton_map;
  relativistic_param relativistic_map;

  Legendre_angle_dist::const_iterator left_data;        // 4 pointers to (eta, probability) data
  Legendre_angle_dist::const_iterator right_data;

  int num_negative;
  bool use_relativistic;

  // where the lower and upper eta values hit the quadrature box
  angle_hit_list lower_hits;
  angle_hit_list upper_hits;
 
  long int quad_count;  // number of 2-d quadratures
  long int Ein_F_count;  // number of calls to Legendre2Body_F::E_quad_F
  long int mu_F_count;  // number of calls to Legendre2Body_F::mu_cm_quad_F

  // quadrature method for center-of-mass cosine
  Quadrature_Method mu_quad_method;
  Interp_Type Ein_interp;

  inline Legendre2d_param(): num_negative( 0 ), quad_count( 0 ),
			     Ein_F_count( 0 ), mu_F_count( 0 ) {}
  inline ~Legendre2d_param() {}

};

//! Class for parameters for the 1-d quadrature based on a set of Legendre coefficients
// ---------------- class Legendre_param ------------------
class Legendre_param : public QuadParamBase
{
public:
  Legendre_coefs coefs;
  Newton_map_param *Newton_map;
  relativistic_param *relativistic_map;
  int num_negative;
  bool flag_set;
  bool use_relativistic;

  Interp_Type Ein_interp;
  inline Legendre_param( ): num_negative( 0 ), flag_set( false )  {}
  inline ~Legendre_param( ) {}

  //! Interpolates between two incident energies
  void interpolate( double Ein, Legendre_angle_dist::const_iterator prev_coef,
		    Legendre_angle_dist::const_iterator next_coef );
};
 
// **************** functions to integrate **********
namespace Legendre2Body_F
{
  // ---------------- mu_cm_quad_F ------------------
  //! Function for the 1-d quadrature over mu_cm
  //! \param mu_cm center-of-mass direction cosine for outgoing particle
  //! \param mu_cm_quad_param parameters for this function
  //! \param value computed contribution to the transfer matrix
  void mu_cm_quad_F( double mu_cm, QuadParamBase *mu_cm_quad_param, coef_vector *value );

  // ---------------- E_quad_F ------------------
  //! Function for the 2-d quadrature over (E_in, mu_cm )
  //! \param E_in inergy of the incident particle
  //! \param e_quad_param parameters for this function
  //! \param value computed contribution to the transfer matrix
  void E_quad_F( double E_in, QuadParamBase *e_quad_param, coef_vector *value );
}

#endif
