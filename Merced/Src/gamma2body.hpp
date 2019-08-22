/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2013-12-10 11:06:56 -0800 (Tue, 10 Dec 2013) $
 * $Author: hedstrom $
 * $Id: gamma2body.hpp 1 2013-12-10 11:06:56 -0800 hedstrom $
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
// define the classes used for Legendre expansions for gammas from capture reactions

#ifndef GAMMA_2BODY_DEF
#define GAMMA_2BODY_DEF

#include "gamma_2_body_map.hpp"
#include "math_util.hpp"
#include "dd_vector.hpp"
#include "box_geom.hpp"
#include "transfer.hpp"

class capture_gamma_Ein_param;  // forward declaration

//! Class to determine the energy of the backward and forward emitted gamma
// ---------------- class red_blue_shift_map ------------------
class red_blue_shift_map : public gamma_2_body_map
{
public:
  //! As a function of incident energy, the energy of the backward gamma 
  //! is decreasing at low incident energy and then increases.
  //! The turn-over energy
  double turnover_T;

  //! The minimum energy of the backward emitted gamma
  double min_E_back;

  //! Default constructor
  inline red_blue_shift_map( ) {}

  //! Default destructor
  inline ~red_blue_shift_map() {}

  //! Gets the turn-over for backward emission
  //! \param below ( T_in_lab, d_red_shift ) below the root
  //! \param above ( T_in_lab, d_red_shift ) above the root
  void get_turnover( dd_entry below, dd_entry above );
};

//! Class for the list of intersections of a curve with an integration box
// ---------------- class gamma_hit_list ------------------
class gamma_hit_list : public hit_list_base
{
private:
  //! Lab-frame gamma energy at the lower incident energy
  dd_entry left_Ein_Egamma;

  //! Lab-frame gamma energy at the higher incident energy
  dd_entry right_Ein_Egamma;

// --------- implement the virtual functions -------------------
  //! Finds the intersections with the bottom of a box
  //! \param lab_E_out the lower desired outgoing energy
  //! \param Ein_hits output: relations between the gamma energy and lab_E_out
  void find_bottom_hits( double lab_E_out, vector< Ein_Eta_Hit > *Ein_hits );

  //! Finds the intersections with the top of a box
  //! \param lab_E_out the higher desired outgoing energy
  //! \param Ein_hits output: relations between the gamma energy and lab_E_out
  void find_top_hits( double lab_E_out, vector< Ein_Eta_Hit > *Ein_hits );

  //! Where do we hit the left-hand side of the box?
  //! \param E_in lower incident energy
  //! \param Eout_ptr desired lower energy of the outgoing particle
  void test_left( double E_in, vector< double >::const_iterator Eout_ptr )
  { test_left_default( E_in, Eout_ptr ); }

  //! Where do we hit the right-hand side of the box?
  //! \param E_in higher incident energy
  //! \param Eout_ptr desired lower energy of the outgoing particle
  void test_right( double E_in, vector< double >::const_iterator Eout_ptr )
  { test_right_default( E_in, Eout_ptr ); }

  //! The code appends the previously calculated left_hit and right_hit
  //! \param Eout_ptr desired lower energy of the outgoing particle
  //! \param E_in_left lower incident energy
  //! \param E_in_right higher incident energy
  void test_sides( vector< double >::const_iterator Eout_ptr,
		   double E_in_left, double E_in_right );

  //! Checks for inconsistencies caused by peculiarities of real arithmetic.
  //! Returns is_OK == true if there are no inconsistencies.
  inline bool consistency( )
  { return consistency_default( ); }
  // --------- end of virtual functions -------------------

  //! Finds the incident energy producing a backward gamma with energy lab_E_out.
  //! Use this routine when the lab-frame gamma energy decreases in incident energy.
  //! \param lab_E_out target energy for the outgoing backward gamma
  //! \param low_Ein pair ( lower incident energy, lab-frame gamma energy )
  //! \param high_Ein pair ( higher incident energy, lab-frame gamma energy )
  double find_decreasing_hit( double lab_E_out,
    const dd_entry &low_Ein, const dd_entry &high_Ein );

  //! Finds the incident energy producing a gamma with energy lab_E_out
  //! Use this routine when the lab-frame gamma energy increases in incident energy.
  //! \param lab_E_out target energy for the outgoing backward gamma
  //! \param low_Ein pair ( lower incident energy, lab-frame gamma energy )
  //! \param high_Ein pair ( higher incident energy, lab-frame gamma energy )
  double find_increasing_hit( double lab_E_out,
    const dd_entry &low_Ein, const dd_entry &high_Ein );

  //! Finds the relation between the outgoing gamma energy and the value lab_Eout
  //! \param lab_Eout the desired lab-frame energy of the outgoing gamma
  //! \param Ein_hits the output pairs ( incident energy, relation to lab_Eout )
  void test_lab_Eout( double lab_Eout, vector< Ein_Eta_Hit > *Ein_hits );

public:
  red_blue_shift_map *map;

  gamma_hit_list( ) : map( 0 ) {}
  ~gamma_hit_list( ) {}

  //! Sets left_Ein_Egamma and right_Ein_Egamma.
  //! \param mu_cm direction cosine of gamma in the center-of-mass frame
  //! \param E_in_left lower end of the incident energy range
  //! \param E_in_right upper end of the incident energy range
  void set_E_range( double mu_cm, double E_in_left, double E_in_right );
};

// ---------------- class capture_gamma ------------------
class capture_gamma : public list< Legendre_coefs >
{
private:
  int first_Ein;  // index of the left-hand end of the first significant energy bin
  int last_Ein;  // index of the right-hand end of the last significant energy bin

  red_blue_shift_map map;
 
  //! Finds the minimum energy of backward gammas
  //! \param transfer the computed transfer matrix
  void get_min_E_back( const T_matrix& transfer );

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

  //! Initializes the quadrature parameters; returns true if the threshold is too high
  //! \param Ein_param the quadrature parameters for integration over incident energy
  void setup_data( capture_gamma_Ein_param *Ein_param );

  //! Sets the range of incident energies for this intergration
  //! \param Ein_param the quadrature parameters for integration over incident energy
  void set_Ein_range( capture_gamma_Ein_param *Ein_param );

  //! Adds to the transfer matrix for all E_out bins for a pair of incident energies.
  //! \param transfer the computed transfer matrix
  //! \param Ein_param the quadrature parameters for integration over incident energy
  void Eout_ladder( T_matrix& transfer, capture_gamma_Ein_param *Ein_param );

  //! Does the integration for one E-E' box between 2 eta = const hyperbolas
  //! \param transfer the computed transfer matrix
  //! \param Eout_count which column of the transfer matrix
  //! \param Ein_param the quadrature parameters for integration over incident energy
  void one_Ebox( T_matrix& transfer, int Eout_count,
     capture_gamma_Ein_param *Ein_param );

  //! Go to the next pair of incident energies.  Returns "true" when finished.
  //! \param E_in the energy of the incident particle
  //! \param Ein_param the quadrature parameters for integration over incident energy
  bool next_ladder( double E_in, capture_gamma_Ein_param *Ein_param );

  //! Adds to an element of transfer the integral between the intersections of 2 eta = const hyperbolas with the E-E' box
  //! \param transfer the computed transfer matrix
  //! \param Eout_count which column of the transfer matrix
  //! \param Ein_param the quadrature parameters for integration over incident energy
  void update_T( T_matrix &transfer, int Eout_count,
    capture_gamma_Ein_param *Ein_param );

public:

  Interp_Type Ein_interp;  // interpolation between incident energies

  capture_gamma( ): Ein_interp( NOTSET ) {}
  ~capture_gamma( ) {}

  particleInfo particle_info;

  double Q;  // the reaction Q value

  //! Reads the data from Python
  //! \param input_file the input data file
  //! \param num_Ein the number of incident energies in the data
  void read_data( data_parser &input_file, int num_Ein );

  //! Calculates the transfer matrix for this particle
  //! \param sigma the cross section data
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param transfer the computed transfer matrix
  void get_T( const dd_vector& sigma, const dd_vector& weight, T_matrix& transfer );

  // for debugging
  void print( );

};

//! Class for parameters for the 2-d quadrature
// ---------------- class capture_gamma_Ein_param ------------------
class capture_gamma_Ein_param : public param_base
{
public:
  // where the forward and backward gamma energies hit the quadrature box
  gamma_hit_list backward_hits;
  gamma_hit_list forward_hits;

  red_blue_shift_map map;
  capture_gamma::const_iterator left_data;        // 4 pointers to (eta, probability) data
  capture_gamma::const_iterator right_data;

  long int quad_count;  // number of 2-d quadratures
  long int Ein_F_count;  // number of calls to Capture_Gamma_F::Ein_F
  long int mu_F_count;  // number of calls to Capture_Gamma_F::mu_F
  int num_negative;  // number of negative sums of Legendre series

  Quadrature_Method mu_quad_method;  // quadrature method for outgoing cosine
  Interp_Type Ein_interp;

  inline capture_gamma_Ein_param(): quad_count( 0 ), Ein_F_count( 0 ),
				    mu_F_count( 0 ), num_negative( 0 ) {}
  inline ~capture_gamma_Ein_param() {}

  //! Sets up the mapping from center-of-mass to lab frame
  //! \param map_ information for the map from center-of-mass to lab frame
  void set_map( red_blue_shift_map &map_ );

};

//! Class for parameters for the 1-d quadrature based on a set of Legendre coefficients
// ---------------- class capture_gamma_mu_param ------------------
class capture_gamma_mu_param : public QuadParamBase
{
public:
  Legendre_coefs coefs;
  red_blue_shift_map *map;
  int num_negative;
  bool flag_set;
  Interp_Type Ein_interp;
  inline capture_gamma_mu_param( ): num_negative( 0 ), flag_set( false )  {}
  inline ~capture_gamma_mu_param( ) {}

  //! Interpolates between two incident energies
  void interpolate( double Ein, capture_gamma::const_iterator prev_coef,
                    capture_gamma::const_iterator next_coef );
};

// ************* functions *************
namespace Capture_Gamma_F
{
  // ------------------ red_shift -----------------------
  //! Returns the energy in lab frame of the backward emitted gamma.
  //! \param T_in_lab kinetic energy of the incident particle in the lab frame
  //! \param params the function parameters
  double red_shift( double T_in_lab, void *params );

  // ------------------ d_red_shift -----------------------
  //! Returns the derivative of the energy of the backward emitted gamma in the
  //! lab frame with respect to the kinetic energy of the incident particle.
  //! \param T_in_lab kinetic energy of the incident particle in the lab frame
  //! \param params the function parameters
  double d_red_shift( double T_in_lab, void *params );

  // ------------------ blue_shift -----------------------
  //! Returns the energy in lab frame of the forward emitted gamma.
  //! \param T_in_lab kinetic energy of the incident particle in the lab frame
  //! \param params the function parameters
  double blue_shift( double T_in_lab, void *params );

  // ---------------- mu_F ------------------
  //! Function for the 1-d quadrature over mu_cm
  //! \param mu_cm center-of-mass direction cosine for outgoing particle
  //! \param mu_cm_quad_param parameters for this function
  //! \param value computed contribution to the transfer matrix
  void mu_F( double mu_cm, QuadParamBase *mu_cm_quad_param, coef_vector *value );

  // ---------------- Ein_F ------------------
  //! Function for the 2-d quadrature over (E_in, mu_cm )
  //! \param E_in inergy of the incident particle
  //! \param eta_quad_param parameters for this function
  //! \param value computed contribution to the transfer matrix
  void Ein_F( double E_in, QuadParamBase *e_quad_param, coef_vector *value );
}

#endif
