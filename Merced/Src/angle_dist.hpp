/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: angle_dist.hpp 1 2006-01 19:06:56Z hedstrom $
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
// define the classes used for ENDL i=1 angular distributions

#ifndef ANGLE_DIST_DEF
#define ANGLE_DIST_DEF

#include "mappings.hpp"
#include "relativistic.hpp"
#include "transfer.hpp"
#include "math_util.hpp"
#include "box_geom.hpp"
#include "param_base.hpp"

class two_body_Ein_param;  // forward declaration

// ---------------- class min_Eout_info ------------------
//! This class holds
//! the negative direction cosine, the minimal outgoing energy in the lab
//! frame, and the corresponding incident energy.
class min_Eout_info
{
public:
  double mu;  // direction cosine in center-of-mass frame (negative)
  double Ein;  // incident energy at the minimum
  double Eout; // minimum outgoing energy in lab frame

  min_Eout_info( ) {}
  ~min_Eout_info( ) {}
};

// ---------------- class min_Eout_info_list ------------------
//! This structure holds a list of
//! the negative direction cosine, the minimal outgoing energy in the lab
//! frame, and the corresponding incident energy.
class min_Eout_info_list : public list< min_Eout_info >
{
public:
  min_Eout_info_list( ) {}
  ~min_Eout_info_list( ) {}
};

//! Class for parameters for the 1-d quadrature based on a pair of (mucm, probability) values.
// ---------------- class two_body_mucm_param ------------------
class two_body_mucm_param : public dd_pair, public QuadParamBase
{
public:
  bool use_relativistic;

  Newton_map_param *Newton_map;
  relativistic_param *relativistic_map;

  inline two_body_mucm_param( ): use_relativistic( false ) {}
  inline ~two_body_mucm_param( ) {}
};

//! Class for the list of intersections of a curve with an integration box
//! In this class the hit_list eta parameter is the center-of-mass direction cosine.
// ---------------- class angle_hit_list ------------------
class angle_hit_list : public hit_list
{
private:
  // Implement the virtual functions
  //! Calculates E_out from E_in, for testing the side of a box
  //! \param E_in energy of the incident particle
  double get_Eout( double E_in );

  //! Finds an intersection with the bottom or top of a box
  //! Returns incident energy corresponding to E_out
  //! \param E_out the level to hit
  //! \param pair_0 ( Ein, Eout ) at lower incident energy
  //! \param pair_1 ( Ein, Eout ) at higher incident energy
  double find_hit( double E_out, const dd_entry &pair_0,
		   const dd_entry &pair_1 );

  //! Finds the intersections with the bottom of a box
  //! \param E_out bottom energy of the outgoing energy bin
  //! \param Ein_hits incident energies which give outgoing energy E_out
  void find_bottom_hits( double E_out, vector< Ein_Eta_Hit > *Ein_hits );

  //! Finds the intersections with the top of a box
  //! \param E_out top energy of the outgoing energy bin
  //! \param Ein_hits incident energies which give outgoing energy E_out
  void find_top_hits( double E_out, vector< Ein_Eta_Hit > *Ein_hits );

public:
  //! data for determination of how this outgoing energy curve hits the quadrature box
  dd_entry flip;  // Ein for minimal Eout, minimal Eout for this mucm 
  dd_entry left_Ein_Eout;  // ( Ein, Eout ) for this mucm and for lower Ein
  dd_entry right_Ein_Eout;  // ( Ein, Eout ) for this mucm and for higher Ein
  bool use_relativistic;

  Newton_map_param Newton_map;
  relativistic_param relativistic_map;

  angle_hit_list( ) : use_relativistic( false ) {}
  ~angle_hit_list( ) {}
};

//! Class for angular distributions
// ---------------- class angle_dist ------------------
class angle_dist : public list< dd_vector >
{
private:
  two_body_map map;
  relativistic_masses relativistic_mass;

  min_Eout_info_list min_Eout_list;  // minimal outgoing energies

  int first_Ein;  // index of the left-hand end of the first significant energy bin
  int last_Ein;  // index of the right-hand end of the last significant energy bin

  //! Accounts for the ENDL convention of isotropic emission at low energies
  void ENDL_kludge( );

  //! Returns the common negative cosines for the relativistic treatment of endothermic reactions
  list< double > negative_mu_list( );

  //! Initializes min_Eout_list, used in the treatment of endothermic reactions
  void init_min_Eout_list( );

  //! Computes min_Eout_list, used in the treatment of endothermic reactions
  void make_min_Eout_list( );

  //! Uses the mass difference to set the threshold
  void set_threshold( );

  //! Initializes the quadrature parameters
  //! \param mu_quad_method method of integration over cosine
  //! \param Ein_param the quadrature parameters for integration over incident energy
  void initialize_param( Quadrature_Method mu_quad_method, two_body_Ein_param
    *Ein_param );

  //! Adds to the transfer matrix for all E_out bins for a pair of incident energies.
  //! \param transfer the computed transfer matrix
  //! \param Ein_param the quadrature parameters for integration over incident energy
  void mucm_ladder( T_matrix& transfer, two_body_Ein_param *Ein_param );

  //! Sets up the map from center-of-mass to laboratory coordinates
  void setup_map( );

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

  //! Initializes the quadrature parameters
  //! \param Ein_param the quadrature parameters for integration over incident energy
  void setup_data( two_body_Ein_param *Ein_param );

  //! Sets the range of incident energies for this intergration
  //! \param Ein_bin the number of this incident energy bin
  //! \param Ein_param the quadrature parameters for integration over incident energy
  void set_Ein_range( int Ein_bin, two_body_Ein_param &Ein_param );

  //! Does the integration for one E-E' box between 2 mucm = const hyperbolas
  //! \param transfer the transfer matrix
  //! \param Eout_count count of the current outgoing energy bin
  //! \param Ein_param the quadrature parameters for integration over incident energy
  void one_Ebox( T_matrix& transfer, int Eout_count, two_body_Ein_param *Ein_param );

  //! Go to the next pair of incident energies.  Returns "true" when finished.
  //! \param E_in the next incident energy
  //! \param Ein_param the quadrature parameters for integration over incident energy
  bool next_ladder( double E_in, two_body_Ein_param *Ein_param );

  //! Adds to an element of transfer the integral between the intersections of 2 mucm = const hyperbolas with the E-E' box
  //! \param transfer the transfer matrix
  //! \param Eout_count count of the current outgoing energy bin
  //! \param Ein_param the quadrature parameters for integration over incident energy
  void update_T( T_matrix &transfer, int Eout_count, two_body_Ein_param *Ein_param );

public:

  double threshold;  // the computed threshold
  double threshold_out;  // outgoing energy at the computed threshold
  two_d_interp Ein_interp;  // interpolation between incident energies
  Interp_qualifier Ein_qualifier;  // for interpolation between incident energies
  Interp_Type mu_interp;

  bool use_relativistic;

  //! Constructor
  inline angle_dist( ) : threshold( -1.0 ),
    	 mu_interp( LINLIN ), use_relativistic( false ) {}

  //! Destructor
  inline ~angle_dist( ) {}

  particleInfo particle_info;

  double Q;  // the reaction Q value

  //! Reads the data from Python
  //! \param input_file input file
  //! \param num_Ein number of incident energies for this reaction
  void read_data( data_parser &input_file, int num_Ein );

  //! Calculates the transfer matrix for this particle
  //! \param sigma the cross section data
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param transfer the computed transfer matrix
  void get_T( const dd_vector& sigma, dd_vector& weight, T_matrix& transfer );

  //! Checks whether an angular probability density is isotropic
  bool isotropic( );

  // for debugging
  void print( );

};

//! Class for parameters for the 2-d quadrature
// ---------------- class two_body_Ein_param ------------------
class two_body_Ein_param : public param_base
{
private:
  //! Interpolates (mu_cm, probability) data to the lower common mu_cm value
  //! \param lower_mu the lower end of the mu_cm range
  void common_low_mucm( double lower_mu );

  //! Interpolates (mu_cm, probability) data to the higher common mu_cm value
  //! \param higher_mu the upper end of the mu_cm range
  void common_high_mucm( double higher_mu );

  //! Returns the min_Eout_info for given center-of-mass outgoing cosine
  //! \param mu_cm the center-of-mass outgoing cosine
  min_Eout_info_list::const_iterator get_min_Eout_info( double mu_cm ) const;

  //! Initializes upper_hits for the value of mucm
  //! \param mucm the center-of-mass direction cosine
  void set_upper_hits( double mucm );

public:
  bool use_relativistic;

  Newton_map_param Newton_map;
  relativistic_param relativistic_map;

  double threshold;  // the computed threshold
  double threshold_out;  // outgoing energy at the computed threshold

  double left_data_Ein;                       // incident energy for left data
  double right_data_Ein;                      // incident energy for right data
  // pointers to the current pair of angular distirbutions
  angle_dist::const_iterator this_mucm_dist;
  angle_dist::const_iterator next_mucm_dist;
 
  // (mu_cm, probability) data interpolated to common mu_cm values
  dd_entry left_data;        // (lower mucm, probability) at low Ein
  dd_entry next_left_data;   // (higher mucm, probability) at low Ein
  dd_entry right_data;       // (lower mucm, probability) at high Ein
  dd_entry next_right_data;  // (higher mucm, probability) at high Ein
  long int quad_count;  // number of 2-d quadratures
  long int Ein_F_count;  // number of calls to angle_dist_F::E_quad_F
  long int mu_F_count;  // number of calls to angle_dist_F::mu_cm_quad_F
  Quadrature_Method mu_quad_method;  // method of integration over cosine

  // where the lower and upper mucm values hit the quadrature box
  angle_hit_list lower_hits;
  angle_hit_list upper_hits;

  // pointers to the (mu_cm, probability) data
  dd_vector::const_iterator left_ptr;
  dd_vector::const_iterator next_left_ptr;
  dd_vector::const_iterator right_ptr;
  dd_vector::const_iterator next_right_ptr;

  // the range of outgoing energies for this data
  double upper_data_max_Eout;
  double lower_data_min_Eout;

  min_Eout_info_list *min_Eout_list;  // minimal outgoing energies

  inline two_body_Ein_param(): use_relativistic( false ), quad_count( 0 ),
      Ein_F_count( 0 ), mu_F_count( 0 ) {}
  inline ~two_body_Ein_param() {}

  //! Copies the relativistic mapping from center-of-mass to lab frame
  //! \param map_ the mapping data to copy
  void set_rel_map( relativistic_masses *map );

  //! Copies the Newtonian mapping from center-of-mass to lab frame
  //! \param map_ the mapping data to copy
  void set_Newton_map( two_body_map *map_ );

  //! Initializes the pointers to the angular probabilities for this E_in
  void reset_start( );

  //! Sets up the next interval of mu_cm values
  //! returns true if there is no more data
  bool next_mucm( );

  //! Initializes lower_hits and upper_hits for the values of mucm
  void reset_hits( );
};

// **************** functions to integrate **********
namespace angle_dist_F
{
  // ---------------- mu_cm_quad_F ------------------
  //! Function for the 1-d quadrature
  //! \param mucm center-of-mass direction cosine for outgoing particle
  //! \param mucm_quad_param parameters for this function
  //! \param value computed contribution to the transfer matrix
  void mu_cm_quad_F( double mucm, QuadParamBase *mucm_quad_param, coef_vector *value );

  // ---------------- E_quad_F ------------------
  //! Function for the 2-d quadrature
  //! \param E_in inergy of the incident particle
  //! \param mucm_quad_param parameters for this function
  //! \param value computed contribution to the transfer matrix
  void E_quad_F( double E_in, QuadParamBase *e_quad_param, coef_vector *value );
}

#endif
