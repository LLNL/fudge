/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2014-04-22 (Tue, Apr 22, 2014) $
 * $Author: hedstrom $
 * $Id: joint_dist.hpp 1 2014-04-22 hedstrom $
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
//header for the classes used on joint energy-angle distributions

#ifndef JOINT_DIST_CLASS
#define JOINT_DIST_CLASS

#include "param_base.hpp"
#include "math_util.hpp"
#include "transfer.hpp"
#include "angle_dist.hpp"
#include "box_geom.hpp"

class joint_dist_param;  // forward declaration

// ----------- class E_mu_P_data  -----------------
//! Class for one set of energy-mu-probability data
class E_mu_P_data
{
public:
  double UB_mu;   // the unit-base direction cosine of the outgoing particle
  double phys_mu;   // the physical direction cosine of the outgoing particle
  double UB_Eout; // the unit-base energy of the outgoing particle
  double phys_Eout; // the physical energy of the outgoing particle
  double Prob;  // the probability density of the outgoing particle

  //! Default constructor
  inline E_mu_P_data( )
  {}

  //! Default destructor
  inline ~E_mu_P_data( )
  {}

  //! Interpolates with respect to unit-base outgoing energy
  //! \param mid_UB_Eout an intermediate unit-base outgoing energy
  //! \param next_data E_mu_P_data at a higher unit-base outgoing energy
  //! \param mid_data E_mu_P_data at unit-base outgoing energy mid_UB_Eout
  void UB_Eout_interp( double mid_UB_Eout, const E_mu_P_data &next_data,
		       E_mu_P_data *mid_data ) const;

  //! Interpolates histogram data with respect to unit-base outgoing energy
  //! \param mid_UB_Eout an intermediate unit-base outgoing energy
  //! \param next_data E_mu_P_data at a higher unit-base outgoing energy
  //! \param mid_data E_mu_P_data at unit-base outgoing energy mid_UB_Eout
  void UB_Eout_histogram( double mid_UB_Eout, const E_mu_P_data &next_data,
		       E_mu_P_data *mid_data ) const;

  //! Interpolates with respect to direction cosine
  //! \param mid_mu an intermediate direction cosine
  //! \param next_data E_mu_P_data at a higher direction cosine
  //! \param mid_data E_mu_P_data at direction cosine
  void mu_interp( double mid_mu, const E_mu_P_data &next_data,
		  E_mu_P_data *mid_data ) const;

  //! Interpolates with respect to incident energy
  //! \param alpha the weight for next_data
  //! \param next_data E_mu_P_data at a higher incident energy
  //! \param mid_data E_mu_P_data at an an intermediate incident energy
  void Ein_interp( double alpha, const E_mu_P_data &next_data,
		   E_mu_P_data *mid_data ) const;

  //! Copies the data
  //! \param to_copy E_mu_P_data to copy
  void copy( const E_mu_P_data &to_copy );
};

// ----------- class joint_dist_hits -----------------
//! Class for the list of intersections of a unit-base curve with an integration box
class joint_dist_hits : public energy_hit_list
{
private:

public:
  // ! where we are in the data
  double Ein_0;  // lower incident energy
  double Ein_1;  // higher incident energy
  E_mu_P_data *Ein0_data;  // lower Ein
  E_mu_P_data *Ein1_data;;  // higher Ein

  //! Constructor
  inline joint_dist_hits( )
  {}

  //! Default destructor
  inline ~joint_dist_hits( )
  {}

  //! Sets the range of integration over incident energy
  //! \param A the lower incident energy for this integral
  //! \param B the higher incident energy for this integral
  void set_incident_range( double A, double B );

  //! Gets the physical outgoing energies at the ends of the incident energy bin
  void get_phys_Eout( );
};

// ----------- class current_data -----------------
//! class for E_mu_P_data at current values of Ein, mu, and Eout
class current_data
{
private:
  double E_in;  // energy of the incident particle

public:
  //! ( mu, outgoing energy, probability density ) data
  E_mu_P_data mu0_Eout0;  // lower mu, lower Eout
  E_mu_P_data mu0_Eout1;  // lower mu, higher Eout
  E_mu_P_data mu1_Eout0;  // higher mu, lower Eout
  E_mu_P_data mu1_Eout1;  // higher mu, higher Eout

  //! range of outgoing energies for the unit-base maps
  unit_base_map mu0_ubase_map;  // lower mu
  unit_base_map mu1_ubase_map;  // higher mu

  //! Constructor
  inline current_data( )
  {}

  //! Default destructor
  inline ~current_data( )
  {}

  //! Sets the incident energy
  //! \param Ein the energy of the incident particle
  inline void set_E_in( double Ein )
  { E_in = Ein; }

  //! Returns the energy of the incident particle
  inline double get_E_in( ) const
  { return E_in; }

  //! Interpolate in incident energy, returns mid_data
  //! \param mid_Ein an intermediate incident energy
  //! \param next_data (E_out, probability) at a higher incident energy
  //! \param mid_data (E_out, probability) at incident energy mid_Ein
  void Ein_interpolate( double mid_Ein, const current_data &next_data,
		    current_data *mid_data ) const;

  //! Interpolate in mu; returns Eout0_data, Eout1_data, and mid_ubase_map
  //! \param mid_mu an intermediate direction cosine
  //! \param Eout0_data (E_out, probability) at mid_mu and lower outgoing energy
  //! \param Eout1_data (E_out, probability) at mid_mu and higher outgoing energy
  //! \param mid_ubase_map the outgoing energy range at mid_mu
  void mu_interpolate( double mid_mu, E_mu_P_data *Eout0_data,
		       E_mu_P_data *Eout1_data, unit_base_map *mid_ubase_map ) const;

};

// ----------- class joint_Eout_param -----------------
//! parameters for 1-d integration in outgoing energy
class joint_Eout_param: public param_base
{
public:
  //! Rule for interpolation in outgoing energy
  Interp_Type Eout_interp;

  //! (Eout, probability) data at the lower outgoing energy
  E_mu_P_data Eout0_data;

  //! (Eout, probability) data at the higher outgoing energy
  E_mu_P_data Eout1_data;

  //! Default constructor
  inline joint_Eout_param( )
  {}

  //! Default destructor
  inline ~joint_Eout_param( )
  {}

};

// ----------- class joint_mu_param -----------------
//! parameters for 2-d integration over mu and E_out
class joint_mu_param: public param_base
{
private:
 
public:
  //! (Eout, probability) data at this incident energy
  current_data this_data;

  //! parameters for integration over outgoing energy
  joint_Eout_param Eout_params;

  //! unit-base maps of the outgoing energy
  unit_base_map mu0_ubase_map;  // outgoing energy range at lower_mu
  unit_base_map mu1_ubase_map;  // outgoing energy range at upper_mu
  unit_base_map mid_ubase_map;  // outgoing energy range at the current mu

  Quadrature_Method Eout_quad_method;  // quadrature method for outgoing energy
  long int Eout_F_count;  // number of calls to joint_dist_F::Eout_F

  //! the outgoing energy bin boundaries
  vector< double >::const_iterator Eout_bottom;
  vector< double >::const_iterator Eout_top;

  //! Default constructor
  inline joint_mu_param( ): Eout_F_count ( 0 ) 
  {}

  //! Default destructor
  inline ~joint_mu_param( )
  {}

  //! determines the geometry for 2d integration over outgoing cosine and energy
  //! returns true if the geometry makes sense
  //! \param this_data (Eout, probability) data at the current incident energy
  bool geometry( const current_data &this_data );
};

// ----------- class one_mu -----------------
//! Class for the energy distribution at given incident energy and angle
class one_mu : public dd_vector
{
private:

public:
  unit_base_map ubase_map;

  //! Default constructor
  inline one_mu( )
  {}

  //! Default destructor
  inline ~one_mu( )
  {}

  //! Use the direction cosine, mu, as the tag
  inline double get_mu( ) const
  {
    return get_tag( );
  }

  //! Sets the direction cosine tag
  //! \param mu the direction cosine of the outgoing particle
  inline void set_mu( double mu )
  {
    set_tag( mu );
  }

  //! We need to make a copy for the one_joint_dist_quad::copy routine
  //! \param to_copy the vector of pairs ( outgoing energy, probability ) to copy
  void copy( const one_mu &to_copy );

  //! Sets up E_mu_P_data at unit-base outgoing energy UB_Eout
  //! \param UB_Eout an intermediate unit-base outgoing energy
  //! \param phys_mu the physical direction cosine
  //! \param prev_data ( unit-base outgoing energy, probability ) at lower E_out
  //! \param next_data ( unit-base outgoing energy, probability ) at higher E_out
  //! \param mid_data E_mu_P_data at unit-base outgoing energy UB_Eout
  void set_E_mu_P_data( double UB_Eout, double phys_mu,
    one_mu::const_iterator prev_data,
    one_mu::const_iterator next_data, E_mu_P_data *mid_data );
};

// ----------- class one_joint_dist -----------------
//! Class for the joint energy-angle distribution at given incident energy
class one_joint_dist : public list< one_mu >
{
private:
  double tag_;

public:
  //! for mapping direction cosines to 0 <= mu <= 1
  unit_base_map mu_ubase_map;

  two_d_interp Ein_interp;
  Interp_Type Eout_interp;
  two_d_interp mu_interp;

  //! Default constructor
  inline one_joint_dist( )
  {}

  //! Default destructor
  inline ~one_joint_dist( )
  {}

  //! Use E_in as the tag
  inline double get_E_in( ) const
  {
    return tag_;
  }

  //! Sets the tag
  //! \param E_in the energy of the incident particle
  inline void set_E_in( double E_in )
  {
    tag_ = E_in;
  }

  //! Makes a copy at the threshold, used by the ENDL_kludge routine
  //! \param to_copy the vectors of pairs ( outgoing energy, probability ) to copy
  void copy( const one_joint_dist &to_copy );

  //! Converts double-differential data from ENDL to ENDF format
  //! \param angles vector of pairs ( outgoing direction cosine, probability )
  void to_ENDF( angle_dist::const_iterator &angles );

  //! Maps the direction cosines to 0 <= mu <= 1
  void to_unit_base( );

  //! Insert energy probability densities for an intermediate mu
  //! This routine is used when ENDL data has (mu, probability) data
  //! at a mu value which is missing from the (mu, Eout, probability) table
  //! \param mu the intermediate direction cosine
  //! \param prev_joint_mu (Eout, probability) data at a lower mu value
  //! \param joint_mu (Eout, probability) data at a higher mu value
  void interpolate_mu( double mu,
    one_joint_dist::iterator &prev_joint_mu,
    one_joint_dist::iterator &joint_mu );
};

// ----------- class joint_dist -----------------
//! Class for joint energy-angle distributions
class joint_dist : public list< one_joint_dist >
{
private:
  int first_Ein;  // index of the left-hand end of the first significant energy bin
  int last_Ein;  // index of the right-hand end of the last significant energy bin

  //! ENDL has the convention that the distribution is assumed independent of energy at low incident energies.
  void ENDL_kludge( );

  // Converts double-differential data from ENDL to ENDF format
  void to_ENDF( );

  // Maps the direction cosines to 0 <= mu <= 1
  void to_unit_base( );

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
  //! \param Ein_param the qudrature parameters
  void setup_param( joint_dist_param *Ein_param );

  //! Sets the range of incident energies for this intergration
  //! \param Ein_param the qudrature parameters
  void set_Ein_range( joint_dist_param *Ein_param );

  //! Handles the ( cosine, Eout, probability ) data for one pair of incident energies
  //! \param transfer the transfer matrix
  //! \param Ein_param the qudrature parameters
  void mu_data_ladder( T_matrix& transfer, joint_dist_param *Ein_param );

  //! Go to the next pair of incident energies.  Returns "true" when finished.
  //! \param E_in energy of the incident particle
  //! \param Ein_param the qudrature parameters
  bool next_Ein_pair( double E_in, joint_dist_param *Ein_param );

  //! Initializes the pointers to the ( Eout, probability ) data for current cosines
  //! \param Ein_param the qudrature parameters
  void start_mu_data( joint_dist_param *Ein_param );

  //! Go to the next pairs of direction cosines.  Returns "true" when finished.
  //! \param Ein_param the qudrature parameters
  bool next_mu_pairs( joint_dist_param *Ein_param );

  //! Adds to the transfer matrix for all E_out bins for given Ein and mu ranges
  //! \param transfer the transfer matrix
  //! \param Ein_param the qudrature parameters
  void Eout_data_ladder( T_matrix& transfer, joint_dist_param *Ein_param );

  //! Starts one staircase of the Eout data
  //! \param Ein_param the qudrature parameters
  void start_Eout_data( joint_dist_param *Ein_param );

  //! go to the next sets of (E_out, probability) pairs
  //! \param Ein_param the qudrature parameters
  bool next_Eout( joint_dist_param *Ein_param );

  //! Adds to the transfer matrix for the current data
  //! \param transfer the transfer matrix
  //! \param Eout_count the current row of the transfer matrix
  //! \param Ein_param the qudrature parameters
  void one_box( T_matrix& transfer, int Eout_count, joint_dist_param *Ein_param );

  //! Increments the transfer matrix
  //! \param transfer the transfer matrix
  //! \param Eout_count the current row of the transfer matrix
  //! \param Ein0_orig the lower incident energy determined by the probability data
  //! \param Ein1_orig the higher incident energy determined by the probability data
  //! \param Ein_param the qudrature parameters
  void update_T( T_matrix& transfer, int Eout_count,
                 double Ein0_orig, double Ein1_orig, joint_dist_param *Ein_param );


public:
  //! The I=1 angular probability densities
  angle_dist angle_data;

  two_d_interp Ein_interp; // interpolation rule for incident energy
  Interp_Type Eout_interp; // interpolation rule for outgoing energy
  two_d_interp mu_interp; // interpolation rule for direction cosine

  bool ENDL_data;

  //! Default constructor
  joint_dist( ): Eout_interp( LINLIN ), ENDL_data( false )
  {}

  //! Default destructor
  inline ~joint_dist( )
  {}

  //! Reads the python data
  //! \param infile input file
  //! \param num_Ein number of incident energies for this reaction
  void read_data( data_parser &inFile, int num_Ein );

  //! Calculates the transfer matrix for this particle.
  //! \param sigma the cross section data
  //! \param mult the outgoing particle multiplicity data
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param transfer the transfer matrix
  void get_T( const dd_vector& sigma, const dd_vector& mult, 
	      const dd_vector& weight, T_matrix& transfer );
};

// ----------- class joint_dist_param -----------------
//! parameters for 3-d integration
class joint_dist_param: public param_base
{
public:
  // where the lower and upper physical Eout values hit the quadrature box
  joint_dist_hits lower_mu0_hits;  // lower mu, lower Eout
  joint_dist_hits lower_mu1_hits;  // higher mu, lower Eout
  joint_dist_hits upper_mu0_hits;  // lower mu, higher Eout
  joint_dist_hits upper_mu1_hits;  // higher mu, higher Eout

  //! where we are in the list, which incident energies
  joint_dist::iterator this_Ein_ptr;
  joint_dist::iterator next_Ein_ptr;

  //! where we are in the data, which mu values
  one_joint_dist::iterator this_Ein_this_mu;  // lower Ein, lower mu
  one_joint_dist::iterator this_Ein_next_mu;  // lower Ein, higher mu
  one_joint_dist::iterator next_Ein_this_mu;  // higher Ein, lower mu
  one_joint_dist::iterator next_Ein_next_mu;  // higher Ein, higher mu

  //! where we are in the (Eout, probability) data
  one_mu::const_iterator Ein0_mu0_Eout0;  // lower Ein, lower mu, lower Eout
  one_mu::const_iterator Ein0_mu0_Eout1;  // lower Ein, lower mu, higher Eout
  one_mu::const_iterator Ein0_mu1_Eout0;  // lower Ein, higher mu, lower Eout
  one_mu::const_iterator Ein0_mu1_Eout1;  // lower Ein, higher mu, higher Eout
  one_mu::const_iterator Ein1_mu0_Eout0;  // higher Ein, lower mu, lower Eout
  one_mu::const_iterator Ein1_mu0_Eout1;  // higher Ein, lower mu, higher Eout
  one_mu::const_iterator Ein1_mu1_Eout0;  // higher Ein, higher mu, lower Eout
  one_mu::const_iterator Ein1_mu1_Eout1;  // higher Ein, higher mu, higher Eout

  //! The corners of the parallelpiped for 3-d interpolation of data
  //! The range of incident energies is as in the data
  double lower_mu;  // the lower direction cosine
  double upper_mu;  // the higher direction cosine
  double lower_Eout;  // the lower outgoing energy
  double upper_Eout;  // the higher outgoing energy

  // parameters for integration over outgoing cosine and energy
  joint_mu_param mu_params;

  //! (Eout, probability) data at the lower incident energy
  current_data Ein0_data;

  //! (Eout, probability) data at the higher incident energy
  current_data Ein1_data;

  //! How the physical lower outgoing energy hits the quadrature box
  energy_hit_list lower_2d_hits;

  //! How the physical higher outgoing energy hits the quadrature box
  energy_hit_list upper_2d_hits;
 
  Quadrature_Method Eout_quad_method;  // quadrature method for outgoing energy
  Quadrature_Method mu_quad_method;  // quadrature method for outgoing cosine

  long int quad_count;  // number of 3-d quadratures
  long int Ein_F_count;  // number of calls to joint_dist_F::Ein_F
  long int Eout_F_count;  // number of calls to joint_dist_F::Eout_F
  long int mu_F_count;  // number of calls to joint_dist_F::mu_F

  //! constructor
  inline joint_dist_param( );

  //! Default destructor
  inline ~joint_dist_param( )
  {}

  //! determines the geometry for 2d integration over outgoing cosine and energy
  //! returns true if the geometry makes sense
  //! \param this_data (Eout, probability) data at the current incident energy
  bool geometry( const current_data &this_data );

  //! Interpolates data to common values of unit-base mu and outgoing energy
  //! \param UB_Eout the desired value of unit-base outgoing energy
  //! \param Ein0_mu0_data computed  E_mu_P_data at lower incident energy, lower_mu
  //! \param Ein0_mu1_data computed  E_mu_P_data at lower incident energy, upper_mu
  //! \param Ein1_mu0_data computed  E_mu_P_data at higher incident energy, lower_mu
  //! \param Ein1_mu1_data computed  E_mu_P_data at higher incident energy, upper_mu
  void common_mu_Eout( double UB_Eout, E_mu_P_data *Ein0_mu0_data,
    E_mu_P_data *Ein0_mu1_data, E_mu_P_data *Ein1_mu0_data,
    E_mu_P_data *Ein1_mu1_data );

};

// ************* functions to integrate ******************
namespace joint_dist_F
{
  // ---------------- joint_dist_F::Eout_F ------------------
  //! Function for the 1-d quadrature over outgoing energy
  //! \param Eout the lab-frame energy of the outgoing particle
  //! \param Eout_quad_param the function parameters
  //! \param value the value of the integrand, a set of Legendre coefficients
  void Eout_F( double Eout, QuadParamBase *Eout_quad_param, coef_vector *value );

  // ---------------- joint_dist_F::mu_F ------------------
  //! Function for the 2-d quadrature over outgoing energy and cosine
  //! \param mu the lab-frame direction cosine of the outgoing particle
  //! \param mu_quad_param the function parameters
  //! \param value the value of the integrand, a set of Legendre coefficients
  void mu_F( double mu, QuadParamBase *mu_quad_param, coef_vector *value );

  // ---------------- joint_dist_F::Ein_F ------------------
  //! Function for the 3-d quadrature over incident energy, cosine, and outgoing energy
  //! \param E_in the energy of the incident particle
  //! \param Ein_quad_param the function parameters
  //! \param value the value of the integrand, a set of Legendre coefficients
  void Ein_F( double E_in, QuadParamBase *Ein_quad_param, coef_vector *value );
}

#endif
