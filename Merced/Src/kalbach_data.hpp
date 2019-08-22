/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2009-09-15  (Tue., Sept. 15, 2009) $
 * $Author: hedstrom $
 * $Id: kalbach_data.hpp 1 2009-09-15 hedstrom $
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
//! Defines the class used for data for the Kalbach model

#ifndef KALBACH_DATA_CLASS
#define KALBACH_DATA_CLASS

#include "dd_vector.hpp"
#include "mappings.hpp"

//! class to identify a nucleon
// ---------------- class nucleon --------------------
class nucleon
{
private:

public:
  double mass;
  double Kalbach_I;  // particle separation energy
  double Kalbach_m;
  int Kalbach_M;
  int Z;
  int A;

  nucleon( ): mass( -1.0 ), A( -1 ) {}
  ~nucleon( ) {}

  //! copies the data
  nucleon& operator=( const nucleon& to_copy );

  //! Sets Z, A, I, M, m
  //! \param ZA the ENDL ZA number: A + 1000*Z
  void set_params( int ZA );

  //! Returns ZA
  inline int get_ZA( )
    {return A + 1000*Z;}

  //! Computes the Kalbach S_a function
  //! \param proj the identifiers for the projectile
  double get_Sa( const nucleon &proj );

  //! Checks for proper initialization
  bool check_data( );
};

//! class for the Kalbach a coefficient
// ---------------- class Kalbach_a ------------------
class Kalbach_a
{
private:
  // Coefficients of the Kalbach polynomial
  double C_1;
  double C_2;
  double C_3;

  double projectile_S;
  double eject_S;

public:
  nucleon target;
  nucleon projectile;
  nucleon eject;
  nucleon compound;
  nucleon residual;
  map_cm_lab *map;

  Kalbach_a( ): C_1( 0.04 ), C_2( 1.8e-6 ), C_3( 6.7e-7 ) {}
  ~Kalbach_a( ) {}

  //! Initializes the mass ratios in map
  void setup_params( );

  //! Computes the S_a and S_b functions
  void set_Sab( );

  //! Computes the Kalbach a function
  //! \param E_in the energy of the incident particle
  //! \param E_out the energy of the outgoing particle
  double get_a( double E_in, double E_out );

  //! Stores the masses
  //! \param particle_info the identities of the particles involved in the reaction
  void copy_masses( const particleInfo &particle_info );
};

//! Class for the current data entries for one incident energy
// ---------------- class kalbach_data -----------------------
class kalbach_data
{
public:
  double E_in;
  // the data entries for this incident energy
  double this_Ecm; // center-of-mass energy
  double this_f0;  // center-of-mass energy probability density
  double this_r;   // r-value
  double next_Ecm;
  double next_f0;
  double next_r;

  two_d_interp Ein_interp; // interpolation with respect to incident energy
  Interp_Type Eout_interp; // interpolation with respect to outgoing energy

  inline kalbach_data( ): Eout_interp( HISTOGRAM ) {}
  inline ~kalbach_data( ) {}

  unit_base_map unit_base;  // unit_base for this incident energy

  //! Does linear interpolation of this data with the next
  //! \param E_in the energy of the incident particle
  //! \param left_data Kalbach data at a lower incident energy
  //! \param right_data Kalbach data at a higher incident energy
  void linlin_interp( double E_in, const kalbach_data& left_data,
    const kalbach_data& right_data );

  //! Does linear interpolation of this unit_base data with the next
  //! \param E_in the energy of the incident particle
  //! \param left_data Kalbach data at a lower incident energy
  //! \param right_data Kalbach data at a higher incident energy
  void unit_base_interp( double E_in, const kalbach_data& left_data,
    const kalbach_data& right_data );

  //! Undoes the unit-base map for one outgoing energy
  //! \param E_unit the unit-base energy of the outgoing particle
  double un_unit_base( double E_unit );

  //! Undoes the unit-base map; used on interpolated data
  void un_unit_base( );

  //! Calculates r and the center-of-mass outgoing probability density
  //! \param Eoutcm center-of-mass outgoing energy
  //! \param computed probability density for this energy
  //! \param computed Kalbach r parameter for this energy
  void get_f0_r( double Eoutcm, double *Ecm_prob, double *r ) const;

};


#endif
