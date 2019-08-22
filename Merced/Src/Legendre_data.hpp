/*
* ******** merced: calculate the transfer matrix *********
* $Revision: 1 $
* $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
* $Author: hedstrom $
* $Id: Legendre_data.hpp 1 2006-02-02 03:06:56Z hedstrom $
*
* ******** merced: calculate the transfer matrix *********
*
* # <<BEGIN-copyright>>
* Copyright (c) 2017, Lawrence Livermore National Security, LLC.
* Produced at the Lawrence Livermore National Laboratory.
* Written by the LLNL Nuclear Data and Theory group
*         (email: mattoon1@llnl.gov)
* LLNL-CODE-725546.
* All rights reserved.
* 
* This file is part of the Merced package, used to generate nuclear reaction
* transfer matrices for deterministic radiation transport.
* 
* 
*     Please also read this link - Our Notice and Modified BSD License
* 
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the disclaimer below.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the disclaimer (as noted below) in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
*       to endorse or promote products derived from this software without specific
*       prior written permission.
* 
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
* THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
* 
* 
* Additional BSD Notice
* 
* 1. This notice is required to be provided under our contract with the U.S.
* Department of Energy (DOE). This work was produced at Lawrence Livermore
* National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
* 
* 2. Neither the United States Government nor Lawrence Livermore National Security,
* LLC nor any of their employees, makes any warranty, express or implied, or assumes
* any liability or responsibility for the accuracy, completeness, or usefulness of any
* information, apparatus, product, or process disclosed, or represents that its use
* would not infringe privately-owned rights.
* 
* 3. Also, reference herein to any specific commercial products, process, or services
* by trade name, trademark, manufacturer or otherwise does not necessarily constitute
* or imply its endorsement, recommendation, or favoring by the United States Government
* or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
* herein do not necessarily state or reflect those of the United States Government or
* Lawrence Livermore National Security, LLC, and shall not be used for advertising or
* product endorsement purposes.
* 
* # <<END-copyright>>
*/

// header for the list of (incident energy, flux)
#ifndef E_FLUX
#define E_FLUX

#include <iostream>
#include <list>
#include <vector>

#include "dd_vector.hpp"
#include "data_parser.hpp"

using namespace std;

// --------------- class Legendre_base -------------------------
//! Class to hold Legendre coefficients j=0, 1, ..., order
class Legendre_base
{
private:
double Energy;

public:
//! Legendre coefficients
double *data;

//! the Legendre order
int order;

//! Constructor
inline Legendre_base( ): Energy(-1.0), order(-1) {}

//! Destructor
inline ~Legendre_base( ) { clean_data( ); }

//! "Energy" could be for the incident or the outgoing particle
inline double get_E_out( ) const { return Energy; }
inline double get_E_in( ) const { return Energy; }

//! Sets the energy for the Legendre data
//! \param Eout the energy for the Legendre data
inline void set_E_out( double Eout ) { Energy = Eout; }
//! \param Ein the energy for the Legendre data
inline void set_E_in( double Ein ) { Energy = Ein; }

void clean_data( );

//! Sets the incident energy and allocates space
//! \param order the Legendre order for the data
void initialize( int order );

//! Access routine for N-th coefficient
//! \param N number of the current Legendre coefficient
double &operator[ ]( int N );

//! Access routine for N-th coefficient; does not change the data
//! \param N number of the current Legendre coefficient
double value( int N ) const;

//! Ignore zero high-order Legendre coefficients
void truncate_zeros( );

//! Scales the vector
//! \param factor multiply all coefficients by this number
Legendre_base& operator*=( double factor );

//! Sums the Legenre series
//! \param mu sum the series at this mu value
double sum_Legendre( double mu );

// For debugging
void print( ) const;
};

// --------------- class Legendre_coefs -------------------------
//! Class to hold Legendre coefficients of flux j=0, 1, ..., order
class Legendre_coefs : public Legendre_base
{
private:
//! Interpolates the flux with weight alpha
//! \param alpha the weight for next_flux
//! \param prev_flux Legendre coefficients at a lower energy
//! \param next_flux Legendre coefficients at a higher energy
void basic_linlin_interp( double alpha, const Legendre_coefs& prev_flux,
const Legendre_coefs& next_flux );

public:
//! Constructor
inline Legendre_coefs( )
{}

//! Constructor
//! \param Order the Legendre order
inline Legendre_coefs( int Order )
{ initialize( Order ); }

//! Destructor
inline ~Legendre_coefs( )
{}

//! Allocates space
//! \param Order the Legendre order
void initialize( int order );

//! Sets all coefficients to zero
void zero_data( );

//! Copies the Legendre coefficients
//! Resets the order, if they are inconsistent
//! \param to_copy the Legendre data to copy
void copy_coef( const Legendre_coefs& to_copy );

//! Only copies the Legendre coefficients, order not reset
//! \param to_copy the Legendre data to copy
void only_copy_coef( const Legendre_coefs& to_copy );

///! Sets the order for interpolated data
//! \param left_order the Legendre order for one set of coefficients
//! \param right_order the Legendre order for another set of coefficients
void set_max_order( int left_order, int right_order );

//! Interpolates the flux at energy E_in
//! \param E_in intermediate energy
//! \param prev_flux Legendre coefficients at a lower energy
//! \param next_flux Legendre coefficients at a higher energy
void linlin_interp( double E_in, const Legendre_coefs& prev_flux,
const Legendre_coefs& next_flux );

//! Interpolates Legendre-coefficient data at energy E_in
//! It is required that left_flux and right_flux be at the same outgoing energy
//! \param E_in intermediate incident energy
//! \param left_Ein lower incident energy
//! \param left_flux Legendre coefficients at incident energy left_Ein
//! \param right_Ein higher incident energy
//! \param right_flux Legendre coefficients at incident energy right_Ein
void Ein_linlin_interp( double E_in, double left_Ein, 
const Legendre_coefs& left_flux, double right_Ein,
const Legendre_coefs& right_flux );

//! Interpolates unit-base Legendre-coefficient data
//! It is required that left_flux and right_flux be at the same outgoing energy
//! \param E_in intermediate incident energy
//! \param alpha the proportionality factor
//! \param left_flux Legendre coefficients at incident energy left_Ein
//! \param right_flux Legendre coefficients at incident energy right_Ein
void unitbase_interp( double E_in, double alpha, 
const Legendre_coefs& left_flux,
const Legendre_coefs& right_flux );

//! Interpolates the flux linearly with respect to the logarithm of the energy
//! \param E_in intermediate energy
//! \param left_flux Legendre coefficients at a lower energy
//! \param right_flux Legendre coefficients at a higher energy
void linlog_interp( double E_in, 
const Legendre_coefs& left_flux,
const Legendre_coefs& right_flux );

};

// --------------- class Legendre_data_range -------------------------
//! Class to hold Legendre coefficient data for a range on outgoing energies
class Legendre_data_range
{
private:
double E_in; // the incident energy

public:
Legendre_coefs prev_data;  // Legendre data for lower outgoing energy
Legendre_coefs next_data;  // Legendre data for upper outgoing energy
unit_base_map ubase_map;  // unit-base map for this data

two_d_interp Ein_interp;
Interp_Type Eout_interp;

inline Legendre_data_range( ) {}
inline ~Legendre_data_range( ) {}

//! Sets the energy of the incident particle
//! \param Ein the energy of the incident particle
inline void set_E_in( double Ein )
{
E_in = Ein;
}

//! Gets the energy of the incident particle
inline double get_E_in( ) const
{
return E_in;
}

//! Sets up a new incident energy
//! \param Ein the energy of the incident particle
//! \param ubasemap the unit-base-map for this data
//! \param Eoutinterp rule for interpolation is outgoing energy
void new_Ein( double Ein, const unit_base_map &ubasemap, Interp_Type Eoutinterp );

//! Sets up prev_data and next_data for a given range of outgoing energies
//! \param prevdata Legendre coefficients at a lower energy
//! \param nextdata Legendre coefficients at a higher energy
//! \param Eout_min bottom of the desired energy range
//! \param Eout_max top of the desired energy range
void set_data( const Legendre_coefs &prevdata, const Legendre_coefs &nextdata,
double Eout_min, double Eout_max );

//! Does unit-base interpolation between incident energies
//! It is required that left_data and right_data be at the same outgoing energy
//! \param E_in an intermediate incident energy
//! \param left_data Legendre coefficients at a lower incident energy
//! \param right_data Legendre coefficients at a higher incident energy
void ubase_interpolate( double E_in,
const Legendre_data_range &left_data, const Legendre_data_range &right_data );

//! Maps from physical variables to unit-base
void to_unit_base( );

//! Maps from unit-base to physical variables
void un_unit_base( );

//! Returns the Legendre coefficients for this outgoing energy
//! \param E_out energy of the outgoing particle
Legendre_coefs Eout_interpolate( double E_out );
};

// --------------- class Legendre_list_base -------------------------
//! Class to hold an array Legendre coefficient data
class Legendre_list_base : public list< Legendre_coefs >
{
private:
//! The energy of the incident particle
double E_in;

//! Finds the total probability
double get_norm( ) const;

public:
int order;
unit_base_map ubase_map;
two_d_interp Ein_interp; // interpolation rule for incident energy
  Interp_Type Eout_interp; // interpolation rule for outgoing energy

  //! Constructor
  inline Legendre_list_base( ): order( -1 ) {}

  //! Destructor
  inline ~Legendre_list_base( ) {}

  //! Returns the energy for the Legendre data
  inline double get_E_in( ) const { return E_in; }

  //! Sets the energy for the Legendre data
  //! \param Ein the energy for the Legendre data
  inline void set_E_in( double Ein ) { E_in = Ein; }

  //! Normalizes the total probability
  void renorm( );

  //! Maps the data to unit base
  void to_unit_base( );

  // For debugging
  void print( );

};

// --------------- class Flux_List -------------------------
//! Class to hold Legendre coefficients of flux j=0, 1, ..., order
class Flux_List : public list< Legendre_coefs >
{
private:

public:
  int order;

  Interp_Type interp;

  //! Constructor
  inline Flux_List( ): order( -1 ) {}

  //! Destructor
  inline ~Flux_List( ) {}

  //! Constructs the list from the python data
  //! \param infile the input file
  //! \param num_Ein the number of incident energies
  void read_flux( data_parser &infile, int num_Ein );

  //! Evaluates the flux at energy E_in with search starting at ptr
  //! \param E_in incident energy requested
  //! \param ptr pointer to the most recent successful search
  Legendre_coefs value( double E_in, Flux_List::const_iterator &ptr ) const;
};

// --------------- class weight_vector -------------------------
//! Class to hold Legendre flux weights for one energy bin
class weight_vector : public Legendre_base
{
private:

public:
  //! Constructor
  inline weight_vector( )
  {}

  //! Destructor
  inline ~weight_vector( )
  {}

  //! Adds the integrals over ( E_left, E_right )
  //! \param e_flux the pairs ( incident energy, Legendre coefficients ) to integrate
  //! \param this_flux pointer to the most recent successful search
  //! \param E_left lower limit of integration
  //! \param E_right upper limit of integration
  void increment( Flux_List &e_flux, 
    Flux_List::const_iterator this_flux, double E_left, double E_right );

  //! Takes the reciprocals
  void invert( );
};
// --------------- class weight_list -------------------------
//! Class to hold all of the Legendre flux weights
class weight_list : public list<  weight_vector >
{
private:

public:
  //! Constructor
  inline weight_list( )
  { }

  //! Destructor
  inline ~weight_list( )
  { }

};

#endif
