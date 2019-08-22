/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2009-08-11 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: Eout_integrals.hpp 1  2009-08-11 03:06:56Z hedstrom $
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

// header for the vector of integrals over Eout/mu
#ifndef EOUT_INTEGRALS
#define EOUT_INTEGRALS

#include <list>

#include "coef_vector.hpp"

using namespace std;

// --------------- Eout_int_param ------------------------
//! functions needed to compute the integrals over Eout and mu
class Eout_int_param
{
public:
  // The range of integration
  double Eout_0;
  double Eout_1;

  inline Eout_int_param( ) {}
  virtual ~Eout_int_param( ) {}

// ************** virtual routines ******************************
  //! Initialize at a given incident energy
  //! \param E_in energy of the incident particle
  virtual void set_Ein( double E_in ) = 0;

  //! Evaluate the integrals over Eout
  //! \param Eout_0 lower outgoing energy limit of integration 
  //! \param Eout_1 upper outgoing energy limit of integration 
  //! \param value computed value of the integrals
  virtual void get_integrals( double Eout_0, double Eout_1, coef_vector &value ) = 0;

  //! Evaluate the integrals over Eout and return the noise in the calculation
  //! \param Eout_0 lower outgoing energy limit of integration 
  //! \param Eout_1 upper outgoing energy limit of integration 
  //! \param value computed value of the integrals
  virtual double tol_get_integrals( double Eout_0, double Eout_1, coef_vector &value ) = 0;
// ************** end of virtual routines ***********************
};

// --------------- Eout_link ---------------------------
//! Integrals over mu and one E_out bin for one incident energy
class Eout_link: public coef_vector
{
public:
  double E_in;

  //! Constructor
  inline Eout_link( ) {}

  //! Destructor
  inline ~Eout_link( ) {}

  //! Linear interpolation between this link and the next
  //! \param Ein the energy to interpolate to
  //! \param next_link data for the next link
  //! \param interp the interpolated data
  void Interpolate( double Ein, const Eout_link& next_link,
		    coef_vector *interp ) const;
};

// --------------- Eout_integrals ---------------------------
class Eout_integrals: public list< Eout_link >
{
public:
  Eout_int_param *Eout_int_params;

  //! Constructor
  inline Eout_integrals( ) {}

  //! Destructor
  inline ~Eout_integrals( ) {}

  //! Adds a new link to the Eout_ints list
  //! \param where the new link is inserted before this one
  //! \param order the Legendre order of the link
  //! \param conserve flag to conserve energy or particle number or both
  //! \param E_in energy of the incident particle
  void new_Eout_int( Eout_integrals::iterator where, int order, 
    Conserve conserve, double E_in );

  //! Adds a new link to the Eout_ints list and returns the noise in the calculation
  //! \param where the new link is inserted before this one
  //! \param order the Legendre order of the link
  //! \param conserve flag to conserve energy or particle number or both
  //! \param E_in energy of the incident particle
  double tol_new_Eout_int( Eout_integrals::iterator where, int order,
    Conserve conserve, double E_in );

  //! Sets up the integrals over E_out for interpolation
  //! \param order the Legendre order of the link
  //! \param conserve flag to conserve energy or particle number or both
  //! \param Ein_0 lower energy of the incident particle
  //! \param Ein_1 higher energy of the incident particle
  void setup_Eout_ints( int order, Conserve conserve, double Ein_0, double Ein_1 );

  //! Check the accuracy of linear interpolation
  //! \param interp_link a link at an intermediate outgoing energy
  //! \param prev_link a link at a lower outgoing energy
  //! \param next_link a link at a higher outgoing energy
  //! \param tol the required tolerance for linear interpolation
  bool interp_OK( const Eout_link& interp_link,
		  Eout_integrals::iterator prev_link,
		  Eout_integrals::iterator next_link, double tol );
};

#endif
