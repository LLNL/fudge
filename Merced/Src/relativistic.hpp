/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2015-05-20 -0800 (Wed, 20 May 2015) $
 * $Author: hedstrom $
 * $Id: relativistic.hpp 1 2015-05-20 03:06:56Z hedstrom $
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

// header for the classes used in relativistic mechanics

#ifndef DEF_RELATIVISTIC
#define DEF_RELATIVISTIC

#include "mappings.hpp"  // for the particleInfo class

// ----------- class relativistic_masses -----------------
//! Class for mass data used in relativistic boosts
class relativistic_masses
{
private:

public:
  particleInfo *rest_masses;
  double M_total;  // sum of all particle masses

  //! the threshold for the reaction
  double threshold;

  double Q_value;

  //! Default constructor
  inline relativistic_masses( ): Q_value(0.0) {}

  //! Default destructor
  inline ~relativistic_masses() {}

  //! Saves the rest masses
  //! \param to_save the rest masses of the particles
  //! \param file_Q the Q value given in the input file
  void setup_masses( particleInfo *to_save, double file_Q );

  //! Calculates the threshold energy
  void get_threshold( );

  //! For exothermic reactions, returns the Newton incident energy for minimal outgoing lab energy
  double Newton_min_Eout( );

  //! Returns the incident energy for zero outgoing lab-frame energy, Newtonian
  double Newton_zero_Eout( );
};
// ----------- class relativistic_param -----------------
//! Class for the parameters in relativistic_F functions
class relativistic_param
{
private:
  //! Gets the momentum from the kinetic energy
  //! \param T particle kinetic energy
  //! \param E0 particle rest mass
  double p_from_T( double T, double E0 );

  //! Gets the kinetic energy from the momentum
  //! \param p particle momentum
  //! \param E0 particle rest mass
  double T_from_p( double p, double E0 );

public:
  relativistic_masses *masses;  // the mass ratios

  double Tin_lab;  // lab kinetic energy of the incident particle
  double pin_lab;  // lab momentum of the incident particle
  double pin_cm;   // center-of-mass momentum of the incident particle
  double Minkowski;  // the Minkowski length, sqrt( E*E - p*p*c*c )
  double T_cm_out;  // cm kinetic energy of the outgoing particle
  double p_cm_out;  // cm momentum of the outgoing particle
  double mu_cm;   // for use in relativistic_F::T_out_lab

  double cosh_chi; // for the boost between frames
  double sinh_chi; // for the boost between frames

  //! Default constructor
  inline relativistic_param( ) {}

  //! Default destructor
  inline ~relativistic_param() {}

  //! Sets up the boost to the lab frame
  //! \param T_in_lab lab frame kinetic energy of the incident particle
  void set_boost( double T_in_lab );

  //! Boosts the ejected particle from the center-of-mass to the lab frame
  //! \param mu_cm direction cosine in the center-of-mass frame
  //! \param T_lab computed lab frame kinetic energy
  //! \param p_lab_parallel computed lab frame parallel component of the momentum
  void boost( double mu_cm, double *T_lab, double *p_lab_parallel );

  //! Calculates the center-of-mass energy and momentum for discrete 2-body reactions
  void get_p_cm_out( );

  //! Calculates the lab-frame kinetic energy for discrete 2-body reactions
  //! \param mu_cm center-of-mass frame direction cosine
  double get_T_lab_out( double mu_cm );

  //! Gets the value of mu_cm given the laboratory energy of outgoing particle.
  //! \param T_lab the desired lab frame kinetic energy
  double get_mu_cm( double T_lab );

  //! Gets the laboratory energy and cosine of outgoing particle
  //! \param mu_cm direction cosine in the center-of-mass frame
  //! \param Tout_lab computed lab frame kinetic energy
  //! \param mu_lab computed lab frame direction cosine
  void get_E_mu_lab( double mu_cm, double *Tout_lab, double *mu_lab );

  //! Returns the incident energy for zero outgoing lab-frame energy, relativistic
  double zero_Eout( );

  //! For negative mu, gets the minimal lab frame outgoing energy and the
  //! corresponding incident energy
  //! \param mu the direction cosine in the center-of-mass frame
  //! \param guess_Ein an initial guess of the required incident energy
  //! \param Eout the computed minimal outgoing energy in the lab frame
  //! \param Ein the computed incident energy for this Eout
  void get_min_Eout( double mu, double guess_Ein, double *Eout, double *Ein );

  //! Finds the incident energies for given T_lab and mu_cm.
  //! Returns the number of solutions: 0, 1, or 2
  //! \param T_lab lab-frame kinetic energy of ejected particle
  //! \param mu_cm center-of-mass direction cosine of ejected particle
  //! \param pair_0 ( Ein, Eout ) at lower incident energy
  //! \param pair_1 ( Ein, Eout ) at higher incident energy
  double find_hit( double E_lab, double mu_cm, const dd_entry &pair_0,
		   const dd_entry &pair_1 );

  //! For Q < 0 and mu < 0 find the incident energy which minimizes Eout.
  //! On exit pair_0 and pair_1 are tighter bounds.
  //! \param mu: the direction cosine in the center-of-mass frame
  //! \param pair_0: the lower bound and its outgoing energy
  //! \param pair_1: the upper bound and its outgoing energy
  //! \param tol: an error tolerance
  dd_entry find_bottom( double mu, dd_entry *pair_0, dd_entry *pair_1,
     double tol );

  //! Find the incident energy which minimizes Eout for mu = -1
  dd_entry find_lowest_bottom(  );
};

// ************* functions *************
namespace relativistic_F
{
  // ------------------ T_out_lab -----------------------
  //! Returns the kinetic energy of the emitted particle in the lab frame
  //! \param T_in_lab kinetic energy of the incident particle in the lab frame
  //! \param params the data for the relativistic boost
  double T_out_lab( double T_in_lab, void *params );

  // ------------------ p_out_lab -----------------------
  //! Returns the parallel momentum of the emitted particle in the lab frame.
  //! Routine is used to find when back scattering gives zero outgoing energy.
  //! \param T_in_lab kinetic energy of the incident particle in the lab frame
  //! \param params the data for the relativistic boost
  double p_out_lab( double T_in_lab, void *params );
}

#endif
