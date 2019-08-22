/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2015-05-20 -0800 (Wed, 20 May 2015) $
 * $Author: hedstrom $
 * $Id: gamma_2_body_map.hpp 1 2015-05-20 03:06:56Z hedstrom $
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

// header for the gamma_2_body_map class

#include "mappings.hpp"  // for particleInfo

#ifndef DEF_GAMMA2BODYMAP
#define DEF_GAMMA2BODYMAP

// ----------- class gamma_2_body_map -----------------
//! Class to handle relativistic gammas from discrete 2-body reactions
class gamma_2_body_map
{
private:

public:
  particleInfo *rest_masses;
  double threshold;     // threshold energy for the reaction

  double Tin_lab;  // lab kinetic energy of the incident particle
  double Minkowski;  // the Minkowski length, sqrt( E*E - p*p*c*c )
  double T_cm_out;  // cm kinetic energy of the outgoing gamma

  double cosh_chi; // for the boost between frames
  double sinh_chi; // for the boost between frames

  //! Do we do the boost relativistically?
  bool relativistic;

  //! Default constructor
  inline gamma_2_body_map( ): relativistic( true ) {}

  //! Default destructor
  inline ~gamma_2_body_map() {}

  //! Saves the rest masses
  void setup_params( particleInfo *to_save );

  //! Sets up the boost to the lab frame
  //! \param T_in_lab lab-frame kinetic energy of incident particle
  void set_boost( double T_in_lab );

  //! Calculates the center-of-mass energy and momentum for discrete 2-body reactions
  void get_T_cm_out( );

  //! Gets the value of mu_cm given the laboratory energy of outgoing gamma.
  //! \param T_lab lab-frame kinetic energy of incident particle
  double get_mu_cm( double T_lab );

  //! Gets the laboratory energy and cosine of outgoing gamma
  //! \param mu_cm center-of-mass direction cosine of ejected particle
  //! \param Tout_lab computed lab-frame kinetic energy of ejected particle
  //! \param mu_lab computed lab-frame direction cosine of ejected particle
  void get_E_mu_lab( double mu_cm, double *Tout_lab, double *mu_lab );

};

#endif
