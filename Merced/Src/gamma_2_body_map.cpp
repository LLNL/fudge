/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2015-05-20 -0800 (Wed, 20 May 2015) $
 * $Author: hedstrom $
 * $Id: gamma_2_body_map.cpp 1 2015-05-20 03:06:56Z hedstrom $
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

// implementation for the gamma_2_body_map class

#include <fstream>  // standard file stream package
#include <iostream>
#include <cmath>

#include "gamma_2_body_map.hpp"
#include "math_util.hpp"
#include "messaging.hpp"
#include "global_params.hpp"

using namespace std;

// ********* class gamma_2_body_map **********************
// ------------- gamma_2_body_map::setup_params --------------------
// Saves the rest masses
void gamma_2_body_map::setup_params( particleInfo *to_save )
{
  rest_masses = to_save;
  double rest_mass_in = rest_masses->mProj + rest_masses->mTarg;
  threshold = ( rest_masses->mRes * rest_masses->mRes - rest_mass_in * rest_mass_in )/
    ( 2.0*rest_masses->mTarg );
  if( threshold < 0.0 )
  {
    // exothermic reaction
    threshold = 0.0;
  }
}
// ------------- gamma_2_body_map::set_boost --------------------
// Sets up the boost to the lab frame
void gamma_2_body_map::set_boost( double T_in_lab )
{
  Tin_lab = T_in_lab;
  double rest_mass_init = rest_masses->mTarg + rest_masses->mProj;
  if( relativistic )
  {
    Minkowski = sqrt( rest_mass_init*rest_mass_init +
          2*rest_masses->mTarg*Tin_lab );
  }
  else
  {
    Minkowski = rest_mass_init + rest_masses->mTarg*Tin_lab/rest_mass_init;
  }
  double numerator;
  if( relativistic )
  {
    numerator = sqrt( Tin_lab*( 2*rest_masses->mProj + Tin_lab ) );
  }
  else
  {
    numerator = sqrt( Tin_lab*2*rest_masses->mProj );
  }
  sinh_chi = numerator/Minkowski;
  cosh_chi = sqrt( 1.0 + sinh_chi*sinh_chi );
  get_T_cm_out( );
}
// ------------- gamma_2_body_map::get_E_mu_lab --------------------
//  Gets the laboratory energy and cosine of outgoing gamma
void gamma_2_body_map::get_E_mu_lab( double mu_cm, double *Tout_lab, double *mu_lab )
{
  double factor = cosh_chi + mu_cm*sinh_chi;
  *Tout_lab = factor*T_cm_out;
  *mu_lab = ( mu_cm*cosh_chi + sinh_chi )/factor;
}
// ------------- gamma_2_body_map::get_T_cm_out --------------------
// Calculates the center-of-mass energy and momentum for discrete 2-body reactions
void gamma_2_body_map::get_T_cm_out( )
{
  // knock-on reaction
  double E_diff = Minkowski - rest_masses->mRes;
  T_cm_out = E_diff * ( Minkowski + rest_masses->mRes ) /( 2*Minkowski );
}
// ------------- gamma_2_body_map::get_mu_cm --------------------
// Gets the value of mu_cm given the laboratory energy of outgoing gamma.
double gamma_2_body_map::get_mu_cm( double T_lab )
{
  // This routine uses the fact that T_lab = (cosh_chi + mu_cm*sinh_chi)*T_cm
  if( sinh_chi <= 0.0 )
  {
    FatalError( "gamma_2_body_map::get_mu_cm", "zero incident energy" );
  }
  double mu_cm = ( T_lab/ T_cm_out - cosh_chi )/sinh_chi;
  if( mu_cm > 1.0 )
  {
    Warning ( "gamma_2_body_map::get_mu_cm", pastenum( "mu_cm too big: ", mu_cm ) );
    mu_cm = 1.0;
  }
  if( mu_cm < -1.0 )
  {
    Warning ( "gamma_2_body_map::get_mu_cm", pastenum( "mu_cm too small: ", mu_cm ) );
    mu_cm = -1.0;
  }
  return mu_cm;
}
