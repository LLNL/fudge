/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2009-09-15  (Tue., Sept. 15, 2009) $
 * $Author: hedstrom $
 * $Id: Ecm_Elab_geom.hpp 1 2009-09-15 hedstrom $
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
//! Defines the classes used for the geometry of the map from center of mass to lab frame



#ifndef ECM_ELAB_GEOM_CLASS
#define ECM_ELAB_GEOM_CLASS

#include "Vcm_Vlab_Hit.hpp"

//! Class for parameters for the 2-d quadrature over cm cosine and Eout_cm
// ---------------- class Ecm_Elab_Ecm_param ------------------
class Ecm_Elab_Ecm_param : public QuadParamBase
{
public:
  map_cm_lab *map;

  // the data entries for this incident energy
  list< Vcm_quadBox_Hit > V_cm_limits;  // a list of values of V_cm which give intersections

  // The V_lab values for this lab E_out bin
  double V_lab_min;
  double V_lab_max;

  double E_in;
  double Ecm_min;  // actual range of integration
  double Ecm_max;
  double data_Ecm_min;  // range of data values
  double data_Ecm_max;
  double min_V_cm;
  double max_V_cm;
  Hit_Corner min_hit_corner;
  Hit_Corner max_hit_corner;

  // lab outgoing energy range
  double lab_Eout_min;
  double lab_Eout_max;

  inline Ecm_Elab_Ecm_param( ) {}
  inline ~Ecm_Elab_Ecm_param( ) {}

  //! Sets up the data for this incident energy
  //! \param Ein energy of the incident particle
  //! \param Eoutmin minimum lab-frame energy for this outgoing energy bin
  //! \param Eoutmax maximum lab-frame energy for this outgoing energy bin
  //! \param Ecmmin minimum outgoing center-of-mass energy for this data
  //! \param Ecmmax maximum outgoing center-of-mass energy for this data
  void setup( double Ein, double Eoutmin, double Eoutmax, double Ecmmin, double Ecmmax);

  //! Identifies the regions of integration over Eout_lab and mu_cm.
  void V_lab_sectors( );

  //! Computes the range of center-of-mass outgoing energies.
  //! Returns the tolerance to use in the quadrature over center-of-mass energy and cosine
  double Ecm_range( );

  //! Computes the center-of-mass outgoing energy for a given Hit_Corner
  //! \param hit_corner identifies the type of corner
  //! \param V_trans the velocity of the center of mass in the lab frame
  //! \param V_bottom velocity corresponding to the bottom of the lab energy bin
  //! \param V_top velocity corresponding to the top of the lab energy bin
  //! \param data_Ecm the center of mass energy of the current data, used if this is not a corner
  double get_Ecm( Hit_Corner hit_corner, double V_trans,
    double V_bottom, double V_top, double data_Ecm );
};

//! Class for parameters for the 1-d quadrature over cm cosine
// ---------------- class Ecm_Elab_mu_param ------------------
class Ecm_Elab_mu_param : public QuadParamBase
{
public:
  double E_in;
  double Eout_cm;
  double Ecm_prob;  // probability density of outgoing energy
  double mu_cm_min;
  double mu_cm_max;
  map_cm_lab *map;

  inline Ecm_Elab_mu_param( ) {}
  inline ~Ecm_Elab_mu_param( ) {}

  //! Sets up the data for this incident energy and this Eout_cm
  //! \param Ein energy of the incident particle
  //! \param Eoutcm center-of-mass energy for this data
  //! \param Ecm_param parameters for quadrature over center-of-mass energy
  void setup( double Ein, double Eoutcm, const Ecm_Elab_Ecm_param& Ecm_param );
};

//! Class for parameters for the 3-d quadrature over Ein, cm cosine, and Eout_cm
// ---------------- class Ecm_Elab_Ein_param ------------------
class Ecm_Elab_Ein_param : public param_base
{
private:

public:
  map_cm_lab *map;

  // parameters for 2-d integration over cm cosine and Eout_cm
  Ecm_Elab_Ecm_param Ecm_params;

  Vcm_quadBox_Hit Vcm_hit_min;  // range of values of V_cm for one quadrature sector
  Vcm_quadBox_Hit Vcm_hit_max;

  // pointer to the lower lab outgoing energy
  vector< double >::const_iterator Eout_ptr;

  Ecm_Elab_Ein_param( ) {}
  ~Ecm_Elab_Ein_param( ) {}

  //! Initializes the quadrature parameters
  //! \param sigma the cross section data
  //! \param multiple the outgoing particle multiplicity data
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param e_flux the approximate flux used to weight the transfer matrix
  //! \param Ein_groups the incident-energy group boundaries
/*  void setup( const dd_vector& sigma, const dd_vector& multiple,
    const dd_vector& weight, const Flux_List& e_flux,
    const Energy_groups& Ein_groups );*/
};

#endif
