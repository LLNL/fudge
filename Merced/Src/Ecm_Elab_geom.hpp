/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2009-09-15  (Tue., Sept. 15, 2009) $
 * $Author: hedstrom $
 * $Id: Ecm_Elab_geom.hpp 1 2009-09-15 hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
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
