/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2015-05-20 -0800 (Wed, 20 May 2015) $
 * $Author: hedstrom $
 * $Id: gamma_2_body_map.hpp 1 2015-05-20 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
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
