/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2008-04-16 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: reaction.hpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// define the classes used for identification of the reaction data

#ifndef REACTION_DEF
#define REACTION_DEF

#include <list>
#include <cmath>
#include <cstdio>
#include <fstream>

#include "param_base.hpp"  // for the quadrature method flag
#include "dd_vector.hpp"  // for the cross section
#include "energy_dist.hpp"
#include "Legendre2Body.hpp"
#include "angle_dist.hpp"
#include "isotropic.hpp"
#include "joint_dist.hpp"
#include "standard_Legendre.hpp"
#include "uncorrelated.hpp"
#include "cm_Legendre.hpp"
#include "cm_doubleDiff.hpp"
#include "lab_doubleDiff.hpp"
#include "Compton.hpp"
#include "coherent.hpp"
#include "evaporation.hpp"
#include "Maxwell.hpp"
#include "Watt.hpp"
#include "MadlandNix.hpp"
#include "kalbach.hpp"
#include "phase_space.hpp"
#include "general_evap.hpp"
#include "two_step.hpp"
#include "transfer.hpp"
#include "Legendre_data.hpp"
#include "data_parser.hpp"

namespace React
{
//! laboratory or center-of-mass frame
enum Frame{ LAB, CM };


//! Class of parameters for integrating multiplicity * flux * sigma * model_weight
//-----------class num_check_param ----------
class num_check_param
{
public:
  Ddvec::dd_pair sigma;  // cross section
  Ddvec::dd_pair mult;   // multiplicity
  Ddvec::dd_pair e_flux; // flux
  Ddvec::dd_pair model_weight;

  inline num_check_param(){}  // constructor
  inline ~num_check_param(){}  // destructor

  //! Computes the product sigma*multiplicity*flux*weight
  //! \param E_in energy of the incident particle
  double value( double E_in );
};

//! Class for reactions
//-----------class reaction ----------
class reaction
{
private:
  Trf::T_matrix transfer;
  Maps::particleInfo particle_info;  // used for center-of-mass data

  std::string comment;

  //! The weight used by the energy_function class
  Ddvec::dd_vector model_weight;

  double version;

  React::Frame frame_in;  // frame for incident particle
  React::Frame frame_out;  // frame for outgoing particle

  //! Interprets the method of integration
  //! Sets the quadrature flags in the Trf::T_matrix, transfer
  //! \param input_file the input file
  Qmeth::Quadrature_Method read_quadrature( Dpar::data_parser &input_file );

public:

  //! cross section
  Ddvec::dd_vector cross_section;
 
  //! angular distributions
  Adist::angle_dist angdist;

  //! angular distributions as Legendre polynomials
  Lg2b::Legendre_angle_dist LegendreAngle;

  //! isotropic energy tables
  Iso::isotropic isotrop;

  //! Legendre moments of energy distributions using quadrature
  Edist::energy_moments energyMoments;

  //! ENDF Legendre moments with unit-base or lin-lin interpolation
  StdLg::standard_Legendre standard_legendre;

  //! joint energy-angle distribution tables
  Jdist::joint_dist joint_dist_table;

  //! data for uncorrelated joint energy-angle distributions using quadrature
  Uncor::uncorrelated uncorr;

  //! data for coherent scattering
  Coh::coherent coherent_model;

  //! data for Compton scattering
  Comp::Compton compton;

  //! data for the evaporation model
  Evap::evaporation evaporation_model;

  //! data for the Maxwell model
  Max::maxwell Maxwell_model;

  //! data for the Watt model
  Wtt::watt Watt_model;

  //! data for the Madland-Nix model
  MNix::MadlandNix MadlandNix_model;

  //! data for the Kalbach model
  Kbach::Kalbach Kalbach_model;

  //! data for the phase_space model
  Phase::phase_space phase_space_model;

  //! data for Legendre expansions in the center-of-mass frame
  CmLeg::cm_Legendre cm_Legendre_model;

  //! data for energy-angle distributions in the center-of-mass frame
  cmDD::joint_dist cm_joint_table;

  //! data for tabular energy-angle distributions in the lab frame
  labDD::lab_joint_dist lab_joint_table;

  //! data for the general evaporation model
  Gevap::general_evap gen_evap_model;

  //! data for a 2-step 2-body reaction
  Tstep::two_step twoStep;

  //! multiplicities
  Ddvec::dd_vector multiple;

  inline reaction() : frame_in( React::LAB ) {}   // constructor
  inline ~reaction(){}  // destructor

  //! process data
  //! \param input_file the input file
  //! \param output_file the output file
  void process_data( Dpar::data_parser &input_file, std::ofstream *output_file );

  //! Reads input data
  //! \param input_file the input file
  void read_input( Dpar::data_parser &input_file );

  //! Reads the data common to all input data; returns true if dataID is found
  //! \param dataID the data identifier
  //! \param input_file the input file
  bool common_input( const std::string &dataID, Dpar::data_parser &input_file );

  //! Reads Global parameters from the input file
  //! \param dataID the data identifier
  //! \param input_file the input file
  bool read_Global( const std::string &dataID, Dpar::data_parser &input_file );

  //! Gets the transfer matrix from discrete 2-body angular data
  //! \param input_file the input file
  void two_body( Dpar::data_parser &input_file );

  //! Gets the transfer matrix from discrete 2-body Legendre polynomial angular data
  //! \param input_file the input file
  void Legendre_two_body( Dpar::data_parser &input_file );

  //! Gets the transfer matrix for a 2-step 2-body reaction
  void do_two_step( Dpar::data_parser &input_file );

//! Reads and processes tables of isotropic energy probability density
  //! \param input_file the input file
  void do_isotropic( Dpar::data_parser &input_file );

  //! Reads and processes Legendre expansions of double differential data
  //! \param input_file the input file
  void Legendre( Dpar::data_parser &input_file );

  //! Reads and processes ENDF Legendre expansions of double differential data
  //! \param input_file the input file
  void do_ENDFLegendre( Dpar::data_parser &input_file );

  //! Reads and processes uncorrelated expansions of double differential data
  //! \param input_file the input file
  void do_uncorr( Dpar::data_parser &input_file );

  //! Reads and processes double_differential tabular data
  //! \param input_file the input file
  void do_joint_dist( Dpar::data_parser &input_file );

  //! Reads and processes double_differential tabular cm data
  //! \param input_file the input file
  void do_energy_angle( Dpar::data_parser &input_file );

  //! Reads and processes scattering factor data for Compton scattering
  //! \param input_file the input file
  void do_Compton( Dpar::data_parser &input_file );

  //! Reads and processes scattering factor data for coherent scattering
  //! \param input_file the input file
  void do_coherent( Dpar::data_parser &input_file );

  //! Reads and processes data for evaporation spectra
  //! \param input_file the input file
  void do_evaporation( Dpar::data_parser &input_file );

  //! Reads and processes data for Maxwell spectra
  //! \param input_file the input file
  void do_Maxwell( Dpar::data_parser &input_file );

  //! Reads and processes data for Watt spectra
  //! \param input_file the input file
  void do_Watt( Dpar::data_parser &input_file );

  //! Reads and processes data for MadlandNix spectra
  //! \param input_file the input file
  void do_MadlandNix( Dpar::data_parser &input_file );

  //! Reads and processes data for Kalbach spectra
  //! \param input_file the input file
  void do_Kalbach( Dpar::data_parser &input_file );

  //! Reads and processes data for phase_space spectra
  //! \param input_file the input file
  void do_phase_space( Dpar::data_parser &input_file );

  //! Reads and processes data for the general evaporation spectrum
  //! \param input_file the input file
  void do_gen_evap( Dpar::data_parser &input_file );

  //! Writes out the transfer matrix to be read by Python
  void write_transfer( );

  //! Writes the cross section for gamma data
  void write_xsec( );

  //! Calculates integrals of multiplicity * sigma * flux, to check row sums
  //! \param threshold, the reaction threshold
  void number_check( double threshold );
};

} // end of namespace React

#endif
