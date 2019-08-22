/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2011-04-29 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: energy_dist.hpp 1 2011-04-29 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// define the classes used for ENDL i=4 energy distributions

#ifndef ENERGY_DIST_CLASS
#define ENERGY_DIST_CLASS

#include <iostream>
#include "data_parser.hpp"
#include "standard_Legendre.hpp"
#include "energy_dist_base.hpp"

using namespace std;

//! Class for energy distributions, one Legendre order
//--------------- class energy_dist ----------------
class energy_dist
{
private:

public:
  //! the ENDL data
  Eprob_vector *EProb_data;

  //! the number of incident energies
  int number_Ein;
  two_d_interp Ein_interp;  // interpolation between incident enrgies
  Interp_Type Eout_interp;  // interpolation between outgoing enrgies

  inline energy_dist( ): number_Ein( 0 ) {}

  ~energy_dist( );

  //! Reads the ENDL data
  //! \param infile input file
  //! \param num_Ein number of incident energies for this reaction
  void read_data( data_parser& infile, int num_Ein );

  //! Maps the data to unit base
  //! \param L_order the Legendre order of the current data
  void unit_base( int L_order );
};

//! list: one link for each Legendre order
//--------------- class energy_moments ----------------
class energy_moments
{
 private:
  int output_order;  // the Legendre order of the output

  //! Converts ENDL data to ENDF format
  void to_ENDF( );

  //! Converts ENDL data to ENDF format for one incident energy
  //! \param Ein_count identifies the current incident energy bin
  void one_Ein_to_ENDF( int Ein_count );

  //! Checks to see that the incident energies are consistent for all Legendre orders
  void check_Ein( ) const;

 protected:
  //! Use the coding for this data in ENDF format
  standard_Legendre ENDF_data;

  //! The original ENDL data
  energy_dist *ENDL_data;

  int data_order;  // the Legendre order of the data
 
  //! Converts isotropic ENDL data to ENDF format
  void zero_order( );

 public:
  two_d_interp Ein_interp;  // interpolation between incident enrgies
  Interp_Type Eout_interp;  // interpolation between outgoing enrgies

  inline energy_moments( ): output_order( -1 ), data_order( -1 ) {}
  ~energy_moments( );

  //! Reads the ENDL data
  //! \param infile input file
  //! \param num_moments number of Legendre moments for this reaction
  void read_data( data_parser& input_file, int num_moments );

   // Calculates the transfer matrix for this particle.
  //! \param sigma the cross section data
  //! \param mult the outgoing particle multiplicity data
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param transfer the transfer matrix
  void get_T( const dd_vector& sigma, const dd_vector& mult,
    const dd_vector& weight, T_matrix& transfer );

  // Prints the lists for debugging
  void print( );
};

#endif
