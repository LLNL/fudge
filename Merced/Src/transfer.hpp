/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: transfer.hpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */
// the class for transfer matrices
#ifndef T_MATRIX_DEF
#define T_MATRIX_DEF

#include <vector>
#include <string>
#include <list>
#include <fstream>

#include "coef_vector.hpp"
#include "Legendre_data.hpp"
#include "data_parser.hpp"
#include "quadrature.hpp"
#include "Energy_groups.hpp"

using namespace std;

//! Class for the transfer matrix
//-----------class T_matrix ----------
class T_matrix
{
public:
  coef_vector *data;
  coef_vector *row_checks;  // used for verification
  coef_vector *row_sums;  // used for verification
  Flux_List e_flux;  // the Legendre coefficients of the incident flux
  int order;   // the Legendre order
  int num_Ein_bins;
  int num_Eout_bins;
  Conserve conserve;
  Quadrature_Method Ein_quad_method;  // quadrature method for incident energy
  Quadrature_Method Eout_quad_method;  // quadrature method for outgoing energy
  Quadrature_Method mu_quad_method;  // quadrature method for outgoing cosine

  bool interpolate_Eout_integrals;  // do we interpolate the integrals over Eout?

  //! The energy groups
  Energy_groups in_groups;
  Energy_groups out_groups;

  //! the weights 1 / \int_(energy group) flux by Legendre order
  weight_list flux_weight;

  //! Average cross section over energy bins, for checking the contribution of an interval
  vector< double > averageCrossSections;

  //! The maximum average cross section
  double maximumCrossSection;
  ofstream *output_file;  // pointer to the output file

  //! Constructor
  T_matrix( );

  //! Destructor
  ~T_matrix( );

  //! Allocates space
  void allocate( );

  //! Operator to grab an entry:
  //! \param Ein_count: incident energy group
  //! \param Eout_count: exit energy group
  coef_vector& operator()( int Ein_count, int Eout_count );

  //! Operator to grab an entry:
  //! \param Ein_count: incident energy group
  //! \param Eout_count: exit energy group
  //! \param ell: the Legendre order of this term
  coef_vector& operator()( int Ein_count, int Eout_count, int ell );

  //! Calculates the weights 1 / \int_(energy group) flux
  void get_flux_weight( );

  //! Applies the weights 1 / \int_(energy group) flux
  void use_weight( );

  //! Scales the weight_E terms by the average E_out, to check gamma output
  void scale_E( );

  // Prints the sums of the rows of the zero-order term.
  // They should agree with the flux-weighted average cross sections
  void check_ell0( );

  // Scales the row check sums by the weights 
  void scale_row_check( );

  //! Computes the net cross section for each energy bin
  //! \param sigma a vector of pairs ( incident energy, cross section )
  void getBinCrossSection( const dd_vector& sigma );

  //! Prints the matrix to the output file
  void write_transfer( );

  //! Prints zeros to the output file for a reaction with high threshold
  void zero_transfer( );
};

#endif
