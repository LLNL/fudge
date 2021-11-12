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
#include "quad_methods.hpp"
#include "Energy_groups.hpp"

namespace Trf
{
//! Class for the transfer matrix
//-----------class T_matrix ----------
class T_matrix
{
public:
  Coef::coef_vector *data;
  Coef::coef_vector *row_checks;  // used for verification
  Coef::coef_vector *row_sums;  // used for verification
  Lgdata::Flux_List e_flux;  // the Legendre coefficients of the incident flux
  int order;   // the Legendre order
  int num_Ein_bins;
  int num_Eout_bins;
  Coef::Conserve conserve;
  Qmeth::Quadrature_Rule Ein_quad_rule;  // quadrature rule for incident energy
  Qmeth::Quadrature_Rule Eout_quad_rule;  // quadrature rule for outgoing energy
  Qmeth::Quadrature_Rule mu_quad_rule;  // quadrature rule for outgoing cosine
  Qmeth::Quadrature_Rule mucm2_quad_rule;  // quadrature rule for second cosine, 2-step reaction
  Qmeth::Quadrature_Rule w_quad_rule;  // quadrature rule for w parameter, 2-step reaction

  //! The energy groups
  Egp::Energy_groups in_groups;
  Egp::Energy_groups out_groups;

  //! the weights 1 / \int_(energy group) flux by Legendre order
  Lgdata::weight_list flux_weight;

  //! Average cross section over energy bins, for checking the contribution of an interval
  std::vector< double > averageCrossSections;

  //! The maximum average cross section
  double maximumCrossSection;
  double threshold;  // reaction threshold
  std::ofstream *output_file;  // pointer to the output file

  //! Constructor
  T_matrix( );

  //! Destructor
  ~T_matrix( );

  //! Allocates space
  void allocate( );

  //! Operator to grab an entry:
  //! \param Ein_count: incident energy group
  //! \param Eout_count: exit energy group
  Coef::coef_vector& operator()( int Ein_count, int Eout_count );

  //! Operator to grab an entry:
  //! \param Ein_count: incident energy group
  //! \param Eout_count: exit energy group
  //! \param ell: the Legendre order of this term
  Coef::coef_vector& operator()( int Ein_count, int Eout_count, int ell );

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
  void getBinCrossSection( const Ddvec::dd_vector& sigma );

  //! Prints the matrix to the output file
  void write_transfer( );

  //! Prints zeros to the output file for a reaction with high threshold
  void zero_transfer( );
};

} // end of namespace Trf

#endif
