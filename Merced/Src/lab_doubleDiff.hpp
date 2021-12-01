/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2021-09-29 $
 * $Author: hedstrom $
 * $Id: lab_doubleDiff.hpp 1 2021-09-21Z hedstrom $
 *
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// declaration of the classes used to handle pointwise energy-angle
// probability density given in the lab frame.

// We read in the data as for center-of-mass data and
// convert the angular tables to Legendre coefficients.

#ifndef LAB_DOUBLE_DIFF
#define LAB_DOUBLE_DIFF

#include "doubleDiff_base.hpp"  // use its classes as input
#include "standard_Legendre.hpp"  // process as Legendre data

namespace labDD
{
  // ----------- class labDD::lab_joint_dist -----------------
  //! Class for joint energy-angle distributions
  class lab_joint_dist : public DDbase::joint_dist_base
  {
  private:
    // use Legendre data
    StdLg::standard_Legendre Legendre_data;

    //! converts the angular table to Legendre coefficients
    void to_Legendre( );

  public:
    
    //! Default constructor
    lab_joint_dist( )
    {}

    //! Default destructor
    inline ~lab_joint_dist( )
    {}

    //! Calculates the transfer matrix for this particle.
    //! \param sigma the cross section data
    //! \param mult the outgoing particle multiplicity data
    //! \param weight the weighting to apply to the transfer matrix entries
    //! \param transfer the transfer matrix
    void get_T( const Ddvec::dd_vector& sigma, const Ddvec::dd_vector& mult,
                const Ddvec::dd_vector& weight, Trf::T_matrix& transfer );

  };

}  // end of namespace labDD

#endif
