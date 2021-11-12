/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2016-02-12 11:30:00 -0800 (Fri, 12 Feb 2016) $
 * $Author: hedstrom $
 * $Id: isotropic.hpp 1 2016-02-12 11:30:00Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// define the classes for isotropic distributions in the lab frame.

#ifndef ISOTROPIC_DEF
#define ISOTROPIC_DEF

#include "energy_dist.hpp"

namespace Iso
{
//! Class for isotropic distributions in the lab frame
// ----------------------- class isotropic --------------------------
class isotropic : public Edist::energy_moments
{
private:

public:
  isotropic( ) {}
  ~isotropic( ) {}

  //! Reads the data
  //! \param infile input file
  //! \param num_Ein number of incident energies for this reaction
  void read_data( Dpar::data_parser& infile, int num_Ein );

  // Calculates the transfer matrix for this particle.
  //! \param sigma the cross section data
  //! \param mult the outgoing particle multiplicity data
  //! \param weight the weighting to apply to the transfer matrix entries
  //! \param transfer the transfer matrix
  void get_T( const Ddvec::dd_vector& sigma, const Ddvec::dd_vector& mult,
	      const Ddvec::dd_vector& weight,
    Trf::T_matrix& transfer );
};

} // end of namespace Iso
#endif
