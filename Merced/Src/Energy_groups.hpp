/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2014-08-05 19:06:56 -0800 (Tue, 05 Aug 2014) $
 * $Author: hedstrom $
 * $Id: Energy_groups.hpp 1 2014-08-05 03:06:56Z hedstrom $
 *
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */
// the class for energy groups
#ifndef EN_GROUP_DEF
#define EN_GROUP_DEF

#include <vector>

#include "data_parser.hpp"

using namespace std;

//! Class for the energy group boundaries
//-----------class Energy_groups ----------
class Energy_groups : public vector< double >
{
public:
  inline Energy_groups( ) {}
  inline ~Energy_groups( ) {}

  //! Reads the energy group boundaries
  //! \param infile input file
  //! \param num_bd number of energy group boundaries
  void read_bd( data_parser &input_file, int num_bd );

  //! Returns the index of the left-hand end of the energy bin containing this_E
  //! \param this_E the given energy
  int first_bin_ID( double this_E ) const;

  //! Returns the index of the right-hand end of the energy bin containing this_E
  //! \param this_E the given energy
  int last_bin_ID( double this_E ) const;
};

#endif
