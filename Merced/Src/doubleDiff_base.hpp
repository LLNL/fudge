/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2021-10-05 $
 * $Author: hedstrom $
 * $Id: doubleDiff_base.hpp 1 2021-10-05Z hedstrom $
 *
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// Declaration of the basic classes for pointwise energy-angle
// probability density

#ifndef DOUBLE_DIFF_BASE
#define DOUBLE_DIFF_BASE

#include <list>

#include "dd_vector.hpp"

namespace DDbase
{
  // ----------- class DDbase::one_Eout -----------------
  //! Class for the energy distribution at given incident energy and angle
  class one_Eout : public Ddvec::dd_vector
  {
  private:

  public:
    // These data all range over -1 <= mu <= 1.
    // Ddvec::unit_base_map ubase_map;

    //! Default constructor
    inline one_Eout( )
    {}

    //! Default destructor
    inline ~one_Eout( )
    {}

    //! Use the direction cosine, Eout, as the tag
    inline double get_Eout( ) const
    {
      return get_tag( );
    }

    //! Sets the direction cosine tag
    //! \param Eout the direction cosine of the outgoing particle
    inline void set_Eout( double Eout )
    {
      set_tag( Eout );
    }

  };

  // ----------- class DDbase::one_joint_dist -----------------
  //! Class for the joint energy-angle distribution at given incident energy
  class one_joint_dist : public std::list< DDbase::one_Eout >
  {
  private:
    double tag_;

  public:
     Ddvec::unit_base_map Eout_ubase_map;

    //! Default constructor
    inline one_joint_dist( )
    {}

    //! Default destructor
    inline ~one_joint_dist( )
    {}

    //! Use E_in as the tag
    inline double get_E_in( ) const
    {
      return tag_;
    }

    //! Sets the tag
    //! \param E_in the energy of the incident particle
    inline void set_E_in( double E_in )
    {
      tag_ = E_in;
    }


    //! Checks the norm of the data
    void check_norm( );

    //! Maps the outgoing energies to [0, 1]
    void to_unit_base( );

  };

  // ----------- class DDbase::joint_dist_base -----------------
  //! Base class for joint energy-angle distributions
  class joint_dist_base : public std::list< DDbase::one_joint_dist >
  {
  private:

    //! Checks the norm of the data
    void check_norm( );

    //! Maps the outgoing energies to 0 <= Eout <= 1
    void to_unit_base( );

  public:
 
    Terp::two_d_interp Ein_interp; // interpolation rule for incident energy
    Terp::two_d_interp Eout_interp; // interpolation rule for outgoing energy
    Terp::Interp_Type mu_interp; // interpolation rule for direction cosine

    //! Default constructor
    joint_dist_base( ): mu_interp( Terp::LINLIN )
    {}

    //! Default destructor
    inline ~joint_dist_base( )
    {}

    //! Reads the python data
    //! \param infile input file
    //! \param num_Ein number of incident energies for this reaction
    void read_data( Dpar::data_parser &inFile, int num_Ein );

  };
  
}  // end of namespace DDbase

#endif

