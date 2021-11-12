/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2018-11-19  (Mon., Nov. 19, 2018) $
 * $Author: hedstrom $
 * $Id: two_step_hit.hpp 1 2018-11-19 hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
//! Defines the classes used for the geometry of the 2-step reaction

#ifndef TWO_STEP_HIT_CLASS
#define TWO_STEP_HIT_CLASS

#include "box_geom.hpp"
#include "mappings.hpp"
#include "relativistic.hpp"

namespace Twohit
{

  //! Class for determination of incident energies which produce changes in geometry

//! Class for the list of intersections of an mu_cm_2 = const. surface with an E_out silo.
//! In this list the "bottom" of the box is *Eout_ptr = const. and the top is *(++Eout_ptr) = const.
// ---------------- class two_step_hit_list ------------------
class two_step_hit_list : public Box::hit_list
{
private:
  
  // **** Implement the virtual functions *****
  //! Calculates E_out from E_in, for testing the side of a box
  //! \param E_in energy of the incident particle
  double get_Eout( double E_in );

  //! Finds an intersection with the bottom or top of a box
  //! Returns incident energy corresponding to E_out
  //! \param E_out the level to hit
  //! \param pair_0 ( Ein, Eout ) at lower incident energy
  //! \param pair_1 ( Ein, Eout ) at higher incident energy
  double find_hit( double E_out, const Ddvec::dd_entry &pair_0,
                   const Ddvec::dd_entry &pair_1 );

  //! Finds the intersections with the bottom of a box
  //! \param E_out_lab bottom boundary of the outgoing energy bin
  //! \param Ein_hits computed incident energies which give intersections
  void find_bottom_hits( double E_out_lab, std::vector< Box::Ein_Eta_Hit > *Ein_hits );

  //! Finds the intersections with the top of a box
  //! \param E_out_lab top boundary of the outgoing energy bin
  //! \param Ein_hits computed incident energies which give intersections
  void find_top_hits( double E_out_lab, std::vector< Box::Ein_Eta_Hit > *Ein_hits );
  // **** end of virtual functions *****

public:
  //! data for determination of how this outgoing energy curve hits the quadrature box
  Ddvec::dd_entry flip;  // Ein for minimal Eout, minimal Eout for this mucm 
  Ddvec::dd_entry left_Ein_Eout;  // ( Ein, Eout ) for this mucm and for lower Ein
  Ddvec::dd_entry right_Ein_Eout;  // ( Ein, Eout ) for this mucm and for higher Ein
  
  bool use_relativistic;
  
  //! For the boost to the lab frame
  Maps::two_step_map_param twoStepMap;
  Rel::relativistic_2_step_param relTwoStepParam;
  
  inline two_step_hit_list( ) : use_relativistic( false ) {}
  inline ~two_step_hit_list( ) {}

  //! Sets the direction cosine for step 1
  //! \param mucm1, the direction cosine for step 1
  void set_mucm_1( double mucm1 );

  //! Sets the direction cosine for step 2
  //! \param mucm2, the direction cosine for step 2
  void set_mucm_2( double mucm2 );

};

//! Class for an array of Twohit::two_step_hit_lists
// ---------------- class two_step_hit_list_array ------------------
class two_step_hit_list_array
{
private:
  int laengd;

public:
  Twohit::two_step_hit_list *hit_lists;

  two_step_hit_list_array( ): laengd( 0 ) {}

  ~two_step_hit_list_array( )
  {
    if( laengd > 0 ) delete [] hit_lists;
  }

  //! allocates space
  //! \param Laengd, the length of the hit_lists array
  void get_array( int Laengd );

  //! Ensures that these hit_lists are defined at the same incident energies.
  void common_hits( );
};

} // end of namespace Twohit

#endif


