/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: box_geom.hpp 1 2006-09-18 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */
// classes used to keep track of the regions of integration

#ifndef BOX_GEOM_DEF
#define BOX_GEOM_DEF

#include <list>
#include <vector>
#include "dd_vector.hpp"  // for Ddvec::dd_pair

namespace Box
{
//! How does a curve eta = const intersect an E-E' quadrature box?
//! center-of-mass data needs the last 4 options later.
// ---------------- Hit_Edge -----------------------------
enum Hit_Edge{ MISSED,     // ENTRY NOT SET
               INSIDE,     // intersects a side or is inside the box
               BOTTOM_IN,  // enters the bottom
               BOTTOM_OUT, // exits the bottom
               TOP_IN,     // enters the top
               TOP_OUT,    // exits the top
               BELOW,      // lies below the box
               ABOVE,      // lies above the box
               ABOVE_FORWARD,     // forward emission above the box, cm data
               BOTTOM_CORNER_IN,  // move from hitting the bottom to inside the box
               BOTTOM_CORNER_OUT, // move from inside the box to hitting the bottom
               TOP_CORNER_IN,     // move from hitting the top to inside the box
               TOP_CORNER_OUT};   // move from inside the box to hitting the top

//! Class for an intersection of eta = const with a quadrature box
// ---------------- class Ein_Eta_Hit ------------------
class Ein_Eta_Hit
{
public:
  double E_in;
  Box::Hit_Edge hit_edge;

  Ein_Eta_Hit( ): E_in(0.0), hit_edge( Box::MISSED ) {}
  Ein_Eta_Hit( const Ein_Eta_Hit& ein_eta_hit );
  ~Ein_Eta_Hit( ){}

  //! Prints a hit
  void print( ) const;

  //! Return the name of the edge.
  std::string edge_name( ) const ;
};

//! General class also used by models in center-of-mass coordinates
//! Class for all intersections of eta = const with a quadrature box
// ---------------- class hit_list_base ------------------
class hit_list_base : public std::list< Box::Ein_Eta_Hit >
{
private:

protected:
  double box_tol;  // When to ignore intersection with the top or bottom of a box

  //! Inserts a link in proper order
  //! \param ein_eta_hit pair: incident energy, value of Box::Hit_Edge
  void stick_in( const Box::Ein_Eta_Hit& ein_eta_hit );

  //! Finds intersections of the curve eta = const with the bottom of the E-E' box
  //! \param Eout_ptr desired lower energy of the outgoing particle
  //! \param E_in_left lower incident energy
  //! \param E_in_right higher incident energy
  void test_bottom( std::vector< double >::const_iterator Eout_ptr,
    double E_in_left, double E_in_right );

  //! Finds intersections of the curve eta = const with the top of the E-E' box
  //! \param Eout_ptr desired higher energy of the outgoing particle
  //! \param E_in_left lower incident energy
  //! \param E_in_right higher incident energy
  void test_top( std::vector< double >::const_iterator Eout_ptr,
    double E_in_left, double E_in_right );

  //! Default check for inconsistencies caused by peculiarities of real arithmetic.
  //! Returns is_OK == true if there are no inconsistencies.
  bool consistency_default( );

  //! Where do we hit the left-hand side of the box?
  //! \param Eout_ptr desired lower energy of the outgoing particle
  //! \param E_in lower incident energy
  void test_left_default( double E_in, std::vector< double >::const_iterator Eout_ptr );

  //! Where do we hit the right-hand side of the box?
  //! \param Eout_ptr desired lower energy of the outgoing particle
  //! \param E_in higher incident energy
  void test_right_default( double E_in, std::vector< double >::const_iterator Eout_ptr );

  // ****** Virtual routines depending on the model *********
  //! Where do we hit the left-hand side of the box?
  //! \param Eout_ptr desired lower energy of the outgoing particle
  //! \param E_in lower incident energy
  virtual void test_left( double E_in, std::vector< double >::const_iterator Eout_ptr ) = 0;

  //! Where do we hit the right-hand side of the box?
  //! \param Eout_ptr desired lower energy of the outgoing particle
  //! \param E_in higher incident energy
  virtual void test_right( double E_in, std::vector< double >::const_iterator Eout_ptr ) = 0;

  //! No top or bottom intersections, so are we above, below, or through the middle?
  //! \param Eout_ptr desired lower energy of the outgoing particle
  //! \param E_in_left lower incident energy
  //! \param E_in_right higher incident energy
  virtual void test_sides( std::vector< double >::const_iterator Eout_ptr,
    double E_in_left, double E_in_right ) = 0;

  //! Checks for inconsistencies caused by peculiarities of real arithmetic.
  //! Returns is_OK == true if there are no inconsistencies.
  virtual bool consistency( ) = 0;

  //! Finds the intersections with the bottom of a box
  //! \param E_out desired lower energy of the outgoing particle
  //! \param Ein_hits list of pairs ( incident energy internal to the box, Box::Hit_Edge )
  virtual void find_bottom_hits( double E_out,
    std::vector< Box::Ein_Eta_Hit > *Ein_hits ) = 0;

  //! Finds the intersections with the top of a box
  //! \param E_out desired higher energy of the outgoing particle
  //! \param Ein_hits list of pairs ( incident energy internal to the box, Box::Hit_Edge )
  virtual void find_top_hits( double E_out,
    std::vector< Box::Ein_Eta_Hit > *Ein_hits ) = 0;
  // ****** End of virtual routines depending on the model ******

 public:
  double eta;

  hit_list_base( )  {}
  virtual ~hit_list_base( ){}

  //! Finds intersections of the curve eta = const with the E-E' box.
  //! Returns is_OK == true if there are no inconsistencies.
  //! \param Eout_ptr desired lower energy of the outgoing particle
  //! \param E_in_left lower incident energy
  //! \param E_in_right higher incident energy
  bool hit_box( double Eta, std::vector< double >::const_iterator Eout_ptr,
    double E_in_left, double E_in_right );

  //! Checks whether this eta = const curve is below the E-E' box
  bool is_below( );

  //! Checks whether this eta = const curve is above the E-E' box
  bool is_above( );

  //! Fills in the list with energies from Ein_hits
  //! \param Ein_hits a vector of incident energies to be inserted
  void fill_in( const std::vector< double >& Ein_hits );

  //! Ensures that these hits and upper_hits are defined at the same incident energies.
  //! \param upper_hits a second list of pairs (incident energy, Box::Hit_Edge)
  void common_hits( hit_list_base &upper_hits );

  //! Prints the list
  void print( );
};

//! Class for all intersections of eta = const with a quadrature box
// ---------------- class hit_list ------------------
class hit_list : public Box::hit_list_base
{
private:

protected:
  // ********** implement some virtual functions *******************
  //! Where do we hit the left-hand side of the box?
  //! \param E_in lower incident energy
  //! \param Eout_ptr desired lower energy of the outgoing particle
  void test_left( double E_in, std::vector< double >::const_iterator Eout_ptr )
  { test_left_default( E_in, Eout_ptr ); }

  //! Where do we hit the right-hand side of the box?
  //! \param E_in higher incident energy
  //! \param Eout_ptr desired lower energy of the outgoing particle
  void test_right( double E_in, std::vector< double >::const_iterator Eout_ptr )
  { test_right_default( E_in, Eout_ptr ); }

  //! No top or bottom intersections, so are we above, below, or through the middle?
  //! \param Eout_ptr desired lower energy of the outgoing particle
  //! \param E_in_left lower incident energy
  //! \param E_in_right higher incident energy
  void test_sides( std::vector< double >::const_iterator Eout_ptr,
		   double E_in_left, double E_in_right );

  //! Checks for inconsistencies caused by peculiarities of real arithmetic.
  //! Returns is_OK == true if there are no inconsistencies.
  inline bool consistency( )
  { return consistency_default( ); }
  // ********* end of implemented virtual functions *******************

  // ****** new virtual routines depending on the model *********
  //! Calculates E_out from E_in, for testing the side of a box
  //! \param E_in the incident energy
  virtual double get_Eout( double E_in ) = 0;
  // ****** End of new virtual routines depending on the model ******

  // ******** Box::hit_list_base virtual functions **************
  //! Finds the intersections with the bottom of a box
  //! \param E_out desired lower energy of the outgoing particle
  //! \param Ein_hits list of pairs ( incident energy internal to the box, Box::Hit_Edge )
  virtual void find_bottom_hits( double E_out,
     std::vector< Box::Ein_Eta_Hit > *Ein_hits ) {}

  //! Finds the intersections with the top of a box
  //! \param E_out desired higher energy of the outgoing particle
  //! \param Ein_hits list of pairs ( incident energy internal to the box, Box::Hit_Edge )
  virtual void find_top_hits( double E_out,
    std::vector< Box::Ein_Eta_Hit > *Ein_hits ) {}
  // ****** End of Box::hit_list_base virtual routines ******

 public:
  hit_list( )  {}
  virtual ~hit_list( ){}

};

//! Class for the list of intersections of a unit-base curve with an integration box
// ---------------- class energy_hit_list ------------------
class energy_hit_list : public Box::hit_list
{
private:
  // Implement the virtual functions
  //! Calculates E_out from E_in, for testing the side of a box
  //! \param E_in energy of the incident particle
  double get_Eout( double E_in );

  //! Finds the intersections with the bottom of a box
  //! \param E_out energy at the bottom of the outgoing energy bin
  //! \param Ein_hits computed incident energies giving intersections with E_out
  void find_bottom_hits( double E_out, std::vector< Box::Ein_Eta_Hit > *Ein_hits );

  //! Finds the intersections with the top of a box
  //! \param E_out energy at the top of the outgoing energy bin
  //! \param Ein_hits computed incident energies giving intersections with E_out
  void find_top_hits( double E_out, std::vector< Box::Ein_Eta_Hit > *Ein_hits );

public:
  //! A pair of values of (E, E') for this eta value; to determine a line.
  //! We want to know where this line has the desired value, E_out.
  Ddvec::dd_pair E_Eout;

  energy_hit_list( ) {}
  ~energy_hit_list( ) {}
};

} // end of namespace Box

#endif
