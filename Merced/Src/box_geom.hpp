/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: box_geom.hpp 1 2006-09-18 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
  Copyright (c) 2017, Lawrence Livermore National Security, LLC.
  Produced at the Lawrence Livermore National Laboratory.
  Written by the LLNL Nuclear Data and Theory group
          (email: mattoon1@llnl.gov)
  LLNL-CODE-725546.
  All rights reserved.
  
  This file is part of the Merced package, used to generate nuclear reaction
  transfer matrices for deterministic radiation transport.
  
  
      Please also read this link - Our Notice and Modified BSD License
  
  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:
      * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the disclaimer below.
      * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the disclaimer (as noted below) in the
        documentation and/or other materials provided with the distribution.
      * Neither the name of LLNS/LLNL nor the names of its contributors may be used
        to endorse or promote products derived from this software without specific
        prior written permission.
  
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
  THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  
  
  Additional BSD Notice
  
  1. This notice is required to be provided under our contract with the U.S.
  Department of Energy (DOE). This work was produced at Lawrence Livermore
  National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
  
  2. Neither the United States Government nor Lawrence Livermore National Security,
  LLC nor any of their employees, makes any warranty, express or implied, or assumes
  any liability or responsibility for the accuracy, completeness, or usefulness of any
  information, apparatus, product, or process disclosed, or represents that its use
  would not infringe privately-owned rights.
  
  3. Also, reference herein to any specific commercial products, process, or services
  by trade name, trademark, manufacturer or otherwise does not necessarily constitute
  or imply its endorsement, recommendation, or favoring by the United States Government
  or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
  herein do not necessarily state or reflect those of the United States Government or
  Lawrence Livermore National Security, LLC, and shall not be used for advertising or
  product endorsement purposes.
  
 * # <<END-copyright>>
 */
// classes used to keep track of the regions of integration

#ifndef BOX_GEOM_DEF
#define BOX_GEOM_DEF

#include <list>
#include <vector>
#include "dd_vector.hpp"  // for dd_pair

using namespace std;

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
  Hit_Edge hit_edge;

  Ein_Eta_Hit( ): E_in(0.0), hit_edge(MISSED) {}
  Ein_Eta_Hit( const Ein_Eta_Hit& ein_eta_hit );
  ~Ein_Eta_Hit( ){}

  //! Prints a hit
  void print( ) const;
};

//! General class also used by models in center-of-mass coordinates
//! Class for all intersections of eta = const with a quadrature box
// ---------------- class hit_list_base ------------------
class hit_list_base : public list< Ein_Eta_Hit >
{
private:

protected:
  double box_tol;  // When to ignore intersection with the top or bottom of a box

  //! Inserts a link in proper order
  //! \param ein_eta_hit pair: incident energy, value of Hit_Edge
  void stick_in( const Ein_Eta_Hit& ein_eta_hit );

  //! Fills in the list with energies from Ein_hits
  //! \param Ein_hits a vector of incident energies to be inserted
  void fill_in( const vector< double >& Ein_hits );

  //! Finds intersections of the curve eta = const with the bottom of the E-E' box
  //! \param Eout_ptr desired lower energy of the outgoing particle
  //! \param E_in_left lower incident energy
  //! \param E_in_right higher incident energy
  void test_bottom( vector< double >::const_iterator Eout_ptr,
    double E_in_left, double E_in_right );

  //! Finds intersections of the curve eta = const with the top of the E-E' box
  //! \param Eout_ptr desired higher energy of the outgoing particle
  //! \param E_in_left lower incident energy
  //! \param E_in_right higher incident energy
  void test_top( vector< double >::const_iterator Eout_ptr,
    double E_in_left, double E_in_right );

  //! Default check for inconsistencies caused by peculiarities of real arithmetic.
  //! Returns is_OK == true if there are no inconsistencies.
  bool consistency_default( );

  //! Where do we hit the left-hand side of the box?
  //! \param Eout_ptr desired lower energy of the outgoing particle
  //! \param E_in lower incident energy
  void test_left_default( double E_in, vector< double >::const_iterator Eout_ptr );

  //! Where do we hit the right-hand side of the box?
  //! \param Eout_ptr desired lower energy of the outgoing particle
  //! \param E_in higher incident energy
  void test_right_default( double E_in, vector< double >::const_iterator Eout_ptr );

  // ****** Virtual routines depending on the model *********
  //! Where do we hit the left-hand side of the box?
  //! \param Eout_ptr desired lower energy of the outgoing particle
  //! \param E_in lower incident energy
  virtual void test_left( double E_in, vector< double >::const_iterator Eout_ptr ) = 0;

  //! Where do we hit the right-hand side of the box?
  //! \param Eout_ptr desired lower energy of the outgoing particle
  //! \param E_in higher incident energy
  virtual void test_right( double E_in, vector< double >::const_iterator Eout_ptr ) = 0;

  //! No top or bottom intersections, so are we above, below, or through the middle?
  //! \param Eout_ptr desired lower energy of the outgoing particle
  //! \param E_in_left lower incident energy
  //! \param E_in_right higher incident energy
  virtual void test_sides( vector< double >::const_iterator Eout_ptr,
    double E_in_left, double E_in_right ) = 0;

  //! Checks for inconsistencies caused by peculiarities of real arithmetic.
  //! Returns is_OK == true if there are no inconsistencies.
  virtual bool consistency( ) = 0;

  //! Finds the intersections with the bottom of a box
  //! \param E_out desired lower energy of the outgoing particle
  //! \param Ein_hits list of pairs ( incident energy internal to the box, Hit_Edge )
  virtual void find_bottom_hits( double E_out,
    vector< Ein_Eta_Hit > *Ein_hits ) = 0;

  //! Finds the intersections with the top of a box
  //! \param E_out desired higher energy of the outgoing particle
  //! \param Ein_hits list of pairs ( incident energy internal to the box, Hit_Edge )
  virtual void find_top_hits( double E_out,
    vector< Ein_Eta_Hit > *Ein_hits ) = 0;
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
  bool hit_box( double Eta, vector< double >::const_iterator Eout_ptr,
    double E_in_left, double E_in_right );

  //! Checks whether this eta = const curve is below the E-E' box
  bool is_below( );

  //! Checks whether this eta = const curve is above the E-E' box
  bool is_above( );

  //! Ensures that these hits and upper_hits are defined at the same incident energies.
  //! \param upper_hits a second list of pairs (incident energy, Hit_Edge)
  void common_hits( hit_list_base &upper_hits );

  //! Prints the list
  void print( );
};

//! Class for all intersections of eta = const with a quadrature box
// ---------------- class hit_list ------------------
class hit_list : public hit_list_base
{
private:

protected:
  // ********** implement some virtual functions *******************
  //! Where do we hit the left-hand side of the box?
  //! \param E_in lower incident energy
  //! \param Eout_ptr desired lower energy of the outgoing particle
  void test_left( double E_in, vector< double >::const_iterator Eout_ptr )
  { test_left_default( E_in, Eout_ptr ); }

  //! Where do we hit the right-hand side of the box?
  //! \param E_in higher incident energy
  //! \param Eout_ptr desired lower energy of the outgoing particle
  void test_right( double E_in, vector< double >::const_iterator Eout_ptr )
  { test_right_default( E_in, Eout_ptr ); }

  //! No top or bottom intersections, so are we above, below, or through the middle?
  //! \param Eout_ptr desired lower energy of the outgoing particle
  //! \param E_in_left lower incident energy
  //! \param E_in_right higher incident energy
  void test_sides( vector< double >::const_iterator Eout_ptr,
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

  // ******** hit_list_base virtual functions **************
  //! Finds the intersections with the bottom of a box
  //! \param E_out desired lower energy of the outgoing particle
  //! \param Ein_hits list of pairs ( incident energy internal to the box, Hit_Edge )
  virtual void find_bottom_hits( double E_out,
     vector< Ein_Eta_Hit > *Ein_hits ) {}

  //! Finds the intersections with the top of a box
  //! \param E_out desired higher energy of the outgoing particle
  //! \param Ein_hits list of pairs ( incident energy internal to the box, Hit_Edge )
  virtual void find_top_hits( double E_out,
    vector< Ein_Eta_Hit > *Ein_hits ) {}
  // ****** End of hit_list_base virtual routines ******

 public:
  hit_list( )  {}
  virtual ~hit_list( ){}

};

//! Class for the list of intersections of a unit-base curve with an integration box
// ---------------- class energy_hit_list ------------------
class energy_hit_list : public hit_list
{
private:
  // Implement the virtual functions
  //! Calculates E_out from E_in, for testing the side of a box
  //! \param E_in energy of the incident particle
  double get_Eout( double E_in );

  //! Finds the intersections with the bottom of a box
  //! \param E_out energy at the bottom of the outgoing energy bin
  //! \param Ein_hits computed incident energies giving intersections with E_out
  void find_bottom_hits( double E_out, vector< Ein_Eta_Hit > *Ein_hits );

  //! Finds the intersections with the top of a box
  //! \param E_out energy at the top of the outgoing energy bin
  //! \param Ein_hits computed incident energies giving intersections with E_out
  void find_top_hits( double E_out, vector< Ein_Eta_Hit > *Ein_hits );

public:
  //! A pair of values of (E, E') for this eta value; to determine a line.
  //! We want to know where this line has the desired value, E_out.
  dd_pair E_Eout;

  energy_hit_list( ) {}
  ~energy_hit_list( ) {}
};

#endif
