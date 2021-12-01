/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2009-09-15  (Tue., Sept. 15, 2009) $
 * $Author: hedstrom $
 * $Id: Vcm_Vlab_Hit.hpp 1 2010-10-05 hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
//! Defines the classes used for the geometry of the cm-to-lab map

#ifndef VCM_VLAB_HIT_CLASS
#define VCM_VLAB_HIT_CLASS

#include "box_geom.hpp"
#include "param_base.hpp"
#include "mappings.hpp"

namespace Vhit
{
//! How does a curve E_cm = const intersect an Eout_lab-mu_lab quadrature box?
// ---------------- Hit_Corner -----------------------------
enum Hit_Corner{ V_INSIDE,        // lies inside the box
                 BOTTOM_FORWARD,  // hits the bottom at mu = 1
                 BOTTOM_BACKWARD, // hits the bottom at mu = -1
                 TOP_FORWARD,     // hits the top at mu = 1
                 TOP_BACKWARD,    // hits the top at mu = -1
                 V_BELOW,         // lies below the box
                 V_ABOVE};        // lies above the box

  
//! Class for an intersection of E_cm = const with an Eout_lab-mu_lab quadrature box
// ---------------- class Vcm_quadBox_Hit ------------------
class Vcm_quadBox_Hit
{
public:
  double V_cm;
  Hit_Corner hit_corner;

  //! default constructor
  inline Vcm_quadBox_Hit( ){}

  //! copy constructor
  //! \param Vcm_quadBox_hit the data to copy
  inline Vcm_quadBox_Hit( const Vcm_quadBox_Hit& Vcm_quadBox_hit )
  {
    V_cm = Vcm_quadBox_hit.V_cm;
    hit_corner = Vcm_quadBox_hit.hit_corner;
  }
  inline ~Vcm_quadBox_Hit( ){}

  //! Identifies the relative locations
  //! \param V_cm center-of-mass velocity of the outgoing particle
  //! \param V_trans velocity of the center of mass in the lab frame
  //! \param V_lab_min minimum velocity for this outgoing energy bin
  //! \param V_lab_max maximum velocity for this outgoing energy bin
  void set_region( double V_cm, double V_trans, double V_lab_min,
                   double V_lab_max );

  //! Prints a hit
  void print( ) const;
};
}    // end of this part of namespace Vhit

namespace Vcm_quadBox_Hit_F
{
  // ---------------- Vcm_quadBox_Hit comparison function ------------------
  //! Order Vcm_quadBox_Hit by its V_cm member
  //! \param first the first Vhit::Vcm_quadBox_Hit to compare
  //! \param second the  second Vhit::Vcm_quadBox_Hit to compare
  bool lessthan_F( Vhit::Vcm_quadBox_Hit first, Vhit::Vcm_quadBox_Hit second );
}

namespace Vhit // more
{
  
//! Class for the list of intersections of an E_cm = const. surface with an E_out silo.
//! In this list the "bottom" of the box is *Eout_ptr = const. and the top is *(++Eout_ptr) = const.
//! The "left side" is mu_lab = -1 and the "right side" is mu_lab = 1.
// ---------------- class Vcm_Vlab_hit_list ------------------
class Vcm_Vlab_hit_list : public Box::hit_list_base
{
private:
  //! The desired outgoing energy in the laboratory frame 
  double E_lab;

  //! Sets the derivative of hit_G with respect to alpha, 
  //! where Ein = (1 - alpha)*Ein_0 + alpha*Ein_1.
  //! \param E_lab the desired outgoing energy in the lab frame
  void set_d_hit_G( double E_lab );
  
  //! Solves hit_G = 0 using Taylor series about Ein_0; returns the number of real roots
  //! \param ein_1 first incident energy giving the desired outgoing energy
  //! \param ein_2 second incident energy giving the desired outgoing energy
  int solve_hit_G0( double *ein_1, double *ein_2 );

  //! Solves hit_G = 0 using Taylor series about Ein_1; returns the number of real roots
  //! \param ein_1 first incident energy giving the desired outgoing energy
  //! \param ein_2 second incident energy giving the desired outgoing energy
  int solve_hit_G1( double *ein_1, double *ein_2 );

  //! Checks the G value for the computed incident energies giving intersections
  //! \param number of incident energies giving the desired outgoing energy
  //! \param ein_1 first incident energy giving the desired outgoing energy
  //! \param ein_2 second incident energy giving the desired outgoing energy
  void do_check_roots( int num_roots, double ein_1, double ein_2 );

  // **** Implement the virtual functions *****
  //! function not used
  //! \param E_in energy of the incident particle
  double get_Eout( double E_in ) {return 0.0;}

  //! Finds the intersections with the bottom of a box
  //! \param E_out_lab bottom boundary of the outgoing energy bin
  //! \param Ein_hits computed incident energies which give intersections
  void find_bottom_hits( double E_out_lab, std::vector< Box::Ein_Eta_Hit > *Ein_hits );

  //! Finds the intersections with the top of a box
  //! \param E_out_lab top boundary of the outgoing energy bin
  //! \param Ein_hits computed incident energies which give intersections
  void find_top_hits( double E_out_lab, std::vector< Box::Ein_Eta_Hit > *Ein_hits );

  //! Where do we hit the left-hand side of the box?
  //! \param dummy not used
  //! \param Eout_ptr bottom boundary of the outgoing energy bin
  void test_left( double dummy, std::vector< double >::const_iterator Eout_ptr );

  //! Where do we hit the right-hand side of the box?
  //! \param dummy not used
  //! \param Eout_ptr bottom boundary of the outgoing energy bin
  void test_right( double dummy, std::vector< double >::const_iterator Eout_ptr );

  //! No top or bottom intersections, so are we above, below, or through the middle?
  //! \param Eout_ptr bottom boundary of the outgoing energy bin
  //! \param dummy_left not used
  //! \param dummy_right not used
  void test_sides( std::vector< double >::const_iterator Eout_ptr,
    double dummy_left, double dummy_right );

  //! Checks for inconsistencies caused by peculiarities of real arithmetic.
  //! Returns is_OK == true if there are no inconsistencies.
  bool consistency( );
  // **** end of virtual functions *****

public:
  //! parameters used by the hit_G function and its derivative
  Maps::Ecm_intersect G0_data;  // at the lower incident energy
  Maps::Ecm_intersect G1_data;  // at the higher incident energy

  inline Vcm_Vlab_hit_list( ) {}
  inline ~Vcm_Vlab_hit_list( ) {}

  //! Calculates hit_G and d_hit_G for given E_out_lab
  //! param E_out_lab desired outgoing energy
  void set_G_data( double E_out_lab );
};

}  // end of namespace Vhit

#endif
