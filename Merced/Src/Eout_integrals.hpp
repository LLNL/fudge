/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2009-08-11 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: Eout_integrals.hpp 1  2009-08-11 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 * 
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */

// header for the vector of integrals over Eout/mu
#ifndef EOUT_INTEGRALS
#define EOUT_INTEGRALS

#include <list>

#include "coef_vector.hpp"

using namespace std;

// --------------- Eout_int_param ------------------------
//! functions needed to compute the integrals over Eout and mu
class Eout_int_param
{
public:
  // The range of integration
  double Eout_0;
  double Eout_1;

  inline Eout_int_param( ) {}
  virtual ~Eout_int_param( ) {}

// ************** virtual routines ******************************
  //! Initialize at a given incident energy
  //! \param E_in energy of the incident particle
  virtual void set_Ein( double E_in ) = 0;

  //! Evaluate the integrals over Eout
  //! \param Eout_0 lower outgoing energy limit of integration 
  //! \param Eout_1 upper outgoing energy limit of integration 
  //! \param value computed value of the integrals
  virtual void get_integrals( double Eout_0, double Eout_1, coef_vector &value ) = 0;

  //! Evaluate the integrals over Eout and return the noise in the calculation
  //! \param Eout_0 lower outgoing energy limit of integration 
  //! \param Eout_1 upper outgoing energy limit of integration 
  //! \param value computed value of the integrals
  virtual double tol_get_integrals( double Eout_0, double Eout_1, coef_vector &value ) = 0;
// ************** end of virtual routines ***********************
};

// --------------- Eout_link ---------------------------
//! Integrals over mu and one E_out bin for one incident energy
class Eout_link: public coef_vector
{
public:
  double E_in;

  //! Constructor
  inline Eout_link( ) {}

  //! Destructor
  inline ~Eout_link( ) {}

  //! Linear interpolation between this link and the next
  //! \param Ein the energy to interpolate to
  //! \param next_link data for the next link
  //! \param interp the interpolated data
  void Interpolate( double Ein, const Eout_link& next_link,
		    coef_vector *interp ) const;
};

// --------------- Eout_integrals ---------------------------
class Eout_integrals: public list< Eout_link >
{
public:
  Eout_int_param *Eout_int_params;

  //! Constructor
  inline Eout_integrals( ) {}

  //! Destructor
  inline ~Eout_integrals( ) {}

  //! Adds a new link to the Eout_ints list
  //! \param where the new link is inserted before this one
  //! \param order the Legendre order of the link
  //! \param conserve flag to conserve energy or particle number or both
  //! \param E_in energy of the incident particle
  void new_Eout_int( Eout_integrals::iterator where, int order, 
    Conserve conserve, double E_in );

  //! Adds a new link to the Eout_ints list and returns the noise in the calculation
  //! \param where the new link is inserted before this one
  //! \param order the Legendre order of the link
  //! \param conserve flag to conserve energy or particle number or both
  //! \param E_in energy of the incident particle
  double tol_new_Eout_int( Eout_integrals::iterator where, int order,
    Conserve conserve, double E_in );

  //! Sets up the integrals over E_out for interpolation
  //! \param order the Legendre order of the link
  //! \param conserve flag to conserve energy or particle number or both
  //! \param Ein_0 lower energy of the incident particle
  //! \param Ein_1 higher energy of the incident particle
  void setup_Eout_ints( int order, Conserve conserve, double Ein_0, double Ein_1 );

  //! Check the accuracy of linear interpolation
  //! \param interp_link a link at an intermediate outgoing energy
  //! \param prev_link a link at a lower outgoing energy
  //! \param next_link a link at a higher outgoing energy
  //! \param tol the required tolerance for linear interpolation
  bool interp_OK( const Eout_link& interp_link,
		  Eout_integrals::iterator prev_link,
		  Eout_integrals::iterator next_link, double tol );
};

#endif
