/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id:param_base.hpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/

//header for param_base, the base class for the quadrature parameters

#ifndef PARAM_BASE_CLASS
#define PARAM_BASE_CLASS

#include "quad_param.hpp"
#include "dd_vector.hpp"
#include "Legendre_data.hpp"
#include "Energy_groups.hpp"

namespace Pbase
{

// ----------- class param_base ------------------
//! Base lass for the parameters for integration over incident energy
class param_base: public Qparam::QuadParamBase
{
private:
  Terp::Interp_Type weight_interp;  // interpolation for the model weight
  
public:
  double data_E_0;     // the lower incident energy for the eta ladder
  double data_E_1;     // the upper incident energy for the eta ladder
  double Ein_0;        // the lower incident energy for quadrature
  double Ein_1;        // the upper incident energy for quadrature
  double Eout_min;     // the bottom of the E-E' quadrature box
  double Eout_max;     // the top of the E-E' quadrature box
  bool use_Eout_min;   // if true, get minimum eta from Eout_min
                       // otherwise get minimum eta from the lower data point
  bool use_Eout_max;   // if true, get maximum eta from Eout_max
                       // otherwise get maximum eta from the upper data point
  
  // pointers to the cross section
  Ddvec::dd_vector::const_iterator this_sigma;
  Ddvec::dd_vector::const_iterator next_sigma;
  Ddvec::dd_vector::const_iterator first_ladder_sigma;  // first sigma for this eta ladder
  Ddvec::dd_vector::const_iterator last_ladder_sigma;   // last sigma for this eta ladder
  Ddvec::dd_vector::const_iterator sigma_end;           // final sigma

  // pointers to the multiplicity
  Ddvec::dd_vector::const_iterator this_mult;
  Ddvec::dd_vector::const_iterator next_mult;
  Ddvec::dd_vector::const_iterator mult_end;

  // pointers to the model weights
  Ddvec::dd_vector::const_iterator this_weight;
  Ddvec::dd_vector::const_iterator next_weight;
  Ddvec::dd_vector::const_iterator weight_end;

  // pointers to the flux data
  Lgdata::Flux_List::const_iterator flux_ptr;
  Lgdata::Flux_List::const_iterator next_flux;
  Lgdata::Flux_List::const_iterator flux_end;

  // pointers to the incident energy boundaries
  Egp::Energy_groups::const_iterator Ein_ptr;
  Egp::Energy_groups::const_iterator next_Ein;
  Egp::Energy_groups::const_iterator Ein_end;
  int Ein_count;

  //! Holds the weight: (cross section) * flux * multiplicity * model weight
  Lgdata::Legendre_coefs current_weight;

  int order;             // The Legendre order of the problem

  //! Default constructor
  inline param_base()
   {}

  //! Default destructor
  inline ~param_base()
    {}

  //!  Gets the range of nontrivial incident energy bins; computes first_Ein and last_Ein
  //! returns true if the threshold is too high for the energy bins
  //! \param sigma_ the cross section data
  //! \param mult_ the outgoing particle multiplicity data
  //! \param weight_ the weighting to apply to the transfer matrix entries
  //! \param e_flux_ the initial approximation to apply to the particle flux
  //! \param Ein_groups the boundaried of the incident energy groups
  //! \param first_Ein the first common incident energy bin
  //! \param last_Ein the last common incident energy bin
  bool get_Ein_range( const Ddvec::dd_vector& sigma_, const Ddvec::dd_vector& mult_,
    const Ddvec::dd_vector& weight_,
    const Lgdata::Flux_List& e_flux_, const Egp::Energy_groups& Ein_groups,
    double *first_Ein, double *last_Ein );

  //!  Sets up the initial quadrature parameters
  //! \param sigma_ the cross section data
  //! \param mult_ the outgoing particle multiplicity data
  //! \param weight_ the weighting to apply to the transfer matrix entries
  //! \param e_flux_ the initial approximation to apply to the particle flux
  //! \param Ein_groups the boundaries of the incident energy groups
   void setup( const Ddvec::dd_vector& sigma_, const Ddvec::dd_vector& mult_,
    const Ddvec::dd_vector& weight_,
    const Lgdata::Flux_List& e_flux_, const Egp::Energy_groups& Ein_groups );

  //! Sets the first common incident energy
  void common_E0( );

  //!  Sets up the initial quadrature parameters for parallel computing by bin
  //! \param Ein_bin_ the number of the incident energy bin
  //! \param sigma_ the cross section data
  //! \param mult_ the outgoing particle multiplicity data
  //! \param weight_ the weighting to apply to the transfer matrix entries
  //! \param e_flux_ the initial approximation to apply to the particle flux
  //! \param Ein_groups the boundaries of the incident energy groups
  void setup_bin( int Ein_bin_, const Ddvec::dd_vector& sigma_, const Ddvec::dd_vector& mult_,
    const Ddvec::dd_vector& weight_,
    const Lgdata::Flux_List& e_flux_, const Egp::Energy_groups& Ein_groups );

  //! Sets the range of integration over incident energy
  void set_Ein_range( );

  //! Sets the data pointers for a new incident energy interval.
  //! Returns "true" if we are finished with the data.
  //! \param E_in the next incident energy
  bool update_pointers( double E_in );

  //! Sets the data pointers for a new incident energy interval in one bin.
  //! Returns "true" if we are finished with the data on this bin.
  //! \param E_in the next incident energy
  bool update_bin_pointers( double E_in );

  //! Sets the range of pointers to cross sections for this set of data
  void set_sigma_range( );

  //! Calculates the weight: (cross section) * flux * multiplicity * model weight
  //! \param E_in the incident energy
  void set_weight( double E_in );

  //! Calculates the weight for gammas: flux * multiplicity * model weight
  //! \param E_in the incident energy
  void flux_weight( double E_in );

};

}  // end of namespace Pbase

#endif
