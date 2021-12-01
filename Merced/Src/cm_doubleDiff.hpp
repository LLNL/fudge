/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 601 $
 * $Date: 2017-12-15 $
 * $Author: hedstrom $
 * $Id: cm_doubleDiff.hpp 601 2017-12-15Z hedstrom $
 *
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// declaration of the classes used to handle pointwise energy-angle
// probability density given in the center-of-mass frame

#ifndef CM_DOUBLE_DIFF
#define CM_DOUBLE_DIFF

#include "doubleDiff_base.hpp"
#include "Ecm_Elab_geom.hpp"
#include "math_util.hpp"
#include "transfer.hpp"

namespace cmDD
{
  // forward reference
  class joint_dist_param;


  // ----------- class joint_Eout_param -----------------
  //! parameters for 2-d integration over Eout and mu
  class joint_Eout_param: public Egeom::Ecm_Elab_Ecm_param
  {
  private:

  public:
 
    Terp::two_d_interp Eout_interp; // interpolation rule for outgoing energy

    double E_in0;  // lower incident energy of the data
    double E_in1;  // higher incident energy of the data

    //! The current incident energy
    double current_E_in;

    //! unit-base map of the outgoing energy
    Ddvec::unit_base_map mid_ubase_map;  // outgoing energy range at the current E_in

    //! where we are in the data, which Eout values
    DDbase::one_joint_dist::const_iterator Ein0_Eout0;  // lower Ein, lower Eout
    DDbase::one_joint_dist::const_iterator Ein0_Eout1;  // lower Ein, higher Eout
    DDbase::one_joint_dist::const_iterator Ein1_Eout0;  // higher Ein, lower Eout
    DDbase::one_joint_dist::const_iterator Ein1_Eout1;  // higher Ein, higher Eout

    Qmeth::Quadrature_Rule mu_quad_rule;  // quadrature rule for outgoing energy
    long int mu_F_count;  // number of calls to joint_dist_F::mu_F

    //! Default constructor
    inline joint_Eout_param( ): mu_F_count ( 0 )
    {}

    //! Default destructor
    inline ~joint_Eout_param( )
    {}

  };

  // ----------- class joint_mu_param -----------------
  //! parameters for 2-d integration over mu and E_out
  class joint_mu_param: public Egeom::Ecm_Elab_mu_param
  {
  private:

    double E_in0;  // lower incident energy of the data
    double E_in1;  // higher incident energy of the data

    //! where we are in the (Ein, Eout, mu, probability) data
    DDbase::one_Eout::const_iterator Ein0_Eout0_mu0;  // low Ein, low Eout, low mu
    DDbase::one_Eout::const_iterator Ein0_Eout1_mu0;  // low Ein, high Eout, low mu
    DDbase::one_Eout::const_iterator Ein0_Eout0_mu1;  // low Ein, low Eout, high mu
    DDbase::one_Eout::const_iterator Ein0_Eout1_mu1;  // low Ein, high Eout, high mu
    DDbase::one_Eout::const_iterator Ein1_Eout0_mu0;  // high Ein, low Eout, low mu
    DDbase::one_Eout::const_iterator Ein1_Eout1_mu0;  // high Ein, high Eout, low mu
    DDbase::one_Eout::const_iterator Ein1_Eout0_mu1;  // high Ein, low Eout, high mu
    DDbase::one_Eout::const_iterator Ein1_Eout1_mu1;  // high Ein, high Eout, high mu

    Qmeth::Quadrature_Rule mu_quad_rule;

    //! Interpolates to the initial mu data
    void start_mu_data( );

    //! Interpolates to the next set of mu data
    //! Returns true when the are no more data
    bool next_mu_data( );

    //! Finds the pointers to the data for the current mu value
    //! \param Mu, the desired mu value
    void find_mu( double Mu );

    //! Finds the interpolated probability density at this mu
    //! \param Mu, the desired mu value
    double get_prob( double Mu );

    //! Sets the interpolated probability density at the higher mu
    void set_high_mu( );

  public:
    //! where we are in the data, which Eout values
    DDbase::one_joint_dist::const_iterator Ein0_Eout0;  // lower Ein, lower Eout
    DDbase::one_joint_dist::const_iterator Ein0_Eout1;  // lower Ein, higher Eout
    DDbase::one_joint_dist::const_iterator Ein1_Eout0;  // higher Ein, lower Eout
    DDbase::one_joint_dist::const_iterator Ein1_Eout1;  // higher Ein, higher Eout

    //! data interpolated to the current incident energy and unit-base outgoing energy
    Ddvec::dd_pair current_data;

    //! interpolated unit-base map of the outgoing energy
    Ddvec::unit_base_map *mid_ubase_map;

    long int mu_F_count;  // number of calls to cmDD::joint_dist_F::mu_F

    //! the outgoing energy bin boundaries
    std::vector< double >::const_iterator Eout_bottom;
    std::vector< double >::const_iterator Eout_top;

    //! Default constructor
    inline joint_mu_param( ): mu_F_count ( 0 )
    {}

    //! Default destructor
    inline ~joint_mu_param( )
    {}

    //! Sets up the data for this incident energy and outgoing energy
    //! \param Ein energy of incident particle
    //! \param Eoutcm center-of-mass energy of outgoing particle
    //! \param Ecm_param parameters for integration over outgoing energy
    void setup( double Ein, double Eoutcm, cmDD::joint_Eout_param& Ecm_param );

    //! Performs the integral over cm mu
    //! Returns true if the interpolation is OK
    //! \param value, the computed integral
    bool scan_data( Coef::coef_vector *value );

  };

  // ----------- class cmDD::joint_dist -----------------
  //! Class for joint energy-angle distributions
  class joint_dist : public DDbase::joint_dist_base
  {
  private:
    //! The smallest incident energy for cross section, multiplicity,
    //! model weight, flux weight, and energy groups
    double E_first;

    int first_Ein;  // index of the left-hand end of the first significant energy bin
    int last_Ein;  // index of the right-hand end of the last significant energy bin

    //! Sets up the map from center-of-mass to laboratory coordinates
    void setup_map( );

    //!  Gets the range of nontrivial incident energy bins; computes first_Ein and last_Ein
    //! returns true if the threshold is too high for the energy bins
    //! \param sigma the cross section data
    //! \param mult the multiplicity of the outgoing particle
    //! \param weight the weighting to apply to the transfer matrix entries
    //! \param e_flux approximate flux used to weight the transfer matrix
    //! \param Ein_groups the boundaries of the incident energy groups
    bool get_Ein_range( const Ddvec::dd_vector& sigma, const Ddvec::dd_vector& mult,
      const Ddvec::dd_vector& weight,
      const Lgdata::Flux_List& e_flux, const Egp::Energy_groups& Ein_groups );

    //! Initializes the quadrature parameters
    //! \param Ein_param the qudrature parameters
    void setup_param( cmDD::joint_dist_param *Ein_param );

    //! Sets the range of incident energies for this intergration
    //! \param Ein_param the qudrature parameters
    void set_Ein_range( cmDD::joint_dist_param *Ein_param );

    //! Go to the next pair of incident energies.  Returns "true" when finished.
    //! \param E_in energy of the incident particle
    //! \param Ein_param the qudrature parameters
    bool next_Ein_pair( double E_in, cmDD::joint_dist_param *Ein_param );

    //! Initializes the pointers to the ( Eout, probability ) data for current cosines
    //! \param Ein_param the qudrature parameters
    void start_mu_data( cmDD::joint_dist_param *Ein_param );

    //! Go to the next pairs of direction cosines.  Returns "true" when finished.
    //! \param Ein_param the qudrature parameters
    bool next_mu_pairs( cmDD::joint_dist_param *Ein_param );

    //! Adds to the transfer matrix for all E_out bins for given Ein and Eout ranges
    //! \param transfer the transfer matrix
    //! \param Ein_param the qudrature parameters
    void Eout_data_ladder( Trf::T_matrix& transfer,
			   cmDD::joint_dist_param *Ein_param );

    //! go to the next sets of (E_out, probability) pairs
    //! \param Ein_param the qudrature parameters
    bool next_Eout( cmDD::joint_dist_param *Ein_param );

    //! Loops through the outgoing energy bins for given cm_Eout data
    //! \param transfer: the entries in the transfer matrix get updated
    //! \param Ein_param the quadrature parameters
    void lab_Eout_ladder( Trf::T_matrix& transfer,
			  cmDD::joint_dist_param *Ein_param );

    //! Adds to the transfer matrix for the current data
    //! \param transfer the transfer matrix
    //! \param Eout_count the current row of the transfer matrix
    //! \param Ein_param the qudrature parameters
    void one_Ebox( Trf::T_matrix& transfer, int Eout_count,
		   cmDD::joint_dist_param *Ein_param );

    //! Increments the transfer matrix
    //! \param transfer the transfer matrix
    //! \param Eout_count the current row of the transfer matrix
    //! \param Ein_param the qudrature parameters
    void update_T( Trf::T_matrix& transfer, int Eout_count,
                   cmDD::joint_dist_param *Ein_param );

    //! Go to the next pair of incident energies.  Returns "true" when finished.
    //! \param  E_in the next incident energy
    //! \param Ein_param the quadrature parameters
    bool next_ladder( double E_in, cmDD::joint_dist_param *Ein_param );

  public:
    Maps::particleInfo particles;
    Maps::map_cm_lab map;

    //! Default constructor
    joint_dist( )
    {}

    //! Default destructor
    inline ~joint_dist( )
    {}

    //! Calculates the transfer matrix for this particle.
    //! \param sigma the cross section data
    //! \param mult the outgoing particle multiplicity data
    //! \param weight the weighting to apply to the transfer matrix entries
    //! \param transfer the transfer matrix
    void get_T( const Ddvec::dd_vector& sigma, const Ddvec::dd_vector& mult,
                const Ddvec::dd_vector& weight, Trf::T_matrix& transfer );
    
  };

  // ----------- class joint_dist_param -----------------
  //! parameters for 3-d integration
  class joint_dist_param: public Egeom::Ecm_Elab_Ein_param
  {
  private:

    //! Sets up the data for unit-base interpolation in incident energy
    void setup_Ein_ubase( );

  public:
    // where the lower and upper eta values hit the quadrature box
    Vhit::Vcm_Vlab_hit_list lower_hits;
    Vhit::Vcm_Vlab_hit_list upper_hits;

    // parameters for integration over outgoing cosine and energy
    cmDD::joint_Eout_param Ecm_params;

    //! where we are in the list, which incident energies
    cmDD::joint_dist::const_iterator Ein0_data;
    cmDD::joint_dist::const_iterator Ein1_data;

    // where we are in the data, which Eout values
    //DDbase::one_joint_dist::const_iterator Ein0_Eout0;  // lower Ein, lower mu
    //DDbase::one_joint_dist::const_iterator Ein0_Eout1;  // lower Ein, higher mu
    //DDbase::one_joint_dist::const_iterator Ein1_Eout0;  // higher Ein, lower mu
    //DDbase::one_joint_dist::const_iterator Ein1_Eout1;  // higher Ein, higher mu

    // current unit-base outgoing energy range
    double prev_data_Eout;
    double next_data_Eout;

    Terp::two_d_interp Ein_interp; // interpolation rule for incident energy
    Terp::two_d_interp Eout_interp; // interpolation rule for outgoing energy

    Qmeth::Quadrature_Rule Eout_quad_rule;  // quadrature rule for outgoing energy
    Qmeth::Quadrature_Rule mu_quad_rule;  // quadrature rule for outgoing cosine

    long int quad_count;  // number of 3-d quadratures
    long int Ein_F_count;  // number of calls to cmDD::joint_dist_F::Ein_F
    long int Eout_F_count;  // number of calls to cmDD::joint_dist_F::Eout_F
    long int mu_F_count;  // number of calls to cmDD::joint_dist_F::mu_F
    int Vcm_hit_count;  // number of calls to quad_F::integrate in cmDD::joint_dist_F::Ein_F

    //! constructor
    inline joint_dist_param( ): quad_count( 0 ), Ein_F_count ( 0 ),
              Eout_F_count ( 0 ), mu_F_count ( 0 ), Vcm_hit_count( 0 ) {}

    //! Default destructor
    inline ~joint_dist_param( )
    {}

    //! Starts one staircase of the Eout_cm data
    void start_Eout_cm( );

    //! Gets the next interval of unit-base outgoing energy
    //! \param Ein_param the quadrature parameters
    bool next_Ecm( );

    //! Sets up the parameters for integration over center-of-mass outgoing energy
    //! \param E_in, the current incident energy
    void set_Ecm_param( double E_in );

  };

}  // end of namespace cmDD

// **************** Functions to integrate *********************
namespace cmDD_F
{
  //! Function for the 1-d quadrature over cm cosine
  //! Returns true if the interpolation is OK.
  //! \param mu the center-of-mass direction cosine of the outgoing particle
  //! \param mu_quad_param the function parameters
  //! \param value the value of the integrand, a set of Legendre coefficients
  bool mu_F( double mu, Qparam::QuadParamBase *mu_quad_param,
	     Coef::coef_vector *value );

  //! Function for the 2-d quadrature over cm cosine and Eout_cm
  //! Returns true if the interpolation is OK.
  //! \param Eout_cm the center-of-mass energy of the outgoing particle
  //! \param Ecm_quad_param the function parameters
  //! \param value the value of the integrand, a set of Legendre coefficients
  bool Ecm_F( double Eout_cm, Qparam::QuadParamBase *Ecm_quad_param,
	      Coef::coef_vector *value );

  //! Function for the 3-d quadrature over E_in, and Eout_cm and cm cosine
  //! Returns true if the interpolation is OK.
  //! The value of Ein_F is itself an integral over Eout_cm and cm cosine.
  //! \param E_in the energy of the incident particle
  //! \param Ein_quad_param the function parameters
  //! \param value the value of the integrand, a set of Legendre coefficients
  bool Ein_F( double E_in, Qparam::QuadParamBase *Ein_quad_param,
	      Coef::coef_vector *value );
}

#endif

