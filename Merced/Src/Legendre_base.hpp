/*
* ******** merced: calculate the transfer matrix *********
* $Revision: 1 $
* $Date: 2021-10-09 $
* $Author: hedstrom $
* $Id: Legendre_base.hpp 1 2021-10-09Z hedstrom $
*
* ******** merced: calculate the transfer matrix *********
*
* # <<BEGIN-copyright>>
* # <<END-copyright>>
*/

// header for the base class for Legendre data

#ifndef LEGENDRE_BASE
#define LEGENDRE_BASE

namespace LgBase
{

  // ------------- class LgBase::Legendre_base ---------------------
  //! Class to hold Legendre coefficients j=0, 1, ..., order
  class Legendre_base
  {
  private:
    double Energy;

  public:
    //! Legendre coefficients
    double *data;

    //! the Legendre order
    int order;

    //! Constructor
    inline Legendre_base( ): Energy(-1.0), order(-1) {}

    //! Destructor
    inline ~Legendre_base( ) { clean_data( ); }

    //! "Energy" could be for the incident or the outgoing particle
    inline double get_E_out( ) const { return Energy; }
    inline double get_E_in( ) const { return Energy; }

    //! Sets the energy for the Legendre data
    //! \param Eout the energy for the Legendre data
    inline void set_E_out( double Eout ) { Energy = Eout; }
    
    //! Sets the energy for the Legendre data
    //! \param Ein the energy for the Legendre data
    inline void set_E_in( double Ein ) { Energy = Ein; }

    void clean_data( );

    //! Allocates space
    //! \param order the Legendre order for the data
    void initialize( int order );

    //! Access routine for N-th coefficient
    //! \param N number of the current Legendre coefficient
    double &operator[ ]( int N );

    //! Access routine for N-th coefficient; does not change the data
    //! \param N number of the current Legendre coefficient
    double value( int N ) const;

    //! Ignore zero high-order Legendre coefficients
    void truncate_zeros( );

    //! Scales the vector
    //! \param factor multiply all coefficients by this number
    Legendre_base& operator*=( double factor );

    //! Sums the Legenre series
    //! \param mu sum the series at this mu value
    double sum_Legendre( double mu );

    // For debugging
    void print( ) const;
  };


}  // end of namespace LgBase

#endif
