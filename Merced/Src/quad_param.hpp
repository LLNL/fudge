/*
* ******** merced: calculate the transfer matrix *********
* $Revision: 1 $
* $Date: 2021-10-09 $
* $Author: hedstrom $
* $Id: quad_param.hpp 1 2021-10-09Z hedstrom $
*
* ******** merced: calculate the transfer matrix *********
*
* # <<BEGIN-copyright>>
* # <<END-copyright>>
*/

// header of the base class for quadrature parameters

#ifndef QUAD_PARAM
#define QUAD_PARAM

namespace Qparam
{

  // ----------- class Qparam::QuadParamBase ------------------
  //! Base lass for the quadrature parameters
  class QuadParamBase
  {
  public:
    int func_count;      // the number of function calls

    //! Default constructor
    inline QuadParamBase(): func_count( 0 )
     {}

    //! Default destructor
    inline ~QuadParamBase() {}

  };

}  // end of namespace Qparam

#endif
