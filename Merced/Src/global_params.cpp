/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: global_params.cpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
//Implementation of the Glb::GlobalParameterClass

#include <cstdlib>
#ifdef _OPENMP
 #include <omp.h>
#endif

#include "global_params.hpp"
#include "data_parser.hpp"         // converting strings to numbers
#include "messaging.hpp"

// *********** Glb::ss_link class ***********************
// ----------- Glb::ss_link::print ---------------
// print the data in a link
void Glb::ss_link::print()
{
  std::cout << "   " << x << "   " << y << std::endl;
}

// *********** Glb::ss_list class ***********************
// ----------- Glb::ss_list::at ---------------
// Function to get the link at "x"
// Returns "false" if the entry is not found.
bool Glb::ss_list::at(std::string x, Glb::ss_list::iterator& link)
{
  for( link = begin(); 
       link != end();
       ++link )
  {
    if ( link->name( ) == x ) 
    {
      return true;
    }
  }
  return false;
}

// ----------- Glb::ss_list::print ---------------
void Glb::ss_list::print()
{
  for( Glb::ss_list::iterator link = begin(); 
       link != end();
       ++link)
  {
    link->print();
  }
  std::cout<<std::endl;
}

// *********** Glb::GlobalParameterClass class ***********************
// ----------- Glb::GlobalParameterClass constructor ---------------
//Default constructor loads the default Parameter list and then reads in the 
//default parameter file and loads them into the Parameter list.
Glb::GlobalParameterClass::GlobalParameterClass()   
{
  // Here is the "standard" list of switches
  set("skip_logging", 0);       //! if == 1, don't log errors in documentation.txt file
  set("message_level", 2);   //!four levels of messages - info, warn, debug, none
                             // 0: only fatal errors
                             // 3: only informational messages and fatal errors
                             // 2: warnings, information, and fatal errors)
                             // 1: debugging information plus all the rest
  set("datafield_precision",8); //! 8 sig-figs in all data fields

  set("input",  "");
  set("output",  "utfil");
  set( "outputLegendreOrder", 3 );  // Legendre order of the transfer matrix
  set( "kinetics", "newtonian" );  // "Newtonian" or "relativistic"
  set( "quad_tol", 0.0001 ); // relative error in the quadrature
  set( "short_interval", 1.0e-12 );  // accept the current estimate for intervals shorter than this relative length
  set( "quad_weight_increase", 2.5 );  // rate of increase of quadrature tolerance for high Legendre orders
  set( "quad_tol_floor", 1.0e-5 );  // floor for the quadrature tolerance for high Legendre orders
  set( "max_divisions", 1000 );  // limit on subdivisions in adaptive quadrature
  set( "flag_adapt_quad_depth", 10 );  // flag when adaptive quadrature hits this depth
  set( "check_row_sum", 0 );  // if > 0, check the sums over each row
  set( "tight_tol", 2.0e-14 );  // a tight tolerance on energy differences
  set( "looser_tol", 1.0e-10 );  // looser tolerance on energy differences
  set( "cum_delta_tol", 1.0e-12 ); // when to use delta-functions for cumulative points
  set( "norm_warn", 0 );  // set > 0 to warn about norm != 1 in the data 
  set( "norm_tol", 2.0e-5 );  // tolerance for warning that norm != 1 in the data
  set( "cum_prob_skip", 1.0e-10 );  // skip intervals of small probability
  set( "sqrt_wt_cutoff", 0.1 ); // use weight sqrt(x) on the interval (0, cutoff)
 
  // *** Warning *** the set( string, double ) command truncates to 6 digits
  // Use set( std::string, std::string )
  set( "m_neutron", "939.565653471" );  // neutron mass in MeV

  set( "max_Eout_ints_size", 5000 );  // a safety check: the maximum number of integrals over E_out---obsolete
  set( "scale_rows", 1 );  // if > 0, scale the rows to match the cross section
  set( "truncate_direct", 1 );  // 0: extrapolate direct interpolation
                                // 1: truncate direct interpolation
  set( "num_threads", 0 );  // number of processors not set
}
// ----------- Glb::GlobalParameterClassss destructor --------
//Default destructor.
Glb::GlobalParameterClass::~GlobalParameterClass()  
{ 
}
// ----------- Glb::GlobalParameterClassss::find ---------------
// Function to get the link at "x"
// Returns "false" if the entry is not found.
bool Glb::GlobalParameterClass::find(std::string x, Glb::ss_list::iterator& link)
{
  return Parameters.at( x, link );
}
// ----------- Glb::GlobalParameterClass::set ---------------
//Sets the parameter values and creates a new parameter if one doesn't exist.
void Glb::GlobalParameterClass::set(std::string gp, double gp_num)
{
  //First thing we need to do is to strip all extraneous blanks from the string gp
  string_F::remove_all_blanks( gp );
  string_F::Tolower( gp );
  std::ostringstream gp_ostr; 
  //  Warning( "Glb::GlobalParameterClass::set",
  //	   "This routine truncates to 6 significant digits" );
  gp_ostr <<  gp_num;
  std::string gp_val = gp_ostr.str();
  //We test for the existence of the parameter name in the list.
  //If it is new, we add a link.  If not, we reset the existing
  //value in the link to gp_val.
  Glb::ss_list::iterator SSlist;

  if ( Parameters.at(gp,SSlist) )
    {
      SSlist->set_value( gp_val );  //reset the parameter value in this link
    }
  else
    {
      Glb::ss_link SSlink;
      SSlink.set_name( gp );  //Load into link
      SSlink.set_value( gp_val );
      Parameters.insert(Parameters.end(), SSlink); //Insert link into list
    }
}
// ----------- Glb::GlobalParameterClass::set ---------------
//Sets the parameter values and creates a new parameter if one doesn't exist.
void Glb::GlobalParameterClass::set(std::string gp, std::string gp_val)
{
  //First thing we need to do is to strip all extraneous blanks from the strings
  string_F::remove_all_blanks( gp );
  string_F::remove_all_blanks( gp_val );
  string_F::Tolower( gp );
  
  //We test for the existence of the parameter name in the list.
  //If it is new, we add a link.  If not, we reset the existing
  //value in the link to gp_val.
  Glb::ss_list::iterator SSlist;

  if ( Parameters.at(gp,SSlist) )
    {
      SSlist->set_value( gp_val );  //reset the parameter value in this link
    }
  else
    {
      Glb::ss_link SSlink;
      SSlink.set_name( gp );  //Load into link
      SSlink.set_value( gp_val );
      Parameters.insert(Parameters.end(), SSlink); //Insert link into list
    }
}
// ----------- Glb::GlobalParameterClass::print ---------------
// Print utility, also saves the parameter values
void Glb::GlobalParameterClass::print( )
{
  std::cout<<"List of global parameters used:"<<std::endl;
  Parameters.print();
}
// ----------- Glb::GlobalParameterClass::read_command_line ---------------
// Function used to parse the command line for parameters and settings
void Glb::GlobalParameterClass::read_command_line(int argc, char* argv[] )
{
    // First pass through arg list, look for the input and "-help", 
    // check formatting and convert it all to strings
    for( int j=1; j<argc; ++j ){
        if ( (static_cast<std::string>(argv[j]) == "-help") || 
             (static_cast<std::string>(argv[j]) == "--help") ||
             (static_cast<std::string>(argv[j]) == "-h") ) get_help();
        if ( argv[j][0] == '-' ) {
            if ( (j+1<argc) || (argv[j+1][0] != '-' ) ) {
                std::string sparam(argv[j]);
                sparam.erase(0,1);
		string_F::Tolower( sparam );
                ++j;
		std::string svalue(argv[j]);
                if( svalue[0] == 't' )
                {
                  svalue = "1";
                }
                else if( svalue[0] == 'f' )
                {
                  svalue = "0";
                }
                set( sparam, svalue );
		set_command( sparam );
            } else {
                std::cerr << "Fatal Error Glb::GlobalParameterClass::read_command_line:"<<std::endl;
                std::cerr<<"    Key '"<<argv[j]<<"' is missing a value!"<<std::endl;
                exit(-10);
            }
        } else {
	  set("input",argv[j]);  // the name of the input data file
        }
    }
}
// ----------- Glb::GlobalParameterClass::get_help ---------------
// Prints usage message, then exits.
void Glb::GlobalParameterClass::get_help( void ){
    std::cout << "\nUsage: merced [-parameter] [value] input_file" <<std::endl;
    std::cout << "\nThis message: merced -help" << std::endl;
    std::cout << "\nExample: merced -E_tol 0.0005 input_file\n"<<std::endl;
    print();
    std::cout << "The specified input is used as a source for the run.\n";
    std::cout << "The output is printed to the screen.\n";
    exit(0);
}
// ----------- Glb::GlobalParameterClass::Value ---------------
// Returns the value of the requested parameter.
double Glb::GlobalParameterClass::Value( std::string ParamName )
{
  Glb::ss_list::iterator XYptr;
  std::string parm = ParamName;
  string_F::Tolower( ParamName );
  if( !Parameters.at(ParamName,XYptr) )
    {
      std::cerr<<"Fatal Error Glb::GlobalParameterClass::Value:"<<
        " no parameter called "<<parm<<" found!"<<std::endl;
      exit(-50);
    }
  return( static_cast<double>( atof(XYptr->value().c_str())));
  //  return( static_cast<double>( stof(XYptr->value()) ) );
}
// ----------- Glb::GlobalParameterClass::Flag ---------------
// Returns the string value of the parameter.
std::string Glb::GlobalParameterClass::Flag( std::string ParamName )
{
  Glb::ss_list::iterator XYptr;
  std::string parm = ParamName;
  string_F::Tolower( ParamName );
  if( !Parameters.at(ParamName,XYptr) )
    {
      std::cerr<<"Fatal Error Glb::GlobalParameterClass::Flag:"<<
        " no parameter called "<<parm<<" found!"<<std::endl;
      exit(-50);
    }
  return( XYptr->value() );
}
// ----------- Glb::GlobalParameterClass::set_command ---------------
// for parameter set in the command line
void Glb::GlobalParameterClass::set_command(std::string gp)
{
  Glb::ss_list::iterator SSlist;

  if ( Parameters.at(gp,SSlist) )
    {
      SSlist->set_in_command = true;
    }
  else
    {
      Msg::FatalError( "Glb::GlobalParameterClass::set_command",
		  gp + " not found" );
    }
}
// ----------- Glb::GlobalParameterClass::get_field_width ---------------
// Gets the datafield precision for output
int Glb::GlobalParameterClass::get_field_width( )
{
  static int field_width = Value( "datafield_precision" ) + 7;
  return field_width;
}
// ----------- Glb::GlobalParameterClass::set_num_threads ---------------
// Sets the number of threads to use in parallel computing
void Glb::GlobalParameterClass::set_num_threads( )
{
  int num_threads;
#ifdef _OPENMP
  num_threads = Value( "num_threads" );
  if( num_threads == 0 )
  { 
    num_threads = omp_get_num_procs( );
  }
  else
  {
    omp_set_num_threads( num_threads );
  }
#else
  num_threads = 1;
#endif
  set( "num_threads", num_threads );
}

