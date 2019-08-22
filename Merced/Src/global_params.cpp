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
//Implementation of the GlobalParameterClass

#include <cstdlib>
#ifdef _OPENMP
 #include <omp.h>
#endif

#include "global_params.hpp"
#include "data_parser.hpp"         // converting strings to numbers
#include "messaging.hpp"

// *********** ss_link class ***********************
// ----------- ss_link::print ---------------
// print the data in a link
void ss_link::print()
{
  cout << "   " << x << "   " << y << endl;
}

// *********** ss_list class ***********************
// ----------- ss_list::at ---------------
// Function to get the link at "x"
// Returns "false" if the entry is not found.
bool ss_list::at(string x, ss_list::iterator& link)
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

// ----------- ss_list::print ---------------
void ss_list::print()
{
  for( ss_list::iterator link = begin(); 
       link != end();
       ++link)
  {
    link->print();
  }
  cout<<endl;
}

// *********** GlobalParameterClass class ***********************
// ----------- GlobalParameterClass constructor ---------------
//Default constructor loads the default Parameter list and then reads in the 
//default parameter file and loads them into the Parameter list.
GlobalParameterClass::GlobalParameterClass()   
{
  // Here is the "standard" list of switches
  set("skip_logging", 0);       //! if == 1, don't log errors in documentation.txt file
  set("message_level",0);       //!three levels of messages - info, warn and severe
                                // info = 0 (all messages written)
                                // warn = 1 (only warnings and errors)
                                // severe >=2 (errors only, these cause exits anyway!)
  set("datafield_precision",8); //! 8 sig-figs in all data fields

  set("input",  "");
  set("output",  "utfil");
  set( "outputLegendreOrder", 3 );  // Legendre order of the transfer matrix
  set( "kinetics", "newtonian" );  // "Newtonian" or "relativistic"
  set( "quad_tol", 0.0001 ); // relative error in the quadrature
  set( "abs_quad_tol", 0.0 ); // absolute error in the quadrature
  set( "max_divisions", 1000 );  // limit on subdivisions in adaptive quadrature
  set( "check_row_sum", 0 );  // if > 0, check the sums over each row
  set( "e_tol", 1.0e-10 );  // tolerance on hitting the corner of an energy box
  set( "abs_tol", 2.0e-14 );  // absolute tolerance if E = 0

  // *** Warning *** the set( string, double ) command truncates to 6 digits
  // Use set( string, string )
  set( "m_neutron", "939.565653471" );  // neutron mass in MeV

  set( "max_Eout_ints_size", 5000 );  // a safety check: the maximum number of integrals over E_out
  set("scale_rows", 1 );  // if > 0, scale the rows to match the cross section
  set( "num_threads", 0 );  // number of processors not set
}
// ----------- GlobalParameterClassss destructor --------
//Default destructor.
GlobalParameterClass::~GlobalParameterClass()  
{ 
}
// ----------- GlobalParameterClassss::find ---------------
// Function to get the link at "x"
// Returns "false" if the entry is not found.
bool GlobalParameterClass::find(string x, ss_list::iterator& link)
{
  return Parameters.at( x, link );
}
// ----------- GlobalParameterClass::set ---------------
//Sets the parameter values and creates a new parameter if one doesn't exist.
void GlobalParameterClass::set(string gp, double gp_num)
{
  //First thing we need to do is to strip all extraneous blanks from the string gp
  string_F::remove_all_blanks( gp );
  string_F::Tolower( gp );
  ostringstream gp_ostr; 
  //  Warning( "GlobalParameterClass::set",
  //	   "This routine truncates to 6 significant digits" );
  gp_ostr <<  gp_num;
  string gp_val = gp_ostr.str();
  //We test for the existence of the parameter name in the list.
  //If it is new, we add a link.  If not, we reset the existing
  //value in the link to gp_val.
  ss_list::iterator SSlist;

  if ( Parameters.at(gp,SSlist) )
    {
      SSlist->set_value( gp_val );  //reset the parameter value in this link
    }
  else
    {
      ss_link SSlink;
      SSlink.set_name( gp );  //Load into link
      SSlink.set_value( gp_val );
      Parameters.insert(Parameters.end(), SSlink); //Insert link into list
    }
}
// ----------- GlobalParameterClass::set ---------------
//Sets the parameter values and creates a new parameter if one doesn't exist.
void GlobalParameterClass::set(string gp, string gp_val)
{
  //First thing we need to do is to strip all extraneous blanks from the strings
  string_F::remove_all_blanks( gp );
  string_F::remove_all_blanks( gp_val );
  string_F::Tolower( gp );
  
  //We test for the existence of the parameter name in the list.
  //If it is new, we add a link.  If not, we reset the existing
  //value in the link to gp_val.
  ss_list::iterator SSlist;

  if ( Parameters.at(gp,SSlist) )
    {
      SSlist->set_value( gp_val );  //reset the parameter value in this link
    }
  else
    {
      ss_link SSlink;
      SSlink.set_name( gp );  //Load into link
      SSlink.set_value( gp_val );
      Parameters.insert(Parameters.end(), SSlink); //Insert link into list
    }
}
// ----------- GlobalParameterClass::print ---------------
// Print utility, also saves the parameter values
void GlobalParameterClass::print( )
{
  cout<<"List of global parameters used:"<<endl;
  Parameters.print();
}
// ----------- GlobalParameterClass::read_command_line ---------------
// Function used to parse the command line for parameters and settings
void GlobalParameterClass::read_command_line(int argc, char* argv[] )
{
    // First pass through arg list, look for the input and "-help", 
    // check formatting and convert it all to strings
    for( int j=1; j<argc; ++j ){
        if ( (static_cast<string>(argv[j]) == "-help") || 
             (static_cast<string>(argv[j]) == "--help") ||
             (static_cast<string>(argv[j]) == "-h") ) get_help();
        if ( argv[j][0] == '-' ) {
            if ( (j+1<argc) || (argv[j+1][0] != '-' ) ) {
                string sparam(argv[j]);
                sparam.erase(0,1);
		string_F::Tolower( sparam );
                ++j;
		string svalue(argv[j]);
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
                cerr << "Fatal Error GlobalParameterClass::read_command_line:"<<endl;
                cerr<<"    Key '"<<argv[j]<<"' is missing a value!"<<endl;
                exit(-10);
            }
        } else {
	  set("input",argv[j]);  // the name of the input data file
        }
    }
}
// ----------- GlobalParameterClass::get_help ---------------
// Prints usage message, then exits.
void GlobalParameterClass::get_help( void ){
    cout << "\nUsage: merced [-parameter] [value] input_file" <<endl;
    cout << "\nThis message: merced -help" << endl;
    cout << "\nExample: merced -E_tol 0.0005 input_file\n"<<endl;
    print();
    cout << "The specified input is used as a source for the run.\n";
    cout << "The output is printed to the screen.\n";
    exit(0);
}
// ----------- GlobalParameterClass::Value ---------------
// Returns the value of the requested parameter.
double GlobalParameterClass::Value( string ParamName )
{
  ss_list::iterator XYptr;
  string parm = ParamName;
  string_F::Tolower( ParamName );
  if( !Parameters.at(ParamName,XYptr) )
    {
      cerr<<"Fatal Error GlobalParameterClass::Value:"<<
        " no parameter called "<<parm<<" found!"<<endl;
      exit(-50);
    }
  return( static_cast<double>( atof(XYptr->value().c_str())));
  //  return( static_cast<double>( stof(XYptr->value()) ) );
}
// ----------- GlobalParameterClass::Flag ---------------
// Returns the string value of the parameter.
string GlobalParameterClass::Flag( string ParamName )
{
  ss_list::iterator XYptr;
  string parm = ParamName;
  string_F::Tolower( ParamName );
  if( !Parameters.at(ParamName,XYptr) )
    {
      cerr<<"Fatal Error GlobalParameterClass::Flag:"<<
        " no parameter called "<<parm<<" found!"<<endl;
      exit(-50);
    }
  return( XYptr->value() );
}
// ----------- GlobalParameterClass::set_command ---------------
// for parameter set in the command line
void GlobalParameterClass::set_command(string gp)
{
  ss_list::iterator SSlist;

  if ( Parameters.at(gp,SSlist) )
    {
      SSlist->set_in_command = true;
    }
  else
    {
      FatalError( "GlobalParameterClass::set_command",
		  gp + " not found" );
    }
}
// ----------- GlobalParameterClass::get_field_width ---------------
// Gets the datafield precision for output
int GlobalParameterClass::get_field_width( )
{
  static int field_width = Value( "datafield_precision" ) + 7;
  return field_width;
}
// ----------- GlobalParameterClass::set_num_threads ---------------
// Sets the number of threads to use in parallel computing
void GlobalParameterClass::set_num_threads( )
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

