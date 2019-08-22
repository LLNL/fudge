/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: main.hpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */
#include <iostream>
#include <string>

#include "messaging.hpp"
#include "global_params.hpp"  // for input parameters
// #include "version.hpp"
// #include "logger.hpp"
#include "reaction.hpp"
#include "data_parser.hpp"

using namespace std;

// Handle the input parameters
GlobalParameterClass Global;

// At the end, we write out the evaluation documentation, version info and
// any processing errors.  This holds that data.
// MessageLogger messageLog;

int main(int argc, char* argv[])
{
  // Read input parameters by user
  Global.read_command_line(argc,argv);

  // Set the number of threads for parallel processing
  Global.set_num_threads( );

  ofstream output_file;

  string input_file_name = Global.Flag( "input" );
  data_parser input_file( input_file_name );
  string output_file_name = Global.Flag( "output" );
  output_file.open( output_file_name.c_str( ) );
  if( !output_file )
  {
    FatalError( "main", "Unable to open output file: " + output_file_name );
  }

  reaction react_data;
  react_data.process_data( input_file, &output_file );

  //  messageLog.write();

  return 0;
}
