/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: data_parser.hpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// declarations for an input file from Python

#ifndef data_parser_DEF
#define data_parser_DEF

#include <cstdio>
#include <fstream>
#include <cstring>

// ------------------- string routines ----------------
namespace string_F
{
  //! Removes extra blanks from a string
  //! \param s the string to modify
  void remove_all_blanks( std::string &s );

  //! Removes final blanks from a string.
  //! \param s the string to modify
  void remove_final_blanks( std::string &s );

  //! A function that lower-cases strings
  //! \param s the string to modify
  void Tolower( std::string& s);
}

namespace Dpar
{
// ------------------- class data_parser ---------------
//! Base class for reading Python file sections
class data_parser
{
private:
  std::ifstream infile;

  //! the current position in next_line
  std::string::size_type pos;

  //! Reads in a new line and deletes trailing blanks, returns "false" at end of file
  bool get_new_line( );

public:
  //! the next line in the file
  std::string next_line;

  //! original line, in case of trouble
  std::string original_line;

  int line_count;

  data_parser( const std::string &inFileName );

  inline ~data_parser( ) { }

  //! Gets the data identifier
  std::string get_dataID( );

  //! Gets the next integer in next_line, starting at pos
  int get_next_int(  );

  //! Gets the next double in next_line, starting at pos
  double get_next_double( );

  //! Extracts and prints a comment
  //! \param output_file the output file
  void get_comment( std::ofstream &output_file );

  //! Extracts text to the end of this line
  std::string get_text( );

  //! Reads the interpolation rule for 2-dimensional tables
  //! \param interpolation: lin-lin, etc.
  //! \param qualifier: direct, unitbase, etc.
  void read_2d_interp( std::string *interpolation, std::string *qualifier );
};

} // end of namespace Dpar

#endif
