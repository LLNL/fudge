/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: data_parser.cpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// coding to read an input file from Python

#include <cstdlib>
#include <cctype>

#include "data_parser.hpp"
#include "messaging.hpp"

// ------------------- string_F::remove_all_blanks -------------------------
// Removes all blanks from a string.
void string_F::remove_all_blanks( string &s )
{
  if( s.size() <= 0 ) return;

  string::size_type nspace = s.find(' ');
  while ( nspace != string::npos )
  {
    s.erase( nspace, 1 );
    nspace = s.find(' ');
  }
}
// ------------------- string_F::remove_final_blanks -------------------------
// Removes final blanks from a string.
void string_F::remove_final_blanks( string &s )
{
  if( s.size() <= 0 ) return;
  {
    string::size_type nspace = s.rfind(' ');
    while ( nspace == s.size()-1 )
    {
      s.erase( s.size()-1,1);
      nspace = s.rfind(' ');
    }
  }
}
// ------------------- string_F::Tolower -------------------------
// A function that lower-cases strings
void string_F::Tolower( string& s)
{
  string::size_type length = s.length();
  for ( string::size_type i = 0; i < length; i++ )
  {
    char this_c = tolower( s[i] );
    s[ i ] = this_c;
  }
}

// ********************* data_parser class ******************
// ------------------- data_parser::data_parser ---------------
// Constructor
data_parser::data_parser( const string &inFileName )
{
  infile.open( inFileName.c_str( ), ios_base::in );
  if( !infile )
  {
    FatalError( "data_parser::data_parser", "No input file: " + inFileName );
  }
  next_line = "";
  pos = string::npos;
  line_count = 0;
}
// ------------------- data_parser::get_new_line ---------------
// Reads in a new line and deletes trailing blanks, returns "false" at end of file
bool data_parser::get_new_line( )
{
  bool got_more = infile.good( );
  while( got_more )
  {
    getline( infile, next_line, '\n' );
    ++line_count;
    if( !infile.good( ) )
    {
      got_more = false;
      break;
    }
    original_line = next_line;
    // Delete '\r' from files uploaded from a Mac.
    string::size_type r_pos = next_line.find_last_of( "\r" );
    if( r_pos != string::npos )
    {
      next_line.erase( r_pos, 1 );
    }
    string::size_type pound_loc = next_line.find( '#' );  // look for comment marker
    if( pound_loc <= next_line.size( ) )
    {
      next_line = next_line.erase( pound_loc, next_line.size( ) );
    }
    string_F::remove_final_blanks( next_line );
    if( next_line.size() > 0 )
    {
      break;
    }
  }
  return got_more;
}

// ------------------- data_parser::get_dataID ---------------
// Gets the data identifier
string data_parser::get_dataID( )
{
  string dataID;
  if( pos == string::npos )
  {
    if( !get_new_line( ) )
    {
      //      Info( "data_parser::get_dataID", "End of input" );
      dataID = "DONE";
    }
    else
    {
      pos = next_line.find_first_of( ':' );
      if( pos == string::npos )
      {
	dataID = "Bad dataID";
      }
      else
      {
        dataID = next_line.substr( 0, pos );
        ++pos;
      }
    }
  }
  return dataID;
}

// ------------------- data_parser::get_next_int ---------------
// Gets the next integer in next_line, starting at pos
int data_parser::get_next_int(  )
{
  if( pos == string::npos )
  {
    if( !get_new_line( ) )
    {
      FatalError( "data_parser::get_next_int", "End of input" );
    }
    pos = 0;
  }
  string numerics( "+-0123456789" );
  pos = next_line.find_first_of( numerics, pos );
  string::size_type next_pos = next_line.find_first_of( " ", pos );
  string number = next_line.substr( pos, next_pos - pos );
  pos = ( next_pos == string::npos ) ? next_pos : next_pos + 1;
  return atoi( number.c_str( ) );
}

// ------------------- data_parser::get_next_double ---------------
// Gets the next double in next_line, starting at pos
double data_parser::get_next_double( )
{
  if( pos == string::npos )
  {
    if( !get_new_line( ) )
    {
      FatalError( "data_parser::get_next_double", "End of input" );
    }
    pos = 0;
  }
  string numerics( "+-0123456789" );
  pos = next_line.find_first_of( numerics, pos );
  string::size_type next_pos = next_line.find_first_of( " ", pos );
  string number = next_line.substr( pos, next_pos - pos );
  pos = ( next_pos == string::npos ) ? next_pos : next_pos + 1;
  return atof( number.c_str( ) );
}

// ------------------- data_parser::get_text ---------------
//! Extracts text to the end of this line
string data_parser::get_text( )
{
  // remove possible single quotes
  string::size_type prev_pos = next_line.find_first_of( "'", pos );
  string::size_type next_pos;
  if( prev_pos == string::npos )
  {
    next_pos = prev_pos;
    prev_pos = pos;
  }
  else
  {
    ++prev_pos;
    next_pos = next_line.find_last_of( "'" );
    if( next_pos < prev_pos )
    {
      Warning( "data_parser::get_text", "unmatched single quotes" );
      next_pos = string::npos;
    }
  }
  // remove initial spaces
  string::size_type length = next_line.length( );
  for ( string::size_type i = prev_pos; (i < length) && isspace( next_line[i] ); i++ )
  {
    ++prev_pos;
  }
  pos = string::npos;
  return next_line.substr( prev_pos, next_pos - prev_pos );
}

// ------------------- data_parser::get_comment ---------------
// Reads a comment
void data_parser::get_comment( ofstream &output_file )
{
  string Comment = next_line.substr( pos, next_line.size( ) );
  output_file << "Comment: " << Comment << endl;
  Info( "Comment", Comment );
  pos = string::npos;
}
// ------------------- data_parser::read_2d_interp ---------------
// Reads the interpolation rule for 2-dimensional tables
void data_parser::read_2d_interp( string *interpolation, string *qualifier )
{
  string rule = get_text( );
  string_F::Tolower( rule );
  string::size_type space_pos = rule.find_first_of( " " );
  if( space_pos == string::npos )
  {
    // histogram in incident energy is special
    if( rule == "flat" )
    {
      *interpolation = rule;
      *qualifier = "direct";
    }
    else
    {
      FatalError( " data_parser::read_2d_interp", "no qualifier found" );
    }
  }
  else
  {
    *interpolation = rule.substr( 0, space_pos );
    *qualifier = rule.substr( space_pos+1, rule.length( ) );
    string_F::remove_all_blanks( *qualifier );
  }
}
