/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: dd_vector.cpp 1 2006-02-02 03:06:56Z hedstrom$
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// implementation of the dd_vector class

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

#include "dd_vector.hpp"
#include "messaging.hpp"
#include "global_params.hpp"

// *********** class Terp::two_d_interp ***********************
// ---------------- Terp::two_d_interp::operator= -----------------
// to copy
Terp::two_d_interp& Terp::two_d_interp::operator=( const Terp::two_d_interp &to_copy )
{
  flag = to_copy.flag;
  qualifier = to_copy.qualifier;

  return *this;
}

// *********** interp_flag_F::read_flag ***********************
// Interprets the interpolation rule
Terp::Interp_Type interp_flag_F::read_flag( const std::string& interp_ID,
  Dpar::data_parser &input_file )
{
  Terp::Interp_Type interp = Terp::NOTSET;

  if( interp_ID == "lin-lin" )
  {
    interp = Terp::LINLIN;
  }
  else if( interp_ID == "flat" )
  {
    interp = Terp::HISTOGRAM;
  }
  else if( interp_ID == "lin-log" )
  {
    interp = Terp::LINLOG;
  }
  else if( interp_ID == "log-lin" )
  {
    interp = Terp::LOGLIN;
  }
  else if( interp_ID == "log-log" )
  {
    interp = Terp::LOGLOG;
  }
  else
  {
    Msg::FatalError( "interp_flag_F::read_flag", 
		Msg::pastenum( "line ", input_file.line_count ) +
              " interpolation flag " +
		input_file.original_line + " undefined." );
  }
  return interp;
}

// *********** interp_flag_F::read_qualifier ***********************
// Interprets the interpolation qualifier
Terp::Interp_qualifier interp_flag_F::read_qualifier( const std::string& qualifier_ID,
  Dpar::data_parser &input_file )
{
  Terp::Interp_qualifier qualifier = Terp::UNSET;

  if( qualifier_ID == "direct" )
  {
    qualifier = Terp::DIRECT;
  }
  else if( qualifier_ID == "unitbase" )
  {
    qualifier = Terp::UNITBASE;
  }
  else if( qualifier_ID == "cumulativepoints" )
  {
    qualifier = Terp::CUMULATIVE_POINTS;
  }
  else if( qualifier_ID == "unitbase-unscaled" )
  {
    qualifier = Terp::UNSCALED_UNITBASE;
  }
  else if( qualifier_ID == "cumulativepoints-unscaled" )
  {
    qualifier = Terp::UNSCALED_CUMULATIVE_POINTS;
  }
  else if( qualifier_ID == "direct-unscaled" )
  {
    qualifier = Terp::UNSCALED_DIRECT;
  }
  else
  {
    Msg::FatalError( "interp_flag_F::read_qualifier", 
		Msg::pastenum( "line ", input_file.line_count ) +
              " interpolation qualifier " +
		input_file.original_line + " undefined." );
  }
  return qualifier;
}

// *********** interp_flag_F::read_1d_interpolation ***********************
// Interprets the interpolation rule
Terp::Interp_Type interp_flag_F::read_1d_interpolation( Dpar::data_parser &input_file )
{
  std::string dataID = input_file.get_dataID( );
  string_F::Tolower( dataID );
  if( dataID != "interpolation" )
  {
    Msg::FatalError( "interp_flag_F::read_1d_interpolation",
                "Expected interpolation rule in " +
		Msg::pastenum( "line ", input_file.line_count ) +
		input_file.original_line );
  }
  std::string interpolation = input_file.get_text( );
  string_F::Tolower( interpolation );

  Terp::Interp_Type interp = interp_flag_F::read_flag( interpolation, input_file );
  return interp;
}

// *********** interp_flag_F::read_2d_interpolation ***********************
// Interprets the interpolation rules for 2d data
void interp_flag_F::read_2d_interpolation( Dpar::data_parser &input_file, 
  Terp::two_d_interp *energy_in_interp, Terp::Interp_Type *data_out_interp )
{
  for( int i = 0; i < 2; ++i )
  {
    std::string dataID = input_file.get_dataID( );
    string_F::Tolower( dataID );
    if( dataID == "incident energy interpolation" )
    {
      // read the interpolation rule and its qualifier
      std::string interpolation;
      std::string qualifier;
      input_file.read_2d_interp( &interpolation, &qualifier );

      // handle the flag
      energy_in_interp->flag = interp_flag_F::read_flag( interpolation,input_file );

      // handle the qualifier
      energy_in_interp->qualifier = interp_flag_F::read_qualifier( qualifier,
         input_file );
    }
    else if( ( dataID == "outgoing energy interpolation" ) ||
	     ( dataID == "outgoing cosine interpolation" ) )
    {
      std::string interpRule = input_file.get_text( );
      *data_out_interp = interp_flag_F::read_flag( interpRule,input_file );
    }
    else
    {
      Msg::FatalError( "interp_flag_F::read_2d_interpolation",
                "Expected interpolation rule in " +
		Msg::pastenum( "line ", input_file.line_count ) +
		input_file.original_line );
    }
  }
}

// *********** interp_flag_F::read_3d_interpolation_ENDL ***********************
// Interprets the interpolation rules for ENDL 3d data
void interp_flag_F::read_3d_interpolation_ENDL( Dpar::data_parser &input_file,
  Terp::two_d_interp *energy_in_interp, Terp::two_d_interp *cosine_interp,
  Terp::Interp_Type *energy_out_interp )
{
  std::string interpolation;
  std::string qualifier;

  for( int i = 0; i < 3; ++i )
  {
    std::string dataID = input_file.get_dataID( );
    string_F::Tolower( dataID );
    if( dataID == "incident energy interpolation" )
    {
      // read the interpolation rule and its qualifier
      input_file.read_2d_interp( &interpolation, &qualifier );

      // handle the flag
      energy_in_interp->flag = interp_flag_F::read_flag( interpolation, input_file );

      // handle the qualifier
      energy_in_interp->qualifier = interp_flag_F::read_qualifier( qualifier,
         input_file );
    }
    else if( dataID == "outgoing cosine interpolation" )
    {
      // read the interpolation rule and its qualifier
      input_file.read_2d_interp( &interpolation, &qualifier );

      // handle the flag
      cosine_interp->flag = interp_flag_F::read_flag( interpolation, input_file );

      // handle the qualifier
      cosine_interp->qualifier = interp_flag_F::read_qualifier( qualifier,
         input_file );
    }
    else if( dataID == "outgoing energy interpolation" )
    {
      std::string interpRule = input_file.get_text( );
      *energy_out_interp = interp_flag_F::read_flag( interpRule, input_file );
    }
    else
    {
      Msg::FatalError( "interp_flag_F::read_3d_interpolation_ENDL",
                "Expected interpolation rule in " +
		Msg::pastenum( "line ", input_file.line_count ) +
		input_file.original_line );
    }
  }
}

// *********** interp_flag_F::read_3d_interpolation_GND ***********************
// Interprets the interpolation rules for GND 3d data
void interp_flag_F::read_3d_interpolation_GND( Dpar::data_parser &input_file,
  Terp::two_d_interp *energy_in_interp, Terp::two_d_interp *energy_out_interp,
  Terp::Interp_Type *cosine_interp )
{
  std::string interpolation;
  std::string qualifier;

  for( int i = 0; i < 3; ++i )
  {
    std::string dataID = input_file.get_dataID( );
    string_F::Tolower( dataID );
    if( dataID == "incident energy interpolation" )
    {
      // read the interpolation rule and its qualifier
      input_file.read_2d_interp( &interpolation, &qualifier );

      // handle the flag
      energy_in_interp->flag = interp_flag_F::read_flag( interpolation, input_file );

      // handle the qualifier
      energy_in_interp->qualifier = interp_flag_F::read_qualifier( qualifier,
         input_file );
    }
    else if( dataID == "outgoing energy interpolation" )
    {
      // read the interpolation rule and its qualifier
      input_file.read_2d_interp( &interpolation, &qualifier );

      // handle the flag
      energy_out_interp->flag = interp_flag_F::read_flag( interpolation, input_file );

      // handle the qualifier
      energy_out_interp->qualifier = interp_flag_F::read_qualifier( qualifier,
         input_file );
    }
    else if( dataID == "outgoing cosine interpolation" )
    {
      std::string interpRule = input_file.get_text( );
      *cosine_interp = interp_flag_F::read_flag( interpRule, input_file );
    }
    else
    {
      Msg::FatalError( "interp_flag_F::read_3d_interpolation_GND",
                "Expected interpolation rule in " +
		Msg::pastenum( "line ", input_file.line_count ) +
		input_file.original_line );
    }
  }
}

// *********** class Ddvec::Ein_vectors ***********************
// -------------- Ddvec::Ein_vectors::copy_head -----------------
// Copies the j-th initial Ddvec::Ein_entries
void Ddvec::Ein_vectors::copy_head( int j, Ddvec::Ein_vector::const_iterator &Einj_ptr,
			     Ddvec::Ein_vector::const_iterator &Eink_ptr,
			     const Ddvec::Ein_vector &Einj, bool start_zero )
{
  Ddvec::Ein_entries new_entries;
  int k = 1 - j;
  new_entries.n[ k ] = 0;
  static double tol = Global.Value( "looser_tol" );
  double limit = Eink_ptr->Ein - tol*std::abs( Eink_ptr->Ein );
  for( ; ( Einj_ptr->Ein < limit ) && ( Einj_ptr != Einj.end( ) ); ++Einj_ptr )
  {
    new_entries.Ein = Einj_ptr->Ein;
    new_entries.n[ j ] = Einj_ptr->n;
    push_back( new_entries );
  }
  if( !start_zero )
  {
    // remember the initial jump
    start_jump[ k ] = end( );
    --start_jump[ k ];
  }
}
// -------------- Ddvec::Ein_vectors::copy_tail -----------------
// Copies the j-th initial Ddvec::Ein_entries
void Ddvec::Ein_vectors::copy_tail( int j, Ddvec::Ein_vector::const_iterator &Einj_ptr,
			     const Ddvec::Ein_vector &Einj, bool end_zero )
{
  Ddvec::Ein_entries new_entries;
  int k = 1 - j;
  if( !end_zero )
  {
    // remember the final jump
    end_jump[ k ] = end( );
    --end_jump[ k ];
  }
  new_entries.n[ k ] = 0;
  for( ; Einj_ptr != Einj.end( ); ++Einj_ptr )
  {
    new_entries.Ein = Einj_ptr->Ein;
    new_entries.n[ j ] = Einj_ptr->n;
    push_back( new_entries );
  }
}
// -------------- Ddvec::Ein_vectors::common_Ein -----------------
// Produces an Ddvec::Ein_vectors with a common set of incident energies
void Ddvec::Ein_vectors::common_Ein( const Ddvec::Ein_vector &Ein_0, const Ddvec::Ein_vector &Ein_1 )
{
  Ddvec::Ein_vector::const_iterator vect0_ptr = Ein_0.begin( );
  Ddvec::Ein_vector::const_iterator vect1_ptr = Ein_1.begin( );
  static double tol = Global.Value( "looser_tol" );
  static double abs_tol = Global.Value( "tight_tol" );
  double next_Ein = vect1_ptr->Ein;
  Ddvec::Ein_entries new_entries;

  //  default: no jumps
  for( int j = 0; j < 2; ++j )
  {
    start_jump[ j ] = end();
    end_jump[ j ] = end();
  }

  // check for different starting energies
  if( vect0_ptr->Ein < next_Ein - tol*std::abs( next_Ein ) )
  {
    copy_head( 0, vect0_ptr, vect1_ptr, Ein_0, Ein_1.start_zero );
  }
  else
  {
    next_Ein = vect0_ptr->Ein;
    if( vect1_ptr->Ein < next_Ein - tol*std::abs( next_Ein ) )
    {
      copy_head( 1, vect1_ptr, vect0_ptr, Ein_1, Ein_0.start_zero  );
    }
  }

  // scan the overlapping energies
  for( ; ( vect0_ptr != Ein_0.end( ) ) && ( vect1_ptr != Ein_1.end( ) ); )
  {
    double E_in = ( vect0_ptr->Ein < vect1_ptr->Ein ) ? vect0_ptr->Ein : vect1_ptr->Ein;
    double x_tol =  ( E_in == 0.0 ) ? abs_tol : tol * std::abs( E_in );
    if( vect0_ptr->Ein < E_in + x_tol )
    {
      if( vect1_ptr->Ein > E_in + x_tol )
      {
	// copy only  *vect0_ptr
	new_entries.Ein = E_in;
	new_entries.n[ 0 ] = vect0_ptr->n;
	new_entries.n[ 1 ] = 0;
	push_back( new_entries );
      }
      else
      {
	// copy both
	new_entries.Ein = E_in;
	new_entries.n[ 0 ] = vect0_ptr->n;
	new_entries.n[ 1 ] = vect1_ptr->n;
	push_back( new_entries );
	++vect1_ptr;
      }
      ++vect0_ptr;
    }
    else
    {
      // copy only  *vect1_ptr
      new_entries.Ein = E_in;
      new_entries.n[ 0 ] = 0;
      new_entries.n[ 1 ] = vect1_ptr->n;
      push_back( new_entries );
      ++vect1_ptr;
    }
  }

  // check for different final energies
  if( vect0_ptr != Ein_0.end( ) )
  {
    copy_tail( 0, vect0_ptr, Ein_0, Ein_1.end_zero );
  }
  else if( vect1_ptr != Ein_1.end( ) )
  {
    copy_tail( 1, vect1_ptr, Ein_1, Ein_0.end_zero );
  }
}

// *********** Ddvec::dd_entry class ***********************
// -------------- Ddvec::dd_entry::linlin_interp -----------------
// Does linear interpolation
double Ddvec::dd_entry::linlin_interp( double E, const Ddvec::dd_entry& next_entry,
				bool *interp_OK ) const
{
  *interp_OK = true;
  double E_diff = next_entry.x - x;

  // check for division by 0
  if(E_diff == 0.0)
  {
    // Msg::Warning("Ddvec::dd_entry::linlin_interp","division by 0");
    Msg::DebugInfo("Ddvec::dd_entry::linlin_interp", "division by 0");
    *interp_OK = false;
    return y;
  }

  double alpha = (E - x)/E_diff;

  static double E_tol = Global.Value( "looser_tol" );
  // check for extrapolation
  if((alpha < 0.0) || (alpha > 1.0 + E_tol))
  {
    // Msg::Warning("Ddvec::dd_entry::linlin_interp","extrapolation");
    Msg::DebugInfo("Ddvec::dd_entry::linlin_interp","extrapolation");
    *interp_OK = false;
  }
  return (alpha*next_entry.y + (1.0 - alpha)*y);
}
// -------------- Ddvec::dd_entry::linlog_interp -----------------
// Does linear interpolation
double Ddvec::dd_entry::linlog_interp( double E, const Ddvec::dd_entry& next_entry,
				bool *interp_OK ) const
{
  *interp_OK = true;
  double E_ratio = next_entry.x / x;

  // check for division by 0
  if(E_ratio == 1.0)
  {
    // Msg::Warning("Ddvec::dd_entry::linlog_interp","division by 0");
    *interp_OK = false;
    return y;
  }

  double alpha = log(E / x)/log( E_ratio );

  // check for extrapolation
  static double E_tol = Global.Value( "looser_tol" );
  if((alpha < 0.0) || (alpha > 1.0 + E_tol) )
  {
    // Msg::Warning("Ddvec::dd_entry::linlog_interp","extrapolation");
    *interp_OK = false;
  }
  return (alpha*next_entry.y + (1.0 - alpha)*y);
}
// -------------- Ddvec::dd_entry::linlin_interp -----------------
// A first step in bilinear interpolation
bool Ddvec::dd_entry::linlin_interp( double E, double left_data_Ein, const Ddvec::dd_entry& left_data,
    double right_data_Ein, const Ddvec::dd_entry& right_data )
{
  // The x-values are supposed to coincide
  if( left_data.x != right_data.x )
  {
    Msg::Warning("Ddvec::dd_entry::linlin_interp2","diagonal interpolation");
    x = ( left_data.x + right_data.x ) / 2;
  }
  else
  {
    x = left_data.x;
  }

  double E_diff = right_data_Ein - left_data_Ein;
  // check for division by 0
  if( E_diff == 0.0 )
  {
    Msg::Warning("Ddvec::dd_entry::linlin_interp2","division by 0");
    y = ( left_data.y + right_data.y ) / 2;
    return false;
  }

  double alpha = ( E - left_data_Ein )/E_diff;

  // check for extrapolation
  static double e_tol = Global.Value( "looser_tol" );
  if((alpha < 0.0) || (alpha > 1.0 + e_tol ) )
  {
    // Msg::Warning("Ddvec::dd_entry::linlin_interp2","extrapolation");
    return false;
  }
  y = alpha*right_data.y + (1.0 - alpha)*left_data.y;
  return true;
}
// -------------- Ddvec::dd_entry::linlog_interp -----------------
// A first step in bilinear-log interpolation
bool Ddvec::dd_entry::linlog_interp( double E, double left_data_Ein, const Ddvec::dd_entry& left_data,
    double right_data_Ein, const Ddvec::dd_entry& right_data )
{
  // The x-values are supposed to coincide
  if( left_data.x != right_data.x )
  {
    Msg::Warning("Ddvec::dd_entry::linlog_interp2","diagonal interpolation");
    x = ( left_data.x + right_data.x ) / 2;
  }
  else
  {
    x = left_data.x;
  }

  double E_ratio = right_data_Ein / left_data_Ein;
  // check for division by 0
  if( E_ratio == 1.0 )
  {
    Msg::Warning("Ddvec::dd_entry::linlog_interp2","division by 0");
    y = ( left_data.y + right_data.y ) / 2;
    return false;
  }

  double alpha = log( E / left_data_Ein )/log( E_ratio );

  // check for extrapolation
  static double e_tol = Global.Value( "looser_tol" );
  if((alpha < 0.0) || (alpha > 1.0 + e_tol ))
  {
    // Msg::Warning("Ddvec::dd_entry::linlog_interp2","extrapolation");
    return false;
  }
  y = alpha*right_data.y + (1.0 - alpha)*left_data.y;
  return true;
}
// -------------- Ddvec::dd_entry::print -----------------
// for debugging
void Ddvec::dd_entry::print( ) const
{
  std::cout.setf( std::ios::scientific, std::ios::floatfield );
  std::cout << x << "  " << y << std::endl;
}

// *********** Ddvec::dd_pair class ***********************

// -------------- Ddvec::dd_pair::set_pair -----------------
void Ddvec::dd_pair::set_pair( const Ddvec::dd_entry& first_,
  const Ddvec::dd_entry& second_ )
{
  first = first_;
  second = second_;
}
// -------------- Ddvec::dd_pair::value -----------------
double Ddvec::dd_pair::value( double eta, bool *interp_OK )
{
  double x_measure = second.x - first.x;
  static double E_tol = Global.Value( "looser_tol" );
  double x_tol = E_tol * x_measure;
  double ans = 0.0;  // Initialize to stop compilers from printing a warning message.

  if( x_tol <= 0.0 )
  {
    Msg::DebugInfo( "Ddvec::dd_pair::value", "data out of order" );
    *interp_OK = false;
    return 0.0;
  }

  if( Eout_interp == Terp::LINLIN )
  {
    double alpha = ( eta - first.x )/( second.x - first.x );
    if( ( eta > second.x + x_tol ) || ( eta < first.x - x_tol ) )
    {
      ans = 0.0;
      *interp_OK = false;
    }
    else
    {
      ans = ( 1.0 - alpha )*first.y + alpha*second.y;
      *interp_OK = true;
    }
  }
  else if( Eout_interp == Terp::HISTOGRAM )
  {
    ans = first.y;
    *interp_OK = true;
  }
  else
  {
    Msg::FatalError( "Ddvec::dd_pair::value",
		     "interpolation type not implemented" );
  }
  return ans;
}
// -------------- Ddvec::dd_pair::inverse -----------------
// Solves for x given y
double Ddvec::dd_pair::inverse( double Eout, bool *interp_OK )
{
  if( second.y == first.y )
  {
    *interp_OK = false;
    return 0.5 * ( first.x + second.x );
  }
  else
  {
    double slope = ( second.y - first.y ) / ( second.x - first.x );
    *interp_OK = true;
    return first.x + ( Eout - first.y ) / slope;
  }
}
// ----------- Ddvec::dd_pair::find_Ein --------------
// Finds the E_in value at which this eta = const line hits E_out
// Values of *flag:
//  -1, inverse fails
//   0, no intersection
//   1, one intersection
//   2, lines overlap
double Ddvec::dd_pair::find_Ein( double E_out, int *flag )
{
  double Ein_hit = -1.0;

  if( second.y == first.y )
  {
    if( first.y == E_out )
    {
      *flag = 2;
      Ein_hit = first.x;
    }
    else
    {
      *flag = 0;
    }
  }
  else
  {
    *flag = 1;
    bool is_OK;
    Ein_hit = inverse( E_out, &is_OK );
    if( !is_OK ) *flag = -1;
  }
  return Ein_hit;
}
// ----------- Ddvec::dd_pair::linlin_interp --------------
// Does linear-linear interpolation with respect to the tag
bool Ddvec::dd_pair::linlin_interp( double E_in, const Ddvec::dd_pair &lower_data,
		      const Ddvec::dd_pair &higher_data )
{
  set_E_in( E_in );
  bool is_OKL = first.linlin_interp( E_in, lower_data.get_E_in( ),
    lower_data.first, higher_data.get_E_in( ), higher_data.first);
  bool is_OKR = second.linlin_interp( E_in, lower_data.get_E_in( ),
    lower_data.second, higher_data.get_E_in( ), higher_data.second);
  return ( is_OKL && is_OKR );
}
// ----------- Ddvec::dd_pair::linlog_interp --------------
// Does linear-log interpolation with respect to the tag
bool Ddvec::dd_pair::linlog_interp( double E_in, const Ddvec::dd_pair &lower_data,
		      const Ddvec::dd_pair &higher_data )
{
  set_E_in( E_in );
  bool is_OKL = first.linlog_interp( E_in, lower_data.get_E_in( ),
    lower_data.first, higher_data.get_E_in( ), higher_data.first);
  bool is_OKR = second.linlog_interp( E_in, lower_data.get_E_in( ),
    lower_data.second, higher_data.get_E_in( ), higher_data.second);
  return ( is_OKL && is_OKR );
}
// ----------- Ddvec::dd_pair::set_data --------------
// Sets up the data for outgoing energies Eout_min and Eout_max
bool Ddvec::dd_pair::set_data( const Ddvec::dd_entry &prev_data, const Ddvec::dd_entry &next_data,
		 double Eout_min, double Eout_max )
{
  if( Eout_min > Eout_max )
  {
    Msg::DebugInfo( "Ddvec::dd_pair::set_data", "improper energy range" );
    return false;
  }
  static double E_tol = Global.Value( "looser_tol" );
  if( Eout_min < ( 1.0 - E_tol )*prev_data.x )
  {
    Msg::DebugInfo( "Ddvec::dd_pair::set_data", "Eout_min too low" );
    return false;
  }
  if( Eout_max > ( 1.0 + E_tol )*next_data.x )
  {
    Msg::DebugInfo( "Ddvec::dd_pair::set_data", "Eout_max too high" );
    return false;
  }

  bool is_OK = true;
  if( Eout_interp == Terp::HISTOGRAM )
  {
    first.x = Eout_min;
    first.y = prev_data.y;
    second.x = Eout_max;
    // second.y is not used for histograms in E_out
    second.y = next_data.y;
  }
  else
  {
    bool is_OKL = true;
    bool is_OKR = true;
    first.x = Eout_min;
    second.x = Eout_max;
    if( Eout_min < ( 1.0 + E_tol )*prev_data.x )
    {
      first.y = prev_data.y;
    }
    else
    {
      first.y = prev_data.linlin_interp( Eout_min, next_data, &is_OKL );
    }
    if( Eout_max > ( 1.0 - E_tol )*next_data.x )
    {
      second.y = next_data.y;
    }
    else
    {
      second.y = prev_data.linlin_interp( Eout_max, next_data, &is_OKR );
    }
    is_OK = is_OKL && is_OKR;
  }

  return is_OK;
}

// ************* class Ddvec::unit_base_map *****************
// ----------- Ddvec::unit_base_map::un_unit_base --------------
// Undoes the unit-base map
// finds Eout for 0 <= eta <= 1
double Ddvec::unit_base_map::un_unit_base( double eta ) const
{
  static double etol = Global.Value( "tight_tol" );
  if( ( eta < 0.0 ) || ( eta > 1.0 + etol ) )
  {
    Msg::FatalError(  "Ddvec::unit_base_map::un_unit_base",
                Msg::pastenum( "eta: ", eta ) +
                " out of range" );
  }


  double E_out = Eout_min + eta * ( Eout_max - Eout_min );
  if( E_out > Eout_max ) E_out = Eout_max;
  return E_out;
}
// ----------- Ddvec::unit_base_map::un_unit_base --------------
// Undoes the unit-base map
// finds Eout for 0 <= eta <= 1
double Ddvec::unit_base_map::un_unit_base( double eta, bool *interp_OK ) const
{
  static double etol = Global.Value( "tight_tol" );
  if( ( eta < 0.0 ) || ( eta > 1.0 + etol ) )
  {
    Msg::DebugInfo(  "Ddvec::unit_base_map::un_unit_base",
                Msg::pastenum( "eta: ", eta ) +
                " out of range" );
    *interp_OK= false;
    if( eta < 0.0 )
    {
      return 0.0;
    }
    else
    {
      return 1.0;
    }   
  }

  *interp_OK = true;
  double E_out = Eout_min + eta * ( Eout_max - Eout_min );
  if( E_out > Eout_max ) E_out = Eout_max;
  return E_out;
}
// ----------- Ddvec::unit_base_map::to_unit_base --------------
// Finds eta for given physical E
double Ddvec::unit_base_map::to_unit_base( double physical_E, bool *interp_OK )
{
  static double etol = Global.Value( "tight_tol" );
  // When used with direction cosines, "Eout_min" could be mu = -1
  if( physical_E - Eout_min < - etol*std::abs( Eout_min ) )
  {
    *interp_OK = false;
    return 0.0;
  }
  else if( physical_E - Eout_max > etol*std::abs( Eout_max ) )
  {
    Msg::DebugInfo(  "Ddvec::unit_base_map::to_unit_base",
                 Msg::pastenum( "physical_E: ", physical_E ) +
                 " out of range" );
    *interp_OK = false;
    return 1.0;
  }
  double eta = ( physical_E - Eout_min ) / ( Eout_max - Eout_min );
  if( eta < 0.0 ) eta = 0.0;
  if( eta > 1.0 ) eta = 1.0;
  *interp_OK = true;
  return eta;
}
// ----------- Ddvec::unit_base_map::copy --------------
// Makes a copy
void Ddvec::unit_base_map::copy( const Ddvec::unit_base_map &to_copy )
{
  Eout_min = to_copy.Eout_min;
  Eout_max = to_copy.Eout_max;
}
// ----------- Ddvec::unit_base_map::interpolate --------------
// Interpolates the energy range
bool Ddvec::unit_base_map::interpolate( double alpha, const Ddvec::unit_base_map &prev,
  const Ddvec::unit_base_map &next )
{
  bool is_OK = ( alpha >= 0.0 ) && ( alpha <= 1.0 );
  Eout_min = ( 1.0 - alpha )*prev.Eout_min + alpha*next.Eout_min;
  Eout_max = ( 1.0 - alpha )*prev.Eout_max + alpha*next.Eout_max;
  return is_OK;
}
// ----------- Ddvec::unit_base_map::too_short --------------
// Tests for a trivial interval
bool Ddvec::unit_base_map::too_short( )
{
  bool Short = ( Eout_min >= Eout_max );
  return Short;
}

// *********** Ddvec::cum_points_pair ***********************
// ------------- Ddvec::cum_points_pair::too_short -------------------
// Tests for a trivial interval
bool Ddvec::cum_points_pair::too_short( )
{
  return ( first.x >= second.x );
}
// -------------- Ddvec::cum_points_pair::to_unit_base -----------------
// Maps the pair to unit base
void Ddvec::cum_points_pair::to_unit_base( )
{
  double Eout_min = first.x;
  double Eout_max = second.x;
  double scale = Eout_max - Eout_min;
  if( scale <= 0.0 )
  {
    Msg::FatalError( "Ddvec::cum_points_pair::to_unit_base",
                "bad energy range" );
  }
  ubase_map.Eout_min = Eout_min;
  ubase_map.Eout_max = Eout_max;

  first.x = 0.0;
  first.y *= scale;
  second.x = 1.0;
  second.y *= scale;
}
// -------------- Ddvec::cum_points_pair::short_to_unit_base -----------------
// Maps the pair to unit base, used for intervals of length zero
void Ddvec::cum_points_pair::short_to_unit_base( double dA 	)
{
  double Eout_min = first.x;
  double Eout_max = second.x;
  double scale = dA;
 
  ubase_map.Eout_min = Eout_min;
  ubase_map.Eout_max = Eout_max;

  first.x = 0.0;
  first.y *= scale;
  second.x = 1.0;
  second.y *= scale;
}
// -------------- Ddvec::cum_points_pair::un_unit_base -----------------
// Undoes the mapping to unit base
void Ddvec::cum_points_pair::un_unit_base( )
{
  double scale = ubase_map.Eout_max - ubase_map.Eout_min;
  if( scale <= 0.0 )
  {
    Msg::FatalError( "Ddvec::cum_points_pair::un_unit_base",
                "bad energy range" );
  }

  first.x = ubase_map.Eout_min;
  first.y *= 1.0/scale;
  second.x = ubase_map.Eout_max;
  second.y *= 1.0/scale;
}


// *********** Ddvec::dd_vector class ***********************
// -------------- Ddvec::dd_vector::dd_vector -----------------
// Copy constructor
Ddvec::dd_vector::dd_vector( const Ddvec::dd_vector &to_copy )
{
  tag = to_copy.tag;
  interp_type = to_copy.interp_type;
  for( Ddvec::dd_vector::const_iterator copy_ptr = to_copy.begin( );
       copy_ptr != to_copy.end( ); ++copy_ptr )
  {
    Ddvec::dd_entry copy_entry( *copy_ptr );
    push_back( copy_entry );
  }
}
// -------------- Ddvec::dd_vector::~dd_vector -----------------
// Destructor
Ddvec::dd_vector::~dd_vector( )
{
}
// -------------- Ddvec::dd_vector::print -----------------
// to print the vector
void Ddvec::dd_vector::print() const
{
  std::cout.setf( std::ios::scientific, std::ios::floatfield);

  // first print the tag
  std::cout << "tag: " << tag << std::endl;

  // loop through the vector
  for(Ddvec::dd_vector::const_iterator entry = begin(); entry != end();
    ++entry)
  {
    std::cout << entry->x << "  "
      << entry->y << std::endl;
  }
  std::cout << std::endl;
}
// --------------- Ddvec::dd_vector::read_data --------------
// Reads the only the data from a file
void Ddvec::dd_vector::read_data( Dpar::data_parser &input_file, int num_sigma )
{
  for( int E_count = 0; E_count < num_sigma; ++E_count )
  {
    Ddvec::dd_entry next_xy;
    next_xy.x = input_file.get_next_double( );
    next_xy.y = input_file.get_next_double( );
    push_back( next_xy );
  }
}
// --------------- Ddvec::dd_vector::read_data_interp --------------
// Reads the data and interpolation type from a file
void Ddvec::dd_vector::read_data_interp( Dpar::data_parser &input_file, int num_sigma )
{
  std::string dataID;
  interp_type = interp_flag_F::read_1d_interpolation( input_file );
  read_data( input_file, num_sigma );
}
// ----------- Ddvec::dd_vector::operator*= ------------
Ddvec::dd_vector& Ddvec::dd_vector::operator*=( double factor )
{
  for( Ddvec::dd_vector::iterator L1 = begin( );
       L1 != end( ); ++L1 )
  {
    L1->y *= factor;
  }

  return *this;
}
// ----------- Ddvec::dd_vector::operator+= ------------
Ddvec::dd_vector& Ddvec::dd_vector::operator+=( Ddvec::dd_vector& vector_2 )
{
  Ddvec::dd_vector::iterator L1 = begin( );

  // add the y values
  for( Ddvec::dd_vector::iterator L2 = vector_2.begin( );
       ( L1 != end( ) ) && ( L2 != vector_2.end( ) );
       ++L1, ++L2 )
  {
    if( L1->x != L2->x ){
      Msg::FatalError( "Ddvec::dd_vector::operator+=", 
		   Msg::pastenum( "different x values: ", L1->x ) +
		   Msg::pastenum( " and: ", L2->x ) );
    }
    L1->y += L2->y;
  }

  return *this;
}
// ----------- Ddvec::dd_vector::operator*= ------------
Ddvec::dd_vector& Ddvec::dd_vector::operator*=( Ddvec::dd_vector& vector_2 )
{
  Ddvec::dd_vector::iterator L1 = begin( );

  // add the y values
  for( Ddvec::dd_vector::iterator L2 = vector_2.begin( );
       ( L1 != end( ) ) && ( L2 !=  vector_2.end( ) );
       ++L1, ++L2 )
  {
    if( L1->x != L2->x ){
      Msg::FatalError( "Ddvec::dd_vector::operator*=", 
		   Msg::pastenum( "different x values: ", L1->x ) +
		   Msg::pastenum( " and: ", L2->x ) );
    }
    L1->y *= L2->y;
  }

  return *this;
}
// ----------- Ddvec::dd_vector::add_entry --------------
// append an ( E_out, Prob ) entry to the vector
void Ddvec::dd_vector::add_entry( double E_out, double Prob )
{
  Ddvec::dd_entry eout_prob_data;  //  the next (E_out, probability) entry

  eout_prob_data.x = E_out;
  eout_prob_data.y = Prob;
  push_back( eout_prob_data );
}
//---------------- Ddvec::dd_vector::make_flat ---------------------
// Makes a constant vector
void Ddvec::dd_vector::make_flat( const Ddvec::dd_vector& sigma, double val )
{
  Ddvec::dd_vector::const_iterator this_sigma = sigma.begin( );
  add_entry( this_sigma->x, val );
  this_sigma = sigma.end( );
  --this_sigma;
  add_entry( this_sigma->x, val );
  interp_type = Terp::HISTOGRAM;
}
//---------------- Ddvec::dd_vector::make_flat ---------------------
// Makes a constant vector
void Ddvec::dd_vector::make_flat( const vector< double >& E_groups, double val )
{
  vector< double >::const_iterator this_group_bd = E_groups.begin( );
  add_entry( *this_group_bd, val );
  this_group_bd = E_groups.end( );
  --this_group_bd;
  add_entry( *this_group_bd, val );
  interp_type = Terp::HISTOGRAM;
}
// ----------- Ddvec::dd_vector::copy ------------
// Copies a vector
void Ddvec::dd_vector::copy( const Ddvec::dd_vector& vector_from )
{
  tag = vector_from.tag;
  interp_type = vector_from.interp_type;
  for( Ddvec::dd_vector::const_iterator L2 = vector_from.begin( );
        L2 != vector_from.end( ); ++L2 )
  {
    add_entry( L2->x, L2->y );
  }
}
// ----------- Ddvec::dd_vector::extrapolate_copy ------------
// Copies a vector with extrapolation
void Ddvec::dd_vector::extrapolate_copy( const Ddvec::dd_vector& vector_from, double min_x,
  double max_x )
{
  tag = vector_from.tag;
  interp_type = vector_from.interp_type;
  static double abs_tol = Global.Value( "tight_tol" );
  Ddvec::dd_vector::const_iterator L2 = vector_from.begin( );
  double x_value = L2->x;
  if( min_x < ( 1.0 - abs_tol )*x_value )
  {
    add_entry( min_x, 0.0 );
    if( L2->y != 0.0 ) add_entry( ( 1.0 - abs_tol )*x_value, 0.0 );
    add_entry( x_value, L2->y );
  }
  else if( min_x < x_value )
  {
    add_entry( min_x, L2->y );
  }
  else
  {
    add_entry( x_value, L2->y );
  }
  for( ++L2; L2 != vector_from.end( ); ++L2 )
  {
    add_entry( L2->x, L2->y );
  }
  Ddvec::dd_vector::iterator Lptr = end( );
  --Lptr;
  x_value = Lptr->x;
  if( max_x > ( 1.0 + abs_tol )*x_value )
  {
    if( Lptr->y != 0.0 ) add_entry( ( 1.0 + abs_tol )*x_value, 0.0 );
    add_entry( max_x, 0.0 );
  }
  else if( max_x > x_value )
  {
    Lptr->x = max_x;
  }
}
// ----------- Ddvec::dd_vector::truncate_copy ------------
// Copies a vector with truncation
bool Ddvec::dd_vector::truncate_copy( const Ddvec::dd_vector& vector_from, double min_x,
			       double max_x, bool do_renorm )
{
  if( min_x >= max_x )
  {
    return false;
  }
  
  tag = vector_from.tag;
  interp_type = vector_from.interp_type;
  static double abs_tol = Global.Value( "tight_tol" );
  Ddvec::dd_vector::const_iterator prev = vector_from.begin( );
  Ddvec::dd_vector::const_iterator next = prev;
  ++next;
  bool interp_OK;
  
  // find the start
  while( next->x < ( 1.0 + abs_tol ) * min_x )
  {
    prev = next;
    ++next;
  }
  if( interp_type == Terp::HISTOGRAM )
  {
    add_entry( min_x, prev->y );
  }
  else  // Terp::LINLIN
  {
    if( prev->x < ( 1.0 - abs_tol ) * min_x )
    {
      double y_value = prev->linlin_interp( min_x, *next, &interp_OK );
      add_entry( min_x, y_value );
    }
    else
    {
      add_entry( min_x, prev->y );
    }
  }
  // do the middle ones
  while( next->x < ( 1.0 - abs_tol ) * max_x )
  {
    add_entry( next->x, next->y );
    prev = next;
    ++next;
  }
  // do the last one
  if( interp_type == Terp::HISTOGRAM )
  {
    add_entry( max_x, prev->y );
  }
  else  // Terp::LINLIN
  {
    if( prev->x < ( 1.0 - abs_tol ) * max_x )
    {
      double y_value = prev->linlin_interp( max_x, *next, &interp_OK );
      add_entry( max_x, y_value );
    }
    else
    {
      add_entry( max_x, prev->y );
    }
  }

  // renorm truncated data?
  bool norm_OK = true;
  if( do_renorm )
  {
    norm_OK = renorm( true );
  }

  return norm_OK;
}
// ----------- Ddvec::dd_vector::prepend ------------
// Prepends x=0 using histogram extrapolation if necessary
void Ddvec::dd_vector::prepend( const Ddvec::dd_vector& to_copy, double E_insert )
{
  tag = to_copy.tag;
  interp_type = to_copy.interp_type;

  // prepend if necessary
  Ddvec::dd_vector::const_iterator this_entry = to_copy.begin( );
  if( this_entry->x > 0.0 )
  {
    add_entry( 0.0, this_entry->y );
  }

  // make sure the result has an entry at E_insert
  if( E_insert > 0.0 )
  {
    static double tol = Global.Value( "looser_tol" );
    // copy up to E_insert
    for( ; this_entry->x < E_insert*( 1.0 - tol ); ++this_entry )
    {
      add_entry( this_entry->x, this_entry->y );
    }
    if( this_entry->x > E_insert*( 1.0 + tol ) )
    {
      Ddvec::dd_vector::const_iterator prev_entry = this_entry;
      --prev_entry;
      bool is_OK;
      double mid_y = prev_entry->linlin_interp( E_insert, *this_entry, &is_OK );
      add_entry( E_insert, mid_y );
    }
  }

  // copy the rest
  for( ; this_entry != to_copy.end( ); ++this_entry )
  {
    add_entry( this_entry->x, this_entry->y );
  }
}
// ----------- Ddvec::dd_vector::parse_vector ------------
// find the incident energies and the jumps
void Ddvec::dd_vector::parse_vector( Ddvec::Ein_vector *E_in_list ) const
{
  static double tol = Global.Value( "looser_tol" );
  Ddvec::Ein_entry new_entry;
  Ddvec::dd_vector::const_iterator vect_ptr = begin( );
  E_in_list->start_zero = ( vect_ptr->y == 0 );
  for( ; vect_ptr != end( ); )
  {
    double E_in = vect_ptr->x;
    double x_tol = tol * std::abs( E_in );
    Ddvec::dd_vector::const_iterator next_ptr = vect_ptr;
    ++next_ptr;
    double next_E_in;
    if ( next_ptr != end( ) )
    {
      next_E_in = next_ptr->x;
    }
    new_entry.Ein = E_in;
    new_entry.n = 1;
    while( ( next_ptr != end( ) ) && ( next_E_in < E_in + x_tol ) )
    {
      ++new_entry.n;
      ++next_ptr;
      next_E_in = next_ptr->x;
    }
    if( new_entry.n > 2 )
    {
      Msg::Warning( "Ddvec::dd_vector::parse_vector",
		    Msg::pastenum( "omitting values for E_in: ", E_in ) );
    }
    E_in_list->push_back( new_entry );
    vect_ptr = next_ptr;
  }
  // is the last entry zero?
  vect_ptr = end( );
  --vect_ptr;
  E_in_list->end_zero = ( vect_ptr->y == 0 );
}
// ----------- Ddvec::dd_vector::fill_with ------------
// compute fill_vector containing entries for each energy in E_in_list
void Ddvec::dd_vector::fill_with( const Ddvec::Ein_vectors &E_in_list, int j,
  Ddvec::dd_vector *fill_vector ) const
{
  if( ( interp_type != Terp::LINLIN ) && ( interp_type != Terp::HISTOGRAM ) )
  {
    Msg::FatalError( "Ddvec::dd_vector::fill_with",
		     "interpolation type not implemented" );
  }

  int k = 1 - j;
  static double tol = Global.Value( "looser_tol" );
  Ddvec::dd_entry new_entry;  // for creating new entries
  // insert zeros at the start of the vector
  Ddvec::Ein_vectors::const_iterator E_in_ptr = E_in_list.begin( );
  double E_in = E_in_ptr->Ein;
  double x_tol = tol * std::abs( E_in );
  Ddvec::dd_vector::const_iterator this_entry = begin( );
  Ddvec::dd_vector::const_iterator prev_entry = this_entry;
  if( E_in < this_entry->x - x_tol )
  {
    for( ; E_in < this_entry->x - x_tol; )
    {
      new_entry.x = E_in_ptr->Ein;
      new_entry.y = 0.0;
      fill_vector->push_back( new_entry );
      // the other vector may have a jump
      if( E_in_ptr->n[ k ] > 1 )
      {
        fill_vector->push_back( new_entry );
      }
      ++E_in_ptr;
      E_in = E_in_ptr->Ein;
      x_tol = tol * std::abs( E_in );
    }
    if( ( E_in_ptr->Ein <= this_entry->x + x_tol ) && ( this_entry->y != 0.0 ) )
    {
      new_entry.x = E_in_ptr->Ein;  // insert a jump
      new_entry.y = 0.0;
      fill_vector->push_back( new_entry );
    }
  }
  else if( E_in_ptr->Ein <= this_entry->x + x_tol )
  {
    new_entry.x = E_in_ptr->Ein;
    new_entry.y = this_entry->y;
    fill_vector->push_back( new_entry );
    ++this_entry;
  }

  // copy the rest of the vector, with insertions according to E_in_list
  for( ++E_in_ptr; E_in_ptr != E_in_list.end( ); ++E_in_ptr )
  {
    E_in = E_in_ptr->Ein;
    x_tol = tol * std::abs( E_in );
    if( this_entry == end( ) )
    {
      // we came to the end of the vector
      break;
    }
    else if( this_entry->x <= E_in + x_tol )
    {
      // there is a node here; just copy and go on to the next
      new_entry.x = E_in;
      new_entry.y = this_entry->y;
      fill_vector->push_back( new_entry );
      if( E_in_ptr->n[ j ] > 1 )  // duplicate incident energies
      {
        if( E_in_ptr->n[ j ] > 2 )  // skip really close entries
        {
          for( int ell = 1; ell < E_in_ptr->n[ j ] - 1; ++ell )
          {
            ++this_entry;
	  }
        }
	++this_entry;
        new_entry.x = E_in;
        new_entry.y = this_entry->y;
        fill_vector->push_back( new_entry );
      }
      prev_entry = this_entry;
      ++this_entry;
    }
    else
    {
      // make a new entry
      new_entry.x = E_in;
      if( interp_type == Terp::LINLIN )
      {
	bool is_OK;
        new_entry.y = prev_entry->linlin_interp( E_in, *this_entry, &is_OK );
      }
      else
      {
        new_entry.y = prev_entry->y;  // histogram data
      }
      fill_vector->push_back( new_entry );
    }
    if( ( E_in_ptr->n[ k ] > 1 ) ||
        ( E_in_ptr == E_in_list.start_jump[ k ] ) ||
        ( E_in_ptr == E_in_list.end_jump[ k ] ) )    // duplicate incident energies
    {
      fill_vector->push_back( new_entry );
    }
  }

  // we may need to insert zeros at the tail
  if( E_in_ptr != E_in_list.end( ) )
  {
    // we may need to add a jump
    if( prev_entry->y != 0.0 )
    {
      new_entry.x = prev_entry->x;  // insert a jump
      new_entry.y = 0.0;
      fill_vector->push_back( new_entry );
    }
    for( ; E_in_ptr != E_in_list.end( ); ++E_in_ptr )
    {
      new_entry.x = E_in_ptr->Ein;
      new_entry.y = 0.0;
      fill_vector->push_back( new_entry );
      if( E_in_ptr->n[ k ] > 1 )
      {
        fill_vector->push_back( new_entry );
      }
    }
  }
}
// ----------- Ddvec::dd_vector::get_norm -----------------
// get the probability for this vector.
// This routine is used for double differential data
double Ddvec::dd_vector::get_norm()
{
  Ddvec::dd_vector::iterator entry = begin();  // the first entry
  if( interp_type == Terp::NOTSET )
  {
    Msg::FatalError( "Ddvec::dd_vector::get_norm",
		     "interpolation type not set" );
  }
  double norm = 0.0;
  double last_E = entry->x;
  double last_P = entry->y;
  double next_E;
  double next_P;

// loop through the vector
  for(++entry; entry != end(); ++entry)
  {
    // get the next energy and probability density
    next_E = entry->x;
    next_P = entry->y;

    // accumulate the probability
    if( interp_type == Terp::LINLIN )
    {
      norm += 0.5*(next_P + last_P)*(next_E - last_E);
    }
    else if( interp_type == Terp::HISTOGRAM )
    {
      norm += last_P*(next_E - last_E);
    }
    else
    {
      Msg::FatalError( "Ddvec::dd_vector::get_norm",
		       "interpolation type not implemented" );
    }

    // save the values for the next interval
    last_E = next_E;
    last_P = next_P;
  }

  return norm;
}
// ----------- Ddvec::dd_vector::renorm ------------
// For probability-density tables set the norm to 1;
// Returns false if we get zero norm from truncated direct interpolation
bool Ddvec::dd_vector::renorm( bool truncated )
{
  double norm = get_norm( );
  Ddvec::dd_vector::iterator entry = begin( );

  // error check
  if( norm == 0.0 )
  {
    if( truncated )
    {
      // direct interpolation with truncation gave norm 0
      return false;
    }
    else
    {
      Msg::Warning("Ddvec::dd_vector::renorm",
        Msg::pastenum( "zero norm in renorm routine for E_in: ", tag ) );
      // fudge something
      entry->y = 1.0;
    }
  }
  
  static int norm_warn = Global.Value( "norm_warn" );
  static double norm_tol = Global.Value( "norm_tol" );
  static bool norm_good = true;
  if( norm_warn && norm_good && ( std::abs( norm - 1.0 ) > norm_tol ) )
  {
    Msg::Warning("Ddvec::dd_vector::renorm",
            Msg::pastenum( "bad norm in renorm routine: ", norm ) );
    norm_good = false;
  }

  // Now, renormalize
  for( ; entry != end( ); ++entry )
  {
    entry->y /= norm;
  }
  
  // this is not truncated data with zero norm
  return true;
}
// ----------- Ddvec::dd_vector::mapto_01 ------------
// Maps a vector to [0, 1] and saves the map; we don't always want to normalize
void Ddvec::dd_vector::mapto_01( Ddvec::unit_base_map *ubase_map )
{
  // find the energy range
  Ddvec::dd_vector::iterator Eprob_ptr = begin( );
  double Eout_min = Eprob_ptr->x;
  ubase_map->Eout_min = Eout_min;
  Ddvec::dd_vector::iterator Eprob_last = end( );
  --Eprob_last;
  double Eout_max = Eprob_last->x;
  ubase_map->Eout_max = Eout_max;
  double scale = Eout_max - Eout_min;
  if( scale <= 0.0 )
  {
    Msg::FatalError( "Ddvec::dd_vector::mapto_01",
                 "improper energy range" );
  }

  for( ; Eprob_ptr != end( ); ++Eprob_ptr )
  {
    Eprob_ptr->x = ( Eprob_ptr->x - Eout_min ) / scale;
  }
}
// ----------- Ddvec::dd_vector::unit_base ------------
// map a vector of probability densities to unit base, and save the map
void Ddvec::dd_vector::unit_base( bool Renorm, Ddvec::unit_base_map *ubase_map )
{
  // First, map the x-values to [0, 1]
  mapto_01( ubase_map );
  // if Renorm, make sure that the norm is 1
  if( Renorm )
  {
    renorm( false );
  }
  else
  {
    // use the computed scale
    double Eout_min = ubase_map->Eout_min;
    double scale = ubase_map->Eout_max - Eout_min;
    *this *= scale;
  }
}
// ----------- Ddvec::dd_vector::un_mapto_01 ------------
// Inverts mapto_01
void Ddvec::dd_vector::un_mapto_01( const Ddvec::unit_base_map *ubase_map )
{
  double Eout_min = ubase_map->Eout_min;
  double scale = ubase_map->Eout_max - Eout_min;
  if( scale <= 0.0 )
  {
    Msg::FatalError( "Ddvec::dd_vector::un_mapto_01",
                 "improper energy range" );
  }

  for( Ddvec::dd_vector::iterator Eprob_ptr = begin( ); Eprob_ptr != end( ); ++Eprob_ptr )
  {
    Eprob_ptr->x = Eout_min + scale*Eprob_ptr->x;
  }
}
// ----------- Ddvec::dd_vector::un_unit_base ------------
// invert the unit-base map of a vector of probability densities
void Ddvec::dd_vector::un_unit_base( bool Renorm, const Ddvec::unit_base_map *ubase_map )
{
  // first, map the x-values
  un_mapto_01( ubase_map );
  // if Renorm, make sure that the norm is 1
  if( Renorm )
  {
    renorm( false );
  }
  else
  {
    // use the computed scale
    double Eout_min = ubase_map->Eout_min;
    double scale = ubase_map->Eout_max - Eout_min;
    *this *= 1.0/scale;
  }
}
// ----------- Ddvec::dd_vector::interpolate ------------
// Do a linear interpolation of 2 vectors
void Ddvec::dd_vector::interpolate( double E_in, const Ddvec::dd_vector &prev_vect, const Ddvec::dd_vector &next_vect )
{
  // work with filled copies
  Ddvec::dd_vector prev_vect_fill;
  Ddvec::dd_vector next_vect_fill;
  dd_vector_F::fill_in_vectors( prev_vect, next_vect, &prev_vect_fill,
     &next_vect_fill, false );
  filled_interpolate( E_in, prev_vect_fill, next_vect_fill );
  interp_type = prev_vect.interp_type;
}
// ----------- Ddvec::dd_vector::filled_interpolate ------------
// Does a linear interpolation of 2 vectors with common x-values
void Ddvec::dd_vector::filled_interpolate( double E_in, const Ddvec::dd_vector &prev_filled,
  const Ddvec::dd_vector &next_filled )
{
  double alpha = ( E_in - prev_filled.get_E_in( ) ) /
    ( next_filled.get_E_in( ) - prev_filled.get_E_in( ) );
  if( ( alpha < 0.0 ) || ( alpha > 1.0 ) )
  {
    Msg::Warning( "Ddvec::dd_vector::filled_interpolate", "extrapolation" );
  }
  Ddvec::dd_entry new_entry;
  Ddvec::dd_vector::const_iterator prev_ptr = prev_filled.begin( );
  Ddvec::dd_vector::const_iterator next_ptr = next_filled.begin( );
  for( ; ( ( prev_ptr != prev_filled.end( ) ) &&
	   ( next_ptr != next_filled.end( ) ) ); ++prev_ptr, ++next_ptr )
  {
    new_entry.x = prev_ptr->x;
    new_entry.y = ( 1.0 - alpha )*prev_ptr->y + alpha*next_ptr->y;
    push_back( new_entry );
  }
  tag = E_in;
}
// ----------- Ddvec::dd_vector::chop_histogram ------------
// Truncates histogram data at the maximum energy
void Ddvec::dd_vector::chop_histogram( double EMax )
{
  Ddvec::dd_vector::iterator last_link = end( );
  --last_link;
  Ddvec::dd_vector::iterator next_last = last_link;
  if( last_link->x > EMax )
  {
    bool chop = false;
    --next_last;
    while( ( next_last->x > EMax ) &&
	   ( next_last != begin( ) ) )
    {
      chop = true;
      last_link = next_last;
      --next_last;
    }
    // do we need to chop?
    next_last = last_link;
    ++next_last;  // the new last link
    if( chop )
    {
      erase( next_last, end( ) );
    }
  }
  last_link->y = 0.0;  // set to zero, whether we chop or not
}
// ----------- Ddvec::dd_vector::scale_E ------------
// Scales the x-component; a kludge to convert ENDF eV to MeV
void Ddvec::dd_vector::scale_E( double eV_to_MeV )
{
  for( Ddvec::dd_vector::iterator this_link = begin( ); this_link != end( );
       ++this_link )
  {
    this_link->x *= eV_to_MeV;
  }
}
// ----------- Ddvec::dd_vector::isotropic ------------
// Checks whether an angular probability density is isotropic
bool Ddvec::dd_vector::isotropic( ) const
{
  bool iso = true;
  static double abs_tol = Global.Value( "tight_tol" );
  for( Ddvec::dd_vector::const_iterator link = begin( ); link != end( );
       ++link )
  {
    if( std::abs( link->y - 0.5 ) > abs_tol )
    {
      iso = false;
      break;
    }
  }
  return iso;
}
// ----------- Ddvec::dd_vector::locate_x ------------
// Finds the interval containing the desired X-value
void Ddvec::dd_vector::locate_x( double X, Ddvec::dd_vector::const_iterator *Prev_ptr,
	       Ddvec::dd_vector::const_iterator *Next_ptr ) const
{
  static double tol = Global.Value( "tight_tol" );

  Ddvec::dd_vector::const_iterator prev_ptr = begin( );
  Ddvec::dd_vector::const_iterator next_ptr = prev_ptr;
  ++next_ptr;

  if( (( prev_ptr->x >= 0.0 ) && ( X < (1-tol) * prev_ptr->x )) ||
      (( prev_ptr->x < 0.0 ) && ( X < (1+tol) * prev_ptr->x )) )
  {
    Msg::FatalError( "Ddvec::dd_vector::locate_x",
		     "X is below the range" );
  }

  while( ( next_ptr->x <= X ) && ( next_ptr != end( ) ) )
  {
    prev_ptr = next_ptr;
    ++next_ptr;
  }

  // Did we overshoot?
  if( next_ptr == end( ) )
  {
    next_ptr = prev_ptr;
    --prev_ptr;
    if( (( X >= 0.0 ) && ( next_ptr->x < (1+tol) * X )) ||
	(( X < 0.0 ) && ( next_ptr->x < (1-tol) * X )) )
    {
      Msg::FatalError( "Ddvec::dd_vector::locate_x",
		       "X is above the range" );
    }
  }

  *Prev_ptr = prev_ptr;
  *Next_ptr = next_ptr;
}
// **************** fill_in_vectors *************************
void dd_vector_F::fill_in_vectors( const Ddvec::dd_vector& vector_0,
  const Ddvec::dd_vector& vector_1, Ddvec::dd_vector *vector_0_fill,
  Ddvec::dd_vector *vector_1_fill, bool re_norm )
// This routine computes vector_0_fill and vector_1_fill with common x-values.
{
  if( vector_0.empty( ) || vector_1.empty( ) )
  {
    std::cout << "vector_0\n";
    vector_0.print( );
    std::cout << "vector_1\n";
    vector_1.print( );
    Msg::FatalError("dd_vector_F::fill_in_vectors",
      "attempt to fill in an empty Ddvec::dd_vector");
  }
  // copy the tags
  vector_0_fill->set_E_in( vector_0.get_E_in( ) );
  vector_1_fill->set_E_in( vector_1.get_E_in( ) );
  if( vector_0.interp_type != vector_1.interp_type )
  {
    Msg::FatalError( "dd_vector_F::fill_in_vectors",
       "inconsistent interpolation types" );
  }
  vector_0_fill->interp_type = vector_0.interp_type;
  vector_1_fill->interp_type = vector_0.interp_type;
  // make a list of the common E_in values
  Ddvec::Ein_vector E_in_list0;
  vector_0.parse_vector( &E_in_list0 );
  Ddvec::Ein_vector E_in_list1;
  vector_1.parse_vector( &E_in_list1 );
  Ddvec::Ein_vectors Ein_common;
  Ein_common.common_Ein( E_in_list0, E_in_list1 );

  // make copies, with jumps
  vector_0.fill_with( Ein_common, 0, vector_0_fill );
  vector_1.fill_with( Ein_common, 1, vector_1_fill );
  if( re_norm )  // this is not done for the I=4 higher Legendre terms
  {
    vector_0_fill->renorm( false );
    vector_1_fill->renorm( false );
  }
  Ddvec::dd_vector::const_iterator v0_ptr = vector_0_fill->begin( );
  Ddvec::dd_vector::const_iterator v1_ptr = vector_1_fill->begin( );
  bool Debug = false;
  if( Debug )
  {
    for( ; ( v0_ptr != vector_0_fill->end( ) ) && ( v0_ptr != vector_0_fill->end( ) );
         ++v0_ptr, ++v1_ptr )
    {
      std::cout << v0_ptr->x << "  " << v0_ptr->y << "  " << v1_ptr->x << "  " << v1_ptr->y << std::endl;
    }
    std::cout << std::endl;
  }
}
