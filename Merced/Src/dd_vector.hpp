/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: dd_vector.hpp 1 2006-02-02 03:06:56Z hedstrom $
 *
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// header for the dd_vector class
#ifndef DD_VECTOR_CLASS
#define DD_VECTOR_CLASS

#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <list>
#include <utility>  // for the pair class

#include "data_parser.hpp"

namespace Terp
{
//! Specifies the interpolation rule
// ------------------------ Interp_Type ---------------
enum Interp_Type{ NOTSET, LINLIN, HISTOGRAM, LINLOG, LOGLIN, LOGLOG };

// ------------------------ Interp_qualifier ---------------
enum Interp_qualifier{ UNSET, DIRECT, UNITBASE, CUMULATIVE_POINTS,
   UNSCALED_DIRECT, UNSCALED_UNITBASE, UNSCALED_CUMULATIVE_POINTS };

// ------------------------ class two_d_interp ----------------------
//! class for 2-dimensional interpolation types
class two_d_interp
{
public:
  Interp_Type flag;  // LINLIN, etc.
  Interp_qualifier qualifier; // UNITBASE, etc.

  two_d_interp( ) {}
  ~two_d_interp( ) {}

  //! to copy
  //! \param to_copy: the data to copy
  two_d_interp& operator=( const two_d_interp &to_copy );
};
} // end of namespace Terp

// ***************** namespace interp_flag_F *************************
//! Routines to read interpolation flags
namespace interp_flag_F
{
  // ------------------------ read_flag ---------------
  //! Interprets the interpolation rule
  //! Returns  HISTOGRAM, LINLIN, etc.
  //! \param interp_ID: "flat", "lin-lin", etc.
  //! \param input_file: to handle fatal errors
  Terp::Interp_Type read_flag( const std::string& interp_ID, Dpar::data_parser &input_file );

  // ------------------------ read_qualifier ---------------
  //! Interprets the interpolation qualifier
  //! Returns  DIRECT, UNITBASE, etc.
  //! \param interp_ID: "direct", "unitbase", etc.
  //! \param input_file: to handle fatal errors
  Terp::Interp_qualifier read_qualifier( const std::string& qualifier_ID,
    Dpar::data_parser &input_file );

  // ------------------------ read_1d_interpolation ---------------
  //! Interprets the interpolation rule and outputs it
  //! \param input_file the input file with data P( Eout | Ein )
  Terp::Interp_Type read_1d_interpolation( Dpar::data_parser &input_file );

  // ------------------------ read_2d_interpolation ---------------
  //! Interprets the interpolation rules for 2d data
  //! \param input_file the input file with data P( Eout | Ein )
  //! \param energy_in_interp interpolation rules with respect to Ein
  //! \param data_out_qualifier interpolation qualifier with respect to Ein: unitbase
  void read_2d_interpolation( Dpar::data_parser &input_file, 
    Terp::two_d_interp *energy_in_interp, Terp::Interp_Type *data_out_interp );

  // ------------------------ read_3d_interpolation_ENDL ---------------
  //! Interprets the interpolation rules for ENDL 3d data
  //! \param input_file the input file with data P( Eout | Ein, cosine )
  //! \param energy_in_interp interpolation with respect to Ein
  //! \param cosine_interp interpolation with respect to cosine
  //! \param energy_out_interp interpolation with respect to Eout
  void read_3d_interpolation_ENDL( Dpar::data_parser &input_file,
    Terp::two_d_interp *energy_in_interp, Terp::two_d_interp *cosine_interp,
    Terp::Interp_Type *energy_out_interp );

  // ------------------------ read_3d_interpolation_GND ---------------
  //! Interprets the interpolation rules for GND 3d data
  //! \param input_file the input file with data P( Eout | Ein, cosine )
  //! \param energy_in_interp interpolation with respect to Ein
  //! \param cosine_interp interpolation with respect to cosine
  //! \param energy_out_interp interpolation with respect to Eout
  void read_3d_interpolation_GND( Dpar::data_parser &input_file,
    Terp::two_d_interp *energy_in_interp, Terp::two_d_interp *cosine_interp,
    Terp::Interp_Type *energy_out_interp );
}
// ****************** end of namespace interp_flag_F *****************

namespace Ddvec
{
// ------------------------ class Ein_entry ---------------
//! A class for one 
class Ein_entry
{
public:
  double Ein;  // the incident energy
  int n;    // the repetition count

  //!inline constructor
  inline Ein_entry( ) {}

  //!inline destructor
  inline ~Ein_entry( ) {}
};

// ------------------------ class Ein_vector ---------------
//! A class used to hold the structure of a dd_vector
class Ein_vector : public std::vector< Ein_entry >
{
public:
  bool start_zero;  // is the first y-value zero?
  bool end_zero;    // is the last y-value zero?

  //!inline constructor
  inline Ein_vector( ) {}

  //!inline destructor
  inline ~Ein_vector( ) {}

  // for debugging
  void print( ) const;
};

// ------------------------ class Ein_entries ---------------
//! A class used in the interpolation of two dd_vectors to a common set of energies
class Ein_entries
{
public:
  double Ein;       // the incident energy
  int n[ 2 ];    // the repetition counts

  //!inline constructor
  inline Ein_entries( ) {}

  //!inline destructor
  inline ~Ein_entries( ) {}
};

// ------------------------ class Ein_vectors ---------------
//! A class used to manage interpolation of two dd_vectors to a common set of energies
class Ein_vectors : public std::vector< Ein_entries >
{
public:
  Ein_vectors::iterator start_jump[ 2 ];  // keep track of an initial jump
  Ein_vectors::iterator end_jump[ 2 ];    // keep track of a final jump

  //!inline constructor
  inline Ein_vectors( ) {}

  //!inline destructor
  inline ~Ein_vectors( ) {}

  //! Copies the j-th initial Ein_entries
  //! \param j = 0,1, which list to copy to
  //! \param Einj_ptr pointer to the first entry to copy
  //! \param Eink_ptr pointer to the first entry not to copy
  //! \param Einj the list to copy
  //! \param start_zero do we save the pointer to the last entry copied?
  void copy_head( int j, Ein_vector::const_iterator &Einj_ptr,
     Ein_vector::const_iterator &Eink_ptr, const Ein_vector &Einj, bool start_zero );

  //! Copies the j-th final Ein_entries
  //! \param j = 0,1, which list to copy to
  //! \param Einj_ptr pointer to the first entry to copy
  //! \param Einj the list to copy
  //! \param end_zero do we save the pointer to the first entry copied?
  void copy_tail( int j, Ein_vector::const_iterator &Einj_ptr,
		  const Ein_vector &Einj, bool end_zero );

  //! Produces Ein_vectors with a common set of incident energies
  //! \param Ein_1 list of Ein_entries, pairs ( energy, number of values )
  //! \param Ein_2 list of Ein_entries, pairs ( energy, number of values )
  void common_Ein( const Ein_vector &Ein_1, const Ein_vector &Ein_2 );
};

// ------------------------ class dd_entry ---------------
//! The basic entry holds a pair of doubles
class dd_entry
{
public:
  double x;
  double y;

  //!Default inline constructor is zero.
  inline dd_entry( )
  {
    x = 0.0;
    y = 0.0;
  }  
  //!Inline copy constructor
  //! param entry data to copy
  inline dd_entry( const dd_entry& entry )
  {
    x = entry.x;
    y = entry.y;
  }
  //!Inline copy constructor
  //! param xx value to copy
  //! param yy value to copy
  inline dd_entry(double xx, double yy)
  {
    x = xx;
    y = yy;
  }

  //!Inline default desctructor is blank.
  inline ~dd_entry() {}

  //! Makes a copy
  //! param to_copy data to copy
  inline dd_entry& operator=( const dd_entry& to_copy )
  {
    x = to_copy.x;
    y = to_copy.y;
    return *this;
  }

  //! Does linear-linear interpolation
  //! param E intermediate energy
  //! param next_entry next ( E, y ) pair
  //! \param interp_OK, true if the interpolation is OK
  double linlin_interp( double E, const dd_entry& next_entry, bool *interp_OK ) const;

  //! Does linear-log interpolation
  //! param E intermediate energy
  //! param next_entry next ( E, y ) pair
  //! \param interp_OK, true if the interpolation is OK
  double linlog_interp( double E, const dd_entry& next_entry, bool *interp_OK ) const;

  //! A first step in bilinear interpolation
  //! Returns true if the interpolation is OK
  //! param E intermediate incident energy
  //! param left_data_Ein lower incident energy, Ein_left
  //! param left_data pair ( Eout, P( Eout | Ein_left ) )
  //! param right_data_Ein lower incident energy, Ein_right
  //! param right_data pair ( Eout, P( Eout | Ein_right ) )
  bool linlin_interp( double E, double left_data_Ein, const dd_entry& left_data,
    double right_data_Ein, const dd_entry& right_data );

  //! A first step in bilinear-log interpolation
  //! Returns true if the interpolation is OK
  //! param E intermediate incident energy
  //! param left_data_Ein lower incident energy, Ein_left
  //! param left_data pair ( Eout, P( Eout | Ein_left ) )
  //! param right_data_Ein lower incident energy, Ein_right
  //! param right_data pair ( Eout, P( Eout | Ein_right ) )
  bool linlog_interp( double E, double left_data_Ein, const dd_entry& left_data,
    double right_data_Ein, const dd_entry& right_data );

  // for debugging
  void print( ) const;
};

// ------------------------ class dd_pair ---------------
//!The dd_pair class is used in interpolation between 2 data points.
class dd_pair : public std::pair< dd_entry, dd_entry >
{
private:
  double tag;

public:
  //! How to interpolate between the two entries
  Terp::Interp_Type Eout_interp;

  //!Default inline constructor is blank.
  inline dd_pair( ): Eout_interp( Terp::LINLIN ) {}

  //!Inline default desctructor is blank.
  inline ~dd_pair() {}

  //! Stores the data entries
  //! \param first_ first pair ( x, y )
  //! \param second_ second pair ( x, y )
  void set_pair( const dd_entry& first_, const dd_entry& second_ );

  //! Do linear interpolation between first and second
  //! \param eta interpolate y for first.x <= eta <= second.x
  //! \param interp_OK, true if the interpolation is OK
  double value( double eta, bool *interp_OK );

  //! Solves for x given y
  //! \param Eout for linear interpolation find x such that y = Eout
  //! \param interp_OK, true if the interpolation is OK
  double inverse( double Eout, bool *interp_OK );

  //! Solves for x given y
  //! \param Eout for linear interpolation find x such that y = Eout
  //! \param flag flag = -1: failure, flag = 0: no solution, flag = 2: many
  double find_Ein( double Eout, int *flag );

  //! The usual tag is the incident energy.
  inline double get_E_in() const
  {
    return tag;
  }

  //! Sets the tag
  //! \param E_in the tag value, usually the energy of the incident particle
  inline void set_E_in( double E_in )
  {
    tag = E_in;
  }

  //! The tag may be direction cosine
  inline double get_mu() const
  {
    return tag;
  }

  //! Sets the tag
  //! \param E_in the tag value when it is direction cosine
  inline void set_mu( double mu )
  {
    tag = mu;
  }

  //! Does linear-linear interpolation with respect to the tag
  //! Returns true if the interpolation is OK
  //! \param E_in the intermediate incident energy
  //! \param lower_data data at lower incident energy
  //! \param higher_data data at higher incident energy
  bool linlin_interp( double E_in, const dd_pair &lower_data,
		      const dd_pair &higher_data );

  //! Does linear-log interpolation with respect to the tag
  //! Returns true if the interpolation is OK
  //! \param E_in the intermediate incident energy
  //! \param lower_data data at lower incident energy
  //! \param higher_data data at higher incident energy
  bool linlog_interp( double E_in, const dd_pair &lower_data,
		      const dd_pair &higher_data );

  //! Sets up the data for outgoing energies Eout_min and Eout_max
  //! Returns true if the interpolation is OK
  //! \param prev_data data at a lower outgoing energy
  //! \param next_data data at a higher outgoing energy
  //! \param Eout_min the desired lower outgoing energy
  //! \param Eout_max the desired higher outgoing energy
  bool set_data( const dd_entry &prev_data, const dd_entry &next_data,
		 double Eout_min, double Eout_max );
};

//! Class for mapping to unit base and back
// ---------------- class unit_base_map ------------------
class unit_base_map
{
public:
  double Eout_min;  // for the unit-base transformation
  double Eout_max;

  inline unit_base_map( ) { }
  inline ~unit_base_map( ) { }

  //! Undoes the unit-base map---finds E' for given eta
  //! For this version the code does a fatal exit for bad eta.
  //! \param eta unit-base energy, 0 <= eta <= 1
  double un_unit_base( double eta ) const;

  //! Undoes the unit-base map---finds E' for given eta
  //! \param eta unit-base energy, 0 <= eta <= 1
  //! \param interp_OK, true if the interpolation is OK
  double un_unit_base( double eta, bool *interp_OK ) const;

  //! Finds eta for given physical E'
  //! \param physical_E the physical energy
  //! \param interp_OK, true if the interpolation is OK
  double to_unit_base( double physical_E, bool *interp_OK );

  //! Makes a copy
  //! \param to_copy the data to copy
  void copy( const unit_base_map &to_copy );

  //! Interpolates the energy range
  //! Returns true if the interpolation is OK
  //! \param alpha interpolation parameter 0 <= alpha <= 1
  //! \param prev the outgoing energy range for the lower incident energy
  //! \param next the outgoing energy range for the higher incident energy
  bool interpolate( double alpha, const unit_base_map &prev,
    const unit_base_map &next );

  //! Tests for a trivial interval
  bool too_short( );
};

// ----------------------- class cum_points_pair -------------------
//! Class for dd_pairs used in interpolation by cumulative points
class cum_points_pair : public dd_pair
{
public:
  unit_base_map ubase_map;

  cum_points_pair( ) {}
  ~cum_points_pair( ) {}
  
  //! Tests for a trivial interval
  bool too_short( );

  //! Maps the pair to unit base
  void to_unit_base( );
  
  //! Maps the pair to unit base, used for intervals of length zero
  //! \param dA, the probability for this interval
  void short_to_unit_base( double dA );

  //! Undoes the mapping to unit base
  void un_unit_base( );
};

// ----------------------- class dd_vector -------------------
//! Vector of dd_entrys.
class dd_vector : public std::vector< dd_entry >
{
private:
  //! A label for the vector.
  double tag;

public:
  Terp::Interp_Type interp_type;

  //! Default constructor
  dd_vector( ): tag( 0.0 ), interp_type( Terp::NOTSET )
  {}

  //! Copy constructor
  //! \param to_copy the data to copy
  dd_vector( const dd_vector &to_copy );

  //! Default destructor
  ~dd_vector( );

  //! Returns the tag for this vector
  inline double get_tag() const
  {
    return tag;
  }
 
  //!Common name used to access the label for this vector.
  inline double get_E_in() const
  {
    return tag;
  }

  //! Sets the tag
  //! \param E_in the incident energy for this data
  inline void set_E_in( double E_in )
  {
    tag = E_in;
  }

  //! Sets the tag
  //! \param Tag the value of the tag for this data
  inline void set_tag( double Tag )
  {
    tag = Tag;
  }
 
  //! Reads the data from a file
  //! \param input_file the input file to read
  //! \param num_sigma the number of pairs ( x, y ) to read
  void read_data( Dpar::data_parser &input_file, int num_sigma );

  //! Reads the data and interpolation type from a file
  //! \param input_file the input file to read
  //! \param num_sigma the number of pairs ( x, y ) to read
  void read_data_interp( Dpar::data_parser &input_file, int num_sigma );

  //! Appends an ( E_out, Prob ) entry to the vector
  //! \param E_out the x-value of the pair ( x, y ) to append to the vector
  //! \param Prob the y-value of the pair ( x, y ) to append to the vector
  void add_entry( double E_out, double Prob );

  //! Makes a constant vector
  //! \param sigma use the energy range of this vector as the x-values for (x, y)
  //! \param val the y-value for the pairs (x, y)
  void make_flat( const dd_vector& sigma, double val );

  //! Makes a constant vector
  //! \param E_groups use the energy range of this vector as the x-values for (x, y)
  //! \param val the y-value for the pairs (x, y)
  void make_flat( const std::vector< double >& E_groups, double val );

  //! Copies a vector
  //! \param vector_from the vector to copy
  void copy( const dd_vector& vector_from );

  //! Copies a vector with extrapolation
  //! \param vector_from the vector to copy
  //! \param min_x if necessary, extrapolate as zero to this lowest energy
  //! \param max_x if necessary, extrapolate as zero to this highest energy
  void extrapolate_copy( const dd_vector& vector_from, double min_x, double max_x );

  //! Copies a vector with truncation
  //! Returns false if truncation gave a vector with zero norm
  //! \param vector,_from the vector to copy
  //! \param min_x, if necessary, chop to this lowest energy
  //! \param max_x, if necessary, chop to this highest energy
  //! \param do_renorm, renormalize if true
  bool truncate_copy( const dd_vector& vector_from, double min_x, double max_x,
		      bool do_renorm );

  //! Prepends x=0 using histogram extrapolation if necessary
  //! \param to_copy, the data to copy
  //! \param E_insert, make sure that the result has an entry at this energy
  void prepend( const dd_vector& to_copy, double E_insert );

  //! Adds one vector to another
  //! \param vector_2 add the y's of this (x, y) vector
  dd_vector& operator+=( dd_vector& vector_2 );

  //! Scales the Y value of this vector by factor
  //! \param factor multiply the y's in the (x, y) pairs by this factor
  dd_vector& operator*=( double factor );

  //! Product of two vectors
  //! \param vector_2 multiply the y's of this (x, y) vector
  dd_vector& operator*=( dd_vector& vector_2 );

  //! Prints the vector
  void print( ) const;

  //! Parses this dd_vector
  
  //! \param E_in_list computed vector of (x-value, multiplicity) of the (x, y) pairs
  void parse_vector( Ein_vector *E_in_list ) const;

  //! Compute fill_vector containing entries for each energy in E_in_list
  //! \param E_in_list 2 vectors of pairs (energy, multiplicity)
  //! \param j which E_in_list to use (j = 0, 1)
  //! \param fill_vector computed vector of interpolated pairs (x, y) from both lists
  void fill_with( const Ein_vectors &E_in_list, int j, dd_vector *fill_vector ) const;

  // get the probability for this vector.
  double get_norm( );

  //! For probability-density tables set the norm to 1
  //! Returns false if we get zero norm from truncated direct interpolation
  // \param truncated, true for data from truncated direct interpolation
  bool renorm( bool truncated );

  //! Maps a vector to [0, 1] and saves the map
  //! \param ubase_map save the original energy range
  void mapto_01( unit_base_map *ubase_map );

  //! Maps a vector of probability densities to unit base, and save the map
  //! \param Renorm whether or not to renorm to 1
  //! \param ubase_map save the original energy range
  void unit_base( bool Renorm, unit_base_map *ubase_map );

  //! Inverts mapto_01
  //! \param ubase_map the original energy range
  void un_mapto_01( const unit_base_map *ubase_map );

  //! Inverts the unit-base map of a vector of probability densities
  //! \param ubase_map the original energy range
  void un_unit_base( bool Renorm, const unit_base_map *ubase_map );

  //! Does a linear interpolation of 2 vectors
  //! \param E_in the intermediate incident energy
  //! \param prev_vect a vector of (x, y) pairs at a lower incident energy
  //! \param next_vect a vector of (x, y) pairs at a higher incident energy
  void interpolate( double E_in, const dd_vector &prev_vect, const dd_vector &next_vect );

  //! Does a linear interpolation of 2 vectors with common x-values
  //! \param E_in the intermediate incident energy
  //! \param prev_filled a vector of (x, y) pairs at a lower incident energy
  //! \param next_filled a vector of (x, y) pairs at a higher incident energy
  void filled_interpolate( double E_in, const dd_vector &prev_filled,
    const dd_vector &next_filled );

  //! Truncates histogram data at the maximum energy
  //! \param EMax chop off all (x, y) pairs with x > EMax
  void chop_histogram( double EMax );

  //! Scales the x-component; a kludge to convert ENDF eV to MeV
  //! \param eV_to_MeV scale factor for the x-values in the pairs (x, y)
  void scale_E( double eV_to_MeV );

  //! Checks whether an angular probability density is isotropic
  bool isotropic( ) const;

  //! Finds the interval containing the desired X-value
  //! \param X, the desired X-value
  //! \param prev_ptr, last link such that prev_ptr->x <= X
  //! \param next_ptr, first link such that next_ptr->x > X
  void locate_x( double X, dd_vector::const_iterator *prev_ptr,
		 dd_vector::const_iterator *next_ptr ) const;
};

} // end of namespace Ddvec

namespace dd_vector_F
{
  // This routine is no longer used.
  
  // -------------------- fill_in_vectors ---------------
  //! This routine computes vector_1_fill and vector_2_fill with common x-values.
  //! \param vector_1 first vector of (x, y) values
  //! \param vector_2 second vector of (x, y) values
  //! \param vector_1_fill first vector of (x, y) values, filled to common x-values
  //! \param vector_2_fill second vector of (x, y) values, filled to common x-values
  //! \param re_norm Do we renormalize?
  void fill_in_vectors( const Ddvec::dd_vector& vector_1, const Ddvec::dd_vector& vector_2,
    Ddvec::dd_vector *vector_1_fill, Ddvec::dd_vector *vector_2_fill, bool re_norm );
}

#endif
