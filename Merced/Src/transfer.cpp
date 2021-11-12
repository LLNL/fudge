/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: transfer.cpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
 */
// Implementation of the transfer matrix class

#include <cstdlib>
#include <iostream>

#include "transfer.hpp"
#include "messaging.hpp"
#include "global_params.hpp"

//****************** class Trf::T_matrix  *****************
//--------------- Trf::T_matrix constructor -----------------
// This makes the code use the lowest cross section energy as the threshold.
Trf::T_matrix::T_matrix( ): order( -1 ), conserve( Coef::BOTH ), threshold( 0.0 )
{
}
//--------------- Trf::T_matrix destructor -----------------
Trf::T_matrix::~T_matrix( )
{
  if( ( conserve != Coef::NOT_SET ) && ( order >= 0 ) )
  {
    delete [] data;
    delete [] row_sums;
    delete [] row_checks;
  }
}
//--------------- Trf::T_matrix::allocate -----------------
void Trf::T_matrix::allocate( )
{
  order = Global.Value( "outputLegendreOrder" );
  if( order > Order::max_output_order )
  {
    Msg::FatalError( "Trf::T_matrix::allocate",
                "increase the size of Order::max_output_order and recompile" );
  }
  e_flux.order = order;

  if( conserve == Coef::NOT_SET )
  {
    Msg::FatalError( "Trf::T_matrix::allocate", 
		 "Set the conservation flag before calling this routine" );
  }
  if( ( num_Ein_bins <= 0 ) || ( num_Eout_bins <= 0 ) )
  {
    Msg::FatalError( "Trf::T_matrix::allocate", 
		 "Set the energy bins before calling this routine" );
  }
  data = new Coef::coef_vector[ num_Ein_bins * num_Eout_bins ];
  // initialize to zero
  for( int i = 0; i < num_Ein_bins * num_Eout_bins; ++i ){
    data[i].set_order( order, conserve );
  }

  // for checking the row sums
  row_checks = new Coef::coef_vector[ num_Ein_bins ];
  row_sums = new Coef::coef_vector[ num_Ein_bins ];
  // initialize to zero
  for( int i = 0; i < num_Ein_bins; ++i )
  {
    row_sums[i].set_order( 0, conserve );
    row_checks[i].set_order( 0, conserve );
  }
}
//--------------- Trf::T_matrix::operator( ) -----------------
Coef::coef_vector& Trf::T_matrix::operator()( int Ein_count, int Eout_count )
{
  int index = Eout_count + num_Eout_bins * Ein_count;
  return data[ index ];
}
//------------ Trf::T_matrix::get_flux_weight -------------------
// Calculates the weights 1 / \int_(energy group) flux
void Trf::T_matrix::get_flux_weight( )
{
  // start the list
  Lgdata::weight_list::iterator next_weight = flux_weight.insert( flux_weight.end( ), Lgdata::weight_vector( ) );
  next_weight->initialize( order );

  // the current trapezoid to add
  double E_left;
  double E_right;
  // pointers to the incident energy bin edges
  Egp::Energy_groups::const_iterator this_Ein = in_groups.begin( );
  Egp::Energy_groups::const_iterator next_Ein = this_Ein;
  ++next_Ein;
  // pointers to the flux data points
  Lgdata::Flux_List::const_iterator this_flux = e_flux.begin( );
  Lgdata::Flux_List::const_iterator next_flux = this_flux;
  ++next_flux;

  // do through the energy bins
  for( ; ; )
  {
    if( next_flux->get_E_in( ) <= *this_Ein )
    {
      // increment the flux data
      this_flux = next_flux;
      ++next_flux;
      if( next_flux == e_flux.end( ) )
      {
        Msg::FatalError( "transfer::get_flux_weight",
                     "flux data ended early" );
      }
    }
    else
    {
      // do the integral for this flux data and energy bin
      E_left = ( this_flux->get_E_in( ) < *this_Ein ) ? *this_Ein :
        this_flux->get_E_in( );
      E_right = ( next_flux->get_E_in( ) > *next_Ein ) ? *next_Ein :
        next_flux->get_E_in( );
      next_weight->increment( e_flux, this_flux, E_left, E_right );
      if( next_flux->get_E_in( ) < *next_Ein )
      {
        // increment the flux data
        this_flux = next_flux;
        ++next_flux;
        if( next_flux == e_flux.end( ) )
        {
          Msg::FatalError( "transfer::get_flux_weight",
                     "Flux data ended early" );
        }
      }
      else
      {
        // finished with this energy bin
        next_weight->invert( );
        // go to the next energy bin
        this_Ein = next_Ein;
        ++next_Ein;
        if( next_Ein == in_groups.end( ) )
        {
          break;  // we are finished
        }
        next_weight = flux_weight.insert( flux_weight.end( ), Lgdata::weight_vector( ) );
        next_weight->initialize( order );
      }
    }
  }
}
//------------ Trf::T_matrix::use_weight -------------------
// Applies the weights 1 / \int_(energy group) flux
void Trf::T_matrix::use_weight( )
{
  Lgdata::weight_list::iterator this_weight = flux_weight.begin( );
  for( int Ein_count = 0; Ein_count < num_Ein_bins; ++Ein_count,
         ++this_weight )
  {
    for( int Eout_count = 0; Eout_count < num_Eout_bins; ++Eout_count )
    {
      ( *this )( Ein_count, Eout_count ) *= *this_weight;
    }
  }
}
//------------ Trf::T_matrix::scale_E -------------------
// Scales the weight_E terms by the average E_out, to check gamma output
void Trf::T_matrix::scale_E( )
{
  Egp::Energy_groups::const_iterator Eout_ptr = out_groups.begin( );
  Egp::Energy_groups::const_iterator next_Eout = Eout_ptr;
  ++next_Eout;
  for( int Eout_count = 0; Eout_count < num_Eout_bins;
       Eout_ptr = next_Eout, ++next_Eout, ++Eout_count )
  {
    double scale_by = 2.0 / ( *Eout_ptr + *next_Eout );   
    for( int Ein_count = 0; Ein_count < num_Ein_bins; ++Ein_count )
    {
      ( *this )( Ein_count, Eout_count ).scale_E( scale_by );
    }
  }
}
//------------ Trf::T_matrix::check_ell0 -------------------
// Prints the sums of the rows of the zero-order term.
// They should agree with the flux-weighted average cross sections
void Trf::T_matrix::check_ell0( )
{
  if( ( conserve == Coef::NUMBER ) || ( conserve == Coef::BOTH ) )
  {
    for( int Ein_count = 0; Ein_count < num_Ein_bins;
         ++Ein_count )
    {
      row_sums[ Ein_count ].weight_1[ 0 ] = 0.0;
      for( int Eout_count = 0; Eout_count < num_Eout_bins;
           ++Eout_count )
      {
	row_sums[ Ein_count ].weight_1[ 0 ] +=
          (*this)( Ein_count, Eout_count ).weight_1[ 0 ];
      }
    }
  }
  if( ( conserve == Coef::ENERGY ) || ( conserve == Coef::BOTH ) )
  {
    for( int Ein_count = 0; Ein_count < num_Ein_bins; 
         ++Ein_count )
    {
      row_sums[ Ein_count ].weight_E[ 0 ] = 0.0;
      for( int Eout_count = 0; Eout_count < num_Eout_bins;
           ++Eout_count )
      {
        row_sums[ Ein_count ].weight_E[ 0 ] += 
          (*this)( Ein_count, Eout_count ).weight_E[ 0 ];
      }
    }
  }
  static int data_precision = Global.Value( "datafield_precision" );
  static int field_width = Global.get_field_width( );
  static int check_row_sum = Global.Value( "check_row_sum" );
  if( ( conserve == Coef::NUMBER ) || ( conserve == Coef::BOTH ) )
  {
    std::cout.setf( std::ios::scientific, std::ios::floatfield);
    if( check_row_sum > 0 )
    {
      std::cout << "Row sums for number\n";
      std::cout << "    integral       row sum       difference       relative\n";
    }
    for( int Ein_count = 0; Ein_count < num_Ein_bins;
         ++Ein_count )
    {
      double row_sum = row_sums[ Ein_count ].weight_1[ 0 ];
      double check_sum = row_checks[ Ein_count ].weight_1[ 0 ];
      double diff = row_sum - check_sum;
      double relative;
      if( check_sum == 0.0 )
      {
	relative = ( row_sum == 0.0 ) ? 0.0 : 1.0;
      }
      else
      {
        relative = diff/check_sum;
      }
      if( check_row_sum > 0 )
      {
        std::cout << std::setw(field_width) <<
            std::setprecision(data_precision) <<
            check_sum << " " <<
          std::setw(field_width) <<
            std::setprecision(data_precision) <<
            row_sum << " " <<
          std::setw(field_width) <<
            diff << " " <<
          std::setw(field_width) <<
            relative << std::endl;
      }
      static int scale_rows = Global.Value( "scale_rows" );
      if( ( scale_rows > 0 ) && ( row_sum > 0.0 ) )
      {
        for( int Eout_count = 0; Eout_count < num_Eout_bins;
             ++Eout_count )
        {
          (*this)( Ein_count, Eout_count ) *= check_sum/row_sum;
	}
      }
    }
  }
/*
 *  if( ( conserve == Coef::ENERGY ) || ( conserve == Coef::BOTH ) )
 *  {
 *    std::cout.setf( std::ios::scientific, std::ios::floatfield);
 *    std::cout << "Row sums for energy\n";
 *    std::cout << "    integral       row sum       difference       relative\n";
 *    for( int Ein_count = 0; Ein_count < num_Ein_bins;
 *         ++Ein_count )
 *    {
 *      double row_sum = row_sums[ Ein_count ].weight_E[ 0 ];
 *      // double check_sum = row_checks[ Ein_count ].weight_E[ 0 ];
 *      double check_sum = row_checks[ Ein_count ].weight_1[ 0 ];
 *      double diff = row_sum - check_sum;
 *      double relative = ( check_sum == 0.0 ) ? 1.0 : diff/check_sum;
 *      std::cout << std::setw(field_width) <<
 *            std::setprecision(data_precision) <<
 *            check_sum << " " <<
 *          std::setw(field_width) <<
 *            std::setprecision(data_precision) <<
 *            row_sum << " " <<
 *          std::setw(field_width) <<
 *            diff << " " <<
 *          std::setw(field_width) <<
 *            relative << " "<< std::endl;
 *    }
 *  }
*/
}
//------------ Trf::T_matrix::scale_row_check -------------------
// Scales the row check sums by the weights
void Trf::T_matrix::scale_row_check( )
{
  Lgdata::weight_list::iterator this_weight = flux_weight.begin( );
  for( int Ein_count = 0; Ein_count < num_Ein_bins; ++Ein_count,
         ++this_weight )
  {
    row_checks[ Ein_count ] *= *this_weight;
  }
}
//------------ Trf::T_matrix::getBinCrossSection -------------------
// Computes the average cross section for each energy bin
void Trf::T_matrix::getBinCrossSection( const Ddvec::dd_vector& sigma )
{
  Egp::Energy_groups::const_iterator nextEinBinPtr = in_groups.begin( );
  double prevEinBin = *nextEinBinPtr;
  ++nextEinBinPtr;
  double nextEinBin = *nextEinBinPtr;
  Ddvec::dd_vector::const_iterator sigmaPtr = sigma.begin( );
  double prevEin = sigmaPtr->x;
  double prevSigma = sigmaPtr->y;
  ++sigmaPtr;
  double nextEin = sigmaPtr->x;
  double nextSigma = sigmaPtr->y;
  double sum = 0.0;
  maximumCrossSection = 0.0;
  double averageSigma;
  for( ;; )
  {
    if( nextEinBin <= prevEin )
    {
      averageCrossSections.push_back( 0.0 );  // We're below threshold
      prevEinBin = nextEinBin;
      ++nextEinBinPtr;
      if( nextEinBinPtr == in_groups.end( ) )
      {
	zero_transfer( );  // all bins are below the threshold
      }
      nextEinBin = *nextEinBinPtr;
    }
    else if( nextEin <= prevEinBin )
    {
      prevSigma = nextSigma;  // we have cross sections below the bottom bin
      prevEin = nextEin;
      ++sigmaPtr;
      if( sigmaPtr == sigma.end( ) )
      {
	break;
      }
      nextEin = sigmaPtr->x;
      nextSigma = sigmaPtr->y;
    }
    else if( nextEin <= nextEinBin )
    {
      sum += ( nextSigma + prevSigma ) * ( nextEin - prevEin );
      prevSigma = nextSigma;
      prevEin = nextEin;
      ++sigmaPtr;
      if( sigmaPtr == sigma.end( ) )
      {
	averageSigma = sum/( 2*( nextEinBin - prevEinBin ) );
	averageCrossSections.push_back( averageSigma );
	if( averageSigma > maximumCrossSection )
	{
	  maximumCrossSection = averageSigma;
	}
	break;
      }
      if( nextEin == nextEinBin )
      {
	prevEinBin = nextEinBin;
        ++nextEinBinPtr;
        if( nextEinBinPtr == in_groups.end( ) )
        {
	  break;
        }
        nextEinBin = *nextEinBinPtr;
        sum = 0.0;
      }
      nextEin = sigmaPtr->x;
      nextSigma = sigmaPtr->y;
    }
    else if( nextEin > nextEinBin )
    {
      double alpha = ( nextEinBin - prevEin ) / ( nextEin - prevEin ); // interpolate
      double midSigma = (1.0 - alpha)*prevSigma + alpha*nextSigma;
      sum += ( midSigma + prevSigma ) * ( nextEinBin - prevEin );
      averageSigma = sum/( 2*( nextEinBin - prevEinBin ) );
      averageCrossSections.push_back( averageSigma );
      if( averageSigma > maximumCrossSection )
      {
	maximumCrossSection = averageSigma;
      }
      sum = 0.0;
      prevSigma = midSigma;
      prevEin = nextEinBin;
      prevEinBin = nextEinBin;
      ++nextEinBinPtr;
      if( nextEinBinPtr == in_groups.end( ) )
      {
	break;
      }
      nextEinBin = *nextEinBinPtr;
    }
  }
}
//------------ Trf::T_matrix::write_transfer -------------------
// Prints the matrix to the output file
void Trf::T_matrix::write_transfer( )
{
  *output_file << "outputLegendreOrder: " << order << std::endl;
  *output_file << "numEinBins: " << num_Ein_bins << std::endl;
  *output_file << "numEoutBins: " << num_Eout_bins << std::endl;

  static int data_precision = Global.Value( "datafield_precision" );
  static int field_width = Global.get_field_width( );
  Coef::coef_vector this_entry( order, conserve );
  // print the weight 1 integrals
  if( ( conserve == Coef::NUMBER ) || ( conserve == Coef::BOTH ) )
  {
    *output_file << "Integrals, weight = 1: numEinBins = " <<
      num_Ein_bins << std::endl;
    output_file->setf( std::ios::scientific, std::ios::floatfield);
    for( int Ein_count = 0; Ein_count < num_Ein_bins;
         ++Ein_count )
    {
      *output_file << "EinBin = " << Ein_count << " : numEoutBins = " <<
        num_Eout_bins << std::endl;
      for( int Eout_count = 0; Eout_count < num_Eout_bins;
         ++Eout_count )
      {
        this_entry = ( *this )( Ein_count, Eout_count );
        for( int t_count = 0; t_count <= order; ++t_count )
        {
          *output_file << std::setw(field_width) <<
            std::setprecision(data_precision) <<
            this_entry.weight_1[t_count] << " ";
        }
        *output_file << std::endl;
      }
    }
  }

  // print the weight E' integrals
  if( ( conserve == Coef::ENERGY ) || ( conserve == Coef::BOTH ) )
  {
    *output_file << "Integrals, weight = E': numEinBins = " <<
      num_Ein_bins << std::endl;
    output_file->setf( std::ios::scientific, std::ios::floatfield);
    for( int Ein_count = 0; Ein_count < num_Ein_bins;
         ++Ein_count )
    {
      *output_file << "EinBin = " << Ein_count << " : numEoutBins = " <<
        num_Eout_bins << std::endl;
      for( int Eout_count = 0; Eout_count < num_Eout_bins;
         ++Eout_count )
      {
        this_entry = ( *this )( Ein_count, Eout_count );
        for( int t_count = 0; t_count <= order; ++t_count )
        {
          *output_file << std::setw(field_width) <<
            std::setprecision(data_precision) <<
            this_entry.weight_E[t_count] <<  " ";
        }
        *output_file << std::endl;
      }
    }
  }
}
//------------ Trf::T_matrix::zero_transfer -------------------
// Prints zeros to the output file for a reaction with high threshold
void Trf::T_matrix::zero_transfer( )
{
  Msg::Warning( "Trf::T_matrix::zero_transfer",
		"omit reaction: all energy data too high" );
  *output_file << "outputLegendreOrder: " << order << std::endl;
  *output_file << "numEinBins: " << num_Ein_bins << std::endl;
  *output_file << "numEoutBins: " << num_Eout_bins << std::endl;

  static int data_precision = Global.Value( "datafield_precision" );
  static int field_width = Global.get_field_width( );
  // print the weight 1 integrals
  if( ( conserve == Coef::NUMBER ) || ( conserve == Coef::BOTH ) )
  {
    *output_file << "Integrals, weight = 1: numEinBins = " <<
      num_Ein_bins << std::endl;
    output_file->setf( std::ios::scientific, std::ios::floatfield);
    for( int Ein_count = 0; Ein_count < num_Ein_bins;
         ++Ein_count )
    {
      *output_file << "EinBin = " << Ein_count << " : numEoutBins = " <<
        num_Eout_bins << std::endl;
      for( int Eout_count = 0; Eout_count < num_Eout_bins;
         ++Eout_count )
      {
        for( int t_count = 0; t_count <= order; ++t_count )
        {
          *output_file << std::setw(field_width) <<
            std::setprecision(data_precision) <<
            0.0 << " ";
        }
        *output_file << std::endl;
      }
    }
  }

  // print the weight E' integrals
  if( ( conserve == Coef::ENERGY ) || ( conserve == Coef::BOTH ) )
  {
    *output_file << "Integrals, weight = E': numEinBins = " <<
      num_Ein_bins << std::endl;
    output_file->setf( std::ios::scientific, std::ios::floatfield);
    for( int Ein_count = 0; Ein_count < num_Ein_bins;
         ++Ein_count )
    {
      *output_file << "EinBin = " << Ein_count << " : numEoutBins = " <<
        num_Eout_bins << std::endl;
      for( int Eout_count = 0; Eout_count < num_Eout_bins;
         ++Eout_count )
      {
        for( int t_count = 0; t_count <= order; ++t_count )
        {
          *output_file << std::setw(field_width) <<
            std::setprecision(data_precision) <<
            0.0 <<  " ";
        }
        *output_file << std::endl;
      }
    }
  }
  // *** Bail out ***
  exit( 0 );
}
