/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2008-04-16 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: reaction.cpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/
// Implement the classes used for identification of the reaction data

#include <iostream>
#include <string>

#include "reaction.hpp"
#include "messaging.hpp"
#include "global_params.hpp"
// #include "logger.hpp"

//****************** class num_check_param *****************
// --------------- num_check_param::value -----------------
// Computes the product sigma*multiplicity*flux*weight
double num_check_param::value( double E_in )
{
  if( E_in < model_weight.first.x )
  {
    return 0.0;
  }
  double vaerde = sigma.value( E_in ) *
    mult.value( E_in ) * e_flux.value( E_in ) * model_weight.value( E_in );

  return vaerde;
}

//****************** class reaction  *****************
// --------------- reaction::read_quadrature -----------------
// Interprets the quadrature rule
Quadrature_Method reaction::read_quadrature( data_parser &input_file )
{
  string quadrature = input_file.get_text( );
  string_F::Tolower( quadrature );
  Quadrature_Method use_quad = ADAPTIVE2;
  if( ( quadrature == "adaptive" ) || ( quadrature == "adaptive2" ) )
  {
    use_quad = ADAPTIVE2;
  }
  else if( quadrature == "adaptive4" )
  {
    use_quad = ADAPTIVE4;
  }
  else if( quadrature == "gauss2" )
  {
    use_quad = GAUSS2;  // Gauss 2-nd order
  }
  else if( quadrature == "gauss4" )
  {
    use_quad = GAUSS4; // Gauss 4-th order
  }
  else if( quadrature == "gauss6" )
  {
    use_quad = GAUSS6; // Gauss 6-th order
  }
  else if( quadrature == "square root" )
  {
    use_quad = ADAPT_HALF; // adaptive Gauss with sqrt(1-x) singularity
  }
  else if( quadrature == "exact" )
  {
    Warning( "reaction::read_quadrature", 
      "Using adaptive---the Omega formula for exact quadrature is disabled." );
    // If you really want this, use svn version 278 of the code.
    use_quad = ADAPTIVE2;
  }
  else
  {
    FatalError( "reaction::read_quadrature", "quadrature rule " +
		quadrature + " undefined." );
  }
  return use_quad;
}
// --------------- reaction::read_Global -----------------
// Reads Global parameters from the input file
bool reaction::read_Global( const string &dataID, data_parser &input_file )
{
  string testCopy = dataID;
  string_F::Tolower( testCopy );
  ss_list::iterator SSlist_ptr;
  string pvalue;

  bool found_it = Global.find( testCopy, SSlist_ptr );
  if( found_it )
  {
    pvalue = input_file.get_text( );
    if( !SSlist_ptr->set_in_command )
    {
      string_F::Tolower( pvalue );
      SSlist_ptr->set_value( pvalue );
    }
  }
  else  // some special flags
  {
    if( testCopy == "inversewavelengthtoenergyfactor" )
    {
      pvalue = input_file.get_text( );
      Global.set( "x_to_energy", pvalue );
      found_it = true;
    }
    else if( testCopy == "electron mass" )
    {
      pvalue = input_file.get_text( );
      Global.set( "m_electron", pvalue );
      found_it = true;
    }
    else if( testCopy == "thompsonscattering" )
    {
      pvalue = input_file.get_text( );
      Global.set( "thompson", pvalue );
      found_it = true; 
    }
  }

  return found_it;
}
// --------------- reaction::process_data -----------------
// process data
void reaction::process_data( data_parser &input_file, ofstream *output_file )
{
  // the default quadrature methods
  transfer.Ein_quad_method = ADAPTIVE2;  // for incident energy
  transfer.Eout_quad_method = ADAPTIVE2;  // for outgoing energy
  transfer.mu_quad_method = ADAPTIVE2;  // for outgoing cosine

  transfer.output_file = output_file;
  string File_type = input_file.get_dataID( );
  // Test the file identity and the version
  if( File_type == "DONE" )
  {
    FatalError( "reaction::process_data", "The input file is empty" );
  }
  else if( ( File_type != "xndfgenTransferMatrix" ) &&
      ( File_type != "xendlTransferMatrix" ) )
  {
    FatalError( "reaction::process_data", 
		pastenum( "line ", input_file.line_count ) +
               ": File type: " + File_type + " not implemented" );
  }
  version = input_file.get_next_double(  );
  if( version != 1.0 )
  {
    FatalError( "reaction::process_data",
      pastenum( "Implement version ", version ) );
  }
  *output_file << "merced: version " << version << endl;

  // Get the type of data
  string processID = input_file.get_dataID( );
  string_F::Tolower( processID );
  while( processID == "comment" )
  {
    input_file.get_comment( *(transfer.output_file) );
    processID = input_file.get_dataID( );
    string_F::Tolower( processID );
  }

  string process_orig = input_file.get_text( );
  string process = process_orig;
  string_F::Tolower( process );
  if( process == "two body transfer matrix" )
  {
    two_body( input_file );
  }
  else if( process == "legendre two body transfer matrix" )
  {
    Legendre_two_body( input_file );
  }
  else if( process == "isotropic table" )
  {
    do_isotropic( input_file );
  }
  else if( process == "legendre eepp data transfer matrix" )
  {
    Legendre( input_file );
  }
  else if( process == "legendre energy-angle data" )
  {
    do_ENDFLegendre( input_file );
  }
  else if( process == "double differential emuepp data transfer matrix" )
  {
    joint_dist_table.ENDL_data = true;
    do_joint_dist( input_file );
  }
  else if( process == "endf double differential emuepp data" )
  {
    joint_dist_table.ENDL_data = false;
    do_joint_dist( input_file );
  }
  else if( process == "uncorrelated energy-angle data transfer matrix" )
  {
    do_uncorr( input_file );
  }
  else if( process == "compton scattering" )
  {
    do_Compton( input_file );
  }
  else if( process == "coherent scattering" )
  {
    do_coherent( input_file );
  }
  else if( process == "evaporation spectrum" )
  {
    do_evaporation( input_file );
  }
  else if( process == "maxwell spectrum" )
  {
    do_Maxwell( input_file );
  }
  else if( process == "watt spectrum" )
  {
    do_Watt( input_file );
  }
  else if( process == "madland-nix spectrum" )
  {
    do_MadlandNix( input_file );
  }
  else if( process == "kalbach spectrum" )
  {
    do_Kalbach( input_file );
  }
  else if( process == "phase space spectrum" )
  {
    do_phase_space( input_file );
  }
  else if( process == "general evaporation" )
  {
    do_gen_evap( input_file );
  }
  else
  {
    FatalError( "reaction::process_data", 
		pastenum( "line ", input_file.line_count ) +
             ": process " + process_orig + " not implemented" );
  }

  if( transfer.order >= 0 )
  {
    // Print the input parameters
    Global.print();
    write_transfer( );
  }
}
// --------------- reaction::common_input -----------------
// Reads and processes input data common to all reactions
bool reaction::common_input( const string &dataID, data_parser &input_file )
{
  bool found_it = false;
  // what type of data is it?
  string tempID = dataID;
  string_F::Tolower( tempID );
  if( tempID == "comment" )
  {
    input_file.get_comment( *(transfer.output_file) );
    found_it = true;
  }
  else if( tempID == "projectile's group boundaries" )
  {
    int num_bd = input_file.get_next_int( );
    transfer.num_Ein_bins = num_bd - 1;
    // read the bin boundaries
    transfer.in_groups.read_bd( input_file, num_bd );
    found_it = true;
  }
  else if( tempID == "product's group boundaries" )
  {
    int num_bd = input_file.get_next_int( );
    transfer.num_Eout_bins = num_bd - 1;
    // read the bin boundaries
    transfer.out_groups.read_bd( input_file, num_bd );
    //! Put all low-energy output in the bottom bin
    transfer.out_groups[ 0 ] = 0.0;
    //! Put all high-energy output in the top bin
    transfer.out_groups[ transfer.num_Eout_bins ] *= 5.0;
    found_it = true;
  }
  else if( tempID == "fluxes" )
  {
    int num_flux = input_file.get_next_int( );
    // read the approximate flux
    transfer.e_flux.read_flux( input_file, num_flux );
    if( transfer.e_flux.interp != LINLIN )
    {
      FatalError( "reaction::common_input", 
		  "only linlear-linear flux is implemented" );
    }
    found_it = true;
  }
  else if( tempID == "cross section" )
  {
    int num_sigma = input_file.get_next_int( );
    // read the cross section
    cross_section.read_data_interp( input_file, num_sigma );
    if( cross_section.interp_type != LINLIN )
    {
      FatalError( "reaction::common_input", 
		  "only linlear-linear cross section is implemented" );
    }
    found_it = true;
  }
  else if( tempID == "multiplicity" )
  {
    int num_mult = input_file.get_next_int( );
    // read the multiplicities
    multiple.read_data_interp( input_file, num_mult );
    if( multiple.interp_type != LINLIN )
    {
      FatalError( "reaction::common_input", 
		  "only linlear-linear multiplicity is implemented" );
    }
    found_it = true;
  }
  else if( tempID == "weight" )
  {
    int num_weight = input_file.get_next_int( );
    // read the weights
    model_weight.read_data_interp( input_file, num_weight );
    if( model_weight.interp_type != HISTOGRAM )
    {
      FatalError( "reaction::common_input", 
		  "only histogram weight is implemented" );
    }
    found_it = true;
  }
  else if( tempID == "projectile's mass" )
  {
    particle_info.mProj = input_file.get_next_double( );
    found_it = true;
  }
  else if( tempID == "target's mass" )
  {
    particle_info.mTarg = input_file.get_next_double( );
    found_it = true;
  }
  else if( tempID == "product's mass" )
  {
    particle_info.mProd = input_file.get_next_double( );
    found_it = true;
  }
  else if( tempID == "residual's mass" )
  {
    particle_info.mRes = input_file.get_next_double( );
    found_it = true;
  }
  else if( tempID == "projectile frame" )
  {
    string frame = input_file.get_text( );
    string_F::Tolower( frame );
    if( frame == "lab" )
    {
      frame_in = LAB;
    }
    else
    {
      FatalError( "reaction::common_input", "invalid incident frame: " + frame );
    }
    found_it = true;
  }
  else if( tempID == "product frame" )
  {
    string frame = input_file.get_text( );
    string_F::Tolower( frame );
    if( frame == "lab" )
    {
      frame_out = LAB;
    }
    else if( frame == "centerofmass" )
    {
      frame_out = CM;
    }
    else
    {
      FatalError( "reaction::common_input", "invalid outgoing frame: " + frame );
    }
    found_it = true;
  }
  else if( tempID == "conserve" )
  {
    string cons = input_file.get_text( );
    string_F::Tolower( cons );
    if( cons == "number" )
    {
      transfer.conserve = NUMBER;
    }
    else if( cons == "energy" )
    {
      transfer.conserve = ENERGY;
    }
    else if( cons == "both" )
    {
      transfer.conserve = BOTH;
    }
    else
    {
      Warning( "reaction::common_input", "invalid conservation type: " + cons );
    }
    found_it = true;
  }
  else if( tempID == "quadrature method" )
  {
    // all quadratures use the same method
    transfer.Ein_quad_method = read_quadrature( input_file );
    transfer.Eout_quad_method = transfer.Ein_quad_method;
    transfer.mu_quad_method = transfer.Ein_quad_method;
    found_it = true;
  }
  else if( tempID == "ein quadrature method" )
  {
    transfer.Ein_quad_method = read_quadrature( input_file );
    found_it = true;
  }
  else if( tempID == "eout quadrature method" )
  {
    transfer.Eout_quad_method = read_quadrature( input_file );
    found_it = true;
  }
  else if( tempID == "mu quadrature method" )
  {
    transfer.mu_quad_method = read_quadrature( input_file );
    found_it = true;
  }
  else if( tempID == "interpolate eout integrals" )
  {
    string interp = input_file.get_text( );
    string_F::Tolower( interp );
    if( interp == "true" )
    {
      transfer.interpolate_Eout_integrals = true;
    }
    else
    {
      transfer.interpolate_Eout_integrals = false;
    }
    found_it = true;
  }
  else if( read_Global( tempID, input_file ) )
  {
    found_it = true;
  }
  else if( dataID == "Bad dataID" )
  {
    FatalError( "reaction::common_input", "Expected dataID in " +
		pastenum( "line ", input_file.line_count ) +
		input_file.original_line );
  }

  return found_it;
}
// --------------- reaction::two_body -----------------
// Reads and processes input data for discrete 2-body reactions
void reaction::two_body( data_parser &input_file )
{
  frame_out = CM;  // default
  for ( string dataID = input_file.get_dataID( ); dataID != "DONE";
    dataID = input_file.get_dataID( ) )
  {
    string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "reaction's q value" )
    {
      angdist.Q =  input_file.get_next_double( );
    }
    else if( tempID == "angular data" )
    {
      // set a temporary threshold based on the cross section data
      // this is reset later, based on the kinetics
      angdist.threshold = cross_section.begin( )->x;
      int num_Ein = input_file.get_next_int( );
      // read the angular probability density
      angdist.read_data( input_file, num_Ein );
    }
    else
    {
      FatalError( "reaction::two_body", 
		pastenum( "line ", input_file.line_count ) +
               ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }
  if( frame_out != CM )
  {
    FatalError( "reaction::two_body",
		"only emission in the center-of-mass frame is supported" );
  }
  // now get the transfer matrix
  transfer.allocate( );
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }
  angdist.particle_info.copy( particle_info );
  angdist.get_T( cross_section, model_weight, transfer );
}
// --------------- reaction::Legendre_two_body ----------------------
// Gets the transfer matrix from discrete 2-body Legendre polynomial angular data
void reaction::Legendre_two_body( data_parser &input_file )
{
  frame_out = CM;  // default
  for ( string dataID = input_file.get_dataID( ); dataID != "DONE";
    dataID = input_file.get_dataID( ) )
  {
    string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "reaction's q value" )
    {
      // not used by gammas from neutron capture
      LegendreAngle.Q =  input_file.get_next_double( );
    }
    else if( tempID == "legendre coefficients" )
    {
      int num_Ein = input_file.get_next_int( );
      // read the angular probability density
      if( particle_info.mProd == 0.0 )
      {
	captureGamma.read_data( input_file, num_Ein );
      }
      else
      {
        LegendreAngle.read_data( input_file, num_Ein );
      }
    }
    else
    {
      FatalError( "reaction::Legendre_two_body", 
		pastenum( "line ", input_file.line_count ) +
               ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }
  if( frame_out != CM )
  {
    FatalError( "reaction::Legendre_two_body",
		"only emission in the center-of-mass frame is supported" );
  }
  // now get the transfer matrix
  transfer.allocate( );
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }
  if( particle_info.mProd == 0.0 )
  {
    captureGamma.particle_info.copy( particle_info );
    captureGamma.get_T( cross_section, model_weight, transfer );
  }
  else
  {
    LegendreAngle.particle_info.copy( particle_info );
    LegendreAngle.get_T( cross_section, model_weight, transfer );
  }
}
// --------------- reaction::Legendre ----------------------
// Reads and processes Legendre expansions of double differential data
void reaction::Legendre( data_parser &input_file )
{
  frame_out = LAB;  // default
  for( string dataID = input_file.get_dataID( ); dataID != "DONE";
       dataID = input_file.get_dataID( ) )
  {
    string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "leeppdata" )
    {
      int num_moments = input_file.get_next_int( );
      // read the Legendre coefficients of the energy probability density
      energyMoments.read_data( input_file, num_moments );
    }
    else
    {
      FatalError( "reaction::Legendre", 
		pastenum( "line ", input_file.line_count ) +
                ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }
  if( frame_out != LAB )
  {
    FatalError( "reaction::Legendre",
		"only emission in the laboratory frame is supported" );
  }
  // now get the transfer matrix
  transfer.allocate( );
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }
  energyMoments.get_T( cross_section, multiple, model_weight, transfer );
}
// --------------- reaction::do_ENDFLegendre ----------------------
// Reads and processes Legendre expansions of double differential data
void reaction::do_ENDFLegendre( data_parser &input_file )
{
  two_d_interp Ein_interp;
  Interp_Type Eout_interp;
  frame_out = LAB;  // default
  for( string dataID = input_file.get_dataID( ); dataID != "DONE";
       dataID = input_file.get_dataID( ) )
  {
    string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "legendre data by incident energy" )
    {
      int num_Ein = input_file.get_next_int( );
      interp_flag_F::read_2d_interpolation( input_file, &Ein_interp, &Eout_interp );
      // read the Legendre coefficients of the energy probability density
      if( frame_out == LAB )
      {
        standard_legendre.Ein_interp = Ein_interp;
        standard_legendre.Eout_interp = Eout_interp;
	standard_legendre.read_data( input_file, num_Ein );
      }
      else
      {
        cm_Legendre_model.Ein_interp = Ein_interp;
        cm_Legendre_model.Eout_interp = Eout_interp;
	cm_Legendre_model.read_data( input_file, num_Ein );
      }
    }
    else
    {
      FatalError( "reaction::do_ENDFLegendre", 
		pastenum( "line ", input_file.line_count ) +
                ": dataID " + input_file.original_line +
                  " not implemented" );
    }
  }
  // now get the transfer matrix
  transfer.allocate( );
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }
  if( frame_out == LAB )
  {
    standard_legendre.get_T( cross_section, multiple, model_weight, transfer );
  }
  else
  {
    cm_Legendre_model.particles.copy( particle_info );
    cm_Legendre_model.get_T( cross_section, multiple, model_weight, transfer );
  }
}
// --------------- reaction::do_isotropic ----------------------
// Reads and processes tables of isotropic energy probability density
void reaction::do_isotropic( data_parser &input_file )
{
  two_d_interp Ein_interp;
  Interp_Type Eout_interp;
  frame_out = LAB;  // default

  for( string dataID = input_file.get_dataID( ); dataID != "DONE";
       dataID = input_file.get_dataID( ) )
  {
    string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "eeppdata" )
    {
      int num_Ein = input_file.get_next_int( );
      interp_flag_F::read_2d_interpolation( input_file, &Ein_interp, &Eout_interp );
      if( frame_out == CM )
      {
        cm_Legendre_model.Ein_interp = Ein_interp;
        cm_Legendre_model.Eout_interp = Eout_interp;
	cm_Legendre_model.read_isotropic( input_file, num_Ein );
      }
      else
      {
        isotrop.Ein_interp = Ein_interp;
        isotrop.Eout_interp = Eout_interp;
        isotrop.read_data( input_file, num_Ein );
      }
    }
    else
    {
      FatalError( "reaction::do_isotropic", 
		pastenum( "line ", input_file.line_count ) +
               ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }

  // now get the transfer matrix
  transfer.allocate( );
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }
  if( frame_out == LAB )
  {
    isotrop.get_T( cross_section, multiple, model_weight, transfer );
  }
  else
  {
    cm_Legendre_model.particles.copy( particle_info );
    cm_Legendre_model.get_T( cross_section, multiple, model_weight, transfer );
  }
}
// --------------- reaction::do_uncorr ----------------------
// Reads and processes uncorrelated expansions of double differential data
void reaction::do_uncorr( data_parser &input_file )
{
  two_d_interp Ein_interp;
  Interp_Type Eout_interp;
  frame_out = LAB;  // default
  bool is_isotropic = false;    // Is the angular data isotropic?
  for( string dataID = input_file.get_dataID( ); dataID != "DONE";
       dataID = input_file.get_dataID( ) )
  {
    string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "angular data" )
    {
      uncorr.mu_table = true;
      int num_Ein = input_file.get_next_int( );
      // read the angular probability density
      angdist.threshold = cross_section.begin( )->x;
      angdist.read_data( input_file, num_Ein );
      is_isotropic = angdist.isotropic( );
    }
    else if( tempID == "legendre coefficients" )
    {
      uncorr.mu_table = false;
      int num_Ein = input_file.get_next_int( );
      // read the Legendre coefficients
      uncorr.read_Legendre( input_file, num_Ein );
    }
    else if( tempID == "eeppdata" )
    {
      int num_Ein = input_file.get_next_int( );
      interp_flag_F::read_2d_interpolation( input_file, &Ein_interp, &Eout_interp );
      // read the energy probability density
      if( is_isotropic && ( frame_out == CM ) )
      {
        cm_Legendre_model.Ein_interp = Ein_interp;
        cm_Legendre_model.Eout_interp = Eout_interp;
	cm_Legendre_model.read_isotropic( input_file, num_Ein );
      }
      else
      {
        uncorr.Ein_interp = Ein_interp;
        uncorr.Eout_interp = Eout_interp;
	uncorr.read_data( input_file, num_Ein, &angdist );
      }
    }
    else
    {
      FatalError( "reaction::do_uncorr", 
		pastenum( "line ", input_file.line_count ) +
               ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }

  // now get the transfer matrix
  transfer.allocate( );
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }
  if( frame_out == LAB )
  {
    uncorr.get_T( cross_section, multiple, model_weight, transfer );
  }
  else
  {
    if( !is_isotropic )
    {
      FatalError( "reaction::do_uncorr",
        "center-of-mass data is implemented only for isotropic data" );
    }
    cm_Legendre_model.particles.copy( particle_info );
    cm_Legendre_model.get_T( cross_section, multiple, model_weight, transfer );
  }
}
// --------------- reaction::do_joint_dist ----------------------
// Reads and processes double-differential tabular data
void reaction::do_joint_dist( data_parser &input_file )
{
  frame_out = LAB;  // default
  for( string dataID = input_file.get_dataID( ); dataID != "DONE";
       dataID = input_file.get_dataID( ) )
  {
    string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "angular data" )
    {
      int num_Ein = input_file.get_next_int( );
      joint_dist_table.angle_data.threshold = cross_section.begin( )->x;
      // read the angular probability density
      joint_dist_table.angle_data.read_data( input_file, num_Ein );
    }
    else if( tempID == "emueppdata" )
    {
      int num_Eout = input_file.get_next_int( );
      // read the double-differential data tables
      joint_dist_table.read_data( input_file, num_Eout );
    }
    else
    {
      FatalError( "reaction::do_joint_dist", 
		pastenum( "line ", input_file.line_count ) +
                ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }
  if( frame_out != LAB )
  {
    FatalError( "reaction::do_joint_dist",
		"only emission in the laboratory frame is supported" );
  }
  // now get the transfer matrix
  transfer.allocate( );
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }
  joint_dist_table.get_T( cross_section, multiple, model_weight, transfer );
}
// --------------- reaction::do_Compton ----------------------
// Reads and processes scattering factor data for Compton scattering
void reaction::do_Compton( data_parser &input_file )
{
  // set defaults
  frame_out = LAB;
  transfer.mu_quad_method = ADAPT_HALF; // adaptive Gauss with sqrt(1-x) singularity
  Global.set( "scale_rows", 0 );   // the check is invalid---sigma not lin-lin

  for( string dataID = input_file.get_dataID( ); dataID != "DONE";
       dataID = input_file.get_dataID( ) )
  {
    string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "scatteringfactordata" )
    {
      int num_factor = input_file.get_next_int( );
      // read the scattering factor data
      compton.file_data.read_data_interp( input_file, num_factor );
    }
    else
    {
      FatalError( "reaction::do_Compton", 
		pastenum( "line ", input_file.line_count ) +
                ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }
  if( frame_out != LAB )
  {
    FatalError( "reaction::do_Compton",
		"only emission in the laboratory frame is supported" );
  }
  // now get the transfer matrix
  transfer.allocate( );
  compton.get_T( transfer, cross_section );
  write_xsec( );
}
// --------------- reaction::do_coherent ----------------------
// Reads and processes scattering factor data for coherent scattering
void reaction::do_coherent( data_parser &input_file )
{
  // set defaults
  frame_out = LAB;
  transfer.mu_quad_method = ADAPT_HALF; // adaptive Gauss with sqrt(1-x) singularity
  Global.set( "scale_rows", 0 );   // the check is invalid---sigma not lin-lin

  for( string dataID = input_file.get_dataID( ); dataID != "DONE";
       dataID = input_file.get_dataID( ) )
  {
    string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "formfactordata" )
    {
      int num_factor = input_file.get_next_int( );
      // read the form factors
      coherent_model.file_data.read_data_interp( input_file, num_factor );
    }
    else if( tempID == "realanomalousfactor" )
    {
      int num_factor = input_file.get_next_int( );
      // read the real anomalous factor
      coherent_model.realAnomalous.read_data_interp( input_file, num_factor );
    }
    else if( tempID == "imaginaryanomalousfactor" )
    {
      int num_factor = input_file.get_next_int( );
      // read the imaginary anomalous factor
      coherent_model.imaginaryAnomalous.read_data_interp( input_file, num_factor );
    }
    else
    {
      FatalError( "reaction::do_coherent", 
		pastenum( "line ", input_file.line_count ) +
                ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }
  if( frame_out != LAB )
  {
    FatalError( "reaction::do_coherent",
		"only emission in the laboratory frame is supported" );
  }
  // now get the transfer matrix
  transfer.allocate( );
  coherent_model.get_T( transfer, cross_section );
  write_xsec( );
}
// --------------- reaction::do_evaporation ----------------------
// Reads and processes data for evaporation spectra
void reaction::do_evaporation( data_parser &input_file )
{
  frame_out = LAB;  // default
  for( string dataID = input_file.get_dataID( ); dataID != "DONE";
       dataID = input_file.get_dataID( ) )
  {
    string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "u" )
    {
      evaporation_model.U = input_file.get_next_double( );
    }
    else if( tempID == "theta" )
    {
      int num_Theta = input_file.get_next_int( );
      // read the Theta values
      evaporation_model.read_data_interp( input_file, num_Theta );
    }
    else
    {
      FatalError( "reaction::do_evaporation", 
		pastenum( "line ", input_file.line_count ) +
                ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }
  if( frame_out != LAB )
  {
    FatalError( "reaction::do_evaporation",
		"only emission in the laboratory frame is supported" );
  }
  // now get the transfer matrix
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }
  transfer.allocate( );
  evaporation_model.get_T( cross_section, multiple, model_weight, transfer );
}
// --------------- reaction::do_Maxwell ----------------------
// Reads and processes data for Maxwell spectra
void reaction::do_Maxwell( data_parser &input_file )
{
  frame_out = LAB;  // default
  for( string dataID = input_file.get_dataID( ); dataID != "DONE";
       dataID = input_file.get_dataID( ) )
  {
    string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "u" )
    {
      Maxwell_model.U = input_file.get_next_double( );
    }
    else if( tempID == "theta" )
    {
      int num_Theta = input_file.get_next_int( );
      // read the Theta values
      Maxwell_model.read_data_interp( input_file, num_Theta );
    }
    else
    {
      FatalError( "reaction::do_Maxwell", 
		pastenum( "line ", input_file.line_count ) +
                ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }
  if( frame_out != LAB )
  {
    FatalError( "reaction::do_Maxwell",
		"only emission in the laboratory frame is supported" );
  }
  // now get the transfer matrix
  transfer.allocate( );
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }
  Maxwell_model.get_T( cross_section, multiple, model_weight, transfer );
}
// --------------- reaction::do_Watt ----------------------
// Reads and processes data for Watt spectra
void reaction::do_Watt( data_parser &input_file )
{
  frame_out = LAB;  // default
  for( string dataID = input_file.get_dataID( ); dataID != "DONE";
       dataID = input_file.get_dataID( ) )
  {
    string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "u" )
    {
      Watt_model.U = input_file.get_next_double( );
    }
    else if( tempID == "b" )
    {
      int num_b = input_file.get_next_int( );
      // read the b values
      Watt_model.b_data.read_data_interp( input_file, num_b );
    }
    else if( ( tempID == "theta" ) || ( tempID == "a" ) )
    {
      int num_Theta = input_file.get_next_int( );
      // read the Theta, i.e., a values
      Watt_model.read_data_interp( input_file, num_Theta );
    }
    else
    {
      FatalError( "reaction::do_Watt", 
		pastenum( "line ", input_file.line_count ) +
                ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }
  if( frame_out != LAB )
  {
    FatalError( "reaction::do_Watt",
		"only emission in the laboratory frame is supported" );
  }
  // now get the transfer matrix
  transfer.allocate( );
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }
  Watt_model.get_T( cross_section, multiple, model_weight, transfer );
}
// --------------- reaction::do_MadlandNix ----------------------
// Reads and processes data for MadlandNix spectra
void reaction::do_MadlandNix( data_parser &input_file )
{
  double Efl = 0.0;  // kinetic energy of the lighter fragment
  double Efh = 0.0;  // kinetic energy of the heavier fragment
  frame_out = LAB;  // default
  // The Madland-Nix model does not specify a maximum outgoing energy
  MadlandNix_model.maxEout = 0.0;  // not set, default
  for( string dataID = input_file.get_dataID( ); dataID != "DONE";
       dataID = input_file.get_dataID( ) )
  {
    string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "maxeout" )
    {
      MadlandNix_model.maxEout = input_file.get_next_double( );
    }
    else if( tempID == "efl" )
    {
      Efl = input_file.get_next_double( );
      MadlandNix_model.Efl = Efl;
    }
    else if( tempID == "efh" )
    {
      Efh = input_file.get_next_double( );
      MadlandNix_model.Efh = Efh;
    }
    else if( tempID == "u" )
    {
      MadlandNix_model.U = input_file.get_next_double( );
      Warning( "reaction::do_MadlandNix",
        "U is not used in the Madland-Nix model." );
    }
    else if( tempID == "tm" )
    {
      int num_TM = input_file.get_next_int( );
      // read the TM values
      MadlandNix_model.read_data_interp( input_file, num_TM );
    }
    else
    {
      FatalError( "reaction::do_MadlandNix", 
		pastenum( "line ", input_file.line_count ) +
                ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }
  if( ( Efl <= 0.0 ) || ( Efh <= 0.0 ) )
  {
    FatalError( "reaction::do_MadlandNix",
      "The kinetic energies of the fission fragments should be positive." );
  }
  if( transfer.conserve != NUMBER )
  {
    Warning( "reaction::do_MadlandNix", 
	     "Energy-conserving transfer matrices not implemented." );
    transfer.conserve = NUMBER;
  }
  if( ( transfer.Ein_quad_method != ADAPTIVE2 ) &&
      ( transfer.Ein_quad_method != ADAPTIVE4 ) )
  {
    Warning( "reaction::do_MadlandNix", 
	     "Gaussian quadrature is inaccurate in the low outgoing energy bins." );
  }
  if( frame_out != LAB )
  {
    FatalError( "reaction::do_MadlandNix",
		"only emission in the laboratory frame is supported" );
  }
  // set up the transfer matrix
  transfer.allocate( );
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }
  MadlandNix_model.get_T( cross_section, multiple, model_weight, transfer );
}
// --------------- reaction::do_Kalbach ----------------------
// Reads and processes data for Kalbach spectra
void reaction::do_Kalbach( data_parser &input_file )
{
  frame_out = CM;  // default
  for( string dataID = input_file.get_dataID( ); dataID != "DONE";
       dataID = input_file.get_dataID( ) )
  {
    string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "projectile's za" )
    {
      Kalbach_model.kalbach_a.projectile.set_params( input_file.get_next_int( ) );
    }
    else if( tempID == "target's za" )
    {
      Kalbach_model.kalbach_a.target.set_params( input_file.get_next_int( ) );
    }
    else if( tempID == "compound's mass" )
    {
      Kalbach_model.kalbach_a.compound.mass = input_file.get_next_double( );
    }
    else if( tempID == "product's za" )
    {
      Kalbach_model.kalbach_a.eject.set_params( input_file.get_next_int( ) );
    }
    else if( tempID == "kalbach probabilities" )
    {
      int num_Ein = input_file.get_next_int( );
      Kalbach_model.read_probability( input_file, num_Ein );
    }
    else if( tempID == "kalbach r parameter" )
    {
      int num_Ein = input_file.get_next_int( );
      Kalbach_model.read_r( input_file, num_Ein );
    }
    else
    {
      FatalError( "reaction::do_Kalbach", 
		pastenum( "line ", input_file.line_count ) +
                ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }
  if( frame_out != CM )
  {
    FatalError( "reaction::do_Kalbach",
		"only emission in the center-of-mass frame is supported" );
  }
  // now get the transfer matrix
  Kalbach_model.kalbach_a.copy_masses( particle_info );
  transfer.allocate( );
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }
  Kalbach_model.get_T( cross_section, multiple, model_weight, transfer );
}
// --------------- reaction::do_phase_space ----------------------
// Reads and processes data for phase_space spectra
void reaction::do_phase_space( data_parser &input_file )
{
  frame_out = CM;  // default
  for( string dataID = input_file.get_dataID( ); dataID != "DONE";
       dataID = input_file.get_dataID( ) )
  {
    string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "total mass" )
    {
      phase_space_model.totalMass = input_file.get_next_double( );
    }
    else if( tempID == "number of particles" )
    {
      phase_space_model.numParticles = input_file.get_next_int( );
    }
    else if( tempID == "q value" )
    {
      phase_space_model.Q_value = input_file.get_next_double( );
    }
    else
    {
      FatalError( "reaction::do_phase_space", 
		pastenum( "line ", input_file.line_count ) +
                ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }
  if( ( transfer.Ein_quad_method != ADAPTIVE2 ) &&
      ( transfer.Ein_quad_method != ADAPTIVE4 ) )
  {
    Warning( "reaction::do_phase_space", 
	     "Gaussian quadrature is inaccurate in the low incident energy bins." );
  }
  if( frame_out != CM )
  {
    FatalError( "reaction::do_phase_space",
		"only emission in the center-of-mass frame is supported" );
  }
  // now get the transfer matrix
  phase_space_model.copy_masses( particle_info );
  transfer.allocate( );
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }
  phase_space_model.get_T( cross_section, multiple, model_weight, transfer );
}
// --------------- reaction::do_gen_evap ----------------------
// Reads and processes data for gen_evap spectra
void reaction::do_gen_evap( data_parser &input_file )
{
  frame_out = LAB;  // default
  for( string dataID = input_file.get_dataID( ); dataID != "DONE";
       dataID = input_file.get_dataID( ) )
  {
    string tempID = dataID;
    string_F::Tolower( tempID );
    // Is this data common to all reactions?
    if( common_input( tempID, input_file ) )
    {
      continue;
    }
    // what special type of data is it?
    else if( tempID == "g" )
    {
      int num_g = input_file.get_next_int( );
      // read the g values
      gen_evap_model.read_data_interp( input_file, num_g );
    }
    else if( tempID == "theta" )
    {
      int num_theta = input_file.get_next_int( );
      // read the theta values
      gen_evap_model.theta.read_data_interp( input_file, num_theta );
    }
    else if( tempID == "u" )
    {
      double U = input_file.get_next_double( );  // not used
      if( U != 0.0 )
      {
	Info( "reaction::do_gen_evap",
	      "The parameter U is not used in this model" );
      }
    }
    else
    {
      FatalError( "reaction::do_gen_evap", 
		pastenum( "line ", input_file.line_count ) +
                ": dataID " + input_file.original_line +
		  " not implemented" );
    }
  }
  if( frame_out != LAB )
  {
    FatalError( "reaction::do_gen_evap",
		"only emission in the laboratory frame is supported" );
  }
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }
  // now get the transfer matrix
  transfer.allocate( );
  gen_evap_model.get_T( cross_section, multiple, model_weight, transfer );
}
// --------------- reaction::write_transfer -----------------
// Writes the output data
void reaction::write_transfer( )
{
  // do the proper scaling
  transfer.get_flux_weight( );
  transfer.use_weight( );

  if( cross_section.size( ) > 0 )
  {
    // check the row sums against the integrals
    number_check( );
    transfer.check_ell0( );
  }
  else
  {
    Warning( "reaction::write_transfer", "Compute the cross section for this reaction" );
  }
  transfer.write_transfer( );
}
// --------------- reaction::write_xsec -----------------
// Writes the cross section for gamma data
void reaction::write_xsec( )
{
  int num_xsec = cross_section.size( );
  *transfer.output_file << "Cross section: n = " << num_xsec << endl;
  *transfer.output_file << "Interpolation: linear-linear" << endl;
  transfer.output_file->setf(ios::scientific,ios::floatfield);
  static int data_precision = Global.Value( "datafield_precision" );
  static int field_width = Global.get_field_width( );
  for( dd_vector::const_iterator xsec = cross_section.begin( ); 
       xsec != cross_section.end( ); ++xsec )
  {
    *transfer.output_file << setw(field_width) <<
            setprecision(data_precision) <<
      xsec->x << " " << xsec->y << endl;
  }
}
// --------------- reaction::number_check -----------------
// Evaluates the integrals of multiplicity * cross section * flux
void reaction::number_check( )
{
  num_check_param num_check;

  // pointers to the cross section
  num_check.sigma.Eout_interp = cross_section.interp_type;
  dd_vector::const_iterator this_sigma = cross_section.begin( );
  dd_vector::const_iterator next_sigma = this_sigma;
  ++next_sigma;

  // pointers to the multiplicity
  if( multiple.size( ) == 0 )
  {
    multiple.make_flat( cross_section, 1.0 );
  }
  num_check.mult.Eout_interp = multiple.interp_type;
  dd_vector::const_iterator this_mult = multiple.begin( );
  dd_vector::const_iterator next_mult = this_mult;
  ++next_mult;

  // pointers to the flux data
  num_check.e_flux.Eout_interp = transfer.e_flux.interp;
  Flux_List::const_iterator flux_ptr = transfer.e_flux.begin( );
  Flux_List::const_iterator next_flux = flux_ptr;
  ++next_flux;

  // pointers to the model weight
  if( model_weight.size( ) == 0 )
  {
    model_weight.make_flat( cross_section, 1.0 );
  }
  num_check.model_weight.Eout_interp = model_weight.interp_type;
  dd_vector::const_iterator this_model_weight = model_weight.begin( );
  dd_vector::const_iterator next_model_weight = this_model_weight;
  ++next_model_weight;

  // count through the incident energy groups
  int Ein_count = 0;
  // which energy group
  vector< double >::const_iterator Ein_ptr = transfer.in_groups.begin( );
  vector< double >::const_iterator next_Ein = Ein_ptr;
  ++next_Ein;

  // synchronize the first energies
  // we may need to increment num_check.sigma
  for( ; next_sigma->x <= *Ein_ptr; 
    this_sigma = next_sigma, ++next_sigma )
  {
    if( next_sigma == cross_section.end( ) )
    {
      Warning( "reaction::number_check",
        "omit reaction: all energy groups too high" );
      transfer.zero_transfer( );
    }
  }
  num_check.sigma.set_pair( *this_sigma, *next_sigma );
  // we may need to increment num_check.mult
  for( ; next_mult->x <= *Ein_ptr; 
    this_mult = next_mult, ++next_mult )
  {
    if( next_mult == multiple.end( ) )
    {
      Warning( "reaction::number_check",
        "omit reaction: no particle multiplicity given" );
      transfer.zero_transfer( );
    }
  }
  num_check.mult.set_pair( *this_mult, *next_mult );
  // we may need to increment the energy bin
  for( ; *next_Ein <= this_sigma->x;
       ++Ein_count, Ein_ptr = next_Ein, ++next_Ein )
  {
    if( next_Ein == transfer.in_groups.end( ) )
    {
      Warning( "reaction::number_check", 
        "omit reaction: threshold too high" );
      transfer.zero_transfer( );
    }
  }
  // we may need to increment the incident flux
  for( ; next_flux->get_E_in( ) <= this_sigma->x;
       flux_ptr = next_flux, ++next_flux )
  {
    if( next_flux == transfer.e_flux.end( ) )
    {
      Warning( "reaction::number_check",
        "omit reaction: ran out of flux data" );
      transfer.zero_transfer( );
    }
  }
  num_check.e_flux.first.x = flux_ptr->get_E_in( );
  num_check.e_flux.first.y = flux_ptr->data[ 0 ];
  num_check.e_flux.second.x = next_flux->get_E_in( );
  num_check.e_flux.second.y = next_flux->data[ 0 ];
  // assume here that the model weight is consistent
  num_check.model_weight.set_pair( *this_model_weight, *next_model_weight );

  // now do the integrals
  for( ; ; )
  {
    // the range of integration
    double E_in0;
    double E_in1;
    E_in0 = ( this_sigma->x < *Ein_ptr ) ? *Ein_ptr :
      this_sigma->x;
    E_in0 = ( this_mult->x < E_in0 ) ? E_in0 :
      this_mult->x;
    E_in0 = ( this_model_weight->x < E_in0 ) ? E_in0 :
      this_model_weight->x;
    E_in0 = ( flux_ptr->get_E_in( ) < E_in0 ) ? E_in0 :
      flux_ptr->get_E_in( );
    E_in1 = ( next_sigma->x > *next_Ein ) ? *next_Ein :
      next_sigma->x;
    E_in1 = ( next_mult->x > E_in1 ) ? E_in1 :
      next_mult->x;
    E_in1 = ( next_model_weight->x > E_in1 ) ? E_in1 :
      next_model_weight->x;
    E_in1 = ( next_flux->get_E_in( ) > E_in1 ) ? E_in1 :
      next_flux->get_E_in( );

    // Use Simpson's rule
    double Simpson = ( ( E_in1 - E_in0 ) / 6.0 ) *
      ( num_check.value( E_in0 ) +
	4.0 * num_check.value( 0.5 * ( E_in0 + E_in1 ) ) +
	num_check.value( E_in1 ) );
    transfer.row_checks[ Ein_count ].weight_1[ 0 ] += Simpson;

    // got to the next interval
    if( E_in1 == next_sigma->x )
    {
      // increment num_check.sigma
      this_sigma = next_sigma;
      ++next_sigma;
      if( next_sigma == cross_section.end( ) )
      {
        break;
      }
      num_check.sigma.set_pair( *this_sigma, *next_sigma );
    }
    if( E_in1 == next_mult->x )
    {
      // increment num_check.mult
      this_mult = next_mult;
      ++next_mult;
      if( next_mult == multiple.end( ) )
      {
        break;
      }
      num_check.mult.set_pair( *this_mult, *next_mult );
    }
    if( E_in1 == next_model_weight->x )
    {
      // increment num_check.model_weight
      this_model_weight = next_model_weight;
      ++next_model_weight;
      if( next_model_weight == model_weight.end( ) )
      {
        break;
      }
      num_check.model_weight.set_pair( *this_model_weight, *next_model_weight );
    }
    if( E_in1 == *next_Ein )
    {
      // go to the next incident energy bin
      ++Ein_count;
      Ein_ptr = next_Ein;
      ++next_Ein;
      if( next_Ein == transfer.in_groups.end( ) )
      {
        break;
      }
    }
    if( E_in1 == next_flux->get_E_in( ) )
    {
      // increment flux
      flux_ptr = next_flux;
      ++next_flux;
      if( next_flux == transfer.e_flux.end( ) )
      {
        break;
      }
      num_check.e_flux.first.x = flux_ptr->get_E_in( );
      num_check.e_flux.first.y = flux_ptr->data[ 0 ];
      num_check.e_flux.second.x = next_flux->get_E_in( );
      num_check.e_flux.second.y = next_flux->data[ 0 ];
    }
  }
  // scale the row sums by the flux weights
  transfer.scale_row_check( );
}
