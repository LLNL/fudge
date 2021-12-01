/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id:messaging.cpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/

#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include "global_params.hpp"
#include "messaging.hpp"


//! Always write out a Fatal Error message, then exit badly
void Msg::FatalError(std::string routine, std::string error){
    std::string msg = routine+": "+error;
    std::cerr << "merced Fatal Error "<<msg<<std::endl;
    //  messageLog.logMessage(msg, 4);
    //  messageLog.write();
    exit(-1);
};

//! Always write out a Severe Error message, then attempt to throw an exception
void Msg::SevereError(std::string routine, std::string error){
    std::string msg = "merced "+routine+": "+error;
    //  messageLog.logMessage(msg, 3);
    //    throw(std::string("Severe Error "+msg));
    std::cerr << "Severe Error "+msg << std::endl;
    exit( -2 );
};

//! Only show a warning if required
void Msg::Warning(std::string routine, std::string warning){
  static bool messagelevel_warn = Global.Value( "message_level" ) >= 2;
    if ( messagelevel_warn) {
        std::string msg = "merced "+routine+": "+warning;
        std::cerr << "Msg::Warning "<<msg<<std::endl;
        //  messageLog.logMessage(msg, 2);
    }
};

//! Only show information if required
void Msg::Info(std::string routine, std::string info){
  static bool messagelevel_info = Global.Value( "message_level" ) >= 3;
    if (messagelevel_info) {
        std::string msg = "merced "+routine+": "+info;
        std::cout << "Msg::Info "<<msg<<std::endl;
        //messageLog.logMessage(msg, 0); // never log Msg::Info statements
    }
};

//! Only show information if required
void Msg::DebugInfo(std::string routine, std::string info){
  static bool messagelevel_info = Global.Value( "message_level" ) == 1;
    if (messagelevel_info) {
        std::string msg = "merced "+routine+": "+info;
        std::cout << "Msg::DebugInfo "<<msg<<std::endl;
        //messageLog.logMessage(msg, 5);
    }
};

//! Only show information if required
void Msg::Unimplemented(std::string routine, std::string info){
  static bool messagelevel_info = Global.Value( "message_level" ) >= 1;
    if (messagelevel_info) {
        std::string msg = "merced "+routine+": "+info;
        std::cout << "Msg::Unimplemented "<<msg<<std::endl;
        //  messageLog.logMessage(msg, 1);
    }
};

//! gizmo to paste doubles on the backs of std::strings
std::string Msg::pastenum(std::string front, double back){
  std::ostringstream tmp; 
    tmp << front << back; 
    return tmp.str();
};

//! gizmo to paste ints on the backs of std::strings
std::string Msg::pastenum(std::string front, int back){
  std::ostringstream tmp; 
    tmp << front << back; 
    return tmp.str();
};
