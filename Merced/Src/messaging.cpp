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
// #include "logger.hpp"

using namespace std;

//! Always write out a Fatal Error message, then exit badly
void FatalError(string routine, string error){
    string msg = routine+": "+error;
    cerr << "merced Fatal Error "<<msg<<endl;
    //  messageLog.logMessage(msg, 4);
    //  messageLog.write();
    exit(-1);
};

//! Always write out a Severe Error message, then attempt to throw an exception
void SevereError(string routine, string error){
    string msg = "merced "+routine+": "+error;
    //  messageLog.logMessage(msg, 3);
    //    throw(string("Severe Error "+msg));
    cerr << "Severe Error "+msg << endl;
    exit( -2 );
};

//! Only show a warning if required
void Warning(string routine, string warning){
  static bool messagelevel_warn = Global.Value( "message_level" ) <= 1;
    if ( messagelevel_warn) {
        string msg = "merced "+routine+": "+warning;
        cerr << "Warning "<<msg<<endl;
        //  messageLog.logMessage(msg, 2);
    }
};

//! Only show information if required
void Info(string routine, string info){
  static bool messagelevel_info = Global.Value( "message_level" ) <= 0;
    if (messagelevel_info) {
        string msg = "merced "+routine+": "+info;
        cout << "Info "<<msg<<endl;
        //messageLog.logMessage(msg, 0); // never log Info statements
    }
};

//! Only show information if required
void Unimplemented(string routine, string info){
  static bool messagelevel_info = Global.Value( "message_level" ) <= 0;
    if (messagelevel_info) {
        string msg = "merced "+routine+": "+info;
        cout << "Unimplemented "<<msg<<endl;
        //  messageLog.logMessage(msg, 1);
    }
};

//! gizmo to paste doubles on the backs of strings
string pastenum(string front, double back){
    ostringstream tmp; 
    tmp << front << back; 
    return tmp.str();
};

//! gizmo to paste ints on the backs of strings
string pastenum(string front, int back){
    ostringstream tmp; 
    tmp << front << back; 
    return tmp.str();
};
