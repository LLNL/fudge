/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id:messaging.cpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
  Copyright (c) 2017, Lawrence Livermore National Security, LLC.
  Produced at the Lawrence Livermore National Laboratory.
  Written by the LLNL Nuclear Data and Theory group
          (email: mattoon1@llnl.gov)
  LLNL-CODE-725546.
  All rights reserved.
  
  This file is part of the Merced package, used to generate nuclear reaction
  transfer matrices for deterministic radiation transport.
  
  
      Please also read this link - Our Notice and Modified BSD License
  
  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:
      * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the disclaimer below.
      * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the disclaimer (as noted below) in the
        documentation and/or other materials provided with the distribution.
      * Neither the name of LLNS/LLNL nor the names of its contributors may be used
        to endorse or promote products derived from this software without specific
        prior written permission.
  
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
  THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  
  
  Additional BSD Notice
  
  1. This notice is required to be provided under our contract with the U.S.
  Department of Energy (DOE). This work was produced at Lawrence Livermore
  National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
  
  2. Neither the United States Government nor Lawrence Livermore National Security,
  LLC nor any of their employees, makes any warranty, express or implied, or assumes
  any liability or responsibility for the accuracy, completeness, or usefulness of any
  information, apparatus, product, or process disclosed, or represents that its use
  would not infringe privately-owned rights.
  
  3. Also, reference herein to any specific commercial products, process, or services
  by trade name, trademark, manufacturer or otherwise does not necessarily constitute
  or imply its endorsement, recommendation, or favoring by the United States Government
  or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
  herein do not necessarily state or reflect those of the United States Government or
  Lawrence Livermore National Security, LLC, and shall not be used for advertising or
  product endorsement purposes.
  
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
