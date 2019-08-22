/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: logger.cpp 1 2006-02-02 03:06:56Z hedstrom $
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
#include <fstream>
#include <string>
#include <sstream>
#include "version.hpp"
#include "logger.hpp"
#include "convert.hpp"
#include "global_params.hpp"

//! output for logEntry
ostream& operator<<(ostream& out, const logEntry& ent){
    switch (ent.severity) {
        case 0: 
            out << rjust("Info",14); 
            break; 
        case 1: 
            out << rjust("Unimplemented",14);  
            break; 
        case 2: 
            out << rjust("Warning",14); 
            break; 
        case 3: 
            out << rjust("Severe Error",14);
            break; 
        case 4: 
            out << rjust("Fatal Error",14); 
    }
    out <<", "<< ent.message;
    return out;
}

//! main output routine for the MessageLogger
bool MessageLogger::write(void){
    bool isEndfDoc(false);
    stringstream oldStuff("");
    // read input file in scope separate from output file's scope so can use same file
    {
        ifstream inStream(inLogFile.c_str());
        char buffer[181];
        if (inStream) {
            while (inStream.getline(buffer,sizeof(buffer))){oldStuff<<buffer;oldStuff<<"\n";}
            isEndfDoc=true;
        } else {
            cout << "Info: Could not read "<<inLogFile<<", will not include ENDF documentation in log file"<<endl;
        }
    }
    ofstream outStream(outLogFile.c_str());
    static int skip_logging = Global.Value( "skip_logging" );
    if (!outStream) return false;
    if (isEndfDoc) outStream << oldStuff.str() << endl << endl;
    if (skip_logging>0) return true;
    //Pretty banner to set off start of program & print the version information
    outStream << 11*string(" ") << "+" << 55*string("-") << "+"<<endl;
    outStream << 11*string(" ") << "|" << center("merced: calculate the transfer matrix",55)<< "|"<<endl;
    outStream << 11*string(" ") << "|" << center("version "+string(MERCED_VERSION),55) << "|" <<endl;
    outStream << 11*string(" ") << "|" << center("revision "+string(MERCED_REVISION),55) << "|" <<endl;
    outStream << 11*string(" ") << "|" << center("last SVN commit "+string(MERCED_VERSION_DATE),55) << "|" <<endl;
    outStream << 11*string(" ") << "+" << 55*string("-") << "+"<<endl;
    outStream << endl;
    if (this->size()>0) {
        outStream << center(60*string("="),62) << endl;
        outStream << center("Problems encountered during energy deposition calculation",62) << endl;
        outStream << center(60*string("="),62) << endl;
        for (MessageLogger::const_iterator it=begin(); it!=end(); ++it){
            stringstream buf;
            buf<<it->first;
            outStream << ljust(buf.str(),35) << it->second << endl;
        }
    } else {
        outStream << center(60*string("="),62) << endl;
        outStream << center("No problems encountered during energy deposition calculation!!!",62) << endl;
        outStream << center(60*string("="),62) << endl;
    }
    return true;
}
