/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id: logger.cpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
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
