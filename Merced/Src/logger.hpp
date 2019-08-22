/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id:logger.hpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/

#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <map>
#include <string>
#include <iostream>
#include "version.hpp"
#include "global_params.hpp"

using namespace std;

//! Severity of the problem and the message that goes with it.  Use this as the value in the map.
class logEntry{
public:
    string message;
    int severity;
    logEntry(string msg, int sev): message(msg), severity(sev){}
    friend ostream& operator<<(ostream& out, const logEntry& ent);
};

//! Map to hold the list of problems.  We append this list to the documentation of the evaluation.
class MessageLogger : public map< logEntry > {
public:
    string inLogFile;
    string outLogFile;
    MessageLogger(string infile="documentation.txt", string outfile="documentation.txt"): 
        inLogFile(infile), outLogFile(outfile){}
    ~MessageLogger(void){}
    bool write(void);
    void logMessage(string msg, int severity)
        {insert(logEntry(msg,severity));}
};

extern MessageLogger messageLog;

#endif
