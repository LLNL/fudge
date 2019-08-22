/*
 * ******** merced: calculate the transfer matrix *********
 * $Revision: 1 $
 * $Date: 2006-02-01 19:06:56 -0800 (Wed, 01 Feb 2006) $
 * $Author: hedstrom $
 * $Id:messaging.hpp 1 2006-02-02 03:06:56Z hedstrom $
 * ******** merced: calculate the transfer matrix *********
 *
 * # <<BEGIN-copyright>>
 * # <<END-copyright>>
*/

// messaging prototypes

#ifndef MESSAGING
#define MESSAGING

#include <string>
using namespace std;

//! Flags a fatal error
//! \param routine the name of the routine in which the error occurs
//! \param error the type of error
void FatalError(string routine, string error);

//! Flags a severe error
//! \param routine the name of the routine in which the error occurs
//! \param error the type of error
void SevereError(string routine, string error);

//! Flags a warning
//! \param routine the name of the routine in which the warning is given
//! \param warning the type of warning
void Warning(string routine, string warning);

//! Flags an information message
//! \param routine the name of the routine in which the information is given
//! \param info the information message
void Info(string routine, string info);

//! Flags an unimplemented routine
//! \param routine the name of the unimplemented routine
//! \param info the information message
void Unimplemented(string routine, string info);

//! Attaches a double to a string
//! \param front the initial string
//! \param back the appended double
string pastenum(string front, double back);

//! Attaches an integer to a string
//! \param front the initial string
//! \param back the appended integer
string pastenum(string front, int back);

#endif
