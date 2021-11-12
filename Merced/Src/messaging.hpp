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

namespace Msg
{

//! Flags a fatal error, always printed
//! \param routine the name of the routine in which the error occurs
//! \param error the type of error
void FatalError(std::string routine, std::string error);

//! Flags a severe error, always printed
//! \param routine the name of the routine in which the error occurs
//! \param error the type of error
void SevereError(std::string routine, std::string error);

//! Flags a debugging information message if message_level >= 1
//! \param routine the name of the routine in which the information is given
//! \param info the information message
void DebugInfo(std::string routine, std::string info);

//! Flags a warning if message_level >= 2
//! \param routine the name of the routine in which the warning is given
//! \param warning the type of warning
void Warning(std::string routine, std::string warning);

//! Flags an information message if message_level >= 3
//! \param routine the name of the routine in which the information is given
//! \param info the information message
void Info(std::string routine, std::string info);

//! Flags an unimplemented routine if message_level >= 2
//! \param routine the name of the unimplemented routine
//! \param info the information message
void Unimplemented(std::string routine, std::string info);

//! Attaches a double to a string
//! \param front the initial string
//! \param back the appended double
std::string pastenum(std::string front, double back);

//! Attaches an integer to a string
//! \param front the initial string
//! \param back the appended integer
std::string pastenum(std::string front, int back);

} // end of namespace Msg

#endif
